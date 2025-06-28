#include "RecoEgamma/EgammaTools/interface/PhotonEnergyCalibrator.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include <CLHEP/Random/RandGaussQ.h>

const EnergyScaleCorrection::ScaleCorrection PhotonEnergyCalibrator::defaultScaleCorr_;
const EnergyScaleCorrection::SmearCorrection PhotonEnergyCalibrator::defaultSmearCorr_;

PhotonEnergyCalibrator::PhotonEnergyCalibrator(const std::string& correctionFile)
    : correctionRetriever_(correctionFile), rng_(nullptr), minEt_(1.0) {}

void PhotonEnergyCalibrator::initPrivateRng(TRandom* rnd) { rng_ = rnd; }

std::array<float, EGEnergySysIndex::kNrSysErrs> PhotonEnergyCalibrator::calibrate(
    reco::Photon& photon,
    const unsigned int runNumber,
    const EcalRecHitCollection* recHits,
    edm::StreamID const& id,
    const PhotonEnergyCalibrator::EventType eventType) const {
  return calibrate(photon, runNumber, recHits, gauss(id), eventType);
}

std::array<float, EGEnergySysIndex::kNrSysErrs> PhotonEnergyCalibrator::calibrate(
    reco::Photon& photon,
    const unsigned int runNumber,
    const EcalRecHitCollection* recHits,
    const float smearNrSigma,
    const PhotonEnergyCalibrator::EventType eventType) const {
  const float scEta = photon.superCluster()->eta();
  const float eta = photon.eta();
  const float et = photon.getCorrectedEnergy(reco::Photon::P4type::regression2) / cosh(eta);

  if (et < minEt_ || edm::isNotFinite(et)) {
    std::array<float, EGEnergySysIndex::kNrSysErrs> retVal;
    retVal.fill(photon.getCorrectedEnergy(reco::Photon::P4type::regression2));
    retVal[EGEnergySysIndex::kScaleValue] = 1.0;
    retVal[EGEnergySysIndex::kScaleUpValue] = 1.0;
    retVal[EGEnergySysIndex::kScaleDownValue] = 1.0;
    retVal[EGEnergySysIndex::kSmearValue] = 0.0;
    retVal[EGEnergySysIndex::kSmearUpValue] = 0.0;
    retVal[EGEnergySysIndex::kSmearDownValue] = 0.0;
    retVal[EGEnergySysIndex::kSmearNrSigma] = smearNrSigma;
    retVal[EGEnergySysIndex::kEcalErrPreCorr] = photon.getCorrectedEnergyError(reco::Photon::P4type::regression2);
    retVal[EGEnergySysIndex::kEcalErrPostCorr] = photon.getCorrectedEnergyError(reco::Photon::P4type::regression2);
    retVal[EGEnergySysIndex::kEcalTrkPreCorr] = 0.;
    retVal[EGEnergySysIndex::kEcalTrkErrPreCorr] = 0.;
    retVal[EGEnergySysIndex::kEcalTrkPostCorr] = 0.;
    retVal[EGEnergySysIndex::kEcalTrkErrPostCorr] = 0.;

    return retVal;
  }
  // Get the seed SC and its gain
  const DetId seedDetId = photon.superCluster()->seed()->seed();
  EcalRecHitCollection::const_iterator seedRecHit = recHits->find(seedDetId);
  unsigned int gainSeedSC = 12;
  if (seedRecHit != recHits->end()) {
    if (seedRecHit->checkFlag(EcalRecHit::kHasSwitchToGain6))
      gainSeedSC = 6;
    if (seedRecHit->checkFlag(EcalRecHit::kHasSwitchToGain1))
      gainSeedSC = 1;
  }

  const EnergyScaleCorrection::ScaleCorrection* scaleCorr =
      correctionRetriever_.getScaleCorr(runNumber, et, scEta, photon.full5x5_r9(), gainSeedSC);
  const EnergyScaleCorrection::SmearCorrection* smearCorr =
      correctionRetriever_.getSmearCorr(et, scEta, photon.full5x5_r9());
  if (scaleCorr == nullptr)
    scaleCorr = &defaultScaleCorr_;
  if (smearCorr == nullptr)
    smearCorr = &defaultSmearCorr_;

  std::array<float, EGEnergySysIndex::kNrSysErrs> uncertainties{};


  //MC central values are not scaled (scale = 1.0), data is not smeared (smearNrSigma = 0)
  //smearing still has a second order effect on data as it enters the E/p combination as an
  //extra uncertainty on the calo energy
  //MC gets all the scale systematics
  if (eventType == EventType::DATA) {
    setEnergyAndSystVarations(scaleCorr->scale(), 0., *scaleCorr, *smearCorr, photon, uncertainties);
  } else if (eventType == EventType::MC) {
    setEnergyAndSystVarations(1.0, smearNrSigma, *scaleCorr, *smearCorr, photon, uncertainties);
  }

  return uncertainties;
}

void PhotonEnergyCalibrator::setEnergyAndSystVarations(
    const float scale,
    const float smearNrSigma,
    const EnergyScaleCorrection::ScaleCorrection& scaleCorr,
    const EnergyScaleCorrection::SmearCorrection& smearCorr,
    reco::Photon& photon,
    std::array<float, EGEnergySysIndex::kNrSysErrs>& energyData) const {

  const float smear = smearCorr.sigma();
  const float smearUp = smearCorr.sigmaUp();
  const float smearDn = smearCorr.sigmaDown();

  const float corr = scale + smear * smearNrSigma;
  const float corrSmearUp = scale + smearUp * smearNrSigma;
  const float corrSmearDn = scale + smearDn * smearNrSigma;
  
  const float corrScaleUp = smearCorr.scaleUp();
  const float corrScaleDn = smearCorr.scaleDown();

  const double oldEcalEnergy = photon.getCorrectedEnergy(reco::Photon::P4type::regression2);
  const double oldEcalEnergyError = photon.getCorrectedEnergyError(reco::Photon::P4type::regression2);

  energyData[EGEnergySysIndex::kEcalPreCorr] = oldEcalEnergy;
  energyData[EGEnergySysIndex::kEcalErrPreCorr] = oldEcalEnergyError;

  const double newEcalEnergy = oldEcalEnergy * corr;
  const double newEcalEnergyError = std::hypot(oldEcalEnergyError * corr, smear * newEcalEnergy);
  photon.setCorrectedEnergy(reco::Photon::P4type::regression2, newEcalEnergy, newEcalEnergyError, true);

  // The total variation
  energyData[EGEnergySysIndex::kScaleUp] = oldEcalEnergy * corrScaleUp;
  energyData[EGEnergySysIndex::kScaleDown] = oldEcalEnergy * corrScaleDn;
  energyData[EGEnergySysIndex::kSmearUp] = oldEcalEnergy * corrSmearUp;
  energyData[EGEnergySysIndex::kSmearDown] = oldEcalEnergy * corrSmearDn;


  energyData[EGEnergySysIndex::kEcalPostCorr] = photon.getCorrectedEnergy(reco::Photon::P4type::regression2);
  energyData[EGEnergySysIndex::kEcalErrPostCorr] = photon.getCorrectedEnergyError(reco::Photon::P4type::regression2);
}

double PhotonEnergyCalibrator::gauss(edm::StreamID const& id) const {
  if (rng_) {
    return rng_->Gaus();
  } else {
    edm::Service<edm::RandomNumberGenerator> rng;
    if (!rng.isAvailable()) {
      throw cms::Exception("Configuration")
          << "XXXXXXX requires the RandomNumberGeneratorService\n"
             "which is not present in the configuration file.  You must add the service\n"
             "in the configuration file or remove the modules that require it.";
    }
    CLHEP::RandGaussQ gaussDistribution(rng->getEngine(id), 0.0, 1.0);
    return gaussDistribution.fire();
  }
}
