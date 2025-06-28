#include "RecoEgamma/EgammaTools/interface/ElectronEnergyCalibrator.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include <CLHEP/Random/RandGaussQ.h>

const EnergyScaleCorrection::ScaleCorrection ElectronEnergyCalibrator::defaultScaleCorr_;
const EnergyScaleCorrection::SmearCorrection ElectronEnergyCalibrator::defaultSmearCorr_;

ElectronEnergyCalibrator::ElectronEnergyCalibrator(const EpCombinationTool& combinator,
                                                   const std::string& correctionFile)
    : correctionRetriever_(correctionFile), epCombinationTool_(&combinator), rng_(nullptr), minEt_(1.0) {}

void ElectronEnergyCalibrator::initPrivateRng(TRandom* rnd) { rng_ = rnd; }

std::array<float, EGEnergySysIndex::kNrSysErrs> ElectronEnergyCalibrator::calibrate(
    reco::GsfElectron& ele,
    const unsigned int runNumber,
    const EcalRecHitCollection* recHits,
    edm::StreamID const& id,
    const ElectronEnergyCalibrator::EventType eventType) const {
  return calibrate(ele, runNumber, recHits, gauss(id), eventType);
}

std::array<float, EGEnergySysIndex::kNrSysErrs> ElectronEnergyCalibrator::calibrate(
    reco::GsfElectron& ele,
    unsigned int runNumber,
    const EcalRecHitCollection* recHits,
    const float smearNrSigma,
    const ElectronEnergyCalibrator::EventType eventType) const {
  const float scEta = ele.superCluster()->eta();
  const float eta = ele.eta();
  const float et = ele.ecalEnergy() / cosh(eta);

  if (et < minEt_ || edm::isNotFinite(et)) {
    std::array<float, EGEnergySysIndex::kNrSysErrs> retVal;
    retVal.fill(ele.energy());
    retVal[EGEnergySysIndex::kScaleValue] = 1.0;
    retVal[EGEnergySysIndex::kScaleUpValue] = 1.0;
    retVal[EGEnergySysIndex::kScaleDownValue] = 1.0;
    retVal[EGEnergySysIndex::kSmearValue] = 0.0;
    retVal[EGEnergySysIndex::kSmearUpValue] = 0.0;
    retVal[EGEnergySysIndex::kSmearDownValue] = 0.0;
    retVal[EGEnergySysIndex::kSmearNrSigma] = smearNrSigma;
    retVal[EGEnergySysIndex::kEcalPreCorr] = ele.ecalEnergy();
    retVal[EGEnergySysIndex::kEcalErrPreCorr] = ele.ecalEnergyError();
    retVal[EGEnergySysIndex::kEcalPostCorr] = ele.ecalEnergy();
    retVal[EGEnergySysIndex::kEcalErrPostCorr] = ele.ecalEnergyError();
    retVal[EGEnergySysIndex::kEcalTrkPreCorr] = ele.energy();
    retVal[EGEnergySysIndex::kEcalTrkErrPreCorr] = ele.corrections().combinedP4Error;
    retVal[EGEnergySysIndex::kEcalTrkPostCorr] = ele.energy();
    retVal[EGEnergySysIndex::kEcalTrkErrPostCorr] = ele.corrections().combinedP4Error;
    return retVal;
  }
  // Get the seed SC and its gain
  const DetId seedDetId = ele.superCluster()->seed()->seed();
  EcalRecHitCollection::const_iterator seedRecHit = recHits->find(seedDetId);
  unsigned int gainSeedSC = 12;
  if (seedRecHit != recHits->end()) {
    if (seedRecHit->checkFlag(EcalRecHit::kHasSwitchToGain6))
      gainSeedSC = 6;
    if (seedRecHit->checkFlag(EcalRecHit::kHasSwitchToGain1))
      gainSeedSC = 1;
  }

  const EnergyScaleCorrection::ScaleCorrection* scaleCorr =
      correctionRetriever_.getScaleCorr(runNumber, et, scEta, ele.full5x5_r9(), gainSeedSC);
  const EnergyScaleCorrection::SmearCorrection* smearCorr =
      correctionRetriever_.getSmearCorr(et, scEta, ele.full5x5_r9());
  if (scaleCorr == nullptr)
    scaleCorr = &defaultScaleCorr_;
  if (smearCorr == nullptr)
    smearCorr = &defaultSmearCorr_;

  std::array<float, EGEnergySysIndex::kNrSysErrs> uncertainties{};

  
  //MC central values are not scaled (scale = 1.0), data is not smeared (smearNrSigma = 0)
  //the smearing (or resolution extra parameter as it might better be called)
  //still has a second order effect on data as it enters the E/p combination as an adjustment
  //to the estimate of the resolution contained in caloEnergyError
  //MC gets all the scale systematics
  if (eventType == EventType::DATA) {
    setEnergyAndSystVarations(scaleCorr->scale(), 0., *scaleCorr, *smearCorr, ele, uncertainties);
  } else if (eventType == EventType::MC) {
    setEnergyAndSystVarations(1.0, smearNrSigma, *scaleCorr, *smearCorr, ele, uncertainties);
  }

  return uncertainties;
}

void ElectronEnergyCalibrator::setEnergyAndSystVarations(
    const float scale,
    const float smearNrSigma,
    const EnergyScaleCorrection::ScaleCorrection& scaleCorr,
    const EnergyScaleCorrection::SmearCorrection& smearCorr,
    reco::GsfElectron& ele,
    std::array<float, EGEnergySysIndex::kNrSysErrs>& energyData) const {

  const float smear = smearCorr.sigma();
  const float smearUp = smearCorr.sigmaUp();
  const float smearDn = smearCorr.sigmaDown();

  const float corr = scale + smear * smearNrSigma;
  const float corrSmearUp = scale + smearUp * smearNrSigma;
  const float corrSmearDn = scale + smearDn * smearNrSigma;
  
  const float corrScaleUp = smearCorr.scaleUp();
  const float corrScaleDn = smearCorr.scaleDown();

  const math::XYZTLorentzVector oldP4 = ele.p4();
  energyData[EGEnergySysIndex::kEcalTrkPreCorr] = ele.energy();
  energyData[EGEnergySysIndex::kEcalTrkErrPreCorr] = ele.corrections().combinedP4Error;
  energyData[EGEnergySysIndex::kEcalPreCorr] = ele.ecalEnergy();
  energyData[EGEnergySysIndex::kEcalErrPreCorr] = ele.ecalEnergyError();

  energyData[EGEnergySysIndex::kScaleUp] = calCombinedMom(ele, corrScaleUp, smear).first;
  energyData[EGEnergySysIndex::kScaleDown] = calCombinedMom(ele, corrScaleDn, smear).first;
  energyData[EGEnergySysIndex::kSmearUp] = calCombinedMom(ele, corrSmearUp, smearUp).first;
  energyData[EGEnergySysIndex::kSmearDown] = calCombinedMom(ele, corrSmearDn, smearDn).first;

  const std::pair<float, float> combinedMomentum = calCombinedMom(ele, corr, smear);
  setEcalEnergy(ele, corr, smear);
  const float energyCorr = combinedMomentum.first / oldP4.t();

  const math::XYZTLorentzVector newP4(
      oldP4.x() * energyCorr, oldP4.y() * energyCorr, oldP4.z() * energyCorr, combinedMomentum.first);

  ele.correctMomentum(newP4, ele.trackMomentumError(), combinedMomentum.second);
  energyData[EGEnergySysIndex::kEcalTrkPostCorr] = ele.energy();
  energyData[EGEnergySysIndex::kEcalTrkErrPostCorr] = ele.corrections().combinedP4Error;

  energyData[EGEnergySysIndex::kEcalPostCorr] = ele.ecalEnergy();
  energyData[EGEnergySysIndex::kEcalErrPostCorr] = ele.ecalEnergyError();
}

void ElectronEnergyCalibrator::setEcalEnergy(reco::GsfElectron& ele, const float scale, const float smear) const {
  const float oldEcalEnergy = ele.ecalEnergy();
  const float oldEcalEnergyErr = ele.ecalEnergyError();
  ele.setCorrectedEcalEnergy(oldEcalEnergy * scale);
  ele.setCorrectedEcalEnergyError(std::hypot(oldEcalEnergyErr * scale, oldEcalEnergy * smear * scale));
}

std::pair<float, float> ElectronEnergyCalibrator::calCombinedMom(reco::GsfElectron& ele,
                                                                 const float scale,
                                                                 const float smear) const {
  const float oldEcalEnergy = ele.ecalEnergy();
  const float oldEcalEnergyErr = ele.ecalEnergyError();

  const auto oldP4 = ele.p4();
  const float oldP4Err = ele.p4Error(reco::GsfElectron::P4_COMBINATION);
  const float oldTrkMomErr = ele.trackMomentumError();

  setEcalEnergy(ele, scale, smear);
  const auto& combinedMomentum = epCombinationTool_->combine(ele, oldEcalEnergyErr * scale);

  ele.setCorrectedEcalEnergy(oldEcalEnergy);
  ele.setCorrectedEcalEnergyError(oldEcalEnergyErr);
  ele.correctMomentum(oldP4, oldTrkMomErr, oldP4Err);

  return combinedMomentum;
}

double ElectronEnergyCalibrator::gauss(edm::StreamID const& id) const {
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
