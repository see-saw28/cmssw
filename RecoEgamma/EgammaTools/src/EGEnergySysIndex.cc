#include "RecoEgamma/EgammaTools/interface/EGEnergySysIndex.h"

namespace {
  std::array<std::string, EGEnergySysIndex::kNrSysErrs> makeEGEnergySysNames() {
    std::array<std::string, EGEnergySysIndex::kNrSysErrs> names;
    names[EGEnergySysIndex::kScaleUp] = "energyScaleUp";
    names[EGEnergySysIndex::kScaleDown] = "energyScaleDown";
    names[EGEnergySysIndex::kScaleUpValue] = "scaleUpValue";
    names[EGEnergySysIndex::kScaleDownValue] = "scaleDownValue";
    names[EGEnergySysIndex::kSmearUp] = "energySigmaUp";
    names[EGEnergySysIndex::kSmearDown] = "energySigmaDown";
    names[EGEnergySysIndex::kSmearUpValue] = "sigmaUpValue";
    names[EGEnergySysIndex::kSmearDownValue] = "sigmaDownValue";
    names[EGEnergySysIndex::kScaleValue] = "scaleValue";
    names[EGEnergySysIndex::kSmearValue] = "sigmaValue";
    names[EGEnergySysIndex::kSmearNrSigma] = "smearNrSigma";
    names[EGEnergySysIndex::kEcalPreCorr] = "ecalEnergyPreCorr";
    names[EGEnergySysIndex::kEcalErrPreCorr] = "ecalEnergyErrPreCorr";
    names[EGEnergySysIndex::kEcalPostCorr] = "ecalEnergyPostCorr";
    names[EGEnergySysIndex::kEcalErrPostCorr] = "ecalEnergyErrPostCorr";
    names[EGEnergySysIndex::kEcalTrkPreCorr] = "ecalTrkEnergyPreCorr";
    names[EGEnergySysIndex::kEcalTrkErrPreCorr] = "ecalTrkEnergyErrPreCorr";
    names[EGEnergySysIndex::kEcalTrkPostCorr] = "ecalTrkEnergyPostCorr";
    names[EGEnergySysIndex::kEcalTrkErrPostCorr] = "ecalTrkEnergyErrPostCorr";
    return names;
  }
}  // namespace

const std::array<std::string, EGEnergySysIndex::kNrSysErrs> EGEnergySysIndex::names_ = makeEGEnergySysNames();
