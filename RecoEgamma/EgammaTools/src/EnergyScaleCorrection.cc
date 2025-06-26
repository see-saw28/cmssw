#include "RecoEgamma/EgammaTools/interface/EnergyScaleCorrection.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <iterator>
#include "correction.h"

EnergyScaleCorrection::EnergyScaleCorrection(const std::string& correctionFileName,
  const std::string& correctionSetScale,
  const std::string& correctionSetSmear,
  unsigned int genSeed)
{
if (!correctionFileName.empty()) {
auto cset = correction::CorrectionSet::from_file(correctionFileName);

corrScale = cset->compound().at(correctionSetScale);
if (!corrScale) {
throw cms::Exception("EnergyScaleCorrection") << "scale correction set not found: " << correctionSetScale;
}

corrSmear = cset->at(correctionSetSmear);
if (!corrSmear) {
throw cms::Exception("EnergyScaleCorrection") << "smearing correction set not found: " << correctionSetSmear;
}
}
}

float EnergyScaleCorrection::scaleCorr(unsigned int runNumber,
                                       double et,
                                       double eta,
                                       double r9,
                                       unsigned int gainSeed
                                      ) const {
  const ScaleCorrection* scaleCorr = getScaleCorr(runNumber, et, eta, r9, gainSeed);
  if (scaleCorr != nullptr)
    return scaleCorr->scale();
  else
    return kDefaultScaleVal_;
}

float EnergyScaleCorrection::scaleCorrUncert(unsigned int runNumber,
                                             double et,
                                             double eta,
                                             double r9,
                                             unsigned int gainSeed
                                            ) const {
  const ScaleCorrection* scaleCorr = getScaleCorr(runNumber, et, eta, r9, gainSeed);
  if (scaleCorr != nullptr)
    return scaleCorr->scaleErr();
  else
    return 0.;
}


float EnergyScaleCorrection::smearingSigma(
    double et, double eta, double r9) const {
  const SmearCorrection* smearCorr = getSmearCorr(et, eta, r9);

  if (smearCorr != nullptr)
    return smearCorr->sigma();
  else
    return kDefaultSmearVal_;
}

const EnergyScaleCorrection::ScaleCorrection* EnergyScaleCorrection::getScaleCorr(
    unsigned int runnr, double et, double ScEta, double r9, unsigned int gainSeed) const {
    
    // build the input vector of the evaluator
    std::vector<correction::Variable::Type> input_vector = {
      std::string("scale"), static_cast<double>(runnr), ScEta, r9, std::abs(ScEta), et, static_cast<double>(gainSeed)
    };

    double scaleValue = corrScale->evaluate(input_vector);
    input_vector[0] = "escale";
    double escaleValue = corrScale->evaluate(input_vector);
    input_vector[0] = "scale_up";
    double scaleUpValue = corrScale->evaluate(input_vector);
    input_vector[0] = "scale_down";
    double scaleDownValue = corrScale->evaluate(input_vector);

    return new ScaleCorrection(scaleValue, escaleValue, scaleUpValue, scaleDownValue);
    
}

const EnergyScaleCorrection::SmearCorrection* EnergyScaleCorrection::getSmearCorr(
    double et, double ScEta, double r9) const {
  
  // build the input vector of the evaluator
  std::vector<correction::Variable::Type> input_vector = {
    std::string("smear"), et, r9, std::abs(ScEta)
  };

  double smearValue = corrSmear->evaluate(input_vector);
  input_vector[0] = "esmear";
  double esmearValue = corrSmear->evaluate(input_vector);
  input_vector[0] = "smear_up";
  double smearUpValue = corrSmear->evaluate(input_vector);
  input_vector[0] = "smear_down";
  double smearDownValue = corrSmear->evaluate(input_vector);
  input_vector[0] = "escale";
  double escaleValue = corrSmear->evaluate(input_vector);
  input_vector[0] = "scale_up";
  double scaleUpValue = corrSmear->evaluate(input_vector);
  input_vector[0] = "scale_down";
  double scaleDownValue = corrSmear->evaluate(input_vector);

  return new SmearCorrection(smearValue, esmearValue, 
                             smearUpValue, smearDownValue,
                             escaleValue, scaleUpValue, scaleDownValue);
}




std::ostream& EnergyScaleCorrection::SmearCorrection::print(std::ostream& os) const {
  os << smear_ << " +/- " << smearErr_ << "\t" << smearUp_ << " +/- " << smearDown_ << "\t" << scaleErr_ << " +/- " << scaleUp_ << " +/- " << scaleDown_;
  return os;
}

