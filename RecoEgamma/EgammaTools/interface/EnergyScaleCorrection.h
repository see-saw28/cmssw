#ifndef RecoEgamma_EgammaTools_EnergyScaleCorrection_h
#define RecoEgamma_EgammaTools_EnergyScaleCorrection_h

//author: Alan Smithee
//description:
//  A port of Shervin Nourbakhsh's EnergyScaleCorrection_class in EgammaAnalysis/ElectronTools
//  this reads the scale & smearing corrections in from a text file for given categories
//  it then allows these values to be accessed

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <bitset>
#include "correction.h"

class EnergyScaleCorrection {

public:
  class ScaleCorrection {
  public:
    ScaleCorrection() : scale_(1.), scaleErr_(0.), scaleUp_(1.), scaleDown_(1.) {}
    ScaleCorrection(float Scale, float ScaleErr, float ScaleUp, float ScaleDown)
        : scale_(Scale), scaleErr_(ScaleErr), scaleUp_(ScaleUp), scaleDown_(ScaleDown) {}

    float scale() const { return scale_; }
    float scaleErr() const { return scaleErr_; };
    float scaleUp() const { return scaleUp_; }
    float scaleDown() const { return scaleDown_; }

    

    friend std::ostream& operator<<(std::ostream& os, const ScaleCorrection& a) { return a.print(os); }
    std::ostream& print(std::ostream& os) const;

  private:
    float scale_, scaleErr_, scaleUp_, scaleDown_;

    };

  struct SmearCorrection {
  public:
    SmearCorrection() : smear_(0.), smearErr_(0.), smearUp_(0.), smearDown_(0.), scaleErr_(0.), scaleUp_(1.), scaleDown_(1.) {}
    SmearCorrection(float Smear, float SmearErr, float SmearUp, float SmearDown,
                    float ScaleErr, float ScaleUp, float ScaleDown)
        : smear_(Smear), smearErr_(SmearErr), smearUp_(SmearUp), smearDown_(SmearDown),
          scaleErr_(ScaleErr), scaleUp_(ScaleUp), scaleDown_(ScaleDown) {}

    friend std::ostream& operator<<(std::ostream& os, const SmearCorrection& a) { return a.print(os); }
    std::ostream& print(std::ostream& os) const;

    float sigma() const { return smear_; }
    float sigmaErr() const { return smearErr_; }
    float sigmaUp() const { return smearUp_; }
    float sigmaDown() const { return smearDown_ < 0 ? 0 : smearDown_; }
    float scaleErr() const { return scaleErr_; }
    float scaleUp() const { return scaleUp_; }
    float scaleDown() const { return scaleDown_; }

  private:
    float smear_, smearErr_, smearUp_, smearDown_;
    float scaleErr_, scaleUp_, scaleDown_;
  };



  EnergyScaleCorrection(const std::string& correctionFileName, unsigned int genSeed = 0);
  EnergyScaleCorrection() {}
  ~EnergyScaleCorrection() {}

  float scaleCorr(unsigned int runnr,
                  double et,
                  double ScEta,
                  double r9,
                  unsigned int gainSeed = 12
                ) const;

  float scaleCorrUncert(unsigned int runnr,
                        double et,
                        double ScEta,
                        double r9,
                        unsigned int gainSeed
                      ) const;

  float smearingSigma(
      double et, double ScEta, double r9) const;


  const ScaleCorrection* getScaleCorr(unsigned int runnr, double et, double ScEta, double r9, unsigned int gainSeed) const;
  const SmearCorrection* getSmearCorr(double et, double ScEta, double r9) const;

private:

  //static data members
  static constexpr float kDefaultScaleVal_ = 1.0;
  static constexpr float kDefaultSmearVal_ = 0.0;

  correction::CompoundCorrection::Ref corrScale;
  correction::Correction::Ref corrSmear;

  template <typename T1, typename T2>
  class Sorter {
  public:
    bool operator()(const std::pair<T1, T2>& lhs, const T1& rhs) const { return lhs.first < rhs; }
    bool operator()(const std::pair<T1, T2>& lhs, const std::pair<T1, T2>& rhs) const { return lhs.first < rhs.first; }
    bool operator()(const T1& lhs, const std::pair<T1, T2>& rhs) const { return lhs < rhs.first; }
    bool operator()(const T1& lhs, const T1& rhs) const { return lhs < rhs; }
  };
};

#endif
