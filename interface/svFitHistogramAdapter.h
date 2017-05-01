#ifndef TauAnalysis_ClassicSVfit_svFitHistogramAdapter_h
#define TauAnalysis_ClassicSVfit_svFitHistogramAdapter_h

#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"

#include <Math/Functor.h>
#include <TH1.h>

namespace classic_svFit
{
  class HistogramAdapter : public ROOT::Math::Functor
  {
   public:
    HistogramAdapter();
    ~HistogramAdapter();

    void setTau1P4(const LorentzVector& tau1P4);
    void setTau2P4(const LorentzVector& tau2P4);

    void bookHistograms(const LorentzVector& vis1P4, const LorentzVector& vis2P4);
    void fillHistograms(const LorentzVector& tau1P4, const LorentzVector& tau2P4) const;

    /// get pT, eta, phi, mass and transverse mass of di-tau system
    double getPt() const;
    double getPtErr() const;
    double getPtLmax() const;
    double getEta() const;
    double getEtaErr() const;
    double getEtaLmax() const;
    double getPhi() const;
    double getPhiErr() const;
    double getPhiLmax() const;
    double getMass() const;
    double getMassErr() const;
    double getMassLmax() const;
    double getTransverseMass() const;
    double getTransverseMassErr() const;
    double getTransverseMassLmax() const;

    void writeHistograms(const std::string& likelihoodFileName) const;

   private:
    virtual double DoEval(const double* x) const;

   protected:
    LorentzVector tau1P4_;
    LorentzVector tau2P4_;

    mutable TH1* histogramPt_;
    mutable TH1* histogramEta_;
    mutable TH1* histogramPhi_;
    mutable TH1* histogramMass_;
    mutable TH1* histogramTransverseMass_;
  };
}

#endif
