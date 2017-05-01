#ifndef TauAnalysis_ClassicSVfit_svFitHistogramAdapter_h
#define TauAnalysis_ClassicSVfit_svFitHistogramAdapter_h

#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"

#include <Math/Functor.h>
#include <TH1.h>

namespace classic_svFit
{

  class SVfitQuantity
  {
   public:
    SVfitQuantity();
    virtual ~SVfitQuantity();

    virtual TH1* CreateHistogram(std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET) const = 0;
    virtual double FitFunction(std::vector<classic_svFit::LorentzVector> const& fittedTauLeptons, std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET) const = 0;

    void SetHistogram(std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET);
    void WriteHistograms() const;

    double Eval(
        std::vector<classic_svFit::LorentzVector> const& fittedTauLeptons,
        std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons,
        classic_svFit::Vector const& measuredMET
    ) const;

    double ExtractValue() const;
    double ExtractUncertainty() const;
    double ExtractLmax() const;

    mutable TH1* histogram_ = nullptr;
  };

  class HiggsPtSVfitQuantity : public SVfitQuantity
  {
   public:
    virtual TH1* CreateHistogram(std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET) const;
    virtual double FitFunction(std::vector<classic_svFit::LorentzVector> const& fittedTauLeptons, std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET) const;
  };
  
  class HiggsEtaSVfitQuantity : public SVfitQuantity
  {
   public:
    virtual TH1* CreateHistogram(std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET) const;
    virtual double FitFunction(std::vector<classic_svFit::LorentzVector> const& fittedTauLeptons, std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET) const;
  };
  
  class HiggsPhiSVfitQuantity : public SVfitQuantity
  {
   public:
    virtual TH1* CreateHistogram(std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET) const;
    virtual double FitFunction(std::vector<classic_svFit::LorentzVector> const& fittedTauLeptons, std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET) const;
  };
  
  class HiggsMassSVfitQuantity : public SVfitQuantity
  {
   public:
    virtual TH1* CreateHistogram(std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET) const;
    virtual double FitFunction(std::vector<classic_svFit::LorentzVector> const& fittedTauLeptons, std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET) const;
  };
  
  class TransverseMassSVfitQuantity : public SVfitQuantity
  {
   public:
    virtual TH1* CreateHistogram(std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET) const;
    virtual double FitFunction(std::vector<classic_svFit::LorentzVector> const& fittedTauLeptons, std::vector<classic_svFit::LorentzVector> const& measuredTauLeptons, classic_svFit::Vector const& measuredMET) const;
  };
  
  class HistogramAdapter : public ROOT::Math::Functor
  {
   public:
    HistogramAdapter();
    ~HistogramAdapter();

    void setMeasurement(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met);
    void setTau1P4(const LorentzVector& tau1P4);
    void setTau2P4(const LorentzVector& tau2P4);

    void bookHistograms(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met);

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
    void fillHistograms(const LorentzVector& tau1P4, const LorentzVector& tau2P4,
                        const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const;

   protected:
    LorentzVector vis1P4_;
    LorentzVector vis2P4_;
    Vector met_;
    
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
