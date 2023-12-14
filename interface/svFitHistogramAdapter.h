#ifndef TauAnalysis_ClassicSVfit_svFitHistogramAdapter_h
#define TauAnalysis_ClassicSVfit_svFitHistogramAdapter_h

#include "TauAnalysis/ClassicSVfit/interface/FittedTauLepton.h"   // FittedTauLepton
#include "TauAnalysis/ClassicSVfit/interface/MeasuredEvent.h"     // MeasuredEvent
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h" // MeasuredTauLepton

#include <Math/Functor.h>                                         // ROOT::Math::Functor
#include <TH1.h>                                                  // TH1

#include <string>                                                 // std::string
#include <vector>                                                 // std::vector<>

namespace classic_svFit
{
  class HistogramTools
  {
   public:
    static
    TH1*
    compHistogramDensity(TH1 const* histogram);
    static
    void
    extractHistogramProperties(
        TH1 const* histogram,
        double& xMaximum,
        double& xMaximum_interpol,
        double& xMean,
        double& xQuantile016,
        double& xQuantile050,
        double& xQuantile084
    );
    static
    double
    extractValue(TH1 const* histogram);
    static
    double
    extractUncertainty(TH1 const* histogram);
    static
    double
    extractLmax(TH1 const* histogram);
    static
    TH1*
    makeHistogram_linBinWidth(const std::string& histogramName, int numBins, double xMin, double xMax);
    static
    TH1*
    makeHistogram_logBinWidth(const std::string& histogramName, double xMin, double xMax, double logBinWidth);
  };

  class SVfitQuantity
  {
   public:
    SVfitQuantity(const std::string& label);
    virtual
    ~SVfitQuantity();

    const TH1*
    getHistogram() const;
    void
    writeHistogram() const;

    void
    fillHistogram(double value, double weight = 1.);

    double
    extractValue() const;
    double
    extractUncertainty() const;
    double
    extractLmax() const;

    bool
    isValidSolution() const;

   protected:
    std::string label_;

    mutable TH1* histogram_ = nullptr;

   private:
    static int nInstances;
   protected:
    std::string uniqueName_;
  };

  class HistogramAdapter : public ROOT::Math::Functor
  {
   public:
    HistogramAdapter(const std::string& label);
    virtual
    ~HistogramAdapter();

    virtual
    HistogramAdapter* clone() const = 0;

    void
    writeHistograms(const std::string& likelihoodFileName) const;

    double
    extractValue(const SVfitQuantity* quantity) const;
    double
    extractUncertainty(const SVfitQuantity* quantity) const;
    double
    extractLmax(const SVfitQuantity* quantity) const;

    bool
    isValidSolution() const;

   protected:
    std::string label_;

    mutable std::vector<SVfitQuantity*> quantities_;
  };

  //-------------------------------------------------------------------------------------------------
  // auxiliary classes to reconstruct pT, eta, phi of single tau leptons
  class SVfitQuantityTau : public SVfitQuantity
  {
   public:
    SVfitQuantityTau(const std::string& label);

    virtual
    TH1*
    createHistogram(const MeasuredTauLepton& measuredTauLepton) const = 0;

    void
    bookHistogram(const MeasuredTauLepton& measuredTauLepton);
  };

  class SVfitQuantityTauPt : public SVfitQuantityTau
  {
   public:
    SVfitQuantityTauPt(const std::string& label);
    virtual
    TH1*
    createHistogram(const MeasuredTauLepton& measuredTauLepton) const;
  };

  class SVfitQuantityTauEta : public SVfitQuantityTau
  {
   public:
    SVfitQuantityTauEta(const std::string& label);
    virtual
    TH1*
    createHistogram(const MeasuredTauLepton& measuredTauLepton) const;
  };

  class SVfitQuantityTauPhi : public SVfitQuantityTau
  {
   public:
    SVfitQuantityTauPhi(const std::string& label);
    virtual
    TH1*
    createHistogram(const MeasuredTauLepton& measuredTauLepton) const;
  };
  
  class HistogramAdapterTau : public HistogramAdapter
  {
   public:
    HistogramAdapterTau(const std::string& label);

    HistogramAdapter*
    clone() const;

    void
    bookHistograms(const MeasuredTauLepton& measuredTauLepton);

    void
    setMeasurement(const MeasuredTauLepton& measuredTauLepton);
    void
    setFittedTauLepton(const FittedTauLepton& fittedTauLepton);

    void
    fillHistograms(double weight = 1.) const;

    /// get pT, eta, phi, mass of tau lepton
    double
    getPt() const;
    double
    getPtErr() const;
    double
    getPtLmax() const;
    double
    getEta() const;
    double
    getEtaErr() const;
    double
    getEtaLmax() const;
    double
    getPhi() const;
    double
    getPhiErr() const;
    double
    getPhiLmax() const;

    /// convenient access to tau lepton four-vector
    LorentzVector
    getP4() const;

  private:
    double
    DoEval(const double* x) const;

   protected:
    MeasuredTauLepton measuredTauLepton_;
    LorentzVector fittedTauP4_;

    SVfitQuantityTauPt* quantity_pt_;
    SVfitQuantityTauEta* quantity_eta_;
    SVfitQuantityTauPhi* quantity_phi_;
  };
  //-------------------------------------------------------------------------------------------------

  //-------------------------------------------------------------------------------------------------
  // auxiliary classes to reconstruct pT, eta, phi, mass, and transverse mass of tau lepton pairs
  class SVfitQuantityDiTau : public SVfitQuantity
  {
   public:
    SVfitQuantityDiTau(const std::string& label);

    virtual
    TH1*
    createHistogram(const MeasuredEvent& measuredEvent) const = 0;

    void
    bookHistogram(const MeasuredEvent& measuredEvent);
  };

  class SVfitQuantityDiTauPt : public SVfitQuantityDiTau
  {
   public:
    SVfitQuantityDiTauPt(const std::string& label);
    virtual
    TH1*
    createHistogram(const MeasuredEvent& measuredEvent) const;
  };

  class SVfitQuantityDiTauEta : public SVfitQuantityDiTau
  {
   public:
    SVfitQuantityDiTauEta(const std::string& label);
    virtual
    TH1*
    createHistogram(const MeasuredEvent& measuredEvent) const;
  };

  class SVfitQuantityDiTauPhi : public SVfitQuantityDiTau
  {
   public:
    SVfitQuantityDiTauPhi(const std::string& label);
    virtual
    TH1*
    createHistogram(const MeasuredEvent& measuredEvent) const;
  };

  class SVfitQuantityDiTauMass : public SVfitQuantityDiTau
  {
   public:
    SVfitQuantityDiTauMass(const std::string& label);
    virtual
    TH1*
    createHistogram(const MeasuredEvent& measuredEvent) const;
  };

  class SVfitQuantityDiTauTransverseMass : public SVfitQuantityDiTau
  {
   public:
    SVfitQuantityDiTauTransverseMass(const std::string& label);
    virtual
    TH1*
    createHistogram(const MeasuredEvent& measuredEvent) const;
  };

  class HistogramAdapterDiTau : public HistogramAdapter
  {
   public:
    HistogramAdapterDiTau(const std::string& label = "ditau");
    virtual
    ~HistogramAdapterDiTau();

    HistogramAdapter*
    clone() const;

    virtual
    void
    bookHistograms(const MeasuredEvent& measuredEvent);

    virtual
    void
    setMeasurement(const MeasuredEvent& measuredEvent);
    virtual
    void
    setFittedTauLeptons(const FittedTauLepton& fittedTauLepton1, const FittedTauLepton& fittedTauLepton2);

    virtual
    void
    fillHistograms(double weight = 1.) const;

    HistogramAdapterTau*
    tau1() const;
    HistogramAdapterTau*
    tau2() const;

    /// get pT, eta, phi, mass and transverse mass of di-tau system
    double
    getPt() const;
    double
    getPtErr() const;
    double
    getPtLmax() const;
    double
    getEta() const;
    double
    getEtaErr() const;
    double
    getEtaLmax() const;
    double
    getPhi() const;
    double
    getPhiErr() const;
    double
    getPhiLmax() const;
    double
    getMass() const;
    double
    getMassErr() const;
    double
    getMassLmax() const;
    double
    getTransverseMass() const;
    double
    getTransverseMassErr() const;
    double
    getTransverseMassLmax() const;

    /// convenient access to four-vector of di-tau system 
    LorentzVector
    getP4() const;

   private:
    double
    DoEval(const double* x) const;

   protected:
    MeasuredEvent measuredEvent_;

    LorentzVector fittedTau1P4_;
    LorentzVector fittedTau2P4_;
    LorentzVector fittedDiTauP4_;

    SVfitQuantityDiTauPt* quantity_pt_;
    SVfitQuantityDiTauEta* quantity_eta_;
    SVfitQuantityDiTauPhi* quantity_phi_;
    SVfitQuantityDiTauMass* quantity_mass_;
    SVfitQuantityDiTauTransverseMass* quantity_transverseMass_;

    HistogramAdapterTau* adapter_tau1_;
    HistogramAdapterTau* adapter_tau2_;
  };
  //-------------------------------------------------------------------------------------------------
}

#endif
