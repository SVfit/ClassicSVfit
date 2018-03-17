#ifndef TauAnalysis_ClassicSVfit_svFitHistogramAdapter_h
#define TauAnalysis_ClassicSVfit_svFitHistogramAdapter_h

#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"

#include <Math/Functor.h>
#include <TH1.h>

namespace classic_svFit
{
  class HistogramTools
  {
   public:
    static TH1* compHistogramDensity(TH1 const* histogram);
    static void extractHistogramProperties(
        TH1 const* histogram,
        double& xMaximum,
        double& xMaximum_interpol,
        double& xMean,
        double& xQuantile016,
        double& xQuantile050,
        double& xQuantile084
    );
    static double extractValue(TH1 const* histogram);
    static double extractUncertainty(TH1 const* histogram);
    static double extractLmax(TH1 const* histogram);
    static TH1* makeHistogram_linBinWidth(const std::string& histogramName, int numBins, double xMin, double xMax);
    static TH1* makeHistogram_logBinWidth(const std::string& histogramName, double xMin, double xMax, double logBinWidth);
  };

  class SVfitQuantity
  {
   public:
    SVfitQuantity(const std::string& label);
    virtual ~SVfitQuantity();

    const TH1* getHistogram() const;
    void writeHistogram() const;

    void fillHistogram(double value);

    double extractValue() const;
    double extractUncertainty() const;
    double extractLmax() const;

    bool isValidSolution() const;

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
    virtual ~HistogramAdapter();

    void writeHistograms(const std::string& likelihoodFileName) const;

    double extractValue(const SVfitQuantity* quantity) const;
    double extractUncertainty(const SVfitQuantity* quantity) const;
    double extractLmax(const SVfitQuantity* quantity) const;

    bool isValidSolution() const;

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

    virtual TH1* createHistogram(const LorentzVector& visP4) const = 0;

    void bookHistogram(const LorentzVector& visP4);
  };

  class SVfitQuantityTauPt : public SVfitQuantityTau
  {
   public:
    SVfitQuantityTauPt(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& visP4) const;
  };

  class SVfitQuantityTauEta : public SVfitQuantityTau
  {
   public:
    SVfitQuantityTauEta(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& visP4) const;
  };

  class SVfitQuantityTauPhi : public SVfitQuantityTau
  {
   public:
    SVfitQuantityTauPhi(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& visP4) const;
  };
  
  class HistogramAdapterTau : public HistogramAdapter
  {
   public:
    HistogramAdapterTau(const std::string& label);

    void bookHistograms(const LorentzVector& visP4);

    void setMeasurement(const LorentzVector& visP4);
    void setTauP4(const LorentzVector& tauP4);

    void fillHistograms(const LorentzVector& tauP4, const LorentzVector& visP4) const;

    /// get pT, eta, phi, mass of tau lepton
    double getPt() const;
    double getPtErr() const;
    double getPtLmax() const;
    double getEta() const;
    double getEtaErr() const;
    double getEtaLmax() const;
    double getPhi() const;
    double getPhiErr() const;
    double getPhiLmax() const;

    /// convenient access to tau lepton four-vector
    LorentzVector getP4() const;

  private:
    double DoEval(const double* x) const;

   protected:
    LorentzVector visP4_;
    LorentzVector tauP4_;

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

    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const = 0;

    void bookHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met);
  };

  class SVfitQuantityDiTauPt : public SVfitQuantityDiTau
  {
   public:
    SVfitQuantityDiTauPt(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const;
  };

  class SVfitQuantityDiTauEta : public SVfitQuantityDiTau
  {
   public:
    SVfitQuantityDiTauEta(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const;
  };

  class SVfitQuantityDiTauPhi : public SVfitQuantityDiTau
  {
   public:
    SVfitQuantityDiTauPhi(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const;
  };

  class SVfitQuantityDiTauMass : public SVfitQuantityDiTau
  {
   public:
    SVfitQuantityDiTauMass(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const;
  };

  class SVfitQuantityDiTauTransverseMass : public SVfitQuantityDiTau
  {
   public:
    SVfitQuantityDiTauTransverseMass(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const;
  };

  class HistogramAdapterDiTau : public HistogramAdapter
  {
   public:
    HistogramAdapterDiTau(const std::string& label = "ditau");
    ~HistogramAdapterDiTau();

    void bookHistograms(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met);

    void setMeasurement(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met);
    void setTau1And2P4(const LorentzVector& tau1P4,  const LorentzVector& tau2P4);

    void fillHistograms(const LorentzVector& tau1P4, const LorentzVector& tau2P4, const LorentzVector& ditauP4,
			const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const;

    HistogramAdapterTau* tau1() const;
    HistogramAdapterTau* tau2() const;

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

    /// convenient access to four-vector of di-tau system 
    LorentzVector getP4() const;

  private:
    double DoEval(const double* x) const;

   protected:
    LorentzVector vis1P4_;
    LorentzVector vis2P4_;
    Vector met_;
    LorentzVector tau1P4_;
    LorentzVector tau2P4_;
    LorentzVector ditauP4_;

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
