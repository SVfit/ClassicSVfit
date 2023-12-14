#ifndef TauAnalysis_ClassicSVfit_HistogramAdapterDiTauSpin_h
#define TauAnalysis_ClassicSVfit_HistogramAdapterDiTauSpin_h

#include "TauAnalysis/ClassicSVfit/interface/BoostToHelicityFrame.h"  // BoostToHelicityFrame
#include "TauAnalysis/ClassicSVfit/interface/FittedTauLepton.h"       // FittedTauLepton
#include "TauAnalysis/ClassicSVfit/interface/HistogramAdapterDiTau.h" // HistogramAdapterDiTau
#include "TauAnalysis/ClassicSVfit/interface/PolarimeterVector.h"     // PolarimeterVector

namespace classic_svFit
{
  //-------------------------------------------------------------------------------------------------
  // auxiliary classes to reconstruct spin polarization vectors Bp and Bm and spin correlation matrix C
  class SVfitQuantityB_i : public SVfitQuantityDiTau
  {
   public:
    SVfitQuantityB_i(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const;
  };

  class SVfitQuantityC_ij : public SVfitQuantityDiTau
  {
   public:
    SVfitQuantityC_ij(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const;
  };
  //-------------------------------------------------------------------------------------------------

  //-------------------------------------------------------------------------------------------------
  class HistogramAdapterDiTauSpin : public HistogramAdapterDiTau
  {
   public:
    HistogramAdapterDiTauSpin(const std::string& label = "ditau");
    ~HistogramAdapterDiTauSpin();

    HistogramAdapter*
    clone() const;

    void
    bookHistograms(const MeasuredEvent& measuredEvent);

    void
    setMeasurement(const MeasuredEvent& measuredEvent);

    void
    setFittedTauLeptons(const FittedTauLepton& fittedTauLepton1, const FittedTauLepton& fittedTauLepton2);

    void
    fillHistograms(double weight) const;

    /// get elements of spin polarization vectors Bp and Bm and of spin correlation matrix C
    double
    getBp_n() const;
    double
    getBp_nErr() const;
    double
    getBp_r() const;
    double
    getBp_rErr() const;
    double
    getBp_k() const;
    double
    getBp_kErr() const;

    double
    getBm_n() const;
    double
    getBm_nErr() const;
    double
    getBm_r() const;
    double
    getBm_rErr() const;
    double
    getBm_k() const;
    double
    getBm_kErr() const;

    double
    getC_nn() const;
    double
    getC_nnErr() const;
    double
    getC_nr() const;
    double
    getC_nrErr() const;
    double
    getC_nk() const;
    double
    getC_nkErr() const;
    double
    getC_rn() const;
    double
    getC_rnErr() const;
    double
    getC_rr() const;
    double
    getC_rrErr() const;
    double
    getC_rk() const;
    double
    getC_rkErr() const;
    double
    getC_kn() const;
    double
    getC_knErr() const;
    double
    getC_kr() const;
    double
    getC_krErr() const;
    double
    getC_kk() const;
    double
    getC_kkErr() const;
   
   private:
    double
    DoEval(const double* x) const;

   protected:
    FittedTauLepton* fittedTauPlus_;
    FittedTauLepton* fittedTauMinus_;

    BoostToHelicityFrame boostToHelicityFrame_;
    PolarimeterVector polarimeterVector_;

    SVfitQuantityB_i* quantity_Bp_n_;
    SVfitQuantityB_i* quantity_Bp_r_;
    SVfitQuantityB_i* quantity_Bp_k_;

    SVfitQuantityB_i* quantity_Bm_n_;
    SVfitQuantityB_i* quantity_Bm_r_;
    SVfitQuantityB_i* quantity_Bm_k_;

    SVfitQuantityC_ij* quantity_C_nn_;
    SVfitQuantityC_ij* quantity_C_nr_;
    SVfitQuantityC_ij* quantity_C_nk_;
    SVfitQuantityC_ij* quantity_C_rn_;
    SVfitQuantityC_ij* quantity_C_rr_;
    SVfitQuantityC_ij* quantity_C_rk_;
    SVfitQuantityC_ij* quantity_C_kn_;
    SVfitQuantityC_ij* quantity_C_kr_;
    SVfitQuantityC_ij* quantity_C_kk_;
  };
  //-------------------------------------------------------------------------------------------------
}
