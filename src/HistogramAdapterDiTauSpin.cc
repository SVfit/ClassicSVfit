#include "TauAnalysis/ClassicSVfit/interface/HistogramAdapterDiTauSpin.h"

namespace classic_svFit
{

//-------------------------------------------------------------------------------------------------
// auxiliary classes to reconstruct spin polarization vectors Bp and Bm and spin correlation matrix C
SVfitQuantityB_i::SVfitQuantityB_i(const std::string& label)
  : SVfitQuantityDiTau(label)
{}

TH1*
SVfitQuantityB_i::createHistogram(const MeasuredEvent& measuredEvent) const
{
  return HistogramTools::makeHistogram_linBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramB_i", 120, -3., +3.);
}

SVfitQuantityC_ij::SVfitQuantityC_ij(const std::string& label)
  : SVfitQuantityDiTau(label)
{}

TH1*
SVfitQuantityC_ij::createHistogram(const MeasuredEvent& measuredEvent) const
{
  return HistogramTools::makeHistogram_linBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramC_ij", 360, -9., +9.);
}

HistogramAdapterDiTauSpin::HistogramAdapterDiTauSpin(const std::string& label)
  : HistogramAdapterDiTau(label)
  , fittedTauPlus_(nullptr)
  , fittedTauMinus_(nullptr)
  , quantity_Bp_n_(nullptr)
  , quantity_Bp_r_(nullptr)
  , quantity_Bp_k_(nullptr)
  , quantity_Bm_n_(nullptr)
  , quantity_Bm_r_(nullptr)
  , quantity_Bm_k_(nullptr)
  , quantity_C_nn_(nullptr)
  , quantity_C_nr_(nullptr)
  , quantity_C_nk_(nullptr)
  , quantity_C_rn_(nullptr)
  , quantity_C_rr_(nullptr)
  , quantity_C_rk_(nullptr)
  , quantity_C_kn_(nullptr)
  , quantity_C_kr_(nullptr)
  , quantity_C_kk_(nullptr)
{
  quantity_Bp_n_ = new SVfitQuantityB_i(label_);
  quantities_.push_back(quantity_Bp_n_);
  quantity_Bp_r_ = new SVfitQuantityB_i(label_);
  quantities_.push_back(quantity_Bp_r_);
  quantity_Bp_k_ = new SVfitQuantityB_i(label_);
  quantities_.push_back(quantity_Bp_k_);

  quantity_Bm_n_ = new SVfitQuantityB_i(label_);
  quantities_.push_back(quantity_Bm_n_);
  quantity_Bm_r_ = new SVfitQuantityB_i(label_);
  quantities_.push_back(quantity_Bm_r_);
  quantity_Bm_k_ = new SVfitQuantityB_i(label_);
  quantities_.push_back(quantity_Bm_k_);

  quantity_C_nn_ = new SVfitQuantityC_ij(label_);
  quantities_.push_back(quantity_C_nn_);
  quantity_C_nr_ = new SVfitQuantityC_ij(label_);
  quantities_.push_back(quantity_C_nr_);
  quantity_C_nk_ = new SVfitQuantityC_ij(label_);
  quantities_.push_back(quantity_C_nk_);
  quantity_C_rn_ = new SVfitQuantityC_ij(label_);
  quantities_.push_back(quantity_C_rn_);
  quantity_C_rr_ = new SVfitQuantityC_ij(label_);
  quantities_.push_back(quantity_C_rr_);
  quantity_C_rk_ = new SVfitQuantityC_ij(label_);
  quantities_.push_back(quantity_C_rk_);
  quantity_C_kn_ = new SVfitQuantityC_ij(label_);
  quantities_.push_back(quantity_C_kn_);
  quantity_C_kr_ = new SVfitQuantityC_ij(label_);
  quantities_.push_back(quantity_C_kr_);
  quantity_C_kk_ = new SVfitQuantityC_ij(label_);
  quantities_.push_back(quantity_C_kk_);
}

HistogramAdapterDiTauSpin::~HistogramAdapterDiTauSpin()
{}

HistogramAdapter*
HistogramAdapterDiTauSpin::clone() const
{
  HistogramAdapterDiTauSpin* retVal = new HistogramAdapterDiTauSpin(this->label_);
  retVal->setMeasurement(this->measuredEvent_);
  retVal->bookHistograms(this->measuredEvent_);
  return retVal;
}

void
HistogramAdapterDiTauSpin::bookHistograms(const MeasuredEvent& measuredEvent)
{
  quantity_Bp_n_->bookHistogram(measuredEvent);
  quantity_Bp_r_->bookHistogram(measuredEvent);
  quantity_Bp_k_->bookHistogram(measuredEvent);
  quantity_Bm_n_->bookHistogram(measuredEvent);
  quantity_Bm_r_->bookHistogram(measuredEvent);
  quantity_Bm_k_->bookHistogram(measuredEvent);
  quantity_C_nn_->bookHistogram(measuredEvent);
  quantity_C_nr_->bookHistogram(measuredEvent);
  quantity_C_nk_->bookHistogram(measuredEvent);
  quantity_C_rn_->bookHistogram(measuredEvent);
  quantity_C_rr_->bookHistogram(measuredEvent);
  quantity_C_rk_->bookHistogram(measuredEvent);
  quantity_C_kn_->bookHistogram(measuredEvent);
  quantity_C_kr_->bookHistogram(measuredEvent);
  quantity_C_kk_->bookHistogram(measuredEvent);
}

void
HistogramAdapterDiTauSpin::setMeasurement(const MeasuredEvent& measuredEvent)
{
  HistogramAdapterDiTau::setMeasurement(measuredEvent);
}

void
HistogramAdapterDiTauSpin::setFittedTauLeptons(const FittedTauLepton& fittedTauLepton1, const FittedTauLepton& fittedTauLepton2)
{
  HistogramAdapterDiTau::setFittedTauLeptons(fittedTauLepton1, fittedTauLepton2);
  boostToHelicityFrame_.setFittedTauLeptons(fittedTauLepton1, fittedTauLepton2);
  std::vector<const FittedTauLepton*> fittedTauLeptons;
  fittedTauLeptons.push_back(&fittedTauLepton1);
  fittedTauLeptons.push_back(&fittedTauLepton2);
  fittedTauPlus_ = nullptr;
  fittedTauMinus_ = nullptr;
  for ( const FittedTauLepton* fittedTauLepton : fittedTauLeptons )
  {
    if      ( fittedTauLepton->getMeasuredTauLepton().charge() > 0 ) fittedTauPlus_ = fittedTauLepton;
    else if ( fittedTauLepton->getMeasuredTauLepton().charge() < 0 ) fittedTauMinus_ = fittedTauLepton;
  }
  assert(fittedTauPlus_ && fittedTauMinus_);
}

void
HistogramAdapterDiTauSpin::fillHistograms(double weight) const
{
  HistogramAdapterDiTau::fillHistograms(weight);
  Vector hPlus = polarimeterVector_(measuredEvent_.measuredTauPlus(), *fittedTauPlus_, BoostToHelicityFrame::kTauPlus, boostToHelicityFrame_);
  double hPlus_r = hPlus.x();
  double hPlus_n = hPlus.y();
  double hPlus_k = hPlus.z();
  Vector hMinus = polarimeterVector_(measuredEvent_.measuredTauMinus(), *fittedTauMinus_, BoostToHelicityFrame::kTauMinus, boostToHelicityFrame_);
  double hMinus_r = hMinus.x();
  double hMinus_n = hMinus.y();
  double hMinus_k = hMinus.z();
  quantity_Bp_n_->fillHistogram(hPlus_n);
  quantity_Bp_r_->fillHistogram(hPlus_r);
  quantity_Bp_k_->fillHistogram(hPlus_k);
  quantity_Bm_n_->fillHistogram(hMinus_n);
  quantity_Bm_r_->fillHistogram(hMinus_r);
  quantity_Bm_k_->fillHistogram(hMinus_k);
  quantity_C_nn_->fillHistogram(hPlus_n*hMinus_n);
  quantity_C_nr_->fillHistogram(hPlus_n*hMinus_r);
  quantity_C_nk_->fillHistogram(hPlus_n*hMinus_k);
  quantity_C_rn_->fillHistogram(hPlus_r*hMinus_n);
  quantity_C_rr_->fillHistogram(hPlus_r*hMinus_r);
  quantity_C_rk_->fillHistogram(hPlus_r*hMinus_k);
  quantity_C_kn_->fillHistogram(hPlus_k*hMinus_n);
  quantity_C_kr_->fillHistogram(hPlus_k*hMinus_r);
  quantity_C_kk_->fillHistogram(hPlus_k*hMinus_k);
}

double HistogramAdapterDiTauSpin::getBp_n() const
{
  return extractValue(quantity_Bp_n_);
}

double HistogramAdapterDiTauSpin::getBp_nErr() const
{
  return extractUncertainty(quantity_Bp_n_);
}

double HistogramAdapterDiTauSpin::getBp_r() const
{
  return extractValue(quantity_Bp_r_);
}

double HistogramAdapterDiTauSpin::getBp_rErr() const
{
  return extractUncertainty(quantity_Bp_r_);
}

double HistogramAdapterDiTauSpin::getBp_k() const
{
  return extractValue(quantity_Bp_k_);
}

double HistogramAdapterDiTauSpin::getBp_kErr() const
{
  return extractUncertainty(quantity_Bp_k_);
}

double HistogramAdapterDiTauSpin::getBm_n() const
{
  return extractValue(quantity_Bm_n_);
}

double HistogramAdapterDiTauSpin::getBm_nErr() const
{
  return extractUncertainty(quantity_Bm_n_);
}

double HistogramAdapterDiTauSpin::getBm_r() const
{
  return extractValue(quantity_Bm_r_);
}

double HistogramAdapterDiTauSpin::getBm_rErr() const
{
  return extractUncertainty(quantity_Bm_r_);
}

double HistogramAdapterDiTauSpin::getBm_k() const
{
  return extractValue(quantity_Bm_k_);
}

double HistogramAdapterDiTauSpin::getBm_kErr() const
{
  return extractUncertainty(quantity_Bm_k_);
}

double HistogramAdapterDiTauSpin::getC_nn() const
{
  return extractValue(quantity_C_nn_);
}

double HistogramAdapterDiTauSpin::getC_nnErr() const
{
  return extractUncertainty(quantity_C_nn_);
}

double HistogramAdapterDiTauSpin::getC_nr() const
{
  return extractValue(quantity_C_nr_);
}

double HistogramAdapterDiTauSpin::getC_nrErr() const
{
  return extractUncertainty(quantity_C_nr_);
}

double HistogramAdapterDiTauSpin::getC_nk() const
{
  return extractValue(quantity_C_nk_);
}

double HistogramAdapterDiTauSpin::getC_nkErr() const
{
  return extractUncertainty(quantity_C_nk_);
}

double HistogramAdapterDiTauSpin::getC_rn() const
{
  return extractValue(quantity_C_rn_);
}

double HistogramAdapterDiTauSpin::getC_rnErr() const
{
  return extractUncertainty(quantity_C_rn_);
}

double HistogramAdapterDiTauSpin::getC_rr() const
{
  return extractValue(quantity_C_rr_);
}

double HistogramAdapterDiTauSpin::getC_rrErr() const
{
  return extractUncertainty(quantity_C_rr_);
}

double HistogramAdapterDiTauSpin::getC_rk() const
{
  return extractValue(quantity_C_rk_);
}

double HistogramAdapterDiTauSpin::getC_rkErr() const
{
  return extractUncertainty(quantity_C_rk_);
}

double HistogramAdapterDiTauSpin::getC_kn() const
{
  return extractValue(quantity_C_kn_);
}

double HistogramAdapterDiTauSpin::getC_knErr() const
{
  return extractUncertainty(quantity_C_kn_);
}

double HistogramAdapterDiTauSpin::getC_kr() const
{
  return extractValue(quantity_C_kr_);
}

double HistogramAdapterDiTauSpin::getC_krErr() const
{
  return extractUncertainty(quantity_C_kr_);
}

double HistogramAdapterDiTauSpin::getC_kk() const
{
  return extractValue(quantity_C_kk_);
}

double HistogramAdapterDiTauSpin::getC_kkErr() const
{
  return extractUncertainty(quantity_C_kk_);
}

double HistogramAdapterDiTauSpin::DoEval(const double* x) const
{
  fillHistograms();
  return 0.;
}
//-------------------------------------------------------------------------------------------------

}
