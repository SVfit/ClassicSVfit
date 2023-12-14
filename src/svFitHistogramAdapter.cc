#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include <TMath.h>
#include <TFile.h>
#include <TObject.h>
#include <TLorentzVector.h>

#include <numeric>

#include <boost/algorithm/string/replace.hpp>

namespace classic_svFit
{

TH1*
HistogramTools::compHistogramDensity(TH1 const* histogram)
{
  TH1* histogram_density = static_cast<TH1*>(histogram->Clone((std::string(histogram->GetName()) + "_density").c_str()));
  histogram_density->Scale(1.0, "width");
  return histogram_density;
}

void HistogramTools::extractHistogramProperties(
    TH1 const* histogram,
    double& xMaximum,
    double& xMaximum_interpol,
    double& xMean,
    double& xQuantile016,
    double& xQuantile050,
    double& xQuantile084
)
{
  // compute median, -1 sigma and +1 sigma limits on reconstructed mass
  if ( histogram->Integral() > 0. ) {
    Double_t q[3];
    Double_t probSum[3];
    probSum[0] = 0.16;
    probSum[1] = 0.50;
    probSum[2] = 0.84;
    (const_cast<TH1*>(histogram))->GetQuantiles(3, q, probSum);
    xQuantile016 = q[0];
    xQuantile050 = q[1];
    xQuantile084 = q[2];
  } else {
    xQuantile016 = 0.;
    xQuantile050 = 0.;
    xQuantile084 = 0.;
  }

  xMean = histogram->GetMean();

  TH1* histogram_density = HistogramTools::compHistogramDensity(histogram);
  if ( histogram_density->Integral() > 0. )
  {
    int binMaximum = histogram_density->GetMaximumBin();
    xMaximum = histogram_density->GetBinCenter(binMaximum);
    double yMaximum = histogram_density->GetBinContent(binMaximum);
    if ( binMaximum > 1 && binMaximum < histogram_density->GetNbinsX() )
    {
      int binLeft       = binMaximum - 1;
      double xLeft      = histogram_density->GetBinCenter(binLeft);
      double yLeft      = histogram_density->GetBinContent(binLeft);

      int binRight      = binMaximum + 1;
      double xRight     = histogram_density->GetBinCenter(binRight);
      double yRight     = histogram_density->GetBinContent(binRight);

      double xMinus     = xLeft - xMaximum;
      double yMinus     = yLeft - yMaximum;
      double xPlus      = xRight - xMaximum;
      double yPlus      = yRight - yMaximum;

      xMaximum_interpol = xMaximum + 0.5*(yPlus*square(xMinus) - yMinus*square(xPlus))/(yPlus*xMinus - yMinus*xPlus);
    }
    else
    {
      xMaximum_interpol = xMaximum;
    }
  }
  else
  {
    xMaximum = 0.;
    xMaximum_interpol = 0.;
  }
  delete histogram_density;
}

double
HistogramTools::extractValue(TH1 const* histogram)
{
  double maximum, maximum_interpol, mean, quantile016, quantile050, quantile084;
  HistogramTools::extractHistogramProperties(histogram, maximum, maximum_interpol, mean, quantile016, quantile050, quantile084);
  double value = maximum;
  return value;
}

double
HistogramTools::extractUncertainty(TH1 const* histogram)
{
  double maximum, maximum_interpol, mean, quantile016, quantile050, quantile084;
  HistogramTools::extractHistogramProperties(histogram, maximum, maximum_interpol, mean, quantile016, quantile050, quantile084);
  double uncertainty = TMath::Sqrt(0.5*(TMath::Power(quantile084 - maximum, 2.) + TMath::Power(maximum - quantile016, 2.)));
  return uncertainty;
}

double
HistogramTools::extractLmax(TH1 const* histogram)
{
  TH1* histogram_density = HistogramTools::compHistogramDensity(histogram);
  double Lmax = histogram_density->GetMaximum();
  delete histogram_density;
  return Lmax;
}

TH1*
HistogramTools::makeHistogram_linBinWidth(const std::string& histogramName, int numBins, double xMin, double xMax)
{
  assert(numBins >= 1);
  TH1* histogram = new TH1D(histogramName.data(), histogramName.data(), numBins, xMin, xMax);
  return histogram;
}

TH1*
HistogramTools::makeHistogram_logBinWidth(const std::string& histogramName, double xMin, double xMax, double logBinWidth)
{
  if ( xMin <= 0. ) xMin = 0.1;
  int numBins = 1 + TMath::Log(xMax/xMin)/TMath::Log(logBinWidth);
  assert(numBins >= 1);
  TArrayF binning(numBins + 1);
  binning[0] = 0.;
  double x = xMin;
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin )
  {
    binning[idxBin] = x;
    x *= logBinWidth;
  }
  TH1* histogram = new TH1D(histogramName.data(), histogramName.data(), numBins, binning.GetArray());
  return histogram;
}

int SVfitQuantity::nInstances = 0;

SVfitQuantity::SVfitQuantity(const std::string& label) 
  : label_(label)
  , histogram_(nullptr)
  , uniqueName_("_SVfitQuantity_" + std::to_string(++SVfitQuantity::nInstances))
{}

SVfitQuantity::~SVfitQuantity()
{
  if ( histogram_ != nullptr ) delete histogram_;
  histogram_ = nullptr;
}

const TH1*
SVfitQuantity::getHistogram() const 
{ 
  return histogram_;
}

void
SVfitQuantity::writeHistogram() const
{
  if ( histogram_ != nullptr )
  {
    std::string histogramName = histogram_->GetName();
    boost::replace_all(histogramName, uniqueName_, "");
    histogram_->Write(histogramName.c_str(), TObject::kWriteDelete);
  }
}

void
SVfitQuantity::fillHistogram(double value, double weight)
{
  histogram_->Fill(value, weight);
}

double
SVfitQuantity::extractValue() const
{
  return HistogramTools::extractValue(histogram_);
}

double
SVfitQuantity::extractUncertainty() const
{
  return HistogramTools::extractUncertainty(histogram_);
}

double
SVfitQuantity::extractLmax() const
{
  return HistogramTools::extractLmax(histogram_);
}

bool
SVfitQuantity::isValidSolution() const
{
  return (extractLmax() > 0.);
}

HistogramAdapter::HistogramAdapter(const std::string& label) 
  : label_(label)
{}

HistogramAdapter::~HistogramAdapter()
{
  for ( SVfitQuantity* quantity : quantities_ )
  {
    delete quantity;
  }
}

void
HistogramAdapter::writeHistograms(const std::string& likelihoodFileName) const
{
  TFile* likelihoodFile = new TFile(likelihoodFileName.data(), "RECREATE");
  likelihoodFile->cd();

  for ( SVfitQuantity* quantity : quantities_ )
  {
    quantity->writeHistogram();
  }

  likelihoodFile->Write();
  likelihoodFile->Close();
  delete likelihoodFile;
}

double
HistogramAdapter::extractValue(const SVfitQuantity* quantity) const
{
  return quantity->extractValue();
}

double
HistogramAdapter::extractUncertainty(const SVfitQuantity* quantity) const
{
  return quantity->extractUncertainty();
}

double
HistogramAdapter::extractLmax(const SVfitQuantity* quantity) const
{
  return quantity->extractLmax();
}

bool
HistogramAdapter::isValidSolution() const
{
  return std::accumulate(quantities_.begin(), quantities_.end(), true,
                         [](bool result, SVfitQuantity* quantity) { return result && quantity->isValidSolution(); });
}

//-------------------------------------------------------------------------------------------------
// auxiliary classes to reconstruct pT, eta, phi of single tau leptons
SVfitQuantityTau::SVfitQuantityTau(const std::string& label)
  : SVfitQuantity(label)
{}

void
SVfitQuantityTau::bookHistogram(const MeasuredTauLepton& measuredTauLepton)
{
  if ( histogram_ != nullptr ) delete histogram_;
  histogram_ = createHistogram(measuredTauLepton);
  histogram_->SetName(std::string(histogram_->GetName() + uniqueName_).c_str());
}

SVfitQuantityTauPt::SVfitQuantityTauPt(const std::string& label)
  : SVfitQuantityTau(label)
{}

TH1* SVfitQuantityTauPt::createHistogram(const MeasuredTauLepton& measuredTauLepton) const
{
  return HistogramTools::makeHistogram_logBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramPt", 1., 1.e+3, 1.025);
}

SVfitQuantityTauEta::SVfitQuantityTauEta(const std::string& label)
  : SVfitQuantityTau(label)
{}

TH1*
SVfitQuantityTauEta::createHistogram(const MeasuredTauLepton& measuredTauLepton) const
{
  return HistogramTools::makeHistogram_linBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramEta", 198, -9.9, +9.9);
}

SVfitQuantityTauPhi::SVfitQuantityTauPhi(const std::string& label)
  : SVfitQuantityTau(label)
{}

TH1*
SVfitQuantityTauPhi::createHistogram(const MeasuredTauLepton& measuredTauLepton) const
{
  return HistogramTools::makeHistogram_linBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramEta", 180, -TMath::Pi(), +TMath::Pi());
}

HistogramAdapterTau::HistogramAdapterTau(const std::string& label)
  : HistogramAdapter(label)
  , quantity_pt_(nullptr)
  , quantity_eta_(nullptr)
  , quantity_phi_(nullptr)
{
  quantity_pt_ = new SVfitQuantityTauPt(label_);
  quantities_.push_back(quantity_pt_);
  quantity_eta_ = new SVfitQuantityTauEta(label_);
  quantities_.push_back(quantity_eta_);
  quantity_phi_ = new SVfitQuantityTauPhi(label_);
  quantities_.push_back(quantity_phi_);
}

HistogramAdapter*
HistogramAdapterTau::clone() const
{
  HistogramAdapterTau* retVal = new HistogramAdapterTau(this->label_);
  retVal->setMeasurement(this->measuredTauLepton_);
  retVal->bookHistograms(this->measuredTauLepton_);
  return retVal;
}

void
HistogramAdapterTau::setMeasurement(const MeasuredTauLepton& measuredTauLepton)
{
  measuredTauLepton_ = measuredTauLepton;
}
 
void
HistogramAdapterTau::setFittedTauLepton(const FittedTauLepton& fittedTauLepton)
{
  fittedTauP4_ = fittedTauLepton.tauP4();
}

void
HistogramAdapterTau::bookHistograms(const MeasuredTauLepton& measuredTauLepton)
{
  quantity_pt_->bookHistogram(measuredTauLepton);
  quantity_eta_->bookHistogram(measuredTauLepton);
  quantity_phi_->bookHistogram(measuredTauLepton);
}

void
HistogramAdapterTau::fillHistograms(double weight) const
{
  quantity_pt_->fillHistogram(fittedTauP4_.pt(), weight);
  quantity_eta_->fillHistogram(fittedTauP4_.eta(), weight);
  quantity_phi_->fillHistogram(fittedTauP4_.phi(), weight);
}

double
HistogramAdapterTau::getPt() const
{
  return extractValue(quantity_pt_);
}

double
HistogramAdapterTau::getPtErr() const
{
  return extractUncertainty(quantity_pt_);
}

double
HistogramAdapterTau::getPtLmax() const
{
  return extractLmax(quantity_pt_);
}

double
HistogramAdapterTau::getEta() const
{
  return extractValue(quantity_eta_);
}

double
HistogramAdapterTau::getEtaErr() const
{
  return extractUncertainty(quantity_eta_);
}

double
HistogramAdapterTau::getEtaLmax() const
{
  return extractLmax(quantity_eta_);
}

double
HistogramAdapterTau::getPhi() const
{
  return extractValue(quantity_phi_);
}

double
HistogramAdapterTau::getPhiErr() const
{
  return extractUncertainty(quantity_phi_);
}

double
HistogramAdapterTau::getPhiLmax() const
{
  return extractLmax(quantity_phi_);
}

LorentzVector
HistogramAdapterTau::getP4() const
{
  TLorentzVector p4;
  p4.SetPtEtaPhiM(this->getPt(), this->getEta(), this->getPhi(), tauLeptonMass);
  return LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E());
}

double
HistogramAdapterTau::DoEval(const double* x) const
{
  fillHistograms();
  return 0.;
}
//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
// auxiliary classes to reconstruct pT, eta, phi, mass, and transverse mass of tau lepton pairs
SVfitQuantityDiTau::SVfitQuantityDiTau(const std::string& label)
  : SVfitQuantity(label)
{}

void
SVfitQuantityDiTau::bookHistogram(const MeasuredEvent& measuredEvent)
{
  if ( histogram_ != nullptr ) delete histogram_;
  histogram_ = createHistogram(measuredEvent);
  histogram_->SetName(std::string(histogram_->GetName() + uniqueName_).c_str());
}

SVfitQuantityDiTauPt::SVfitQuantityDiTauPt(const std::string& label)
  : SVfitQuantityDiTau(label)
{}

TH1*
SVfitQuantityDiTauPt::createHistogram(const MeasuredEvent& measuredEvent) const
{
  return HistogramTools::makeHistogram_logBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramPt", 1., 1.e+3, 1.025);
}

SVfitQuantityDiTauEta::SVfitQuantityDiTauEta(const std::string& label)
  : SVfitQuantityDiTau(label)
{}

TH1*
SVfitQuantityDiTauEta::createHistogram(const MeasuredEvent& measuredEvent) const
{
  return HistogramTools::makeHistogram_linBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramEta", 198, -9.9, +9.9);
}

SVfitQuantityDiTauPhi::SVfitQuantityDiTauPhi(const std::string& label)
  : SVfitQuantityDiTau(label)
{}

TH1* SVfitQuantityDiTauPhi::createHistogram(const MeasuredEvent& measuredEvent) const
{
  return HistogramTools::makeHistogram_linBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramPhi", 180, -TMath::Pi(), +TMath::Pi());
}

SVfitQuantityDiTauMass::SVfitQuantityDiTauMass(const std::string& label)
  : SVfitQuantityDiTau(label)
{}

TH1*
SVfitQuantityDiTauMass::createHistogram(const MeasuredEvent& measuredEvent) const
{
  const std::vector<MeasuredTauLepton>& measuredTauLeptons = measuredEvent.measuredTauLeptons();
  assert(measuredTauLeptons.size() == 2);
  const LorentzVector& vis1P4 = measuredTauLeptons[0].p4();
  const LorentzVector& vis2P4 = measuredTauLeptons[1].p4();
  double visMass = (vis1P4 + vis2P4).mass();
  double minMass = visMass/1.0125;
  double maxMass = TMath::Max(1.e+4, 1.e+1*minMass);
  return HistogramTools::makeHistogram_logBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramMass", minMass, maxMass, 1.025);
}

SVfitQuantityDiTauTransverseMass::SVfitQuantityDiTauTransverseMass(const std::string& label)
  : SVfitQuantityDiTau(label)
{}

TH1*
SVfitQuantityDiTauTransverseMass::createHistogram(const MeasuredEvent& measuredEvent) const
{
  const std::vector<MeasuredTauLepton>& measuredTauLeptons = measuredEvent.measuredTauLeptons();
  assert(measuredTauLeptons.size() == 2);
  const LorentzVector& vis1P4 = measuredTauLeptons[0].p4();
  const LorentzVector& vis2P4 = measuredTauLeptons[1].p4();
  LorentzVector visDiTauP4 = vis1P4 + vis2P4;
  double visTransverseMass2 = square(vis1P4.Et() + vis2P4.Et()) - (square(visDiTauP4.px()) + square(visDiTauP4.py()));
  double visTransverseMass = TMath::Sqrt(TMath::Max(1., visTransverseMass2));
  double minTransverseMass = visTransverseMass/1.0125;
  double maxTransverseMass = TMath::Max(1.e+4, 1.e+1*minTransverseMass);
  return HistogramTools::makeHistogram_logBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramTransverseMass", minTransverseMass, maxTransverseMass, 1.025);
}
    
HistogramAdapterDiTau::HistogramAdapterDiTau(const std::string& label)
  : HistogramAdapter(label)
  , quantity_pt_(nullptr)
  , quantity_eta_(nullptr)
  , quantity_phi_(nullptr)
  , quantity_mass_(nullptr)
  , quantity_transverseMass_(nullptr)
  , adapter_tau1_(nullptr)
  , adapter_tau2_(nullptr)
{
  quantity_pt_ = new SVfitQuantityDiTauPt(label_);
  quantities_.push_back(quantity_pt_);
  quantity_eta_ = new SVfitQuantityDiTauEta(label_);
  quantities_.push_back(quantity_eta_);
  quantity_phi_ = new SVfitQuantityDiTauPhi(label_);
  quantities_.push_back(quantity_phi_);
  quantity_mass_ = new SVfitQuantityDiTauMass(label_);
  quantities_.push_back(quantity_mass_);
  quantity_transverseMass_ = new SVfitQuantityDiTauTransverseMass(label_);
  quantities_.push_back(quantity_transverseMass_);

  adapter_tau1_ = new HistogramAdapterTau(label_ + "_tau1");
  adapter_tau2_ = new HistogramAdapterTau(label_ + "_tau2");
}

HistogramAdapter*
HistogramAdapterDiTau::clone() const
{
  HistogramAdapterDiTau* retVal = new HistogramAdapterDiTau(this->label_);
  retVal->setMeasurement(this->measuredEvent_);
  retVal->bookHistograms(this->measuredEvent_);
  return retVal;
}

HistogramAdapterDiTau::~HistogramAdapterDiTau()
{
  delete adapter_tau1_;
  delete adapter_tau2_;
}

void
HistogramAdapterDiTau::setMeasurement(const MeasuredEvent& measuredEvent)
{
  measuredEvent_ = measuredEvent;
  const std::vector<MeasuredTauLepton>& measuredTauLeptons = measuredEvent.measuredTauLeptons();
  assert(measuredTauLeptons.size() == 2);
  adapter_tau1_->setMeasurement(measuredTauLeptons[0]);
  adapter_tau2_->setMeasurement(measuredTauLeptons[1]);
}

void
HistogramAdapterDiTau::setFittedTauLeptons(const FittedTauLepton& fittedTauLepton1,  const FittedTauLepton& fittedTauLepton2)
{
  fittedTau1P4_ = fittedTauLepton1.tauP4();
  fittedTau2P4_ = fittedTauLepton2.tauP4();
  fittedDiTauP4_ = fittedTau1P4_ + fittedTau2P4_;
  adapter_tau1_->setFittedTauLepton(fittedTauLepton1);
  adapter_tau2_->setFittedTauLepton(fittedTauLepton2);
}

void
HistogramAdapterDiTau::bookHistograms(const MeasuredEvent& measuredEvent)
{
  quantity_pt_->bookHistogram(measuredEvent);
  quantity_eta_->bookHistogram(measuredEvent);
  quantity_phi_->bookHistogram(measuredEvent);
  quantity_mass_->bookHistogram(measuredEvent);
  quantity_transverseMass_->bookHistogram(measuredEvent);
  measuredEvent_ = measuredEvent;
  const std::vector<MeasuredTauLepton>& measuredTauLeptons = measuredEvent.measuredTauLeptons();
  assert(measuredTauLeptons.size() == 2);
  adapter_tau1_->bookHistograms(measuredTauLeptons[0]);
  adapter_tau2_->bookHistograms(measuredTauLeptons[1]);
}

void
HistogramAdapterDiTau::fillHistograms(double weight) const
{
  quantity_pt_->fillHistogram(fittedDiTauP4_.pt(), weight);
  quantity_eta_->fillHistogram(fittedDiTauP4_.eta(), weight);
  quantity_phi_->fillHistogram(fittedDiTauP4_.phi(), weight);
  quantity_mass_->fillHistogram(fittedDiTauP4_.mass(), weight);
  double transverseMass2 = square(fittedTau1P4_.Et() + fittedTau2P4_.Et()) - (square(fittedDiTauP4_.px()) + square(fittedDiTauP4_.py()));
  quantity_transverseMass_->fillHistogram(TMath::Sqrt(TMath::Max(1., transverseMass2)), weight);
  adapter_tau1_->fillHistograms(weight);
  adapter_tau2_->fillHistograms(weight);
}

HistogramAdapterTau*
HistogramAdapterDiTau::tau1() const 
{ 
  return adapter_tau1_; 
}
 
HistogramAdapterTau*
HistogramAdapterDiTau::tau2() const 
{
  return adapter_tau2_; 
}

double
HistogramAdapterDiTau::getPt() const
{
  return extractValue(quantity_pt_);
}

double
HistogramAdapterDiTau::getPtErr() const
{
  return extractUncertainty(quantity_pt_);
}

double
HistogramAdapterDiTau::getPtLmax() const
{
  return extractLmax(quantity_pt_);
}

double
HistogramAdapterDiTau::getEta() const
{
  return extractValue(quantity_eta_);
}

double
HistogramAdapterDiTau::getEtaErr() const
{
  return extractUncertainty(quantity_eta_);
}

double
HistogramAdapterDiTau::getEtaLmax() const
{
  return extractLmax(quantity_eta_);
}

double
HistogramAdapterDiTau::getPhi() const
{
  return extractValue(quantity_phi_);
}

double
HistogramAdapterDiTau::getPhiErr() const
{
  return extractUncertainty(quantity_phi_);
}

double
HistogramAdapterDiTau::getPhiLmax() const
{
  return extractLmax(quantity_phi_);
}

double
HistogramAdapterDiTau::getMass() const
{
  return extractValue(quantity_mass_);
}

double
HistogramAdapterDiTau::getMassErr() const
{
  return extractUncertainty(quantity_mass_);
}

double
HistogramAdapterDiTau::getMassLmax() const
{
  return extractLmax(quantity_mass_);
}

double
HistogramAdapterDiTau::getTransverseMass() const
{
  return extractValue(quantity_transverseMass_);
}

double
HistogramAdapterDiTau::getTransverseMassErr() const
{
  return extractUncertainty(quantity_transverseMass_);
}

double
HistogramAdapterDiTau::getTransverseMassLmax() const
{
  return extractLmax(quantity_transverseMass_);
}

LorentzVector
HistogramAdapterDiTau::getP4() const
{
  TLorentzVector p4;
  p4.SetPtEtaPhiM(this->getPt(), this->getEta(), this->getPhi(), this->getMass());
  return LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E());
}

double
HistogramAdapterDiTau::DoEval(const double* x) const
{
  fillHistograms();
  return 0.;
}
//-------------------------------------------------------------------------------------------------

}
