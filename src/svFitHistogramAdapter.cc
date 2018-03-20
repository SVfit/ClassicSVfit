#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include <TMath.h>
#include <TFile.h>
#include <TObject.h>
#include <TLorentzVector.h>

#include <numeric>

#include <boost/algorithm/string/replace.hpp>

using namespace classic_svFit;

TH1* HistogramTools::compHistogramDensity(TH1 const* histogram)
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
  if ( histogram_density->Integral() > 0. ) {
    int binMaximum = histogram_density->GetMaximumBin();
    xMaximum = histogram_density->GetBinCenter(binMaximum);
    double yMaximum = histogram_density->GetBinContent(binMaximum);
    if ( binMaximum > 1 && binMaximum < histogram_density->GetNbinsX() ) {
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
    } else {
      xMaximum_interpol = xMaximum;
    }
  } else {
    xMaximum = 0.;
    xMaximum_interpol = 0.;
  }
  delete histogram_density;
}

double HistogramTools::extractValue(TH1 const* histogram)
{
  double maximum, maximum_interpol, mean, quantile016, quantile050, quantile084;
  HistogramTools::extractHistogramProperties(histogram, maximum, maximum_interpol, mean, quantile016, quantile050, quantile084);
  double value = maximum;
  return value;
}

double HistogramTools::extractUncertainty(TH1 const* histogram)
{
  double maximum, maximum_interpol, mean, quantile016, quantile050, quantile084;
  HistogramTools::extractHistogramProperties(histogram, maximum, maximum_interpol, mean, quantile016, quantile050, quantile084);
  double uncertainty = TMath::Sqrt(0.5*(TMath::Power(quantile084 - maximum, 2.) + TMath::Power(maximum - quantile016, 2.)));
  return uncertainty;
}

double HistogramTools::extractLmax(TH1 const* histogram)
{
  TH1* histogram_density = HistogramTools::compHistogramDensity(histogram);
  double Lmax = histogram_density->GetMaximum();
  delete histogram_density;
  return Lmax;
}

TH1* HistogramTools::makeHistogram_linBinWidth(const std::string& histogramName, int numBins, double xMin, double xMax)
{
  TH1* histogram = new TH1D(histogramName.data(), histogramName.data(), numBins, xMin, xMax);
  return histogram;
}

TH1* HistogramTools::makeHistogram_logBinWidth(const std::string& histogramName, double xMin, double xMax, double logBinWidth)
{
  if ( xMin <= 0. ) xMin = 0.1;
  int numBins = 1 + TMath::Log(xMax/xMin)/TMath::Log(logBinWidth);
  TArrayF binning(numBins + 1);
  binning[0] = 0.;
  double x = xMin;
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) {
    binning[idxBin] = x;
    x *= logBinWidth;
  }
  TH1* histogram = new TH1D(histogramName.data(), histogramName.data(), numBins, binning.GetArray());
  return histogram;
}

int SVfitQuantity::nInstances = 0;

SVfitQuantity::SVfitQuantity(const std::string& label) 
  : label_(label)
  , uniqueName_("_SVfitQuantity_" + std::to_string(++SVfitQuantity::nInstances))
{}

SVfitQuantity::~SVfitQuantity()
{
  if ( histogram_ != nullptr ) delete histogram_;
}

const TH1* SVfitQuantity::getHistogram() const 
{ 
  return histogram_;
}

void SVfitQuantity::writeHistogram() const
{
  if ( histogram_ != nullptr ) {
    std::string histogramName = histogram_->GetName();
    boost::replace_all(histogramName, uniqueName_, "");
    histogram_->Write(histogramName.c_str(), TObject::kWriteDelete);
  }
}

void SVfitQuantity::fillHistogram(double value)
{
  histogram_->Fill(value);
}

double SVfitQuantity::extractValue() const
{
  return HistogramTools::extractValue(histogram_);
}

double SVfitQuantity::extractUncertainty() const
{
  return HistogramTools::extractUncertainty(histogram_);
}

double SVfitQuantity::extractLmax() const
{
  return HistogramTools::extractLmax(histogram_);
}

bool SVfitQuantity::isValidSolution() const
{
  return (extractLmax() > 0.);
}

HistogramAdapter::HistogramAdapter(const std::string& label) 
  : label_(label)
{}

HistogramAdapter::~HistogramAdapter()
{
  for ( std::vector<SVfitQuantity*>::iterator quantity = quantities_.begin(); 
	quantity != quantities_.end(); ++quantity ) {
    delete (*quantity);
  }
}

void HistogramAdapter::writeHistograms(const std::string& likelihoodFileName) const
{
  TFile* likelihoodFile = new TFile(likelihoodFileName.data(), "RECREATE");
  likelihoodFile->cd();

  for (std::vector<SVfitQuantity*>::iterator quantity = quantities_.begin(); 
       quantity != quantities_.end(); ++quantity ) {
    (*quantity)->writeHistogram();
  }

  likelihoodFile->Write();
  likelihoodFile->Close();
  delete likelihoodFile;
}

double HistogramAdapter::extractValue(const SVfitQuantity* quantity) const
{
  return quantity->extractValue();
}

double HistogramAdapter::extractUncertainty(const SVfitQuantity* quantity) const
{
  return quantity->extractUncertainty();
}

double HistogramAdapter::extractLmax(const SVfitQuantity* quantity) const
{
  return quantity->extractLmax();
}

bool HistogramAdapter::isValidSolution() const
{
  return std::accumulate(quantities_.begin(), quantities_.end(), true,
                         [](bool result, SVfitQuantity* quantity) { return result && quantity->isValidSolution(); });
}

//-------------------------------------------------------------------------------------------------
// auxiliary classes to reconstruct pT, eta, phi of single tau leptons
SVfitQuantityTau::SVfitQuantityTau(const std::string& label)
  : SVfitQuantity(label)
{}

void SVfitQuantityTau::bookHistogram(const LorentzVector& visP4)
{
  if ( histogram_ != nullptr ) delete histogram_;
  histogram_ = createHistogram(visP4);
  histogram_->SetName(std::string(histogram_->GetName() + uniqueName_).c_str());
}

SVfitQuantityTauPt::SVfitQuantityTauPt(const std::string& label)
  : SVfitQuantityTau(label)
{}

TH1* SVfitQuantityTauPt::createHistogram(const LorentzVector& visP4) const
{
  return HistogramTools::makeHistogram_logBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramPt", 1., 1.e+3, 1.025);
}

SVfitQuantityTauEta::SVfitQuantityTauEta(const std::string& label)
  : SVfitQuantityTau(label)
{}

TH1* SVfitQuantityTauEta::createHistogram(const LorentzVector& visP4) const
{
  return HistogramTools::makeHistogram_linBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramEta", 198, -9.9, +9.9);
}

SVfitQuantityTauPhi::SVfitQuantityTauPhi(const std::string& label)
  : SVfitQuantityTau(label)
{}

TH1* SVfitQuantityTauPhi::createHistogram(const LorentzVector& visP4) const
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

void HistogramAdapterTau::setMeasurement(const LorentzVector& visP4)
{
  visP4_ = visP4;
}
 
void HistogramAdapterTau::setTauP4(const LorentzVector& tauP4)
{
  tauP4_ = tauP4;
}

void HistogramAdapterTau::bookHistograms(const LorentzVector& visP4)
{
  quantity_pt_->bookHistogram(visP4);
  quantity_eta_->bookHistogram(visP4);
  quantity_phi_->bookHistogram(visP4);
}

void HistogramAdapterTau::fillHistograms(const LorentzVector& tauP4, const LorentzVector& visP4) const
{
  quantity_pt_->fillHistogram(tauP4.pt());
  quantity_eta_->fillHistogram(tauP4.eta());
  quantity_phi_->fillHistogram(tauP4.phi());
}

double HistogramAdapterTau::getPt() const
{
  return extractValue(quantity_pt_);
}

double HistogramAdapterTau::getPtErr() const
{
  return extractUncertainty(quantity_pt_);
}

double HistogramAdapterTau::getPtLmax() const
{
  return extractLmax(quantity_pt_);
}

double HistogramAdapterTau::getEta() const
{
  return extractValue(quantity_eta_);
}

double HistogramAdapterTau::getEtaErr() const
{
  return extractUncertainty(quantity_eta_);
}

double HistogramAdapterTau::getEtaLmax() const
{
  return extractLmax(quantity_eta_);
}

double HistogramAdapterTau::getPhi() const
{
  return extractValue(quantity_phi_);
}

double HistogramAdapterTau::getPhiErr() const
{
  return extractUncertainty(quantity_phi_);
}

double HistogramAdapterTau::getPhiLmax() const
{
  return extractLmax(quantity_phi_);
}

classic_svFit::LorentzVector HistogramAdapterTau::getP4() const
{
  TLorentzVector p4;
  p4.SetPtEtaPhiM(this->getPt(), this->getEta(), this->getPhi(), tauLeptonMass);
  return classic_svFit::LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E());
}

double HistogramAdapterTau::DoEval(const double* x) const
{
  fillHistograms(tauP4_, visP4_);
  return 0.;
}
//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
// auxiliary classes to reconstruct pT, eta, phi, mass, and transverse mass of tau lepton pairs
SVfitQuantityDiTau::SVfitQuantityDiTau(const std::string& label)
  : SVfitQuantity(label)
{}

void SVfitQuantityDiTau::bookHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met)
{
  if ( histogram_ != nullptr ) delete histogram_;
  histogram_ = createHistogram(vis1P4, vis2P4, met);
  histogram_->SetName(std::string(histogram_->GetName() + uniqueName_).c_str());
}

SVfitQuantityDiTauPt::SVfitQuantityDiTauPt(const std::string& label)
  : SVfitQuantityDiTau(label)
{}

TH1* SVfitQuantityDiTauPt::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  return HistogramTools::makeHistogram_logBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramPt", 1., 1.e+3, 1.025);
}

SVfitQuantityDiTauEta::SVfitQuantityDiTauEta(const std::string& label)
  : SVfitQuantityDiTau(label)
{}

TH1* SVfitQuantityDiTauEta::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  return HistogramTools::makeHistogram_linBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramEta", 198, -9.9, +9.9);
}

SVfitQuantityDiTauPhi::SVfitQuantityDiTauPhi(const std::string& label)
  : SVfitQuantityDiTau(label)
{}

TH1* SVfitQuantityDiTauPhi::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  return HistogramTools::makeHistogram_linBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramPhi", 180, -TMath::Pi(), +TMath::Pi());
}

SVfitQuantityDiTauMass::SVfitQuantityDiTauMass(const std::string& label)
  : SVfitQuantityDiTau(label)
{}

TH1* SVfitQuantityDiTauMass::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  double visMass = (vis1P4 + vis2P4).mass();
  double minMass = visMass/1.0125;
  double maxMass = TMath::Max(1.e+4, 1.e+1*minMass);
  return HistogramTools::makeHistogram_logBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramMass", minMass, maxMass, 1.025);
}

SVfitQuantityDiTauTransverseMass::SVfitQuantityDiTauTransverseMass(const std::string& label)
  : SVfitQuantityDiTau(label)
{}

TH1* SVfitQuantityDiTauTransverseMass::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  classic_svFit::LorentzVector measuredDiTauSystem = vis1P4 + vis2P4;
  double visTransverseMass2 = square(vis1P4.Et() + vis2P4.Et()) - (square(measuredDiTauSystem.px()) + square(measuredDiTauSystem.py()));
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

HistogramAdapterDiTau::~HistogramAdapterDiTau()
{
  delete adapter_tau1_;
  delete adapter_tau2_;
}

void HistogramAdapterDiTau::setMeasurement(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met)
{
  vis1P4_ = vis1P4;
  vis2P4_ = vis2P4;
  met_ = met;
  adapter_tau1_->setMeasurement(vis1P4);
  adapter_tau2_->setMeasurement(vis2P4);
}

void HistogramAdapterDiTau::setTau1And2P4(const LorentzVector& tau1P4, const LorentzVector& tau2P4) 
{
  tau1P4_ = tau1P4;
  tau2P4_ = tau2P4;
  ditauP4_ = tau1P4_ + tau2P4_;
  adapter_tau1_->setTauP4(tau1P4);
  adapter_tau2_->setTauP4(tau2P4);
}

void HistogramAdapterDiTau::bookHistograms(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met)
{
  quantity_pt_->bookHistogram(vis1P4, vis2P4, met);
  quantity_eta_->bookHistogram(vis1P4, vis2P4, met);
  quantity_phi_->bookHistogram(vis1P4, vis2P4, met);
  quantity_mass_->bookHistogram(vis1P4, vis2P4, met);
  quantity_transverseMass_->bookHistogram(vis1P4, vis2P4, met);
  adapter_tau1_->bookHistograms(vis1P4);
  adapter_tau2_->bookHistograms(vis2P4);
}

void HistogramAdapterDiTau::fillHistograms(const LorentzVector& tau1P4, const LorentzVector& tau2P4, const LorentzVector& ditauP4,
					   const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  quantity_pt_->fillHistogram(ditauP4.pt());
  quantity_eta_->fillHistogram(ditauP4.eta());
  quantity_phi_->fillHistogram(ditauP4.phi());
  quantity_mass_->fillHistogram(ditauP4.mass());
  double transverseMass2 = square(tau1P4.Et() + tau2P4.Et()) - (square(ditauP4.px()) + square(ditauP4.py()));
  quantity_transverseMass_->fillHistogram(TMath::Sqrt(TMath::Max(1., transverseMass2)));
  adapter_tau1_->fillHistograms(tau1P4, vis1P4);
  adapter_tau2_->fillHistograms(tau2P4, vis2P4);
}

HistogramAdapterTau* HistogramAdapterDiTau::tau1() const 
{ 
  return adapter_tau1_; 
}
 
HistogramAdapterTau* HistogramAdapterDiTau::tau2() const 
{
  return adapter_tau2_; 
}

double HistogramAdapterDiTau::getPt() const
{
  return extractValue(quantity_pt_);
}

double HistogramAdapterDiTau::getPtErr() const
{
  return extractUncertainty(quantity_pt_);
}

double HistogramAdapterDiTau::getPtLmax() const
{
  return extractLmax(quantity_pt_);
}

double HistogramAdapterDiTau::getEta() const
{
  return extractValue(quantity_eta_);
}

double HistogramAdapterDiTau::getEtaErr() const
{
  return extractUncertainty(quantity_eta_);
}

double HistogramAdapterDiTau::getEtaLmax() const
{
  return extractLmax(quantity_eta_);
}

double HistogramAdapterDiTau::getPhi() const
{
  return extractValue(quantity_phi_);
}

double HistogramAdapterDiTau::getPhiErr() const
{
  return extractUncertainty(quantity_phi_);
}

double HistogramAdapterDiTau::getPhiLmax() const
{
  return extractLmax(quantity_phi_);
}

double HistogramAdapterDiTau::getMass() const
{
  return extractValue(quantity_mass_);
}

double HistogramAdapterDiTau::getMassErr() const
{
  return extractUncertainty(quantity_mass_);
}

double HistogramAdapterDiTau::getMassLmax() const
{
  return extractLmax(quantity_mass_);
}

double HistogramAdapterDiTau::getTransverseMass() const
{
  return extractValue(quantity_transverseMass_);
}

double HistogramAdapterDiTau::getTransverseMassErr() const
{
  return extractUncertainty(quantity_transverseMass_);
}

double HistogramAdapterDiTau::getTransverseMassLmax() const
{
  return extractLmax(quantity_transverseMass_);
}

classic_svFit::LorentzVector HistogramAdapterDiTau::getP4() const
{
  TLorentzVector p4;
  p4.SetPtEtaPhiM(this->getPt(), this->getEta(), this->getPhi(), this->getMass());
  return classic_svFit::LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E());
}

double HistogramAdapterDiTau::DoEval(const double* x) const
{
  fillHistograms(tau1P4_, tau2P4_, ditauP4_, vis1P4_, vis2P4_, met_);
  return 0.;
}
//-------------------------------------------------------------------------------------------------
