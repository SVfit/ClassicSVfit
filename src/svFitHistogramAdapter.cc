#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include <TMath.h>
#include <TFile.h>

using namespace classic_svFit;

HistogramAdapter::HistogramAdapter()
  : histogramPt_(0),
    histogramPt_density_(0),
    histogramEta_(0),
    histogramEta_density_(0),
    histogramPhi_(0),
    histogramPhi_density_(0),
    histogramMass_(0),
    histogramMass_density_(0),
    histogramTransverseMass_(0),
    histogramTransverseMass_density_(0)
{}

HistogramAdapter::~HistogramAdapter()
{
  delete histogramPt_;
  delete histogramPt_density_;
  delete histogramEta_;
  delete histogramEta_density_;
  delete histogramPhi_;
  delete histogramPhi_density_;
  delete histogramMass_;
  delete histogramMass_density_;
  delete histogramTransverseMass_;
  delete histogramTransverseMass_density_;
}

namespace
{
  TH1* makeHistogram(const std::string& histogramName, double xMin, double xMax, double logBinWidth)
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
}

void HistogramAdapter::bookHistograms(const LorentzVector& vis1P4, const LorentzVector& vis2P4)
{
  // CV: book histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
  LorentzVector visDiTauP4 = vis1P4 + vis2P4;
  delete histogramPt_;
  histogramPt_ = makeHistogram("ClassicSVfitIntegrand_histogramPt", 1., 1.e+3, 1.025);
  delete histogramPt_density_;
  histogramPt_density_ = (TH1*)histogramPt_->Clone(Form("%s_density", histogramPt_->GetName()));
  delete histogramEta_;
  histogramEta_ = new TH1D("ClassicSVfitIntegrand_histogramEta", "ClassicSVfitIntegrand_histogramEta", 198, -9.9, +9.9);
  delete histogramEta_density_;
  histogramEta_density_ = (TH1*)histogramEta_->Clone(Form("%s_density", histogramEta_->GetName()));
  delete histogramPhi_;
  histogramPhi_ = new TH1D("ClassicSVfitIntegrand_histogramPhi", "ClassicSVfitIntegrand_histogramPhi", 180, -TMath::Pi(), +TMath::Pi());
  delete histogramPhi_density_;
  histogramPhi_density_ = (TH1*)histogramPhi_->Clone(Form("%s_density", histogramPhi_->GetName()));
  double mVis_measured = visDiTauP4.mass();
  double minMass = mVis_measured/1.0125;
  double maxMass = TMath::Max(1.e+4, 1.e+1*minMass);
  delete histogramMass_;
  histogramMass_ = makeHistogram("ClassicSVfitIntegrand_histogramMass", minMass, maxMass, 1.025);
  delete histogramMass_density_;
  histogramMass_density_ = (TH1*)histogramMass_->Clone(Form("%s_density", histogramMass_->GetName()));
  double mTvis2_measured = square(vis1P4.Et() + vis2P4.Et()) - (square(visDiTauP4.px()) + square(visDiTauP4.py()));
  double mTvis_measured = TMath::Sqrt(TMath::Max(1., mTvis2_measured));
  double minTransverseMass = mTvis_measured/1.0125;
  double maxTransverseMass = TMath::Max(1.e+4, 1.e+1*minTransverseMass);
  delete histogramTransverseMass_;
  histogramTransverseMass_ = makeHistogram("ClassicSVfitIntegrand_histogramTransverseMass", minTransverseMass, maxTransverseMass, 1.025);
  delete histogramTransverseMass_density_;
  histogramTransverseMass_density_ = (TH1*)histogramTransverseMass_->Clone(Form("%s_density", histogramTransverseMass_->GetName()));
}

void HistogramAdapter::fillHistograms(const LorentzVector& tau1P4, const LorentzVector& tau2P4) const
{
  // CV: fill histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
  LorentzVector fittedDiTauP4 = tau1P4 + tau2P4;
  histogramPt_->Fill(fittedDiTauP4.pt());
  histogramEta_->Fill(fittedDiTauP4.eta());
  histogramPhi_->Fill(fittedDiTauP4.phi());
  histogramMass_->Fill(fittedDiTauP4.mass());
  double transverseMass2 = square(tau1P4.Et() + tau2P4.Et()) - (square(fittedDiTauP4.px()) + square(fittedDiTauP4.py()));
  double transverseMass = TMath::Sqrt(TMath::Max(1., transverseMass2));
  histogramTransverseMass_->Fill(transverseMass);
}

void HistogramAdapter::writeHistograms(const std::string& likelihoodFileName) const
{
  TFile* likelihoodFile = new TFile(likelihoodFileName.data(), "RECREATE");
  histogramPt_->Write();
  histogramPt_density_->Write();
  histogramEta_->Write();
  histogramEta_density_->Write();
  histogramPhi_->Write();
  histogramPhi_density_->Write();
  histogramMass_->Write();
  histogramMass_density_->Write();
  histogramTransverseMass_->Write();
  histogramTransverseMass_density_->Write();
  delete likelihoodFile;
}

namespace
{
  void compHistogramDensity(const TH1* histogram, TH1* histogram_density)
  {
    for ( int idxBin = 1; idxBin <= histogram->GetNbinsX(); ++idxBin ) {
      double binContent = histogram->GetBinContent(idxBin);
      double binError = histogram->GetBinError(idxBin);
      double binWidth = histogram->GetBinWidth(idxBin);
      assert(binWidth > 0.);
      histogram_density->SetBinContent(idxBin, binContent/binWidth);
      histogram_density->SetBinError(idxBin, binError/binWidth);
    }
  }

  void extractHistogramProperties(const TH1* histogram, const TH1* histogram_density,
                                  double& xMaximum, double& xMaximum_interpol,
                                  double& xMean,
                                  double& xQuantile016, double& xQuantile050, double& xQuantile084)
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
  }

  double extractValue(const TH1* histogram, TH1* histogram_density)
  {
    double maximum, maximum_interpol, mean, quantile016, quantile050, quantile084;
    compHistogramDensity(histogram, histogram_density);
    extractHistogramProperties(
      histogram, histogram_density,
      maximum, maximum_interpol, mean, quantile016, quantile050, quantile084);
    double value = maximum;
    return value;
  }

  double extractUncertainty(const TH1* histogram, TH1* histogram_density)
  {
    double maximum, maximum_interpol, mean, quantile016, quantile050, quantile084;
    compHistogramDensity(histogram, histogram_density);
    extractHistogramProperties(
      histogram, histogram_density,
      maximum, maximum_interpol, mean, quantile016, quantile050, quantile084);
    double uncertainty = TMath::Sqrt(0.5*(TMath::Power(quantile084 - maximum, 2.) + TMath::Power(maximum - quantile016, 2.)));
    return uncertainty;
  }

  double extractLmax(const TH1* histogram, TH1* histogram_density)
  {
    compHistogramDensity(histogram, histogram_density);
    double Lmax = histogram_density->GetMaximum();
    return Lmax;
  }
}

double HistogramAdapter::getPt() const
{
  return extractValue(histogramPt_, histogramPt_density_);
}

double HistogramAdapter::getPtErr() const
{
  return extractUncertainty(histogramPt_, histogramPt_density_);
}

double HistogramAdapter::getPtLmax() const
{
  return extractLmax(histogramPt_, histogramPt_density_);
}

double HistogramAdapter::getEta() const
{
  return extractValue(histogramEta_, histogramEta_density_);
}

double HistogramAdapter::getEtaErr() const
{
  return extractUncertainty(histogramEta_, histogramEta_density_);
}

double HistogramAdapter::getEtaLmax() const
{
  return extractLmax(histogramEta_, histogramEta_density_);
}

double HistogramAdapter::getPhi() const
{
  return extractValue(histogramPhi_, histogramPhi_density_);
}

double HistogramAdapter::getPhiErr() const
{
  return extractUncertainty(histogramPhi_, histogramPhi_density_);
}

double HistogramAdapter::getPhiLmax() const
{
  return extractLmax(histogramPhi_, histogramPhi_density_);
}

double HistogramAdapter::getMass() const
{
  return extractValue(histogramMass_, histogramMass_density_);
}

double HistogramAdapter::getMassErr() const
{
  return extractUncertainty(histogramMass_, histogramMass_density_);
}

double HistogramAdapter::getMassLmax() const
{
  return extractLmax(histogramMass_, histogramMass_density_);
}

double HistogramAdapter::getTransverseMass() const
{
  return extractValue(histogramTransverseMass_, histogramTransverseMass_density_);
}

double HistogramAdapter::getTransverseMassErr() const
{
  return extractUncertainty(histogramTransverseMass_, histogramTransverseMass_density_);
}

double HistogramAdapter::getTransverseMassLmax() const
{
  return extractLmax(histogramTransverseMass_, histogramTransverseMass_density_);
}
