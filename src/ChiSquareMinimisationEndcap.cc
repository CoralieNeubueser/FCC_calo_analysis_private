#include "ChiSquareMinimisationEndcap.h"

#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenVertexCollection.h"
#include "datamodel/PositionedCaloHitCollection.h"

#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

#include "TVector3.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

// STL
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <bitset>

std::vector<int> EEndcap;
std::vector<double> vecEemEndcap;
std::vector<double> vecEhadEndcap;
std::vector<double> vecEfirstEndcap;
std::vector<double> vecElastEndcap;

Double_t chiSquareFitEndcap(const Double_t *par) {
  Double_t fitvalue = 0.;
  Double_t E0 = 0.;
  
  //loop over all events
  for(int i=0; i<vecEemEndcap.size(); i++){
    double Eem = vecEemEndcap.at(i);
    double Ehad = vecEhadEndcap.at(i);
    double Efirst = vecEfirstEndcap.at(i);
    double Elast = vecElastEndcap.at(i);
    
    E0 = Eem*par[0] + Ehad*par[1] + par[2]*Efirst*par[0] + par[3]*pow(Eem*par[0],2);
    
    fitvalue += pow(E0 - EEndcap.at(i),2)/(EEndcap.at(i)*par[4]); 
  }
  return fitvalue;    
}

ChiSquareMinimisationEndcap::ChiSquareMinimisationEndcap(double aEnergy, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal):
m_energy(aEnergy), ecal_decoder(aBitfieldEcal), hcal_decoder(aBitfieldHcal) {
  Initialize_histos();
}

ChiSquareMinimisationEndcap::~ChiSquareMinimisationEndcap() {}

void ChiSquareMinimisationEndcap::Initialize_histos()
{
 
  h_parameter = new TH1F("h_parameter","", 5, 0,5);
  m_histograms.push_back(h_parameter);

  h_energyECAL = new TH1F("h_energyECAL","", 100, 0, m_energy);
  m_histograms.push_back(h_energyECAL);

  h_energyHCAL = new TH1F("h_energyHCAL","", 100, 0, m_energy);
  m_histograms.push_back(h_energyHCAL);

  h_benchmark = new TH1F("h_benchmark","", 1000, 0, m_energy+0.5*m_energy);
  m_histograms.push_back(h_benchmark);
  
  h_emScale = new TH1F("h_emScale", "", 1000, 0, m_energy+0.5*m_energy);
  m_histograms.push_back(h_emScale);

}

void ChiSquareMinimisationEndcap::processEvent(podio::EventStore& aStore, int aEventId, int aEnergy, bool aVerbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::PositionedCaloHitCollection*    colECalPositions(nullptr);
  const fcc::PositionedCaloHitCollection*    colCalPositions(nullptr);

  bool colMCParticlesOK   = aStore.get("GenParticles", colMCParticles);
  bool colECalPositionsOK = aStore.get("cellECalPositions" , colECalPositions);
  bool colCalPositionsOK  = aStore.get("cellHCalPositions" , colCalPositions);

  // Energy depsitied in 1st HCAL layer
  E_firstLayer = 0.;
  E_lastLayer = 0.;
  //Total cell energy per event in ECal
  SumE_cell = 0.;
  //Total hit energy per event in HCal
  SumH_cell = 0.;

  EEndcap.push_back(aEnergy);

  int entriesLastEcalLayer=0;

  //PositionedHits collection
  if (colECalPositionsOK && colCalPositionsOK) {
    if (aVerbose) {
      std::cout << " Collections:            " << std::endl;
      std::cout << " -> #ECalPositions:      " << colECalPositions->size()   << std::endl;
      std::cout << " -> #HCalPositionedHits: " << colCalPositions->size()    << std::endl;
    }

    // loop over ECAL entries)
    for (auto& iecl=colECalPositions->begin(); iecl!=colECalPositions->end(); ++iecl)
      {
	double cellEnergy = iecl->core().energy;
	SumE_cell += cellEnergy;
	TVector3 cellPosition(iecl->position().x,iecl->position().y,iecl->position().z);
	double r = sqrt(pow(iecl->position().x,2)+pow(iecl->position().y,2));
     	int layerId = ecal_decoder.value("layer",iecl->core().cellId);
	if (aVerbose) std::cout << "ECAL cell layer: " << layerId << std::endl;
	if ( layerId == 7 ){	
     	  E_lastLayer += cellEnergy;
	  entriesLastEcalLayer ++;
	}
      }

    // loop over HCal entries
    for (auto& iecl=colCalPositions->begin(); iecl!=colCalPositions->end(); ++iecl)
      {
	if (iecl->core().energy > m_thr){
	  double hitEnergy = iecl->core().energy;
	  TVector3 hitPosition(iecl->position().x,iecl->position().y,iecl->position().z);
	  double r = sqrt(pow(iecl->position().x,2)+pow(iecl->position().y,2));
	  int layerId = hcal_decoder.value("layer",iecl->core().cellId);
	  // Add energy in first HCal layer
	  if ( layerId == 0 ){	
	    E_firstLayer += hitEnergy;
	  }
	  SumH_cell += hitEnergy;
	}
      } 
    // Fill vectors
    vecEemEndcap.push_back(SumE_cell); 
    vecEhadEndcap.push_back(SumH_cell); 
    vecEfirstEndcap.push_back(E_firstLayer);
    vecElastEndcap.push_back(E_lastLayer);

    // Fill histograms
    h_energyECAL ->Fill(SumE_cell);
    h_energyHCAL ->Fill(SumH_cell);

    if (aVerbose) 
      std::cout << "Total cell energy, ECAL (GeV): " << SumE_cell << " total hit energy, HCAL (GeV): " << SumH_cell << std::endl;
    if (entriesLastEcalLayer==0) std::cout << "Attention!!! no ecal layer entries!!!" << std::endl;

    // Calibration to EM scale
    h_emScale->Fill(SumE_cell + SumH_cell);  

    // Benchmark reconstruction
    double E0_rec = SumE_cell*a + SumH_cell + c*sqrt(abs(E_lastLayer*a*E_firstLayer)) + d*pow(SumE_cell*a,2);
    h_benchmark ->Fill(E0_rec);
  }
  else {
    if (aVerbose) {
      std::cout << "No CaloPositionedHits Collection!!!!!" << std::endl;
    }
  }
}


void ChiSquareMinimisationEndcap::finishLoop(int aNumEvents, bool aVerbose) {
    
  std::cout << "Running minimisation for " << EEndcap.size() << " #events! \n";
  if (EEndcap.size() != vecEhadEndcap.size()){
    std::cout << "Something's wrong!!" << std::endl; 
  } 
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  
  //  ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );  
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
    
  // create funciton wrapper for minmizer
  // a IMultiGenFunction type 
  ROOT::Math::Functor f(&chiSquareFitEndcap,5); 
  Double_t steps[5] = {0.01, 0.01, 0.01, 0.0001, 0.01};
  Double_t step = {0.00001};
  // starting point
  double variable[5] = {1.,1.,0.74,-0.0038, 1.};
  min->SetFunction(f);
	  
  // Set the free variables to be minimized!
  //  min->SetFixedVariable(0,"a", variable[0] ); 
  min->SetVariable(0,"a", variable[0], steps[0] ); 
  min->SetVariable(1,"b", variable[1], steps[1] ); 
  min->SetVariable(2,"c", variable[2], steps[2] ); 
  min->SetVariable(3,"d", variable[3], steps[3] );
  min->SetVariable(4,"e", variable[4], steps[4] );
  min->FixVariable(1);
  min->FixVariable(4);
  // Set parameter Limits 
  min->SetVariableLimits(0, 0, 1e6 ); 
  min->SetVariableLimits(2, 0, 1e6 ); 
  min->SetVariableLimits(3,-1e6, 1e6 ); 
	  
  // do the minimization
  min->Minimize(); 
	  
  const Double_t *xs = min->X();
  const Double_t *ys = min->Errors();
  std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] <<  "," << xs[3] << "," << xs[4] << "): " 
	    << min->MinValue()  << std::endl;
  std::cout << "Minimum: #delta f(" << ys[0] << "," << ys[1] << "," << ys[2] << "," << ys[3] << "," << ys[4] << std::endl;

  h_parameter->SetBinContent(0, xs[0]);
  h_parameter->SetBinContent(1, xs[1]);
  h_parameter->SetBinContent(2, xs[2]);
  h_parameter->SetBinContent(3, xs[3]);
  h_parameter->SetBinContent(4, xs[4]);

  h_parameter->SetBinError(0, ys[0]);
  h_parameter->SetBinError(1, ys[1]);
  h_parameter->SetBinError(2, ys[2]);
  h_parameter->SetBinError(3, ys[3]);
  h_parameter->SetBinError(4, ys[4]);
}


