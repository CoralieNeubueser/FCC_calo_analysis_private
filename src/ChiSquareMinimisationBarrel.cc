#include "ChiSquareMinimisationBarrel.h"

#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenVertexCollection.h"
#include "datamodel/CaloHitCollection.h"

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

std::vector<double> vecEBarrel;
std::vector<double> vecEemBarrel;
std::vector<double> vecEhadBarrel;
std::vector<double> vecEHCAL_firstBarrel;
std::vector<double> vecEECAL_firstBarrel;
std::vector<double> vecEECAL_lastBarrel;

// Sampling fracions for EM scale                                                                                                                                  
double pion_hcal = .024;
double sf_hcal = .0249;

Double_t chiSquareFitBarrel(const Double_t *par) {
  Double_t fitvalue = 0.;
 
  //loop over all events
  for(int i=0; i<vecEemBarrel.size(); i++){

    double E_Barrel = vecEBarrel.at(i);
    double Eem_Barrel = vecEemBarrel.at(i);
    // Calibrate HCal cells to EM scale
    double Ehad_Barrel = vecEhadBarrel.at(i);// *pion_hcal/sf_hcal;
    double EHfirst_Barrel = vecEHCAL_firstBarrel.at(i); // *pion_hcal/sf_hcal;
    double EEfirst_Barrel = vecEECAL_firstBarrel.at(i);
    double EElast_Barrel = vecEECAL_lastBarrel.at(i);
    
    Double_t E0_Barrel = Eem_Barrel*par[0] + Ehad_Barrel*par[1] + par[2]*sqrt(abs(EElast_Barrel*par[0]*EHfirst_Barrel*par[1])) + par[3]*pow(Eem_Barrel*par[0],2);
    
    fitvalue += pow((E_Barrel-E0_Barrel),2)/E_Barrel; 
  }
  return fitvalue;    
}

ChiSquareMinimisationBarrel::ChiSquareMinimisationBarrel(double aEnergy, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal):
m_energy(aEnergy), ecal_decoder(aBitfieldEcal), hcal_decoder(aBitfieldHcal) {
  Initialize_histos();
}

ChiSquareMinimisationBarrel::~ChiSquareMinimisationBarrel() {}

void ChiSquareMinimisationBarrel::Initialize_histos()
{
  vecEBarrel.clear();
  vecEemBarrel.clear();
  vecEhadBarrel.clear();
  vecEHCAL_firstBarrel.clear();
  vecEECAL_firstBarrel.clear();
  vecEECAL_lastBarrel.clear();

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

void ChiSquareMinimisationBarrel::processEvent(podio::EventStore& aStore, int aEventId, double aEnergy, bool aVerbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::CaloHitCollection*    colECalPositions(nullptr);
  const fcc::CaloHitCollection*    colCalPositions(nullptr);

  bool colMCParticlesOK   = aStore.get("GenParticles", colMCParticles);
  bool colECalPositionsOK = aStore.get("ECalBarrelCells" , colECalPositions);
  bool colCalPositionsOK  = aStore.get("HCalBarrelCells" , colCalPositions);

  // Energy depsitied in 1st HCAL layer
  EHCAL_firstLayer = 0.;
  EECAL_firstLayer = 0.;
  EECAL_lastLayer = 0.;
  //Total cell energy per event in ECal
  SumEcorr_cell = 0.;
  SumE_cell = 0.;
  //Total hit energy per event in HCal
  SumH_cell = 0.;

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
	SumE_cell += iecl->core().energy;
	//TVector3 cellPosition(iecl->position().x,iecl->position().y,iecl->position().z);
	//double r = sqrt(pow(iecl->position().x,2)+pow(iecl->position().y,2));
     	int layerId = ecal_decoder.value("layer",iecl->core().cellId);
	if (aVerbose) std::cout << "ECAL cell layer: " << layerId << std::endl;
	if ( layerId == 0 ){	
     	  EECAL_firstLayer += iecl->core().energy;
	}
	else if ( layerId == 7 ){	
     	  EECAL_lastLayer += iecl->core().energy;
	  entriesLastEcalLayer ++;
	}
      }

    // loop over HCal entries
    for (auto& iecl=colCalPositions->begin(); iecl!=colCalPositions->end(); ++iecl)
      {
	if (iecl->core().energy > m_thr){
	  //TVector3 hitPosition(iecl->position().x,iecl->position().y,iecl->position().z);
	  //double r = sqrt(pow(iecl->position().x,2)+pow(iecl->position().y,2));
	  int layerId = hcal_decoder.value("layer",iecl->core().cellId);
	  // Add energy in first HCal layer
	  if ( layerId == 0 ){	
	    EHCAL_firstLayer += iecl->core().energy;
	  }
	  SumH_cell += iecl->core().energy;
	}
      } 
         
    // Fill vectors
    vecEBarrel.push_back(aEnergy);
    vecEemBarrel.push_back(SumE_cell); 
    vecEhadBarrel.push_back(SumH_cell); 
    vecEHCAL_firstBarrel.push_back(EHCAL_firstLayer);
    vecEECAL_firstBarrel.push_back(EECAL_firstLayer);
    vecEECAL_lastBarrel.push_back(EECAL_lastLayer);

    // Fill histograms
    h_energyECAL ->Fill(SumE_cell);
    h_energyHCAL ->Fill(SumH_cell);

    if (aVerbose) 
      std::cout << "Total cell energy, ECAL (GeV): " << SumE_cell << " total hit energy, HCAL (GeV): " << SumH_cell << std::endl;
    if (entriesLastEcalLayer==0) std::cout << "Attention!!! no ecal layer entries!!!" << std::endl;

    // Calibration to EM scale
    h_emScale->Fill(SumE_cell + SumH_cell);  

    // Benchmark reconstruction
    double E0_rec = SumE_cell*a + SumH_cell*b + (c*sqrt(fabs(EECAL_lastLayer*a*EHCAL_firstLayer*b))) + d*pow(SumE_cell*a,2);
    h_benchmark ->Fill(E0_rec);
  }
  else {
    if (aVerbose) {
      std::cout << "No CaloPositionedHits Collection!!!!!" << std::endl;
    }
  }
}


void ChiSquareMinimisationBarrel::finishLoop(int aNumEvents, bool aVerbose) {
    
  std::cout << "Running minimisation for " << vecEBarrel.size() << " #events! \n";
  if (vecEBarrel.size() != vecEhadBarrel.size()){
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
  ROOT::Math::Functor f(&chiSquareFitBarrel,4); 
  Double_t steps[4] = {0.001, 0.001, 0.001, 0.001};
  // starting point
  // use correct EM scale from 2.4->2.48% --> fit does not converge, stick to hadron scale 2.4%  
  double variable[4] = {1., 1., .5, .01};
  min->SetFunction(f);
	  
  // Set the free variables to be minimized!
  //  min->SetFixedVariable(0,"a", variable[0] ); 
  min->SetVariable(0,"a", variable[0], steps[0] ); 
  min->SetVariable(1,"b", variable[1], steps[1] ); 
  min->SetVariable(2,"c", variable[2], steps[2] ); 
  min->SetVariable(3,"d", variable[3], steps[3] );
  // min->SetVariable(4,"e", variable[4], steps[4] );
  min->FixVariable(1); 
  //  min->FixVariable(4); 	  
  //  min->SetVariableLimits(2, 0, 1e6);
   
  // do the minimization
  min->Minimize(); 
	  
  const Double_t *xs = min->X();
  const Double_t *ys = min->Errors();
  std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] <<  "," << xs[3] << "," << xs[4] << "): " 
	    << min->MinValue()  << std::endl;
  std::cout << "Minimum: #delta f(" << ys[0] << "," << ys[1] << "," << ys[2] << "," << ys[3] << "," << xs[4] <<std::endl;

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


