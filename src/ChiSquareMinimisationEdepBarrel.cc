#include "ChiSquareMinimisationEdepBarrel.h"

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

std::vector<double> vecEBarrelEdep;
std::vector<double> vecEemBarrelEdep;
std::vector<double> vecEhadBarrelEdep;
std::vector<double> vecEHCAL_firstBarrelEdep;
std::vector<double> vecEECAL_firstBarrelEdep;
std::vector<double> vecEECAL_lastBarrelEdep;

Double_t chiSquareFitBarrelEdep(const Double_t *par) {
  Double_t fitvalue = 0.;
 
  //loop over all events
  for(int i=0; i<vecEemBarrelEdep.size(); i++){

    double E_BarrelEdep = vecEBarrelEdep.at(i);
    double Eem_BarrelEdep = vecEemBarrelEdep.at(i);
    double Ehad_BarrelEdep = vecEhadBarrelEdep.at(i);
    double EHfirst_BarrelEdep = vecEHCAL_firstBarrelEdep.at(i);
    double EEfirst_BarrelEdep = vecEECAL_firstBarrelEdep.at(i);
    double EElast_BarrelEdep = vecEECAL_lastBarrelEdep.at(i);
    
    double Etot = Eem_BarrelEdep + Ehad_BarrelEdep; 
    double parA =  par[0] + par[1]*Etot + par[2]/sqrt(Etot);
    double parB =  par[4] + pow(Etot,par[5]) + par[6]*log(Etot);
    //double parC = par[6]*Etot+par[7]/pow(Etot,par[8]);

    Double_t E0_BarrelEdep = Eem_BarrelEdep*parA + Ehad_BarrelEdep*par[3] + parB*sqrt(abs(EElast_BarrelEdep*parA*EHfirst_BarrelEdep*par[3])) + par[7]*pow(Eem_BarrelEdep*parA,2);
    
    fitvalue += pow((E_BarrelEdep-E0_BarrelEdep),2)/E_BarrelEdep; 
  }
  return fitvalue;    
}

ChiSquareMinimisationEdepBarrel::ChiSquareMinimisationEdepBarrel(double aEnergy, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal):
m_energy(aEnergy), ecal_decoder(aBitfieldEcal), hcal_decoder(aBitfieldHcal) {
  Initialize_histos();
}

ChiSquareMinimisationEdepBarrel::~ChiSquareMinimisationEdepBarrel() {}

void ChiSquareMinimisationEdepBarrel::Initialize_histos()
{
  vecEBarrelEdep.clear();
  vecEemBarrelEdep.clear();
  vecEhadBarrelEdep.clear();
  vecEHCAL_firstBarrelEdep.clear();
  vecEECAL_firstBarrelEdep.clear();
  vecEECAL_lastBarrelEdep.clear();

  h_parameter = new TH1F("h_parameter","", 9, 0,9);
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

void ChiSquareMinimisationEdepBarrel::processEvent(podio::EventStore& aStore, int aEventId, double aEnergy, bool aVerbose) {
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
	  int layerId = hcal_decoder.value("layer",iecl->core().cellId);
	  // Add energy in first HCal layer
	  if ( layerId == 0 ){	
	    EHCAL_firstLayer += iecl->core().energy;
	  }
	  SumH_cell += iecl->core().energy;
	}
      } 
         
    // Fill vectors
    vecEBarrelEdep.push_back(aEnergy);
    vecEemBarrelEdep.push_back(SumE_cell); 
    vecEhadBarrelEdep.push_back(SumH_cell); 
    vecEHCAL_firstBarrelEdep.push_back(EHCAL_firstLayer);
    vecEECAL_firstBarrelEdep.push_back(EECAL_firstLayer);
    vecEECAL_lastBarrelEdep.push_back(EECAL_lastLayer);

    // Fill histograms
    h_energyECAL ->Fill(SumE_cell);
    h_energyHCAL ->Fill(SumH_cell);

    if (aVerbose) 
      std::cout << "Total cell energy, ECAL (GeV): " << SumE_cell << " total hit energy, HCAL (GeV): " << SumH_cell << std::endl;
    if (entriesLastEcalLayer==0) std::cout << "Attention!!! no ecal layer entries!!!" << std::endl;

    // Calibration to EM scale
    h_emScale->Fill(SumE_cell + SumH_cell);  

    // Benchmark reconstruction
    double E0_rec = SumE_cell*a + SumH_cell*b + (c*sqrt(abs(EECAL_lastLayer*a*EHCAL_firstLayer*b))) + d*pow(SumE_cell*a,2);
    h_benchmark ->Fill(E0_rec);
  }
  else {
    if (aVerbose) {
      std::cout << "No CaloPositionedHits Collection!!!!!" << std::endl;
    }
  }
}


void ChiSquareMinimisationEdepBarrel::finishLoop(int aNumEvents, bool aVerbose) {
    
  std::cout << "Running minimisation for " << vecEBarrelEdep.size() << " #events! \n";
  if (vecEBarrelEdep.size() != vecEhadBarrelEdep.size()){
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
  ROOT::Math::Functor f(&chiSquareFitBarrelEdep,8); 
  Double_t steps = 0.001;
  // starting point
  //  double variable[9] = {1.,-0.0001,2.,-0.1,-0.00001,0.1,0.00001,-0.5,1.};
  double variable[8] = {0.957055625701,3.77259715365e-08,1.02754379428,1.,0.24268958826,0.0969499656497,-0.216931357087,.0};
  min->SetFunction(f);
	  
  // Set the free variables to be minimized!
  //  min->SetFixedVariable(0,"a", variable[0] ); 
  min->SetVariable(0,"a1", variable[0], steps ); 
  min->SetVariable(1,"a2", variable[1], steps );
  min->SetVariable(2,"a3", variable[2], steps );
  min->SetVariable(3,"b", variable[3], steps ); 
  min->SetVariable(4,"c1", variable[4], steps ); 
  min->SetVariable(5,"c2", variable[5], steps );
  min->SetVariable(6,"c3", variable[6], steps );
  min->SetVariable(7,"d", variable[7], steps );
//  min->SetVariable(6,"c0", variable[6], steps );
//  min->SetVariable(7,"c1", variable[7], steps );
//  min->SetVariable(8,"c2", variable[8], steps );
  min->SetVariableLimits(3, 0.5, 1.5 );
  min->FixVariable(3);
  min->FixVariable(7);

  // do the minimization
  min->Minimize(); 
	  
  const Double_t *xs = min->X();
  const Double_t *ys = min->Errors();
  std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] <<  "," << xs[3] << "," << xs[4] << "," << xs[5] << "," << xs[6] << "," << xs[7] << "," << xs[8] <<"): " 
	    << min->MinValue()  << std::endl;
  std::cout << "Minimum: #delta f(" << ys[0] << "," << ys[1] << "," << ys[2] << "," << ys[3] << "," << ys[4] << "," << ys[5] << "," << ys[6] << "," << ys[7] << "," << ys[8] << std::endl;

  h_parameter->SetBinContent(0, xs[0]);
  h_parameter->SetBinContent(1, xs[1]);
  h_parameter->SetBinContent(2, xs[2]);
  h_parameter->SetBinContent(3, xs[3]);
  h_parameter->SetBinContent(4, xs[4]);
  h_parameter->SetBinContent(5, xs[5]);
  h_parameter->SetBinContent(6, xs[6]);
  h_parameter->SetBinContent(7, xs[7]);
//  h_parameter->SetBinContent(8, xs[8]);

  h_parameter->SetBinError(0, ys[0]);
  h_parameter->SetBinError(1, ys[1]);
  h_parameter->SetBinError(2, ys[2]);
  h_parameter->SetBinError(3, ys[3]);
  h_parameter->SetBinError(4, ys[4]);
  h_parameter->SetBinError(5, ys[5]);
  h_parameter->SetBinError(6, ys[6]);
  h_parameter->SetBinError(7, ys[7]);
//  h_parameter->SetBinError(8, ys[8]);
}


