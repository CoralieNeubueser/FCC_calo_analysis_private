#include "ChiSquareMinimisationClusterBarrel.h"

#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenVertexCollection.h"
#include "datamodel/CaloClusterCollection.h"

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

bool eDep = false;

std::vector<double> vecEClusterBarrel;

std::vector<double> vecEMcluster;
std::vector<double> vecHADcluster;

std::vector<double> vecEMlikeClusterECal;
std::vector<double> vecEMlikeClusterHCal;
std::vector<double> vecEMlikeClusterlastECal;
std::vector<double> vecEMlikeClusterfirstHCal;

std::vector<double> vecHADlikeClusterECal;
std::vector<double> vecHADlikeClusterHCal;
std::vector<double> vecHADlikeClusterlastECal;
std::vector<double> vecHADlikeClusterfirstHCal;

Double_t chiSquareFitClusterBarrel(const Double_t *par) {
  Double_t fitvalue = 0.;
 
  //loop over all events
  for(uint i=0; i<vecEClusterBarrel.size(); i++){

    double E_Cluster = vecEClusterBarrel.at(i);

    double E_EMcluster = vecEMcluster.at(i);
    double E_HADcluster = vecHADcluster.at(i);

    double E_EMlikeCluster_ecal = vecEMlikeClusterECal.at(i);
    double E_EMlikeCluster_hcal = vecEMlikeClusterHCal.at(i);
    double E_EMlikeCluster_lastECal = vecEMlikeClusterlastECal.at(i);
    double E_EMlikeCluster_firstHCal = vecEMlikeClusterfirstHCal.at(i);

    double E_HADlikeCluster_ecal = vecHADlikeClusterECal.at(i);
    double E_HADlikeCluster_hcal = vecHADlikeClusterHCal.at(i);
    double E_HADlikeCluster_lastECal = vecHADlikeClusterlastECal.at(i);
    double E_HADlikeCluster_firstHCal = vecHADlikeClusterfirstHCal.at(i);

    double Etot_EM = E_EMlikeCluster_ecal + E_EMlikeCluster_hcal / 1.; // calibrate HCal cells to EM scale
    double Etot_HAD = E_HADlikeCluster_ecal + E_HADlikeCluster_hcal;

    if (Etot_EM + Etot_HAD + E_HADlikeCluster_hcal + E_EMlikeCluster_ecal == 0)
      continue;
    double ParA1 = 0;
    double ParB1 = 0;
    double ParA2 = 0;
    double ParB2 = 0;
    double ParC = 0;
    
    if (Etot_EM != 0){
      ParA1 = par[0] + par[1]*Etot_EM + par[2]/sqrt(Etot_EM);
      ParB1 = par[3] + pow(Etot_EM,par[4]) + par[5]*log(Etot_EM);
    }
    if (Etot_HAD != 0 && eDep){
      ParA2 = par[0] + par[1]*Etot_HAD + par[2]/sqrt(Etot_HAD);
      ParB2 = par[3] + pow(Etot_HAD,par[4]) + par[5]*log(Etot_HAD); // pow(Etot_HAD,par[4]) + par[5]*log(Etot_HAD);
      //      ParC = par[6]*Etot_HAD + par[7]/pow(Etot_HAD, par[8]); ///exp(Etot_HAD) + par[7]/pow(Etot_HAD,par[8]);
    }

//    // Minimise 3 parameter for each energy seperately:
//    ParA1 = par[0];
//    ParB1 = par[1];
//    ParA2 = par[0];
//    ParB2 = par[1];
//    ParC = par[2];

    // correction on EM like cluster
    Double_t E0_EMlikeCluster = E_EMlikeCluster_ecal*ParA1 + E_EMlikeCluster_hcal + ParB1*sqrt(fabs(E_EMlikeCluster_lastECal*ParA1*E_EMlikeCluster_firstHCal));
    // correction on hadron like cluster
    Double_t E0_HADlikeCluster =  E_HADlikeCluster_ecal*ParA2 + E_HADlikeCluster_hcal + ParB2*sqrt(fabs(E_HADlikeCluster_lastECal*ParA2*E_HADlikeCluster_firstHCal)); // + ParC*pow(E_HADlikeCluster_ecal*ParA2, 2);
    
    if (!eDep){
      E0_HADlikeCluster = E_HADlikeCluster_ecal*par[0] + E_HADlikeCluster_hcal + par[1]*sqrt(fabs(E_HADlikeCluster_lastECal*par[0]*E_HADlikeCluster_firstHCal)) + par[2]*pow(E_HADlikeCluster_ecal*par[0], 2);
      E0_EMlikeCluster = 0.;
      E_EMcluster = 0.;
      E_HADcluster = 0.;
    }

//    std::cout << "True Energy                : " << E_Cluster << std::endl;
//    std::cout << "Energy in EM cluster       : " << E_EMcluster << std::endl;
//    std::cout << "Energy in HAD cluster      : " << E_HADcluster << std::endl;
//    std::cout << "Energy in EM like cluster  : " << E0_EMlikeCluster << std::endl;
//    std::cout << "Energy in HAD like cluster : " << E0_HADlikeCluster << std::endl;
    Double_t Etot_cluster = E_EMcluster + E_HADcluster + E0_EMlikeCluster + E0_HADlikeCluster;
    
    fitvalue += pow((E_Cluster - Etot_cluster),2)/E_Cluster; 
  }
  return fitvalue;    
}

ChiSquareMinimisationClusterBarrel::ChiSquareMinimisationClusterBarrel(double aEnergy, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal, bool aEnergyDependence):
  m_energy(aEnergy), m_decoderECal(aBitfieldEcal), m_decoderHCal(aBitfieldHcal), m_eDep(aEnergyDependence) {
  Initialize_histos();
}

ChiSquareMinimisationClusterBarrel::~ChiSquareMinimisationClusterBarrel() {}

void ChiSquareMinimisationClusterBarrel::Initialize_histos()
{
  eDep = m_eDep;

  vecEClusterBarrel.clear();

  vecEMcluster.clear();
  vecHADcluster.clear();
  
  vecEMlikeClusterECal.clear();
  vecEMlikeClusterHCal.clear();
  vecEMlikeClusterlastECal.clear();
  vecEMlikeClusterfirstHCal.clear();
  
  vecHADlikeClusterECal.clear();
  vecHADlikeClusterHCal.clear();
  vecHADlikeClusterlastECal.clear();
  vecHADlikeClusterfirstHCal.clear();

  h_parameter = new TH1F("h_parameter","", 15, 0,16);
  m_histograms.push_back(h_parameter);

  h_energyECAL = new TH1F("h_energyECAL","", 15000, 0, m_energy+0.5*m_energy);
  m_histograms.push_back(h_energyECAL);

  h_energyHCAL = new TH1F("h_energyHCAL","", 15000, 0, m_energy+0.5*m_energy);
  m_histograms.push_back(h_energyHCAL);

  h_energyHADcluster = new TH1F("h_energyHADcluster","", 15000, 0, m_energy+0.5*m_energy);
  m_histograms.push_back(h_energyHADcluster);

  h_energyEMcluster = new TH1F("h_energyEMcluster","", 15000, 0, m_energy+0.5*m_energy);
  m_histograms.push_back(h_energyEMcluster);

  h_benchmark = new TH1F("h_benchmark","", 150000, 0, m_energy+0.5*m_energy);
  m_histograms.push_back(h_benchmark);
  
  h_emScale = new TH1F("h_emScale", "", 150000, 0, m_energy+0.5*m_energy);
  m_histograms.push_back(h_emScale);

}

void ChiSquareMinimisationClusterBarrel::processEvent(podio::EventStore& aStore, int aEventId, double aEnergy, bool aVerbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::CaloClusterCollection* clusters(nullptr);

  bool colMCParticlesOK   = aStore.get("GenParticles", colMCParticles);
  bool colClustersOK = aStore.get("caloClustersBarrelNoise" , clusters);

  double EMcluster = 0.;
  double HADcluster = 0.;
  
  double EMlikeClusterECal = 0.;
  double EMlikeClusterHCal = 0.;
  double EMlikeClusterlastECal = 0.;
  double EMlikeClusterfirstHCal = 0.;
  
  double HADlikeClusterECal = 0.;
  double HADlikeClusterHCal = 0.;
  double HADlikeClusterlastECal = 0.;
  double HADlikeClusterfirstHCal = 0.;
 
  //PositionedHits collection
  if (colClustersOK) {
    if (aVerbose) {
      std::cout << " Collections:            " << std::endl;
      std::cout << " -> #BarrelCluster:      " << clusters->size()   << std::endl;
    }
    // loop over Cluster
    for (auto& iecl=clusters->begin(); iecl!=clusters->end(); ++iecl)
      {
	auto cluster = *iecl;
	bool cellsInBoth = false;
	std::map<uint,double> energyBoth;
	// sum energies in E and HCal layers for benchmark correction
	double energyLastECal = 0.;
	double energyFirstHCal = 0.;
	
	// Loop over cluster cells 
	for (uint it = 0; it < cluster.hits_size(); it++){
	  auto cellId = cluster.hits(it).core().cellId;
	  auto cellEnergy = cluster.hits(it).core().energy;
	  m_decoder->setValue(cellId);
	  uint systemId = m_decoder->value("system",cellId);
	  int layerId;
	  if (systemId == m_systemIdECal){
	    m_decoderECal.setValue(cellId);
	    layerId = m_decoderECal.value("layer",cellId);
	    if (aVerbose) {
	      std::cout << "ECal layer: " << layerId << std::endl;
	    }	    
	    if( layerId == m_lastECalLayer ) 
	      energyLastECal += cellEnergy;
	  }
	  else if  (systemId == m_systemIdHCal){
	    m_decoderHCal.setValue(cellId);
	    layerId = m_decoderHCal.value("layer",cellId);
	    if (aVerbose) {
	      std::cout << "HCal layer: " << layerId << std::endl;
	    }
	    if ( layerId == m_firstHCalLayer )
	      energyFirstHCal += cellEnergy;
	  }
	  energyBoth[systemId] += cellEnergy;
	}
	
	if (energyBoth[m_systemIdHCal] > 1e-3 &&  energyBoth[m_systemIdECal] > 1e-3)
	  cellsInBoth = true;
	
	// check if cluster energy is equal to sum over cells
	if (static_cast<int>(iecl->core().energy*100.0) != static_cast<int>((energyBoth[m_systemIdECal] + energyBoth[m_systemIdHCal])*100.0))
	  if (aVerbose) {
	    std::cout << "The cluster energy is not equal to sum over cell energy: " << iecl->core().energy << ", " << (energyBoth[m_systemIdECal] + energyBoth[m_systemIdHCal]) << std::endl;
	  }
	// 2. Calibrate the cluster if it contains cells in both systems
	if(cellsInBoth || !m_eDep) {
	  // Calculate the fraction of energy in ECal
	  auto energyFraction = energyBoth[m_systemIdECal] / iecl->core().energy;
	  if (aVerbose) {
	    std::cout << "Energy fraction in ECal : " << energyFraction << std::endl;
	  }
	  bool calibECal = false;
	  if (energyFraction >= m_fractionECal && m_eDep) {
	    // calibrate HCal cells to EM scale
	    // assuming HCal cells are calibrated to hadron scale
	    EMlikeClusterECal += energyBoth[m_systemIdECal];
	    EMlikeClusterHCal += energyBoth[m_systemIdHCal] * m_ehHCal;
	    EMlikeClusterlastECal += energyLastECal;
	    EMlikeClusterfirstHCal += energyFirstHCal;
	  }
	  else {
	    // calibrate ECal cells to hadron scale
	    // assuming ECal cells are calibrated to EM scale
	    HADlikeClusterECal += energyBoth[m_systemIdECal];
	    HADlikeClusterHCal += energyBoth[m_systemIdHCal];
	    HADlikeClusterlastECal += energyLastECal;
	    HADlikeClusterfirstHCal += energyFirstHCal;
	  }
	}
	else if (energyBoth[m_systemIdECal] > 1e-3 ){ // Fill the unchanged cluster in output collection
	  EMcluster += energyBoth[m_systemIdECal];
	}
	else if (energyBoth[m_systemIdHCal] > 1e-3 ){ // Fill the unchanged cluster in output collection
	  HADcluster += energyBoth[m_systemIdHCal];
	}
      } 
    // Fill vectors
    vecEClusterBarrel.push_back(aEnergy);
    
    vecEMcluster.push_back( EMcluster );
    vecHADcluster.push_back( HADcluster );
    
    vecEMlikeClusterECal.push_back( EMlikeClusterECal );
    vecEMlikeClusterHCal.push_back( EMlikeClusterHCal );
    vecEMlikeClusterlastECal.push_back( EMlikeClusterlastECal );
    vecEMlikeClusterfirstHCal.push_back( EMlikeClusterfirstHCal );
    
    vecHADlikeClusterECal.push_back( HADlikeClusterECal );
    vecHADlikeClusterHCal.push_back( HADlikeClusterHCal );
    vecHADlikeClusterlastECal.push_back( HADlikeClusterlastECal );
    vecHADlikeClusterfirstHCal.push_back( HADlikeClusterfirstHCal );
    
    // fill histograms for 4 types of clusters
    h_energyEMcluster->Fill(EMlikeClusterECal + EMlikeClusterHCal/1.);
    h_energyHADcluster->Fill(HADlikeClusterECal + HADlikeClusterHCal);
    h_energyECAL->Fill(EMcluster);
    h_energyHCAL->Fill(HADcluster);

    if (aVerbose) 
      std::cout << "Total cluster energy, ECAL (GeV): " << EMcluster << " total hit energy, HCAL (GeV): " << HADcluster << std::endl;
    if ((EMlikeClusterlastECal+HADlikeClusterlastECal)==0) std::cout << "Attention!!! no ecal layer entries!!!" << std::endl;
    
  }
  else {
    if (aVerbose) {
      std::cout << "No Cluster Collection!!!!!" << std::endl;
    }
  }
}


void ChiSquareMinimisationClusterBarrel::finishLoop(int aNumEvents, bool aVerbose) {
    
  int numPar = 3;
  if (m_eDep)
    numPar = 6;

  std::cout << "Running minimisation for " << vecEClusterBarrel.size() << " #events! \n";
  if (vecEClusterBarrel.size() != vecEMcluster.size()){
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
  ROOT::Math::Functor f(&chiSquareFitClusterBarrel, numPar); 
  Double_t steps[9] = {0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
  // starting point
  double variable[9] = {0.957055625701,3.77259715365e-08,1.02754379428,0.24268958826,0.0969499656497,-0.216931357087,-1.54642630749e-10,-4519735.41454,7.02561585394};
  min->SetFunction(f);
	  
  // Set the free variables to be minimized!
  //  min->SetFixedVariable(0,"a", variable[0] ); 
  min->SetVariable(0,"a1", variable[0], steps[0] ); 
  min->SetVariable(1,"a2", variable[1], steps[1] ); 
  min->SetVariable(2,"a3", variable[2], steps[2] ); 

  if (numPar>3){
    min->SetVariable(3,"b1", variable[3], steps[3] ); 
    min->SetVariable(4,"b2", variable[4], steps[4] ); 
    min->SetVariable(5,"b3", variable[5], steps[5] ); 
  }
//  min->SetVariable(6,"c1", variable[6], steps[6] ); 
//  min->SetVariable(7,"c2", variable[7], steps[7] ); 
//  min->SetVariable(8,"c3", variable[8], steps[8] ); 

//  min->SetVariable(9,"a1em", variable[0], steps[0] ); 
//  min->SetVariable(10,"a2em", variable[1], steps[1] ); 
//  min->SetVariable(11,"a3em", variable[2], steps[2] ); 
//  min->SetVariable(12,"b1em", variable[3], steps[3] ); 
//  min->SetVariable(13,"b2em", variable[4], steps[4] ); 
//  min->SetVariable(14,"b3em", variable[5], steps[5] ); 
	   
  // do the minimization
  min->Minimize(); 
	  
  const Double_t *xs = min->X();
  const Double_t *ys = min->Errors();
  std::cout << "Minimum: f("; 
  for (uint i=0;i<numPar;i++){   
    std::cout << "," << xs[i] << "\n";
  }
  std::cout << "Minimum: #delta f("; 
  for (uint i=0;i<numPar;i++){
    std::cout << "," << ys[i] << "\n";
  }
    
  h_parameter->SetBinContent(0, xs[0]);
  h_parameter->SetBinContent(1, xs[1]);
  h_parameter->SetBinContent(2, xs[2]);
  if(m_eDep){
    h_parameter->SetBinContent(3, xs[3]);
    h_parameter->SetBinContent(4, xs[4]);
    h_parameter->SetBinContent(5, xs[5]);
  }
//  h_parameter->SetBinContent(6, xs[6]);
//  h_parameter->SetBinContent(7, xs[7]);
//  h_parameter->SetBinContent(8, xs[8]);

//  h_parameter->SetBinContent(9, xs[9]);
//  h_parameter->SetBinContent(10, xs[10]);
//  h_parameter->SetBinContent(11, xs[11]);
//  h_parameter->SetBinContent(12, xs[12]);
//  h_parameter->SetBinContent(13, xs[13]);
//  h_parameter->SetBinContent(14, xs[14]);
 
  h_parameter->SetBinError(0, ys[0]);
  h_parameter->SetBinError(1, ys[1]);
  h_parameter->SetBinError(2, ys[2]);
  if(m_eDep){
    h_parameter->SetBinError(3, ys[3]);
    h_parameter->SetBinError(4, ys[4]);
    h_parameter->SetBinError(5, ys[5]);
  }
//  h_parameter->SetBinError(6, ys[6]);
//  h_parameter->SetBinError(7, ys[7]);
//  h_parameter->SetBinError(8, ys[8]);

//  h_parameter->SetBinError(9, ys[9]);
//  h_parameter->SetBinError(10, ys[10]);
//  h_parameter->SetBinError(11, ys[11]);
//  h_parameter->SetBinError(12, ys[12]);
//  h_parameter->SetBinError(13, ys[13]);
//  h_parameter->SetBinError(14, ys[14]);
}


