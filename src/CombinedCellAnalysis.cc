#include "CombinedCellAnalysis.h"

// FCC-EDM
#include "datamodel/MCParticleCollection.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/PositionedCaloHitCollection.h"

// PODIO
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

// ROOT
#include "TVector3.h"

// STL
#include <iostream>

CombinedCellAnalysis::CombinedCellAnalysis(double aEnergy, double aA1, double aA2, double aA3, double aB, double aC1, double aC2, double aC3, double aD, double aLinCorrA, double aLinCorrB, double aLinCorrC, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal) : m_energy(aEnergy), m_a1(aA1), m_a2(aA2), m_a3(aA3), m_b(aB), m_c1(aC1), m_c2(aC2), m_c3(aC3), m_d(aD), m_linCorrA(aLinCorrA), m_linCorrB(aLinCorrB), m_linCorrC(aLinCorrC), ecal_decoder(aBitfieldEcal), hcal_decoder(aBitfieldHcal), SumE_ecal(0), SumE_hcal(0), h_cellEnergy(nullptr) {
  Initialize_histos();
}

CombinedCellAnalysis::~CombinedCellAnalysis() {}


void CombinedCellAnalysis::Initialize_histos() {
  double sigma = sqrt(pow(0.6*sqrt(m_energy)+0.1*(m_energy),2));
  h_benchmark = new TH1F("h_benchmark","", 1000,0,2);
  m_histograms.push_back(h_benchmark);

  h_benchmarkLinCorr = new TH1F("h_benchmarkLinCorr","", 1000,0,2);
  m_histograms.push_back(h_benchmarkLinCorr);

  h_cellEnergy = new TH1F("h_cellEnergy","", 1000,0,2);
  m_histograms.push_back(h_cellEnergy);

  h_numCells = new TH1F("h_numCells","", 1000, 0, 1000000);
  m_histograms.push_back(h_numCells);

  h_cellEnergy_ecal = new TH1F("h_cellEnergy_ecal","",1000,0,2);
  m_histograms.push_back(h_cellEnergy_ecal);

  h_cellEnergy_hcal = new TH1F("h_cellEnergy_hcal","",1000,0,2);
  m_histograms.push_back(h_cellEnergy_hcal);

}

void CombinedCellAnalysis::processEvent(podio::EventStore& aStore, int aEventId, bool verbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::CaloHitCollection*     colECal(nullptr);
  const fcc::CaloHitCollection*     colHCal(nullptr);

  bool colMCParticlesOK       = aStore.get("GenParticles", colMCParticles);
  bool colECalOK     = aStore.get("ECalBarrelCells" , colECal);
  bool colHCalOK     = aStore.get("HCalBarrelCells" , colHCal);

  double etaVertex = 0;
  double phiVertex = 0;
  double thetaVertex = 0;

  //MCParticle and Vertices collection
  if (colMCParticlesOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #MCTruthParticles:     " << colMCParticles->size()    << std::endl;
    }
    //Loop through the collection
    for (auto& iparticle=colMCParticles->begin(); iparticle!=colMCParticles->end(); ++iparticle) {
      TVector3 particle(iparticle->core().p4.px,iparticle->core().p4.py,iparticle->core().p4.pz);
      etaVertex = particle.Eta();
      phiVertex = particle.Phi();
      thetaVertex = 2 * atan( exp(-etaVertex) );
      //Fill histograms
      //      h_ptGen->Fill( sqrt( pow(iparticle->core().p4.px,2)+pow(iparticle->core().p4.py,2) ) );
    }
  }
  else {
    if (verbose) {
      std::cout << "No MCTruth info available" << std::endl;
    }
  }

  //Total hit energy per event
  SumE_ecal = 0.;
  SumE_hcal = 0.;
  double phiE_ecal = 0.;
  double phiE_hcal = 0.;
  double etaE_ecal = 0.;
  double etaE_hcal = 0.;
  double thetaE_ecal = 0.;
  double zE_ecal = 0.;
  double thetaE_hcal = 0.;
  double zE_hcal = 0.;
  E_firstLayer = 0.;
  E_lastLayer = 0.;

  //Cells in ECal collection
  if (colECalOK && colHCalOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #newCaloPositions:     " << colECal->size()    << std::endl;;
      std::cout << " -> #newCaloPositions:     " << colHCal->size()    << std::endl;;
    }
    for (auto& iecl=colECal->begin(); iecl!=colECal->end(); ++iecl){
      SumE_ecal += iecl->core().energy;
      auto layer = ecal_decoder.value("layer",iecl->core().cellId);
      if (layer == 7)
	E_lastLayer += iecl->core().energy;
    }
    for (auto& iecl=colHCal->begin(); iecl!=colHCal->end(); ++iecl){
      SumE_hcal += iecl->core().energy;
      auto layer = hcal_decoder.value("layer",iecl->core().cellId);
      if (layer == 0)
	E_firstLayer += iecl->core().energy;
    }
   
    if (verbose) {
      std::cout << "Total cell energy:                 " << SumE_ecal + SumE_hcal << std::endl;
    }


    //Fill histograms
    h_cellEnergy      ->Fill((SumE_ecal+SumE_hcal)/m_energy);
    h_cellEnergy_ecal ->Fill(SumE_ecal/m_energy);
    h_cellEnergy_hcal ->Fill(SumE_hcal/m_energy);
    h_numCells ->Fill(colECal->size() + colHCal->size());
    
    // Benchmark reconstruction
    double parA = m_a1 + m_a2*(SumE_ecal+SumE_hcal) + m_a3/sqrt(SumE_ecal+SumE_hcal);
    if (m_a2 == 0 && m_a3 == 0 )
      parA = m_a1;
    double parC = m_c1 + pow((SumE_ecal+SumE_hcal), m_c2) + m_c3*log(SumE_ecal+SumE_hcal); 
    if (m_c2 == 0 && m_c3 == 0 )
      parC = m_c1;

    double E0_rec = SumE_ecal*parA + SumE_hcal*m_b + parC*sqrt(fabs(E_firstLayer*parA*E_lastLayer*m_b)) + m_d*pow(SumE_ecal*parA,2);
    double Ecorr_rec = pow((E0_rec - m_linCorrA)/m_linCorrB, 1/m_linCorrC);
    h_benchmark          ->Fill(E0_rec/m_energy);
    h_benchmarkLinCorr   ->Fill(Ecorr_rec/m_energy);
  }
  else {
    std::cout << "No colECalPositions Collection!!!!!" << std::endl;
    std::cout << "No colHCalPositions Collection!!!!!" << std::endl;
  } 
}
void CombinedCellAnalysis::finishLoop(int aNumEvents, bool aVerbose) {
  std::cout << "Total energy: " << h_cellEnergy->GetMean() << std::endl;
  aNumEvents = h_cellEnergy->GetEntries();
}
