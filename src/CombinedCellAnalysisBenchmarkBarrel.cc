#include "CombinedCellAnalysisBenchmarkBarrel.h"

// FCC-EDM
#include "datamodel/MCParticleCollection.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/PositionedCaloHitCollection.h"

// PODIO
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

// ROOT
#include "TMath.h"
#include "TVector3.h"

// STL
#include <iostream>
#include <string>
#include <sstream>
#include <bitset>

CombinedCellAnalysisBenchmarkBarrel::CombinedCellAnalysisBenchmarkBarrel(double aEnergy, double aThr, double aA, double aB, double aC, double aLinA, double aLinB, double aLinC, double aA0, double aA1, double aA2, double aB0, double aB1, double aB2, double aC0, double aC1, double aC2, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal) : m_energy(aEnergy), m_thr(aThr), m_a(aA), m_b(aB), m_c(aC), m_aLin(aLinA), m_bLin(aLinB), m_cLin(aLinC), m_a0(aA0), m_a1(aA1), m_a2(aA2), m_b0(aB0), m_b1(aB1), m_b2(aB2), m_c0(aC0), m_c1(aC1), m_c2(aC2), ecal_decoder(aBitfieldEcal), hcal_decoder(aBitfieldHcal) {
  Initialize_histos();
}

CombinedCellAnalysisBenchmarkBarrel::~CombinedCellAnalysisBenchmarkBarrel() {}

void CombinedCellAnalysisBenchmarkBarrel::Initialize_histos() {
  double sigma = sqrt(pow(0.6*sqrt(m_energy)+0.02*(m_energy),2));

  h_eTot = new TH1D("h_eTot","", 4*std::ceil(4*sigma)+std::ceil(3*sigma), std::floor(m_energy-4*sigma), std::ceil(m_energy+3*sigma));
  m_histograms.push_back(h_eTot);
 
  h_classicBenchmark = new TH1D("h_classicBenchmark","", 4*std::ceil(4*sigma)+std::ceil(3*sigma), std::floor(m_energy-4*sigma), std::ceil(m_energy+3*sigma));
  m_histograms.push_back(h_classicBenchmark);

  h_eDepBenchmark = new TH1D("h_eDepBenchmark","", 4*std::ceil(4*sigma)+std::ceil(3*sigma), std::floor(m_energy-4*sigma), std::ceil(m_energy+3*sigma));
  m_histograms.push_back(h_eDepBenchmark);
}

void CombinedCellAnalysisBenchmarkBarrel::processEvent(podio::EventStore& aStore, int aEventId, bool verbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::PositionedCaloHitCollection*     colECalPositions(nullptr);
  const fcc::PositionedCaloHitCollection*     colHCalPositions(nullptr);


  bool colMCParticlesOK       = aStore.get("GenParticles", colMCParticles);
  // for Pb spacer 
  bool colECalPositionsOK     = aStore.get("cellECalPositions" , colECalPositions);
  bool colHCalPositionsOK     = aStore.get("cellHCalPositions" , colHCalPositions);

  //Total hit energy per event
  SumE_ecal = 0.;
  SumE_hcal = 0.;
  E_firstLayer = 0.;
  E_lastLayer = 0.;
  double cryoE = 0.;

  //Cells in ECal collection
  if (colECalPositionsOK && colHCalPositionsOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #newCaloPositions:     " << colECalPositions->size()    << std::endl;;
      std::cout << " -> #newCaloPositions:     " << colHCalPositions->size()    << std::endl;;
    }
    for (auto& iecl=colECalPositions->begin(); iecl!=colECalPositions->end(); ++iecl){
      SumE_ecal += iecl->core().energy;
      TVector3 hitPosition(iecl->position().x,iecl->position().y,iecl->position().z);
      double r = sqrt(pow(iecl->position().x,2)+pow(iecl->position().y,2));
      TVector3 vec(iecl->position().x,iecl->position().y,iecl->position().z);
      double phi = atan2( iecl->position().y, iecl->position().x );
      double eta = vec.Eta();
      int layerId = ecal_decoder.value("layer",iecl->core().cellId);
      int cryoId = ecal_decoder.value("cryo",iecl->core().cellId);
      
      // Add energy in last ECal layer
      if ( layerId == 7 ){	
	E_lastLayer += iecl->core().energy;
      }
      if (verbose){
	std::cout << " ECAL :  \n";
	std::cout << " r " << r << ", eta " << eta << ", phi " << phi << ", energy " << iecl->core().energy << "\n";
	std::cout << " layer " << layerId << std::endl;
      }
      
      if (cryoId==1){
	// Fill histogram for correlation plot
	// Re-calibrate the enegery in the cryostat.. due to layerId0 -> calibrated to EM scale of first ecal layer
	cryoE += iecl->core().energy * 0.12125;
	continue;
      } 
    }
    for (auto& iecl=colHCalPositions->begin(); iecl!=colHCalPositions->end(); ++iecl){
      if (iecl->core().energy > m_thr){
	SumE_hcal += iecl->core().energy;
	TVector3 hitPosition(iecl->position().x,iecl->position().y,iecl->position().z);
	double r = sqrt(pow(iecl->position().x,2)+pow(iecl->position().y,2));
	TVector3 vec(iecl->position().x,iecl->position().y,iecl->position().z);
	double phi = atan2( iecl->position().y, iecl->position().x );
	double eta = vec.Eta();
	int layerId = hcal_decoder.value("layer",iecl->core().cellId);

	// Add energy in first HCal layer
	if ( layerId == 0 ){	
	  E_firstLayer += iecl->core().energy;
	}
	if (verbose){
	  std::cout << "HCAL :  \n";
	  std::cout << " r " << r << ", eta " << eta << ", phi " << phi << ", energy " << iecl->core().energy << "\n";
	  std::cout << " layer " << layerId << std::endl;
	}
      }
    }
   
    if (verbose) {
      std::cout << "Total cell energy:                 " << SumE_ecal + SumE_hcal << std::endl;
    }


    //Fill histograms
    h_eTot            ->Fill(SumE_ecal+SumE_hcal);
    
    // Etot
    double Etot = SumE_ecal + SumE_hcal;
    // Benchmark reconstruction
    double E_rec0 = SumE_ecal*m_a + SumE_hcal + m_b*sqrt(abs(E_firstLayer*m_a*E_lastLayer)) + m_c*pow(SumE_ecal*m_a,2);
    // energy dependent Benchmark reconstruction
    double parA = m_a0 + m_a1*Etot + m_a2/sqrt(Etot);
    double parB = m_b0 + m_b1*Etot + m_b2*log(Etot);
    double parC = m_c0*Etot+m_c1/pow(Etot,m_c2);
    double E_rec1 = SumE_ecal*parA + SumE_hcal + parB*sqrt(abs(E_firstLayer*parA*E_lastLayer)) + parC*pow(SumE_ecal*parA,2);
    
    h_eTot             ->Fill( Etot );
    h_classicBenchmark ->Fill( pow((E_rec0 - m_aLin)/m_bLin, 1/m_cLin) );
    h_eDepBenchmark    ->Fill( E_rec1 );
  }
  else {
    std::cout << "No colECalPositions Collection!!!!!" << std::endl;
    std::cout << "No colHCalPositions Collection!!!!!" << std::endl;
  } 
}
void CombinedCellAnalysisBenchmarkBarrel::finishLoop(int aNumEvents, bool aVerbose) {
  std::cout << "Total energy: " << h_eTot->GetMean() << ", above threshold " << m_thr << std::endl;  
}
