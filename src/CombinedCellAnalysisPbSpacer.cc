#include "CombinedCellAnalysisPbSpacer.h"

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

CombinedCellAnalysisPbSpacer::CombinedCellAnalysisPbSpacer(double aEnergy, double aThr, double aA, double aB, double aC) : m_energy(aEnergy), ecal_decoder(""), hcal_decoder(""), m_thr(aThr), m_a(aA), m_b(aB), m_c(aC) {
  Initialize_histos();
}
CombinedCellAnalysisPbSpacer::CombinedCellAnalysisPbSpacer(double aEnergy, double aThr, double aA, double aB, double aC, double aLinA, double aLinB, double aLinC, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal) : m_energy(aEnergy), m_thr(aThr), m_a(aA), m_b(aB), m_c(aC),m_aLin(aLinA), m_bLin(aLinB), m_cLin(aLinC), ecal_decoder(aBitfieldEcal), hcal_decoder(aBitfieldHcal) {
  Initialize_histos();
}

CombinedCellAnalysisPbSpacer::~CombinedCellAnalysisPbSpacer() {}


void CombinedCellAnalysisPbSpacer::Initialize_histos() {
  double sigma = sqrt(pow(0.8*sqrt(m_energy)+0.02*(m_energy),2));

  h_benchmark = new TH1D("h_benchmark","", 4*std::ceil(4*sigma)+std::ceil(3*sigma), std::floor(m_energy-4*sigma), std::ceil(m_energy+3*sigma));
  m_histograms.push_back(h_benchmark);

  h_benchmark2nd = new TH1D("h_benchmark2nd","", 4*std::ceil(4*sigma)+std::ceil(3*sigma), std::floor(m_energy-4*sigma), std::ceil(m_energy+3*sigma));
  m_histograms.push_back(h_benchmark2nd);

  h_cellEnergy = new TH1D("h_cellEnergy","", 4*std::ceil(4*sigma)+std::ceil(3*sigma), std::floor(m_energy-4*sigma), std::ceil(m_energy+3*sigma));
  m_histograms.push_back(h_cellEnergy);

  h_cellEnergy_ecal = new TH1D("h_cellEnergy_ecal","",100,0,m_energy);
  m_histograms.push_back(h_cellEnergy_ecal);

  h_cellEnergy_hcal = new TH1D("h_cellEnergy_hcal","",100,0,m_energy);
  m_histograms.push_back(h_cellEnergy_hcal);

  h_cellEnergy_ecalP = new TH1D("h_cellEnergy_ecalP","",100,0,m_energy);
  m_histograms.push_back(h_cellEnergy_ecalP);

  h_cellEnergy_hcalP = new TH1D("h_cellEnergy_hcalP","",100,0,m_energy);
  m_histograms.push_back(h_cellEnergy_hcalP);

  h_phiRec_ecal = new TH1D("h_phiRec_ecal","",1000,-.5,.5);
  m_histograms.push_back(h_phiRec_ecal);

  h_phiRec_hcal = new TH1D("h_phiRec_hcal","",1000,-.5,.5);
  m_histograms.push_back(h_phiRec_hcal);

  h_phiRec = new TH1D("h_phiRec","",1000,-.5,.5);
  m_histograms.push_back(h_phiRec);

  h_etaRec = new TH1D("h_etaRec","",1000,-.5,.5);
  m_histograms.push_back(h_etaRec);

  h_etaRec_ecal = new TH1D("h_etaRec_ecal","",1000,-.5,.5);
  m_histograms.push_back(h_etaRec_ecal);

  h_etaRec_hcal = new TH1D("h_etaRec_hcal","",1000,-.5,.5);
  m_histograms.push_back(h_etaRec_hcal);

  h_thetaRec_ecal = new TH1D("h_thetaRec_ecal","",1000,-.5,.5);
  m_histograms.push_back(h_thetaRec_ecal);

  h_thetaRec_hcal = new TH1D("h_thetaRec_hcal","",1000,-.5,.5);
  m_histograms.push_back(h_thetaRec_hcal);

  h_thetaRec = new TH1D("h_thetaRec","",1000,-.5,.5);
  m_histograms.push_back(h_thetaRec);

  h_cellEnergy_first = new TH1D("h_cellEnergy_first","Energy in first HCAL layer",100,0,m_energy/2);
  m_histograms.push_back(h_cellEnergy_first);

  h_cellEnergy_last = new TH1D("h_cellEnergy_last","Energy in last ECAL layer",100,0,m_energy/2);
  m_histograms.push_back(h_cellEnergy_last);

  h_cellId = new TH1D("h_cellId","", 1000, 0,5000e6);
  m_histograms.push_back(h_cellId);

  h_lostECorr = new TH2D("h_lostECorr","", 50,0,50, 50,0,10);
  h_lostEMultiCorr = new TH2D("h_lostEMultiCorr","", 25,0,25, 25,0,10);

  h_etaphi1 = new TH2D("h_etaphi1","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi2 = new TH2D("h_etaphi2","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi3 = new TH2D("h_etaphi3","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi4 = new TH2D("h_etaphi4","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi5 = new TH2D("h_etaphi5","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi6 = new TH2D("h_etaphi6","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi7 = new TH2D("h_etaphi7","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi8 = new TH2D("h_etaphi8","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi9 = new TH2D("h_etaphi9","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi10 = new TH2D("h_etaphi10","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi11 = new TH2D("h_etaphi11","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi12 = new TH2D("h_etaphi12","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi13 = new TH2D("h_etaphi13","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi14 = new TH2D("h_etaphi14","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi15 = new TH2D("h_etaphi15","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi16 = new TH2D("h_etaphi16","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi17 = new TH2D("h_etaphi17","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi18 = new TH2D("h_etaphi18","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  m_histograms.push_back(h_etaphi1);
  m_histograms.push_back(h_etaphi2);
  m_histograms.push_back(h_etaphi3);
  m_histograms.push_back(h_etaphi4);
  m_histograms.push_back(h_etaphi5);
  m_histograms.push_back(h_etaphi6);
  m_histograms.push_back(h_etaphi7);
  m_histograms.push_back(h_etaphi8);
  m_histograms.push_back(h_etaphi9);
  m_histograms.push_back(h_etaphi10);
  m_histograms.push_back(h_etaphi11);
  m_histograms.push_back(h_etaphi12);
  m_histograms.push_back(h_etaphi13);
  m_histograms.push_back(h_etaphi14);
  m_histograms.push_back(h_etaphi15);
  m_histograms.push_back(h_etaphi16);
  m_histograms.push_back(h_etaphi17);
  m_histograms.push_back(h_etaphi18);
 
  h_ene_eta = new TH1D("h_ene_eta","", int(4./0.025), -2.0,2.0);
  m_histograms.push_back(h_ene_eta);

  h_ene_phi = new TH1D("h_ene_phi","", 128, -3.1416 ,3.1416);
  m_histograms.push_back(h_ene_phi);

  h_ene_x = new TH1D("h_ene_x","", 10000, -5000.,5000.);
  m_histograms.push_back(h_ene_x);

  h_ene_y = new TH1D("h_ene_y","", 10000, -5000.,5000.);
  m_histograms.push_back(h_ene_y);

  h_ene_z = new TH1D("h_ene_z","", 510, -4600.,4600.);
  m_histograms.push_back(h_ene_z);

  h_ene_r = new TH1D("h_ene_r","", 20, -0.5, 19.5);
  m_histograms.push_back(h_ene_r);

  h_lambdaEcal = new TH1D("h_lambdaEcal","", 13, 0, 13);
  m_histograms.push_back(h_lambdaEcal);
  h_lambdaHcal = new TH1D("h_lambdaHcal","", 13, 0, 13);
  m_histograms.push_back(h_lambdaHcal);
}

void CombinedCellAnalysisPbSpacer::processEvent(podio::EventStore& aStore, int aEventId, bool verbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::PositionedCaloHitCollection*     colECalPositions(nullptr);
  const fcc::PositionedCaloHitCollection*     colHCalPositions(nullptr);


  bool colMCParticlesOK       = aStore.get("GenParticles", colMCParticles);
  // for Pb spacer 
  bool colECalPositionsOK     = aStore.get("ECalBarrelCellPositions" , colECalPositions);
  bool colHCalPositionsOK     = aStore.get("HCalBarrelCellPositions" , colHCalPositions);

  double etaVertex = 0;
  double phiVertex = 0;
  double thetaVertex = 0;
  //Direction of gen. particle
  TVector3 directionParticle(0.,0.,0.);

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
      //unit vector
      directionParticle = particle.Unit();
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
  double cryoE = 0.;

  //Cells in ECal collection
  if (colECalPositionsOK && colHCalPositionsOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #newCaloPositions:     " << colECalPositions->size()    << std::endl;;
      std::cout << " -> #newCaloPositions:     " << colHCalPositions->size()    << std::endl;;
    }
    for (auto& iecl=colECalPositions->begin(); iecl!=colECalPositions->end(); ++iecl){
      SumE_ecal += iecl->core().energy*m_a;
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

      h_ene_x->Fill(iecl->position().x,iecl->core().energy*m_a);
      h_ene_y->Fill(iecl->position().y,iecl->core().energy*m_a);
      h_ene_z->Fill(iecl->position().z,iecl->core().energy*m_a);
      
      zE_ecal += iecl->position().z*iecl->core().energy*m_a;
      phiE_ecal += phi*iecl->core().energy*m_a;
      etaE_ecal += eta*iecl->core().energy*m_a;
      thetaE_ecal += (2. * atan( exp(-eta) ))*iecl->core().energy*m_a;
      
      h_ene_r           ->Fill(layerId, iecl->core().energy*m_a);
      h_ene_phi         ->Fill(phi, iecl->core().energy*m_a);
      h_ene_eta         ->Fill(eta, iecl->core().energy*m_a);
      if ((layerId +1) == 1)
	h_etaphi1->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*m_a);
      else if ((layerId +1) == 2)
	h_etaphi2->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*m_a);
      else if ((layerId +1) == 3)
	h_etaphi3->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*m_a);
      else if ((layerId +1) == 4)
	h_etaphi4->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*m_a);
      else if ((layerId +1) == 5)
	h_etaphi5->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*m_a);
      else if ((layerId +1) == 6)
	h_etaphi6->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*m_a);
      else if ((layerId +1) == 7)
	h_etaphi7->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*m_a);
      else if ((layerId +1) == 8)
	h_etaphi8->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*m_a);

      TVector3 showerStart = directionParticle*(rMinEcal/float(sin(thetaVertex)));
      TVector3 hitVector = hitPosition-showerStart;
      double hitLong = directionParticle*hitVector;
      // multiplied by *m_a, correction from benchmark
      h_lambdaEcal->Fill((lambdaOffset/float(cos(directionParticle.Eta())) + hitLong/lambdaEcal), iecl->core().energy*m_a);
    }
    for (auto& iecl=colHCalPositions->begin(); iecl!=colHCalPositions->end(); ++iecl){
      if (iecl->core().energy > m_thr){
	SumE_hcal += iecl->core().energy*sf_hcalCorr;
	TVector3 hitPosition(iecl->position().x,iecl->position().y,iecl->position().z);
	double r = sqrt(pow(iecl->position().x,2)+pow(iecl->position().y,2));
	TVector3 vec(iecl->position().x,iecl->position().y,iecl->position().z);
	double phi = atan2( iecl->position().y, iecl->position().x );
	double eta = vec.Eta();
	int layerId = hcal_decoder.value("layer",iecl->core().cellId);

	// Add energy in first HCal layer
	if ( layerId == 0 ){	
	  E_firstLayer += iecl->core().energy*sf_hcalCorr;
	}
	if (verbose){
	  std::cout << "HCAL :  \n";
	  std::cout << " r " << r << ", eta " << eta << ", phi " << phi << ", energy " << iecl->core().energy << "\n";
	  std::cout << " layer " << layerId << std::endl;
	}
	h_ene_x->Fill(iecl->position().x,iecl->core().energy*sf_hcalCorr);
	h_ene_y->Fill(iecl->position().y,iecl->core().energy*sf_hcalCorr);
	h_ene_z->Fill(iecl->position().z,iecl->core().energy*sf_hcalCorr);

	zE_hcal += iecl->position().z*iecl->core().energy*sf_hcalCorr;
	phiE_hcal += phi*iecl->core().energy*sf_hcalCorr;
	etaE_hcal += eta*iecl->core().energy*sf_hcalCorr;
	thetaE_hcal += (2. * atan( exp(-eta) ))*iecl->core().energy*sf_hcalCorr;
	
	h_ene_r           ->Fill((layerId+8), iecl->core().energy*sf_hcalCorr);
	h_ene_phi         ->Fill(phi, iecl->core().energy*sf_hcalCorr);
	h_ene_eta         ->Fill(eta, iecl->core().energy*sf_hcalCorr);
	if ((layerId +9) == 9)
	  h_etaphi9->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*sf_hcalCorr);
	else if ((layerId +9) == 10)
	  h_etaphi10->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*sf_hcalCorr);
	else if ((layerId +9) == 11)
	  h_etaphi11->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*sf_hcalCorr);
	else if ((layerId +9) == 12)
	  h_etaphi12->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*sf_hcalCorr);
	else if ((layerId +9) == 13)
	  h_etaphi13->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*sf_hcalCorr);
	else if ((layerId +9) == 14)
	  h_etaphi14->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*sf_hcalCorr);
	else if ((layerId +9) == 15)
	  h_etaphi15->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*sf_hcalCorr);
	else if ((layerId +9) == 16)
	  h_etaphi16->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*sf_hcalCorr);
	else if ((layerId +9) == 17)
	  h_etaphi17->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*sf_hcalCorr);
	else if ((layerId +9) == 18)
	  h_etaphi18->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy*sf_hcalCorr);
      
	TVector3 showerStart = directionParticle*(rMinHcal/float(sin(thetaVertex)));
	TVector3 hitVector = hitPosition-showerStart;
	double hitLong = directionParticle*hitVector;
	// multiplied by *m_b, correction from benchmark
	h_lambdaHcal->Fill((lambdaOffsetHcal/float(cos(directionParticle.Eta())) + hitLong/lambdaHcal), iecl->core().energy * m_b);
      }
    }
   
    if (verbose) {
      std::cout << "Total cell energy:                 " << SumE_ecal + SumE_hcal << std::endl;
    }


    //Fill histograms
    h_cellEnergy      ->Fill(SumE_ecal+SumE_hcal);
    h_cellEnergy_ecal ->Fill(SumE_ecal);
    h_cellEnergy_hcal ->Fill(SumE_hcal);
    h_cellEnergy_ecalP->Fill(SumE_ecal*m_a);
    h_cellEnergy_hcalP->Fill(SumE_hcal);
    h_phiRec_ecal     ->Fill(phiE_ecal/SumE_ecal - phiVertex);
    h_etaRec_ecal     ->Fill(etaE_ecal/SumE_ecal - etaVertex);
    h_thetaRec_ecal   ->Fill(thetaE_ecal/SumE_ecal - thetaVertex);
    h_phiRec_hcal     ->Fill(phiE_hcal/SumE_hcal - phiVertex);
    h_etaRec_hcal     ->Fill(etaE_hcal/SumE_hcal - etaVertex);
    h_thetaRec_hcal   ->Fill(thetaE_hcal/SumE_hcal - thetaVertex);
    h_phiRec          ->Fill((phiE_hcal+phiE_ecal)/(SumE_hcal+SumE_ecal) - phiVertex);
    h_etaRec          ->Fill((etaE_hcal+etaE_ecal)/(SumE_hcal+SumE_ecal) - etaVertex);
    h_thetaRec        ->Fill((thetaE_hcal+thetaE_ecal)/(SumE_hcal+SumE_ecal) - thetaVertex);
    double barycenterR_ecal = (zE_ecal/SumE_ecal) * tan((thetaE_ecal/SumE_ecal));
    double barycenterR_hcal = (zE_hcal/SumE_hcal) * tan((thetaE_hcal/SumE_hcal));

    // Benchmark reconstruction
    double E0_rec = SumE_ecal*m_a + SumE_hcal + m_b*sqrt(abs(E_firstLayer*m_a*E_lastLayer)) + m_c*pow(SumE_ecal*m_a,2);
    //    std::cout << " Energy benchmark : " << E0_rec << "\n";
    h_benchmark        ->Fill( E0_rec );
    h_benchmark2nd     ->Fill( pow((E0_rec - m_aLin)/m_bLin, 1/m_cLin) );
    h_cellEnergy_first ->Fill( E_firstLayer );
    h_cellEnergy_last  ->Fill( E_lastLayer );
    h_lostECorr        ->Fill((E_firstLayer+E_lastLayer), cryoE);
    h_lostEMultiCorr        ->Fill(sqrt(abs(E_firstLayer*E_lastLayer)), cryoE);
  }
  else {
    std::cout << "No colECalPositions Collection!!!!!" << std::endl;
    std::cout << "No colHCalPositions Collection!!!!!" << std::endl;
  } 
}
void CombinedCellAnalysisPbSpacer::finishLoop(int aNumEvents, bool aVerbose) {
  std::cout << "Total energy: " << h_cellEnergy->GetMean() << ", above threshold " << m_thr << std::endl;
  aNumEvents = h_ene_r->GetEntries();
  
  //h_ene_r  ->Scale(1./h_ene_r->Integral());
  h_ene_phi->Scale(1./h_ene_phi->Integral());
  h_ene_eta->Scale(1./h_ene_eta->Integral());
  
  std::cout << "Integral phi " << h_ene_phi->Integral() << std::endl;
}
