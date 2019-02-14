#include "CombinedCellAnalysisFwd.h"

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

CombinedCellAnalysisFwd::CombinedCellAnalysisFwd(double aEnergy, double aThr, double aA, double aB, double aC, double aD) : m_energy(aEnergy), ecal_decoder(""), hcal_decoder(""), m_thr(aThr), m_a(aA), m_b(aB), m_c(aC), m_d(aD), SumE_ecal(0), SumE_hcal(0), h_cellEnergy(nullptr) {
  Initialize_histos();
}
CombinedCellAnalysisFwd::CombinedCellAnalysisFwd(double aEnergy, double aThr, double aA, double aB, double aC, double aD, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal, const std::string& aEcalCollection, const std::string& aHcalCollection) : m_energy(aEnergy), m_thr(aThr), m_a(aA), m_b(aB), m_c(aC), m_d(aD), SumE_ecal(0), SumE_hcal(0), h_cellEnergy(nullptr), ecal_decoder(aBitfieldEcal), hcal_decoder(aBitfieldHcal), ecal_collection(aEcalCollection), hcal_collection(aHcalCollection) {
  Initialize_histos();
}

CombinedCellAnalysisFwd::~CombinedCellAnalysisFwd() {}


void CombinedCellAnalysisFwd::Initialize_histos() {
  h_benchmark = new TH1F("h_benchmark","",1500,0,15000);
  m_histograms.push_back(h_benchmark);

  h_cellEnergy = new TH1F("h_cellEnergy","",1500,0,15000);
  m_histograms.push_back(h_cellEnergy);

  double Et = m_energy / cosh(3);

  h_ET = new TH1F("h_ET","",100,Et-Et/2,Et+Et/2);
  m_histograms.push_back(h_ET);

  h_cellEnergy_ecal = new TH1F("h_cellEnergy_ecal","",1000,0,m_energy/2);
  m_histograms.push_back(h_cellEnergy_ecal);

  h_cellEnergy_hcal = new TH1F("h_cellEnergy_hcal","",1000,0,m_energy/2);
  m_histograms.push_back(h_cellEnergy_hcal);

  h_sf_ecal = new TH1F("h_sf_ecal","",100,0,.01);
  m_histograms.push_back(h_sf_ecal);

  h_sf_hcal = new TH1F("h_sf_hcal","",100,0,.01);
  m_histograms.push_back(h_sf_hcal);

  h_ET_ecal = new TH1F("h_ET_ecal","",100,0,Et);
  m_histograms.push_back(h_ET_ecal);

  h_ET_hcal = new TH1F("h_ET_hcal","",100,0,Et);
  m_histograms.push_back(h_ET_hcal);

  h_cellEnergy_ecalP = new TH1F("h_cellEnergy_ecalP","",100,0,m_energy);
  m_histograms.push_back(h_cellEnergy_ecalP);

  h_cellEnergy_hcalP = new TH1F("h_cellEnergy_hcalP","",100,0,m_energy);
  m_histograms.push_back(h_cellEnergy_hcalP);

  h_phiRec_ecal = new TH1F("h_phiRec_ecal","",1000,-.4,.4);
  m_histograms.push_back(h_phiRec_ecal);

  h_phiRec_hcal = new TH1F("h_phiRec_hcal","",1000,-.4,.4);
  m_histograms.push_back(h_phiRec_hcal);

  h_phiRec = new TH1F("h_phiRec","",1000,-.4,.4);
  m_histograms.push_back(h_phiRec);

  h_etaRec = new TH1F("h_etaRec","",1000,-.1,.1);
  m_histograms.push_back(h_etaRec);

  h_etaRec_ecal = new TH1F("h_etaRec_ecal","",1000,-.1,.1);
  m_histograms.push_back(h_etaRec_ecal);

  h_etaRec_hcal = new TH1F("h_etaRec_hcal","",1000,-.1,.1);
  m_histograms.push_back(h_etaRec_hcal);

  h_thetaRec_ecal = new TH1F("h_thetaRec_ecal","",1000,-.1,.1);
  m_histograms.push_back(h_thetaRec_ecal);

  h_thetaRec_hcal = new TH1F("h_thetaRec_hcal","",1000,-.1,.1);
  m_histograms.push_back(h_thetaRec_hcal);

  h_thetaRec = new TH1F("h_thetaRec","",1000,-.1,.1);
  m_histograms.push_back(h_thetaRec);

  h_cellEnergy_first = new TH1F("h_cellEnergy_first","Energy in first HCAL layer",100,0,m_energy/2);
  m_histograms.push_back(h_cellEnergy_first);

  h_cellEnergy_last = new TH1F("h_cellEnergy_last","Energy in last ECAL layer",100,0,m_energy/2);
  m_histograms.push_back(h_cellEnergy_last);

  h_cellId = new TH1F("h_cellId","", 1000, 0,5000e6);
  m_histograms.push_back(h_cellId);

  h_lostECorr = new TH2F("h_lostECorr","", 50,0,50, 50,0,50);

  h_etaphi1 = new TH2F("h_etaphi1","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi2 = new TH2F("h_etaphi2","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi3 = new TH2F("h_etaphi3","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi4 = new TH2F("h_etaphi4","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi5 = new TH2F("h_etaphi5","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi6 = new TH2F("h_etaphi6","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi7 = new TH2F("h_etaphi7","", 160,-4,4, 128, -TMath::Pi(), TMath::Pi());
  h_etaphi8 = new TH2F("h_etaphi8","", 160,-4,4, 128, -TMath::Pi(), TMath::Pi());
  h_etaphi9 = new TH2F("h_etaphi9","", 160,-4,4, 128, -TMath::Pi(), TMath::Pi());
  h_etaphi10 = new TH2F("h_etaphi10","", 160,-4,4, 128, -TMath::Pi(), TMath::Pi());
  h_etaphi11 = new TH2F("h_etaphi11","", 160,-4,4, 128, -TMath::Pi(), TMath::Pi());
  h_etaphi12 = new TH2F("h_etaphi12","", 160,-4,4, 128, -TMath::Pi(), TMath::Pi());
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
 
  h_ene_eta = new TH1F("h_ene_eta","", 160,-4,4);
  m_histograms.push_back(h_ene_eta);

  h_ene_phi = new TH1F("h_ene_phi","", 128, -TMath::Pi(), TMath::Pi());
  m_histograms.push_back(h_ene_phi);

  h_ene_x = new TH1F("h_ene_x","", 10000, -5000.,5000.);
  m_histograms.push_back(h_ene_x);

  h_ene_y = new TH1F("h_ene_y","", 10000, -5000.,5000.);
  m_histograms.push_back(h_ene_y);

  h_ene_z = new TH1F("h_ene_z","", 1000, -20000., 20000.);
  m_histograms.push_back(h_ene_z);

  h_ene_r = new TH1F("h_ene_r","", 95, -0.5, 94.5);
  m_histograms.push_back(h_ene_r);

  h_lambda = new TH1F("h_lambda","", 17, 0, 17);
  m_histograms.push_back(h_lambda);
}

void CombinedCellAnalysisFwd::processEvent(podio::EventStore& aStore, int aEventId, bool verbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::PositionedCaloHitCollection*     colECalPositions(nullptr);
  const fcc::PositionedCaloHitCollection*     colHCalPositions(nullptr);


  bool colMCParticlesOK       = aStore.get("GenParticles", colMCParticles);
  // for Pb spacer 
  bool colECalPositionsOK     = aStore.get(ecal_collection , colECalPositions);
  bool colHCalPositionsOK     = aStore.get(hcal_collection , colHCalPositions);

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
      thetaVertex = particle.Theta(); //2 * atan( exp(-etaVertex) );
      //unit vector
      directionParticle = particle.Unit();
      //Fill histograms
      //      h_ptGen->Fill( sqrt( pow(iparticle->core().p4.px,2)+pow(iparticle->core().p4.py,2) ) );
    }
  }
  else {
    std::cout << "No MCTruth info available" << std::endl;
  }

  //Total hit energy per event
  SumE_ecal = 0.;
  SumE_hcal = 0.;
  SumET_ecal = 0.;
  SumET_hcal = 0.;
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
  double absorberE_ecal = 0.;
  double absorberE_hcal = 0.;

  //Cells in ECal collection
  if (colECalPositionsOK && colHCalPositionsOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #newCaloPositions:     " << colECalPositions->size()    << std::endl;;
      std::cout << " -> #newCaloPositions:     " << colHCalPositions->size()    << std::endl;;
    }
    for (auto& iecl=colECalPositions->begin(); iecl!=colECalPositions->end(); ++iecl){
      TVector3 hitPosition(iecl->position().x,iecl->position().y,iecl->position().z);
      double r = hitPosition.Perp();
      double phi = atan2( iecl->position().y, iecl->position().x );
      double eta = hitPosition.Eta();
      double ET =  iecl->core().energy / cosh(eta);
      int layerId = ecal_decoder.value("layer",iecl->core().cellId);
      int cryoId = ecal_decoder.value("cryo",iecl->core().cellId);
      int typeId = ecal_decoder.value("type",iecl->core().cellId);
      if (cryoId==1){
	// Fill histogram for correlation plot
	// Re-calibrate the enegery in the cryostat.. due to layerId0 -> calibrated to EM scale of first ecal layer
	cryoE += iecl->core().energy;
	continue;
      }
      // Energy in absorber+readout (passive material)
      if (typeId > 0) {
	absorberE_ecal += iecl->core().energy; 
      }
      // Energy in sensitive
      else {
	SumE_ecal += iecl->core().energy;
	SumET_ecal += ET;
	// Add energy in FIRST ECal layer
	if ( layerId > -1 && layerId < 7 ){	
	  E_lastLayer += iecl->core().energy;
	}
	if (verbose){
	  std::cout << " ECAL :  \n";
	  std::cout << " r " << r << ", eta " << eta << ", phi " << phi << ", energy " << iecl->core().energy << "\n";
	  std::cout << " layer " << layerId << std::endl;
	}
      
	h_ene_x->Fill(iecl->position().x,iecl->core().energy);
	h_ene_y->Fill(iecl->position().y,iecl->core().energy);
	h_ene_z->Fill(iecl->position().z,iecl->core().energy);
	
	zE_ecal += iecl->position().z*iecl->core().energy;
	phiE_ecal += phi*iecl->core().energy;
	etaE_ecal += eta*iecl->core().energy;
	thetaE_ecal += (2. * atan( exp(-eta) ))*iecl->core().energy;
      
	h_ene_r           ->Fill(r,iecl->core().energy);
	h_ene_phi         ->Fill(phi,iecl->core().energy);
	h_ene_eta         ->Fill(eta,iecl->core().energy);
	if (layerId > -1 && layerId < 7)
	  h_etaphi1->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy);
	else if (layerId > 6 && layerId < 14)
	  h_etaphi2->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy);
	else if (layerId > 13 && layerId < 21)
	  h_etaphi3->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy);
	else if (layerId > 20 && layerId < 28)
	  h_etaphi4->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy);
	else if (layerId > 27 && layerId < 35)
	  h_etaphi5->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy);
	else if (layerId > 34)
	  h_etaphi6->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy);
	
	TVector3 showerStart = directionParticle*(zMinEcal/float(cos(thetaVertex)));
	TVector3 hitVector = hitPosition-showerStart;
	double hitLong = directionParticle*hitVector;
	h_lambda->Fill((lambdaOffset + hitLong/lambdaEcal), iecl->core().energy);
      }
    }
    for (auto& iecl=colHCalPositions->begin(); iecl!=colHCalPositions->end(); ++iecl){
      if (iecl->core().energy > m_thr){
	TVector3 hitPosition(iecl->position().x,iecl->position().y,iecl->position().z);
	double r = hitPosition.Perp();
	TVector3 vec(iecl->position().x,iecl->position().y,iecl->position().z);
	double phi = atan2( iecl->position().y, iecl->position().x );
	double eta = vec.Eta();
	double ET =  iecl->core().energy / cosh(eta);
	int layerId = hcal_decoder.value("layer",iecl->core().cellId);
	int typeId = hcal_decoder.value("type",iecl->core().cellId);
	// Energy in absorber+readout (passive material)
	if (typeId > 0) {
	  absorberE_hcal += iecl->core().energy; 
	}
	// Energy in sensitive
	else {
	  SumE_hcal += iecl->core().energy * hcalCorr_sf;
	  SumET_hcal += ET;
	  // Add energy in first HCal layer
	  if ( layerId > -1 && layerId < 8 ){	
	    E_firstLayer += iecl->core().energy * hcalCorr_sf;
	  }
	  if (verbose){
	    std::cout << "HCAL :  \n";
	    std::cout << " r " << r << ", eta " << eta << ", phi " << phi << ", energy " << iecl->core().energy * hcalCorr_sf << "\n";
	    std::cout << " layer " << layerId << std::endl;
	  }
	  h_ene_x->Fill(iecl->position().x,iecl->core().energy* hcalCorr_sf);
	  h_ene_y->Fill(iecl->position().y,iecl->core().energy* hcalCorr_sf);
	  h_ene_z->Fill(iecl->position().z,iecl->core().energy* hcalCorr_sf);

	  zE_hcal += iecl->position().z*iecl->core().energy* hcalCorr_sf;
	  phiE_hcal += phi*iecl->core().energy* hcalCorr_sf;
	  etaE_hcal += eta*iecl->core().energy* hcalCorr_sf;
	  thetaE_hcal += (2. * atan( exp(-eta) ))*iecl->core().energy* hcalCorr_sf;
	
	  h_ene_r           ->Fill(r,iecl->core().energy* hcalCorr_sf);
	  h_ene_phi         ->Fill(phi,iecl->core().energy* hcalCorr_sf);
	  h_ene_eta         ->Fill(eta,iecl->core().energy)* hcalCorr_sf;
	  if (layerId > -1 && layerId < 8)
	    h_etaphi7->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy* hcalCorr_sf);
	  else if (layerId > 7 && layerId < 17)
	    h_etaphi8->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy* hcalCorr_sf);
	  else if (layerId > 16 && layerId < 26)
	    h_etaphi9->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy* hcalCorr_sf);
	  else if (layerId > 25 && layerId < 35)
	    h_etaphi10->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy* hcalCorr_sf);
	  else if (layerId > 34 && layerId < 44)
	    h_etaphi11->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy* hcalCorr_sf);
	  else if (layerId > 43 && layerId < 53)
	    h_etaphi12->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy* hcalCorr_sf);
	  
	  TVector3 showerStart = directionParticle*(zMinHcal/float(cos(thetaVertex)));
	  TVector3 hitVector = hitPosition-showerStart;
	  double hitLong = directionParticle*hitVector;
	  h_lambda->Fill((lambdaOffsetHcal + hitLong/lambdaHcal), iecl->core().energy* hcalCorr_sf);
	}
      }
    }
   
    if (verbose) {
      std::cout << "Total cell energy:                 " << SumE_ecal + SumE_hcal << std::endl;
      std::cout << "Total cell energy, diff fsampling: " << SumE_ecal + SumE_hcal << std::endl;
    }


    //Fill histograms
    h_cellEnergy      ->Fill(SumE_ecal+SumE_hcal);
    h_cellEnergy_ecal ->Fill(SumE_ecal);
    h_cellEnergy_hcal ->Fill(SumE_hcal);
    h_ET              ->Fill(SumET_ecal+SumET_hcal);
    h_ET_ecal          ->Fill(SumET_ecal);
    h_ET_hcal          ->Fill(SumET_hcal);
    h_cellEnergy_ecalP->Fill(SumE_ecal*m_a);
    h_cellEnergy_hcalP->Fill(SumE_hcal*m_b);
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
    //    h_thetaRec        ->Fill( atan( (barycenterR_hcal-barycenterR_ecal)/(zE_hcal/SumE_hcal - zE_ecal/SumE_ecal) ) - thetaVertex);

    // Benchmark reconstruction
    double E0_rec = SumE_ecal*m_a + SumE_hcal + m_c*sqrt(abs(E_firstLayer*m_a*E_lastLayer)) + m_d*pow(SumE_ecal*m_a,2);
    h_benchmark        ->Fill( E0_rec );
    h_sf_hcal -> Fill( SumE_hcal/float(SumE_hcal+absorberE_hcal) );
    h_sf_ecal -> Fill( SumE_ecal/float(SumE_ecal+absorberE_ecal) );

    h_cellEnergy_first ->Fill( E_firstLayer );
    h_cellEnergy_last  ->Fill( E_lastLayer );
    h_lostECorr        ->Fill((E_firstLayer+E_lastLayer), cryoE);
  }
  else {
    std::cout << "No colECalPositions Collection!!!!!" << std::endl;
    std::cout << "No colHCalPositions Collection!!!!!" << std::endl;
  } 
}
void CombinedCellAnalysisFwd::finishLoop(int aNumEvents, bool aVerbose) {
  std::cout << "Total energy: " << h_cellEnergy->GetMean() << ", above threshold " << m_thr << std::endl;
  aNumEvents = h_ene_r->GetEntries();
  
  // h_lambda->Scale(1./(double)aNumEvents);
  h_ene_r  ->Scale(1./h_ene_r->Integral());
  h_ene_phi->Scale(1./h_ene_phi->Integral());
  h_ene_eta->Scale(1./h_ene_eta->Integral());
  
  std::cout << "Integral phi " << h_ene_phi->Integral() << std::endl;
}
