#include "CombinedCellAnalysisEndcap.h"

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

CombinedCellAnalysisEndcap::CombinedCellAnalysisEndcap(double aEnergy, double aThr, double aA, double aB, double aC, double aD) : m_energy(aEnergy), ecal_decoder(""), hcal_decoder(""), hcalEB_decoder(""), m_thr(aThr), m_a(aA), m_b(aB), m_c(aC), m_d(aD), SumE_ecal(0), SumE_hcal(0), h_cellEnergy(nullptr) {
  Initialize_histos();
}
CombinedCellAnalysisEndcap::CombinedCellAnalysisEndcap(double aEnergy, double aThr, double aA, double aB, double aC, double aD, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal, const std::string& aBitfieldHcalEB, const std::string& aEcalCollection, const std::string& aHcalCollection, const std::string& aHcalEBCollection) : m_energy(aEnergy), m_thr(aThr), m_a(aA), m_b(aB), m_c(aC), m_d(aD), SumE_ecal(0), SumE_hcal(0), h_cellEnergy(nullptr), ecal_decoder(aBitfieldEcal), hcal_decoder(aBitfieldHcal), hcalEB_decoder(aBitfieldHcalEB), ecal_collection(aEcalCollection), hcal_collection(aHcalCollection), hcalEB_collection(aHcalEBCollection) {
  Initialize_histos();
}

CombinedCellAnalysisEndcap::~CombinedCellAnalysisEndcap() {}


void CombinedCellAnalysisEndcap::Initialize_histos() {
  double sigma = sqrt(pow(0.8*sqrt(m_energy)+0.02*(m_energy),2));
  h_benchmark = new TH1D("h_benchmark","", 2*std::ceil(4*sigma)+std::ceil(3*sigma), std::floor(m_energy-4*sigma), std::ceil(m_energy+3*sigma));
  m_histograms.push_back(h_benchmark);

  h_cellEnergy = new TH1D("h_cellEnergy","", 2*std::ceil(4*sigma)+std::ceil(3*sigma), std::floor(m_energy-4*sigma), std::ceil(m_energy+3*sigma));
  m_histograms.push_back(h_cellEnergy);

  double Et = m_energy / cosh(2);

  h_ET = new TH1D("h_ET","",100,Et-Et/2,Et+Et/2);
  m_histograms.push_back(h_ET);

  h_cellEnergy_ecal = new TH1D("h_cellEnergy_ecal","",1000,0,m_energy/2);
  m_histograms.push_back(h_cellEnergy_ecal);

  h_cellEnergy_hcal = new TH1D("h_cellEnergy_hcal","",1000,0,m_energy/2);
  m_histograms.push_back(h_cellEnergy_hcal);

  h_cellEnergy_hcalEB = new TH1D("h_cellEnergy_hcalEB","",1000,0,m_energy/2);
  m_histograms.push_back(h_cellEnergy_hcalEB);

  h_sf_ecal = new TH1D("h_sf_ecal","",100,0,.1);
  m_histograms.push_back(h_sf_ecal);

  h_sf_hcal = new TH1D("h_sf_hcal","",100,0,.1);
  m_histograms.push_back(h_sf_hcal);

  h_ET_ecal = new TH1D("h_ET_ecal","",100,0,Et);
  m_histograms.push_back(h_ET_ecal);

  h_ET_hcal = new TH1D("h_ET_hcal","",100,0,Et);
  m_histograms.push_back(h_ET_hcal);

  h_ET_hcalEB = new TH1D("h_ET_hcalEB","",100,0,Et);
  m_histograms.push_back(h_ET_hcalEB);

  h_cellEnergy_ecalP = new TH1D("h_cellEnergy_ecalP","",1000,0,m_energy);
  m_histograms.push_back(h_cellEnergy_ecalP);

  h_cellEnergy_hcalP = new TH1D("h_cellEnergy_hcalP","",1000,0,m_energy);
  m_histograms.push_back(h_cellEnergy_hcalP);

  h_phiRec_ecal = new TH1D("h_phiRec_ecal","",1000,-.4,.4);
  m_histograms.push_back(h_phiRec_ecal);

  h_phiRec_hcal = new TH1D("h_phiRec_hcal","",1000,-.4,.4);
  m_histograms.push_back(h_phiRec_hcal);

  h_phiRec = new TH1D("h_phiRec","",1000,-.4,.4);
  m_histograms.push_back(h_phiRec);

  h_etaRec = new TH1D("h_etaRec","",1000,-.1,.1);
  m_histograms.push_back(h_etaRec);

  h_etaRec_ecal = new TH1D("h_etaRec_ecal","",1000,-.1,.1);
  m_histograms.push_back(h_etaRec_ecal);

  h_etaRec_hcal = new TH1D("h_etaRec_hcal","",1000,-.1,.1);
  m_histograms.push_back(h_etaRec_hcal);

  h_thetaRec_ecal = new TH1D("h_thetaRec_ecal","",1000,-.1,.1);
  m_histograms.push_back(h_thetaRec_ecal);

  h_thetaRec_hcal = new TH1D("h_thetaRec_hcal","",1000,-.1,.1);
  m_histograms.push_back(h_thetaRec_hcal);

  h_thetaRec = new TH1D("h_thetaRec","",1000,-.1,.1);
  m_histograms.push_back(h_thetaRec);

  h_cellEnergy_first = new TH1D("h_cellEnergy_first","Energy in first HCAL layer",100,0,m_energy/2);
  m_histograms.push_back(h_cellEnergy_first);

  h_cellEnergy_last = new TH1D("h_cellEnergy_last","Energy in last ECAL layer",100,0,m_energy/2);
  m_histograms.push_back(h_cellEnergy_last);

  h_cellId = new TH1D("h_cellId","", 1000, 0,5000e6);
  m_histograms.push_back(h_cellId);

  h_lostECorr = new TH2D("h_lostECorr","", 50,0,50, 50,0,50);

  h_etaphi1 = new TH2D("h_etaphi1","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi2 = new TH2D("h_etaphi2","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi3 = new TH2D("h_etaphi3","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi4 = new TH2D("h_etaphi4","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi5 = new TH2D("h_etaphi5","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi6 = new TH2D("h_etaphi6","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi7 = new TH2D("h_etaphi7","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi8 = new TH2D("h_etaphi8","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
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
  h_etaphi19 = new TH2D("h_etaphi19","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi20 = new TH2D("h_etaphi20","", 320,-4,4, 256, -TMath::Pi(), TMath::Pi());
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
  m_histograms.push_back(h_etaphi19);
  m_histograms.push_back(h_etaphi20);
 
  h_ene_eta = new TH1D("h_ene_eta","", 320, -4.0,4.0);
  m_histograms.push_back(h_ene_eta);

  h_ene_phi = new TH1D("h_ene_phi","", 256, -TMath::Pi(), TMath::Pi());
  m_histograms.push_back(h_ene_phi);

  h_ene_x = new TH1D("h_ene_x","", 10000, -5000.,5000.);
  m_histograms.push_back(h_ene_x);

  h_ene_y = new TH1D("h_ene_y","", 10000, -5000.,5000.);
  m_histograms.push_back(h_ene_y);

  h_ene_z = new TH1D("h_ene_z","", 900, -9000.,9000.);
  m_histograms.push_back(h_ene_z);

  h_ene_r = new TH1D("h_ene_r","", 239, -0.5, 238.5);
  m_histograms.push_back(h_ene_r);

  h_lambdaEcal = new TH1D("h_lambdaEcal","", 13, 0, 13);
  m_histograms.push_back(h_lambdaEcal);
  h_lambdaHcal = new TH1D("h_lambdaHcal","", 13, 0, 13);
  m_histograms.push_back(h_lambdaHcal);
  h_lambdaHcalEB = new TH1D("h_lambdaHcalEB","", 13, 0, 13);
  m_histograms.push_back(h_lambdaHcalEB);
}

void CombinedCellAnalysisEndcap::processEvent(podio::EventStore& aStore, int aEventId, bool verbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::PositionedCaloHitCollection*     colECalPositions(nullptr);
  const fcc::PositionedCaloHitCollection*     colHCalPositions(nullptr);
  const fcc::PositionedCaloHitCollection*     colHCalEBPositions(nullptr);


  bool colMCParticlesOK       = aStore.get("GenParticles", colMCParticles);
  // for Pb spacer 
  bool colECalPositionsOK     = aStore.get(ecal_collection , colECalPositions);
  bool colHCalPositionsOK     = aStore.get(hcal_collection , colHCalPositions);
  bool colHCalEBPositionsOK   = aStore.get(hcalEB_collection , colHCalEBPositions);

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
    }
  }
  else {
    std::cout << "No MCTruth info available" << std::endl;
  }

  //Total hit energy per event
  SumE_ecal = 0.;
  SumE_hcal = 0.;
  SumE_hcalEB = 0.;
  SumET_ecal = 0.;
  SumET_hcal = 0.;
  SumET_hcalEB = 0.;
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
	if ( layerId > -1 && layerId < 26 ){	
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
	
	h_ene_r           ->Fill(layerId,iecl->core().energy);
	h_ene_phi         ->Fill(phi,iecl->core().energy);
	h_ene_eta         ->Fill(eta,iecl->core().energy);
	if (layerId > -1 && layerId < 26)
	  h_etaphi1->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy);
	else if (layerId > 25 && layerId < 52)
	  h_etaphi2->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy);
	else if (layerId > 51 && layerId < 78)
	  h_etaphi3->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy);
	else if (layerId > 77 && layerId < 103)
	  h_etaphi4->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy);
	else if (layerId > 102 && layerId < 129)
	  h_etaphi5->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy);
	else if (layerId > 128)
	  h_etaphi6->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy);
	
	TVector3 showerStart = directionParticle*(zMinEcal/float(cos(thetaVertex)));
	TVector3 hitVector = hitPosition-showerStart;
	double hitLong = directionParticle*hitVector;
	h_lambdaEcal->Fill((lambdaOffset + hitLong/lambdaEcal), iecl->core().energy);
      }
    }
    for (auto& iecl=colHCalPositions->begin(); iecl!=colHCalPositions->end(); ++iecl){
      if (iecl->core().energy > m_thr){
	TVector3 hitPosition(iecl->position().x,iecl->position().y,iecl->position().z);
	double r = hitPosition.Perp();
	double phi = atan2( iecl->position().y, iecl->position().x );
	double eta = hitPosition.Eta();
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
	  if ( layerId > -1 && layerId < 13 ){	
	    E_firstLayer += iecl->core().energy;
	  }
	  if (verbose){
	    std::cout << "HCAL :  \n";
	    std::cout << " r " << r << ", eta " << eta << ", phi " << phi << ", energy " << iecl->core().energy  * hcalCorr_sf<< "\n";
	    std::cout << " layer " << layerId << std::endl;
	  }
	  h_ene_x->Fill(iecl->position().x,iecl->core().energy * hcalCorr_sf);
	  h_ene_y->Fill(iecl->position().y,iecl->core().energy * hcalCorr_sf);
	  h_ene_z->Fill(iecl->position().z,iecl->core().energy * hcalCorr_sf);
	  
	  zE_hcal += iecl->position().z*iecl->core().energy * hcalCorr_sf;
	  phiE_hcal += phi*iecl->core().energy * hcalCorr_sf;
	  etaE_hcal += eta*iecl->core().energy * hcalCorr_sf;
	  thetaE_hcal += (2. * atan( exp(-eta) ))*iecl->core().energy * hcalCorr_sf;
	  
	  h_ene_r           ->Fill((layerId + 156),iecl->core().energy * hcalCorr_sf);
	  h_ene_phi         ->Fill(phi,iecl->core().energy * hcalCorr_sf);
	  h_ene_eta         ->Fill(eta,iecl->core().energy * hcalCorr_sf);
	  if (layerId > -1 && layerId < 13)
	    h_etaphi7->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy * hcalCorr_sf);
	  else if (layerId > 12 && layerId < 27)
	    h_etaphi8->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy * hcalCorr_sf);
	  else if (layerId > 26 && layerId < 41)
	    h_etaphi9->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy * hcalCorr_sf);
	  else if (layerId > 40 && layerId < 55)
	    h_etaphi10->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy * hcalCorr_sf);
	  else if (layerId > 54 && layerId < 69)
	    h_etaphi11->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy * hcalCorr_sf);
	  else if (layerId > 68 && layerId < 83)
	    h_etaphi12->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy * hcalCorr_sf);
	  
	  TVector3 showerStart = directionParticle*(zMinHcal/float(cos(thetaVertex)));
	  TVector3 hitVector = hitPosition-showerStart;
	  double hitLong = directionParticle*hitVector;
	  h_lambdaHcal->Fill((lambdaOffsetHcal + hitLong/lambdaHcal), iecl->core().energy* hcalCorr_sf);
	}
      }
    }
  }
  else {
    std::cout << "No colECalPositions Collection!!!!!" << std::endl;
    std::cout << "No colHCalPositions Collection!!!!!" << std::endl;
  } 
  if (colHCalEBPositionsOK){
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #HCalEBPositions:     " << colHCalEBPositions->size()    << std::endl;;
    }
    for (auto& iecl=colHCalEBPositions->begin(); iecl!=colHCalEBPositions->end(); ++iecl){
      if (iecl->core().energy > m_thr){
	TVector3 hitPosition(iecl->position().x,iecl->position().y,iecl->position().z);
	double r = hitPosition.Perp();
	double phi = atan2( iecl->position().y, iecl->position().x );
	double eta = hitPosition.Eta();
	double ET =  iecl->core().energy / cosh(eta);
	int layerId = hcalEB_decoder.value("layer",iecl->core().cellId);
	int typeId = hcalEB_decoder.value("type",iecl->core().cellId);
	
	// type == 0,2, first module type
	// type == 1,3, second module type
	// fill layer maps for module 2 types
	if (typeId == 1 || typeId == 3) {
	  if (layerId > -1 && layerId < 1)
	    h_etaphi13->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy );
	  else if (layerId > 0 && layerId < 2)
	    h_etaphi14->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy );
	  else if (layerId > 1 && layerId < 3)
	    h_etaphi15->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy );
	  else if (layerId > 2 && layerId < 4)
	    h_etaphi16->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy );
	  else if (layerId > 3 && layerId < 5)
	    h_etaphi17->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy );
	  else if (layerId > 4 && layerId < 6)
	    h_etaphi18->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy );
	  else if (layerId > 5 && layerId < 7)
	    h_etaphi19->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy );
	  else if (layerId > 6 && layerId < 8)
	    h_etaphi20->Fill(eta-etaVertex, phi-phiVertex, iecl->core().energy );
	}
	SumE_hcalEB += iecl->core().energy;
	SumET_hcalEB += ET;
      
	if (iecl->position().z < zMinHcal){
	  TVector3 showerStart = directionParticle*(zMinEcal/float(cos(thetaVertex)));
	  TVector3 hitVector = hitPosition-showerStart;
	  double hitLong = directionParticle*hitVector;
	  h_lambdaHcalEB->Fill((lambdaOffset + hitLong/lambdaEcal), iecl->core().energy);
	}
	else {
	  TVector3 showerStart = directionParticle*(zMinHcal/float(cos(thetaVertex)));
	  TVector3 hitVector = hitPosition-showerStart;
	  double hitLong = directionParticle*hitVector;
	  h_lambdaHcalEB->Fill((lambdaOffsetHcal + hitLong/lambdaHcal), iecl->core().energy* hcalCorr_sf);
	}
      }
    }
  }
  else {
    std::cout << "No colHCalEBPositions Collection!!!!!" << std::endl;
  }

  
  if (verbose) {
    std::cout << "Total cell energy:                 " << SumE_ecal + SumE_hcal + SumE_hcalEB << std::endl;
    std::cout << "Total cell energy, diff fsampling: " << SumE_ecal + SumE_hcal + SumE_hcalEB << std::endl;
  }

    //Fill histograms
    h_cellEnergy        ->Fill(SumE_ecal + SumE_hcal + SumE_hcalEB);
    h_cellEnergy_ecal   ->Fill(SumE_ecal);
    h_cellEnergy_hcal   ->Fill(SumE_hcal);
    h_cellEnergy_hcalEB ->Fill(SumE_hcalEB);
    h_ET                ->Fill(SumET_ecal + SumET_hcal + SumET_hcalEB);
    h_ET_ecal           ->Fill(SumET_ecal);
    h_ET_hcal           ->Fill(SumET_hcal);
    h_ET_hcalEB         ->Fill(SumET_hcalEB);
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
    double E0_rec = SumE_ecal*m_a + SumE_hcal*m_b + SumET_hcalEB + m_c*sqrt(abs(E_firstLayer*m_a)) + m_d*pow(SumE_ecal*m_a,2);
    h_benchmark        ->Fill( E0_rec );
    //    std::cout << " Energy benchmark : " << E0_rec << "\n";
    h_sf_hcal -> Fill( SumE_hcal/float(SumE_hcal+absorberE_hcal) );
    h_sf_ecal -> Fill( SumE_ecal/float(SumE_ecal+absorberE_ecal) );

    h_cellEnergy_first ->Fill( E_firstLayer );
    h_lostECorr        ->Fill( E_firstLayer, absorberE_ecal);

    if (verbose) {
      std::cout << " Sampling fractions :     \n";
      std::cout << " Energy (in passive) ECAL : " << absorberE_ecal << "\n";
      std::cout << " Energy (in passive) HCAL : " << absorberE_hcal << "\n";
      std::cout << " Sampling fractions HCAL:   " << SumE_hcal/float(SumE_hcal+absorberE_hcal) << "\n";
      std::cout << " Sampling fractions ECAL:   " << SumE_ecal/float(SumE_ecal+absorberE_ecal) << std::endl;
    }
    
}
void CombinedCellAnalysisEndcap::finishLoop(int aNumEvents, bool aVerbose) {
  std::cout << "Total energy: " << h_cellEnergy->GetMean() << ", above threshold " << m_thr << std::endl;
  aNumEvents = h_ene_r->GetEntries();
 
  h_ene_r  ->Scale(1./h_ene_r->Integral());
  h_ene_phi->Scale(1./h_ene_phi->Integral());
  h_ene_eta->Scale(1./h_ene_eta->Integral());
  
  std::cout << "Integral phi " << h_ene_phi->Integral() << std::endl;
}
