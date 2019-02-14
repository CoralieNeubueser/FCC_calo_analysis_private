#include "CellAnalysis.h"

// FCC-EDM
#include "datamodel/CaloHitCollection.h"
#include "datamodel/PositionedCaloHitCollection.h"

// PODIO
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

// ROOT
#include "TVector3.h"

// STL
#include <iostream>

CellAnalysis::CellAnalysis(double aEnergy, double aSf, double aThr) : m_energy(aEnergy), m_sf(aSf), m_thr(aThr), SumE_cell(0), h_cellEnergy(nullptr) {
  Initialize_histos();
}

CellAnalysis::~CellAnalysis() {}


void CellAnalysis::Initialize_histos() {
  h_cellEnergy = new TH1F("h_cellEnergy","",1000,0,2);
  m_histograms.push_back(h_cellEnergy);
  h_cellEnergy_check = new TH1F("h_cellEnergy_check","",100000,0,10000);
  m_histograms.push_back(h_cellEnergy_check);

  h_cellId = new TH1F("h_cellId","", 1000, 0,5000e6);
  m_histograms.push_back(h_cellId);

  h_ene_eta = new TH1F("h_ene_eta","", int(4./0.025), -2.0,2.0);
  m_histograms.push_back(h_ene_eta);

  h_ene_phi = new TH1F("h_ene_phi","", 128, -3.1416 ,3.1416);
  m_histograms.push_back(h_ene_phi);

  h_ene_x = new TH1F("h_ene_x","", 10000, -5000.,5000.);
  m_histograms.push_back(h_ene_x);

  h_ene_y = new TH1F("h_ene_y","", 10000, -5000.,5000.);
  m_histograms.push_back(h_ene_y);

  h_ene_z = new TH1F("h_ene_z","", 510, -4600.,4600.);
  m_histograms.push_back(h_ene_z);

  h_ene_r = new TH1F("h_ene_r","", 10, 2850.,4650.);
  m_histograms.push_back(h_ene_r);

  h_ene_eta_check = new TH1F("h_ene_eta_check","", 400, -2.0,2.0);
  m_histograms.push_back(h_ene_eta_check);

  h_ene_phi_check = new TH1F("h_ene_phi_check","", 600, -3.1416 ,3.1416);
  m_histograms.push_back(h_ene_phi_check);

  h_ene_r_check = new TH1F("h_ene_r_check","", 40, 2700.,4700.);
  m_histograms.push_back(h_ene_r_check);
}

void CellAnalysis::processEvent(podio::EventStore& aStore, int aEventId, bool verbose) {
  //Get the collections
  const fcc::CaloHitCollection*               colHCalCell(nullptr);
  const fcc::PositionedCaloHitCollection*     colHCalPositions_new(nullptr);
  const fcc::PositionedCaloHitCollection*     colHCalPositionedHits_old(nullptr);

  bool colHCalCellOK                   = aStore.get("HCalBarrelCells" , colHCalCell);
  bool colHCalPositions_newOK          = aStore.get("HCalPositions" , colHCalPositions_new);
  bool colHCalPositionedHits_oldOK     = aStore.get("HCalPositionedHits" , colHCalPositionedHits_old);

  //Total hit energy per event
  SumE_cell = 0.;

  int cellsBelowThreshold=0;
  double lostEnergy=0;
  //Cell collection
  if (colHCalCellOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #caloCells:     " << colHCalCell->size()    << std::endl;;
    }
    //Loop through the collection
    for (auto& iecl=colHCalCell->begin(); iecl!=colHCalCell->end(); ++iecl)
      {
	//if (verbose) std::cout << "HCal cell energy " << iecl->core().energy << std::endl;
	if (iecl->core().energy > m_thr){
	  SumE_cell += iecl->core().energy;
	  h_cellId->Fill(iecl->core().cellId);
	  
 	}
	else{
	  cellsBelowThreshold++;
	  lostEnergy += iecl->core().energy;
	}
      }
    //Fill histograms
    h_cellEnergy->Fill(SumE_cell/m_sf/m_energy);
  }
  else {
    std::cout << "No CaloHit Collection!!!!!" << std::endl;
  }
  if (verbose) {
    std::cout << "Number of cells below threshold:  " << cellsBelowThreshold << std::endl;
    std::cout << "Energy of cells below threshold:  " << lostEnergy << std::endl;
    std::cout << "Total energy 1: " << SumE_cell/m_sf << std::endl;
  }
  SumE_cell = 0.;

  //Cells collection
  if (colHCalPositions_newOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #newCaloPositions:     " << colHCalPositions_new->size()    << std::endl;;
    }
    for (auto& iecl=colHCalPositions_new->begin(); iecl!=colHCalPositions_new->end(); ++iecl)
      {
	if (iecl->core().energy > m_thr){
	  SumE_cell += iecl->core().energy;
	  double r = sqrt(pow(iecl->position().x,2)+pow(iecl->position().y,2));
	  TVector3 vec(iecl->position().x,iecl->position().y,iecl->position().z);
	  double phi = atan2( iecl->position().y, iecl->position().x );
	  double eta = vec.Eta();
	  if (verbose){
	    std::cout << " r " << r << ", eta " << eta << ", phi " << phi << ", energy " << iecl->core().energy << std::endl;
	  }
	  h_ene_x->Fill(iecl->position().x,iecl->core().energy/m_sf);
	  h_ene_y->Fill(iecl->position().y,iecl->core().energy/m_sf);
	  h_ene_z->Fill(iecl->position().z,iecl->core().energy/m_sf);
	  
	  h_ene_r->Fill(r,iecl->core().energy/m_sf);
	  h_ene_phi->Fill(phi,iecl->core().energy/m_sf);
	  h_ene_eta->Fill(eta,iecl->core().energy/m_sf);
	}
      }
    if (verbose) {
      std::cout << "Total cell energy: " << SumE_cell/m_sf << std::endl;
    }
    //Fill histograms
    h_cellEnergy_check->Fill(SumE_cell/m_sf);
  }
  else {
    //    std::cout << "No colHCalPositions_new Collection!!!!!" << std::endl;
  }

  SumE_cell = 0.;
  //Hits collection                                                                                           
  if (colHCalPositionedHits_oldOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #newCaloPositionedHits:     " << colHCalPositionedHits_old->size()    << std::endl;;
    }
    for (auto& iecl=colHCalPositionedHits_old->begin(); iecl!=colHCalPositionedHits_old->end(); ++iecl)
      {
        //if (verbose) std::cout << "HCal cell energy " << iecl->core().energy << std::endl; 
        SumE_cell += iecl->core().energy;
        //if (iecl->core().energy>0.005) {                                                                                     
        double r = sqrt(pow(iecl->position().x,2)+pow(iecl->position().y,2));
        //if (verbose) std::cout << " x " << iecl->position().x << " y " << iecl->position().y << std::endl;
        TVector3 vec(iecl->position().x,iecl->position().y,iecl->position().z);
        double phi = atan2( iecl->position().y, iecl->position().x );
        double eta = vec.Eta();
        if (verbose){
	  std::cout << " r " << r << ", eta " << eta << ", phi " << phi << ", energy " << iecl->core().energy << std::endl;
        }
        h_ene_r_check->Fill(r,iecl->core().energy/m_sf);
        h_ene_phi_check->Fill(phi,iecl->core().energy/m_sf);
        h_ene_eta_check->Fill(eta,iecl->core().energy/m_sf);
      }

    if (verbose) {
      std::cout << "Total energy 2: " << SumE_cell/m_sf << std::endl;
    }
    //Fill histograms                                                                                                                                                                                                                        
    h_cellEnergy_check->Fill(SumE_cell/m_sf);
  }
  else {
    // std::cout << "No colHCalPositionedHits_old Collection!!!!!" << std::endl;
  }
}

void CellAnalysis::finishLoop(int aNumEvents, bool aVerbose) {
  std::cout << "Total energy: " << h_cellEnergy->GetMean() << ", above threshold " << m_thr << ", check " <<  h_cellEnergy_check->GetMean() << std::endl;
  aNumEvents = h_ene_r->GetEntries();
  int aSecNumEvents = h_ene_r_check->GetEntries();

  h_ene_r->Scale(1./h_ene_r->Integral());
  h_ene_phi->Scale(1./h_ene_phi->Integral());
  h_ene_eta->Scale(1./h_ene_eta->Integral());
  
  h_ene_r_check->Scale(1./h_ene_r_check->Integral());
  h_ene_phi_check->Scale(1./h_ene_phi_check->Integral());
  h_ene_eta_check->Scale(1./h_ene_eta_check->Integral());

  std::cout << "Integral phi " << h_ene_phi->Integral() << " check " << h_ene_phi_check->Integral() << std::endl;
}
