#include "SimpleAnalysis.h"

#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenVertexCollection.h"
#include "datamodel/PositionedCaloHitCollection.h"

// ROOT
#include "TH1F.h"

// STL
#include <iostream>

SimpleAnalysis::SimpleAnalysis(double aEnergy, double aSf):
  m_energy(aEnergy), m_sf(aSf), hHitEnergy(nullptr), hCellEnergy(nullptr), hGenPt(nullptr)  {
  Initialize_histos();
}
SimpleAnalysis::~SimpleAnalysis(){}


void SimpleAnalysis::Initialize_histos() {
  hHitEnergy = new TH1F("hHitEnergy","", 10000, 0, 100);
  m_histograms.push_back(hHitEnergy);

  hCellEnergy = new TH1F("hCellenergy","", 100000, 0, 10000);
  m_histograms.push_back(hCellEnergy);

  hGenPt = new TH1F("hGenPt","", 100, m_energy-0.2*m_energy, m_energy+0.2*m_energy);
  m_histograms.push_back(hGenPt);
}

void SimpleAnalysis::processEvent(podio::EventStore& store, int aEventId, bool verbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::PositionedCaloHitCollection*     colHCalPositionedHits(nullptr);

  bool colMCParticlesOK = store.get("GenParticles", colMCParticles);
  bool colHCalPositionedHitsOK     = store.get("HCalPositionedHits" , colHCalPositionedHits);

  //Total hit energy per event
  double SumE_hit_hcal = 0.;

  //PositionedHits collection
  if (colHCalPositionedHitsOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #HCalPositionedHits:     " << colHCalPositionedHits->size()    << std::endl;
    }
    //Loop through the collection
    for (auto& iecl=colHCalPositionedHits->begin(); iecl!=colHCalPositionedHits->end(); ++iecl) {
      //      if (verbose) std::cout << "HCal hit energy " << iecl->core().energy << std::endl;
	//if (verbose) std::cout << "used HCal hit energy " << iecl->core().energy << std::endl;
	SumE_hit_hcal += iecl->core().energy;
    }

    if (verbose) std::cout << "Total hit energy (GeV): " << SumE_hit_hcal << " total cell energy (GeV): " << SumE_hit_hcal/m_sf*100 << " hit collection size: " << colHCalPositionedHits->size() << std::endl;

    //Fill histograms
    hHitEnergy->Fill(SumE_hit_hcal);
    hCellEnergy->Fill(SumE_hit_hcal/m_sf*100);
  }
  else {
    if (verbose) {
      std::cout << "No CaloPositionedHits Collection!!!!!" << std::endl;
    }
  }

  //MCParticle and Vertices collection
  if (colMCParticlesOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #MCTruthParticles:     " << colMCParticles->size()    << std::endl;
    }
    //Loop through the collection
    for (auto& iparticle=colMCParticles->begin(); iparticle!=colMCParticles->end(); ++iparticle) {
      //Fill histogram
      hGenPt->Fill( sqrt( pow(iparticle->core().p4.px,2)+
        pow(iparticle->core().p4.py,2) ) );
    }
  }
  else {
    if (verbose) {
      std::cout << "No MCTruth info available" << std::endl;
    }
  }
}

void SimpleAnalysis::finishLoop(int aNumEvents, bool aVerbose) {
  std::cout << "Total energy: " << hCellEnergy->GetMean() << std::endl;
}
