#include "CellSignificances.h"

// FCC-EDM
#include "datamodel/MCParticleCollection.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/PositionedCaloHitCollection.h"
#include "datamodel/CaloHit.h"

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

CellSignificances::CellSignificances(double aEnergy, const std::string& aBitfieldEcal, const std::string& aHitCollection, bool aElectrNoiseAdded) : m_energy(aEnergy), ecal_decoder(aBitfieldEcal),hit_collection(aHitCollection), elecNoise(aElectrNoiseAdded) {
  Initialize_histos();
}

CellSignificances::~CellSignificances() {}


void CellSignificances::Initialize_histos() {

  // Get noise per cellID
  TFile file("/afs/cern.ch/work/c/cneubuse/public/FCChh/cellNoise_map_segHcal_electronicsNoiseLevel.root", "READ");
  TTree* tree = nullptr;
  file.GetObject("noisyCells", tree);
  ULong64_t readCellId;
  double readNoisyCells;
  double readNoisyCellsOffset;
  tree->SetBranchAddress("cellId", &readCellId);
  tree->SetBranchAddress("noiseLevel", &readNoisyCells);
  tree->SetBranchAddress("noiseOffset", &readNoisyCellsOffset);
  for (uint i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    m_map.insert(std::pair<uint64_t, std::pair<double, double>>(readCellId, std::make_pair(readNoisyCells, readNoisyCellsOffset)));
  }

  double sigma = sqrt(pow(0.6*sqrt(m_energy)+0.1*(m_energy),2));

  h_energy_EMlayer1 = new TH1F("h_energy_EMlayer1","",100,-10,10);
  m_histograms.push_back(h_energy_EMlayer1);

  h_energy_EMlayer2 = new TH1F("h_energy_EMlayer2","",100,-10,10);
  m_histograms.push_back(h_energy_EMlayer2);

  h_energy_EMlayer3 = new TH1F("h_energy_EMlayer3","",100,-10,10);
  m_histograms.push_back(h_energy_EMlayer3);

  h_energy_EMlayer8 = new TH1F("h_energy_EMlayer8","",100,-10,10);
  m_histograms.push_back(h_energy_EMlayer8);

  h_energy_HadLayer1 = new TH1F("h_energy_HadLayer1","",100,-10,10);
  m_histograms.push_back(h_energy_HadLayer1);

  h_energy_HadLayer2 = new TH1F("h_energy_HadLayer2","",100,-10,10);
  m_histograms.push_back(h_energy_HadLayer2);

  h_energy_HadLayer6 = new TH1F("h_energy_HadLayer6","",100,-10,10);
  m_histograms.push_back(h_energy_HadLayer6);

  h_energy_HadLayer10 = new TH1F("h_energy_HadLayer10","",100,-10,10);
  m_histograms.push_back(h_energy_HadLayer10);
}

void CellSignificances::processEvent(podio::EventStore& aStore, int aEventId, bool verbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::CaloHitCollection*     colCal(nullptr);
  
  bool colMCParticlesOK       = aStore.get("GenParticles", colMCParticles);
  // Cell Collection
  bool colOK     = aStore.get(hit_collection , colCal);

  if (colOK) {
    if (verbose) {
      std::cout << " Collections:              " << std::endl;
      std::cout << " -> #CellCollection:           " << colCal->size()    << std::endl;
    }

    for (auto const& icl=colCal->begin(); icl!=colCal->end(); ++icl){
      fcc::CaloHit cell = *icl;
      int systemId = ecal_decoder.value("system",cell.core().cellId);
      double thr = m_map[cell.core().cellId].first;
      int layer;
      if (systemId == 5){
	layer = ecal_decoder.value("layer",cell.core().cellId) +1;
      }
      else{
	layer = hcal_decoder->value("layer",cell.core().cellId) +1;
      }
      double cellSig = cell.core().energy / thr;

      if (systemId == 5){
	if(layer == 1)
	  h_energy_EMlayer1->Fill(cellSig);
	else if(layer == 2)
	  h_energy_EMlayer2->Fill(cellSig);
	else if(layer == 3)
	  h_energy_EMlayer3->Fill(cellSig);
	else if(layer == 8)
	  h_energy_EMlayer8->Fill(cellSig);
      }
      else{
	if(layer == 1)
	  h_energy_HadLayer1->Fill(cellSig);
	else if(layer == 2)
	  h_energy_HadLayer2->Fill(cellSig);
	else if(layer == 6)
	  h_energy_HadLayer6->Fill(cellSig);
	else if(layer == 10)
	  h_energy_HadLayer10->Fill(cellSig);
      }
    
    }
    
  }
  else {
    std::cout << "No cell Collection!!!!!" << std::endl;
  } 
}
void CellSignificances::finishLoop(int aNumEvents, bool aVerbose) {
  std::cout << "Mean energy first HCal layer: " << h_energy_HadLayer1->GetMean() << std::endl;
}
