#include "ClusterAnalysisTree.h"

// FCC-EDM
#include "datamodel/MCParticleCollection.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/PositionedCaloHitCollection.h"
#include "datamodel/CaloClusterCollection.h"
#include "datamodel/CaloCluster.h"

// PODIO
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

// ROOT
#include "TFile.h"
#include "TMath.h"
#include "TBranch.h"
#include "TVector3.h"

// STL
#include <iostream>
#include <string>
#include <sstream>
#include <bitset>

ClusterAnalysisTree::ClusterAnalysisTree(double aEnergy, const std::string& aBitfieldEcal, const std::string& aClusterCollection, const std::string& aEcalCollection, const std::string& aRootFileName) : m_energy(aEnergy), ecal_decoder(aBitfieldEcal), clusterCollection(aClusterCollection), ecal_collection(aEcalCollection), m_rootFileName(aRootFileName) {
  clearVariables();
  Initialize_tree();
}

ClusterAnalysisTree::~ClusterAnalysisTree() {}


void ClusterAnalysisTree::Initialize_tree() {
  std::cout << "Open root file :  " << m_rootFileName.c_str() << std::endl;
  m_file = TFile::Open(m_rootFileName.c_str(), "recreate");
  m_tree = new TTree("tree","tree");
  std::cout << " initialising tree .. "   << std::endl;
  m_tree->Branch("totEnergy", &m_Energy);
  m_tree->Branch("cluster", &m_cluster);
  m_tree->Branch("clusterCells", &m_clusterCells);
  m_tree->Branch("clusterEnergy", &m_clusterEnergy);
  m_tree->Branch("cellEta", &m_etaCell);
  m_tree->Branch("cellTheta", &m_thetaCell);
  m_tree->Branch("cellPhi", &m_phiCell);
  m_tree->Branch("cellLayer", &m_layerCell);
  m_tree->Branch("cellType", &m_typeCell);
  m_tree->Branch("cellEnergy", &m_energyCell);
  std::cout << " tree initialised. "   << std::endl;
}

void ClusterAnalysisTree::processEvent(podio::EventStore& aStore, int aEventId, bool verbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::PositionedCaloHitCollection*     colECal(nullptr);
  const fcc::CaloClusterCollection*     clusters(nullptr);

  bool colMCParticlesOK       = aStore.get("GenParticles", colMCParticles);
  // Cell Collection
  bool colECalOK     = aStore.get(ecal_collection , colECal);
  // Cluster Collection
  bool clusterOK     = aStore.get(clusterCollection , clusters);

  clearVariables();

  double etaVertex = 0;
  double phiVertex = 0;
  double thetaVertex = 0;
  //Direction of gen. particle
  TVector3 directionParticle(0.,0.,0.);
  SumECluster = 0;

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
    std::cout << "No MCTruth info available" << std::endl;
  }
  
  if (clusterOK && colECalOK) {
    if (verbose) {
      std::cout << " Collections:              " << std::endl;
      std::cout << " -> #newCluster:           " << clusters->size()    << std::endl;;
      std::cout << " -> #newECalPositions:     " << colECal->size()    << std::endl;;
    }
    // Loop over Clusters
    m_cluster = clusters->size();
    for (auto const& icl=clusters->begin(); icl!=clusters->end(); ++icl){
      fcc::CaloCluster cluster = *icl;
      TVector3 hitPosition(cluster.position().x,cluster.position().y,cluster.position().z);
      double r = sqrt(pow(hitPosition.X(),2) + pow(hitPosition.Y(),2));
      double eta = hitPosition.Eta();
      double theta = hitPosition.Theta();
      double ET =  cluster.core().energy / cosh(eta);
      SumECluster += cluster.core().energy;
      SumET += ET;

      m_clusterEnergy.push_back( cluster.core().energy );
      m_clusterCells.push_back( cluster.hits_size() );

      std::vector<double> etaVec;
      std::vector<double> thetaVec;
      std::vector<double> phiVec;
      std::vector<int> layerVec;
      std::vector<int> typeVec;
      std::vector<double> energyVec;

      if (verbose) 
	std::cout << "Found " << cluster.hits_size() << " cells in cluster." << std::endl;
    
     //auto it=cluster.hits_begin();
      for (uint it = 0; it < cluster.hits_size(); it++){
	// find cell in positionedCellCollection
	for (auto& iecl=colECal->begin(); iecl!=colECal->end(); ++iecl){
	  auto cell = *iecl;
	  // std::cout << "CellID in positionedCaloHitCollection " << cell.core().cellId << std::endl;
	  if (cluster.hits(it).core().cellId == cell.core().cellId && cell.core().cellId != 0) {
	    if (verbose) {
	      std::cout << "Found cluster cell in CellCollection." << std::endl;
	    }
	    TVector3 hitPosition(cell.position().x,cell.position().y,cell.position().z);
	    double r = sqrt(pow(cell.position().x,2)+pow(cell.position().y,2));
	    TVector3 vec(iecl->position().x,cell.position().y,cell.position().z);
	    double phi = atan2( cell.position().y, cell.position().x );
	    double eta = vec.Eta();
	    double theta = vec.Theta();
	    int layerId = ecal_decoder.value("layer",cell.core().cellId);
	    
	    etaVec.push_back(eta - etaVertex);
	    phiVec.push_back(phi - phiVertex);
	    thetaVec.push_back(theta - thetaVertex);
	    layerVec.push_back(layerId);
            typeVec.push_back(cell.core().bits);
            energyVec.push_back(cell.core().energy);
	  }
	  else 
	    continue;
	}
      }
      m_etaCell.push_back(etaVec);
      m_thetaCell.push_back(thetaVec);
      m_phiCell.push_back(phiVec);
      m_layerCell.push_back(layerVec);
      m_typeCell.push_back(typeVec);
      m_energyCell.push_back(energyVec);
      // }
    }      
    if (verbose) {
      std::cout << "Total cluster energy:                 " << SumECluster << std::endl;
    }
        
    m_Energy = SumECluster;

    m_tree->Fill();
  }
  else {
    std::cout << "No cluster Collection!!!!!" << std::endl;
  } 
}

void ClusterAnalysisTree::finishLoop(int aNumEvents, bool aVerbose) {
  std::cout << "Total number events: " << m_tree->GetEntries() << std::endl;
  aNumEvents = m_tree->GetEntries();
  
  if (m_file and m_file->IsOpen()) {
    //  m_file->cd();
    //  m_tree->Write();
    m_file->Write();
    delete m_file;
  }
}

void ClusterAnalysisTree::clearVariables() {

  m_Energy = 0.;
  m_cluster = 0;
  m_clusterEnergy.clear();
  m_clusterCells.clear();
  m_etaCell.clear();
  m_thetaCell.clear();
  m_phiCell.clear();
  m_layerCell.clear();
  m_typeCell.clear();
  m_energyCell.clear();

}
