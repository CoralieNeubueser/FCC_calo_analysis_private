#include "BasicClusterAnalysis.h"

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
#include "TMath.h"
#include "TVector3.h"

// STL
#include <iostream>
#include <string>
#include <sstream>
#include <bitset>

BasicClusterAnalysis::BasicClusterAnalysis(double aEnergy, const std::string& aBitfieldEcal, const std::string& aClusterCollection, const std::string& aEcalCollection) : m_energy(aEnergy), ecal_decoder(aBitfieldEcal), clusterCollection(aClusterCollection), ecal_collection(aEcalCollection) {
  Initialize_histos();
}

BasicClusterAnalysis::~BasicClusterAnalysis() {}


void BasicClusterAnalysis::Initialize_histos() {
  double sigma = sqrt(pow(0.6*sqrt(m_energy)+0.1*(m_energy),2));

  h_clusters = new TH1F("h_clusters","", 1000,0,2);
  m_histograms.push_back(h_clusters);

  double Et = m_energy / cosh(0);

  h_ET = new TH1F("h_ET","",1000,0,2);
  m_histograms.push_back(h_ET);

  h_nClusters = new TH1F("h_nClusters","",100000,0,20000);
  m_histograms.push_back(h_nClusters);

  h_clusterEnergy = new TH1F("h_clusterEnergy","",100000,0,20000);
  m_histograms.push_back(h_clusterEnergy);

  h_clusterCells = new TH1F("h_clusterCells","",100000,0,20000);
  m_histograms.push_back(h_clusterCells);

  h_clusterPt = new TH1F("h_clusterPt","", 100000,0,20000);
  m_histograms.push_back(h_clusterPt);

  h_corrEnergyNumber = new TH2F("h_corrEnergyNumber","",100,0,m_energy+m_energy/2.,1000,0,10000);
  m_histograms.push_back(h_corrEnergyNumber);

  h_phiRec = new TH1F("h_phiRec","",1000,-.4,.4);
  m_histograms.push_back(h_phiRec);

  h_etaRec = new TH1F("h_etaRec","",1000,-.1,.1);
  m_histograms.push_back(h_etaRec);

  h_thetaRec = new TH1F("h_thetaRec","",1000,-.1,.1);
  m_histograms.push_back(h_thetaRec);

  h_cellId = new TH1F("h_cellId","", 1000, 0,5000e6);
  m_histograms.push_back(h_cellId);

  h_ene_eta = new TH1F("h_ene_eta","", int(4./0.025),-4,4);
  m_histograms.push_back(h_ene_eta);

  h_ene_phi = new TH1F("h_ene_phi","", 128, -TMath::Pi(), TMath::Pi());
  m_histograms.push_back(h_ene_phi);

  h_ene_x = new TH1F("h_ene_x","", 10000, -5000.,5000.);
  m_histograms.push_back(h_ene_x);

  h_ene_y = new TH1F("h_ene_y","", 10000, -5000.,5000.);
  m_histograms.push_back(h_ene_y);

  h_ene_z = new TH1F("h_ene_z","",  510, -4600.,4600.);
  m_histograms.push_back(h_ene_z);

  h_ene_r = new TH1F("h_ene_r","", 1000, -2500.,2500.);
  m_histograms.push_back(h_ene_r);

  h_lambda = new TH1F("h_lambda","", 13, 0, 13);
  m_histograms.push_back(h_lambda);
}

void BasicClusterAnalysis::processEvent(podio::EventStore& aStore, int aEventId, bool verbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::PositionedCaloHitCollection*     colECal(nullptr);
  const fcc::CaloClusterCollection*     clusters(nullptr);

  bool colMCParticlesOK       = aStore.get("GenParticles", colMCParticles);
  // Cell Collection
  bool colECalOK     = aStore.get(ecal_collection , colECal);
  // Cluster Collection
  bool clusterOK     = aStore.get(clusterCollection , clusters);

  double etaVertex = 0;
  double phiVertex = 0;
  double thetaVertex = 0;
  //Direction of gen. particle
  TVector3 directionParticle(0.,0.,0.);
  SumECluster = 0;
  double Et = m_energy / cosh(0);
  double pTGen = 0;
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
      pTGen  = m_energy*particle.Unit().Perp();
     //unit vector
      directionParticle = particle.Unit();
    }
  }
  else {
    //if (verbose) {
      std::cout << "No MCTruth info available" << std::endl;
      //}
  }
  
  if (clusterOK) {
    if (verbose) {
      std::cout << " Collections:              " << std::endl;
      std::cout << " -> #newCluster:           " << clusters->size()    << std::endl;;
    }
    // std::cout << " Clusters per event:              " << clusters->size() << std::endl;
    h_nClusters -> Fill(clusters->size());
    int clusterId =0;
    for (auto const& icl=clusters->begin(); icl!=clusters->end(); ++icl){
      fcc::CaloCluster cluster = *icl;
      TVector3 hitPosition(cluster.position().x,cluster.position().y,cluster.position().z);
      h_ene_x->Fill(hitPosition.X());
      h_ene_y->Fill(hitPosition.Y());
      h_ene_z->Fill(hitPosition.Z());
      h_ene_r->Fill(hitPosition.Perp());
      h_ene_eta->Fill(hitPosition.Eta());      
      h_ene_phi->Fill(hitPosition.Phi());      
      double r   = sqrt(pow(hitPosition.X(),2) + pow(hitPosition.Y(),2));
      double eta = hitPosition.Eta();
      double ET  = cluster.core().energy / cosh(eta);
      double pT  = cluster.core().energy*hitPosition.Unit().Perp();
      SumECluster += cluster.core().energy;
      SumET += ET;

      h_clusterEnergy->Fill(cluster.core().energy);
      h_clusterPt->Fill(pT/pTGen);
      h_clusterCells->Fill(cluster.hits_size());
      h_corrEnergyNumber->Fill(cluster.core().energy,cluster.hits_size());
      if (verbose) 
	std::cout << "Found " << cluster.hits_size() << " cells in cluster." << std::endl;
      clusterId++;
    }
    
    if (verbose) {
      std::cout << "Total cluster energy:                 " << SumECluster << std::endl;
    }
    
    
    //Fill histograms
    h_clusters      ->Fill(SumECluster/m_energy);
    h_ET            ->Fill(SumET/Et);
 
    
  }
  else {
    std::cout << "No cluster Collection!!!!!" << std::endl;
  } 
}
void BasicClusterAnalysis::finishLoop(int aNumEvents, bool aVerbose) {
  std::cout << "Total energy: " << h_clusters->GetMean() << std::endl;
  aNumEvents = h_ene_r->GetEntries();
  
  // h_lambda->Scale(1./(double)aNumEvents);
  h_ene_r  ->Scale(1./h_ene_r->Integral());
  h_ene_phi->Scale(1./h_ene_phi->Integral());
  h_ene_eta->Scale(1./h_ene_eta->Integral());
  
  std::cout << "Integral phi " << h_ene_phi->Integral() << std::endl;
}
