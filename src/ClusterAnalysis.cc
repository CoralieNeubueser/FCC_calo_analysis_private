#include "ClusterAnalysis.h"

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
#include "TVector2.h"
#include "TVector3.h"

// STL
#include <iostream>
#include <string>
#include <sstream>
#include <bitset>

ClusterAnalysis::ClusterAnalysis(double aEnergy, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal, const std::string& aClusterCollection, const std::string& aEcalCollection, bool normalizeToMCParticle, double deltaRCut) : m_energy(aEnergy), ecal_decoder(aBitfieldEcal), hcal_decoder(aBitfieldHcal), clusterCollection(aClusterCollection), ecal_collection(aEcalCollection), m_norm(normalizeToMCParticle), m_deltaRCut(deltaRCut){
  Initialize_histos();
}

ClusterAnalysis::~ClusterAnalysis() {}


void ClusterAnalysis::Initialize_histos() {
  double sigma = sqrt(pow(0.8*sqrt(m_energy)+0.02*(m_energy),2));

  h_clusters = new TH1F("h_clusters","", 4*std::ceil(4*sigma)+std::ceil(3*sigma), std::floor(m_energy-4*sigma), std::ceil(m_energy+3*sigma));
  m_histograms.push_back(h_clusters);

  double Et = m_energy / cosh(0);

  h_ET = new TH1F("h_ET","",100,Et-Et/2,Et+Et/2);
  m_histograms.push_back(h_ET);

  h_nClusters = new TH1F("h_nClusters","",100,0,100);
  m_histograms.push_back(h_nClusters);

  h_sharedClusters = new TH1F("h_sharedClusters","",100,0,100);
  m_histograms.push_back(h_sharedClusters);

  h_clusterEnergy = new TH1F("h_clusterEnergy","",100,0,m_energy+m_energy/2.);
  m_histograms.push_back(h_clusterEnergy);

  h_clusterCells = new TH1F("h_clusterCells","",1000,0,10000);
  m_histograms.push_back(h_clusterCells);

  h_numCells = new TH1F("h_numCells","", 1000, 0, 1000000);
  m_histograms.push_back(h_numCells);

  h_clusterCellTypes = new TH1F("h_clusterCellTypes","",5,-0.5,4.5);
  m_histograms.push_back(h_clusterCellTypes);

  h_phiRec = new TH1F("h_phiRec","",1000,-.4,.4);
  m_histograms.push_back(h_phiRec);

  h_etaRec = new TH1F("h_etaRec","",1000,-.1,.1);
  m_histograms.push_back(h_etaRec);

  h_thetaRec = new TH1F("h_thetaRec","",1000,-.1,.1);
  m_histograms.push_back(h_thetaRec);

  h_cellId = new TH1F("h_cellId","", 1000, 0,5000e6);
  m_histograms.push_back(h_cellId);

  int etaBinsHCal = 20;
  int etaBinsECal = 40;
  double etaMin = .25;
  double etaMinECal = .2;

  h_etaphi1 = new TH2F("h_etaphi1","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi2 = new TH2F("h_etaphi2","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi3 = new TH2F("h_etaphi3","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi4 = new TH2F("h_etaphi4","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi5 = new TH2F("h_etaphi5","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi6 = new TH2F("h_etaphi6","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi7 = new TH2F("h_etaphi7","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi8 = new TH2F("h_etaphi8","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphi9 = new TH2F("h_etaphi9","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi10 = new TH2F("h_etaphi10","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi11 = new TH2F("h_etaphi11","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi12 = new TH2F("h_etaphi12","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi13 = new TH2F("h_etaphi13","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi14 = new TH2F("h_etaphi14","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi15 = new TH2F("h_etaphi15","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi16 = new TH2F("h_etaphi16","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi17 = new TH2F("h_etaphi17","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphi18 = new TH2F("h_etaphi18","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
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

  h_etaphiSharedCluster1 = new TH2F("h_etaphiSharedCluster1","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster2 = new TH2F("h_etaphiSharedCluster2","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster3 = new TH2F("h_etaphiSharedCluster3","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster4 = new TH2F("h_etaphiSharedCluster4","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster5 = new TH2F("h_etaphiSharedCluster5","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster6 = new TH2F("h_etaphiSharedCluster6","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster7 = new TH2F("h_etaphiSharedCluster7","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster8 = new TH2F("h_etaphiSharedCluster8","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster9 = new TH2F("h_etaphiSharedCluster9","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster10 = new TH2F("h_etaphiSharedCluster10","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster11 = new TH2F("h_etaphiSharedCluster11","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster12 = new TH2F("h_etaphiSharedCluster12","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster13 = new TH2F("h_etaphiSharedCluster13","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster14 = new TH2F("h_etaphiSharedCluster14","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster15 = new TH2F("h_etaphiSharedCluster15","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster16 = new TH2F("h_etaphiSharedCluster16","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster17 = new TH2F("h_etaphiSharedCluster17","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiSharedCluster18 = new TH2F("h_etaphiSharedCluster18","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  m_histograms.push_back(h_etaphiSharedCluster1);
  m_histograms.push_back(h_etaphiSharedCluster2);
  m_histograms.push_back(h_etaphiSharedCluster3);
  m_histograms.push_back(h_etaphiSharedCluster4);
  m_histograms.push_back(h_etaphiSharedCluster5);
  m_histograms.push_back(h_etaphiSharedCluster6);
  m_histograms.push_back(h_etaphiSharedCluster7);
  m_histograms.push_back(h_etaphiSharedCluster8);
  m_histograms.push_back(h_etaphiSharedCluster9);
  m_histograms.push_back(h_etaphiSharedCluster10);
  m_histograms.push_back(h_etaphiSharedCluster11);
  m_histograms.push_back(h_etaphiSharedCluster12);
  m_histograms.push_back(h_etaphiSharedCluster13);
  m_histograms.push_back(h_etaphiSharedCluster14);
  m_histograms.push_back(h_etaphiSharedCluster15);
  m_histograms.push_back(h_etaphiSharedCluster16);
  m_histograms.push_back(h_etaphiSharedCluster17);
  m_histograms.push_back(h_etaphiSharedCluster18);

  h_etaphiType1 = new TH2F("h_etaphiType1","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiType2 = new TH2F("h_etaphiType2","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiType3 = new TH2F("h_etaphiType3","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiType4 = new TH2F("h_etaphiType4","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiType5 = new TH2F("h_etaphiType5","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiType6 = new TH2F("h_etaphiType6","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiType7 = new TH2F("h_etaphiType7","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiType8 = new TH2F("h_etaphiType8","", 800,-4,4, 704, -TMath::Pi(), TMath::Pi());
  h_etaphiType9 = new TH2F("h_etaphiType9","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiType10 = new TH2F("h_etaphiType10","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiType11 = new TH2F("h_etaphiType11","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiType12 = new TH2F("h_etaphiType12","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiType13 = new TH2F("h_etaphiType13","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiType14 = new TH2F("h_etaphiType14","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiType15 = new TH2F("h_etaphiType15","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiType16 = new TH2F("h_etaphiType16","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiType17 = new TH2F("h_etaphiType17","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  h_etaphiType18 = new TH2F("h_etaphiType18","", etaBinsHCal,-etaMin,etaMin, 256, -TMath::Pi(), TMath::Pi());
  m_histograms.push_back(h_etaphiType1);
  m_histograms.push_back(h_etaphiType2);
  m_histograms.push_back(h_etaphiType3);
  m_histograms.push_back(h_etaphiType4);
  m_histograms.push_back(h_etaphiType5);
  m_histograms.push_back(h_etaphiType6);
  m_histograms.push_back(h_etaphiType7);
  m_histograms.push_back(h_etaphiType8);
  m_histograms.push_back(h_etaphiType9);
  m_histograms.push_back(h_etaphiType10);
  m_histograms.push_back(h_etaphiType11);
  m_histograms.push_back(h_etaphiType12);
  m_histograms.push_back(h_etaphiType13);
  m_histograms.push_back(h_etaphiType14);
  m_histograms.push_back(h_etaphiType15);
  m_histograms.push_back(h_etaphiType16);
  m_histograms.push_back(h_etaphiType17);
  m_histograms.push_back(h_etaphiType18);

  h_etaphiEnergy1 = new TH2F("h_etaphiEnergy1","", etaBinsECal,-etaMinECal,etaMinECal, etaBinsECal, -etaMinECal, etaMinECal);
  h_etaphiEnergy2 = new TH2F("h_etaphiEnergy2","", etaBinsECal,-etaMinECal,etaMinECal, etaBinsECal, -etaMinECal, etaMinECal);
  h_etaphiEnergy3 = new TH2F("h_etaphiEnergy3","", etaBinsECal,-etaMinECal,etaMinECal, etaBinsECal, -etaMinECal, etaMinECal);
  h_etaphiEnergy4 = new TH2F("h_etaphiEnergy4","", etaBinsECal,-etaMinECal,etaMinECal, etaBinsECal, -etaMinECal, etaMinECal);
  h_etaphiEnergy5 = new TH2F("h_etaphiEnergy5","", etaBinsECal,-etaMinECal,etaMinECal, etaBinsECal, -etaMinECal, etaMinECal);
  h_etaphiEnergy6 = new TH2F("h_etaphiEnergy6","", etaBinsECal,-etaMinECal,etaMinECal, etaBinsECal, -etaMinECal, etaMinECal);
  h_etaphiEnergy7 = new TH2F("h_etaphiEnergy7","", etaBinsECal,-etaMinECal,etaMinECal, etaBinsECal, -etaMinECal, etaMinECal);
  h_etaphiEnergy8 = new TH2F("h_etaphiEnergy8","", etaBinsECal,-etaMinECal,etaMinECal, etaBinsECal, -etaMinECal, etaMinECal);
  h_etaphiEnergy9 = new TH2F("h_etaphiEnergy9","", etaBinsHCal,-etaMin,etaMin, etaBinsHCal, -etaMin, etaMin);
  h_etaphiEnergy10 = new TH2F("h_etaphiEnergy10","", etaBinsHCal,-etaMin,etaMin, etaBinsHCal, -etaMin, etaMin);
  h_etaphiEnergy11 = new TH2F("h_etaphiEnergy11","", etaBinsHCal,-etaMin,etaMin, etaBinsHCal, -etaMin, etaMin);
  h_etaphiEnergy12 = new TH2F("h_etaphiEnergy12","", etaBinsHCal,-etaMin,etaMin, etaBinsHCal, -etaMin, etaMin);
  h_etaphiEnergy13 = new TH2F("h_etaphiEnergy13","", etaBinsHCal,-etaMin,etaMin, etaBinsHCal, -etaMin, etaMin);
  h_etaphiEnergy14 = new TH2F("h_etaphiEnergy14","", etaBinsHCal,-etaMin,etaMin, etaBinsHCal, -etaMin, etaMin);
  h_etaphiEnergy15 = new TH2F("h_etaphiEnergy15","", etaBinsHCal,-etaMin,etaMin, etaBinsHCal, -etaMin, etaMin);
  h_etaphiEnergy16 = new TH2F("h_etaphiEnergy16","", etaBinsHCal,-etaMin,etaMin, etaBinsHCal, -etaMin, etaMin);
  h_etaphiEnergy17 = new TH2F("h_etaphiEnergy17","", etaBinsHCal,-etaMin,etaMin, etaBinsHCal, -etaMin, etaMin);
  h_etaphiEnergy18 = new TH2F("h_etaphiEnergy18","", etaBinsHCal,-etaMin,etaMin, etaBinsHCal, -etaMin, etaMin);
  m_histograms.push_back(h_etaphiEnergy1);
  m_histograms.push_back(h_etaphiEnergy2);
  m_histograms.push_back(h_etaphiEnergy3);
  m_histograms.push_back(h_etaphiEnergy4);
  m_histograms.push_back(h_etaphiEnergy5);
  m_histograms.push_back(h_etaphiEnergy6);
  m_histograms.push_back(h_etaphiEnergy7);
  m_histograms.push_back(h_etaphiEnergy8);
  m_histograms.push_back(h_etaphiEnergy9);
  m_histograms.push_back(h_etaphiEnergy10);
  m_histograms.push_back(h_etaphiEnergy11);
  m_histograms.push_back(h_etaphiEnergy12);
  m_histograms.push_back(h_etaphiEnergy13);
  m_histograms.push_back(h_etaphiEnergy14);
  m_histograms.push_back(h_etaphiEnergy15);
  m_histograms.push_back(h_etaphiEnergy16);
  m_histograms.push_back(h_etaphiEnergy17);
  m_histograms.push_back(h_etaphiEnergy18);
 
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

  s_EPhi = new THStack("s_EPhi","");
  m_stacks.push_back(s_EPhi);

}

void ClusterAnalysis::processEvent(podio::EventStore& aStore, int aEventId, bool verbose) {

  std::cout << "Event number " << aEventId << std::endl;
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

  //MCParticle and Vertices collection
  if (colMCParticlesOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #MCTruthParticles:     " << colMCParticles->size()    << std::endl;
    }
    if (m_norm){
      //Loop through the collection
      for (auto& iparticle=colMCParticles->begin(); iparticle!=colMCParticles->end(); ++iparticle) {
	TVector3 particle(iparticle->core().p4.px,iparticle->core().p4.py,iparticle->core().p4.pz);
	etaVertex = particle.Eta();
	phiVertex = particle.Phi();
	thetaVertex = particle.Theta();
	//unit vector
	directionParticle = particle.Unit();
      }
    }
  }
  else {
    //if (verbose) {
    std::cout << "No MCTruth info available" << std::endl;
    //}
  }
  
  if (clusterOK && colECalOK) {
    //    if (verbose) {
    std::cout << " Collections:              " << std::endl;
    std::cout << " -> #newCluster:           " << clusters->size()    << std::endl;;
    std::cout << " -> #newCellPositions:     " << colECal->size()    << std::endl;;
    //}
    // std::cout << " Clusters per event:              " << clusters->size() << std::endl;
    h_nClusters -> Fill(clusters->size());
    int clusterId =0;
    int numClusterCells =0;
    int sharedClusters = 0;

    std::vector<uint64_t> checkCellIds;

    for (auto const& icl=clusters->begin(); icl!=clusters->end(); ++icl){
      fcc::CaloCluster cluster = *icl;
      TVector3 hitPosition(cluster.position().x,cluster.position().y,cluster.position().z);
      h_ene_x->Fill(hitPosition.X());
      h_ene_y->Fill(hitPosition.Y());
      h_ene_z->Fill(hitPosition.Z());
      h_ene_r->Fill(hitPosition.Perp());
      h_ene_eta->Fill(hitPosition.Eta());      
      h_ene_phi->Fill(hitPosition.Phi());      

      double deltaR = sqrt(pow(hitPosition.Eta()-etaVertex,2) + pow(hitPosition.Phi()-phiVertex,2));
      double r = sqrt(pow(hitPosition.X(),2) + pow(hitPosition.Y(),2));
      double eta = hitPosition.Eta();
      double ET =  cluster.core().energy / cosh(eta);
      SumECluster += cluster.core().energy;
      SumET += ET;
      clusterId += 1;

      float sumEnergyClusterCells = 0;

      float cellSystem = 0;
      bool cellsInBoth = false;

      h_clusterEnergy->Fill(cluster.core().energy);
      h_clusterCells->Fill(cluster.hits_size());
      numClusterCells += cluster.hits_size();
      if (verbose) 
	std::cout << "Found " << cluster.hits_size() << " cells in cluster." << std::endl;
      
      if (deltaR > m_deltaRCut)
	continue;

      TH1D* e_phi = new TH1D("e_phi","",256, -TMath::Pi(), TMath::Pi());
      for (uint it = 0; it < cluster.hits_size(); it++){
	// find cell in positionedCellCollection
	for (auto& iecl=colECal->begin(); iecl!=colECal->end(); ++iecl){
	  auto cell = *iecl;
	  if (cluster.hits(it).core().cellId == cell.core().cellId && cell.core().cellId != 0) {
	    if (verbose) {
	      std::cout << "Found cluster cell in CellCollection." << std::endl;
	    }
	    TVector3 hitPosition(cell.position().x,cell.position().y,cell.position().z);
	    double r = sqrt(pow(cell.position().x,2)+pow(cell.position().y,2));
	    TVector3 vec(iecl->position().x,cell.position().y,cell.position().z);
	    double phi = vec.Phi();
	    double eta = vec.Eta();

	    if(cluster.hits(it).core().energy !=  cell.core().energy)
	      std::cout << "Cell energy in cluster : " << cluster.hits(it).core().energy << ", cell energy in caloHitCollection : " << cell.core().energy << "\n"; 

	    checkCellIds.push_back(cell.core().cellId);
	    //std::cout << "Cell types : " << cell.core().bits << std::endl;
	    h_clusterCellTypes->Fill( cell.core().bits );

	    sumEnergyClusterCells += cell.core().energy;

	    e_phi->Fill(phi, cell.core().energy);
	    
	    int systemId = ecal_decoder.value("system",cell.core().cellId);
	    //std::cout << "System id : " << systemId << std::endl;
	    if (systemId != cellSystem && cellSystem!=0)
	      cellsInBoth = true;

	    cellSystem = systemId;
	    
	    int layerId = 0;
	    if ( systemId == 5 ){
	      layerId = ecal_decoder.value("layer",cell.core().cellId);
	    }
	    else if ( systemId == 8 ) {
	      layerId = hcal_decoder.value("layer",cell.core().cellId) + 8;
	    }
	    //std::cout << "Layer id : " << layerId << std::endl;

	    double deltaPhi = TVector2::Phi_mpi_pi( phi-phiVertex );
	    
	    if ((layerId +1) == 1){
	      h_etaphi1->Fill(eta-etaVertex , deltaPhi, clusterId);//cell.core().energy);
	      h_etaphiType1->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);
	      h_etaphiEnergy1->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	    } else if ((layerId +1) == 2) {
	      h_etaphi2->Fill(eta-etaVertex , deltaPhi , clusterId);//cell.core().energy);
	      h_etaphiType2->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);            
              h_etaphiEnergy2->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	    } else if ((layerId +1) == 3) {
	      h_etaphi3->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	      h_etaphiType3->Fill(eta-etaVertex ,deltaPhi , cell.core().bits);//cell.core().energy);               
	      h_etaphiEnergy3->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	    } else if ((layerId +1) == 4) {
	      h_etaphiType4->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy4->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	      h_etaphi4->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	    } else if ((layerId +1) == 5) {
	      h_etaphiType5->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy5->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	      h_etaphi5->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	    } else if ((layerId +1) == 6) {
	      h_etaphiType6->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy6->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	      h_etaphi6->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	    } else if ((layerId +1) == 7) {
	      h_etaphiType7->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy7->Fill(eta-etaVertex , deltaPhi, cell.core().energy);
	      h_etaphi7->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	    } else if ((layerId +1) == 8) {
	      h_etaphiType8->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy8->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	      h_etaphi8->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	    } else if ((layerId +1) == 9) {
	      h_etaphiType9->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy9->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	      h_etaphi9->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	    } else if ((layerId +1) == 10) {
	      h_etaphiType10->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy10->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	      h_etaphi10->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	    } else if ((layerId +1) == 11) {
	      h_etaphiType11->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy11->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	      h_etaphi11->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	    } else if ((layerId +1) == 12) {
	      h_etaphiType12->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy12->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	      h_etaphi12->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	    } else if ((layerId +1) == 13) {
	      h_etaphiType13->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy13->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	      h_etaphi13->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	    } else if ((layerId +1) == 14) {
	      h_etaphi14->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	      h_etaphiType14->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy14->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	    } else if ((layerId +1) == 15) {
	      h_etaphi15->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	      h_etaphiType15->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy15->Fill(eta-etaVertex ,deltaPhi , cell.core().energy);
	    } else if ((layerId +1) == 16) {
	      h_etaphi16->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	      h_etaphiType16->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy16->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	    } else if ((layerId +1) == 17) {
	      h_etaphi17->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	      h_etaphiType17->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy17->Fill(eta-etaVertex , deltaPhi , cell.core().energy);
	    } else if ((layerId +1) == 18) {
	      h_etaphi18->Fill(eta-etaVertex , deltaPhi , clusterId);// cell.core().energy);
	      h_etaphiType18->Fill(eta-etaVertex , deltaPhi , cell.core().bits);//cell.core().energy);                  
              h_etaphiEnergy18->Fill(eta-etaVertex , deltaPhi , cell.core().energy);	    }
	  }
	}
      }
      
      if(cellsInBoth){
	sharedClusters++;
	for (uint it = 0; it < cluster.hits_size(); it++){
	  // find cell in positionedCellCollection
	  for (auto& iecl=colECal->begin(); iecl!=colECal->end(); ++iecl){
	    auto cell = *iecl;
	    if (cluster.hits(it).core().cellId == cell.core().cellId && cell.core().cellId != 0) {
	      if (verbose) {
		std::cout << "Found cluster cell in CellCollection." << std::endl;
	      }
	      TVector3 hitPosition(cell.position().x,cell.position().y,cell.position().z);
	      double r = sqrt(pow(cell.position().x,2)+pow(cell.position().y,2));
	      TVector3 vec(iecl->position().x,cell.position().y,cell.position().z);
	      double phi = vec.Phi();
	      double eta = vec.Eta();
	      int systemId = ecal_decoder.value("system",cell.core().cellId);
	      int layerId = 0;
	      if ( systemId == 5 ){
		layerId = ecal_decoder.value("layer",cell.core().cellId);
	      }
	      else if ( systemId == 8 ) {
		layerId = hcal_decoder.value("layer",cell.core().cellId) + 8;
	      }
	      if ((layerId +1) == 1){
		h_etaphiSharedCluster1->Fill(eta , phi, clusterId);//cell.core().energy);
	      } else if ((layerId +1) == 2) {
		h_etaphiSharedCluster2->Fill(eta , phi , clusterId);//cell.core().energy);
	      } else if ((layerId +1) == 3) {
		h_etaphiSharedCluster3->Fill(eta , phi , clusterId);// cell.core().energy);
	      } else if ((layerId +1) == 4) {
		h_etaphiSharedCluster4->Fill(eta , phi , clusterId);//cell.core().energy);                  
	      } else if ((layerId +1) == 5) {
		h_etaphiSharedCluster5->Fill(eta , phi , clusterId);//cell.core().energy);                  
	      } else if ((layerId +1) == 6) {
		h_etaphiSharedCluster6->Fill(eta , phi , clusterId);//cell.core().energy);                  
	      } else if ((layerId +1) == 7) {
		h_etaphiSharedCluster7->Fill(eta , phi , clusterId);//cell.core().energy);                  
	      } else if ((layerId +1) == 8) {
		h_etaphiSharedCluster8->Fill(eta , phi , clusterId);//cell.core().energy);                  
	      } else if ((layerId +1) == 9) {
		h_etaphiSharedCluster9->Fill(eta , phi , clusterId);// cell.core().energy);
	      } else if ((layerId +1) == 10) {
		h_etaphiSharedCluster10->Fill(eta , phi , clusterId);// cell.core().energy);
	      } else if ((layerId +1) == 11) {
		h_etaphiSharedCluster11->Fill(eta , phi , clusterId);// cell.core().energy);
	      } else if ((layerId +1) == 12) {
		h_etaphiSharedCluster12->Fill(eta , phi , clusterId);// cell.core().energy);
	      } else if ((layerId +1) == 13) {
		h_etaphiSharedCluster13->Fill(eta , phi , clusterId);// cell.core().energy);
	      } else if ((layerId +1) == 14) {
		h_etaphiSharedCluster14->Fill(eta , phi , clusterId);// cell.core().energy);
	      } else if ((layerId +1) == 15) {
		h_etaphiSharedCluster15->Fill(eta , phi , clusterId);// cell.core().energy);
	      } else if ((layerId +1) == 16) {
		h_etaphiSharedCluster16->Fill(eta , phi , clusterId);// cell.core().energy);
	      } else if ((layerId +1) == 17) {
		h_etaphiSharedCluster17->Fill(eta , phi , clusterId);// cell.core().energy);
	      } else if ((layerId +1) == 18) {
		h_etaphiSharedCluster18->Fill(eta , phi , clusterId);// cell.core().energy);
	      }
	    }
	  }
	}
      }
      s_EPhi->Add(e_phi);
      if (fabs(sumEnergyClusterCells - cluster.core().energy) > 0.1)
	std::cout << "Energy sum of cells != cluster energy!!!  " << cluster.core().energy << ", " << sumEnergyClusterCells << std::endl;
	
    }

    std::sort(checkCellIds.begin(), checkCellIds.end());
    auto it = std::adjacent_find(checkCellIds.begin(), checkCellIds.end());
    if (it!=checkCellIds.end())
      std::cout << "Found doublicated cellids " << *it << std::endl; 
    if (verbose) {
      std::cout << "Total cluster energy:                 " << SumECluster << std::endl;
    }
    
    std::cout << "Total number of cells in cluster: "<< checkCellIds.size() << std::endl;
    //Fill histograms
    h_clusters    ->Fill(SumECluster);
    h_ET          ->Fill(SumET);
    h_numCells    ->Fill(numClusterCells);
    h_sharedClusters->Fill(sharedClusters);
    
  }
  else {
    std::cout << "No cluster Collection!!!!!" << std::endl;
  } 
}
void ClusterAnalysis::finishLoop(int aNumEvents, bool aVerbose) {
  std::cout << "Total energy: " << h_clusters->GetMean() << std::endl;
  aNumEvents = h_ene_r->GetEntries();
  
  // h_lambda->Scale(1./(double)aNumEvents);
  h_ene_r  ->Scale(1./h_ene_r->Integral());
  h_ene_phi->Scale(1./h_ene_phi->Integral());
  h_ene_eta->Scale(1./h_ene_eta->Integral());
  
  std::cout << "Integral phi " << h_ene_phi->Integral() << std::endl;
}
