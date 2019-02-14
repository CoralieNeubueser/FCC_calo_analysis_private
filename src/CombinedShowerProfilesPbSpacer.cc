#include "CombinedShowerProfilesPbSpacer.h"

#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenVertexCollection.h"
#include "datamodel/PositionedCaloHitCollection.h"

#include "TVector3.h"
#include "TMath.h"

// STL
#include <vector>
#include <iostream>

CombinedShowerProfilesPbSpacer::CombinedShowerProfilesPbSpacer(double aEnergy):
  m_energy(aEnergy) {
  Initialize_histos();
}

CombinedShowerProfilesPbSpacer::~CombinedShowerProfilesPbSpacer() {}


void CombinedShowerProfilesPbSpacer::Initialize_histos()
{

  double r_x[20] = {1920, 1940, 2030, 2120, 2210, 2300, 2390, 2480, 2570, 2850, 2950, 3050, 3200, 3350, 3500, 3650, 3900, 4150, 4400, 4650}; 
  double x0_x[20] = {0}; 
  x0_x[0] = .3; 
  for (int i=0; i<8; i++){
    x0_x[i+1] = x0_x[i] + layerThickness_ecal[i]/X0_ecal;
  }
  x0_x[9] = 2.7;
  for (int i=9; i<19; i++){
    x0_x[i+1] = x0_x[i] + layerThickness_hcal[i-9]/X0_hcal;
  }

  h_corrEnergy = new TH1F("h_corrEnergy","", 100, m_energy-0.5*m_energy, m_energy+0.5*m_energy);
  m_histograms.push_back(h_corrEnergy);

  h_cellEnergy = new TH1F("h_cellEnergy","", 100, m_energy-0.5*m_energy, m_energy+0.5*m_energy);
  m_histograms.push_back(h_cellEnergy);

  h_radialProfile = new TH1F("h_radialProfile","", 100, 0, 20);
  m_histograms.push_back(h_radialProfile);

  h_longProfile = new TH1F("h_longProfile","", 15, 0, 15);
  m_histograms.push_back(h_longProfile);

  h_longProfile_m = new TH1F("h_longProfile_m","", 100, 0, 5000);
  m_histograms.push_back(h_longProfile_m);

  h_radialProfile_particle = new TH1F("h_radialProfile_particle","", 100, 0, 20);
  m_histograms.push_back(h_radialProfile_particle);

  // long. distance ECal in X0
  h_longProfile_particle = new TH1F("h_longProfile_particle","", 15, 0, 15);
  m_histograms.push_back(h_longProfile_particle);

  h_longProfile_particle_m = new TH1F("h_longProfile_particle_m","", 100, 0, 5000);
  m_histograms.push_back(h_longProfile_particle_m);

  h_ptGen = new TH1F("h_ptGen","", 100, m_energy-0.2*m_energy, m_energy+0.2*m_energy);
  m_histograms.push_back(h_ptGen);

  h_r = new TH1F("h_r","", 19, r_x);
  m_histograms.push_back(h_r);

  h_firstLayer = new TH1F("h_firstLayer","", 100,0,m_energy);
  m_histograms.push_back(h_firstLayer);

  h_firstLayerVsHcalTot = new TH2F("h_firstLayerVsHcalTot","", 100,0,0.2*m_energy,100,0,0.5*m_energy); 
  m_histograms.push_back(h_firstLayerVsHcalTot);

  h_r_hcal = new TH1F("h_r_hcal","", 19, r_x);
  m_histograms.push_back(h_r_hcal);

  //  h_longProfile_particle_hcal = new TH1F("h_longProfile_particle_hcal","", 15, X0_offset, 180);
  h_longProfile_particle_hcal = new TH1F("h_longProfile_particle_hcal","", 15, 0, 15);//19, x0_x);
  m_histograms.push_back(h_longProfile_particle_hcal);

  h_eta_particle = new TH1F("h_eta_particle","", 1000,0,1);
  m_histograms.push_back(h_eta_particle);
}


void CombinedShowerProfilesPbSpacer::processEvent(podio::EventStore& aStore, int aEventId, bool aVerbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::PositionedCaloHitCollection*    colECalPositionedHits(nullptr);
  // const fcc::PositionedCaloHitCollection*    colCalPositionedHits(nullptr);
  const fcc::PositionedCaloHitCollection*    colCalPositions(nullptr);

  bool colMCParticlesOK = aStore.get("GenParticles", colMCParticles);
  bool colECalPositionedHitsOK   = aStore.get("cellECalPositions" , colECalPositionedHits);
  //bool colCalPositionedHitsOK    = aStore.get("HCalPositionedHits" , colCalPositionedHits);
  bool colCalPositionsOK         = aStore.get("cellHCalPositions" , colCalPositions);

  // Energy depsitied in 1st HCAL layer
  E_firstLayer = 0.;
  //Total hit energy per event in ECal
  SumE_hit = 0.;
  //Total hit energy per event in HCal
  SumH_hit = 0.;
  //EM shower axis - assuming single shower per event!!!!
  //Direction of gen. particle
  TVector3 directionParticle(0.,0.,0.);
  //Direction of hits in the first layer
  TVector3 directionHits(0.,0.,0.);

  //MCParticle and Vertices collection
  if (colMCParticlesOK) {
    if (aVerbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #MCTruthParticles:     " << colMCParticles->size()    << std::endl;
    }
    //Loop through the collection
    for (auto& iparticle=colMCParticles->begin(); iparticle!=colMCParticles->end(); ++iparticle) {
      TVector3 particle(iparticle->core().p4.px,iparticle->core().p4.py,iparticle->core().p4.pz);
      //unit vector
      directionParticle = particle.Unit();
      //Fill histograms
      h_ptGen->Fill( sqrt( pow(iparticle->core().p4.px,2)+pow(iparticle->core().p4.py,2) ) );
    }
  }
  else {
    if (aVerbose) {
      std::cout << "No MCTruth info available" << std::endl;
    }
  }


  //Positioned Cells collection
  if (colECalPositionedHitsOK && colCalPositionsOK) {
    if (aVerbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #ECalPositionedHits:    " << colECalPositionedHits->size()    << std::endl;;
      std::cout << " -> #HCalPositionedHits:    " << colCalPositions->size()    << std::endl;;
    }
    TH1F *h_hitR   = new TH1F("h_hitR", "h_hitR", 1000, RcaloMin_ecal, RcaloMin_ecal+1.5*layerThickness);
    TH1F *h_hitPhi = new TH1F("h_hitPhi", "h_hitPhi", 100000, -TMath::Pi(), TMath::Pi());
    TH1F *h_hitEta = new TH1F("h_hitEta", "h_hitEta", 1000000, -EtaMax, EtaMax);

    //Loop through the ECal collection
    //First: find the mean position in the first layer, calculate total energy
    for (auto& iecl=colECalPositionedHits->begin(); iecl!=colECalPositionedHits->end(); ++iecl)
      {
	double hitEnergy = iecl->core().energy;
	SumE_hit += hitEnergy;
	
	TVector3 hit_position(iecl->position().x,iecl->position().y,iecl->position().z);
	double R_firstLayer = RcaloMin_ecal+layerThickness;
	if (hit_position.Pt()<RcaloMin_ecal) {
	  std::cout <<"Hit before the calorimeter????? Please check the value of RcaloMin!"<< std::endl;
	}
	if (hit_position.Pt()<R_firstLayer) {
	  h_hitR->Fill(hit_position.Perp(), hitEnergy);
	  h_hitPhi->Fill(hit_position.Phi(), hitEnergy);
	  h_hitEta->Fill(hit_position.Eta(), hitEnergy);
	}
      }
    TVector3 meanFirstLayer_vector;
    //std::cout << "Low edge " << h_hitEta->GetBinLowEdge(1) << " maximum " << h_hitEta->GetMaximumBin() << " mean " << h_hitEta->GetMean() << std::endl;
    double meanR   = h_hitR  ->GetBinLowEdge(1)+h_hitR  ->GetMaximumBin()*h_hitR  ->GetBinWidth(1)+h_hitR  ->GetBinWidth(1)*0.5;
    double meanPhi = h_hitPhi->GetBinLowEdge(1)+h_hitPhi->GetMaximumBin()*h_hitPhi->GetBinWidth(1)+h_hitPhi->GetBinWidth(1)*0.5;
    double meanEta = h_hitEta->GetBinLowEdge(1)+h_hitEta->GetMaximumBin()*h_hitEta->GetBinWidth(1)+h_hitEta->GetBinWidth(1)*0.5;

    meanFirstLayer_vector.SetPtEtaPhi(meanR, meanEta, meanPhi);
    //direction from hits in the first calorimeter layer
    directionHits = meanFirstLayer_vector.Unit();
    h_eta_particle->Fill( directionParticle.Eta() );

    delete h_hitR;
    delete h_hitPhi;
    delete h_hitEta;

    //Second loop: radial & longitudinal profiles (ECal entries)
    for (auto& iecl=colECalPositionedHits->begin(); iecl!=colECalPositionedHits->end(); ++iecl)
      {
	double hitEnergy = iecl->core().energy;
	TVector3 hitPosition(iecl->position().x,iecl->position().y,iecl->position().z);
	double r = sqrt(pow(iecl->position().x,2)+pow(iecl->position().y,2));
     	// Add energy in last ECal layer
	if ( r > (RcaloMax_ecal - layerThickness)){	
     	  E_firstLayer += hitEnergy/sf_ecal;
	}

	//	std::cout << "radii ecal  = " << r << std::endl;

	//gen. particle direction
	TVector3 showerStart_particle = directionParticle*RcaloMin_ecal;
	TVector3 hitVector_particle = hitPosition-showerStart_particle;
	double hitLong_particle = directionParticle*hitVector_particle;
	double hitRadial_particle = hitVector_particle.Perp(directionParticle);

	//hits in first layer direction
	//Start of the hits in the first active layer
	TVector3 showerStart = directionHits*RcaloMin_ecal;//meanFirstLayer_vector.Mag();
	//TVector3 showerStart = directionHits*RcaloMin;
	TVector3 hitVector = hitPosition-showerStart;
	double hitLong = directionHits*hitVector;
	double hitRadial = hitVector.Perp(directionHits);

	//Fill longitudinal and radial profile histograms
	h_radialProfile_particle->Fill(hitRadial_particle/X0_ecal, hitEnergy);
	h_longProfile_particle  ->Fill(X0_zeroOffset/float(cos(directionParticle.Eta())) + hitLong_particle/X0_ecal, hitEnergy);
	h_longProfile_particle_m->Fill(hitLong_particle, hitEnergy);
	h_radialProfile         ->Fill(hitRadial/X0_ecal, hitEnergy);
	h_longProfile           ->Fill(X0_zeroOffset/float(cos(directionParticle.Eta())) + hitLong/X0_ecal, hitEnergy);
	h_longProfile_m         ->Fill(hitLong, hitEnergy);

	h_r           ->Fill( r, hitEnergy );

      }

    //Second loop: radial & longitudinal profiles (HCal entries)
    for (auto& iecl=colCalPositions->begin(); iecl!=colCalPositions->end(); ++iecl)
      {
	double hitEnergy = iecl->core().energy;
	TVector3 hitPosition(iecl->position().x,iecl->position().y,iecl->position().z);
	double r = sqrt(pow(iecl->position().x,2)+pow(iecl->position().y,2));
     	// Add energy in first HCal layer
     	if ( r < (RcaloMin_hcal+layerThickness)){	
     	  E_firstLayer += hitEnergy/sf_hcal;
	}
	//gen. particle direction
	TVector3 showerStart_particle = directionParticle*RcaloMin_ecal;
	TVector3 showerStartFromHcal_particle = directionParticle*RcaloMin_hcal;
	TVector3 hitVector_particle = hitPosition-showerStart_particle;
	TVector3 hitVectorFromHcal_particle = hitPosition-showerStartFromHcal_particle;
	
	double hitLong_particle = directionParticle*hitVector_particle;
	double hitLongFromHcal_particle = directionParticle*hitVectorFromHcal_particle;
	double hitRadial_particle = hitVector_particle.Perp(directionParticle);

	//hits in first layer direction
	//Start of the hits in the first active layer
	TVector3 showerStart     = directionHits*RcaloMin_hcal;
	//TVector3 showerStart = directionHits*RcaloMin;
	TVector3 hitVector = hitPosition-showerStart;
	double hitLong = directionHits*hitVector;
	double hitRadial = hitVector.Perp(directionHits);
       
	//	std::cout << "radii hcal  = " << r << std::endl;
	//std::cout << "cos(eta) hcal  = " << float(cos(directionParticle.Eta())) << std::endl;
	
	double positionLong     = hitLong_particle/X0_hcal;
	double positionLongFromHcal     = hitLongFromHcal_particle/X0_hcal;
	double positionLongHits = hitLong/X0_hcal;
	
	//Fill longitudinal and radial profile histograms
	h_radialProfile_particle   ->Fill(hitRadial_particle/X0_hcal, hitEnergy);
	h_longProfile_particle_hcal->Fill(X0_offset/float(cos(directionParticle.Eta())) + positionLongFromHcal, hitEnergy);
	h_longProfile_particle_m   ->Fill(hitLong_particle, hitEnergy);
	h_radialProfile            ->Fill(hitRadial/X0_hcal, hitEnergy);
	h_longProfile              ->Fill(positionLongHits, hitEnergy);
	h_longProfile_m            ->Fill(hitLong, hitEnergy);

	h_r_hcal      ->Fill( r, hitEnergy );

	SumH_hit += hitEnergy;
      }
    
    // correction dependent on last ECAL and first HCAL layer energy
    h_firstLayer->Fill(E_firstLayer);

    h_firstLayerVsHcalTot->Fill(E_firstLayer, (m_energy-SumH_hit/sf_hcal-SumE_hit/sf_ecal));
    // Simple energy sum    
    Ecorr = SumE_hit + SumH_hit;
    
    //Fill histograms
    h_corrEnergy ->Fill(Ecorr);
    // on EM scale
    h_cellEnergy ->Fill(SumE_hit + SumH_hit*0.029/sf_hcal);

    // // normalise to pion energy
    // h_r                     ->Scale(1./(double)Ecorr);
    // h_r_hcal                ->Scale(1./(double)Ecorr);
    // h_longProfile           ->Scale(1./(double)Ecorr);
    // h_longProfile_m         ->Scale(1./(double)Ecorr);
    // h_longProfile_particle  ->Scale(1./(double)Ecorr);
    // h_longProfile_particle_hcal->Scale(1./(double)Ecorr);
    // h_longProfile_particle_m->Scale(1./(double)Ecorr);

    if (aVerbose) std::cout << "Total cell energy, ECAL (GeV): " << SumE_hit/0.1915 << " total cell energy, HCAL (GeV): " << SumH_hit/sf_hcal << ", corrected total energy : " << Ecorr << std::endl;

  }
  else {
    if (aVerbose) {
      std::cout << "No CaloPositionedHits Collection!!!!!" << std::endl;
    }
  }
}

void CombinedShowerProfilesPbSpacer::finishLoop(int aNumEvents, bool aVerbose) {
  h_r                     ->Scale(1./(double)aNumEvents);
  h_r_hcal                ->Scale(1./(double)aNumEvents);
  h_eta_particle          ->Scale(1./(double)aNumEvents);
  h_radialProfile         ->Scale(1./(double)aNumEvents);
  h_longProfile           ->Scale(1./(double)aNumEvents);
  h_longProfile_m         ->Scale(1./(double)aNumEvents);
  h_radialProfile_particle->Scale(1./(double)aNumEvents);
  h_longProfile_particle  ->Scale(1./(double)aNumEvents);
  h_longProfile_particle_hcal->Scale(1./(double)aNumEvents);
  h_longProfile_particle_m->Scale(1./(double)aNumEvents);

  h_longProfile_particle_m->GetXaxis()->SetRangeUser(0,3000);
  h_longProfile_particle_m->GetXaxis()->SetTitle("longitudinal distance from shower start [mm]");
  h_longProfile_particle_m->GetYaxis()->SetTitle("#LTE_{rec}#GT [GeV]");
    
    std::cout << "BinsMin " << h_r->GetBinCenter(h_r->GetNbinsX()) << std::endl;
  //  std::cout << "BinsMax " << h_r->GetMaximumBin() << std::endl;
  //  std::cout << "End of loop" << std::endl;
}
