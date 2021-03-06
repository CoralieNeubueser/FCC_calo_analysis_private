#ifndef COMBINEDSHOWERPROFILES_H
#define COMBINEDSHOWERPROFILES_H

#include "BaseAnalysis.h"

#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class CombinedShowerProfiles: public BaseAnalysis {

 public:
  CombinedShowerProfiles(double aEnergy);
  ~CombinedShowerProfiles();

  virtual void processEvent(podio::EventStore& store, int aEventId, bool verbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;

  /// Energy corrected for loss between the Calorimeters
  TH1F* h_corrEnergy;
  /// Calibrated hit energy (SF)
  TH1F* h_cellEnergy;
  /// Radial profile, direction of the shower from hits in the first layer
  TH1F* h_radialProfile;
  /// Long. profile, direction of the shower from hits in the first layer
  TH1F* h_longProfile;
  /// Long. profile in mm, direction of the shower from hits in the first layer
  TH1F* h_longProfile_m;
  /// Radial profile, direction of the shower from the generated particle
  TH1F* h_radialProfile_particle;
  /// Long. profile, direction of the shower from the generated particle
  TH1F* h_longProfile_particle;
  /// Long. profile HCAL, direction of the shower from the generated particle
  TH1F* h_longProfile_particle_hcal;
  /// Long. profile in mm, direction of the shower from the generated particle
  TH1F* h_longProfile_particle_m;
  /// pt of the generated particle
  TH1F* h_ptGen;
  /// eta from hits' direction
  TH1F* h_eta;
  /// eta from particles' direction
  TH1F* h_eta_particle;
  /// radius
  TH1F* h_r;
  /// radius HCAL
  TH1F* h_r_hcal;
  /// energy in first HCAL layer
  TH1F* h_firstLayer;
  /// energy in first HCAL layer
  TH2F* h_firstLayerVsHcalTot;

 private:
  void Initialize_histos();
  
  // Parameters of the correction function 2.67125 + x* 0.802343 ECal last layer+ HCAl first layer for whole energy range
  double a = 2.67; //in GeV
  double b = 0.80; //

  // from correcting response -4.32945 + x*0.916021
  double a1 = -4.32945; //in GeV
  double b1 = 0.916021; //

  // Sampling fracions for EM scale
  double sf_ecal = .188; 
  double sf_hcal = .032;

  double Ecorr;
  double E_firstLayer;
  double m_energy;
  double SumE_hit;    // Total hit energy per event in ECAL
  double SumH_hit;    // Total hit energy per event in HCAL

  double layerThickness_ecal[8] = {20,90,90,90,90,90,90,90};
  // first HCAL layer for 1.8X0 from ECAL cryo
  double layerThickness_hcal[10] = {100,100,150,150,150,150,250,250,250,250};
  double cellVol_ecal[8] = {4.315, 20.16, 21.53, 22.787, 24.374, 25.66, 26.68, 28.06};
  double cellVol_hcal[10] = {87087, 90049, 140626, 147288, 154000, 160709, 282645, 301218, 319787, 338351};

  double r_x[20] = {1920, 1940, 2030, 2120, 2210, 2300, 2390, 2480, 2570, 2860.5, 2960.5, 3060.5, 3210.5, 3360.5, 3510.5, 3660.5, 3910.5, 4160.5, 4410.5, 4660.5}; 
  double x0_x[20] = {0}; 

  const double RcaloMin_ecal = 1920.; // min ECAL, Cryostat 5cm + LAr bath 9cm increased radius from 1780 to 1920
  const double RcaloMax_ecal = 2570.; // max ECAL, Cryostat 10cm + 5cm LAr bath decreased radius from 2720 to 2570
  const double RcaloMin_hcal = 2860.5; // min HCAL
  const double layerThickness = 100; // mainly used as first layer thickness, 6 ECAL
  const double EtaMax = 10.0;
  const double X0_ecal = 294.1; // in mm for ECAL, eta=0
  const double X0_hcal = 200.7; // in mm for HCAL, eta=0
  const double X0_zeroOffset = 1.8;  // X0 before the ECAL starts
  const double X0_offset = 32.5;  // X0 before the HCAL starts
  //  const double long_offset = 1018;  // mm from start in ECal to beginning of HCal with eta=0.36
  const double long_offset = 930; // mm from start in ECal to beginning of HCal with eta=0.36 
};

#endif /* COMBINEDSHOWERPROFILES_H */
