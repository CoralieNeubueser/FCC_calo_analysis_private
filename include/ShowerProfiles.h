#ifndef SHOWERPROFILES_H
#define SHOWERPROFILES_H

#include "BaseAnalysis.h"

#include "TObject.h"
#include "TH1F.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class ShowerProfiles: public BaseAnalysis {

 public:
  ShowerProfiles(double aEnergy, double aSf);
  ~ShowerProfiles();

  virtual void processEvent(podio::EventStore& store, int aEventId, bool verbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;

  /// Hit energy
  TH1F* h_hitEnergy;
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
  /// Long. profile in mm, direction of the shower from the generated particle
  TH1F* h_longProfile_particle_m;
  /// pt of the generated particle
  TH1F* h_ptGen;
  /// eta from hits' direction
  TH1F* h_eta;
  /// eta from particles' direction
  TH1F* h_eta_particle;

 private:
  void Initialize_histos();
  double m_energy;
  double m_sf;
  double SumE_hit;    // Total hit energy per event

  const double RcaloMin = 2850.; // 1950 ECAL
  const double RcaloThickness = 1908.; // !eta=0.36! 1800*1.06 800 ECAL
  const double layerThickness = 106; // mainly used as first layer thickness, 6 ECAL
  const double EtaMax = 10.0;
  const double X0 = 21.73; // in mm for HCAL, 15.586 old average X0 for 4 mm Lar + 2 mm Pb
  const double Lambda = 9; // average Lambda

};

#endif /* SHOWERPROFILES_H */
