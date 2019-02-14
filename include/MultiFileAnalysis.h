#ifndef MULTIFILEANALYSIS_H
#define MULTIFILEANALYSIS_H

#include "TObject.h"
#include "TH1.h"
#include "TProfile.h"
#include <string>

namespace podio {
  class EventStore;
  class ROOTReader;
}

class MultiFileAnalysis {
public:
  MultiFileAnalysis();
  virtual ~MultiFileAnalysis();

  void loop(const int numFiles, const char** aFilenames, const double *aEnergies, bool aVerbose = false);  //Open the files in the reader and loop through all events
 protected:
  std::vector<TH1*> m_histograms;
  std::vector<TProfile*> m_profiles;
private:
  virtual void processEvent(podio::EventStore& aStore, int aEventId, double aEnergy, bool aVerbose = false) = 0;
  virtual void finishLoop(int aNumEvents, bool aVerbose) = 0;
  virtual void Delete_histos();
  virtual void Reset_histos();
};

#endif /* MULTIFILEANALYSIS_H */
