#ifndef BASETWOFILEANALYSIS_H
#define BASETWOFILEANALYSIS_H

#include "TObject.h"
#include "TH1.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class BaseTwoFileAnalysis {
public:
  BaseTwoFileAnalysis();
  virtual ~BaseTwoFileAnalysis();

  void loop(const std::string& aFilenameSim ,const std::string& aFilenameRec, bool aVerbose = false);  //Open the file in the reader and loop through the events
protected:
  std::vector<TH1*> m_histograms;
private:
  virtual void processEvent(podio::EventStore& aStoreSim, podio::EventStore& aStoreRec, int aEventId, bool aVerbose = false) = 0;
  virtual void finishLoop(int aNumEvents, bool aVerbose) = 0;
  virtual void Delete_histos();
  virtual void Reset_histos();
};

#endif /* BASETWOFILEANALYSIS_H */
