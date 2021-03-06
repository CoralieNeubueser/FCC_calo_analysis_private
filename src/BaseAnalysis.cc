#include "BaseAnalysis.h"

// podio specific includes
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

// STL
#include <vector>
#include <iostream>
#include <string>

BaseAnalysis::BaseAnalysis() {
  TH1::AddDirectory(kFALSE);
}

BaseAnalysis::~BaseAnalysis() {
  Delete_histos();
}

void BaseAnalysis::loop(const std::string& aFilename, bool aVerbose, bool aOneEvent) {
  //Reset histograms
  Reset_histos();

  //Open file in the reader
  auto reader = podio::ROOTReader();
  auto store = podio::EventStore();
  try {
    reader.openFile(aFilename);
    std::cout << "BaseAnalysis opening file " << aFilename << std::endl;
  }
  catch(std::runtime_error& err) {
    std::cerr<<err.what()<<". Quitting."<<std::endl;
    exit(1);
  }
  store.setReader(&reader);

  //Loop over all events
  unsigned nEvents = reader.getEntries();
  if (aOneEvent)
    nEvents = 1;
  std::cout << "Number of events: " << nEvents << std::endl;
  for(unsigned i=0; i<nEvents; ++i) {
    //    std::cout<<"reading event "<<i<<std::endl;
    if(i>11) aVerbose = false;

    processEvent(store, i, aVerbose);

    store.clear();
    reader.endOfEvent();
  }

  finishLoop(nEvents, aVerbose);

  std::cout << "End of loop" << std::endl;
  reader.closeFile();
  return;
}

void BaseAnalysis::Reset_histos() {
  for (auto iterator=m_histograms.begin(); iterator<m_histograms.end(); iterator++) {
    (*iterator)->Reset();
    (*iterator)->Sumw2();
  }
  for (auto iterator=m_profiles.begin(); iterator<m_profiles.end(); iterator++) {
    (*iterator)->Reset();
    (*iterator)->Sumw2();
  }
}

void BaseAnalysis::Delete_histos() {
  for (auto iterator=m_histograms.begin(); iterator<m_histograms.end(); iterator++) {
    delete (*iterator);
  }
  for (auto iterator=m_profiles.begin(); iterator<m_profiles.end(); iterator++) {
    delete (*iterator);
  }
  for (auto iterator=m_stacks.begin(); iterator<m_stacks.end(); iterator++) {
    delete (*iterator);
  }
  m_histograms.clear();
  m_profiles.clear();
  m_stacks.clear();
}
