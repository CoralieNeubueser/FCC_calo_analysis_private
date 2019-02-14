#include "TreeAnalysis.h"

// podio specific includes
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

// STL
#include <vector>
#include <iostream>
#include <string>

TreeAnalysis::TreeAnalysis() {
}

TreeAnalysis::~TreeAnalysis() {
}

void TreeAnalysis::loop(const std::string& aFilename, bool aVerbose) {
  //Open file in the reader
  auto reader = podio::ROOTReader();
  auto store = podio::EventStore();
  try {
    reader.openFile(aFilename);
    std::cout << "TreeAnalysis opening file " << aFilename << std::endl;
  }
  catch(std::runtime_error& err) {
    std::cerr<<err.what()<<". Quitting."<<std::endl;
    exit(1);
  }
  store.setReader(&reader);

  //Loop over all events
  unsigned nEvents = reader.getEntries();
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
  return;
}

