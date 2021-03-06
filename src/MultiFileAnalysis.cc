#include "MultiFileAnalysis.h"

// podio specific includes
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

// STL
#include <vector>
#include <iostream>
#include <string>

MultiFileAnalysis::MultiFileAnalysis() {
  TH1::AddDirectory(kFALSE);
}

MultiFileAnalysis::~MultiFileAnalysis() {
  Delete_histos();
}

void MultiFileAnalysis::loop(const int numFiles, const char** aFilenames, const double *aEnergies, bool aVerbose) {
  //Reset histograms
  Reset_histos();
  unsigned totEvents =0;  
  auto store = podio::EventStore();

  std::cout << "Number of input files : " << numFiles << std::endl;
  for (uint ifile=0; ifile<numFiles; ifile++) {
    //Open file in the reader
    unsigned nEvents =0;  
    auto reader = podio::ROOTReader();
    auto store = podio::EventStore();
    try {
      reader.openFile(aFilenames[ifile]);
      std::cout << "MultiFileAnalysis opening file " << aFilenames[ifile] << std::endl;
      std::cout << "Energy                         " << aEnergies[ifile] << std::endl;
    }
    catch(std::runtime_error& err) {
      std::cerr<<err.what()<<". Quitting."<<std::endl;
      exit(1);
    }
    store.setReader(&reader);
    
    //Loop over all events
    nEvents += reader.getEntries();
    totEvents += reader.getEntries();
    std::cout << "Number of events: " << nEvents << std::endl;
    for(unsigned i=0; i<nEvents; ++i) {
      //    std::cout<<"reading event "<<i<<std::endl;
      if(i>11) aVerbose = false;
      
      processEvent(store, i, aEnergies[ifile], aVerbose);

      store.clear();
      reader.endOfEvent();
    }
  }
  std::cout << "Number of total events: " << totEvents << std::endl;
  
  finishLoop(totEvents, aVerbose);
  
  std::cout << "End of loop" << std::endl;
  return;
}

void MultiFileAnalysis::Reset_histos() {
  for (auto iterator=m_histograms.begin(); iterator<m_histograms.end(); iterator++) {
    (*iterator)->Reset();
    (*iterator)->Sumw2();
  }
  for (auto iterator=m_profiles.begin(); iterator<m_profiles.end(); iterator++) {
    (*iterator)->Reset();
    (*iterator)->Sumw2();
  }
}

void MultiFileAnalysis::Delete_histos() {
  for (auto iterator=m_histograms.begin(); iterator<m_histograms.end(); iterator++) {
    delete (*iterator);
  }
  for (auto iterator=m_profiles.begin(); iterator<m_profiles.end(); iterator++) {
    delete (*iterator);
  }
  m_histograms.clear();
  m_profiles.clear();
}
