// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//classes included to extract electron information
#include "DataFormats/PhotonReco/interface/Photon.h"
#include "DataFormats/PhotonReco/interface/PhotonFwd.h" 

//classes included to extract tracking for the electrons
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//additional classes for storage, containers and operations
#include<vector>
#include<string>
#include "TFile.h"
#include "TTree.h"
#include <stdlib.h>
#include<iostream>
#include<fstream>
#include<sstream>



//
// class declaration
//

class PhotonObjectInfoExtractorToCsv : public edm::EDAnalyzer {
   public:
      explicit PhotonObjectInfoExtractorToCsv(const edm::ParameterSet&);
      ~PhotonObjectInfoExtractorToCsv();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
 
 //declare a function to do the photon analysis
      void analyzePhotons(const edm::Event& iEvent, const edm::Handle<reco::PhotonCollection> &photons);
  //function to store info in csv
  void dumpPhotonsToCsv();
  //declare the input tag for the photons collection to be used (read from cofiguration)
  edm::InputTag photonsInput;
  int maxNumObjt;

  //Declare some variables for storage
  std::ofstream myfile;
  int maxpart;
  std::ostringstream oss;
  std::string theHeader;


  //and declare variable that will go into the root tree
  int runno; //run number
  int evtno; //event number
  int nphoton;//number of photons in the event
  std::string photon_partype; //type of particle
  std::vector<float> photon_e;
  std::vector<float> photon_pt;
  std::vector<float> photon_px;
  std::vector<float> photon_py;
  std::vector<float> photon_pz;
  std::vector<float> photon_eta;
  std::vector<float> photon_phi;
  std::vector<float> photon_ch;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PhotonObjectInfoExtractorToCsv::PhotonObjectInfoExtractorToCsv(const edm::ParameterSet& iConfig)

{
  //This should match the configuration in the corresponding python file
  photonsInput = iConfig.getParameter<edm::InputTag>("InputCollection");
  maxNumObjt = iConfig.getUntrackedParameter<int>("maxNumberPhotons",5);

}


PhotonObjectInfoExtractorToCsv::~PhotonObjectInfoExtractorToCsv()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhotonObjectInfoExtractorToCsv::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //get the global information first
   runno = iEvent.id().run();
   evtno  = iEvent.id().event();

   //Declare a container (or handle) where to store your photons.
   //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideDataFormatRecoPhoton
   Handle<reco::PhotonCollection> myphotons;

   //This is where the information gets extracted from the EDM file
   //Essentially, this corresponds to the information stored in a specific
   //branch within the EDM files.  For example, for recoPhotons one could get
   //just "photons", which would be the most used reco photons collection,
   //but also "photonsFromCosmics", which could be used to study fakes.
   //If you explore the branches in the EDM file with a ROOT TBrowser, 
   //you can indeed find a
   //"recoPhotons_photons__RECO" branch and a "recoPhotons_photonsFromCosmics__RECO" one.
   //Therefore, following the documentation
   //(https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMGetDataFromEvent?rev=20),
   //one could simply write "photons" instead of 
   //the photonsInput variable, which is extracted from
   //the configuration above.  However, using such a configuration variable
   //allows one to access a different branch or type of photons, some cosmic photons
   //for example, without having to recompile the code.
   iEvent.getByLabel(photonsInput, myphotons); 

   //Now, to keep it orderly, pass the collection to a subroutine that extracts
   //some of  the photon information
   //We do need to pass the event.  We could have also passed
   //the event setup if it were needed.
   analyzePhotons(iEvent,myphotons);
   dumpPhotonsToCsv();
   return;

}

// ------------ function to analyze photons
void 
PhotonObjectInfoExtractorToCsv::analyzePhotons(const edm::Event& iEvent, const edm::Handle<reco::PhotonCollection> &photons)
{
  //clear the storage containers for this objects in this event
  nphoton=0;
  photon_e.clear();
  photon_px.clear();
  photon_py.clear();
  photon_pz.clear();
  photon_pt.clear();
  photon_eta.clear();
  photon_phi.clear();
  photon_ch.clear();

  //check if the collection is valid
  if(photons.isValid()){
    //get the number of photons in the event
	//loop over all the photons in this event
    int idx = 0;
	for (reco::PhotonCollection::const_iterator recoPhoton = photons->begin(); recoPhoton!=photons->end(); ++recoPhoton){
      //find only globlal photons for this specific example
      //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPhotonAnalysis?rev=88
	  if(recoPhoton->isGlobalPhoton()) {
	    photon_partype = "G"; 
	    photon_e.push_back(recoPhoton->energy());
	    photon_pt.push_back(recoPhoton->pt());
	    photon_px.push_back(recoPhoton->px());
	    photon_py.push_back(recoPhoton->py());
	    photon_pz.push_back(recoPhoton->pz());
	    photon_eta.push_back(recoPhoton->eta());
	    photon_phi.push_back(recoPhoton->phi());
	    photon_ch.push_back(recoPhoton->charge());
	    ++idx;
	  }
	}
	nphoton = idx;
  }
  
}

// ------------ function to analyze photon
void PhotonObjectInfoExtractorToCsv::dumpPhotonsToCsv()
{
  unsigned int maxnumobjt = maxNumObjt;
  if(nphoton>0){
  oss.str("");oss.clear();oss<<runno;
  myfile<<oss.str();
  oss.str("");oss.clear();oss<<evtno;
  myfile<<","<<oss.str();
    for (unsigned int j=0;j<maxnumobjt;j++){
      oss.str("");oss.clear();oss<<photon_partype;
      myfile<<","<<oss.str();
      oss.str("");oss.clear();oss<<photon_e[j];
      j<photon_e.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      //      std::cout<<maxnumobjt<<"\t"<<nmu<<"\t"<<photon_e.size()<<"\t"<<photon_e[j]<<"\t"<<oss.str()<<std::endl;
      oss.str("");oss.clear();oss<<photon_px[j];
      j<photon_px.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<photon_py[j];
      j<photon_py.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<photon_pz[j];
      j<photon_pz.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<photon_pt[j];
      j<photon_pt.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<photon_eta[j];
      j<photon_eta.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<photon_phi[j];
      j<photon_phi.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<photon_ch[j];
      j<photon_ch.size() ? myfile<<","<<oss.str():myfile<<",0.0";
    }
  myfile<<"\n";
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
PhotonObjectInfoExtractorToCsv::beginJob()
{
  //Define storage
  myfile.open("PhotonObjectInfo.csv");
  //Write the header.
  //create the header string accordingly
  theHeader = "Run,Event";
  for(int j =1;j<maxNumObjt+1;j++){
    oss.str(""); oss<<j;
    std::string idxstr = oss.str();
    theHeader += ",type"+idxstr+",E"+idxstr+",px"+idxstr+",py"+idxstr+",pz"+idxstr+",pt"+idxstr+",eta"+idxstr+",phi"+idxstr+",Q"+idxstr;
  }
  
  myfile<<theHeader<<"\n";

}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonObjectInfoExtractorToCsv::endJob() 
{

  //save file
  myfile.close();

}

// ------------ method called when starting to processes a run  ------------
void 
PhotonObjectInfoExtractorToCsv::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PhotonObjectInfoExtractorToCsv::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PhotonObjectInfoExtractorToCsv::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PhotonObjectInfoExtractorToCsv::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonObjectInfoExtractorToCsv::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonObjectInfoExtractorToCsv);std::ostringstream oss;
