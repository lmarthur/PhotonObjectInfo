#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//classes included to extract photon information
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

//classes included to extract tracking for the photons
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//additional classes for storage, containers and operations
#include<vector>
#include<string>
#include "TFile.h"
#include "TTree.h"
#include <stdlib.h>



//
// class declaration
//

class PhotonObjectInfoExtractor : public edm::EDAnalyzer {
   public:
      explicit PhotonObjectInfoExtractor(const edm::ParameterSet&);
      ~PhotonObjectInfoExtractor();

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
      void analyzePhotons(const edm::Event& iEvent);
  //declare the input tag for the photons collection to be used (read from cofiguration)
  edm::InputTag photonsInput;
  
  //These variables will be global

  //Declare some variables for storage
  TFile* myfile;//root file
  TTree* mytree;//root tree

  //and declare variable that will go into the root tree
  int runno; //run number
  int evtno; //event number
  int nphotons; //number of photons in the event
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
PhotonObjectInfoExtractor::PhotonObjectInfoExtractor(const edm::ParameterSet& iConfig)

{
  //This should match the configuration in the corresponding python file
  photonsInput = iConfig.getParameter<edm::InputTag>("InputCollection");

}


PhotonObjectInfoExtractor::~PhotonObjectInfoExtractor()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhotonObjectInfoExtractor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //get the global information first
   runno = iEvent.id().run();
   evtno  = iEvent.id().event();
   
   //Now, to keep it orderly, pass the collection to a subroutine that extracts
   //some of  the photon information
   //We do need to pass the event.  We could have also passed
   //the event setup if it were needed.  For example, if you need to 
   //store the trigger information you will need to follow
   //this example (https://github.com/cms-opendata-analyses/trigger_examples/tree/master/TriggerInfo/TriggerInfoAnalyzer) and check how to do it.
   analyzePhotons(iEvent);

   //Here, if one were to write a more general PhysicsObjectsInfoExtractor.cc
   //code, this is where the rest of the objects extraction will be, for exmaple:

   //analyzeElectrons(iEvent);
   //analyzeJets(iEvent, iSetup);
   //analyzeMet (iEvent, iSetup);
   //......

  

   //fill the root tree
   mytree->Fill();
   return;

}

// ------------ function to analyze photons
void 
PhotonObjectInfoExtractor::analyzePhotons(const edm::Event& iEvent)
{
  //clear the storage containers for this objects in this event
  //these were declared above and are global
  nphoton=0;
  photon_e.clear();
  photon_pt.clear();
  photon_px.clear();
  photon_py.clear();
  photon_pz.clear();
  photon_eta.clear();
  photon_phi.clear();
  photon_ch.clear()

  edm::Handle<reco::PhotonCollection> myphotons;

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
   iEvent.getByLabel(photonsInput, myphoton); 


  //check if the collection is valid
  if(myphoton.isValid()){
      //get the number of photons in the event
      nphotons=(*myphoton).size();
	//loop over all the photons in this event
	for (reco::PhotonCollection::const_iterator recoPhoton = myphotons->begin(); recoPhoton!=myphotons->end(); ++recoPhoton){
      //find only globlal photons for this specific example
      //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPhotonAnalysis?rev=88
	  //Note that this would be already a selection cut, i.e.
	  //requiring it to be global is a constrain on what kind of photon it is
      if(recoPhoton->isGlobalPhoton()) {
	  photon_e.push_back(recoPhoton->energy());
	  photon_pt.push_back(recoPhoton->pt());
	  photon_px.push_back(recoPhoton->px());
	  photon_py.push_back(recoPhoton->py());
	  photon_pz.push_back(recoPhoton->pz());
	  photon_eta.push_back(recoPhoton->eta());
	  photon_phi.push_back(recoPhoton->phi());
	  photon_ch.push_back(recoPhoton->ch());

	  // get the track combinig the information from both the Tracker and the Spectrometer
	  reco::TrackRef recoCombinedGlbTrack = recoPhoton->combinedPhoton();
//	  photon_glbtrk_pt.push_back(recoCombinedGlbTrack->pt());
//  photon_glbtrk_eta.push_back(recoCombinedGlbTrack->eta());
//  photon_glbtrk_phi.push_back(recoCombinedGlbTrack->phi());

	  //here one could apply some identification
	  //cuts to show how to do particle id, and store
	  //refined variables
	}
      else{
	//Here I put default values for those photons that are not global
	//so the containers do not show up as empty. One could do
	//this in a smarter way though.
	photon_e.push_back(-999);
	photon_pt.push_back(-999);
	photon_px.push_back(-999);
	photon_py.push_back(-999);
	photon_pz.push_back(-999);
	photon_eta.push_back(-999);
	photon_phi.push_back(-999);
	photon_ch.push_back(-999);
	//photon_glbtrk_pt.push_back(-999);
	//photon_glbtrk_eta.push_back(-999);
	//photon_glbtrk_phi.push_back(-999);
      }
      }
    }
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
PhotonObjectInfoExtractor::beginJob()
{
  //Define storage variables
  myfile = new TFile("PhotonObjectInfo.root","RECREATE");
  mytree = new TTree("mytree","Rootuple with object information");
  //point root branches to the right place
  //this is a typical ROOT way of doing it
  //one could of course try to store the information in a different format
  //for exmaple json (nested) format or plain csv 
  //(that would be something nice to implement).
  mytree->Branch("runno",&runno,"runno/I");
  mytree->Branch("evtno",&evtno,"evtno/I");
  mytree->Branch("nphoton",&nphoton,"nphoton/I");
  mytree->Branch("photon_e",&photon_e);
  mytree->Branch("photon_pt",&photon_pt);
  mytree->Branch("photon_px",&photon_px);
  mytree->Branch("photon_py",&photon_py);
  mytree->Branch("photon_pz",&photon_pz);
  mytree->Branch("photon_eta",&photon_eta);
  mytree->Branch("photon_phi",&photon_phi);
  mytree->Branch("photon_ch",&photon_ch);
  //mytree->Branch("photon_glbtrk_pt",&photon_glbtrk_pt);
  //mytree->Branch("photon_glbtrk_eta",&photon_glbtrk_eta);
  //mytree->Branch("photon_glbtrk_phi",&photon_glbtrk_phi);

  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonObjectInfoExtractor::endJob() 
{

  //save file
  myfile->Write();

}

// ------------ method called when starting to processes a run  ------------
void 
PhotonObjectInfoExtractor::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PhotonObjectInfoExtractor::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PhotonObjectInfoExtractor::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PhotonObjectInfoExtractor::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonObjectInfoExtractor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonObjectInfoExtractor);
