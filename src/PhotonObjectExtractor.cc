#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//classes included to extract electron information
#include "DataFormats/ElectronReco/interface/Electron.h"
#include "DataFormats/ElectronReco/interface/ElectronFwd.h" 

//classes included to extract tracking for the electrons
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

class ElectronObjectInfoExtractor : public edm::EDAnalyzer {
   public:
      explicit ElectronObjectInfoExtractor(const edm::ParameterSet&);
      ~ElectronObjectInfoExtractor();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
 
 //declare a function to do the electron analysis
      void analyzeElectrons(const edm::Event& iEvent);
  //declare the input tag for the electrons collection to be used (read from cofiguration)
  edm::InputTag electronsInput;
  
  //These variables will be global

  //Declare some variables for storage
  TFile* myfile;//root file
  TTree* mytree;//root tree

  //and declare variable that will go into the root tree
  int runno; //run number
  int evtno; //event number
  int nelectron; //number of electrons in the event
  std::vector<float> electron_e;
  std::vector<float> electron_pt;
  std::vector<float> electron_px;
  std::vector<float> electron_py;
  std::vector<float> electron_pz;
  std::vector<float> electron_eta;
  std::vector<float> electron_phi;
  std::vector<float> electron_ch;

  

  
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
ElectronObjectInfoExtractor::ElectronObjectInfoExtractor(const edm::ParameterSet& iConfig)

{
  //This should match the configuration in the corresponding python file
  electronsInput = iConfig.getParameter<edm::InputTag>("InputCollection");

}


ElectronObjectInfoExtractor::~ElectronObjectInfoExtractor()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronObjectInfoExtractor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //get the global information first
   runno = iEvent.id().run();
   evtno  = iEvent.id().event();
   
   //Now, to keep it orderly, pass the collection to a subroutine that extracts
   //some of  the electron information
   //We do need to pass the event.  We could have also passed
   //the event setup if it were needed.  For example, if you need to 
   //store the trigger information you will need to follow
   //this example (https://github.com/cms-opendata-analyses/trigger_examples/tree/master/TriggerInfo/TriggerInfoAnalyzer) and check how to do it.
   analyzeElectrons(iEvent);

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

// ------------ function to analyze electrons
void 
ElectronObjectInfoExtractor::analyzeElectrons(const edm::Event& iEvent)
{
  //clear the storage containers for this objects in this event
  //these were declared above and are global
  nelectron=0;
  electron_e.clear();
  electron_pt.clear();
  electron_px.clear();
  electron_py.clear();
  electron_pz.clear();
  electron_eta.clear();
  electron_phi.clear();
  electron_ch.clear()

  edm::Handle<reco::ElectronCollection> myelectrons;

   //This is where the information gets extracted from the EDM file
   //Essentially, this corresponds to the information stored in a specific
   //branch within the EDM files.  For example, for recoElectrons one could get
   //just "electrons", which would be the most used reco electrons collection,
   //but also "electronsFromCosmics", which could be used to study fakes.
   //If you explore the branches in the EDM file with a ROOT TBrowser, 
   //you can indeed find a
   //"recoElectrons_electrons__RECO" branch and a "recoElectrons_electronsFromCosmics__RECO" one.
   //Therefore, following the documentation
   //(https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMGetDataFromEvent?rev=20),
   //one could simply write "electrons" instead of 
   //the electronsInput variable, which is extracted from
   //the configuration above.  However, using such a configuration variable
   //allows one to access a different branch or type of electrons, some cosmic electrons
   //for example, without having to recompile the code.
   iEvent.getByLabel(electronsInput, myelectrons); 


  //check if the collection is valid
  if(myelectrons.isValid()){
      //get the number of electrons in the event
      nelectrons=(*myelectrons).size();
	//loop over all the electrons in this event
	for (reco::ElectronCollection::const_iterator recoElectron = myelectrons->begin(); recoElectron!=myelectrons->end(); ++recoElectron){
      //find only globlal electrons for this specific example
      //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookElectronAnalysis?rev=88
	  //Note that this would be already a selection cut, i.e.
	  //requiring it to be global is a constrain on what kind of electron it is
      if(recoElectron->isGlobalElectron()) {
	  electron_e.push_back(recoElectron->energy());
	  electron_pt.push_back(recoElectron->pt());
	  electron_px.push_back(recoElectron->px());
	  electron_py.push_back(recoElectron->py());
	  electron_pz.push_back(recoElectron->pz());
	  electron_eta.push_back(recoElectron->eta());
	  electron_phi.push_back(recoElectron->phi());
	  electron_ch.push_back(recoElectron->ch());

	  // get the track combinig the information from both the Tracker and the Spectrometer
	  reco::TrackRef recoCombinedGlbTrack = recoElectron->combinedElectron();
//	  electron_glbtrk_pt.push_back(recoCombinedGlbTrack->pt());
//  electron_glbtrk_eta.push_back(recoCombinedGlbTrack->eta());
//  electron_glbtrk_phi.push_back(recoCombinedGlbTrack->phi());

	  //here one could apply some identification
	  //cuts to show how to do particle id, and store
	  //refined variables
	}
      else{
	//Here I put default values for those electrons that are not global
	//so the containers do not show up as empty. One could do
	//this in a smarter way though.
	electron_e.push_back(-999);
	electron_pt.push_back(-999);
	electron_px.push_back(-999);
	electron_py.push_back(-999);
	electron_pz.push_back(-999);
	electron_eta.push_back(-999);
	electron_phi.push_back(-999);
	electron_ch.push_back(-999);
	//electron_glbtrk_pt.push_back(-999);
	//electron_glbtrk_eta.push_back(-999);
	//electron_glbtrk_phi.push_back(-999);
      }
      }
    }
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
ElectronObjectInfoExtractor::beginJob()
{
  //Define storage variables
  myfile = new TFile("ElectronObjectInfo.root","RECREATE");
  mytree = new TTree("mytree","Rootuple with object information");
  //point root branches to the right place
  //this is a typical ROOT way of doing it
  //one could of course try to store the information in a different format
  //for exmaple json (nested) format or plain csv 
  //(that would be something nice to implement).
  mytree->Branch("runno",&runno,"runno/I");
  mytree->Branch("evtno",&evtno,"evtno/I");
  mytree->Branch("nelectron",&nelectron,"nelectron/I");
  mytree->Branch("electron_e",&electron_e);
  mytree->Branch("electron_pt",&electron_pt);
  mytree->Branch("electron_px",&electron_px);
  mytree->Branch("electron_py",&electron_py);
  mytree->Branch("electron_pz",&electron_pz);
  mytree->Branch("electron_eta",&electron_eta);
  mytree->Branch("electron_phi",&electron_phi);
  mytree->Branch("electron_ch",&electron_ch);
  //mytree->Branch("electron_glbtrk_pt",&electron_glbtrk_pt);
  //mytree->Branch("electron_glbtrk_eta",&electron_glbtrk_eta);
  //mytree->Branch("electron_glbtrk_phi",&electron_glbtrk_phi);


  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronObjectInfoExtractor::endJob() 
{

  //save file
  myfile->Write();

}

// ------------ method called when starting to processes a run  ------------
void 
ElectronObjectInfoExtractor::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ElectronObjectInfoExtractor::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ElectronObjectInfoExtractor::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ElectronObjectInfoExtractor::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronObjectInfoExtractor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronObjectInfoExtractor);
