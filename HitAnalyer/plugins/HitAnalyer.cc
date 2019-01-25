// -*- C++ -*-
//
// Package:    RPC+CSCTrigger/HitAnalyer
// Class:      HitAnalyer
// 
/**\class HitAnalyer HitAnalyer.cc RPC+CSCTrigger/HitAnalyer/plugins/HitAnalyer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ms. Choi Jieun
//         Created:  Fri, 18 Jan 2019 04:00:35 GMT
//
//


// system include files
#include <memory>

#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "TH1.h"
#include "TTree.h"


//CSCCorrelatedLCTDigi
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using namespace edm;

class HitAnalyer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit HitAnalyer(const edm::ParameterSet&);
      ~HitAnalyer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
      TTree *tree;
      TH1D *EventInfo;

      unsigned int b_EVENT, b_RUN, b_LUMI;

      edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> corrlctsToken_;
//
//       int b_NCSCCorrLCTDigi_ME31;
//      int b_NCSCCorrLCTDigi_ME41;

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
HitAnalyer::HitAnalyer(const edm::ParameterSet& iConfig)

{
//   MuonDigiRPC_( consumes<CSCCorrelatedLCTDigiCollection>( iConfig.getParameter<edm::InputTag> ( "MuonDigiRPC" )));
//   corrlctsToken_ = consumes<CSCCorrelatedLCTDigiCollection>(corrlctsTag_);
//   corrlctsToken_ = consumes<CSCCorrelatedLCTDigiCollection>(iConfig.getParameter<edm::InputTag> ( "MuonDigiRPC"));
  auto MuonDigiRPCLabel = iConfig.getParameter<edm::InputTag>("MuonDigiRPC");
   corrlctsToken_ = consumes<CSCCorrelatedLCTDigiCollection>(edm::InputTag(MuonDigiRPCLabel.label(), "MPCSORTED" ));
//   corrlctsToken_ = consumes<CSCCorrelatedLCTDigiCollection>(iConfig.getParameter<edm::InputTag> ( "MuonDigiRPC" ));
//   MuonDigiCSC_( consumes< somethig >( iConfig.getParameter<edm::InputTag> ( "MuonDigiCSC" )));

   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("tree", "Tree for RPC+CSCTrigger");

   EventInfo = fs->make<TH1D>("EventInfo","Event Information",2,0,2);
   EventInfo->GetXaxis()->SetBinLabel(1,"Total Number of Events");
   EventInfo->GetXaxis()->SetBinLabel(2,"Selected Number of Events");

}

HitAnalyer::~HitAnalyer()
{


   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
HitAnalyer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
//   using namespace csccorrelatedlctdigis;

   EventInfo->Fill(0.5);

   edm::Handle<CSCCorrelatedLCTDigiCollection> corrlcts;
   iEvent.getByToken(corrlctsToken_, corrlcts);


   for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc=corrlcts.product()->begin(); csc!=corrlcts.product()->end(); csc++)
     {
       CSCCorrelatedLCTDigiCollection::Range range1 = corrlcts.product()->get((*csc).first);
       for(CSCCorrelatedLCTDigiCollection::const_iterator lct=range1.first; lct!=range1.second; lct++)
         {
	   std::cout << lct->getTrknmb() << std::endl;
	}
	}

/*
   Handle<> ;
   iEvent.getByLabel("simMuonRPCDigis", );

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
*/

   b_EVENT = b_RUN = b_LUMI = 0;
   

   //Fill Branch

    b_EVENT  = iEvent.id().event();
    b_RUN    = iEvent.id().run();
    b_LUMI   = iEvent.id().luminosityBlock();

    tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
HitAnalyer::beginJob()
{

   tree->Branch("EVENT", &b_EVENT, "EVENT/i");
   tree->Branch("RUN"  , &b_RUN  , "RUN/i");
   tree->Branch("LUMI" , &b_LUMI , "LUMI/i");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HitAnalyer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HitAnalyer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HitAnalyer);
