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


#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"
//#include "RecoLocalMuon/CSCValidation/macros/triggerCscIdSector.h"

#include "L1Trigger/CSCTriggerPrimitives/src/CSCMotherboardME3141RPC.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeomManager.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeometry.h"

#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

//#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
//#include "Geometry/CSCGeometry/interface/CSCGeometry.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using namespace edm;
using namespace std;

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
    
      TH1D *Ndigis;

      unsigned int b_EVENT, b_RUN, b_LUMI;
      unsigned int b_numberofDigis;

      int b_Trknmb;
      int b_BX;
      int b_cscId;

      int b_endcap;
      int b_station;
      int b_sector;
      int b_subsector;
      int b_ring;
      int b_strip;
      int b_keyWire;

      int b_test;

//      edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> corrlctsToken_;
//      edm::EDGetTokenT<RPCDigiCollection> rpcDigiToken_;

      edm::EDGetTokenT<MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>> corrlctsToken_;
      edm::EDGetTokenT<MuonDigiCollection<RPCDetId,RPCDigi>> rpcDigiToken_;

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
//   corrlctsToken_ = consumes<CSCCorrelatedLCTDigiCollection>(iConfig.getParameter<edm::InputTag> ( "MuonDigiRPC"));

//   auto corrlctsDigiLabel = iConfig.getParameter<edm::InputTag>("simCSCTriggerpreDigis");
//   corrlctsToken_ = consumes<CSCCorrelatedLCTDigiCollection>(edm::InputTag(corrlctsDigiLabel.label(), "MPCSORTED" ));
//   auto RPCDigiLabel = iConfig.getParameter<edm::InputTag>("simMuonRPCDigis");
//   rpcDigiToken_ = consumes<RPCDigiCollection>(edm::InputTag(RPCDigiLabel.label(), "" ));

   auto corrlctsDigiLabel = iConfig.getParameter<edm::InputTag>("simCSCTriggerpreDigis");
   corrlctsToken_ = consumes<MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>>(edm::InputTag(corrlctsDigiLabel.label(), "MPCSORTED" ));
   auto RPCDigiLabel = iConfig.getParameter<edm::InputTag>("simMuonRPCDigis");
   rpcDigiToken_ = consumes<MuonDigiCollection<RPCDetId,RPCDigi>>(edm::InputTag(RPCDigiLabel.label(), "" ));

//   cscSegments = consumes<CSCSegmentCollection>(iConfig.getUntrackedParameter < edm::InputTag > ("cscSegments"));
//   tracks = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter < edm::InputTag > ("tracks"))

   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("tree", "Tree for RPC+CSCTrigger");

   EventInfo = fs->make<TH1D>("EventInfo","Event Information",2,0,2);
   EventInfo->GetXaxis()->SetBinLabel(1,"Total Number of Events");
   EventInfo->GetXaxis()->SetBinLabel(2,"Selected Number of Events");

   Ndigis = fs->make<TH1D>("Ndigis", "number of digis per chamber", 10, 0, 10);
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

//   edm::Handle<CSCCorrelatedLCTDigiCollection> corrlcts;
//   iEvent.getByToken(corrlctsToken_, corrlcts);
//   edm::Handle<RPCDigiCollection> rpcdigis;
//   iEvent.getByToken(rpcDigiToken_, rpcdigis);

   edm::Handle<MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>> corrlcts;
   iEvent.getByToken(corrlctsToken_, corrlcts);
   edm::Handle<MuonDigiCollection<RPCDetId,RPCDigi>> rpcdigis;
   iEvent.getByToken(rpcDigiToken_, rpcdigis);


   //csctrack
/*
   edm::Handle<CSCSegmentCollection> allCSCSegments;
   iEvent.getByToken(cscSegments, allCSCSegments);
   
   edm::Handle<reco::TrackCollection> alltracks;
   iEvent.getByToken(tracks,alltracks);

   reco::TrackCollection::const_iterator Track;
   CSCSegmentCollection::const_iterator segmentCSC;
   CSCSegment segmentCSCToStore;

   auto selectedCscSegments = std::make_unique<CSCSegmentCollection>();
   std::vector<CSCDetId> chamberIdCSC;
  */ 

   if (!corrlcts.isValid()) {
     edm::LogInfo("DataNotFound") << "can't find CSCCorrleatedLCTDigiCollection with label "<< corrlcts << std::endl;
     return;
   }

   if (!rpcdigis.isValid()) {
     edm::LogInfo("DataNotFound") << "can't find RPCDigiCollection with label "<< rpcdigis << std::endl;
     return;
   }

   b_EVENT = b_RUN = b_LUMI = 0;

   b_Trknmb = 0;
   b_BX = 0;
   b_cscId = 0;
   b_numberofDigis = 0;

   b_endcap = b_station = b_subsector = b_sector = b_ring = b_strip = b_keyWire = 0;
   b_test = 0;

   b_EVENT  = iEvent.id().event();
   b_RUN    = iEvent.id().run();
   b_LUMI   = iEvent.id().luminosityBlock();
 
   std::cout << "New event\n" << endl;
   for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc=corrlcts.product()->begin(); csc!=corrlcts.product()->end(); csc++)

      {  
       CSCCorrelatedLCTDigiCollection::Range range1 = corrlcts.product()->get((*csc).first);
       cout << "\nwe are in ranage 1 " << endl;
       b_test = 0;
       b_numberofDigis = 0;
       for(CSCCorrelatedLCTDigiCollection::const_iterator lct=range1.first; lct!=range1.second; lct++)
         {

           if ( !rpcdigis.isValid() ) continue;
           if ( !corrlcts.isValid() ) continue;

           // cout << "b_test = " << b_test << endl;

           b_endcap = (*csc).first.endcap()-1;
           b_station = (*csc).first.station()-1;
           b_sector  = (*csc).first.triggerSector()-1;
           b_subsector = CSCTriggerNumbering::triggerSubSectorFromLabels((*csc).first);
           b_strip = lct->getStrip();
           b_keyWire = lct->getKeyWG();
           
           b_cscId = lct->getCSCID();
	   b_Trknmb = lct->getTrknmb();
	   b_BX = lct->getBX();
           b_numberofDigis++;
           b_test++;
           Ndigis->Fill(b_numberofDigis);

	   cout << "getCSCID() = " << b_cscId << endl;
	   cout << "I'm here:: endcap: " << b_endcap << " station: " << b_station << " sector: " << b_sector << " subsector: " << b_subsector << " strip: " << b_strip << " wire: " << b_keyWire << endl;
       
           cout << "number of digis: " << b_numberofDigis << endl;

	 }
       }


   tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
HitAnalyer::beginJob()
{

   tree->Branch("EVENT", &b_EVENT, "EVENT/i");
   tree->Branch("RUN"  , &b_RUN  , "RUN/i");
   tree->Branch("LUMI" , &b_LUMI , "LUMI/i");

   tree->Branch("Trknmb" , &b_Trknmb , "Trknmb/i");
   tree->Branch("BX" , &b_BX , "BX/i");
   tree->Branch("CSCID" , &b_cscId , "CSCID/i");
   tree->Branch("numberofDigis" , &b_numberofDigis , "numberofDigis/i");

   tree->Branch("endcap" , &b_endcap , "endcap/i");
   tree->Branch("station" , &b_station , "station/i");
   tree->Branch("sector" , &b_sector , "sector/i");
   tree->Branch("subsector" , &b_subsector , "subsector/i");
   tree->Branch("strip" , &b_strip , "strip/i");
   tree->Branch("keyWire" , &b_keyWire , "keyWire/i");
   tree->Branch("test" , &b_test , "test/i");

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
