// -*- C++ -*-
//
// Package:    RPC+CSCTrigger/CSCExtrapoltoRPC
// Class:      CSCExtrapoltoRPC
// 
/**\class CSCExtrapoltoRPC CSCExtrapoltoRPC.cc RPC+CSCTrigger/CSCExtrapoltoRPC/plugins/CSCExtrapoltoRPC.cc

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
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TH1.h"
#include "TTree.h"

#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeomManager.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeometry.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/RecSegment.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "RecoLocalMuon/RPCRecHit/src/CSCObjectMap.h"
#include "RecoLocalMuon/RPCRecHit/interface/CSCSegtoRPC.h"
#include "RecoLocalMuon/RPCRecHit/src/CSCStationIndex.h"

#include "L1Trigger/CSCTriggerPrimitives/src/CSCMotherboardME3141RPC.h"
#include <L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeometry.h>
#include <Geometry/RPCGeometry/interface/RPCRollSpecs.h>
#include <DataFormats/MuonDetId/interface/CSCTriggerNumbering.h>
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "L1Trigger/CSCTriggerPrimitives/src/CSCMotherboard.h"
#include <L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeomManager.h>

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

class CSCExtrapoltoRPC : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
//      explicit CSCExtrapoltoRPC(const edm::ParameterSet& iConfig, unsigned endcap, unsigned station, unsigned sector, unsigned subsector, unsigned chamber);
      explicit CSCExtrapoltoRPC(const edm::ParameterSet&)
;
      ~CSCExtrapoltoRPC();

      std::unique_ptr<RPCRecHitCollection> && thePoints(){ return std::move(_ThePoints); }

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      std::unique_ptr<RPCRecHitCollection> _ThePoints; 

      // ----------member data ---------------------------
      
      TTree *tree;
      TH1D *EventInfo;
    
      TH1D *Ndigis;
      TH1D *fNdigis;
      TH1D *bNdigis;
      TH1D *ME31NDigis;
      TH1D *ME41NDigis;

      unsigned int b_EVENT, b_RUN, b_LUMI;
      unsigned int b_numberofDigis;

      unsigned int b_fNDigis;
      unsigned int b_bNDigis;
      unsigned int b_ME31NDigis;
      unsigned int b_ME41NDigis;

      int b_Trknmb;
      int b_cscBX;
      int b_cscId;

      int b_CSCendcap;
      int b_CSCstation;
      int b_CSCsector;
      int b_CSCsubsector;
      int b_CSCstrip;
      int b_CSCkeyWire;

      edm::EDGetTokenT<MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>> corrlctsToken_;
      edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;
/*
      const unsigned theEndcap;
      const unsigned theStation;
      const unsigned theSector;
      const unsigned theSubsector;
      const unsigned theTrigChamber;
*/

};


//CSCExtrapoltoRPC::CSCExtrapoltoRPC(const edm::ParameterSet& iConfig, unsigned endcap, unsigned station, unsigned sector, unsigned subsector, unsigned chamber) : theEndcap(endcap), theStation(station), theSector(sector), theSubsector(subsector), theTrigChamber(chamber)
CSCExtrapoltoRPC::CSCExtrapoltoRPC(const edm::ParameterSet& iConfig)
{
   auto corrlctsDigiLabel = iConfig.getParameter<edm::InputTag>("simCSCTriggerpreDigis");
   corrlctsToken_ = consumes<MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>>(edm::InputTag(corrlctsDigiLabel.label(), "MPCSORTED" ));
   auto RPCDigiLabel = iConfig.getParameter<edm::InputTag>("simMuonRPCDigis");
   rpcRecHitsToken_ = consumes<RPCRecHitCollection>(edm::InputTag(RPCDigiLabel.label(), "" ));

   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("tree", "Tree for RPC+CSCTrigger");

   EventInfo = fs->make<TH1D>("EventInfo","Event Information",2,0,2);
   EventInfo->GetXaxis()->SetBinLabel(1,"Total Number of Events");
   EventInfo->GetXaxis()->SetBinLabel(2,"Selected Number of Events");

   Ndigis = fs->make<TH1D>("Ndigis", "", 10, 0, 10);
   Ndigis->GetXaxis()->SetTitle("Number of digis per chamber");
   Ndigis->GetYaxis()->SetTitle("Number of chamber");

   fNdigis = fs->make<TH1D>("fNdigis", "number of digis per chamber (forward)", 10, 0, 10);
   bNdigis = fs->make<TH1D>("bNdigis", "number of digis per chamber (backward)", 10, 0, 10);

   ME31NDigis = fs->make<TH1D>("ME31NDigis", "number of digis per chamber (ME3/1)", 10, 0, 10);
   ME31NDigis->GetXaxis()->SetTitle("Number of digis per chamber");
   ME31NDigis->GetYaxis()->SetTitle("Number of chamber");

   ME41NDigis = fs->make<TH1D>("ME41NDigis", "number of digis per chamber (ME4/1)", 10, 0, 10);
   ME41NDigis->GetXaxis()->SetTitle("Number of digis per chamber");
   ME41NDigis->GetYaxis()->SetTitle("Number of chamber");
}

CSCExtrapoltoRPC::~CSCExtrapoltoRPC()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//
// ------------ method called for each event  ------------
void
CSCExtrapoltoRPC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
   cout << "CSCDetId::minStationId(): " << CSCDetId::minStationId() << endl;

   cout << "maxTriggerCscId(): " << CSCTriggerNumbering::maxTriggerCscId() << endl;
   cout << "minTriggerCscId(): " << CSCTriggerNumbering::minTriggerCscId()  << endl;
   cout << "maxTriggerSectorId(): " << CSCTriggerNumbering::maxTriggerSectorId() << endl;
   cout << "minTriggerSectorId(): " << CSCTriggerNumbering::minTriggerSectorId() << endl;
   cout << "maxTriggerSubSectorId(): " << CSCTriggerNumbering::maxTriggerSubSectorId() << endl;
   cout << "minTriggerSubSectorId(): " << CSCTriggerNumbering::minTriggerSubSectorId() << endl;



   using namespace edm;

   EventInfo->Fill(0.5);

   edm::Handle<MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>> corrlcts;
   iEvent.getByToken(corrlctsToken_, corrlcts);

   edm::Handle<RPCRecHitCollection> rpcRecHits;
   iEvent.getByToken(rpcRecHitsToken_, rpcRecHits);

   if (!corrlcts.isValid()) {
     edm::LogInfo("DataNotFound") << "can't find CSCCorrleatedLCTDigiCollection with label "<< corrlcts << std::endl;
     return;
   }

   if (!rpcRecHits.isValid()) {
     edm::LogInfo("DataNotFound") << "can't find RPCRecHitDigiCollection with label "<< rpcRecHits << std::endl;
     return;
   }

   _ThePoints = std::make_unique<RPCRecHitCollection>();

   b_EVENT = b_RUN = b_LUMI = 0;

   b_Trknmb = 0;
   b_cscBX = 0;
   b_cscId = 0;

   b_CSCendcap = b_CSCstation = b_CSCsubsector = b_CSCsector = b_CSCstrip = b_CSCkeyWire = 0;

   b_EVENT  = iEvent.id().event();
   b_RUN    = iEvent.id().run();
   b_LUMI   = iEvent.id().luminosityBlock();

/*
    CSCTriggerGeomManager* geo_manager(CSCTriggerGeometry::get());
    const CSCChamber* cscChamber(geo_manager->chamber(theEndcap, theStation, theSector, theSubsector, theTrigChamber));
    const CSCLayer* keyLayer(cscChamber->layer(3));
    const CSCLayerGeometry* keyLayerGeometry(keyLayer->geometry());

    const LocalPoint lpCSC(keyLayerGeometry->topology()->localPosition(i));
*/

   std::cout << "New event\n" << endl;
   for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc=corrlcts.product()->begin(); csc!=corrlcts.product()->end(); csc++){  
       CSCCorrelatedLCTDigiCollection::Range range1 = corrlcts.product()->get((*csc).first);
//       cout << "\nwe are in ranage 1 " << endl;
       b_numberofDigis = b_fNDigis = b_bNDigis = b_ME31NDigis = b_ME41NDigis = 0;

       for(CSCCorrelatedLCTDigiCollection::const_iterator lct=range1.first; lct!=range1.second; lct++){

           const CSCDetId csc_id((*csc).first.rawId());
           cout << "DetId: " << csc_id << endl;

           b_CSCendcap = (*csc).first.endcap()-1;
           b_CSCstation = (*csc).first.station()-1;
           b_CSCsector  = (*csc).first.triggerSector()-1;
           b_CSCsubsector = CSCTriggerNumbering::triggerSubSectorFromLabels((*csc).first);
           b_CSCstrip = lct->getStrip();
           b_CSCkeyWire = lct->getKeyWG();
           
           b_cscId = lct->getCSCID();
	   b_Trknmb = lct->getTrknmb();
	   b_cscBX = lct->getBX();

	   cout << "getCSCID() = " << b_cscId << endl;
	   cout << "I'm here in CSC:: endcap: " << b_CSCendcap << " station: " << b_CSCstation << " sector: " << b_CSCsector << " subsector: " << b_CSCsubsector << " strip: " << b_CSCstrip << " wire: " << b_CSCkeyWire << endl;

	   cout << "I'm here in CSC:: endcap: " << csc_id.endcap() << " station: " << csc_id.station() << endl;


//           CSCTriggerGeomManager* geo_manager(CSCTriggerGeometry::get());
//           const CSCChamber* cscChamber(geo_manager->chamber(theEndcap, theStation, theSector, theSubsector, theTrigChamber));
           //for test
/*
           int newStation = 0;
           int newEndcap = 0;

           if (b_CSCstation == 0) newStation = 1;
           else newStation = 2;
                     if (b_CSCendcap == 0) newEndcap = 1;
           if (b_CSCendcap == 0) newEndcap = 1;
           else newEndcap = -1;
*/ 
//           const CSCChamber* cscChamber(geo_manager->chamber(b_CSCendcap, b_CSCstation, b_CSCsector, b_CSCsubsector, b_cscId));
//           const CSCLayer* keyLayer(cscChamber->layer(3));
//           const CSCLayerGeometry* keyLayerGeometry(keyLayer->geometry());
//           const LocalPoint lpCSC(keyLayerGeometry->topology()->localPosition(b_CSCstrip));

           //to check forward and backward endcap
           if ( b_CSCendcap == 0 ) b_fNDigis++; 
           if ( b_CSCendcap == 1 ) b_bNDigis++; 

           if ( b_CSCstation == 2 ){
             if ( b_cscId == 1 || b_cscId == 2 || b_cscId == 3 ) b_ME31NDigis++;
             }
           if ( b_CSCstation == 3 ){
             if ( b_cscId == 1 || b_cscId == 2 || b_cscId == 3 ) b_ME41NDigis++;
             }
       }

       cout << "number of digis: " << b_numberofDigis << endl;   
       cout << "ME31 NDigis: " << b_ME31NDigis << " ME41 NDigis: " << b_ME41NDigis << endl;   

       //fill histo
       Ndigis->Fill(b_numberofDigis);
       if (b_fNDigis != 0) fNdigis->Fill(b_fNDigis);
       if (b_bNDigis != 0) bNdigis->Fill(b_bNDigis);
       if (b_ME31NDigis != 0) ME31NDigis->Fill(b_ME31NDigis);
       if (b_ME41NDigis != 0) ME41NDigis->Fill(b_ME41NDigis);

   }

   tree->Fill();
   EventInfo->Fill(1.5);

}


// ------------ method called once each job just before starting event loop  ------------
void 
CSCExtrapoltoRPC::beginJob()
{

   tree->Branch("EVENT", &b_EVENT, "EVENT/i");
   tree->Branch("RUN"  , &b_RUN  , "RUN/i");
   tree->Branch("LUMI" , &b_LUMI , "LUMI/i");

   tree->Branch("Trknmb" , &b_Trknmb , "Trknmb/i");
   tree->Branch("cscBX" , &b_cscBX , "cscBX/i");
   tree->Branch("CSCID" , &b_cscId , "CSCID/i");
   tree->Branch("numberofDigis" , &b_numberofDigis , "numberofDigis/i");

   tree->Branch("CSCendcap" , &b_CSCendcap , "CSCendcap/i");
   tree->Branch("CSCstation" , &b_CSCstation , "CSCstation/i");
   tree->Branch("CSCsector" , &b_CSCsector , "CSCsector/i");
   tree->Branch("CSCsubsector" , &b_CSCsubsector , "CSCsubsector/i");
   tree->Branch("CSCstrip" , &b_CSCstrip , "CSCstrip/i");
   tree->Branch("CSCkeyWire" , &b_CSCkeyWire , "CSCkeyWire/i");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
CSCExtrapoltoRPC::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CSCExtrapoltoRPC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CSCExtrapoltoRPC);
