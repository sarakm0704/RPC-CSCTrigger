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
#include "Geometry/Records/interface/MuonGeometryRecord.h"
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
      unsigned int b_cscId;

      int b_CSCendcap;
      int b_CSCstation;
      unsigned int b_CSCsector;
      unsigned int b_CSCsubsector;
      int b_CSCstrip;
      int b_CSCkeyWire;

      edm::EDGetTokenT<MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>> corrlctsToken_;
      edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;

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

   edm::ESHandle<CSCGeometry> cscGeo;
   iSetup.get<MuonGeometryRecord>().get( cscGeo );     

   edm::ESHandle<RPCGeometry> rpcGeo;
   iSetup.get<MuonGeometryRecord>().get( rpcGeo );     
 
   cout << "CSCDetId::minRingId(): " << CSCDetId::minRingId() << endl;
   cout << "CSCDetId::maxRingId(): " << CSCDetId::maxRingId() << endl;

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
           b_numberofDigis++;

	   cout << "getCSCID() = " << b_cscId << endl;
	   cout << "I'm here in CSC:: endcap: " << b_CSCendcap << " station: " << b_CSCstation << " sector: " << b_CSCsector << " subsector: " << b_CSCsubsector << " strip: " << b_CSCstrip << " wire: " << b_CSCkeyWire << endl;

	   cout << "I'm here in CSC:: endcap: " << csc_id.endcap() << " station: " << csc_id.station() << endl;

           //CSCTriggerGeomManager* geo_manager(CSCTriggerGeometry::get());
           //const CSCChamber* cscChamber(geo_manager->chamber(newEndcap, newStation, b_CSCsector, b_CSCsubsector, b_cscId));
          
         //  const RPCChamber* rpcChamber=rpcGeo->chamber(rpc_id);
	  // auto randRoll(rpcChamber->roll(2));
                     
           const CSCChamber* cscChamber=cscGeo->chamber(csc_id);
           const CSCLayer* keyLayer1(cscChamber->layer(1));
           const CSCLayer* keyLayer2(cscChamber->layer(2));
           const CSCLayer* keyLayer3(cscChamber->layer(3));
           const CSCLayer* keyLayer4(cscChamber->layer(4));
           const CSCLayer* keyLayer5(cscChamber->layer(5));
           const CSCLayer* keyLayer6(cscChamber->layer(6));
           const CSCLayerGeometry* keyLayerGeometry1(keyLayer1->geometry());
           const CSCLayerGeometry* keyLayerGeometry2(keyLayer2->geometry());
           const CSCLayerGeometry* keyLayerGeometry3(keyLayer3->geometry());
           const CSCLayerGeometry* keyLayerGeometry4(keyLayer4->geometry());
           const CSCLayerGeometry* keyLayerGeometry5(keyLayer5->geometry());
           const CSCLayerGeometry* keyLayerGeometry6(keyLayer6->geometry());
           const LocalPoint lpCSC_strip1(keyLayerGeometry1->topology()->localPosition(b_CSCstrip));
           const LocalPoint lpCSC_strip2(keyLayerGeometry2->topology()->localPosition(b_CSCstrip));
           const LocalPoint lpCSC_strip3(keyLayerGeometry3->topology()->localPosition(b_CSCstrip));
           const LocalPoint lpCSC_strip4(keyLayerGeometry4->topology()->localPosition(b_CSCstrip));
           const LocalPoint lpCSC_strip5(keyLayerGeometry5->topology()->localPosition(b_CSCstrip));
           const LocalPoint lpCSC_strip6(keyLayerGeometry6->topology()->localPosition(b_CSCstrip));
           const LocalPoint lpCSC_wg1(keyLayerGeometry1->topology()->localPosition(b_CSCkeyWire));
           const LocalPoint lpCSC_wg2(keyLayerGeometry2->topology()->localPosition(b_CSCkeyWire));
           const LocalPoint lpCSC_wg3(keyLayerGeometry3->topology()->localPosition(b_CSCkeyWire));
           const LocalPoint lpCSC_wg4(keyLayerGeometry4->topology()->localPosition(b_CSCkeyWire));
           const LocalPoint lpCSC_wg5(keyLayerGeometry5->topology()->localPosition(b_CSCkeyWire));
           const LocalPoint lpCSC_wg6(keyLayerGeometry6->topology()->localPosition(b_CSCkeyWire));
           const GlobalPoint gpst(keyLayer3->toGlobal(lpCSC_strip3));
           const GlobalPoint gpwg(keyLayer3->toGlobal(lpCSC_wg3));

 	   float X0=lpCSC_strip3.x();
 	   float Y0=lpCSC_strip3.y();
	   float Z0=lpCSC_strip3.z();

          // const LocalPoint lpRPC(randRoll->toLocal(gp));

           cout << "LocalPoint1 strip " << lpCSC_strip1 << endl;
           cout << "x" << X0 << "y" << Y0 << "z" << Z0<< endl;

           cout << "\nLocalPoint2 strip " << lpCSC_strip2 << endl;
           cout << "LocalPoint3 strip " << lpCSC_strip3 << endl;
           cout << "LocalPoint4 strip " << lpCSC_strip4 << endl;
           cout << "LocalPoint5 strip " << lpCSC_strip5 << endl;
           cout << "LocalPoint6 strip " << lpCSC_strip6 << endl;
           cout << "LocalPoint wiregroup " << lpCSC_wg1 << endl;
           cout << "LocalPoint wiregroup " << lpCSC_wg2 << endl;
           cout << "LocalPoint wiregroup " << lpCSC_wg3 << endl;
           cout << "LocalPoint wiregroup " << lpCSC_wg4 << endl;
           cout << "LocalPoint wiregroup " << lpCSC_wg5 << endl;
           cout << "LocalPoint wiregroup " << lpCSC_wg6 << endl;
           cout << "GlobalPoint st3" << gpst << endl;
           cout << "GlobalPoint wg3" << gpwg << endl;
      


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

       cout << "\n\nnumber of digis: " << b_numberofDigis << endl;   
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
