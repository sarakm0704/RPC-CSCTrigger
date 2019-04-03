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
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TH1.h"
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/RecSegment.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "Geometry/RPCGeometry/interface/RPCRollSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "RecoLocalMuon/RPCRecHit/interface/CSCSegtoRPC.h"
#include "RecoLocalMuon/RPCRecHit/src/CSCObjectMap.h"
#include "RecoLocalMuon/RPCRecHit/src/CSCStationIndex.h"

#include "L1Trigger/CSCTriggerPrimitives/src/CSCMotherboard.h"
#include "L1Trigger/CSCTriggerPrimitives/src/CSCMotherboardME3141RPC.h"

using namespace edm;
using namespace std;

class CSCExtrapoltoRPC : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit CSCExtrapoltoRPC(const edm::ParameterSet&);
    ~CSCExtrapoltoRPC();

    std::unique_ptr<RPCRecHitCollection> && thePoints(){ return std::move(_ThePoints); }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    //https://github.com/cms-sw/cmssw/blob/02d4198c0b6615287fd88e9a8ff650aea994412e/L1Trigger/L1TMuon/interface/GeometryTranslator.h
    const CSCGeometry& getCSCGeometry() const { return *cscGeo; }

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    TTree *tree;
    TH1D *EventInfo;

    TH1D *Ndigis;
    TH1D *fNdigis;
    TH1D *bNdigis;
    TH1D *ME31NDigis;
    TH1D *ME41NDigis;
    TH1D *Nrechit;

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

    std::unique_ptr<RPCRecHitCollection> _ThePoints;
    edm::EDGetTokenT<MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>> corrlctsToken_;
    edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;

    edm::ESHandle<CSCGeometry> cscGeo;
    edm::ESHandle<RPCGeometry> rpcGeo;

    //GlobalPoint getGlobalPosition(unsigned int rawId, const CSCCorrelatedLCTDigi& lct) const;
    GlobalPoint getGlobalPosition(CSCDetId rawId, const CSCCorrelatedLCTDigi& lct) const;
    const CSCGeometry* cscGeomtery_;
};


GlobalPoint
//CSCExtrapoltoRPC::getGlobalPosition(unsigned int rawId, const CSCCorrelatedLCTDigi& lct) const{ //105x
CSCExtrapoltoRPC::getGlobalPosition(CSCDetId rawId, const CSCCorrelatedLCTDigi& lct) const{

  CSCDetId cscId = CSCDetId(rawId);
  CSCDetId key_id(cscId.endcap(), cscId.station(), cscId.ring(), cscId.chamber(), CSCConstants::KEY_CLCT_LAYER);

  const auto& cscChamber = getCSCGeometry().chamber(cscId);
  float fractional_strip = lct.getFractionalStrip();
  const auto& layer_geo = cscChamber->layer(CSCConstants::KEY_CLCT_LAYER)->geometry();

  // LCT::getKeyWG() also starts from 0
  float wire = layer_geo->middleWireOfGroup(lct.getKeyWG() + 1);
  const LocalPoint& csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, wire);
  const GlobalPoint& csc_gp = cscGeo->idToDet(key_id)->surface().toGlobal(csc_intersect);

  return csc_gp;

}

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
  Ndigis->GetYaxis()->SetTitle("Number of digis");

  fNdigis = fs->make<TH1D>("fNdigis", "number of digis per chamber (forward)", 10, 0, 10);
  bNdigis = fs->make<TH1D>("bNdigis", "number of digis per chamber (backward)", 10, 0, 10);

  ME31NDigis = fs->make<TH1D>("ME31NDigis", "number of digis per chamber (ME3/1)", 10, 0, 10);
  ME31NDigis->GetXaxis()->SetTitle("Number of digis per chamber");
  ME31NDigis->GetYaxis()->SetTitle("Number of digis");

  ME41NDigis = fs->make<TH1D>("ME41NDigis", "number of digis per chamber (ME4/1)", 10, 0, 10);
  ME41NDigis->GetXaxis()->SetTitle("Number of digis per chamber");
  ME41NDigis->GetYaxis()->SetTitle("Number of digis");

  Nrechit = fs->make<TH1D>("Nrechit", "", 10, 0, 10);
  Nrechit->GetXaxis()->SetTitle("Number of rechit per chamber");
  Nrechit->GetYaxis()->SetTitle("Number of digis");
   
}

CSCExtrapoltoRPC::~CSCExtrapoltoRPC()
{

}

void
CSCExtrapoltoRPC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  iSetup.get<MuonGeometryRecord>().get( cscGeo );     
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
  b_CSCendcap = b_CSCstation = b_CSCsubsector = b_CSCsector = b_CSCstrip = b_CSCkeyWire = 0;

  b_Trknmb = 0;
  b_cscBX = 0;
  b_cscId = 0;

  b_EVENT  = iEvent.id().event();
  b_RUN    = iEvent.id().run();
  b_LUMI   = iEvent.id().luminosityBlock();

  std::cout << "New event\n" << endl;
    for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc=corrlcts.product()->begin(); csc!=corrlcts.product()->end(); csc++){  
      CSCCorrelatedLCTDigiCollection::Range range1 = corrlcts.product()->get((*csc).first);
//       cout << "\nwe are in ranage 1 " << endl;
      b_numberofDigis = b_fNDigis = b_bNDigis = b_ME31NDigis = b_ME41NDigis = 0;

      for(CSCCorrelatedLCTDigiCollection::const_iterator lct=range1.first; lct!=range1.second; lct++){

        cout << "\nnew iterator" << endl;
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
        cout << "I'm here in CSC:: endcap: "    << b_CSCendcap <<
                                " station: "    << b_CSCstation <<
                                " sector: "     << b_CSCsector <<
                                " subsector: "  << b_CSCsubsector <<
                                " strip: "      << b_CSCstrip <<
                                " wire: "       << b_CSCkeyWire << endl;

        cout << "I'm here in CSC:: endcap: " << csc_id.endcap() <<
                                " station: " << csc_id.station() << endl;

        //CSCTriggerGeomManager* geo_manager(CSCTriggerGeometry::get());
        //const CSCChamber* cscChamber(geo_manager->chamber(newEndcap, newStation, b_CSCsector, b_CSCsubsector, b_cscId));

        //const RPCChamber* rpcChamber=rpcGeo->chamber(rpc_id);
        //const LocalPoint lpRPC(randRoll->toLocal(gp));
        //auto randRoll(rpcChamber->roll(2));
                     
        const CSCChamber* cscChamber=cscGeo->chamber(csc_id);
        const CSCLayer* keyLayer3(cscChamber->layer(3));
        const CSCLayerGeometry* keyLayerGeometry3(keyLayer3->geometry());
        const LocalPoint lpCSC_strip3(keyLayerGeometry3->topology()->localPosition(b_CSCstrip));
        const GlobalPoint gpst(keyLayer3->toGlobal(lpCSC_strip3));
        //const LocalPoint lpCSC_wg3(keyLayerGeometry3->topology()->localPosition(b_CSCkeyWire));
        //const GlobalPoint gpwg(keyLayer3->toGlobal(lpCSC_wg3));

        float fractional_strip = 0.5 * (lct->getStrip() + 1) - 0.25;
        auto layer_geo = cscChamber->layer(CSCConstants::KEY_CLCT_LAYER)->geometry();
        float wire = layer_geo->middleWireOfGroup(lct->getKeyWG() + 1);
        LocalPoint csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, wire);
        GlobalPoint gp_intersect(keyLayer3->toGlobal(csc_intersect));

        cout << "LocalPoint3 strip " << lpCSC_strip3 << endl;
        cout << "GlobalPoint st3" << gpst << endl;
        cout << "intersection" << csc_intersect << endl;     
        //cout << "LocalPoint wiregroup " << lpCSC_wg3 << endl;
        //cout << "GlobalPoint wg3" << gpwg << endl;

        GlobalPoint gp_test(getGlobalPosition(csc_id, *lct));
        cout << "Globalfunction intersect test: " << gp_intersect << ":" << gp_test << endl;

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

      cout << "\nnumber of digis: " << b_numberofDigis << endl;   
      //cout << "ME31 NDigis: " << b_ME31NDigis << " ME41 NDigis: " << b_ME41NDigis << endl;   

      //fill histo
      Ndigis->Fill(b_numberofDigis);
      if (b_fNDigis != 0) fNdigis->Fill(b_fNDigis);
      if (b_bNDigis != 0) bNdigis->Fill(b_bNDigis);
      if (b_ME31NDigis != 0) ME31NDigis->Fill(b_ME31NDigis);
      if (b_ME41NDigis != 0) ME41NDigis->Fill(b_ME41NDigis);

  }

  int nRPC = 0;
  RPCRecHitCollection::const_iterator rpcIt;
  for (rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {
    cout << "I'm new rechits" << endl;

    nRPC++;
    RPCDetId rpcid = (RPCDetId)(*rpcIt).rpcId();

    int kRoll  = rpcid.roll();
    int kSubsector  = rpcid.subsector();
    int kRegion  = rpcid.region();
    int kStation = rpcid.station();
    int kRing = rpcid.ring();
    int kSector = rpcid.sector();
    int kLayer = rpcid.layer();
    int bx = (*rpcIt).BunchX();
    int clSize = (*rpcIt).clusterSize();

    cout << "Region/Station/Ring/Sector/SubSector/Roll/Layer: " <<
            kRegion << " / " << kStation << " / " << kRing << " / " <<
            kSector << " / " << kSubsector <<  " / " <<
            kRoll << " / " << kLayer <<
            "\tbx,clSize: " << bx << ", " << clSize << endl;    
    cout << "number of RPCRecHits: " << nRPC << endl;
    Nrechit->Fill(nRPC);
  }

  tree->Fill();
  EventInfo->Fill(1.5);

}

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

void 
CSCExtrapoltoRPC::endJob() 
{

}

void
CSCExtrapoltoRPC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(CSCExtrapoltoRPC);
