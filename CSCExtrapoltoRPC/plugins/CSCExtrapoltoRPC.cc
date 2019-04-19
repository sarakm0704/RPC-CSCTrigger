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
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
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
//#include "L1Trigger/CSCTriggerPrimitives/src/CSCMotherboardME3141RPC.h"
#include <L1Trigger/CSCCommonTrigger/interface/CSCConstants.h>

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

    TH1D *h_ME31NDigis;
    TH1D *h_ME41NDigis;
    TH1D *h_ME31NDigis0;
    TH1D *h_ME41NDigis0;

    TH1D *h_ME31NDigis_a;
    TH1D *h_ME41NDigis_a;
    TH1D *h_ME31NDigis0_a;
    TH1D *h_ME41NDigis0_a;

    TH1D *NRecHits;
    TH1D *h_RE31NRecHits;
    TH1D *h_RE41NRecHits;

    TH1D *h_S3NDigis;
    TH1D *h_S4NDigis;
    TH1D *h_S3NRecHits;
    TH1D *h_S4NRecHits;

    unsigned int b_EVENT, b_RUN, b_LUMI;
    unsigned int b_numberofDigis;

    unsigned int b_S3NDigis;
    unsigned int b_S4NDigis;

    unsigned int b_ME11NDigis;
    unsigned int b_ME21NDigis;

    double b_ME31NDigis;
    double b_ME41NDigis;
    int bx_ME31NDigis;
    int bx_ME41NDigis;

    int b_ME31NDigis_Total;
    int b_ME41NDigis_Total;

    unsigned int a_ME31NDigis;
    unsigned int a_ME41NDigis;

    int nRPC;
    unsigned int b_S3NRecHits;
    unsigned int b_S4NRecHits;
    unsigned int b_RE31NRecHits;
    unsigned int b_RE41NRecHits;

    int b_Trknmb;
    int b_cscId;

    int b_rpcBX;
    int b_cscBX;

    TH1D *h_xNMatchedME31;
    TH1D *h_xNMatchedME41;

    TH1D *h_yNMatchedME31;
    TH1D *h_yNMatchedME41;
  
    TH2D *h_MatchedME31;
    TH2D *h_MatchedME41;

    double xME3115;
    double xME3114;
    double xME3113;
    double xME3112;
    double xME3111;
    double xME3110;
    double xME3109;
    double xME3108;
    double xME3107;
    double xME3106;
    double xME3105;
    double xME3104;
    double xME3103;
    double xME3102;
    double xME3101;

    double xME4115;
    double xME4114;
    double xME4113;
    double xME4112;
    double xME4111;
    double xME4110;
    double xME4109;
    double xME4108;
    double xME4107;
    double xME4106;
    double xME4105;
    double xME4104;
    double xME4103;
    double xME4102;
    double xME4101;

    double yME3115;
    double yME3114;
    double yME3113;
    double yME3112;
    double yME3111;
    double yME3110;
    double yME3109;
    double yME3108;
    double yME3107;
    double yME3106;
    double yME3105;
    double yME3104;
    double yME3103;
    double yME3102;
    double yME3101;

    double yME4115;
    double yME4114;
    double yME4113;
    double yME4112;
    double yME4111;
    double yME4110;
    double yME4109;
    double yME4108;
    double yME4107;
    double yME4106;
    double yME4105;
    double yME4104;
    double yME4103;
    double yME4102;
    double yME4101;

    double ME31[15][15];
    double ME41[15][15];

    std::unique_ptr<RPCRecHitCollection> _ThePoints;
    edm::EDGetTokenT<MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>> corrlctsToken_;
    edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;

    edm::ESHandle<CSCGeometry> cscGeo;
    edm::ESHandle<RPCGeometry> rpcGeo;

    //GlobalPoint getGlobalPosition(unsigned int rawId, const CSCCorrelatedLCTDigi& lct) const;
    GlobalPoint getCSCGlobalPosition(CSCDetId rawId, const CSCCorrelatedLCTDigi& lct) const;
    GlobalPoint getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& rpcIt) const;
    const CSCGeometry* cscGeomtery_;
};


GlobalPoint
//CSCExtrapoltoRPC::getGlobalPosition(unsigned int rawId, const CSCCorrelatedLCTDigi& lct) const{ //105x
CSCExtrapoltoRPC::getCSCGlobalPosition(CSCDetId rawId, const CSCCorrelatedLCTDigi& lct) const{

  CSCDetId cscId = CSCDetId(rawId);
  CSCDetId key_id(cscId.endcap(), cscId.station(), cscId.ring(), cscId.chamber(), CSCConstants::KEY_CLCT_LAYER);

  const auto& cscChamber = getCSCGeometry().chamber(cscId);
  float fractional_strip = lct.getFractionalStrip();
  //CSCs have 6 layers. The key (refernce) layer is the third layer (CSCConstant)
  const auto& layer_geo = cscChamber->layer(CSCConstants::KEY_CLCT_LAYER)->geometry();

  // LCT::getKeyWG() also starts from 0
  float wire = layer_geo->middleWireOfGroup(lct.getKeyWG() + 1);
  const LocalPoint& csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, wire);
  const GlobalPoint& csc_gp = cscGeo->idToDet(key_id)->surface().toGlobal(csc_intersect);

  return csc_gp;

}

GlobalPoint
CSCExtrapoltoRPC::getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& rpcIt) const{

  RPCDetId rpcid = RPCDetId(rpcId);
  const LocalPoint& rpc_lp = rpcIt.localPosition();
  const GlobalPoint& rpc_gp = rpcGeo->idToDet(rpcid)->surface().toGlobal(rpc_lp);
 
  return rpc_gp;

}

/*
RPCRecHit
CSCExtrapoltoRPC::matchingRPC(const CSCCorrelatedLCTDigi& lct, RPCBXwindow[:]) const{

  for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {

    if(matched == 0) return nullptr;
    if(matched > 0) return rpcIt;

  }

}
*/


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

  h_xNMatchedME31 = fs->make<TH1D>("h_xNMatchedME31", "", 15, 0, 15);
  h_xNMatchedME31->GetXaxis()->SetTitle("X cutoff (cm)");
  h_xNMatchedME31->GetYaxis()->SetTitle("Matched (%)");

  h_xNMatchedME41 = fs->make<TH1D>("h_xNMatchedME41", "", 15, 0, 15);
  h_xNMatchedME41->GetXaxis()->SetTitle("X cutoff (cm)");
  h_xNMatchedME41->GetYaxis()->SetTitle("Matched (%)");

  h_yNMatchedME31 = fs->make<TH1D>("h_yNMatchedME31", "", 15, 0, 15);
  h_yNMatchedME31->GetXaxis()->SetTitle("Y cutoff (cm)");
  h_yNMatchedME31->GetYaxis()->SetTitle("Matched (%)");

  h_yNMatchedME41 = fs->make<TH1D>("h_yNMatchedME41", "", 15, 0, 15);
  h_yNMatchedME41->GetXaxis()->SetTitle("Y cutoff (cm)");
  h_yNMatchedME41->GetYaxis()->SetTitle("Matched (%)");

  h_ME31NDigis = fs->make<TH1D>("ME31NDigis_before", "number of digi per chamber (ME3/1)", 10, 0, 10);
  h_ME31NDigis->GetXaxis()->SetTitle("Number of digi per chamber");
  h_ME31NDigis->GetYaxis()->SetTitle("Number of chamber");

  h_ME41NDigis = fs->make<TH1D>("ME41NDigis_before", "number of digi per chamber (ME4/1)", 10, 0, 10);
  h_ME41NDigis->GetXaxis()->SetTitle("Number of digi per chamber");
  h_ME41NDigis->GetYaxis()->SetTitle("Number of chamber");

  h_ME31NDigis0 = fs->make<TH1D>("ME31NDigis0_before", "number of digi per chamber (ME3/1)", 10, 0, 10);
  h_ME31NDigis0->GetXaxis()->SetTitle("Number of digi per chamber");
  h_ME31NDigis0->GetYaxis()->SetTitle("Number of chamber");

  h_ME41NDigis0 = fs->make<TH1D>("ME41NDigis0_before", "number of digi per chamber (ME4/1)", 10, 0, 10);
  h_ME41NDigis0->GetXaxis()->SetTitle("Number of digi per chamber");
  h_ME41NDigis0->GetYaxis()->SetTitle("Number of chamber");

  h_ME31NDigis_a = fs->make<TH1D>("ME31NDigis_after", "number of digi per chamber (ME3/1)", 10, 0, 10);
  h_ME31NDigis_a->GetXaxis()->SetTitle("Number of digi per chamber");
  h_ME31NDigis_a->GetYaxis()->SetTitle("Number of chamber");

  h_ME41NDigis_a = fs->make<TH1D>("ME41NDigis_after", "number of digi per chamber (ME4/1)", 10, 0, 10);
  h_ME41NDigis_a->GetXaxis()->SetTitle("Number of digi per chamber");
  h_ME41NDigis_a->GetYaxis()->SetTitle("Number of chamber");

  h_ME31NDigis0_a = fs->make<TH1D>("ME31NDigis0_after", "number of digi per chamber (ME3/1)", 10, 0, 10);
  h_ME31NDigis0_a->GetXaxis()->SetTitle("Number of digi per chamber");
  h_ME31NDigis0_a->GetYaxis()->SetTitle("Number of chamber");

  h_ME41NDigis0_a = fs->make<TH1D>("ME41NDigis0_after", "number of digi per chamber (ME4/1)", 10, 0, 10);
  h_ME41NDigis0_a->GetXaxis()->SetTitle("Number of digi per chamber");
  h_ME41NDigis0_a->GetYaxis()->SetTitle("Number of chamber");

  NRecHits = fs->make<TH1D>("NRecHits", "", 10, 0, 10);
  NRecHits->GetXaxis()->SetTitle("Number of rechit per chamber");
  NRecHits->GetYaxis()->SetTitle("Number of chamber");

  h_RE31NRecHits = fs->make<TH1D>("h_RE31NRecHits", "number of rechit per chamber (RE3/1)", 20, 0, 20);
  h_RE31NRecHits->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_RE31NRecHits->GetYaxis()->SetTitle("Number of chamber");
   
  h_RE41NRecHits = fs->make<TH1D>("h_RE41NRecHits", "number of rechit per chamber (RE4/1)", 20, 0, 20);
  h_RE41NRecHits->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_RE41NRecHits->GetYaxis()->SetTitle("Number of chamber");

  h_S3NDigis = fs->make<TH1D>("h_S3NDigis", "number of digi in station 3", 10, 0, 10);
  h_S3NDigis->GetXaxis()->SetTitle("Number of digi");
  h_S3NDigis->GetYaxis()->SetTitle("Number of chamber");
  
  h_S4NDigis = fs->make<TH1D>("h_S4NDigis", "number of digi in station 4", 10, 0, 10);
  h_S4NDigis->GetXaxis()->SetTitle("Number of digi");
  h_S4NDigis->GetYaxis()->SetTitle("Number of chamber");

  h_S3NRecHits = fs->make<TH1D>("h_S3NRecHits", "number of rechit in station 3", 20, 0, 20);
  h_S3NRecHits->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_S3NRecHits->GetYaxis()->SetTitle("Number of chamber");
  
  h_S4NRecHits = fs->make<TH1D>("h_S4NRecHits", "number of rechit in station 4", 20, 0, 20);
  h_S4NRecHits->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_S4NRecHits->GetYaxis()->SetTitle("Number of chamber");

  h_MatchedME31 = fs->make<TH2D>("h_MatchedME31", "", 15, 0, 15, 15, 0, 15);
  h_MatchedME31->GetXaxis()->SetTitle("X cutoff (cm)");
  h_MatchedME31->GetYaxis()->SetTitle("Y cutoff (cm)");

  h_MatchedME41 = fs->make<TH2D>("h_MatchedME41", "", 15, 0, 15, 15, 0, 15);
  h_MatchedME41->GetXaxis()->SetTitle("X cutoff (cm)");
  h_MatchedME41->GetYaxis()->SetTitle("Y cutoff (cm)");

}

CSCExtrapoltoRPC::~CSCExtrapoltoRPC()
{

}

void
CSCExtrapoltoRPC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  iSetup.get<MuonGeometryRecord>().get( cscGeo );     
  iSetup.get<MuonGeometryRecord>().get( rpcGeo );     

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

  b_EVENT  = iEvent.id().event();
  b_RUN    = iEvent.id().run();
  b_LUMI   = iEvent.id().luminosityBlock();

  std::vector<GlobalPoint> rpcMatched_gp;

  bool isMatchxME3115 = false;
  bool isMatchxME3114 = false;
  bool isMatchxME3113 = false;
  bool isMatchxME3112 = false;
  bool isMatchxME3111 = false;
  bool isMatchxME3110 = false;
  bool isMatchxME3109 = false;
  bool isMatchxME3108 = false;
  bool isMatchxME3107 = false;
  bool isMatchxME3106 = false;
  bool isMatchxME3105 = false;
  bool isMatchxME3104 = false;
  bool isMatchxME3103 = false;
  bool isMatchxME3102 = false;
  bool isMatchxME3101 = false;
  bool isMatchyME3115 = false;
  bool isMatchyME3114 = false;
  bool isMatchyME3113 = false;
  bool isMatchyME3112 = false;
  bool isMatchyME3111 = false;
  bool isMatchyME3110 = false;
  bool isMatchyME3109 = false;
  bool isMatchyME3108 = false;
  bool isMatchyME3107 = false;
  bool isMatchyME3106 = false;
  bool isMatchyME3105 = false;
  bool isMatchyME3104 = false;
  bool isMatchyME3103 = false;
  bool isMatchyME3102 = false;
  bool isMatchyME3101 = false;

  bool isMatchxME4115 = false;
  bool isMatchxME4114 = false;
  bool isMatchxME4113 = false;
  bool isMatchxME4112 = false;
  bool isMatchxME4111 = false;
  bool isMatchxME4110 = false;
  bool isMatchxME4109 = false;
  bool isMatchxME4108 = false;
  bool isMatchxME4107 = false;
  bool isMatchxME4106 = false;
  bool isMatchxME4105 = false;
  bool isMatchxME4104 = false;
  bool isMatchxME4103 = false;
  bool isMatchxME4102 = false;
  bool isMatchxME4101 = false;
  bool isMatchyME4115 = false;
  bool isMatchyME4114 = false;
  bool isMatchyME4113 = false;
  bool isMatchyME4112 = false;
  bool isMatchyME4111 = false;
  bool isMatchyME4110 = false;
  bool isMatchyME4109 = false;
  bool isMatchyME4108 = false;
  bool isMatchyME4107 = false;
  bool isMatchyME4106 = false;
  bool isMatchyME4105 = false;
  bool isMatchyME4104 = false;
  bool isMatchyME4103 = false;
  bool isMatchyME4102 = false;
  bool isMatchyME4101 = false;

  bool isMatchME31[15][15];
  bool isMatchME41[15][15];

  cout << "\nNew event" << endl;
  for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc=corrlcts.product()->begin(); csc!=corrlcts.product()->end(); csc++){  
    CSCCorrelatedLCTDigiCollection::Range range1 = corrlcts.product()->get((*csc).first);
    b_numberofDigis = b_ME11NDigis  = b_ME21NDigis = b_ME31NDigis = b_ME41NDigis = 0;
    bx_ME31NDigis = bx_ME41NDigis =0;
    b_cscBX = b_rpcBX = 0;

    a_ME31NDigis = a_ME41NDigis = 0;

    const CSCDetId csctest_id((*csc).first.rawId());

    int ismatched = 0;

    b_S3NDigis = b_S4NDigis = 0;
    b_S3NRecHits = b_S4NRecHits = 0;

    rpcMatched_gp.clear();
    GlobalPoint gp_cscint;

    nRPC = b_RE31NRecHits = b_RE41NRecHits = 0;
    //this is only for the number of station per CSCchamber
    for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {

      RPCDetId rpcid = (RPCDetId)(*rpcIt).rpcId();
      nRPC++;
      if(rpcid.station() == 3) b_S3NRecHits++;
      if(rpcid.station() == 4) b_S4NRecHits++;

      if(rpcid.station() == 3 && rpcid.ring() == 1) b_RE31NRecHits++;
      if(rpcid.station() == 4 && rpcid.ring() == 1) b_RE41NRecHits++;
    }


    for(CSCCorrelatedLCTDigiCollection::const_iterator lct=range1.first; lct!=range1.second; lct++){
      const CSCDetId csc_id((*csc).first.rawId());
      b_cscBX = lct->getBX();
      int strip = lct->getStrip();
      int wg = lct->getKeyWG();

      cout << "strip / wg" << strip << " / " << wg << endl;

      gp_cscint = GlobalPoint(0.0,0.0,0.0);
      gp_cscint = getCSCGlobalPosition(csc_id, *lct);

      b_numberofDigis++;

      if (csc_id.station() == 3) b_S3NDigis++;
      if (csc_id.station() == 4) b_S4NDigis++;

      if (csc_id.station() == 3 && csc_id.ring() == 1){

        b_ME31NDigis = b_ME31NDigis + 1;
        if (b_cscBX == 6) b_ME31NDigis_Total = b_ME31NDigis_Total + 1;
        if (b_cscBX == 6) bx_ME31NDigis++;

        GlobalPoint gp_ME31(getCSCGlobalPosition(csc_id, *lct));
      }
      else if (csc_id.station() == 4 && csc_id.ring() == 1){

        b_ME41NDigis = b_ME41NDigis + 1;
        if (b_cscBX == 6) b_ME41NDigis_Total = b_ME41NDigis_Total + 1;
        if (b_cscBX == 6) bx_ME41NDigis++;

        GlobalPoint gp_ME41(getCSCGlobalPosition(csc_id, *lct));
      }

      int tmp_matched = ismatched;

      for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {

        RPCDetId rpcid = (RPCDetId)(*rpcIt).rpcId();
        b_rpcBX = (*rpcIt).BunchX();

        GlobalPoint gp_rpc(0.0,0.0,0.0);
        gp_rpc = getRPCGlobalPosition(rpcid, *rpcIt);

        if (rpcid.region() == 0) continue; //skip the barrels

        if (gp_rpc.x() == 0 && gp_rpc.y() == 0 && gp_rpc.z() == 0 ) continue;
        if (gp_cscint.x() == 0 && gp_cscint.y() == 0 && gp_cscint.z() == 0 ) continue;
        
        float Rx = gp_rpc.x();
        float Ry = gp_rpc.y();

        float Cx = gp_cscint.x();
        float Cy = gp_cscint.y();

        float extrapol2D = sqrt((Rx-Cx)*(Rx-Cx)+(Ry-Cy)*(Ry-Cy));

        if (extrapol2D < 10){

          bool dupl = false;

          for (std::vector<GlobalPoint>::const_iterator rpcMatching = rpcMatched_gp.begin(); rpcMatching != rpcMatched_gp.end(); rpcMatching++){         
            GlobalPoint a = *rpcMatching; 
            float tmp_dist = sqrt( (gp_rpc.x()-a.x())*(gp_rpc.x()-a.x()) + (gp_rpc.y()-a.y())*(gp_rpc.y()-a.y()) );
            if ( tmp_dist < 5 ) dupl = true;
          }

          if ( dupl != true ){
            if (csc_id.endcap() == 1 && csc_id.station() == 3 && csc_id.ring() == 1 &&
                rpcid.region() == 1 && rpcid.station() == 3 && rpcid.ring() == 1){
              ismatched++;
              rpcMatched_gp.push_back(gp_rpc);
            }
            if (csc_id.endcap() == 1 && csc_id.station() == 4 && csc_id.ring() == 1 &&
                rpcid.region() == 1 && rpcid.station() == 4 && rpcid.ring() == 1){
              ismatched++;
              rpcMatched_gp.push_back(gp_rpc);
            }
            if (csc_id.endcap() == 2 && csc_id.station() == 3 && csc_id.ring() == 1 &&
                rpcid.region() == -1 && rpcid.station() == 3 && rpcid.ring() == 1){
              ismatched++;
              rpcMatched_gp.push_back(gp_rpc);
            }
            if (csc_id.endcap() == 2 && csc_id.station() == 4 && csc_id.ring() == 1 &&
                rpcid.region() == -1 && rpcid.station() == 4 && rpcid.ring() == 1){
              ismatched++;
              rpcMatched_gp.push_back(gp_rpc);
            }
          }
        }
        if ( ismatched != tmp_matched ) break; //move on to next lct

/*
        int kRoll  = rpcid.roll();
        int kSubsector  = rpcid.subsector();
        int kRegion  = rpcid.region();
        int kStation = rpcid.station();
        int kRing = rpcid.ring();
        int kSector = rpcid.sector();
        int kLayer = rpcid.layer();
        int bx = (*rpcIt).BunchX();
        int clSize = (*rpcIt).clusterSize();
    
        cout << "I'm here in RPC Region: " << kRegion <<
                              " Station: " << kStation <<
                              " Ring: "    << kRing <<
                              " Sector: "  << kSector <<
                            " Subsector: " << kSubsector <<
                                 " Roll: " << kRoll <<
                                " Layer: " << kLayer <<
                           "\tbx,clSize: " << bx << ", " << clSize << endl;    
        cout << "number of RPCRecHits: " << nRPC << endl;
    
        int nrechitRE31 = 0;
        int nrechitRE41 = 0;
        if (kRegion != 0){
          if (kStation == 3 && kRing == 1) nrechitRE31++;    
          else if(kStation == 4 && kRing == 1) nrechitRE41++;
    
        }
    
    */
      }//RPCRecHit loop

      for (int i=0; i<15; i++){
        for (int j=0; j<15; j++){
          isMatchME31[i][j] = false;
          isMatchME41[i][j] = false;
        }
      }
      isMatchxME3115 = false;
      isMatchxME3114 = false;
      isMatchxME3113 = false;
      isMatchxME3112 = false;
      isMatchxME3111 = false;
      isMatchxME3110 = false;
      isMatchxME3109 = false;
      isMatchxME3108 = false;
      isMatchxME3107 = false;
      isMatchxME3106 = false;
      isMatchxME3105 = false;
      isMatchxME3104 = false;
      isMatchxME3103 = false;
      isMatchxME3102 = false;
      isMatchxME3101 = false;
      isMatchxME4115 = false;
      isMatchxME4114 = false;
      isMatchxME4113 = false;
      isMatchxME4112 = false;
      isMatchxME4111 = false;
      isMatchxME4110 = false;
      isMatchxME4109 = false;
      isMatchxME4108 = false;
      isMatchxME4107 = false;
      isMatchxME4106 = false;
      isMatchxME4105 = false;
      isMatchxME4104 = false;
      isMatchxME4103 = false;
      isMatchxME4102 = false;
      isMatchxME4101 = false;

      isMatchyME3115 = false;
      isMatchyME3114 = false;
      isMatchyME3113 = false;
      isMatchyME3112 = false;
      isMatchyME3111 = false;
      isMatchyME3110 = false;
      isMatchyME3109 = false;
      isMatchyME3108 = false;
      isMatchyME3107 = false;
      isMatchyME3106 = false;
      isMatchyME3105 = false;
      isMatchyME3104 = false;
      isMatchyME3103 = false;
      isMatchyME3102 = false;
      isMatchyME3101 = false;
      isMatchyME4115 = false;
      isMatchyME4114 = false;
      isMatchyME4113 = false;
      isMatchyME4112 = false;
      isMatchyME4111 = false;
      isMatchyME4110 = false;
      isMatchyME4109 = false;
      isMatchyME4108 = false;
      isMatchyME4107 = false;
      isMatchyME4106 = false;
      isMatchyME4105 = false;
      isMatchyME4104 = false;
      isMatchyME4103 = false;
      isMatchyME4102 = false;
      isMatchyME4101 = false;

      for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {

        RPCDetId rpcid = (RPCDetId)(*rpcIt).rpcId();

        GlobalPoint gp_rpc(0.0,0.0,0.0);
        gp_rpc = getRPCGlobalPosition(rpcid, *rpcIt);

        if (rpcid.region() == 0) continue; //skip the barrels

        if (gp_rpc.x() == 0 && gp_rpc.y() == 0 && gp_rpc.z() == 0 ) continue;
        if (gp_cscint.x() == 0 && gp_cscint.y() == 0 && gp_cscint.z() == 0 ) continue;

        float Dx = abs(gp_rpc.x()-gp_cscint.x());
        float Dy = abs(gp_rpc.y()-gp_cscint.y());

        if (csc_id.station() == 3 && csc_id.ring() == 1 && rpcid.station() == 3 && rpcid.ring() == 1){
          for (int i = 0; i < 15; i++){
            for (int j = 0; j < 15; j++){
              if (Dx < i+1 && Dy < j+1) isMatchME31[i][j] = true;
            }
          }
          if (Dx < 15) isMatchxME3115 = true;
          if (Dx < 14) isMatchxME3114 = true;
          if (Dx < 13) isMatchxME3113 = true;
          if (Dx < 12) isMatchxME3112 = true;
          if (Dx < 11) isMatchxME3111 = true;
          if (Dx < 10) isMatchxME3110 = true;
          if (Dx < 9) isMatchxME3109 = true;
          if (Dx < 8) isMatchxME3108 = true;
          if (Dx < 7) isMatchxME3107 = true;
          if (Dx < 6) isMatchxME3106 = true;
          if (Dx < 5) isMatchxME3105 = true;
          if (Dx < 4) isMatchxME3104 = true;
          if (Dx < 3) isMatchxME3103 = true;
          if (Dx < 2) isMatchxME3102 = true;
          if (Dx < 1) isMatchxME3101 = true;
          if (Dy < 15) isMatchyME3115 = true;
          if (Dy < 14) isMatchyME3114 = true;
          if (Dy < 13) isMatchyME3113 = true;
          if (Dy < 12) isMatchyME3112 = true;
          if (Dy < 11) isMatchyME3111 = true;
          if (Dy < 10) isMatchyME3110 = true;
          if (Dy < 9) isMatchyME3109 = true;
          if (Dy < 8) isMatchyME3108 = true;
          if (Dy < 7) isMatchyME3107 = true;
          if (Dy < 6) isMatchyME3106 = true;
          if (Dy < 5) isMatchyME3105 = true;
          if (Dy < 4) isMatchyME3104 = true;
          if (Dy < 3) isMatchyME3103 = true;
          if (Dy < 2) isMatchyME3102 = true;
          if (Dy < 1) isMatchyME3101 = true;
        }
        if (csc_id.station() == 4 && csc_id.ring() == 1 && rpcid.station() == 4 && rpcid.ring() == 1){
          for (int i = 0; i < 15; i++){
            for (int j = 0; j < 15; j++){
              if (Dx < i+1 && Dy < j+1) isMatchME41[i][j] = true;
            }
          }
          if (Dx < 15) isMatchxME4115 = true;
          if (Dx < 14) isMatchxME4114 = true;
          if (Dx < 13) isMatchxME4113 = true;
          if (Dx < 12) isMatchxME4112 = true;
          if (Dx < 11) isMatchxME4111 = true;
          if (Dx < 10) isMatchxME4110 = true;
          if (Dx < 9) isMatchxME4109 = true;
          if (Dx < 8) isMatchxME4108 = true;
          if (Dx < 7) isMatchxME4107 = true;
          if (Dx < 6) isMatchxME4106 = true;
          if (Dx < 5) isMatchxME4105 = true;
          if (Dx < 4) isMatchxME4104 = true;
          if (Dx < 3) isMatchxME4103 = true;
          if (Dx < 2) isMatchxME4102 = true;
          if (Dx < 1) isMatchxME4101 = true;
          if (Dy < 15) isMatchyME4115 = true;
          if (Dy < 14) isMatchyME4114 = true;
          if (Dy < 13) isMatchyME4113 = true;
          if (Dy < 12) isMatchyME4112 = true;
          if (Dy < 11) isMatchyME4111 = true;
          if (Dy < 10) isMatchyME4110 = true;
          if (Dy < 9) isMatchyME4109 = true;
          if (Dy < 8) isMatchyME4108 = true;
          if (Dy < 7) isMatchyME4107 = true;
          if (Dy < 6) isMatchyME4106 = true;
          if (Dy < 5) isMatchyME4105 = true;
          if (Dy < 4) isMatchyME4104 = true;
          if (Dy < 3) isMatchyME4103 = true;
          if (Dy < 2) isMatchyME4102 = true;
          if (Dy < 1) isMatchyME4101 = true;
        }

/*
        if (Dx < 10){
          bool dupl = false;
          for (std::vector<GlobalPoint>::const_iterator rpcMatching = rpcMatched_gp.begin(); rpcMatching != rpcMatched_gp.end(); rpcMatching = true){         
            GlobalPoint a = *rpcMatching; 
            if (gp_rpc.x() == a.x() && gp_rpc.y() == a.y() && gp_rpc.z() == a.z()) dupl = true;
          }
          if ( dupl != true ){
            if (csc_id.station() == 3 && csc_id.ring() == 1 && rpcid.station() == 3 && rpcid.ring() == 1){
              MatchedME31++;
              rpcMatched_gp.push_back(gp_rpc);
            }
            if (csc_id.station() == 4 && csc_id.ring() == 1 && rpcid.station() == 4 && rpcid.ring() == 1){
              MatchedME41++;
              rpcMatched_gp.push_back(gp_rpc);
            }
          }
        }

        if ( ismatched != tmp_matched ) break; //move on to next lct


        int kRoll  = rpcid.roll();
        int kSubsector  = rpcid.subsector();
        int kRegion  = rpcid.region();
        int kStation = rpcid.station();
        int kRing = rpcid.ring();
        int kSector = rpcid.sector();
        int kLayer = rpcid.layer();
        int bx = (*rpcIt).BunchX();
        int clSize = (*rpcIt).clusterSize();
    
        cout << "I'm here in RPC Region: " << kRegion <<
                              " Station: " << kStation <<
                              " Ring: "    << kRing <<
                              " Sector: "  << kSector <<
                            " Subsector: " << kSubsector <<
                                 " Roll: " << kRoll <<
                                " Layer: " << kLayer <<
                           "\tbx,clSize: " << bx << ", " << clSize << endl;    
        cout << "number of RPCRecHits: " << nRPC << endl;
    
        int nrechitRE31 = 0;
        int nrechitRE41 = 0;
        if (kRegion != 0){
          if (kStation == 3 && kRing == 1) nrechitRE31++;    
          else if(kStation == 4 && kRing == 1) nrechitRE41++;
    
        }
    
    */
      }//RPCRecHit loop

      for (int i = 0; i < 15; i++){
        for (int j = 0; j < 15; j++){
         
          if (b_cscBX == 6){
            if(isMatchME31[i][j]) ME31[i][j]++;
            if(isMatchME41[i][j]) ME41[i][j]++;
          }
        }
      }
        
      if (isMatchxME3115) xME3115++;
      if (isMatchxME3114) xME3114++;
      if (isMatchxME3113) xME3113++;
      if (isMatchxME3112) xME3112++;
      if (isMatchxME3111) xME3111++;
      if (isMatchxME3110) xME3110++;
      if (isMatchxME3109) xME3109++;
      if (isMatchxME3108) xME3108++;
      if (isMatchxME3107) xME3107++;
      if (isMatchxME3106) xME3106++;
      if (isMatchxME3105) xME3105++;
      if (isMatchxME3104) xME3104++;
      if (isMatchxME3103) xME3103++;
      if (isMatchxME3102) xME3102++;
      if (isMatchxME3101) xME3101++;

      if (isMatchyME3115) yME3115++;
      if (isMatchyME3114) yME3114++;
      if (isMatchyME3113) yME3113++;
      if (isMatchyME3112) yME3112++;
      if (isMatchyME3111) yME3111++;
      if (isMatchyME3110) yME3110++;
      if (isMatchyME3109) yME3109++;
      if (isMatchyME3108) yME3108++;
      if (isMatchyME3107) yME3107++;
      if (isMatchyME3106) yME3106++;
      if (isMatchyME3105) yME3105++;
      if (isMatchyME3104) yME3104++;
      if (isMatchyME3103) yME3103++;
      if (isMatchyME3102) yME3102++;
      if (isMatchyME3101) yME3101++;

      if (isMatchxME4115) xME4115++;
      if (isMatchxME4114) xME4114++;
      if (isMatchxME4113) xME4113++;
      if (isMatchxME4112) xME4112++;
      if (isMatchxME4111) xME4111++;
      if (isMatchxME4110) xME4110++;
      if (isMatchxME4109) xME4109++;
      if (isMatchxME4108) xME4108++;
      if (isMatchxME4107) xME4107++;
      if (isMatchxME4106) xME4106++;
      if (isMatchxME4105) xME4105++;
      if (isMatchxME4104) xME4104++;
      if (isMatchxME4103) xME4103++;
      if (isMatchxME4102) xME4102++;
      if (isMatchxME4101) xME4101++;

      if (isMatchyME4115) yME4115++;
      if (isMatchyME4114) yME4114++;
      if (isMatchyME4113) yME4113++;
      if (isMatchyME4112) yME4112++;
      if (isMatchyME4111) yME4111++;
      if (isMatchyME4110) yME4110++;
      if (isMatchyME4109) yME4109++;
      if (isMatchyME4108) yME4108++;
      if (isMatchyME4107) yME4107++;
      if (isMatchyME4106) yME4106++;
      if (isMatchyME4105) yME4105++;
      if (isMatchyME4104) yME4104++;
      if (isMatchyME4103) yME4103++;
      if (isMatchyME4102) yME4102++;
      if (isMatchyME4101) yME4101++;

      a_ME31NDigis = b_ME31NDigis;
      a_ME41NDigis = b_ME41NDigis;

    }//CSCLCT loop

    if (b_ME31NDigis == 4 && ismatched == 2) a_ME31NDigis=a_ME31NDigis-2;
    if (b_ME41NDigis == 4 && ismatched == 2) a_ME41NDigis=a_ME41NDigis-2;

    if ( b_ME31NDigis != 0 ) h_ME31NDigis->Fill(b_ME31NDigis);
    if ( b_ME41NDigis != 0 ) h_ME41NDigis->Fill(b_ME41NDigis);

    if (bx_ME31NDigis != 0) h_ME31NDigis0->Fill(bx_ME31NDigis);
    if (bx_ME41NDigis != 0) h_ME41NDigis0->Fill(bx_ME41NDigis);

    if ( a_ME31NDigis != 0 ) h_ME31NDigis_a->Fill(a_ME31NDigis);
    if ( a_ME41NDigis != 0 ) h_ME41NDigis_a->Fill(a_ME41NDigis);

    h_ME31NDigis0_a->Fill(a_ME31NDigis);
    h_ME41NDigis0_a->Fill(a_ME41NDigis);

    h_S3NDigis->Fill(b_S3NDigis);
    h_S4NDigis->Fill(b_S4NDigis);

    NRecHits->Fill(nRPC);
    h_S3NRecHits->Fill(b_S3NRecHits);
    h_S4NRecHits->Fill(b_S4NRecHits);
    h_RE31NRecHits->Fill(b_RE31NRecHits);
    h_RE41NRecHits->Fill(b_RE41NRecHits);

/*
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
      cout << "I'm here in CSC: endcap: "    << b_CSCendcap <<
                             " station: "    << b_CSCstation <<
                             " sector: "     << b_CSCsector <<
                             " subsector: "  << b_CSCsubsector <<
                             " strip: "      << b_CSCstrip <<
                             " wire: "       << b_CSCkeyWire << endl;
      cout << "BX = " << b_cscBX << endl;
      cout << "and I'm with CSCDetId: endcap: "    << csc_id.endcap() <<
                                   " station: "    << csc_id.station() <<
                                   " ring: "     << csc_id.ring() <<
                                   " chamber: "  << csc_id.chamber() <<
                                   " layer: "      << csc_id.layer() << endl;
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
    cout << "ME31 NDigis: " << b_ME31NDigis << " ME41 NDigis: " << b_ME41NDigis << endl;   
*/
    //fill histo
//      Ndigis->Fill(b_numberofDigis);
//      if (b_fNDigis != 0) fNdigis->Fill(b_fNDigis);
//      if (b_bNDigis != 0) bNdigis->Fill(b_bNDigis);
//      if (b_ME31NDigis != 0) ME31NDigis->Fill(b_ME31NDigis);
//      if (b_ME41NDigis != 0) ME41NDigis->Fill(b_ME41NDigis);

  }//CSCChamberloop

  tree->Fill();
  EventInfo->Fill(1.5);

}

void 
CSCExtrapoltoRPC::beginJob()
{

  tree->Branch("EVENT", &b_EVENT, "EVENT/i");
  tree->Branch("RUN"  , &b_RUN  , "RUN/i");
  tree->Branch("LUMI" , &b_LUMI , "LUMI/i");

  tree->Branch("numberofDigis" , &b_numberofDigis , "numberofDigis/i");
  tree->Branch("ME31NDigis" , &b_ME31NDigis , "ME31NDigis/i");
  tree->Branch("ME41NDigis" , &b_ME41NDigis , "ME41NDigis/i");

  tree->Branch("RE31NRecHits" , &b_RE31NRecHits , "RE31NRecHits/i");
  tree->Branch("RE41NRecHits" , &b_RE41NRecHits , "RE41NRecHits/i");

  tree->Branch("cscBX" , &b_cscBX , "cscBX/i");
  tree->Branch("rpcBX" , &b_rpcBX , "rpcBX/i");

  xME3115 = xME3114 = xME3113 = xME3112 = xME3111 = xME3110 = xME3109 = xME3108 = xME3107 = xME3106 = xME3105 = xME3104 = xME3103 = xME3102 = xME3101 = 0;
  xME4115 = xME4114 = xME4113 = xME4112 = xME4111 = xME4110 = xME4109 = xME4108 = xME4107 = xME4106 = xME4105 = xME4104 = xME4103 = xME4102 = xME4101 = 0;
  yME3115 = yME3114 = yME3113 = yME3112 = yME3111 = yME3110 = yME3109 = yME3108 = yME3107 = yME3106 = yME3105 = yME3104 = yME3103 = yME3102 = yME3101 = 0;
  yME4115 = yME4114 = yME4113 = yME4112 = yME4111 = yME4110 = yME4109 = yME4108 = yME4107 = yME4106 = yME4105 = yME4104 = yME4103 = yME4102 = yME4101 = 0;

  b_ME31NDigis_Total = b_ME41NDigis_Total = 0;

  for (int i=0; i<15; i++){
    for (int j=0; j<15; j++){
      ME31[i][i] = 0;
      ME41[i][i] = 0;
    }
  }

}

void 
CSCExtrapoltoRPC::endJob() 
{
  if (b_ME31NDigis_Total != 0 ){
    for (int i=0; i<15; i++){
      for (int j=0; j< 15; j++){
        h_xNMatchedME31->SetBinContent(i+1, ME31[i][14]/b_ME31NDigis_Total*100);
        h_yNMatchedME31->SetBinContent(j+1, ME31[14][j]/b_ME31NDigis_Total*100);
      
        h_xNMatchedME41->SetBinContent(i+1, ME41[i][14]/b_ME41NDigis_Total*100);
        h_yNMatchedME41->SetBinContent(j+1, ME41[14][j]/b_ME41NDigis_Total*100);

        h_MatchedME31->SetBinContent(i+1,j+1,ME31[i][j]/b_ME31NDigis_Total*100);
        h_MatchedME41->SetBinContent(i+1,j+1,ME41[i][j]/b_ME41NDigis_Total*100);

      }
    }
  }

  

/*
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(15, ME31[]/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(14, xME3114/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(13, xME3113/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(12, xME3112/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(11, xME3111/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(10, xME3110/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(9, xME3109/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(8, xME3108/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(7, xME3107/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(6, xME3106/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(5, xME3105/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(4, xME3104/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(3, xME3103/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(2, xME3102/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_xNMatchedME31->SetBinContent(1, xME3101/b_ME31NDigis_Total*100);

  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(15, xME4115/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(14, xME4114/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(13, xME4113/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(12, xME4112/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(11, xME4111/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(10, xME4110/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(9, xME4109/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(8, xME4108/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(7, xME4107/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(6, xME4106/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(5, xME4105/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(4, xME4104/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(3, xME4103/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(2, xME4102/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_xNMatchedME41->SetBinContent(1, xME4101/b_ME41NDigis_Total*100);

  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(15, yME3115/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(14, yME3114/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(13, yME3113/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(12, yME3112/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(11, yME3111/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(10, yME3110/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(9, yME3109/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(8, yME3108/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(7, yME3107/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(6, yME3106/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(5, yME3105/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(4, yME3104/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(3, yME3103/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(2, yME3102/b_ME31NDigis_Total*100);
  if (b_ME31NDigis_Total != 0 ) h_yNMatchedME31->SetBinContent(1, yME3101/b_ME31NDigis_Total*100);

  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(15, yME4115/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(14, yME4114/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(13, yME4113/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(12, yME4112/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(11, yME4111/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(10, yME4110/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(9, yME4109/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(8, yME4108/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(7, yME4107/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(6, yME4106/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(5, yME4105/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(4, yME4104/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(3, yME4103/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(2, yME4102/b_ME41NDigis_Total*100);
  if (b_ME41NDigis_Total != 0 ) h_yNMatchedME41->SetBinContent(1, yME4101/b_ME41NDigis_Total*100);
*/
//  if (b_ME31NDigis_Total != 0 && b_ME41NDigis_Total != 0) h_Matched->()

  tree->Fill();  

}

void
CSCExtrapoltoRPC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(CSCExtrapoltoRPC);
