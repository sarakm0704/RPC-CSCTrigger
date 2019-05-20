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
#include <L1Trigger/CSCCommonTrigger/interface/CSCConstants.h>

using namespace edm;
using namespace std;

class CSCExtrapoltoRPC : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit CSCExtrapoltoRPC(const edm::ParameterSet&);
    ~CSCExtrapoltoRPC();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    //https://github.com/cms-sw/cmssw/blob/02d4124c0b6690287fd88e9a8ff650aea254412e/L1Trigger/L1TMuon/interface/GeometryTranslator.h
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
    TH1D *h_RE31NRecHits0;
    TH1D *h_RE41NRecHits0;

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
    int bx_RE31NRecHits;
    int bx_RE41NRecHits;

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

    double ME31[25][25];
    double ME41[25][25];

    bool isMatchME31[25][25];
    bool isMatchME41[25][25];

    int EventNum;

    edm::EDGetTokenT<MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>> corrlctsToken_;
    edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;

    edm::ESHandle<CSCGeometry> cscGeo;
    edm::ESHandle<RPCGeometry> rpcGeo;

    edm::Handle<RPCRecHitCollection> rpcRecHits;

    //GlobalPoint getGlobalPosition(unsigned int rawId, const CSCCorrelatedLCTDigi& lct) const;
    GlobalPoint getCSCGlobalPosition(CSCDetId rawId, const CSCCorrelatedLCTDigi& lct) const;
    GlobalPoint getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& rpcIt) const;
    std::pair<RPCRecHitCollection::const_iterator, float*> matchingRPC(CSCDetId rawId, const CSCCorrelatedLCTDigi& lct, float dx_cutoff, float dy_cutoff) const;
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

//RPCRecHitCollection::const_iterator
std::pair<RPCRecHitCollection::const_iterator, float*>
CSCExtrapoltoRPC::matchingRPC(CSCDetId rawId, const CSCCorrelatedLCTDigi& lct, float dx_cutoff, float dy_cutoff) const{
//const float & makes this faster?

  GlobalPoint gp_cscint(0.0,0.0,0.0);
  gp_cscint = getCSCGlobalPosition(rawId, lct);

  float min_distance = std::numeric_limits<float>::max();
  RPCRecHitCollection::const_iterator rpc_hit_matched;
  float min_DxDy[2] = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
  
  for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {
    RPCDetId rpcid = (RPCDetId)(*rpcIt).rpcId();

    GlobalPoint gp_rpc(0.0,0.0,0.0);
    gp_rpc = getRPCGlobalPosition(rpcid, *rpcIt);

    if (rpcid.region() == 0) continue; //skip the barrels

    float Dx = abs(gp_rpc.x()-gp_cscint.x());
    float Dy = abs(gp_rpc.y()-gp_cscint.y());
    float distance = sqrt(Dx*Dx + Dy*Dy);

    if (rawId.station() == 3 && rawId.ring() == 1 && rpcid.station() == 3 && rpcid.ring() == 1){
      if (Dx < dx_cutoff && Dy < dy_cutoff && distance < min_distance){
        min_DxDy[0]=Dx;
        min_DxDy[1]=Dy;
        min_distance = distance;
        rpc_hit_matched = rpcIt;
      }
    }
  }
//  return rpc_hit_matched;
  return std::make_pair(rpc_hit_matched,min_DxDy);
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
  Ndigis->GetYaxis()->SetTitle("Number of chamber");

  h_xNMatchedME31 = fs->make<TH1D>("h_xNMatchedME31", "Matching Efficiency in ME3/1", 25, 0, 25);
  h_xNMatchedME31->GetXaxis()->SetTitle("X cutoff (cm)");
  h_xNMatchedME31->GetYaxis()->SetTitle("Matched (%)");

  h_xNMatchedME41 = fs->make<TH1D>("h_xNMatchedME41", "Matching Efficiency in ME4/1", 25, 0, 25);
  h_xNMatchedME41->GetXaxis()->SetTitle("X cutoff (cm)");
  h_xNMatchedME41->GetYaxis()->SetTitle("Matched (%)");

  h_yNMatchedME31 = fs->make<TH1D>("h_yNMatchedME31", "Matching Efficiency in ME3/1", 25, 0, 25);
  h_yNMatchedME31->GetXaxis()->SetTitle("Y cutoff (cm)");
  h_yNMatchedME31->GetYaxis()->SetTitle("Matched (%)");

  h_yNMatchedME41 = fs->make<TH1D>("h_yNMatchedME41", "Matching Efficiency in ME4/1", 25, 0, 25);
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

  h_RE31NRecHits0 = fs->make<TH1D>("h_RE31NRecHits0", "number of rechit per chamber (RE3/1)", 20, 0, 20);
  h_RE31NRecHits0->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_RE31NRecHits0->GetYaxis()->SetTitle("Number of chamber");
   
  h_RE41NRecHits0 = fs->make<TH1D>("h_RE41NRecHits0", "number of rechit per chamber (RE4/1)", 20, 0, 20);
  h_RE41NRecHits0->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_RE41NRecHits0->GetYaxis()->SetTitle("Number of chamber");

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

  h_MatchedME31 = fs->make<TH2D>("h_MatchedME31", "Matching efficiency in ME31", 25, 0, 25, 25, 0, 25);
  h_MatchedME31->GetXaxis()->SetTitle("X cutoff (cm)");
  h_MatchedME31->GetYaxis()->SetTitle("Y cutoff (cm)");

  h_MatchedME41 = fs->make<TH2D>("h_MatchedME41", "Matching efficiency in ME41", 25, 0, 25, 25, 0, 25);
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

//  edm::Handle<RPCRecHitCollection> rpcRecHits;
  iEvent.getByToken(rpcRecHitsToken_, rpcRecHits);

  if (!corrlcts.isValid()) {
    edm::LogInfo("DataNotFound") << "can't find CSCCorrleatedLCTDigiCollection with label "<< corrlcts << std::endl;
    return;
  }
  if (!rpcRecHits.isValid()) {
    edm::LogInfo("DataNotFound") << "can't find RPCRecHitDigiCollection with label "<< rpcRecHits << std::endl;
    return;
  }

  b_EVENT = b_RUN = b_LUMI = 0;

  b_EVENT  = iEvent.id().event();
  b_RUN    = iEvent.id().run();
  b_LUMI   = iEvent.id().luminosityBlock();

  std::vector<GlobalPoint> rpcMatched_gp;

  EventNum++;

  cout << "\nNew event" << endl;
  cout << "Event " << EventNum << endl;

  b_S3NRecHits = b_S4NRecHits = 0;
  bx_RE31NRecHits = bx_RE41NRecHits =0;
  nRPC = b_RE31NRecHits = b_RE41NRecHits = 0;

  for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {

    RPCDetId rpcid = (RPCDetId)(*rpcIt).rpcId();
    nRPC++;
    if(rpcid.station() == 3) b_S3NRecHits++;
    if(rpcid.station() == 4) b_S4NRecHits++;

    if(rpcid.station() == 3 && rpcid.ring() == 1) b_RE31NRecHits++;
    if(rpcid.station() == 4 && rpcid.ring() == 1) b_RE41NRecHits++;
    if((*rpcIt).BunchX() == 0 && rpcid.station() == 3 && rpcid.ring() == 1) bx_RE31NRecHits++;
    if((*rpcIt).BunchX() == 0 && rpcid.station() == 4 && rpcid.ring() == 1) bx_RE41NRecHits++;

    //print globalposition
    GlobalPoint gp_rpcsample(0.0,0.0,0.0);
    gp_rpcsample = getRPCGlobalPosition(rpcid, *rpcIt);
//    cout << "RPCGlobalPosition" << gp_rpcsample << endl; 
 
  }


  for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc=corrlcts.product()->begin(); csc!=corrlcts.product()->end(); csc++){  
    CSCCorrelatedLCTDigiCollection::Range range1 = corrlcts.product()->get((*csc).first);
    b_numberofDigis = b_ME11NDigis  = b_ME21NDigis = b_ME31NDigis = b_ME41NDigis = 0;
    bx_ME31NDigis = bx_ME41NDigis =0;
    b_cscBX = b_rpcBX = 0;

    a_ME31NDigis = a_ME41NDigis = 0;
    int ismatched = 0;

    b_S3NDigis = b_S4NDigis = 0;

    rpcMatched_gp.clear();
    GlobalPoint gp_cscint;

    range1.first++;
//    range1.first++;
//    range1.first++;
//    range1.first++;
    if (range1.first != range1.second) continue; // check that there are two digis in the chamber, there is probably a better way but it works...
//    range1.first--;
//    range1.first--;
//    range1.first--;
    range1.first--;
    int ndigi = 0;

    for(CSCCorrelatedLCTDigiCollection::const_iterator lct=range1.first; lct!=range1.second; lct++){
      ndigi++;
      const CSCDetId csc_id((*csc).first.rawId());

      b_cscBX = lct->getBX();

      gp_cscint = GlobalPoint(0.0,0.0,0.0);
      gp_cscint = getCSCGlobalPosition(csc_id, *lct);
      cout << "CSCGlobalPoint" << gp_cscint << endl;
      double xslope = gp_cscint.x()/gp_cscint.z();
      double yslope = gp_cscint.y()/gp_cscint.z();

      b_numberofDigis++;
      if(abs(gp_cscint.eta())<1.9) continue;

      if (csc_id.station() == 3) b_S3NDigis++;
      if (csc_id.station() == 4) b_S4NDigis++;

      if (csc_id.station() == 3 && csc_id.ring() == 1){

        b_ME31NDigis = b_ME31NDigis + 1;
        b_ME31NDigis_Total = b_ME31NDigis_Total + 1;
        if (b_cscBX == 6) bx_ME31NDigis++;
      }
      else if (csc_id.station() == 4 && csc_id.ring() == 1){

        b_ME41NDigis = b_ME41NDigis + 1;
        b_ME41NDigis_Total = b_ME41NDigis_Total + 1;
        if (b_cscBX == 6) bx_ME41NDigis++;
      }

      for (int i=0; i<25; i++){
        for (int j=0; j<25; j++){
          isMatchME31[i][j] = false;
          isMatchME41[i][j] = false;
        }
      }

      for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {

        RPCDetId rpcid = (RPCDetId)(*rpcIt).rpcId();

        GlobalPoint gp_rpc(0.0,0.0,0.0);
        gp_rpc = getRPCGlobalPosition(rpcid, *rpcIt);
        cout << "RPCGlobalPoint" << gp_rpc << endl;

        if (gp_rpc.z() * gp_cscint.z() < 0 ) continue;

        double dz = gp_rpc.z() - gp_cscint.z();

        double dx = dz*xslope;
        double dy = dz*yslope;

        GlobalPoint gp_transcsc(gp_cscint.x()+dx, gp_cscint.y()+dy, gp_rpc.z());
        cout << "translated CSCGlobalPoint" << gp_transcsc << endl;

        LocalPoint lp_extrapol = rpcGeo->idToDet(rpcid)->surface().toLocal(gp_transcsc);

        //phi relation
//        float rpcphi=0;
//        float cscphi=0;
      
//        (CenterPointRollGlobal.barePhi()<0)? 
//          rpcphi = 2*3.149092+CenterPointRollGlobal.barePhi():rpcphi=CenterPointRollGlobal.barePhi();
      
//        (CenterPointCSCGlobal.barePhi()<0)?
//          cscphi = 2*3.1490926536+CenterPointCSCGlobal.barePhi():cscphi=CenterPointCSCGlobal.barePhi();

//        cout << "dphi" << abs(rpcphi-cscphi) << endl;
//        if (abs(rpcphi-cscphi) > 0.4) continue;

        if (rpcid.region() == 0) continue; //skip the barrels

        if (gp_rpc.x() == 0 && gp_rpc.y() == 0 && gp_rpc.z() == 0 ) continue;
        if (gp_cscint.x() == 0 && gp_cscint.y() == 0 && gp_cscint.z() == 0 ) continue;

        //global distance
//        float Dx = abs(gp_rpc.x()-gp_cscint.x());
//        float Dy = abs(gp_rpc.y()-gp_cscint.y());

        //local distance
        LocalPoint lp_rpc(0.0,0.0,0.0);
        lp_rpc = (*rpcIt).localPosition();
        float Dx = abs(lp_rpc.x()-lp_extrapol.x());
        float Dy = abs(lp_rpc.y()-lp_extrapol.y());

//        cout << "Dx= " << Dx << "Dy= " << Dy << endl;

        if (csc_id.station() == 3 && csc_id.ring() == 1 && rpcid.station() == 3 && rpcid.ring() == 1){
          for (int i = 0; i < 25; i++){
            for (int j = 0; j < 25; j++){
              if (Dx < i+1 && Dy < j+1){
                isMatchME31[i][j] = true;
              }
            }
          }
        }
        if (csc_id.station() == 4 && csc_id.ring() == 1 && rpcid.station() == 4 && rpcid.ring() == 1){
          for (int i = 0; i < 25; i++){
            for (int j = 0; j < 25; j++){
              if (Dx < i+1 && Dy < j+1){
                isMatchME41[i][j] = true;
              }
            }
          }
        }

      }//RPCRecHit loop

      for (int i = 0; i < 25; i++){
        for (int j = 0; j < 25; j++){
         
          if(isMatchME31[i][j]) ME31[i][j]++;
          if(isMatchME41[i][j]) ME41[i][j]++;
          
        }
      }

      a_ME31NDigis = b_ME31NDigis;
      a_ME41NDigis = b_ME41NDigis;

    }//CSCLCT loop
//    cout << ndigi << endl;

    if (b_ME31NDigis == 4 && ismatched == 2) a_ME31NDigis=a_ME31NDigis-2;
    if (b_ME41NDigis == 4 && ismatched == 2) a_ME41NDigis=a_ME41NDigis-2;

    if (b_ME31NDigis != 0) h_ME31NDigis->Fill(b_ME31NDigis);
    if (b_ME41NDigis != 0) h_ME41NDigis->Fill(b_ME41NDigis);

    if (bx_ME31NDigis != 0) h_ME31NDigis0->Fill(bx_ME31NDigis);
    if (bx_ME41NDigis != 0) h_ME41NDigis0->Fill(bx_ME41NDigis);

    if ( a_ME31NDigis != 0 ) h_ME31NDigis_a->Fill(a_ME31NDigis);
    if ( a_ME41NDigis != 0 ) h_ME41NDigis_a->Fill(a_ME41NDigis);

    h_ME31NDigis0_a->Fill(a_ME31NDigis);
    h_ME41NDigis0_a->Fill(a_ME41NDigis);

    if (b_S3NDigis !=0 ) h_S3NDigis->Fill(b_S3NDigis);
    if (b_S4NDigis !=0 ) h_S4NDigis->Fill(b_S4NDigis);

  }//CSCChamber loop

  NRecHits->Fill(nRPC);
  h_S3NRecHits->Fill(b_S3NRecHits);
  h_S4NRecHits->Fill(b_S4NRecHits);
  h_RE31NRecHits->Fill(b_RE31NRecHits);
  h_RE41NRecHits->Fill(b_RE41NRecHits);
  h_RE31NRecHits0->Fill(bx_RE31NRecHits);
  h_RE41NRecHits0->Fill(bx_RE41NRecHits);

  tree->Fill();
  EventInfo->Fill(1.5);

}

void 
CSCExtrapoltoRPC::beginJob()
{

  EventNum = 0;

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

  b_ME31NDigis_Total = b_ME41NDigis_Total = 0;

  for (int i=0; i<25; i++){
    for (int j=0; j<25; j++){
      ME31[i][j] = 0;
      ME41[i][j] = 0;
    }
  }

}

void 
CSCExtrapoltoRPC::endJob() 
{

  if (b_ME31NDigis_Total != 0){
    for (int i=0; i<25; i++){
      h_xNMatchedME31->SetBinContent(i+1, ME31[i][24]/b_ME31NDigis_Total*100);
      h_xNMatchedME41->SetBinContent(i+1, ME41[i][24]/b_ME41NDigis_Total*100);
    }
  }
  if (b_ME41NDigis_Total != 0){
    for (int j=0; j< 25; j++){
      h_yNMatchedME31->SetBinContent(j+1, ME31[24][j]/b_ME31NDigis_Total*100);
      h_yNMatchedME41->SetBinContent(j+1, ME41[24][j]/b_ME41NDigis_Total*100);
    }
  }
    
  if (b_ME31NDigis_Total != 0 && b_ME41NDigis_Total != 0){
    for (int i=0; i<25; i++){
      for (int j=0; j<25; j++){
        h_MatchedME31->SetBinContent(i+1,j+1,ME31[i][j]/b_ME31NDigis_Total*100);
        h_MatchedME41->SetBinContent(i+1,j+1,ME41[i][j]/b_ME41NDigis_Total*100);  
      }
    } 
  }

  tree->Fill();  

  cout << "matched ratio at 12cm * 16cm ME31 " << ME31[12][16]/b_ME31NDigis_Total*100 << endl;
  cout << "matched ratio at 12cm * 16cm ME41 " << ME41[12][16]/b_ME41NDigis_Total*100 << endl;

}

void
CSCExtrapoltoRPC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(CSCExtrapoltoRPC);
