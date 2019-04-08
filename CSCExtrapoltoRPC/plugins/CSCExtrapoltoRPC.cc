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
    TH1D *RE31Nrechit;
    TH1D *RE41Nrechit;

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
//CSCExtrapoltoRPC::getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& rpcIt) const{
CSCExtrapoltoRPC::getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& recHit) const{

  RPCDetId rpcid = RPCDetId(rpcId);
//  const LocalPoint& rpc_lp = rpcIt.localPosition();
  const LocalPoint& rpc_lp = recHit.localPosition();
  const GlobalPoint& rpc_gp = rpcGeo->idToDet(rpcid)->surface().toGlobal(rpc_lp);
 
  return rpc_gp;

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
  Nrechit->GetYaxis()->SetTitle("Number of rechits");

  RE31Nrechit = fs->make<TH1D>("RE31Nrechit", "number of rechits per chamber (RE3/1)", 10, 0, 10);
  RE31Nrechit->GetXaxis()->SetTitle("Number of rechit per chamber");
  RE31Nrechit->GetYaxis()->SetTitle("Number of rechits");
   
  RE41Nrechit = fs->make<TH1D>("RE41Nrechit", "number of rechits per chamber (RE4/1)", 10, 0, 10);
  RE41Nrechit->GetXaxis()->SetTitle("Number of rechit per chamber");
  RE41Nrechit->GetYaxis()->SetTitle("Number of rechits");
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
  b_CSCendcap = b_CSCstation = b_CSCsubsector = b_CSCsector = b_CSCstrip = b_CSCkeyWire = 0;

  b_Trknmb = b_cscBX = b_cscId = 0;

  b_EVENT  = iEvent.id().event();
  b_RUN    = iEvent.id().run();
  b_LUMI   = iEvent.id().luminosityBlock();

  std::cout << "\nNew event" << endl;
    for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc=corrlcts.product()->begin(); csc!=corrlcts.product()->end(); csc++){  
      CSCCorrelatedLCTDigiCollection::Range range1 = corrlcts.product()->get((*csc).first);
//       cout << "\nwe are in ranage 1 " << endl;
      b_numberofDigis = b_fNDigis = b_bNDigis = b_ME31NDigis = b_ME41NDigis = 0;

      for(CSCCorrelatedLCTDigiCollection::const_iterator lct=range1.first; lct!=range1.second; lct++){

        const CSCDetId csc_id((*csc).first.rawId());

        GlobalPoint gp_cscint(getCSCGlobalPosition(csc_id, *lct));  
        if (csc_id.station() == 3 && csc_id.ring() == 1){
          b_ME31NDigis++;
          cout << "CSCGlobalposition in ME3/1" << gp_cscint << "BX: " << lct->getBX() << endl;
        }
        else if (csc_id.station() == 4 && csc_id.ring() == 1){
          b_ME41NDigis++;
          cout << "CSCGlobalposition in ME4/1" << gp_cscint << "BX: " << lct->getBX() << endl;
        }
//        cout << "Globalfunction intersect: " << gp_cscint << endl;
      }
      
      Ndigis->Fill(b_numberofDigis);
      ME31NDigis->Fill(b_ME31NDigis);
      ME41NDigis->Fill(b_ME41NDigis);
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

      //fill histo
      Ndigis->Fill(b_numberofDigis);
      if (b_fNDigis != 0) fNdigis->Fill(b_fNDigis);
      if (b_bNDigis != 0) bNdigis->Fill(b_bNDigis);
      if (b_ME31NDigis != 0) ME31NDigis->Fill(b_ME31NDigis);
      if (b_ME41NDigis != 0) ME41NDigis->Fill(b_ME41NDigis);
*/
  }

  int numberofRecHits = 0;
  int nrechitRE31 = 0;
  int nrechitRE41 = 0;
  int nRPC = 0;
  for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {
    RPCDetId rpcId = (*rpcIt).rpcId();
    typedef std::pair<RPCRecHitCollection::const_iterator, RPCRecHitCollection::const_iterator> rangeRecHits;
    rangeRecHits recHitCollection = rpcRecHits->get(rpcId);
    RPCRecHitCollection::const_iterator recHit;

    numberofRecHits = nrechitRE31 = nrechitRE41 = 0;
    
    for (recHit = recHitCollection.first; recHit != recHitCollection.second; recHit++){
//    bool isNewChamber = true;
      nRPC++;
      numberofRecHits++;
      RPCDetId rpcid = (*rpcIt).rpcId();
//      GlobalPoint gp_rpc(getRPCGlobalPosition(rpcid, *rpcIt));
      GlobalPoint gp_rpc(getRPCGlobalPosition(rpcid, *recHit));

      if (rpcid.region() == -1 || rpcid.region() == 1){
        if (rpcid.station() == 3 && rpcid.ring() == 1){
          nrechitRE31++;
          cout << "RPCGlobalposition in ME3/1" << gp_rpc << "BX/clSize: " << (*rpcIt).BunchX()  << "/" << (*rpcIt).clusterSize() << endl;
        }
        else if (rpcid.station() == 4 && rpcid.ring() == 1){
          nrechitRE41++;
          cout << "RPCGlobalposition in ME4/1" << gp_rpc << "BX/clSize: " << (*rpcIt).BunchX() << "/" << (*rpcIt).clusterSize() << endl;
        }
      }

    }
    
    Nrechit->Fill(numberofRecHits);
    RE31Nrechit->Fill(nrechitRE31);
    RE41Nrechit->Fill(nrechitRE41);
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

    cout << "RE31: " << nrechitRE31 << " RE41: " << nrechitRE41 << endl;
    if (nrechitRE31 != 0) RE31Nrechit->Fill(nrechitRE31);
    if (nrechitRE41 != 0) RE41Nrechit->Fill(nrechitRE41);
*/


   
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
