// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author everyone

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <vector>
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "PHOSBase/Geometry.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"

#include <iostream>
#include <cmath>

using namespace o2;
using namespace o2::framework;

struct MCTrackSel {

  Configurable<int> nITSClst{"nITSClst",3,"N of clusters on ITS"};
  Configurable<int> nTPCClst1{"nTPCClst1",30,"First border N of clusters on TPC"};
  Configurable<int> nTPCClst2{"nTPCClst2",50,"Second border N of clusters on TPC"};
  Configurable<int> nTPCClst3{"nTPCClst3",70,"Third border N of clusters on TPC"};
  Configurable<int> PtUpper{"PtUpper",100,"Upper limit for Traverse momentum"};


  TrackSelection trackselOptions()
  {
    TrackSelection selectedTracks;
    // selectedTracks.SetTrackType(aod::track::Run2GlobalTrack);
    selectedTracks.SetPtRange(0.1f, 1e10f);
    selectedTracks.SetEtaRange(-0.8f, 0.8f);
    selectedTracks.SetRequireITSRefit(false);  //
    selectedTracks.SetRequireTPCRefit(true);   // Not Sure what these are
    selectedTracks.SetRequireGoldenChi2(true); //
    selectedTracks.SetMinNCrossedRowsTPC(70);
    selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
    selectedTracks.SetMaxChi2PerClusterTPC(4.f);
    selectedTracks.SetRequireHitsInITSLayers(1, {0, 1}); // one hit in any SPD
    selectedTracks.SetMaxChi2PerClusterITS(36.f);
    selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / pow(pt, 1.1f); }); //there is a function that relates the DCA to the PT. Need to figure this out
    selectedTracks.SetMaxDcaZ(2.f);
    return selectedTracks;
  }
  // Histogram registry: an object to hold your histograms
  TrackSelection mySelection; //Initialisation???

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    mySelection = trackselOptions();

    const AxisSpec axisEta{250, 0, 1, "Eta"};
    const AxisSpec axisDCA{250, 0, 0.3, "axisDCA"};
    const AxisSpec axisNoHitDCA{250, 0, 0.3, "axisNoHitDCA"};
    const AxisSpec ITSaxisCounter{250, 0, 11, "Number of Clusters"};
    const AxisSpec TPCaxisCounter{250, 0, 160, "Number of Clusters"};
    // PHI
    const AxisSpec axisPhi{250,0,2*M_PI,"Phi"};
    // ETA
    const AxisSpec axiseta{250,-1.5,1.5,"eta"};

    const AxisSpec axisTPCChi{250, 0, 5, "TPCChi"};
    const AxisSpec axisTPCCR{250, 0, 300, "TPCCr"};
    const AxisSpec axisRatio{250, 0, 2, "Ratio"};
    const AxisSpec axisVertex{100, -0.3, +0.3};

    histos.add("SelectionHistogram", "SelectionHistogram", kTH1F, {axisEta});
    histos.add("DCAHistogram", "DCAHistogram", kTH1F, {axisDCA});
    histos.add("TPCChiHistogram", "TPCChiHistogram", kTH1F, {axisTPCChi});
    histos.add("TPCCrHistogram", "TPCCrHistogram", kTH1F, {axisTPCCR});
    histos.add("RatioHistogram", "RatioHistogram", kTH1F, {axisRatio});
    // ETA
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axiseta});
    histos.add("etaHistogramLOW", "etaHistogramLOW", kTH1F, {axiseta});
    histos.add("etaHistogramMED", "etaHistogramMED", kTH1F, {axiseta});
    histos.add("etaHistogramHIGH", "etaHistogramHIGH", kTH1F, {axiseta});
    
    //PHI
    histos.add("phiHistogram", "phiHistogram", kTH1F, {axisPhi});
    histos.add("phiHistogramLOW", "phiHistogramLOW", kTH1F, {axisPhi});
    histos.add("phiHistogramMED", "phiHistogramMED", kTH1F, {axisPhi});
    histos.add("phiHistogramHIGH", "phiHistogramHIGH", kTH1F, {axisPhi});

    const AxisSpec axisITSRefitChi{250, 0, 5, "ITSRefit"};
    histos.add("ITSRefitChi", "ITSRefitChi", kTH1F, {axisITSRefitChi});

    // TPC crossed rows Chi2
    const AxisSpec axis035TPCChi{250, 0, 5};
    histos.add("035TPCChi", "035TPCChi", kTH1F, {axis035TPCChi});
    histos.add("3570TPCChi", "3570TPCChi", kTH1F, {axis035TPCChi});
    histos.add("70TPCChi", "70TPCChi", kTH1F, {axis035TPCChi});

    // ITS Chi2
    const AxisSpec axis13HitITS{250, 0, 5};

    histos.add("13HitITS", "13HitITS", kTH1F, {axis13HitITS});
    histos.add("NoHitITS", "NoHitITS", kTH1F, {axis13HitITS});
    histos.add("1AllHitITS", "1AllHitITS", kTH1F, {axis13HitITS});
    histos.add("H1370HitITS", "H1370HitITS", kTH1F, {axis13HitITS});
    histos.add("ITSChi", "ITSChi", kTH1F, {axis13HitITS});
    histos.add("ChiLOW","ChiLOW",kTH1F, {axisTPCChi});
    histos.add("ChiMED","ChiMED",kTH1F, {axisTPCChi});
    histos.add("ChiHIGH","ChiHIGH",kTH1F, {axisTPCChi});
    histos.add("ChiBOTH","ChiBOTH",kTH1F, {axisTPCChi});
    histos.add("ChiFIRST","ChiFIRST",kTH1F, {axisTPCChi});
    histos.add("ChiSECOND","ChiSECOND",kTH1F, {axisTPCChi});
    histos.add("ChiNONE","ChiNONE",kTH1F, {axisTPCChi});
    //histos.add("ITSChi","ITSChi",kTH1F, {axisTPCChi});

    // DCA
    const AxisSpec axis13HitITSDCA{250, 0, 0.3};
    histos.add("13HitITSDCA", "13HitITSDCA", kTH1F, {axis13HitITSDCA});
    histos.add("NoHitITSDCA", "NoHitITSDCA", kTH1F, {axis13HitITSDCA});
    histos.add("1AllHitITSDCA", "1AllHitITSDCA", kTH1F, {axis13HitITSDCA});
    histos.add("H1370HitITSDCA", "H1370HitITSDCA", kTH1F, {axis13HitITSDCA});
    histos.add("DCAlow","DCAlow", kTH1F, {axisDCA});
    histos.add("DCAmed","DCAmed", kTH1F, {axisDCA});
    histos.add("DCAhigh","DCAhigh", kTH1F, {axisDCA});
    histos.add("DCABOTH","DCABOTH", kTH1F, {axisDCA});
    histos.add("DCAFIRST","DCAFIRST", kTH1F, {axisDCA});
    histos.add("DCASECOND","DCASECOND", kTH1F, {axisDCA});
    histos.add("DCANONE","DCANONE", kTH1F, {axisDCA});

    // pt
    const AxisSpec axis13HitITSPT{250, 0, 15};
    histos.add("13HitITSPT", "13HitITSPT", kTH1F, {axis13HitITSPT});
    histos.add("NoHitITSPT", "NoHitITSPT", kTH1F, {axis13HitITSPT});
    histos.add("1AllHitITSPT", "1AllHitITSPT", kTH1F, {axis13HitITSPT});
    histos.add("H1370HitITSPT", "H1370HitITSPT", kTH1F, {axis13HitITSPT});
    histos.add("PtHistogram", "PtHistogram", kTH1F, {axis13HitITSPT});
    histos.add("PtHistogramLOW", "PtHistogramLOW", kTH1F, {axis13HitITSPT});
    histos.add("PtHistogramMED", "PtHistogramMED", kTH1F, {axis13HitITSPT});
    histos.add("PtHistogramHIGH", "PtHistogramHIGH", kTH1F, {axis13HitITSPT});
    histos.add("PtHistogramTOTAL", "PtHistogramTOTAL", kTH1F, {axis13HitITSPT});
    histos.add("PtHistogramBOTH", "PtHistogramBOTH", kTH1F, {axis13HitITSPT});
    histos.add("PtHistogramFIRST", "PtHistogramFIRST", kTH1F, {axis13HitITSPT});
    histos.add("PtHistogramSECOND", "PtHistogramSECOND", kTH1F, {axis13HitITSPT});
    histos.add("PtHistogramNONE", "PtHistogramNONE", kTH1F, {axis13HitITSPT});


    // ITSclst vs Pt
    histos.add("VertexPoint", "VertexPoint", kTH2F, {axisVertex,axisVertex});
    histos.add("ITSclstvsPt", "ITSclstvsPt", kTH2F, {ITSaxisCounter,axis13HitITSPT});
    histos.add("ITSchivsPt", "ITSchivsPt", kTH2F, {axis13HitITS,axis13HitITSPT});
    histos.add("TPCchivsPt", "TPCchivsPt", kTH2F, {axisTPCChi,axis13HitITSPT});
    histos.add("PtvsTPCclst", "PtvsTPCclst", kTH2F, {axis13HitITSPT,TPCaxisCounter});
    // DCA vs PT
    histos.add("DCAvsPt", "DCAvsPt", kTH2F, {axisDCA, axis13HitITSPT});
    histos.add("NoHitDCAvsPt", "NoHitDCAvsPt", kTH2F, {axisNoHitDCA, axis13HitITSPT});
    histos.add("TrackSelectDCAvsPt", "TrackSelectDCAvsPt", kTH2F, {axisDCA, axis13HitITSPT});
    histos.add("ITSClst","ITSClst",kTH1F,{ITSaxisCounter});
    histos.add("TPCClst","TPCClst",kTH1F,{TPCaxisCounter});
    histos.add("TPCClst1","TPCClst1",kTH1F,{TPCaxisCounter});
    histos.add("TPCClst2","TPCClst2",kTH1F,{TPCaxisCounter});
    histos.add("TPCClst3","TPCClst3",kTH1F,{TPCaxisCounter});

  }

void process(aod::Collisions const& colls, soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra_001, aod::TracksDCA, aod::TrackSelectionExtension> const& tracks)
  {
    for (auto& coll : colls){
      histos.fill(HIST("VertexPoint"), coll.posX(), coll.posY());
    }
    for (auto& track : tracks) {
      //bool isSelected = mySelection.IsSelected(track);
      histos.fill(HIST("PtHistogramTOTAL"), track.pt());
      //histos.fill(HIST("ITSClst"),track.itsNCls());
      //histos.fill(HIST("TPCClst"),track.tpcNClsCrossedRows());
      
      if (abs(track.eta()) > 0.8 || track.pt() < 0.1 || track.pt() > PtUpper || track.itsChi2NCl() < 0.5 
      || track.itsChi2NCl()>2.0 || track.tpcChi2NCl() < 0.5 || track.tpcChi2NCl() > 2.0)
        continue;

      // TPC Chi2
      histos.fill(HIST("TPCClst"),track.tpcNClsCrossedRows());
      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("phiHistogram"), track.phi());
      //histos.fill(HIST("PtHistogram"), track.pt());
      histos.fill(HIST("ITSChi"), track.itsChi2NCl());
      histos.fill(HIST("ITSClst"),track.itsNCls());
      histos.fill(HIST("ITSclstvsPt") , track.pt(), track.itsNCls());
      histos.fill(HIST("ITSchivsPt") , track.pt(), track.itsChi2NCl());
      histos.fill(HIST("TPCchivsPt") , track.pt(), track.tpcChi2NCl());
      histos.fill(HIST("PtvsTPCclst"), track.pt(), track.tpcNClsCrossedRows());
      //histos.fill(HIST("ITSclstvsPt") , track.itsNCls(), track.pt());
            
      if ((track.itsClsSizeInLayer(0) != 0 || track.itsClsSizeInLayer(1) != 0) && track.itsNCls() >= nITSClst) {

        //histos.fill(HIST("ITSClst"),track.itsNCls());
        histos.fill(HIST("PtHistogram"), track.pt());
        
        
        /*
        if (track.tpcNClsCrossedRows() < 130 && track.tpcNClsCrossedRows() >= nTPCClst1) {
          histos.fill(HIST("035TPCChi"), track.tpcChi2NCl());
          histos.fill(HIST("TPCClst1"),track.tpcNClsCrossedRows());
          histos.fill(HIST("PtHistogramLOW"),track.pt());
          histos.fill(HIST("DCAlow"),track.dcaXY());
          histos.fill(HIST("etaHistogramLOW"),track.eta());
          histos.fill(HIST("phiHistogramLOW"),track.phi());
          histos.fill(HIST("ChiLOW"),track.tpcChi2NCl());

        }
        */

        if (track.tpcNClsCrossedRows() < 130 && track.tpcNClsCrossedRows() >= nTPCClst2) {
          histos.fill(HIST("3570TPCChi"), track.tpcChi2NCl());
          histos.fill(HIST("TPCClst2"),track.tpcNClsCrossedRows());
          histos.fill(HIST("PtHistogramMED"),track.pt());
          histos.fill(HIST("DCAmed"),track.dcaXY());
          histos.fill(HIST("etaHistogramMED"),track.eta());
          histos.fill(HIST("phiHistogramMED"),track.phi());
          histos.fill(HIST("ChiMED"),track.tpcChi2NCl());
        }

        /*
        if (track.tpcNClsCrossedRows() < 130 && track.tpcNClsCrossedRows() >= nTPCClst3) {
          histos.fill(HIST("70TPCChi"), track.tpcChi2NCl());
          histos.fill(HIST("TPCClst3"),track.tpcNClsCrossedRows());
          histos.fill(HIST("PtHistogramHIGH"),track.pt());
          histos.fill(HIST("DCAhigh"),track.dcaXY());
          histos.fill(HIST("etaHistogramHIGH"),track.eta());
          histos.fill(HIST("phiHistogramHIGH"),track.phi());
          histos.fill(HIST("ChiHIGH"),track.tpcChi2NCl());
        }
        */
        /**
        if (track.tpcNClsCrossedRows() >= 20 && track.tpcNClsCrossedRows() <= 70)
        {
          histos.fill(HIST("etaHistogramLOWHALF"), track.eta());
          histos.fill(HIST("phiHistogramLOWHALF"), track.phi());
          histos.fill(HIST("PtHistogramLOWHALF"), track.pt());
          histos.fill(HIST("DCALOWHALF"),track.dcaXY());
          histos.fill(HIST("ChiLOWHALF"),track.tpcChi2NCl());
        }

        if (track.tpcNClsCrossedRows() >= 70 && track.tpcNClsCrossedRows() <= 150)
        {
          histos.fill(HIST("etaHistogramHIGHHALF"), track.eta());
          histos.fill(HIST("phiHistogramHIGHHALF"), track.phi());
          histos.fill(HIST("PtHistogramHIGHHALF"), track.pt());
          histos.fill(HIST("DCAHIGHHALF"),track.dcaXY());
          histos.fill(HIST("ChiHIGHHALF"),track.tpcChi2NCl());
        }
        
        histos.fill(HIST("TPCChiHistogram"), track.tpcChi2NCl());

        histos.fill(HIST("TPCCrHistogram"), track.tpcNClsCrossedRows());

        histos.fill(HIST("RatioHistogram"), track.tpcCrossedRowsOverFindableCls());
        */
      }
      // ITS BOTH hit on first 2 layers
      if ((track.itsClsSizeInLayer(0) != 0 && track.itsClsSizeInLayer(1) != 0) && track.itsNCls() >= nITSClst) {
        
        if (track.tpcNClsCrossedRows() < 130 && track.tpcNClsCrossedRows() >= nTPCClst2) {
          histos.fill(HIST("PtHistogramBOTH"), track.pt());
          histos.fill(HIST("DCABOTH"), track.dcaXY());
          histos.fill(HIST("ChiBOTH"), track.tpcChi2NCl());
        }
      
      }
      // ITS NO hit on NONE layer
      if (track.itsNCls() >= nITSClst) {
        
        if (track.tpcNClsCrossedRows() < 130 && track.tpcNClsCrossedRows() >= nTPCClst2) {
          histos.fill(HIST("PtHistogramNONE"), track.pt());
          histos.fill(HIST("DCANONE"), track.dcaXY());
          histos.fill(HIST("ChiNONE"), track.tpcChi2NCl());
        }
      
      }
      // ITS ONLY 1st hit
      if (track.itsClsSizeInLayer(0) != 0 && track.itsNCls() >= nITSClst) {
        
        if (track.tpcNClsCrossedRows() < 130 && track.tpcNClsCrossedRows() >= nTPCClst2) {
          histos.fill(HIST("PtHistogramFIRST"), track.pt());
          histos.fill(HIST("DCAFIRST"), track.dcaXY());
          histos.fill(HIST("ChiFIRST"), track.tpcChi2NCl());
        }
      
      }
      // ITS ONLY 2nd hit
      if (track.itsClsSizeInLayer(1) != 0 && track.itsNCls() >= nITSClst) {
        
        if (track.tpcNClsCrossedRows() < 130 && track.tpcNClsCrossedRows() >= nTPCClst2) {
          histos.fill(HIST("PtHistogramSECOND"), track.pt());
          histos.fill(HIST("DCASECOND"), track.dcaXY());
          histos.fill(HIST("ChiSECOND"), track.tpcChi2NCl());
        }
      
      }

      /*
      // ITS Chi2, DCA, PT
      if (track.itsClsSizeInLayer(0) == 0 && track.itsNCls() >= nITSClst) {

        histos.fill(HIST("NoHitITS"), track.itsChi2NCl());
        histos.fill(HIST("NoHitITSDCA"), track.dcaXY());
        histos.fill(HIST("NoHitITSPT"), track.pt());
        histos.fill(HIST("NoHitDCAvsPt"), track.dcaXY(), track.pt());

      }

      if (track.itsClsSizeInLayer(0) != 0 && track.itsNCls() >= nITSClst && track.tpcNClsCrossedRows() > nTPCClst2) {

        histos.fill(HIST("H1370HitITS"), track.itsChi2NCl());
        histos.fill(HIST("H1370HitITSDCA"), track.dcaXY());
        histos.fill(HIST("H1370HitITSPT"), track.pt());
        histos.fill(HIST("DCAvsPt"), track.dcaXY(), track.pt());

      }

      if (track.itsNClsInnerBarrel() != 0) { //this means of course that the first one could also be empty I'm sure because the DCA worsened. 

        histos.fill(HIST("13HitITS"), track.itsChi2NCl());
        histos.fill(HIST("13HitITSDCA"), track.dcaXY());
        histos.fill(HIST("13HitITSPT"), track.pt());
      }

      if (track.itsClsSizeInLayer(0) != 0) {

        histos.fill(HIST("1AllHitITS"), track.itsChi2NCl());
        histos.fill(HIST("1AllHitITSDCA"), track.dcaXY());
        histos.fill(HIST("1AllHitITSPT"), track.pt());
      }



      if (mySelection.IsSelected(track)) {
       //if (track.isQualityTrack()){ 
        histos.fill(HIST("PtHistogram"), track.pt());
        histos.fill(HIST("SelectionHistogram"), track.tpcChi2NCl());
        histos.fill(HIST("DCAHistogram"), track.dcaXY());
        histos.fill(HIST("TrackSelectDCAvsPt"), track.dcaXY(), track.pt());

        //histos.fill(HIST("ITSChiHistogram"), track.itsChi2NCl());
        //histos.fill(HIST("ITSRefitChi"), track.itsChi2NCl());
        histos.fill(HIST("ITSChi"), track.itsChi2NCl());

        if (track.passedITSRefit()) {
        }
        */
      
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MCTrackSel>(cfgc)};
}