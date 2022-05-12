//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 23 13:02:03 2021 by ROOT version 6.22/06
// from TTree EligantSelector/Energy Station
// found on file: run1005_03.root
//////////////////////////////////////////////////////////

#pragma once

#ifndef EligantSelector_h
#define EligantSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH2.h>
#include <TString.h>
#include <TObjString.h>

#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>
#include <fstream>
#include <limits>
#include <csignal>
#include <ctime>
#include <TEntryList.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TVector3.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <TCutG.h>
#include <TList.h>
#include <TKey.h>
#include <TSysEvtHandler.h>
#include <TSystem.h>
#include <TApplication.h>
#include <deque>
#include "TTreeIndex.h"
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* getenv */
#include <map>
#include <functional>
#include <queue>
#include <vector>

//#include "TObjString.h."
// Headers needed by this particular selector


class EligantSelector : public TSelector {
public :
    TTreeReader     fReader;  //!the tree reader
    TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
    TFile          *foutFile;
    TTree          *outputTree;


//   UChar_t	 uMod; 
//   UChar_t	 uChannel; 

 class TDelilaEvent {
 public:
    UChar_t	        fMod; 
    UChar_t	        fChannel;    
//     ULong64_t	    fTimeStamp;
//     double_t	    fTimeStampFS;//FineTS
    UShort_t	    fEnergy;//ChargeLong
    UShort_t	    fEnergyShort;//ChargeShot  
    UShort_t        det_def;//0 - nothing; 1 - core/single HPge; 2 - segment; 3 - CeBr; 4 - CsI; 5 - BGO1; 6 - BGO2; 7 -BGO - 3; 8 - solar cell; 9 - pulser
    float	        EnergyCal;
    float           EnergyDC;
    UShort_t        domain;
    UShort_t        cs_domain;
    UShort_t        channel;//ch daq
    UShort_t        CS;//0 - no; 1 - yes
    double_t        Time;
    double_t        TimeTrg;
    double_t        TimeBunch;
    Float_t         theta;
    Float_t	        phi;
    ULong64_t       trg;
    UShort_t        bunch;
    UShort_t        fold;    
    TDelilaEvent(): domain(-1),channel(-1),fEnergy(-1),CS(0),cs_domain(0),Time(0),trg(0),bunch(0),fold(0){};
 };
 
    ULong64_t	    fTimeStamp;
    double_t	    fTimeStampFS;//FineTS
    UShort_t	    fEnergyShort;//ChargeShot  

  class TDelilaDetector { 
  public:
    Int_t	 dom;
    int	 ch;//ch daq
//     Int_t	 serial;
    TString	 serial;
    Float_t  theta;
    Float_t	 phi;  
    UShort_t detType;//0 - nothing; 1 - core; 2 - segment; 3 - CeBr; 4 - CsI; 5 - BGO1; 6 - BGO2; 9 - pulser
    Int_t	 TimeOffset;
    double_t bgo_time_offset;
    Int_t 	 threshold; 
    Int_t	 pol_order;
    Int_t    cs_dom;
    std::vector<float> calibE;
    TDelilaDetector(): dom(-1),phi(-1),theta(-1),TimeOffset(0),calibE(0),threshold(-1),ch(-1),pol_order(-1){};
 };
 
  std::deque<TDelilaEvent> delilaQu;
  
//   std::map<unsigned int, TDelilaDetector > LUT_DELILA;
  std::map<int, TDelilaDetector > LUT_ELIGANT;    
  std::map<int, int > LUT_TA;
  std::map<int, double_t > LUT_TA_TRG;

  TDelilaEvent DelilaEvent;  
  TDelilaEvent DelilaEventTreated ;
  TDelilaEvent lastDelilaEvent;  
  TDelilaEvent lastEliadeZeroEvent;
  TDelilaEvent LastTriggerEvent;
  TDelilaEvent LastBunchEvent;    
  
  TDelilaEvent PulserEvent;
  TDelilaEvent DomainZeroEvent;    
  
  TDelilaEvent startEventCore;  
  TDelilaEvent startEventSegment;  

  TBranch *b_channel;
  TBranch *b_tstmp;
  TBranch *b_tstmp_fine;
  TBranch *b_energ;  //ChargeLong
  TBranch *b_energ_short;  //ChargeShot
  TBranch *b_mod;  
  
  Long64_t nb_entries;
  
  Long64_t bunch_length;
  Long64_t bunch_reset;
 
  TH1F* hChannelHit;
  TH1F* hDomainHit;
  TH1F* hSegmentHit;
  TH1F* hDetTypeHit;

  std::map<int, TH1F*> hDelila;
  std::map<int, TH1F*> hDelilaCS;
  std::map<int, TH1F*> hDelilaDC;
  std::map<int, TH1F*> hDelilaCS_DC;
  
  std::map<int, TH1F*> hDelila_long;
  std::map<int, TH1F*> hDelilaCS_long;
  std::map<int, TH1F*> hDelilaDC_long;
  std::map<int, TH1F*> hDelilaCS_DC_long;

  std::map<UInt_t, std::string> detector_name;

  TH1F* hTriggerTrigger;
  TH1F* hBunchBunch;
  TH1F* hNBunch;
  TH1F* hBunchFold;
  TH1F* hTimeInBunch;
  TH1F* hEventsPerTrigger;
  TH2F* mFoldEnergy;
  
  TH1F* hCoincID;
  TH1F* hTriggerDomain;
  
  TH1F* hNN_time_diff;
  
  TH1F* hdelilaQu_size;
  
  TH2F* mDelila;
  TH2F* mDelila_raw;
  TH2F* mDelilaCS;
  TH2F* mDelilaDC;//keV
  TH2F* mDelilaCS_DC;//keV
  
//   TH2F* mDelila_long;
//   TH2F* mDelila_raw_long;
//   TH2F* mDelilaCS_long;
//   TH2F* mDelilaDC_long;//keV
//   TH2F* mDelilaCS_DC_long;//keV

  TH2F* mGammaGammaDC;
  TH2F* mGammaGammaCS_DC;
  
  TH2F* mEnergyTimeDiff_trigger;
  TH2F* mDomainTimeDiff_trigger;
  TH2F* mDomainTimeDiff_bunch;
  
  TH3F* cGGTheta;
  TH3F* cGGG;
  
  std::map<int, TH2F*> mGG;
  std::map<int, TH2F*> mGG_theta;
  std::map<int, TH2F*> mGG_CS;
  std::map<int, TH2F*> mGG_DC;
  std::map<int, TH2F*> mGG_CS_DC;
  std::map<int, TH2F*> mGG_time_diff;
  
  std::map<int, TH2F*> mGG_long;
  std::map<int, TH2F*> mGG_CS_long;
  std::map<int, TH2F*> mGG_DC_long;
  std::map<int, TH2F*> mGG_CS_DC_long;
  
//   std::map<int, TH2F*> mTimeDiff;
  std::map<int, TH1F*> hMult;

  std::map<UInt_t, std::string> gg_coinc_id;
  std::map<UInt_t, Float_t> coinc_gates;//in ps

  
  std::map<int, std::string> domain_list;
  std::map<int, TH2F*> mEnergy_time_diff;
  
  TH2F* mTimeDiffCS;
  
  TH2F* mThetaPhi; 
  TH2F* mGammaGamma;
  TH2F* mTimeDiff_gg;
  TH1F* hMult_gg;
  TH2F* mGammaGammaCS;
  TH2F* mTimeDiff_gg_CS;
  TH1F* hMult_gg_CS;
  
  TH2F* mLaBr_LabBr_time_diff;

  
  TH2F *mPulser0TimeDiff;
  
  TH1F* hTimeDiffPulser;
  TH2F* mPulserPulser;
  
  TH1F* hTimeZero;
  TH1F* hTimeSort;
    
  TH2F* mTimeCalibDomain0;
  TH2F* mTimeCalib;
  TH2F* mTimeCalibBGO;
  TH2F* mTimeCalibBGO_cs_dom;
  
  TH2F* mNeutron;
  TH2F* mShortLong;
//   TH2F* mTimeCalibDomain0;
  
    
  std::clock_t start;
  double duration;
 
  ULong64_t nevents;
//   ULong64_t nevents_reset;
  int reset_counter;
 
  ULong64_t lastTime;
  
  TObjArray *toks;

   EligantSelector(TTree * /*tree*/ =0) { }
   virtual ~EligantSelector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   virtual Long64_t GetEntries() { return fChain ? fChain->GetEntries() : 0;}

   virtual void  Read_ELIADE_LookUpTable();
   virtual void  Read_TimeAlignment_LookUpTable();
   virtual void  Read_TimeAlignment_Trigger();
   virtual void  Read_Confs();
   virtual void  Print_ELIADE_LookUpTable();
   virtual void  Print_TimeAlignment_LookUpTable();
   virtual void  Print_TimeAlignment_Trigger_LookUpTable();

   virtual float CalibDet(float,int);
   virtual int GetCoincTimeCorrection(int dom1, int dom2);
      
   virtual void TreatDelilaEvent();
   virtual int GetCoincID(int dom1, int dom2);
   virtual int GetCoinc_det_def(int det_def1, int det_def2);
   virtual void CheckPulserAllignement(int zero_dom);
   virtual void PrintDelilaEvent(TDelilaEvent &ev_);
   virtual void SetUpNewTrigger();
   virtual void FillOutputTree();
   virtual bool TriggerDecision();
   virtual void TimeAlignementTrigger();
   
   virtual void TreatLaBrSingle();
   virtual void TreatHpGeSingle();
   virtual void TreatNeutronSingle();
   
   virtual void TreaNeutronMultiplicity();

   ClassDef(EligantSelector,0);
   
};

#endif

#ifdef EligantSelector_cxx
void EligantSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
   
   
  foutFile->cd();
  outputTree = new TTree("SelectedDelila","SelectedDelila");
//   outputTree->Branch("fTEventTS",&DelilaEventTreated .fTimeStamp,"TimeStamp/l");
//   outputTree->Branch("fTEventFS",&DelilaEventTreated .fTimeStampFS,"TimeStamp/D");
  outputTree->Branch("fTimeBunch",&DelilaEventTreated .TimeBunch,"TimeBunch/D");
  outputTree->Branch("fTime",&DelilaEventTreated .Time,"Time/D");
  outputTree->Branch("fEnergy",&DelilaEventTreated .fEnergy,"Energy/F");
  outputTree->Branch("fEnergy_kev",&DelilaEventTreated .EnergyCal,"Energy_kev/F");
  outputTree->Branch("fEnergyDC",&DelilaEventTreated .EnergyDC,"EnergyDC/F");
  outputTree->Branch("fDomain",&DelilaEventTreated .domain,"Domain/b");
  outputTree->Branch("fDetType",&DelilaEventTreated .det_def,"def/b");
  outputTree->Branch("fCS",&DelilaEventTreated .CS,"CS/b");
  outputTree->Branch("fTRG",&DelilaEventTreated .trg,"Trigger/b");
  outputTree->Branch("fFold",&DelilaEventTreated .fold,"Fold/b");
}

Bool_t EligantSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef EligantSelector_cxx
