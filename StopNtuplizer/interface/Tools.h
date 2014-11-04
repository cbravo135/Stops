#ifndef Tools_H
#define Tools_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h" 
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h" 
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"



#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"

//Root Classes

#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLegend.h"

//Standard C++ classes
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <ostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <memory>
#include <iomanip>

using namespace std;

namespace tools {

  struct event_info
  {
    std::vector<int> el,mu;
 
    float met;
    float ht;
    int channel;
    float weight;
    ULong64_t event;
    int  lumi;
    uint run;
    TString comment;
    event_info():  met(-1.), ht(-1.), channel(0), weight(1.), event(0), lumi(0), 
		   run(0), comment("") {}
    
  } ;

 
  float mass(float pt1 , float pt2, float eta1 , float eta2, float phi1, float phi2);
  
   void ERR( edm::InputTag& IT );


   struct lepton_pair 
   {
     int i; //first lepton index                                                                                      
     int j; //second lepton index                                                                                     
     float sumPt;
     int channel ;  //mumu =1 , emu =2, ee =3                                                                         
    };

/////////////////////////  Selectors //////////////////
std::vector<const pat::Muon* > MuonSelector(const std::vector<pat::Muon>  & thePatMuons, 
double v_muon_pt,	 
double v_muon_eta,  
reco::Vertex::Point PV,
double v_muon_d0,	 
double v_muon_reliso,
bool usePFiso);


 std::vector<const pat::Muon* > MuonEWKinoSelector(const std::vector<pat::Muon>  & thePatMuons,
					     double v_muon_pt,
					     double v_muon_eta,
                        		     reco::Vertex::Point PV);

 std::vector<const pat::Muon* > MuonEWKinoAntiSelector(const std::vector<pat::Muon>  & thePatMuons,
						   double v_muon_pt,
						   double v_muon_eta,
						   reco::Vertex::Point PV);

 std::vector<const pat::Muon* > MuonEWKinoIsoSelector(const std::vector<pat::Muon>  & thePatMuons,
						       double v_muon_pt,
						       double v_muon_eta,
						       reco::Vertex::Point PV);

 std::vector<const pat::Muon* > MuonSSSelector(const std::vector<pat::Muon>  & thePatMuons,
						       double v_muon_pt,
						       double v_muon_eta,
						       reco::Vertex::Point PV);

 std::vector<const pat::Muon* > MuonAntiSSSelector(const std::vector<pat::Muon>  & thePatMuons,
					       double v_muon_pt,
					       double v_muon_eta,
					       reco::Vertex::Point PV);

 std::vector<const pat::Muon* > MuonNoIsoSSSelector(const std::vector<pat::Muon>  & thePatMuons,
					       double v_muon_pt,
					       double v_muon_eta,
					       reco::Vertex::Point PV);
 std::vector<const pat::Muon* > MuonIsoSSSelector(const std::vector<pat::Muon>  & thePatMuons,
						    double v_muon_pt,
						    double v_muon_eta,
						    reco::Vertex::Point PV);


 std::vector<const pat::Muon* > MuonAntiSelector(const std::vector<pat::Muon>  & thePatMuons,
						      double v_muon_pt,
						      double v_muon_eta,
						      reco::Vertex::Point PV);




std::vector<const pat::Electron* > ElectronSelector(const std::vector<pat::Electron>  & thePatElectrons, 
double v_electron_pt, 
reco::Vertex::Point PV,
edm::Handle< std::vector<reco::Conversion> > &theConversions,
reco::BeamSpot::Point BS);

std::vector<const pat::Electron* > ElectronEWKinoSelector(const std::vector<pat::Electron>  & thePatElectrons,
						     double v_electron_pt,
						     double v_electron_eta,
						     double myRho,
						     reco::Vertex::Point PV,
	                                             edm::Handle< std::vector<reco::Conversion> > &theConversions,
						     reco::BeamSpot::Point BS);


 std::vector<const pat::Electron* > ElectronEWKinoConvSelector(const std::vector<pat::Electron>  & thePatElectrons,
							   double v_electron_pt,
							   double v_electron_eta,
							   double myRho,
							   reco::Vertex::Point PV,
							   edm::Handle< std::vector<reco::Conversion> > &theConversions,
							   reco::BeamSpot::Point BS);


 std::vector<const pat::Electron* > ElectronEWKinoAntiSelector(const std::vector<pat::Electron>  & thePatElectrons,
							   double v_electron_pt,
							   double v_electron_eta,
							   double myRho,
							   reco::Vertex::Point PV,
							   edm::Handle< std::vector<reco::Conversion> > &theConversions,
							   reco::BeamSpot::Point BS);

 std::vector<const pat::Electron* > ElectronEWKinoIsoSelector(const std::vector<pat::Electron>  & thePatElectrons,
							       double v_electron_pt,
							       double v_electron_eta,
							       double myRho,
							       reco::Vertex::Point PV,
							       edm::Handle< std::vector<reco::Conversion> > &theConversions,
							       reco::BeamSpot::Point BS);


 std::vector<const pat::Electron* > ElectronAntiSelector(const std::vector<pat::Electron>  & thePatElectrons,
							      double v_electron_pt,
							      double v_electron_eta,
							      double myRho,
							      reco::Vertex::Point PV,
							      edm::Handle< std::vector<reco::Conversion> > &theConversions,
							      reco::BeamSpot::Point BS);


 std::vector<const pat::Electron* > ElectronSSSelector(const std::vector<pat::Electron>  & thePatElectrons,
							       double v_electron_pt,
							       double v_electron_eta,
							       double myRho,
							       reco::Vertex::Point PV,
							       edm::Handle< std::vector<reco::Conversion> > &theConversions,
							       reco::BeamSpot::Point BS);

 std::vector<const pat::Electron* > ElectronAntiSSSelector(const std::vector<pat::Electron>  & thePatElectrons,
						       double v_electron_pt,
						       double v_electron_eta,
						       double myRho,
						       reco::Vertex::Point PV,
						       edm::Handle< std::vector<reco::Conversion> > &theConversions,
						       reco::BeamSpot::Point BS);

 std::vector<const pat::Electron* > ElectronNoIsoSSSelector(const std::vector<pat::Electron>  & thePatElectrons,
						       double v_electron_pt,
						       double v_electron_eta,
						       double myRho,
						       reco::Vertex::Point PV,
						       edm::Handle< std::vector<reco::Conversion> > &theConversions,
						       reco::BeamSpot::Point BS);

 std::vector<const pat::Electron* > ElectronIsoSSSelector(const std::vector<pat::Electron>  & thePatElectrons,
							    double v_electron_pt,
							    double v_electron_eta,
							    double myRho,
							    reco::Vertex::Point PV,
							    edm::Handle< std::vector<reco::Conversion> > &theConversions,
							    reco::BeamSpot::Point BS);




 std::vector<const reco::PFTau * > TauSelector(edm::Handle<reco::PFTauCollection> & PFTaus,
							   double v_tau_pt,
							   double v_tau_eta,
							   edm::Handle<reco::PFTauDiscriminator> & electron,
							   edm::Handle<reco::PFTauDiscriminator> & muon,
							   edm::Handle<reco::PFTauDiscriminator> & iso,
							   edm::Handle<reco::PFTauDiscriminator> & decay);

 std::vector<const reco::PFTau * > TauAntiSelector(edm::Handle<reco::PFTauCollection> & PFTaus,
					       double v_tau_pt,
					       double v_tau_eta,
					       edm::Handle<reco::PFTauDiscriminator> & electron,
					       edm::Handle<reco::PFTauDiscriminator> & muon,
					       edm::Handle<reco::PFTauDiscriminator> & iso,
					       edm::Handle<reco::PFTauDiscriminator> & decay);


std::vector<const pat::Jet* > JetSelector(const std::vector<pat::Jet>  & thePatJets, 
double  value_jet_et,  
double  value_jet_eta, 
bool    bool_jet_id, 
bool    jetLeptonCleaning, 
double  value_jet_leptonVetoDR,
std::vector<const pat::Electron*> & vElectrons,
std::vector<const pat::Muon*> & vMuons,
std::vector<const reco::PFTau * > & vTaus);
    

 double TriLeptonMass( std::vector<const pat::Electron*>  & thePatElectrons, std::vector<const pat::Muon*>  & thePatMuons,std::vector<const reco::PFTau * > & vTaus);
 double ZMass( std::vector<const pat::Electron*>  & thePatElectrons, std::vector<const pat::Muon*>  & thePatMuons,std::vector<const reco::PFTau * > & vTaus, int OSSF);
 double MTcalc( std::vector<const pat::Electron*>  & thePatElectrons, std::vector<const pat::Muon*>  & thePatMuons, double pfMET_px, double pfMET_py, int OSSF);
 double MT(TLorentzVector lep3,double pfMET_px, double pfMET_py, int OSSF);
 int DY( std::vector<const pat::Electron*>  & thePatElectrons, std::vector<const pat::Muon*>  & thePatMuons);
 int OSSF( std::vector<const pat::Electron*>  & thePatElectrons, std::vector<const pat::Muon*>  & thePatMuons,std::vector<const reco::PFTau * > & vTaus);
 bool PassEventAcceptance(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons);
 bool PassEventAcceptance2(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons);
 double MCTPerp1(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons, double pfMET_px, double pfMET_py, int OSSF);
 double MCTPerp2(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons, double pfMET_px, double pfMET_py, int OSSF);
 std::vector<TLorentzVector> LeptonSelector(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons, std::vector<const reco::PFTau * > & vTaus,int OSSF);
 double JER (double eta) ;
    
double getZMassForDiEle(const std::vector<pat::Electron>  & thePatElectrons,                                                                                                                                 
    edm::Handle< std::vector<reco::Conversion> > &theConversions,                                                                                                                                                         
    reco::BeamSpot::Point BS,                                                                                                                                                                    
    reco::Vertex::Point PV,  
    double  Rho,       
    const pat::Electron *myEle );  
    
    double getZMassForDiMu(const std::vector<pat::Muon>  & thePatMuons,   
                                reco::Vertex::Point PV,
                                const pat::Muon *myMu );

float getPtUncertainty(float pt, float eta);

}  


#endif
