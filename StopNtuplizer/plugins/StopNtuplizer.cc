// -*- C++ -*-
//
// Package:    Stops/StopNtuplizer
// Class:      StopNtuplizer
// 
/**\class StopNtuplizer StopNtuplizer.cc Stops/StopNtuplizer/plugins/StopNtuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Cameron Bravo
//         Created:  Fri, 03 Oct 2014 22:46:38 GMT
//
//


#include <iostream>

#include "StopNtuplizer.h"

using namespace edm;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
StopNtuplizer::StopNtuplizer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

   

}


StopNtuplizer::~StopNtuplizer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
StopNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
    nEv_h->Fill(0.0);
    edm::Handle<edm::View<pat::Muon> > muons;
    iEvent.getByLabel("slimmedMuons",muons);
    edm::View<pat::Muon> patMuons = *muons;

    edm::Handle<edm::View<pat::Electron> > electrons;
    iEvent.getByLabel("slimmedElectrons",electrons);
    edm::View<pat::Electron> patElectrons = *electrons;

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel("offlineSlimmedPrimaryVertices",vertices);
    
    reco::VertexCollection vertexV = *(vertices.product());

    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByLabel("slimmedJets",jets);
    edm::View<pat::Jet> patJets = *jets;

    edm::Handle<pat::METCollection> metcoll;
    iEvent.getByLabel("slimmedMETs", metcoll);
    if(!metcoll.isValid()) return;

    const pat::MET *met=NULL;
    met=&(metcoll->front());

    muons_pt.clear();
    muons_phi.clear();
    muons_eta.clear();
    muons_E.clear();
    muons_ID.clear();
    bool tightMu = false;
    Nmuons = 0;

    elecs_pt.clear();
    elecs_phi.clear();
    elecs_eta.clear();
    elecs_E.clear();
    bool tightE = false;
    Nelecs = 0;

    jets_pt.clear();
    jets_phi.clear();
    jets_eta.clear();
    jets_E.clear();
    Njets = 0;

    tight_elec = -1;
    tight_muon = -1;

    NPV = vertexV.size();
    
    if(muons.isValid())
    {
        for(edm::View<pat::Muon>::const_iterator i_muon = patMuons.begin();i_muon != patMuons.end(); ++i_muon) 
        {
            if(!(i_muon->isLooseMuon())) continue;
            muons_pt.push_back(i_muon->pt());
            muons_eta.push_back(i_muon->eta());
            muons_phi.push_back(i_muon->phi());
            muons_E.push_back(i_muon->energy());
            Nmuons++;
            if(i_muon->isTightMuon(vertexV.at(0)))
            {
                if(!tightMu)
                {
                    tightMu = true;
                    tight_muon = Nmuons - 1;
                }
                muons_ID.push_back(1);
            }
            else
            {
                muons_ID.push_back(0);
            }
        }
    }

    if(electrons.isValid())
    {
        for(edm::View<pat::Electron>::const_iterator i_electron = patElectrons.begin();i_electron != patElectrons.end(); ++i_electron) 
        {
            TE_h->Fill(i_electron->electronID("eidTight"));
            LE_h->Fill(i_electron->electronID("eidLoose"));
            if(!tightE)
            {
                if(i_electron->electronID("eidTight") == 7)
                {
                    elecs_pt.push_back(i_electron->pt());
                    elecs_eta.push_back(i_electron->eta());
                    elecs_phi.push_back(i_electron->phi());
                    elecs_E.push_back(i_electron->energy());
                    tightE = true;
                    tight_elec = Nelecs;
                    Nelecs++;
                }
                else if(i_electron->electronID("eidLoose") == 7)
                {
                    elecs_pt.push_back(i_electron->pt());
                    elecs_eta.push_back(i_electron->eta());
                    elecs_phi.push_back(i_electron->phi());
                    elecs_E.push_back(i_electron->energy());
                    Nelecs++;
                }
            }
            else
            {
                if(i_electron->electronID("eidLoose") == 7)
                {
                    elecs_pt.push_back(i_electron->pt());
                    elecs_eta.push_back(i_electron->eta());
                    elecs_phi.push_back(i_electron->phi());
                    elecs_E.push_back(i_electron->energy());
                    Nelecs++;
                }
            }
        }
    }
    
    if(jets.isValid())
    {
        for(edm::View<pat::Jet>::const_iterator i_jet = patJets.begin();i_jet != patJets.end(); ++i_jet) 
        {
                jets_pt.push_back(i_jet->pt());
                jets_eta.push_back(i_jet->eta());
                jets_phi.push_back(i_jet->phi());
                jets_E.push_back(i_jet->energy());
                Njets++;
            
        }
    }

    met_pt = met->pt();
    met_eta = met->eta();
    met_phi = met->phi();
    met_E = met->energy();

    if(tightMu) tree->Fill();


}


// ------------ method called once each job just before starting event loop  ------------
void 
StopNtuplizer::beginJob()
{
    //gROOT->ProcessLine("#include <vector>");
    tree = fs->make<TTree>("tree","tree");
    nEv_h = fs->make<TH1D>("nEv_h","nEv_h",20,-2.0,5.0);
    TE_h = fs->make<TH1I>("TE_h","TE_h",21,-1,20);
    LE_h = fs->make<TH1I>("LE_h","LE_h",21,-1,20);
    tree->Branch("muons_pt",&muons_pt);
    tree->Branch("muons_eta",&muons_eta);
    tree->Branch("muons_phi",&muons_phi);
    tree->Branch("muons_E",&muons_E);
    tree->Branch("muons_ID",&muons_ID);
    tree->Branch("tight_muon",&tight_muon);
    tree->Branch("Nmuons",&Nmuons);
    tree->Branch("elecs_pt",&elecs_pt);
    tree->Branch("elecs_eta",&elecs_eta);
    tree->Branch("elecs_phi",&elecs_phi);
    tree->Branch("elecs_E",&elecs_E);
    tree->Branch("tight_elec",&tight_elec);
    tree->Branch("Nelecs",&Nelecs);
    tree->Branch("jets_pt",&jets_pt);
    tree->Branch("jets_eta",&jets_eta);
    tree->Branch("jets_phi",&jets_phi);
    tree->Branch("jets_E",&jets_E);
    tree->Branch("Njets",&Njets);
    tree->Branch("met_pt",&met_pt);
    tree->Branch("met_eta",&met_eta);
    tree->Branch("met_phi",&met_phi);
    tree->Branch("met_E",&met_E);
    tree->Branch("NPV",&NPV);
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
StopNtuplizer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
StopNtuplizer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
StopNtuplizer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
StopNtuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
StopNtuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
StopNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(StopNtuplizer);
