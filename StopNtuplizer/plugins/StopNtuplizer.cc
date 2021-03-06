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

    parameters = iConfig;

    convCollectionLabel_ = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("convLabel"));


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

    edm::Handle< reco::BeamSpot > theBeamSpot;
    iEvent.getByLabel("offlineBeamSpot", theBeamSpot );
    reco::BeamSpot::Point  BS= theBeamSpot->position();

    edm::Handle< reco::GenParticleCollection > theGenParts;
    iEvent.getByLabel("prunedGenParticles", theGenParts );

    edm::Handle< reco::ConversionCollection > theConversions;
    iEvent.getByToken(convCollectionLabel_, theConversions);

    const pat::MET *met=NULL;
    met=&(metcoll->front());

    muons_pt.clear();
    muons_phi.clear();
    muons_eta.clear();
    muons_E.clear();
    muons_ID.clear();
    m_muons.clear();
    bool tightMu = false;
    Nmuons = 0;

    Ngen_muons = 0;
    gen_muons_pt.clear();
    gen_muons_eta.clear();
    gen_muons_phi.clear();
    gen_muons_E.clear();

    Ngen_elecs = 0;
    gen_elecs_pt.clear();
    gen_elecs_eta.clear();
    gen_elecs_phi.clear();
    gen_elecs_E.clear();


    elecs_pt.clear();
    elecs_phi.clear();
    elecs_eta.clear();
    elecs_E.clear();
    m_elecs.clear();
    //bool tightE = false;
    Nelecs = 0;

    jets_pt.clear();
    jets_phi.clear();
    jets_eta.clear();
    jets_E.clear();
    Njets = 0;

    tight_elec = -1;
    tight_muon = -1;

    NPV = vertexV.size();
    reco::Vertex::Point PV1p = vertexV.begin()->position();
    
    if(theGenParts.isValid())
    {
        for(unsigned int i = 0; i < theGenParts->size() ; i++) 
        {
            const reco::GenParticle &genPart = (*theGenParts)[i];
            /*std::cout << "pdgID: " << genPart.pdgId() << "\tstatus: " << genPart.status() << "\tmotherN: " << genPart.numberOfMothers();
            for(unsigned int ii = 0; ii < genPart.numberOfMothers(); ii++)
            {
                std::cout << "\tmotherID: " << genPart.motherRef(ii)->pdgId();
            }
            std::cout << std::endl;*/
            if(genPart.numberOfMothers() == 0) continue;
            //if(!((TMath::Abs(genPart.pdgId()) == 11 || TMath::Abs(genPart.pdgId()) == 13) && genPart.status() == 1 && TMath::Abs(genPart.motherRef(0)->pdgId()) == 24)) continue;
            if(!((TMath::Abs(genPart.pdgId()) == 11 || TMath::Abs(genPart.pdgId()) == 13) && TMath::Abs(genPart.motherRef(0)->pdgId()) == 24)) continue;
            if(TMath::Abs(genPart.pdgId()) == 11)
            {
                gen_elecs_pt.push_back(genPart.pt());
                gen_elecs_eta.push_back(genPart.eta());
                gen_elecs_phi.push_back(genPart.phi());
                gen_elecs_E.push_back(genPart.energy());
                Ngen_elecs++;
            }
            else
            {
                gen_muons_pt.push_back(genPart.pt());
                gen_muons_eta.push_back(genPart.eta());
                gen_muons_phi.push_back(genPart.phi());
                gen_muons_E.push_back(genPart.energy());
                Ngen_muons++;
            }
        }
    }

    if(muons.isValid())
    {
        for(edm::View<pat::Muon>::const_iterator i_muon = patMuons.begin();i_muon != patMuons.end(); ++i_muon) 
        {
            if(!(i_muon->isLooseMuon())) continue;
            if(i_muon->pt() < 30) continue;
            double chargedHadronIso = i_muon->pfIsolationR04().sumChargedHadronPt;
            double neutralHadronIso = i_muon->pfIsolationR04().sumNeutralHadronEt;
            double photonIso = i_muon->pfIsolationR04().sumPhotonEt;
            double beta = i_muon->pfIsolationR04().sumPUPt;
            double pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/i_muon->pt() ;
            if(pfRelIso>0.15) continue;

            muons_pt.push_back(i_muon->pt());
            muons_eta.push_back(i_muon->eta());
            muons_phi.push_back(i_muon->phi());
            muons_E.push_back(i_muon->energy());
            Nmuons++;
            TLorentzVector muon4V;
            muon4V.SetPtEtaPhiE(i_muon->pt(),i_muon->eta(),i_muon->phi(),i_muon->energy());
            bool matched = false;
            for(int i =0; i < Ngen_muons; i++)
            {
                TLorentzVector genmV;
                genmV.SetPtEtaPhiE(gen_muons_pt[i],gen_muons_eta[i],gen_muons_phi[i],gen_muons_E[i]);
                if(genmV.DeltaR(muon4V) < 0.05 && !matched)
                {
                    matched = true;
                    m_muons.push_back(i);
                }
            }
            if(!matched) m_muons.push_back(-1);
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

    std::vector<pat::Electron> ElecV;

    if(electrons.isValid())
    {
        for(edm::View<pat::Electron>::const_iterator i_E = patElectrons.begin();i_E != patElectrons.end(); ++i_E) 
        {
            ElecV.push_back(*i_E);
        }
        std::vector<const pat::Electron* > Es = tools::ElectronSelector(ElecV,10.0,PV1p,theConversions,BS);
        int Ne = Es.size();
        for(int i = 0;i < Ne; ++i)
        {
            if(Es[i]->pt() < 30) continue;
            elecs_pt.push_back(Es[i]->pt());
            elecs_eta.push_back(Es[i]->eta());
            elecs_phi.push_back(Es[i]->phi());
            elecs_E.push_back(Es[i]->energy());
            Nelecs++;
            TLorentzVector elec4V;
            elec4V.SetPtEtaPhiE(Es[i]->pt(),Es[i]->eta(),Es[i]->phi(),Es[i]->energy());
            bool matched = false;
            for(int i =0; i < Ngen_elecs; i++)
            {
                TLorentzVector geneV;
                geneV.SetPtEtaPhiE(gen_elecs_pt[i],gen_elecs_eta[i],gen_elecs_phi[i],gen_elecs_E[i]);
                if(geneV.DeltaR(elec4V) < 0.05 && !matched)
                {
                    matched = true;
                    m_elecs.push_back(i);
                }
            }
            if(!matched) m_elecs.push_back(-1);

        }

    }
    
    if(jets.isValid())
    {
        for(edm::View<pat::Jet>::const_iterator i_jet = patJets.begin();i_jet != patJets.end(); ++i_jet) 
        {
            if(i_jet->pt() < 30) continue;
            bool keepJet = true;
            TLorentzVector jet4v;
            jet4v.SetPtEtaPhiE(i_jet->pt(),i_jet->eta(),i_jet->phi(),i_jet->energy());
            for(int i = 0;i < Nelecs; ++i) 
            {
                TLorentzVector e4v;
                e4v.SetPtEtaPhiE(elecs_pt.at(i),elecs_eta.at(i),elecs_phi.at(i),elecs_E.at(i));
                if(e4v.DeltaR(jet4v) < 0.3) keepJet = false;
            }
            for(int i = 0;i < Nmuons; ++i) 
            {
                TLorentzVector m4v;
                m4v.SetPtEtaPhiE(muons_pt.at(i),muons_eta.at(i),muons_phi.at(i),muons_E.at(i));
                if(m4v.DeltaR(jet4v) < 0.3) keepJet = false;
            }
            if(!keepJet) continue;
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

    //if(tightMu || Nelecs > 0) tree->Fill();
    tree->Fill();


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
    tree->Branch("m_muons",&m_muons);
    tree->Branch("elecs_pt",&elecs_pt);
    tree->Branch("elecs_eta",&elecs_eta);
    tree->Branch("elecs_phi",&elecs_phi);
    tree->Branch("elecs_E",&elecs_E);
    tree->Branch("tight_elec",&tight_elec);
    tree->Branch("Nelecs",&Nelecs);
    tree->Branch("m_elecs",&m_elecs);
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
    tree->Branch("Ngen_muons",&Ngen_muons);
    tree->Branch("gen_muons_pt",&gen_muons_pt);
    tree->Branch("gen_muons_eta",&gen_muons_eta);
    tree->Branch("gen_muons_phi",&gen_muons_phi);
    tree->Branch("gen_muons_E",&gen_muons_E);
    tree->Branch("Ngen_elecs",&Ngen_elecs);
    tree->Branch("gen_elecs_pt",&gen_elecs_pt);
    tree->Branch("gen_elecs_eta",&gen_elecs_eta);
    tree->Branch("gen_elecs_phi",&gen_elecs_phi);
    tree->Branch("gen_elecs_E",&gen_elecs_E);
  
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
