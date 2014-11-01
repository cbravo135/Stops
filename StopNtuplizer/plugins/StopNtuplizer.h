//system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Stops/StopNtuplizer/interface/Tools.h"

#include "TTree.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TROOT.h"

#include "math.h"

#include <vector>


class StopNtuplizer : public edm::EDAnalyzer {
    public:
        explicit StopNtuplizer(const edm::ParameterSet&);
        ~StopNtuplizer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::Service<TFileService> fs;

        TTree * tree;
        TH1D * nEv_h;
        TH1I * TE_h;
        TH1I * LE_h;

        //Tree Variables

        std::vector<double> muons_pt;
        std::vector<double> muons_eta;
        std::vector<double> muons_phi;
        std::vector<double> muons_E;
        std::vector<int> muons_ID;
        int tight_muon;
        int Nmuons;

        std::vector<double> elecs_pt;
        std::vector<double> elecs_eta;
        std::vector<double> elecs_phi;
        std::vector<double> elecs_E;
        int tight_elec;
        int Nelecs;

        std::vector<double> jets_pt;
        std::vector<double> jets_eta;
        std::vector<double> jets_phi;
        std::vector<double> jets_E;
        int Njets;

        double met_pt;
        double met_eta;
        double met_phi;
        double met_E;

        int NPV;

};
