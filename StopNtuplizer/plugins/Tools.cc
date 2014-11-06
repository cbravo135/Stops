
 #include "Stops/StopNtuplizer/interface/Tools.h"
 
 



  //Electron Selector 

std::vector<const pat::Electron* > tools::ElectronSelector(const std::vector<pat::Electron>  & thePatElectrons, 
double v_electron_pt, 
reco::Vertex::Point PV,
edm::Handle< std::vector<reco::Conversion> > &theConversions,
reco::BeamSpot::Point BS)
{



    std::vector<const pat::Electron* > vElectrons;  
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) 
	{
	    //std::cout<<"Electron "<<el->pt()<<" "<<el->eta()<<" "<<el->phi()<<std::endl;
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
	    if( el->pt() < v_electron_pt ) continue;
	    //std::cout<<"gap"<<std::endl;
	    if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
	    //std::cout<<"conversion"<<std::endl;
        bool vtxFitConversion = ConversionTools::hasMatchedConversion( reco::GsfElectron (*el) , theConversions, BS);
        if( vtxFitConversion )  continue; 
	    //std::cout<<"missing hits"<<std::endl;
	    if( el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > 1 ) continue; 

        reco::GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
        float absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
        double relIsoWithDBeta = absiso/el->pt();
	  
	    if( TMath::Abs( el->eta() ) <= 1.479 )
	    {
	        //std::cout<<"deltaphi"<<std::endl;
	        if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.0323 ) continue;
	        //std::cout<<"deltaeta"<<std::endl;
	        if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.0106 ) continue;
	        //std::cout<<"sigietaieta"<<std::endl;
	        if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.0107 ) continue;
	        //std::cout<<"hovere"<<std::endl;
	        if( TMath::Abs(el->hadronicOverEm()) > 0.067 ) continue;  //recommended is 0.12 but HLT applies 0.1
	        //std::cout<<"d0"<<std::endl;
	        if( TMath::Abs(el->gsfTrack()->dxy(PV)) > 0.0131  )  continue; 
	        //std::cout<<"dz "<< el->gsfTrack()->dz(PV)<<" PV "<<PV.X()<<" "<<PV.Y()<<" "<<PV.Z() <<std::endl;
            if( TMath::Abs(el->gsfTrack()->dz(PV))  > 0.22310  ) continue;  
	        //std::cout<<"E/p"<<std::endl;
            if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.1043 ) continue;
            if(relIsoWithDBeta > 0.2179) continue;
	    }
	    else 
	    {
	        //std::cout<<"deltaphi"<<std::endl;
	        if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.0455 ) continue;
	        //std::cout<<"deltaeta"<<std::endl;
	        if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.0108 ) continue;
	        //std::cout<<"sigietaieta"<<std::endl;
	        if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.0318 ) continue;
	        //std::cout<<"hovere"<<std::endl;
	        if( TMath::Abs(el->hadronicOverEm()) > 0.097 ) continue;   /// at the HLT 0.075  recommended is 0.1
	        //std::cout<<"d0"<<std::endl;
	        if( TMath::Abs(el->gsfTrack()->dxy(PV)) > 0.0845  )  continue; 
	        //std::cout<<"dz "<< el->gsfTrack()->dz(PV)<<" PV "<<PV.X()<<" "<<PV.Y()<<" "<<PV.Z() <<std::endl;
            if( TMath::Abs(el->gsfTrack()->dz(PV))  > 0.7523  ) continue;  
	        //std::cout<<"E/p"<<std::endl;
            if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.1201 ) continue;
            if(relIsoWithDBeta > 0.254) continue;
	    }

	 
	    //std::cout<<"Electron pass"<<std::endl;
	    vElectrons.push_back(&*el );
	}
    return vElectrons;
}
