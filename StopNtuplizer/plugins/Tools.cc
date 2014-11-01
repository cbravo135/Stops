
 #include "Stops/StopNtuplizer/interface/Tools.h"
 
 



  //Electron Selector 

std::vector<const pat::Electron* > tools::ElectronSelector(const std::vector<pat::Electron>  & thePatElectrons, 
double v_electron_pt, 
double v_electron_eta, 
reco::Vertex::Point PV,
double v_electron_d0, 
double v_electron_reliso,
bool usePFiso, 
bool bool_electron_chargeConsistency,
edm::Handle< std::vector<reco::Conversion> > &theConversions,
reco::BeamSpot::Point BS)
{


v_electron_d0= 0.02;
double v_electron_dz = 0.2;
//bool bool_electron_ecalDriven= true;

  std::vector<const pat::Electron* > vElectrons;  
        for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) 
	{
	  std::cout<<"Electron "<<el->pt()<<" "<<el->eta()<<" "<<el->phi()<<std::endl;
          const reco::GsfTrackRef gsfTrack = el->gsfTrack();
         
	  if( el->pt() < v_electron_pt ) continue;
	  std::cout<<"electron eta"<<std::endl;
	  if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
	  // should we use SC eta ?  ele->superCluster()->eta();
	  std::cout<<"gap"<<std::endl;
	  if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
	  //	  std::cout<<"ecaldriven"<<std::endl;
	  //	  if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
	  std::cout<<"d0"<<std::endl;
	  if( TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  continue; 
	  std::cout<<"dz "<< el->gsfTrack()->dz(PV)<<" PV "<<PV.X()<<" "<<PV.Y()<<" "<<PV.Z() <<std::endl;
          if( TMath::Abs(el->gsfTrack()->dz(PV))  > v_electron_dz  ) continue;  

	  //	  if(!el->isGsfCtfScPixChargeConsistent() )  continue;
	       
	  // customized selection from RA5
	  /*	  if( el->pt() < 20. )
	    {
	     if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	    }
	  */
          
          //if( TMath::Fabs(( 1./ el->trackMomentumAtVtx().p() ) - ( 1./ el->ecalEnergy() ) ) > 0.05 ) continue;
	  std::cout<<"E/p"<<std::endl;
           if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
	   std::cout<<"conversion"<<std::endl;
	   //            bool vtxFitConversion = ConversionTools::hasMatchedConversion(*el, theConversions, BS);
            bool vtxFitConversion = ConversionTools::hasMatchedConversion( reco::GsfElectron (*el) , theConversions, BS);
           if( vtxFitConversion )  continue; 
	   std::cout<<"missing hits"<<std::endl;
	    if( el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > 1 ) continue; 
	  
	  if( TMath::Abs( el->eta() ) < 1.5 )
	    {
	      std::cout<<"deltaphi"<<std::endl;
	      if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.15 ) continue;
	      std::cout<<"deltaeta"<<std::endl;
	      if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
	      std::cout<<"sigietaieta"<<std::endl;
	      if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
	      std::cout<<"hovere"<<std::endl;
	      if( TMath::Abs(el->hadronicOverEm()) > 0.12 ) continue;  //recommended is 0.12 but HLT applies 0.1
	    }
	  else 
	    {
	      std::cout<<"deltaphi"<<std::endl;
	      if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.10 ) continue;
	      std::cout<<"deltaeta"<<std::endl;
	      if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.009 ) continue;
	      std::cout<<"sigietaieta"<<std::endl;
	      if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
	      std::cout<<"hovere"<<std::endl;
	      if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) continue;   /// at the HLT 0.075  recommended is 0.1
	    }

	  //float ecalIso = TMath::Abs(el->eta()) > 1.47 ? el->dr03EcalRecHitSumEt() : TMath::Max(el->dr03EcalRecHitSumEt()-1.,0.); 
	  //float detRelIso = (el->dr03TkSumPt() + el->dr03HcalTowerSumEt() + ecalIso ) / el->pt() ;
           //puChargedHadronIso(), particleIso()
          //float pfRelIso  = (el->chargedHadronIso() + el->neutralHadronIso() + el->photonIso() ) /el->pt() ;
          //float RelIso = detRelIso; 
	  //if(usePFiso) RelIso = pfRelIso;  
	  //	  if( RelIso  > v_electron_reliso ) continue;
	 
	  std::cout<<"Electron pass"<<std::endl;
	  vElectrons.push_back(&*el );
	}
      return vElectrons;
}
