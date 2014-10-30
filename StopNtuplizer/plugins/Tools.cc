
 #include "Stops/StopNtuplizer/interface/Tools.h"
 
 
  float tools::mass(float pt1 , float pt2, float eta1 , float eta2, float phi1, float phi2)
  {
    TLorentzVector v1,v2;
    v1.SetPtEtaPhiM(pt1, eta1, phi1, 0.);
    v2.SetPtEtaPhiM(pt2,eta2, phi2, 0.);
    v1 = v1 + v2;  
    return v1.Mag();
  }

   void tools::ERR( edm::InputTag& IT )
  {
    cerr << "[ERROR] : " << IT << " is not a valid input label for this event.  SKIPPING EVENT " << endl;
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double tools::MTcalc(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons, double pfMET_px, double pfMET_py, int OSSF){
  TLorentzVector lep1, lep2, dilep, lep3;  // lepton 4-vectors                                                                                                                                                                             
  double difference=10000;
  double pfMET_pt=sqrt(pfMET_px*pfMET_px+pfMET_py*pfMET_py);
  if(OSSF<1) return -1;
  if(OSSF==1){
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
      const pat::Electron *ele1 = thePatElectrons[i];
      for(unsigned int j = i+1 ; j < thePatElectrons.size() ;j++ ) {
	const pat::Electron *ele2 = thePatElectrons[j];
	if(ele2->charge()==ele1->charge()) continue;
	lep1.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
	lep2.SetPtEtaPhiE(ele2->pt(),ele2->eta(),ele2->phi(),ele2->energy());
	dilep=lep1+lep2;
	if (fabs(dilep.M()-91)<difference) {
	  difference=fabs(dilep.M()-91);
	  if(thePatElectrons.size()==2){
	    const pat::Muon *mu1 = thePatMuons[0];
	    lep3.SetPtEtaPhiE(mu1->pt(),mu1->eta(),mu1->phi(),mu1->energy());
	  }
	  else{
	    int electron=3-i-j;
	    const pat::Electron *ele3 = thePatElectrons[electron];
	    lep3.SetPtEtaPhiE(ele3->pt(),ele3->eta(),ele3->phi(),ele3->energy());
	  }
	}

      }
    }
    for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
      const pat::Muon *mu1 = thePatMuons[i];
      for(unsigned int j = i+1 ; j < thePatMuons.size() ;j++ ) {
	const pat::Muon *mu2 = thePatMuons[j];
	if(mu2->charge()==mu1->charge()) continue;
	lep1.SetPtEtaPhiE(mu1->pt(),mu1->eta(),mu1->phi(),mu1->energy());
	lep2.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
	dilep=lep1+lep2;
	if (fabs(dilep.M()-91)<difference) {
	  difference=fabs(dilep.M()-91);
	  if(thePatMuons.size()==2){
	    const pat::Electron *ele1 = thePatElectrons[0];
	    lep3.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
	  }
	  else{
	    int muon=3-i-j;
	    const pat::Muon *mu3 = thePatMuons[muon];
	    lep3.SetPtEtaPhiE(mu3->pt(),mu3->eta(),mu3->phi(),mu3->energy());
	  }
	}
      }
    }
  }
  else if(OSSF==2){
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
      const pat::Electron *ele1 = thePatElectrons[i];
      for(unsigned int j = 0 ; j < thePatMuons.size() ;j++ ) {
	const pat::Muon *mu2 = thePatMuons[j];
	if(mu2->charge()==ele1->charge()) continue;
	lep1.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
	lep2.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
	dilep=lep1+lep2;
	if(fabs(dilep.M()-50)<difference){
	  difference=(fabs(dilep.M()-50));
	  if(thePatElectrons.size()>1){
	    int electron=1-i;
	    const pat::Electron *ele2 = thePatElectrons[electron];
	    lep3.SetPtEtaPhiE(ele2->pt(),ele2->eta(),ele2->phi(),ele2->energy());
	  }
	  if(thePatMuons.size()>1){
	    int muon=1-i;
	    const pat::Muon *mu2 = thePatMuons[muon];
	    lep3.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
	  }


	}
      }
    }
  }



  if(difference<10000) return sqrt((pfMET_pt+lep3.Pt())*(pfMET_pt+lep3.Pt())-(pfMET_px+lep3.Px())*(pfMET_px+lep3.Px())-(pfMET_py+lep3.Py())*(pfMET_py+lep3.Py()));
  else return -1;
}


  //Muon Selector 
std::vector<const pat::Muon* > tools::MuonSelector(const std::vector<pat::Muon>  & thePatMuons, 
double v_muon_pt,	 
double v_muon_eta,  
reco::Vertex::Point PV,
double v_muon_d0,	 
double v_muon_reliso,
bool usePFiso)
{  

v_muon_d0= 0.02;
double  v_muon_dz = 0.5;
int v_muon_numberOfMatchedStations =2; 
int v_muon_nPixelValidHits =1.;
int v_muon_numberOfValidMuonHits = 1;
double v_muon_chi2Norm= 10;
int v_muon_nValidHits= 6; 
double v_muon_hadVetoEt= 6.;
double v_muon_emVetoEt=4.; 

   std::vector<const pat::Muon* > vMuons;
      for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) 
	{
	   std::cout<<"Muon "<<mu->pt()<<" "<<mu->eta()<<" "<<mu->phi()<<std::endl;
	  if ( mu->pt()  < v_muon_pt ) continue;
	  std::cout<<"Eta cut"<<std::endl;
	  if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
//	  if ( !mu->isTrackerMuon() ) continue; 
	  std::cout<<"Global Muon"<<std::endl;
	  if ( !mu->isGlobalMuon()  ) continue;
	  std::cout<<"PF Muon"<<std::endl;
	  if ( !mu->isPFMuon() ) continue; 
	  std::cout<<"#stations"<<std::endl;
	  if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
	  const reco::TrackRef innerTrack = mu->innerTrack();
	  if( innerTrack.isNull() ) continue;
	  std::cout<<"valid_hits"<<std::endl;
	  if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue; 
	  std::cout<<" pixel hits"<<std::endl;
	  if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
	  
	  const reco::TrackRef globalTrack = mu->globalTrack() ;
	  if( globalTrack.isNull() ) continue;
	  std::cout<<"chi2"<<std::endl;
	  if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
	  std::cout<<"valid_muon_hits"<<std::endl;
	  if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
	  std::cout<<"d0"<<std::endl;
          if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;   
	  std::cout<<"dz"<<std::endl;
  	  if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;   

	  //	  if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
	  //	  if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;

          float pfRelIso  = (mu->pfIsolationR03().sumChargedHadronPt + mu->pfIsolationR03().sumNeutralHadronEt + mu->pfIsolationR03().sumPhotonEt) /mu->pt() ;
	  float detRelIso = (mu->isolationR03().emEt + mu->isolationR03().hadEt + mu->isolationR03().sumPt) / mu->pt() ;
          float RelIso = detRelIso; 
	  if(usePFiso) RelIso = pfRelIso;  
	  //	  if( RelIso > v_muon_reliso ) continue;
	  std::cout<<"Muon pass"<<std::endl;
	  vMuons.push_back(&*mu);
       }
       
       return vMuons;      

	  //float dxy = ( (mu->vy() - PV.y())*mu->px() - (mu->vx() - PV.x()) * mu->py()) / mu->pt();   /// impact Para w.r.t PV
	  //if(TMath::Abs(dxy )   > v_muon_d0  ) continue;
	  //TVector3 momPerp(0,0,0);
  	  //momPerp.SetPtEtaPhi(mu->pt(),mu->eta(),mu->phi());
  	  //TVector3 posPerp(mu->vx()-PV.x(), mu->vy() - PV.y(), 0);
 	  //float dzcorr = mu->vz() - PV.z() - posPerp.Dot(momPerp)/mu->pt() * (mu->pz()/mu->pt());
	  //if(TMath::Abs(dzcorr )   > v_muon_dz  ) continue;
}



std::vector<const pat::Muon* > tools::MuonEWKinoSelector(const std::vector<pat::Muon>  & thePatMuons,
						   double v_muon_pt,
						   double v_muon_eta,
						   reco::Vertex::Point PV)
{

     double v_muon_d0= 0.02;
  //  double v_muon_d0=10.2;
   double  v_muon_dz = 0.5;
  //double  v_muon_dz = 20.5;   

  int v_muon_numberOfMatchedStations =2;
  int v_muon_nPixelValidHits =1.;
  int v_muon_numberOfValidMuonHits = 1;
  double v_muon_chi2Norm= 10;
  int v_muon_nValidHits= 6;

  std::vector<const pat::Muon* > vMuons;
  for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
    {
       if ( mu->pt()  < v_muon_pt ) continue;

       if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
       if ( !mu->isGlobalMuon()  ) continue;
       if ( !mu->isPFMuon() ) continue;
      if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim                                                                                                                     
      const reco::TrackRef innerTrack = mu->innerTrack();
      if( innerTrack.isNull() ) continue;
      if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
      if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
      if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;

      const reco::TrackRef globalTrack = mu->globalTrack() ;
      if( globalTrack.isNull() ) continue;
      if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
      if(globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
      if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
      if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
                                                                                                                                                          

      //      double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
      //      double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
      //      double photonIso = mu->pfIsolationR03().sumPhotonEt;
      //      double beta = mu->pfIsolationR03().sumPUPt;
      //      double pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
      //      if(pfRelIso>0.15) continue;
      double chargedHadronIso = mu->pfIsolationR04().sumChargedHadronPt;
      double neutralHadronIso = mu->pfIsolationR04().sumNeutralHadronEt;
      double photonIso = mu->pfIsolationR04().sumPhotonEt;
      double beta = mu->pfIsolationR04().sumPUPt;
      double pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
      if(pfRelIso>0.15) continue;


      vMuons.push_back(&*mu);

    }

  return vMuons;
                                                                                                                      
}


std::vector<const pat::Muon* > tools::MuonSSSelector(const std::vector<pat::Muon>  & thePatMuons,
							 double v_muon_pt,
							 double v_muon_eta,
							 reco::Vertex::Point PV)
{

  double v_muon_d0= 0.005;
  double  v_muon_dz = 0.1;
  int v_muon_numberOfMatchedStations =2;
  int v_muon_nPixelValidHits =1.;
  int v_muon_numberOfValidMuonHits = 1;
  double v_muon_chi2Norm= 10;
  int v_muon_nValidHits= 6;
  double v_muon_hadVetoEt= 6.;
  double v_muon_emVetoEt=4.;

  std::vector<const pat::Muon* > vMuons;
  for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
    {
      if ( mu->pt()  < v_muon_pt ) continue;

      if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
      if ( !mu->isGlobalMuon()  ) continue;
      if ( !mu->isPFMuon() ) continue;
      if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim                                                                                                                         
      const reco::TrackRef innerTrack = mu->innerTrack();
      if( innerTrack.isNull() ) continue;
      if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
      if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;

      const reco::TrackRef globalTrack = mu->globalTrack() ;
      if( globalTrack.isNull() ) continue;
      if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
      if(globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
      if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
      if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;

      if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
      if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;



      double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
      double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
      double photonIso = mu->pfIsolationR03().sumPhotonEt;
      double beta = mu->pfIsolationR03().sumPUPt;
      double pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
      if(pfRelIso>0.10) continue;
      vMuons.push_back(&*mu);

    }

  return vMuons;

}
std::vector<const pat::Muon* > tools::MuonIsoSSSelector(const std::vector<pat::Muon>  & thePatMuons,
						     double v_muon_pt,
						     double v_muon_eta,
						     reco::Vertex::Point PV)
{

  double v_muon_d0= 0.005;
  double  v_muon_dz = 0.1;
  int v_muon_numberOfMatchedStations =2;
  int v_muon_nPixelValidHits =1.;
  int v_muon_numberOfValidMuonHits = 1;
  double v_muon_chi2Norm= 10;
  int v_muon_nValidHits= 6;
  double v_muon_hadVetoEt= 6.;
  double v_muon_emVetoEt=4.;

  std::vector<const pat::Muon* > vMuons;
  for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
    {
      if ( mu->pt()  < v_muon_pt ) continue;
      int failID=0;

      if ( TMath::Abs( mu->eta() ) > v_muon_eta ) failID++;
      if ( !mu->isGlobalMuon()  ) failID++;
      if ( !mu->isPFMuon() ) failID++;
      if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) failID++;   //we should add this to skim                                                                                                                        \
                                                                                                                                                                                                                                            
	const reco::TrackRef innerTrack = mu->innerTrack();
	if( innerTrack.isNull() ) continue;
	if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) failID++;
	if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) failID++;

	const reco::TrackRef globalTrack = mu->globalTrack() ;
	if( globalTrack.isNull() ) continue;
	if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) failID++;
	if(globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) failID++;
	if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) failID++;
	if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) failID++;

	if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) failID++;
	if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) failID++;



	double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
	double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
	double photonIso = mu->pfIsolationR03().sumPhotonEt;
	double beta = mu->pfIsolationR03().sumPUPt;
	double pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
	//	if(pfRelIso<0.15) continue;
	if(failID==0) continue;
	vMuons.push_back(&*mu);

    }

return vMuons;

}

std::vector<const pat::Muon* > tools::MuonAntiSSSelector(const std::vector<pat::Muon>  & thePatMuons,
							 double v_muon_pt,
							 double v_muon_eta,
							 reco::Vertex::Point PV)
{

  double v_muon_d0= 0.005;
  double  v_muon_dz = 0.1;
  int v_muon_numberOfMatchedStations =2;
  int v_muon_nPixelValidHits =1.;
  int v_muon_numberOfValidMuonHits = 1;
  double v_muon_chi2Norm= 10;
  int v_muon_nValidHits= 6;
  double v_muon_hadVetoEt= 6.;
  double v_muon_emVetoEt=4.;

  std::vector<const pat::Muon* > vMuons;
  for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
    {
      if ( mu->pt()  < v_muon_pt ) continue;

      if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
      if ( !mu->isGlobalMuon()  ) continue;
      if ( !mu->isPFMuon() ) continue;
      if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim              

      const reco::TrackRef innerTrack = mu->innerTrack();
      if( innerTrack.isNull() ) continue;
      if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
      if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;

      const reco::TrackRef globalTrack = mu->globalTrack() ;
      if( globalTrack.isNull() ) continue;
      if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
      if(globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
      if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
      if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;

      if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
      if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;



      double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
      double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
      double photonIso = mu->pfIsolationR03().sumPhotonEt;
      double beta = mu->pfIsolationR03().sumPUPt;
      double pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
      if(pfRelIso<0.1) continue;
      vMuons.push_back(&*mu);

    }


return vMuons;

}


std::vector<const reco::PFTau *> tools::TauSelector(edm::Handle<reco::PFTauCollection>  & PFTaus,
							  double v_tau_pt,
							  double v_tau_eta,
							  edm::Handle<reco::PFTauDiscriminator> & electron,
							  edm::Handle<reco::PFTauDiscriminator> & muon,
                                                          edm::Handle<reco::PFTauDiscriminator> & iso,
                                                          edm::Handle<reco::PFTauDiscriminator> & decay){
  std::vector<const reco::PFTau *> vTaus;

   for( unsigned  i=0; i<PFTaus->size(); i++ ) {
     //     if((*PFTaus)[i].pt()>10) std::cout<<"Tau "<<(*PFTaus)[i].pt()<<" "<< (*PFTaus)[i].eta()<<" "<< (*PFTaus)[i].phi()<<std::endl;
     //        std::cout<<"pt threshold"<<std::endl;
	if((*PFTaus)[i].pt()<v_tau_pt) continue;
	//	std::cout<<"eta cut"<<std::endl;
	if(TMath::Abs((*PFTaus)[i].eta())>v_tau_eta) continue;
	//	std::cout<<"Input tau "<<(*PFTaus)[i].pt()<<" "<< (*PFTaus)[i].eta()<<" "<< (*PFTaus)[i].phi()<<std::endl;
	reco::PFTauRef tauCandidate(PFTaus, i);
	//	std::cout<<"electron discr."<<std::endl;
        if( (*electron)[tauCandidate] < 0.5 ) continue;
	//	std::cout<<"muon discr."<<std::endl;
        if( (*muon)[tauCandidate] < 0.5 ) continue;
	//	std::cout<<"iso discr."<<std::endl;
        if( (*iso)[tauCandidate] < 0.5 ) continue;
	//	std::cout<<"decay discr."<<std::endl;
        if( (*decay)[tauCandidate] < 0.5 ) continue;
	//	std::cout<<"tau pass"<<std::endl;
	vTaus.push_back(&((*PFTaus)[i]));
   }
  return vTaus;
}

std::vector<const reco::PFTau *> tools::TauAntiSelector(edm::Handle<reco::PFTauCollection>  & PFTaus,
						    double v_tau_pt,
						    double v_tau_eta,
						    edm::Handle<reco::PFTauDiscriminator> & electron,
						    edm::Handle<reco::PFTauDiscriminator> & muon,
						    edm::Handle<reco::PFTauDiscriminator> & iso,
						    edm::Handle<reco::PFTauDiscriminator> & decay){
  std::vector<const reco::PFTau *> vTaus;

  for( unsigned  i=0; i<PFTaus->size(); i++ ) {
  if((*PFTaus)[i].pt()<v_tau_pt) continue;
  if(TMath::Abs((*PFTaus)[i].eta())>v_tau_eta) continue;
  //  std::cout<<"Anti-selected Tau "<<(*PFTaus)[i].pt()<<" "<< (*PFTaus)[i].eta()<<" "<< (*PFTaus)[i].phi()<<std::endl;
  reco::PFTauRef tauCandidate(PFTaus, i);
  int failID=0;
  //  std::cout<<"electron discr."<<std::endl;
  if( (*electron)[tauCandidate] < 0.5 ) failID++;
  //  std::cout<<"muon discr."<<std::endl;
  if( (*muon)[tauCandidate] < 0.5 ) failID++;
  if( (*iso)[tauCandidate] < 0.5 ) failID++;
  if( (*decay)[tauCandidate] < 0.5 ) failID++;
  if(failID==0) continue;
  vTaus.push_back(&((*PFTaus)[i]));
 }
return vTaus;
}


std::vector<const pat::Muon* > tools::MuonEWKinoAntiSelector(const std::vector<pat::Muon>  & thePatMuons,
							 double v_muon_pt,
							 double v_muon_eta,
							 reco::Vertex::Point PV)
{

  //  double v_muon_d0= 0.02;
  double v_muon_d0= 0.2;

  double  v_muon_dz = 0.5;
  int v_muon_numberOfMatchedStations =2;
  int v_muon_nPixelValidHits =1.;
  int v_muon_numberOfValidMuonHits = 1;
  double v_muon_chi2Norm= 10;
  int v_muon_nValidHits= 6;

  std::vector<const pat::Muon* > vMuons;
  for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
    {
      if ( mu->pt()  < v_muon_pt ) continue;

      if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
      if ( !mu->isGlobalMuon()  ) continue;
      if ( !mu->isPFMuon() ) continue;
      if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim                                                                                                                         
      const reco::TrackRef innerTrack = mu->innerTrack();
      if( innerTrack.isNull() ) continue;
      if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
      if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
      if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;

      const reco::TrackRef globalTrack = mu->globalTrack() ;
      if( globalTrack.isNull() ) continue;
      if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
      if(globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
      if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
      if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;


      //      double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
      //      double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
      //      double photonIso = mu->pfIsolationR03().sumPhotonEt;
      //      double beta = mu->pfIsolationR03().sumPUPt;
      //      double pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
      //      if(pfRelIso<0.15) continue;

      double chargedHadronIso = mu->pfIsolationR04().sumChargedHadronPt;
      double neutralHadronIso = mu->pfIsolationR04().sumNeutralHadronEt;
      double photonIso = mu->pfIsolationR04().sumPhotonEt;
      double beta = mu->pfIsolationR04().sumPUPt;
      double pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
      if(pfRelIso<0.2) continue;


      // if(pfRelIso<0.1) continue;  
      vMuons.push_back(&*mu);

    }

  return vMuons;

}

std::vector<const pat::Muon* > tools::MuonEWKinoIsoSelector(const std::vector<pat::Muon>  & thePatMuons,
							     double v_muon_pt,
							     double v_muon_eta,
							     reco::Vertex::Point PV)
{

  //  double v_muon_d0= 0.02;
  double v_muon_d0= 0.2;
  double  v_muon_dz = 0.5;
  int v_muon_numberOfMatchedStations =2;
  int v_muon_nPixelValidHits =1.;
  int v_muon_numberOfValidMuonHits = 1;
  double v_muon_chi2Norm= 10;
  int v_muon_nValidHits= 6;

  std::vector<const pat::Muon* > vMuons;
  for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
    {
      if(!mu->isGlobalMuon() && !mu->isTrackerMuon()) continue;
      if ( mu->pt()  < v_muon_pt ) continue;
      if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
      int failID=0;
      if ( !mu->isGlobalMuon()  ) failID++;
      if ( !mu->isPFMuon() ) failID++;
      if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) failID++;   //we should add this to skim                                                                                                
      const reco::TrackRef innerTrack = mu->innerTrack();
      if( innerTrack.isNull() ) failID++;
      else if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) failID++;
      else if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) failID++;
      else if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) failID++;
      else if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) failID++;
      else if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) failID++;

      const reco::TrackRef globalTrack = mu->globalTrack() ;
      if( globalTrack.isNull() ) failID++;
      else if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) failID++;
      else if(globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) failID++;


      double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
      double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
      double photonIso = mu->pfIsolationR03().sumPhotonEt;
      double beta = mu->pfIsolationR03().sumPUPt;
      double pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
      //      if(pfRelIso>0.15) continue;
      if(failID==0) continue;
      // if(pfRelIso<0.1) continue;                                                                                                                                                                                                          
      vMuons.push_back(&*mu);

    }

  return vMuons;

}

std::vector<const pat::Muon* > tools::MuonAntiSelector(const std::vector<pat::Muon>  & thePatMuons,
							    double v_muon_pt,
							    double v_muon_eta,
							    reco::Vertex::Point PV)
{

  double v_muon_d0= 0.02;
  double  v_muon_dz = 0.5;
  int v_muon_numberOfMatchedStations =2;
  int v_muon_nPixelValidHits =1.;
  int v_muon_numberOfValidMuonHits = 1;
  double v_muon_chi2Norm= 10;
  int v_muon_nValidHits= 6;

  std::vector<const pat::Muon* > vMuons;
  for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
    {
      if ( mu->pt()  < v_muon_pt ) continue;
      if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
      int failID=0;
      if ( !mu->isGlobalMuon()  ) failID++;
      if ( !mu->isPFMuon() ) failID++;
      if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) failID++;   //we should add this to skim                                                                                                                         
      const reco::TrackRef innerTrack = mu->innerTrack();
      if( innerTrack.isNull() ) failID++;
      else if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) failID++;
      else if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) failID++;
      else if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) failID++;
      else if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) failID++;
      else if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) failID++;
      const reco::TrackRef globalTrack = mu->globalTrack() ;
      if( globalTrack.isNull() ) failID++;
      else if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) failID++;
      else if(globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) failID++;


      double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
      double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
      double photonIso = mu->pfIsolationR03().sumPhotonEt;
      double beta = mu->pfIsolationR03().sumPUPt;
      double pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
      if(pfRelIso<0.15 && failID==0) continue;
      vMuons.push_back(&*mu);

    }

  return vMuons;

}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
bool bool_electron_ecalDriven= true;

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

	  float ecalIso = TMath::Abs(el->eta()) > 1.47 ? el->dr03EcalRecHitSumEt() : TMath::Max(el->dr03EcalRecHitSumEt()-1.,0.); 
	  float detRelIso = (el->dr03TkSumPt() + el->dr03HcalTowerSumEt() + ecalIso ) / el->pt() ;
           //puChargedHadronIso(), particleIso()
          //float pfRelIso  = (el->chargedHadronIso() + el->neutralHadronIso() + el->photonIso() ) /el->pt() ;
          float RelIso = detRelIso; 
	  //if(usePFiso) RelIso = pfRelIso;  
	  //	  if( RelIso  > v_electron_reliso ) continue;
	 
	  std::cout<<"Electron pass"<<std::endl;
	  vElectrons.push_back(&*el );
	}
      return vElectrons;
}
std::vector<const pat::Electron* > tools::ElectronEWKinoSelector(const std::vector<pat::Electron>  & thePatElectrons,
							   double v_electron_pt,
							   double v_electron_eta,
							   double myRho,
							   reco::Vertex::Point PV,
							   edm::Handle< std::vector<reco::Conversion> > &theConversions,
							   reco::BeamSpot::Point BS)
{

     double v_electron_d0= 0.02;
  //  double v_electron_d0= 10.02;
  double v_electron_dz = 0.2;
  //  double v_electron_dz = 20.2;

  std::vector<const pat::Electron* > vElectrons;
  for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
    {
      const reco::GsfTrackRef gsfTrack = el->gsfTrack();
      if( el->pt() < v_electron_pt ) continue;
      if( TMath::Abs(el->eta()) > v_electron_eta ) continue;                           

      if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
      if( TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  continue;
      if( TMath::Abs(el->gsfTrack()->dz(PV))  > v_electron_dz  ) continue;

      if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
      bool vtxFitConversion = ConversionTools::hasMatchedConversion( reco::GsfElectron (*el) , theConversions, BS);
      if( vtxFitConversion )  continue;
      if( el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > 1 ) continue;

      if( TMath::Abs( el->eta() ) < 1.5 )
	{
	  if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.15 ) continue;
	  if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
          if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
	  if( TMath::Abs(el->hadronicOverEm()) > 0.12 ) continue;  //recommended is 0.12 but HLT applies 0.1                                                                                                                             
	}
      else
	{
	  if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.10 ) continue;
	  if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.009 ) continue;
	  if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
	  if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) continue;   /// at the HLT 0.075  recommended is 0.1                                                                                                                              
	}
      //      double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13 };
      double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14 };


      double CorrectedTerm;
      if( TMath::Abs( el->superCluster()->eta() ) < 1.0)    CorrectedTerm = myRho * Aeff[ 0 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.0 && TMath::Abs( el->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.479 && TMath::Abs( el->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.0 && TMath::Abs( el->superCluster()->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.2 && TMath::Abs( el->superCluster()->eta() ) < 2.3  )     CorrectedTerm = myRho * Aeff[ 4 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.3 && TMath::Abs( el->superCluster()->eta() ) < 2.4  )     CorrectedTerm = myRho * Aeff[ 5 ];
      else CorrectedTerm = myRho * Aeff[ 6 ];

      float pfRelIso = (el->chargedHadronIso() + TMath::Max(0.0, el->neutralHadronIso() + el->photonIso() - CorrectedTerm ) ) /el->pt() ;
      if(pfRelIso>0.15) continue;
      //      if(pfRelIso>0.1) continue;
      vElectrons.push_back(&*el );
    }
  return vElectrons;
}


std::vector<const pat::Electron* > tools::ElectronEWKinoConvSelector(const std::vector<pat::Electron>  & thePatElectrons,
								 double v_electron_pt,
								 double v_electron_eta,
								 double myRho,
								 reco::Vertex::Point PV,
								 edm::Handle< std::vector<reco::Conversion> > &theConversions,
								 reco::BeamSpot::Point BS)
{


  double v_electron_d0= 0.02;
  double v_electron_dz = 0.2;

  std::vector<const pat::Electron* > vElectrons;
  for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
    {
      const reco::GsfTrackRef gsfTrack = el->gsfTrack();
      if( el->pt() < v_electron_pt ) continue;
      if( TMath::Abs(el->superCluster()->eta()) > v_electron_eta ) continue;

      if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
      if( TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  continue;
      if( TMath::Abs(el->gsfTrack()->dz(PV))  > v_electron_dz  ) continue;

      if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
      bool vtxFitConversion = ConversionTools::hasMatchedConversion( reco::GsfElectron (*el) , theConversions, BS);
      if( !vtxFitConversion && el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() <2 ) continue;

      if( TMath::Abs( el->eta() ) < 1.5 )
        {
          if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.15 ) continue;
          if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
          if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
          if( TMath::Abs(el->hadronicOverEm()) > 0.12 ) continue;  //recommended is 0.12 but HLT applies 0.1                                                                                                                                 
        }
      else
        {
          if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.10 ) continue;
          if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.009 ) continue;
          if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
          if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) continue;   /// at the HLT 0.075  recommended is 0.1                                                                                                                                  
        }
//      double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13 };
      double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14 };

      double CorrectedTerm;
      if( TMath::Abs( el->superCluster()->eta() ) < 1.0)    CorrectedTerm = myRho * Aeff[ 0 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.0 && TMath::Abs( el->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.479 && TMath::Abs( el->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.0 && TMath::Abs( el->superCluster()->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.2 && TMath::Abs( el->superCluster()->eta() ) < 2.3  )     CorrectedTerm = myRho * Aeff[ 4 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.3 && TMath::Abs( el->superCluster()->eta() ) < 2.4  )     CorrectedTerm = myRho * Aeff[ 5 ];
      else CorrectedTerm = myRho * Aeff[ 6 ];

      float pfRelIso = (el->chargedHadronIso() + TMath::Max(0.0, el->neutralHadronIso() + el->photonIso() - CorrectedTerm ) ) /el->pt() ;
      if(pfRelIso>0.15) continue;
      //      if(pfRelIso>0.1) continue;                                                                                                                                                                                                     
      vElectrons.push_back(&*el );
    }
  return vElectrons;
}


std::vector<const pat::Electron* > tools::ElectronSSSelector(const std::vector<pat::Electron>  & thePatElectrons,
								 double v_electron_pt,
								 double v_electron_eta,
								 double myRho,
								 reco::Vertex::Point PV,
								 edm::Handle< std::vector<reco::Conversion> > &theConversions,
								 reco::BeamSpot::Point BS)
{

  double v_electron_d0= 0.01;
  double v_electron_dz = 0.1;

  std::vector<const pat::Electron* > vElectrons;
  for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
    {
      const reco::GsfTrackRef gsfTrack = el->gsfTrack();
      if( el->pt() < v_electron_pt ) continue;
      if( TMath::Abs(el->superCluster()->eta()) > v_electron_eta ) continue;

      if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
      if( TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  continue;
      if( TMath::Abs(el->gsfTrack()->dz(PV))  > v_electron_dz  ) continue;

      if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
      bool vtxFitConversion = ConversionTools::hasMatchedConversion( reco::GsfElectron (*el) , theConversions, BS);
      if( vtxFitConversion )  continue;
      //      if(!el->isGsfCtfScPixChargeConsistent() )  continue;            

      if( el->pt() < 20. )
	{
	  if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	}

      if( el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;

      if( TMath::Abs( el->superCluster()->eta() ) < 1.4442 )
	{
	  if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06 ) continue;
	  if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
	  if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
	  if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) continue;  //recommended is 0.12 but HLT applies 0.1                                                                                                                              
	}
      else if( TMath::Abs( el->superCluster()->eta() ) < 2.4 )
	{
	  if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
	  if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
	  if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
	  if( TMath::Abs(el->hadronicOverEm()) > 0.075 ) continue;   /// at the HLT 0.075  recommended is 0.1                                                                                                                            
	}


      //      double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13 };
      double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14 };

      double CorrectedTerm;
      if( TMath::Abs( el->eta() ) < 1.0)    CorrectedTerm = myRho * Aeff[ 0 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.0 && TMath::Abs( el->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.479 && TMath::Abs( el->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.0 && TMath::Abs( el->superCluster()->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.2 && TMath::Abs( el->superCluster()->eta() ) < 2.3  )     CorrectedTerm = myRho * Aeff[ 4 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.3 && TMath::Abs( el->superCluster()->eta() ) < 2.4  )     CorrectedTerm = myRho * Aeff[ 5 ];
      else CorrectedTerm = myRho * Aeff[ 6 ];

      float pfRelIso = (el->chargedHadronIso() + TMath::Max(0.0, el->neutralHadronIso() + el->photonIso() - CorrectedTerm ) ) /el->pt() ;
      //       if(pfRelIso>0.15) continue;
       if(pfRelIso>0.09) continue;                                                                                                                                                                                                     
      vElectrons.push_back(&*el );
    }
  return vElectrons;
}

std::vector<const pat::Electron* > tools::ElectronIsoSSSelector(const std::vector<pat::Electron>  & thePatElectrons,
							     double v_electron_pt,
							     double v_electron_eta,
							     double myRho,
							     reco::Vertex::Point PV,
							     edm::Handle< std::vector<reco::Conversion> > &theConversions,
							     reco::BeamSpot::Point BS)
{

  double v_electron_d0= 0.01;
  double v_electron_dz = 0.1;

  std::vector<const pat::Electron* > vElectrons;
  for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
    {
      const reco::GsfTrackRef gsfTrack = el->gsfTrack();
      if( el->pt() < v_electron_pt ) continue;
      if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
      int failID=0;

      if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
      if( TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  failID++;
      if( TMath::Abs(el->gsfTrack()->dz(PV))  > v_electron_dz  ) failID++;

      if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) failID++;
      bool vtxFitConversion = ConversionTools::hasMatchedConversion( reco::GsfElectron (*el) , theConversions, BS);
      if( vtxFitConversion )  failID++;
      //      if(!el->isGsfCtfScPixChargeConsistent() )  failID++;

      if( el->pt() < 20. )
	{
          if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) failID++;
        }

      if( el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > 0 ) failID++;

      if( TMath::Abs( el->eta() ) < 1.4442 )
        {
          if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06 ) failID++;
          if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) failID++;
          if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) failID++;
          if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) failID++;  //recommended is 0.12 but HLT applies 0.1                                                                                                                                  
        }
      else if( TMath::Abs( el->eta() ) < 2.4 )
        {
          if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) failID++;
          if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) failID++;
          if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) failID++;
          if( TMath::Abs(el->hadronicOverEm()) > 0.075 ) failID++;   /// at the HLT 0.075  recommended is 0.1                                                                                                                                
        }


      //      double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13 };
      double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14 };

      double CorrectedTerm;
      if( TMath::Abs( el->eta() ) < 1.0)    CorrectedTerm = myRho * Aeff[ 0 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.0 && TMath::Abs( el->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.479 && TMath::Abs( el->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.0 && TMath::Abs( el->superCluster()->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.2 && TMath::Abs( el->superCluster()->eta() ) < 2.3  )     CorrectedTerm = myRho * Aeff[ 4 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.3 && TMath::Abs( el->superCluster()->eta() ) < 2.4  )     CorrectedTerm = myRho * Aeff[ 5 ];

      float pfRelIso = (el->chargedHadronIso() + TMath::Max(0.0, el->neutralHadronIso() + el->photonIso() - CorrectedTerm ) ) /el->pt() ;
      if(failID==0) continue;
      //     if(pfRelIso>0.15) continue;                                                                                                                                                                                                     
      //if(pfRelIso<0.15) continue;
      vElectrons.push_back(&*el );
    }
  return vElectrons;
}

std::vector<const pat::Electron* > tools::ElectronAntiSSSelector(const std::vector<pat::Electron>  & thePatElectrons,
								 double v_electron_pt,
								 double v_electron_eta,
								 double myRho,
								 reco::Vertex::Point PV,
								 edm::Handle< std::vector<reco::Conversion> > &theConversions,
								 reco::BeamSpot::Point BS)
{

  double v_electron_d0= 0.01;
  double v_electron_dz = 0.1;

  std::vector<const pat::Electron* > vElectrons;
  for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
    {
      const reco::GsfTrackRef gsfTrack = el->gsfTrack();
      if( el->pt() < v_electron_pt ) continue;
      if( TMath::Abs(el->eta()) > v_electron_eta ) continue;

      if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
      if( TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  continue;
      if( TMath::Abs(el->gsfTrack()->dz(PV))  > v_electron_dz  ) continue;

      if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
      bool vtxFitConversion = ConversionTools::hasMatchedConversion( reco::GsfElectron (*el) , theConversions, BS);
      if( vtxFitConversion )  continue;
      if(!el->isGsfCtfScPixChargeConsistent() )  continue;

      if( el->pt() < 20. )
        {
          if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	  }

	  if( el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;

	  if( TMath::Abs( el->eta() ) < 1.4442 )
	    {
	      if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06 ) continue;
	      if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
	      if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
	      if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) continue;  //recommended is 0.12 but HLT applies 0.1                                                                                                                                 \
                                                                                                                                                                                                                                             
		}
	  else if( TMath::Abs( el->eta() ) < 2.4 )
	    {
	      if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
	      if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
	      if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
	      if( TMath::Abs(el->hadronicOverEm()) > 0.075 ) continue;   /// at the HLT 0.075  recommended is 0.1                                                                                                                               \
                                                                                                                                                                                                                                             
		}


	  //	  double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13 };
	  double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14 };

	  double CorrectedTerm;
	  if( TMath::Abs( el->eta() ) < 1.0)    CorrectedTerm = myRho * Aeff[ 0 ];
	  else if( TMath::Abs( el->superCluster()->eta() ) > 1.0 && TMath::Abs( el->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
	  else if( TMath::Abs( el->superCluster()->eta() ) > 1.479 && TMath::Abs( el->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
	  else if( TMath::Abs( el->superCluster()->eta() ) > 2.0 && TMath::Abs( el->superCluster()->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
	  else if( TMath::Abs( el->superCluster()->eta() ) > 2.2 && TMath::Abs( el->superCluster()->eta() ) < 2.3  )     CorrectedTerm = myRho * Aeff[ 4 ];
	  else if( TMath::Abs( el->superCluster()->eta() ) > 2.3 && TMath::Abs( el->superCluster()->eta() ) < 2.4  )     CorrectedTerm = myRho * Aeff[ 5 ];
	  else CorrectedTerm = myRho * Aeff[ 6 ];

	  float pfRelIso = (el->chargedHadronIso() + TMath::Max(0.0, el->neutralHadronIso() + el->photonIso() - CorrectedTerm ) ) /el->pt() ;
	  //     if(pfRelIso>0.15) continue;                                                                                                                                                                                                    \
                                                                                                                                                                                                                                             
	  if(pfRelIso<0.09) continue;
	  vElectrons.push_back(&*el );
    }
     return vElectrons;
}


std::vector<const pat::Electron* > tools::ElectronEWKinoAntiSelector(const std::vector<pat::Electron>  & thePatElectrons,
								 double v_electron_pt,
								 double v_electron_eta,
								 double myRho,
								 reco::Vertex::Point PV,
								 edm::Handle< std::vector<reco::Conversion> > &theConversions,
								 reco::BeamSpot::Point BS)
{

  double v_electron_d0= 0.02;
  double v_electron_dz = 0.2;

  std::vector<const pat::Electron* > vElectrons;
  for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
    {
      const reco::GsfTrackRef gsfTrack = el->gsfTrack();
      if( el->pt() < v_electron_pt ) continue;
      if( TMath::Abs(el->eta()) > v_electron_eta ) continue;

      if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
      if( TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  continue;
      if( TMath::Abs(el->gsfTrack()->dz(PV))  > v_electron_dz  ) continue;

      if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
      bool vtxFitConversion = ConversionTools::hasMatchedConversion( reco::GsfElectron (*el) , theConversions, BS);
      if( vtxFitConversion )  continue;
      if( el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > 1 ) continue;
      if( TMath::Abs( el->eta() ) < 1.5 )
        {
          if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.15 ) continue;
          if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
          if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
          if( TMath::Abs(el->hadronicOverEm()) > 0.12 ) continue;  //recommended is 0.12 but HLT applies 0.1                                                                                                                                 
        }
      else
        {
          if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.10 ) continue;
          if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.009 ) continue;
          if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
          if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) continue;   /// at the HLT 0.075  recommended is 0.1                                                                                                                                  
        }
      //      double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13 };
      double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14 };

      double CorrectedTerm;
      if( TMath::Abs( el->superCluster()->eta() ) < 1.0)    CorrectedTerm = myRho * Aeff[ 0 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.0 && TMath::Abs( el->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.479 && TMath::Abs( el->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.0 && TMath::Abs( el->superCluster()->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.2 && TMath::Abs( el->superCluster()->eta() ) < 2.3  )     CorrectedTerm = myRho * Aeff[ 4 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.3 && TMath::Abs( el->superCluster()->eta() ) < 2.4  )     CorrectedTerm = myRho * Aeff[ 5 ];
      else CorrectedTerm = myRho * Aeff[ 6 ];

      float pfRelIso = (el->chargedHadronIso() + TMath::Max(0.0, el->neutralHadronIso() + el->photonIso() - CorrectedTerm ) ) /el->pt() ;
      if(pfRelIso<0.15) continue;

      vElectrons.push_back(&*el );
    }
  return vElectrons;
}

std::vector<const pat::Electron* > tools::ElectronEWKinoIsoSelector(const std::vector<pat::Electron>  & thePatElectrons,
								     double v_electron_pt,
								     double v_electron_eta,
								     double myRho,
								     reco::Vertex::Point PV,
								     edm::Handle< std::vector<reco::Conversion> > &theConversions,
								     reco::BeamSpot::Point BS)
{

  double v_electron_d0= 0.02;
  double v_electron_dz = 0.2;

  std::vector<const pat::Electron* > vElectrons;
  for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
    {
      const reco::GsfTrackRef gsfTrack = el->gsfTrack();
      if( el->pt() < v_electron_pt ) continue;
      if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
     
      if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;

      int failID=0;
      if( TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  failID++;
      if( TMath::Abs(el->gsfTrack()->dz(PV))  > v_electron_dz  ) failID++;

      if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) failID++;
      bool vtxFitConversion = ConversionTools::hasMatchedConversion( reco::GsfElectron (*el) , theConversions, BS);
      if( vtxFitConversion )  failID++;
      if( el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > 1 ) failID++;
      if( TMath::Abs( el->eta() ) < 1.5 )
        {
          if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.15 ) failID++;
          if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) failID++;
          if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) failID++;
          if( TMath::Abs(el->hadronicOverEm()) > 0.12 ) failID++;  //recommended is 0.12 but HLT applies 0.1                                                                                                                              
        }
      else
        {
          if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.10 ) failID++;
          if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.009 ) failID++;
          if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) failID++;
          if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) failID++;   /// at the HLT 0.075  recommended is 0.1                                                                                                                               
	}

      //      double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13 };
      double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14 };

      double CorrectedTerm;
      if( TMath::Abs( el->superCluster()->eta() ) < 1.0)    CorrectedTerm = myRho * Aeff[ 0 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.0 && TMath::Abs( el->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.479 && TMath::Abs( el->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.0 && TMath::Abs( el->superCluster()->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.2 && TMath::Abs( el->superCluster()->eta() ) < 2.3  )     CorrectedTerm = myRho * Aeff[ 4 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.3 && TMath::Abs( el->superCluster()->eta() ) < 2.4  )     CorrectedTerm = myRho * Aeff[ 5 ];
      else CorrectedTerm = myRho * Aeff[ 6 ];

      float pfRelIso = (el->chargedHadronIso() + TMath::Max(0.0, el->neutralHadronIso() + el->photonIso() - CorrectedTerm ) ) /el->pt() ;
      //      if(pfRelIso>0.15) continue;
      if(failID==0) continue;

      vElectrons.push_back(&*el );
    }
  return vElectrons;
}

std::vector<const pat::Electron* > tools::ElectronAntiSelector(const std::vector<pat::Electron>  & thePatElectrons,
								    double v_electron_pt,
								    double v_electron_eta,
								    double myRho,
								    reco::Vertex::Point PV,
								    edm::Handle< std::vector<reco::Conversion> > &theConversions,
								    reco::BeamSpot::Point BS)
{

  double v_electron_d0= 0.02;
  double v_electron_dz = 0.2;

  std::vector<const pat::Electron* > vElectrons;
  for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
    {
      const reco::GsfTrackRef gsfTrack = el->gsfTrack();
      if( el->pt() < v_electron_pt ) continue;
      if( TMath::Abs(el->eta()) > v_electron_eta ) continue;

      if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;

      int failID=0;
      if( TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  failID++;
      if( TMath::Abs(el->gsfTrack()->dz(PV))  > v_electron_dz  ) failID++;

      if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) failID++;
      bool vtxFitConversion = ConversionTools::hasMatchedConversion( reco::GsfElectron (*el) , theConversions, BS);
      if( vtxFitConversion )  failID++;
      if( el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > 1 ) failID++;
      if( TMath::Abs( el->eta() ) < 1.5 )
        {
          if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.15 ) failID++;
          if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) failID++;
          if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) failID++;
          if( TMath::Abs(el->hadronicOverEm()) > 0.12 ) failID++;  //recommended is 0.12 but HLT applies 0.1                                                                                                                                 
        }
      else
        {
          if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.10 ) failID++;
          if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.009 ) failID++;
          if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) failID++;
          if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) failID++;   /// at the HLT 0.075  recommended is 0.1                                                                                                                                  
        }

      //      double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13 };                                                                                                                                                             
      double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14 };

      double CorrectedTerm;
      if( TMath::Abs( el->superCluster()->eta() ) < 1.0)    CorrectedTerm = myRho * Aeff[ 0 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.0 && TMath::Abs( el->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 1.479 && TMath::Abs( el->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.0 && TMath::Abs( el->superCluster()->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.2 && TMath::Abs( el->superCluster()->eta() ) < 2.3  )     CorrectedTerm = myRho * Aeff[ 4 ];
      else if( TMath::Abs( el->superCluster()->eta() ) > 2.3 && TMath::Abs( el->superCluster()->eta() ) < 2.4  )     CorrectedTerm = myRho * Aeff[ 5 ];
      else CorrectedTerm = myRho * Aeff[ 6 ];

      float pfRelIso = (el->chargedHadronIso() + TMath::Max(0.0, el->neutralHadronIso() + el->photonIso() - CorrectedTerm ) ) /el->pt() ;
      if(pfRelIso<0.15 && failID==0) continue;

      vElectrons.push_back(&*el );
    }
  return vElectrons;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Jet* > tools::JetSelector(const std::vector<pat::Jet>  & thePatJets, 
double  v_jet_et,  
double  v_jet_eta, 
bool    bool_jet_id, 
bool    jetLeptonCleaning, 
double  v_jet_leptonVetoDR,
std::vector<const pat::Electron*> & vElectrons,
std::vector<const pat::Muon*> & vMuons,
std::vector<const reco::PFTau *> & vTaus
)
{


     std::vector< const pat::Jet* > vJets;

      for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ ) 
	{
	  //	  if(jet->pt()>10)   std::cout<<"Jet "<<jet->pt()<<" "<<jet->eta()<<" "<<jet->phi()<<" raw "<<jet->correctedJet("Uncorrected").pt()<<std::endl;

	  if( jet->pt() < v_jet_et )continue;
          if( TMath::Abs( jet->eta() ) > v_jet_eta ) continue;
	  //  std::cout<<"Jet ID"<<std::endl;
	  if( bool_jet_id )
	    {
	      //  std::cout<<"Neutral Hadronic energy Fraction"<<std::endl;
	      if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
	      //  std::cout<<"Neutral EM energy Fraction"<<std::endl;
	      if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
	      //    std::cout<<"Multiplicity"<<std::endl;
	      if(jet->nConstituents() < 2 ) continue;
	      if( TMath::Abs( jet->eta() ) < 2.4 ) 
		{
		  //		  std::cout<<"Charged Hadronic energy Fraction"<<std::endl;
		  if( jet->chargedHadronEnergyFraction() == 0. ) continue;
		  //	  std::cout<<"Charged EM energy Fraction"<<std::endl;
		  if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
		  //	  std::cout<<"Charged Multiplicity"<<std::endl;
		  if( jet->chargedMultiplicity() == 0 ) continue;
		}
	    }

	  bool vetoJet = false;
          for(unsigned int i = 0 ; i < vMuons.size() ;i++ ) 
            {
             const pat::Muon *mu = vMuons[i]; 
	      float dphi = TMath::ACos( TMath::Cos( mu->phi()-jet->phi() ));
	      float deta = mu->eta()-jet->eta();
	      float dr = TMath::Sqrt( dphi*dphi + deta*deta) ;
      	    
	      if(jetLeptonCleaning && dr < v_jet_leptonVetoDR)
		{
		  vetoJet = true;
		  break;
		}
	    }
          for(unsigned int i = 0 ; i < vElectrons.size() ;i++ ) 
            {
             const pat::Electron *el = vElectrons[i]; 
	      float dphi = TMath::ACos( TMath::Cos( el->phi()-jet->phi() ) );
	      float deta = el->eta() - jet->eta();
	      float dr = TMath::Sqrt( dphi*dphi + deta*deta );
	      
	      if( jetLeptonCleaning && dr < v_jet_leptonVetoDR) 
		{
		  vetoJet = true;
		  break;
		}
	    }

          for(unsigned int i = 0 ; i < vTaus.size() ;i++ )
            {
	      const reco::PFTau *tau = vTaus[i];
              float dphi = TMath::ACos( TMath::Cos( tau->phi()-jet->phi() ) );
              float deta = tau->eta() - jet->eta();
              float dr = TMath::Sqrt( dphi*dphi + deta*deta );

              if( jetLeptonCleaning && dr < v_jet_leptonVetoDR)
                {
                  vetoJet = true;
                  break;
                }
            }

	  // std::cout<<"Jet cleaning"<<std::endl;
	  if( vetoJet ) continue;	  
	  //  std::cout<<"Jet passes"<<std::endl;
	  vJets.push_back( &*jet );
      }
	  return vJets;
}
bool tools::PassEventAcceptance2(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons){
  int n20=0;
  int n10=0;
  for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
    const pat::Electron *ele1 = thePatElectrons[i];
    if(ele1->pt()>20){
      n20++;
      n10++;
    }
    else if(ele1->pt()>10) n10++;
  }
  for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
    const pat::Muon *mu1 = thePatMuons[i];
    if(mu1->pt()>20){
      n20++;
      n10++;
    }
    else if(mu1->pt()>10) n10++;
  }
  if(n20>0) return kTRUE;
  else return kFALSE;
}


bool tools::PassEventAcceptance(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons){
  int n20=0;
  int n10=0;
  for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
    const pat::Electron *ele1 = thePatElectrons[i];
    if(ele1->pt()>20){
      n20++;
      n10++;
    }
    else if(ele1->pt()>10) n10++;
  }
  for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
    const pat::Muon *mu1 = thePatMuons[i];
    if(mu1->pt()>20){
      n20++;
      n10++;
    }
    else if(mu1->pt()>10) n10++;
  }
  if(n20>0 && n10>1) return kTRUE;
  else return kFALSE;
}


int tools::OSSF(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons, std::vector<const reco::PFTau *> & vTaus){
  TLorentzVector lep1, lep2, dilep, lep3;  // lepton 4-vectors      
  int nele=thePatElectrons.size();
  int nmu=thePatMuons.size() ;
  int ntau=vTaus.size();
    int nelecharge=0;
    int nmucharge=0;
    int ntaucharge=0;
    int nlepcharge=0;
      for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
	const pat::Electron *ele1 = thePatElectrons[i];
	nelecharge+=ele1->charge();
        nlepcharge+=ele1->charge();
      }
      for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
	const pat::Muon *mu1 = thePatMuons[i];
        nmucharge+=mu1->charge();
        nlepcharge+=mu1->charge();
      }
      for(unsigned int i = 0 ; i < vTaus.size() ;i++ ) {
        const reco::PFTau *tau1 = vTaus[i];
        ntaucharge+=tau1->charge();
        nlepcharge+=tau1->charge();
      }
      if(ntau==1 && fabs(nlepcharge)<(nele+nmu+ntau)-1 && fabs(nelecharge+nmucharge)==(nele+nmu)) return 3;
      else if((nele>1 && fabs(nelecharge)<nele-1) || (nmu>1 && fabs(nmucharge)<nmu-1)) return 1;
      else if(nele+nmu>1 && fabs(nelecharge+nmucharge)<(nele+nmu)-1) return 2;
      else return 0;
}

double tools::TriLeptonMass(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons, std::vector<const reco::PFTau *> & vTaus){
  TLorentzVector lep1, lep2, trilep, lep3;  // lepton 4-vectors                                                                                                                                                                               
  for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
    const pat::Electron *ele1 = thePatElectrons[i];
    lep1.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
    trilep=trilep+lep1;
  }
  for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
    const pat::Muon *mu1 = thePatMuons[i];
    lep1.SetPtEtaPhiE(mu1->pt(),mu1->eta(),mu1->phi(),mu1->energy());
    trilep=trilep+lep1;  
  }
  for(unsigned int i = 0 ; i < vTaus.size() ;i++ ) {
    const reco::PFTau *tau1 = vTaus[i];
    lep1.SetPtEtaPhiE(tau1->pt(),tau1->eta(),tau1->phi(),tau1->energy());
    trilep=trilep+lep1;
  }

  return trilep.M();
}


double tools::ZMass(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons, std::vector<const reco::PFTau *> & vTaus,int OSSF){
  TLorentzVector lep1, lep2, dilep, lep3;  // lepton 4-vectors                                                                                                                                                                               
  double Zmass=10000;
  if(OSSF==1){
  for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
      const pat::Electron *ele1 = thePatElectrons[i];
      for(unsigned int j = i+1 ; j < thePatElectrons.size() ;j++ ) {
	const pat::Electron *ele2 = thePatElectrons[j];
	if(ele2->charge()==ele1->charge()) continue;
	lep1.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
        lep2.SetPtEtaPhiE(ele2->pt(),ele2->eta(),ele2->phi(),ele2->energy());
	dilep=lep1+lep2;
	if(fabs(dilep.M()-91)<fabs(Zmass-91)){
	  Zmass=dilep.M();
	}
      }
  }

  for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
    const pat::Muon *mu1 = thePatMuons[i];
    for(unsigned int j = i+1 ; j < thePatMuons.size() ;j++ ) {
      const pat::Muon *mu2 = thePatMuons[j];
      if(mu2->charge()==mu1->charge()) continue;
      lep1.SetPtEtaPhiE(mu1->pt(),mu1->eta(),mu1->phi(),mu1->energy());
      lep2.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
      dilep=lep1+lep2;
      if(fabs(dilep.M()-91)<fabs(Zmass-91)){
	Zmass=dilep.M();
      }
    }
  }
  }
  else if(OSSF==2){
   double difference=10000;
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
      const pat::Electron *ele1 = thePatElectrons[i];
      for(unsigned int j = 0 ; j < thePatMuons.size() ;j++ ) {
	const pat::Muon *mu2 = thePatMuons[j];
	if(mu2->charge()==ele1->charge()) continue;
	lep1.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
	lep2.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
	dilep=lep1+lep2;
	if(fabs(dilep.M()-50)<difference){
          difference=fabs(dilep.M()-50);
	  Zmass=dilep.M();
	}
      }
    }
    if(vTaus.size()==1){
	const reco::PFTau *tau = vTaus[0];
	lep1.SetPtEtaPhiE(tau->pt(),tau->eta(),tau->phi(),tau->energy());
        for(unsigned int i=0;i<thePatElectrons.size();i++){
	  const pat::Electron *ele1 = thePatElectrons[i];
	  if(tau->charge()==ele1->charge()) continue;
	  lep2.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
	  dilep=lep1+lep2;
         if(fabs(dilep.M()-60)<difference){
          difference=fabs(dilep.M()-60);
          Zmass=dilep.M();
        }
	}
     for(unsigned int j = 0 ;j < thePatMuons.size();j++){
       const pat::Muon * mu2=thePatMuons[j];
        if(mu2->charge()==tau->charge()) continue;
        lep2.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
        dilep=lep1+lep2;
        if(fabs(dilep.M()-60)<difference){
	  difference=fabs(dilep.M()-60);
          Zmass=dilep.M();
	}
    }

  }

  }
    else if(OSSF==3){
      const reco::PFTau * tau=vTaus[0];
      for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
	const pat::Electron *ele1 = thePatElectrons[i];
	if(tau->charge()==ele1->charge()) continue;
	lep1.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
	lep2.SetPtEtaPhiE(tau->pt(),tau->eta(),tau->phi(),tau->energy());
	dilep=lep1+lep2;
	if(fabs(dilep.M()-60)<fabs(Zmass-60)){
	  Zmass=dilep.M();
	}
      }
      for(unsigned int j = 0 ; j < thePatMuons.size() ;j++ ) {
	const pat::Muon *mu2 = thePatMuons[j];
	 if(mu2->charge()==tau->charge()) continue;
	  lep1.SetPtEtaPhiE(tau->pt(),tau->eta(),tau->phi(),tau->energy());
	  lep2.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
	  dilep=lep1+lep2;
	  if(fabs(dilep.M()-60)<fabs(Zmass-60)){
	    Zmass=dilep.M();
	  }
	}
    }

  return Zmass;
}

std::vector<TLorentzVector> tools::LeptonSelector(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons, std::vector<const reco::PFTau *> & vTaus, int OSSF){
  TLorentzVector lep1, lep2, dilep, lep3;  // lepton 4-vectors                                                   
  TLorentzVector lep1saved, lep2saved, lep3saved;  // lepton 4-vectors 
  std::vector<TLorentzVector> SelectedLeptons;                 
  double difference=10000;

  if(OSSF==0) return SelectedLeptons;
  else if(OSSF==1){
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
      const pat::Electron *ele1 = thePatElectrons[i];
      for(unsigned int j = i+1 ; j < thePatElectrons.size() ;j++ ) {
	const pat::Electron *ele2 = thePatElectrons[j];
	if(ele2->charge()==ele1->charge()) continue;
	lep1.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
	lep2.SetPtEtaPhiE(ele2->pt(),ele2->eta(),ele2->phi(),ele2->energy());
	dilep=lep1+lep2;
	if (fabs(dilep.M()-91)<difference) {
	  difference=fabs(dilep.M()-91);
	  if(thePatElectrons.size()==2 && thePatMuons.size()==1){
	    const pat::Muon *mu1 = thePatMuons[0];
	    lep3.SetPtEtaPhiE(mu1->pt(),mu1->eta(),mu1->phi(),mu1->energy());
	    lep1saved=lep1;
            lep2saved=lep2;
            lep3saved=lep3;
	  }
	  else if(thePatElectrons.size()==2 && vTaus.size()==1){
            const reco::PFTau *tau = vTaus[0];
            lep3.SetPtEtaPhiE(tau->pt(),tau->eta(),tau->phi(),tau->energy());
            lep1saved=lep1;
            lep2saved=lep2;
            lep3saved=lep3;
          }

	  else {
	    int electron=3-i-j;
	    const pat::Electron *ele3 = thePatElectrons[electron];
	    lep3.SetPtEtaPhiE(ele3->pt(),ele3->eta(),ele3->phi(),ele3->energy());
            lep1saved=lep1;
            lep2saved=lep2;
            lep3saved=lep3;
	  }
	}
      }
    }
    for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
      const pat::Muon *mu1 = thePatMuons[i];
      for(unsigned int j = i+1 ; j < thePatMuons.size() ;j++ ) {
	const pat::Muon *mu2 = thePatMuons[j];
	if(mu2->charge()==mu1->charge()) continue;
	lep1.SetPtEtaPhiE(mu1->pt(),mu1->eta(),mu1->phi(),mu1->energy());
	lep2.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
	dilep=lep1+lep2;
	if (fabs(dilep.M()-91)<difference) {
	  difference=fabs(dilep.M()-91);
	  if(thePatMuons.size()==2 && thePatElectrons.size()==1 ){
	    const pat::Electron *ele1 = thePatElectrons[0];
	    lep3.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
            lep1saved=lep1;
            lep2saved=lep2;
            lep3saved=lep3;
	  }
          else if(thePatMuons.size()==2 && vTaus.size()==1){
            const reco::PFTau *tau = vTaus[0];
            lep3.SetPtEtaPhiE(tau->pt(),tau->eta(),tau->phi(),tau->energy());
            lep1saved=lep1;
            lep2saved=lep2;
            lep3saved=lep3;
          }
	  else{
	    int muon=3-i-j;
	    const pat::Muon *mu3 = thePatMuons[muon];
	    lep3.SetPtEtaPhiE(mu3->pt(),mu3->eta(),mu3->phi(),mu3->energy());
            lep1saved=lep1;
            lep2saved=lep2;
            lep3saved=lep3;

	  }
	}
      }
    }
  }
  else if(OSSF==2){
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
      const pat::Electron *ele1 = thePatElectrons[i];
      for(unsigned int j = 0 ; j < thePatMuons.size() ;j++ ) {
        const pat::Muon *mu2 = thePatMuons[j];
        if(mu2->charge()==ele1->charge()) continue;
        lep1.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
        lep2.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
        dilep=lep1+lep2;
        if(fabs(dilep.M()-50)<difference){
	  lep1saved=lep1;
	  lep2saved=lep2;
	  difference=(fabs(dilep.M()-50));
	  if(thePatElectrons.size()>1){
	    int electron=1-i;
	    const pat::Electron *ele2 = thePatElectrons[electron];
	    lep3.SetPtEtaPhiE(ele2->pt(),ele2->eta(),ele2->phi(),ele2->energy());
            lep3saved=lep3;
	  }
          if(thePatMuons.size()>1){
            int muon=1-j;
            const pat::Muon *mu2 = thePatMuons[muon];
            lep3.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
            lep3saved=lep3;
          }
          else if(vTaus.size()==1){
            const reco::PFTau *tau = vTaus[0];
            lep3.SetPtEtaPhiE(tau->pt(),tau->eta(),tau->phi(),tau->energy());
            lep3saved=lep3;
	    TLorentzVector dilep2;
	    if(tau->charge()==ele1->charge()){
	      dilep2=lep2+lep3;
	      if(fabs(dilep2.M()-60)<difference){
		lep1saved=lep2;
                lep2saved=lep3;
                lep3saved=lep1;
	      }
	    }
	    else {
	      dilep2=lep1+lep3;
              if(fabs(dilep2.M()-60)<difference){
                lep1saved=lep1;
                lep2saved=lep3;
                lep3saved=lep2;
	      }

	    }
          }

	}
      }
    }
  }


  else if(OSSF==3){
    const reco::PFTau *tau = vTaus[0];
    lep1.SetPtEtaPhiE(tau->pt(),tau->eta(),tau->phi(),tau->energy());
    lep1saved=lep1;
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
      const pat::Electron *ele1 = thePatElectrons[i];
      if(tau->charge()==ele1->charge()) continue;
      lep2.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
      dilep=lep1+lep2;
      if (fabs(dilep.M()-60)<difference) {
	difference=fabs(dilep.M()-60);
	if(thePatElectrons.size()==1 && thePatMuons.size()==1){
	  const pat::Muon *mu1 = thePatMuons[0];
	  lep3.SetPtEtaPhiE(mu1->pt(),mu1->eta(),mu1->phi(),mu1->energy());
	  lep1saved=lep1;
	  lep2saved=lep2;
	  lep3saved=lep3;
	}
 	else {
	  int electron=1-i;
	  const pat::Electron *ele3 = thePatElectrons[electron];
	  lep3.SetPtEtaPhiE(ele3->pt(),ele3->eta(),ele3->phi(),ele3->energy());
	  lep1saved=lep1;
	  lep2saved=lep2;
	  lep3saved=lep3;
	}
      }
    }
      for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
	const pat::Muon *mu1 = thePatMuons[i];
	if(tau->charge()==mu1->charge()) continue;
	lep2.SetPtEtaPhiE(mu1->pt(),mu1->eta(),mu1->phi(),mu1->energy());
	dilep=lep1+lep2;
	if (fabs(dilep.M()-60)<difference) {
	  difference=fabs(dilep.M()-60);
	  if(thePatMuons.size()==1 && thePatElectrons.size()==1){
	    const pat::Electron *ele1 = thePatElectrons[0];
	    lep3.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
	    lep1saved=lep1;
	    lep2saved=lep2;
	    lep3saved=lep3;
	  }
	  else {
	    int muon=1-i;
	    const pat::Muon *mu3 = thePatMuons[muon];
	    lep3.SetPtEtaPhiE(mu3->pt(),mu3->eta(),mu3->phi(),mu3->energy());
	    lep1saved=lep1;
	    lep2saved=lep2;
	    lep3saved=lep3;
	  }
	}
      }
    
  }

  SelectedLeptons.push_back(lep1saved);
  SelectedLeptons.push_back(lep2saved);
  SelectedLeptons.push_back(lep3saved);
  return SelectedLeptons;

}


  double tools::MT(TLorentzVector lep3, double pfMET_px, double pfMET_py, int OSSF){
      double pfMET_pt=sqrt(pfMET_px*pfMET_px+pfMET_py*pfMET_py);
      if(OSSF<1) return -1;
      return sqrt((pfMET_pt+lep3.Pt())*(pfMET_pt+lep3.Pt())-(pfMET_px+lep3.Px())*(pfMET_px+lep3.Px())-(pfMET_py+lep3.Py())*(pfMET_py+lep3.Py()));
  }

  int tools::DY(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons){
    TLorentzVector lep1, lep2, dilep, lep3;  // lepton 4-vectors                                                                                         
      for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
	  const pat::Electron *ele1 = thePatElectrons[i];
	  for(unsigned int j = i+1 ; j < thePatElectrons.size() ;j++ ) {
	    const pat::Electron *ele2 = thePatElectrons[j];
	    if(ele2->charge()==ele1->charge()) continue;
	    lep1.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
	    lep2.SetPtEtaPhiE(ele2->pt(),ele2->eta(),ele2->phi(),ele2->energy());
	    dilep=lep1+lep2;
	    if(dilep.M()<12) return 1;
	  }
      }
	
      for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
	  const pat::Muon *mu1 = thePatMuons[i];
	  for(unsigned int j = i+1 ; j < thePatMuons.size() ;j++ ) {
	    const pat::Muon *mu2 = thePatMuons[j];
	    if(mu2->charge()==mu1->charge()) continue;
	    lep1.SetPtEtaPhiE(mu1->pt(),mu1->eta(),mu1->phi(),mu1->energy());
	    lep2.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
	    dilep=lep1+lep2;
	    if(dilep.M()<12) return 1;
	  }
      }
      for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
	const pat::Electron *ele1 = thePatElectrons[i];
	for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
        const pat::Muon *mu1 = thePatMuons[i];
	if(mu1->charge()==ele1->charge()) continue;
	lep1.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
	lep2.SetPtEtaPhiE(mu1->pt(),mu1->eta(),mu1->phi(),mu1->energy());
	 dilep=lep1+lep2;
	 if(dilep.M()<12) return 1;
	}
      }

      return 0;
  }

double tools::MCTPerp1(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons, double pfMET_px, double pfMET_py, int OSSF){

  double mctperp=1000000;
  TVector3 MET(pfMET_px,pfMET_py,0);
  TVector3 UT;
  if(OSSF<1) return -1;
  if(OSSF==1){
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
      const pat::Electron *ele1 = thePatElectrons[i];
      for(unsigned int j = i+1 ; j < thePatElectrons.size() ;j++ ) {
	const pat::Electron *ele2 = thePatElectrons[j];
	if(ele2->charge()==ele1->charge()) continue;
	TVector3 lep1(ele1->px(),ele1->py(),0);
	TVector3 lep2(ele2->px(),ele2->py(),0);
	UT=-lep1-lep2-MET;
	TVector3 lep1_perp=lep1-UT.Dot(lep1)/(UT.Mag2())*UT;
        TVector3 lep2_perp=lep2-UT.Dot(lep2)/(UT.Mag2())*UT;
	double mctperptest=sqrt(2*(lep1_perp.Pt()*lep2_perp.Pt()+lep1_perp.Dot(lep2_perp)));
	if(mctperptest<mctperp) mctperp=mctperptest;
      }
    }
    for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
      const pat::Muon *mu1 = thePatMuons[i];
      for(unsigned int j = i+1 ; j < thePatMuons.size() ;j++ ) {
	const pat::Muon *mu2 = thePatMuons[j];
	if(mu2->charge()==mu1->charge()) continue;
        TVector3 lep1(mu1->px(),mu1->py(),0);
        TVector3 lep2(mu2->px(),mu2->py(),0);
        UT=-lep1-lep2-MET;
        TVector3 lep1_perp=lep1-UT.Dot(lep1)/(UT.Mag2())*UT;
        TVector3 lep2_perp=lep2-UT.Dot(lep2)/(UT.Mag2())*UT;
        double mctperptest=sqrt(2*(lep1_perp.Pt()*lep2_perp.Pt()+lep1_perp.Dot(lep2_perp)));
        if(mctperptest<mctperp) mctperp=mctperptest;
      }
    }
  }
  else if(OSSF==2){
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
      const pat::Electron *ele1 = thePatElectrons[i];
      for(unsigned int j = 0 ; j < thePatMuons.size() ;j++ ) {
	const pat::Muon *mu2 = thePatMuons[j];
	if(mu2->charge()==ele1->charge()) continue;
        TVector3 lep1(ele1->px(),ele1->py(),0);
        TVector3 lep2(mu2->px(),mu2->py(),0);
        UT=-lep1-lep2-MET;
       TVector3 lep1_perp=lep1-UT.Dot(lep1)/(UT.Mag2())*UT;
        TVector3 lep2_perp=lep2-UT.Dot(lep2)/(UT.Mag2())*UT;
        double mctperptest=sqrt(2*(lep1_perp.Pt()*lep2_perp.Pt()+lep1_perp.Dot(lep2_perp)));
        if(mctperptest<mctperp) mctperp=mctperptest;

      }
    }
  }
  return mctperp;
}


double tools::MCTPerp2(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons, double pfMET_px, double pfMET_py, int OSSF){

  double mctperp=-1;
  TVector3 MET(pfMET_px,pfMET_py,0);
  TVector3 UT;
  if(OSSF<1) return -1;
  if(OSSF==1){
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
      const pat::Electron *ele1 = thePatElectrons[i];
      for(unsigned int j = i+1 ; j < thePatElectrons.size() ;j++ ) {
        const pat::Electron *ele2 = thePatElectrons[j];
        if(ele2->charge()==ele1->charge()) continue;
        TVector3 lep1(ele1->px(),ele1->py(),0);
        TVector3 lep2(ele2->px(),ele2->py(),0);
        UT=-lep1-lep2-MET;
        TVector3 lep1_perp=lep1-UT.Dot(lep1)/(UT.Mag2())*UT;
        TVector3 lep2_perp=lep2-UT.Dot(lep2)/(UT.Mag2())*UT;
        double mctperptest=sqrt(2*(lep1_perp.Pt()*lep2_perp.Pt()+lep1_perp.Dot(lep2_perp)));
        if(mctperptest>mctperp) mctperp=mctperptest;
      }
    }
    for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
      const pat::Muon *mu1 = thePatMuons[i];
      for(unsigned int j = i+1 ; j < thePatMuons.size() ;j++ ) {
        const pat::Muon *mu2 = thePatMuons[j];
        if(mu2->charge()==mu1->charge()) continue;
        TVector3 lep1(mu1->px(),mu1->py(),0);
        TVector3 lep2(mu2->px(),mu2->py(),0);
        UT=-lep1-lep2-MET;
        TVector3 lep1_perp=lep1-UT.Dot(lep1)/(UT.Mag2())*UT;
        TVector3 lep2_perp=lep2-UT.Dot(lep2)/(UT.Mag2())*UT;
        double mctperptest=sqrt(2*(lep1_perp.Pt()*lep2_perp.Pt()+lep1_perp.Dot(lep2_perp)));
	if(mctperptest>mctperp) mctperp=mctperptest;
      }
    }
  }
  else if(OSSF==2){
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
      const pat::Electron *ele1 = thePatElectrons[i];
      for(unsigned int j = 0 ; j < thePatMuons.size() ;j++ ) {
        const pat::Muon *mu2 = thePatMuons[j];
        if(mu2->charge()==ele1->charge()) continue;
        TVector3 lep1(ele1->px(),ele1->py(),0);
        TVector3 lep2(mu2->px(),mu2->py(),0);
        UT=-lep1-lep2-MET;
        TVector3 lep1_perp=lep1-UT.Dot(lep1)/(UT.Mag2())*UT;
        TVector3 lep2_perp=lep2-UT.Dot(lep2)/(UT.Mag2())*UT;
        double mctperptest=sqrt(2*(lep1_perp.Pt()*lep2_perp.Pt()+lep1_perp.Dot(lep2_perp)));
        if(mctperptest>mctperp) mctperp=mctperptest;

      }
    }
  }
  return mctperp;
}



double tools::getZMassForDiEle(const std::vector<pat::Electron>  & thePatElectrons,                                                                                                                                 
     edm::Handle< std::vector<reco::Conversion> > &theConversions,                                                                                                                                     
     reco::BeamSpot::Point BS,                                                                                                                                                                       
     reco::Vertex::Point PV,  
     double  Rho,       
     const pat::Electron* myEle) {
    
    
  //    double Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13  };
    double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14 };

    double ZMass = 0;
    
    
    bool bool_pfIsolation = true;

    
    double v_electron_pt = 10.;
    double v_electron_eta = 2.4 ; 
    double v_electron_d0 = 0.04; 
    double v_electron_reliso = 0.2;
    double v_electron_dz = 0.2;
    bool bool_electron_ecalDriven = true;
    bool usePFiso =false;
    bool usePFisoCorr = false;
    
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {

    const reco::GsfTrackRef gsfTrack = el->gsfTrack();
    
    if( el->pt() < v_electron_pt ) continue;
    if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
  
    if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
  
    if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
    
    if(TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  continue;
    
    if(TMath::Abs(el->gsfTrack()->dz(PV)) > v_electron_dz  ) continue;
    
    

    
    if( TMath::Abs( el->eta() ) < 1.4442 )
    {
        if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.8 ) continue;
        if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
        if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
        if( TMath::Abs(el->hadronicOverEm()) > 0.15 ) continue;  //recommended is 0.12 but HLT applies 0.1
    }
    else if( TMath::Abs( el->eta() ) < 2.4 ) 
    {
        if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.7 ) continue;
        if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.01 ) continue;
        if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
    }
    
 
       
    double CorrectedTerm = 0.; 
    if( TMath::Abs( el->eta() ) < 1.0                                             )     CorrectedTerm = Rho * Aeff[ 0 ];
    else if( TMath::Abs( el->eta() ) > 1.0   && TMath::Abs( el->eta() ) < 1.479)     CorrectedTerm = Rho * Aeff[ 1 ];
    else if( TMath::Abs( el->eta() ) > 1.479 && TMath::Abs( el->eta() ) < 2.0  )     CorrectedTerm = Rho * Aeff[ 2 ];
    else if( TMath::Abs( el->eta() ) > 2.0   && TMath::Abs( el->eta() ) < 2.2  )     CorrectedTerm = Rho * Aeff[ 3 ];
    else if( TMath::Abs( el->eta() ) > 2.2   && TMath::Abs( el->eta() ) < 2.3  )     CorrectedTerm = Rho * Aeff[ 4 ];
    else if( TMath::Abs( el->eta() ) > 2.3   && TMath::Abs( el->eta() ) < 2.4  )     CorrectedTerm = Rho * Aeff[ 5 ];
    
    float pfRelIso = (el->chargedHadronIso() + TMath::Max(0.0, el->neutralHadronIso() + el->photonIso() - CorrectedTerm ) ) /el->pt() ;
    
    
    
    // Detector Isolation  
    float ecalIso = TMath::Abs(el->eta()) > 1.47 ? el->dr03TkSumPt() : TMath::Max(el->dr03TkSumPt()-1.,0.); 
    float detRelIso = (el->dr03EcalRecHitSumEt() + el->dr03HcalTowerSumEt() + ecalIso ) / el->pt() ;
    
    
    if( bool_pfIsolation ){
        if( pfRelIso > v_electron_reliso ) continue;
    }
    else{
        
        if( detRelIso > v_electron_reliso ) continue;
    }
    
    
    if( myEle->charge() == el->charge()  ) continue;
    
    TLorentzVector P1; P1.SetPtEtaPhiE( myEle->pt(), myEle->eta(), myEle->phi(), myEle->energy() );
    TLorentzVector P2; P2.SetPtEtaPhiE( el->pt(), el->eta(), el->phi(), el->energy() );
    
    float Mass = ( P1 + P2 ).M();
    
    
    if( Mass < 106. && Mass > 76. )   ZMass = Mass; 
 
    
    }

return ZMass;
}




double tools::getZMassForDiMu(const std::vector<pat::Muon>  & thePatMuons,   
                            reco::Vertex::Point PV,
                            const pat::Muon *myMu ){
    
    
    bool bool_pfIsolation = true;
    float v_muon_pt  = 10.;
    float v_muon_eta = 2.4;
    float v_muon_dz  = 999.;
    float v_muon_iso = 0.2;
  //  bool ZMass = false; 
    double ZMass = 0;

    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        if ( mu->isTrackerMuon() || mu->isGlobalMuon() ) {
        if ( !mu->isPFMuon() ) continue; 
        
    //    TVector3 momPerp(0,0,0);
    //    momPerp.SetPtEtaPhi(mu->pt(),mu->eta(),mu->phi());
    //    TVector3 posPerp(mu->vx()-PV.x(), mu->vy() - PV.y(), 0);
    //    float dzcorr = mu->vz() - PV.z() - posPerp.Dot(momPerp)/mu->pt() * (mu->pz()/mu->pt());
     //   if(TMath::Abs(dzcorr ) > v_muon_dz  ) continue;
        
        // PF Isolation  
        float chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
        float neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
        float photonIso = mu->pfIsolationR03().sumPhotonEt;
        float beta = mu->pfIsolationR03().sumPUPt;
        float pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
        
        // Det Isolation 
        float detRelIso = (mu->isolationR03().emEt + mu->isolationR03().hadEt + mu->isolationR03().sumPt) / mu->pt() ;
        
        if( bool_pfIsolation ){
        
	    if( pfRelIso > v_muon_iso ) continue;
        }
        else{
            
            if( detRelIso > v_muon_iso ) continue;
        }
        
        
        if( mu->charge() == myMu->charge() ) continue;  // We need opposite sign muons
        
        TLorentzVector P1; P1.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
        TLorentzVector P2; P2.SetPtEtaPhiE( myMu->pt(), myMu->eta(), myMu->phi(), myMu->energy() );
        
        float Mass = ( P1 + P2 ).M();
        
        if( Mass > 76. && Mass < 106. ) ZMass = Mass;
        
        
        
    }
}
    
return ZMass;
}





float tools::getPtUncertainty(float pt, float eta) {
	// returns the relative uncertainty on a jet pt, as a function of pt and eta
	// DISCLAIMER this function is only valid for the Same-Sign 2012 HPA analysis!!
	
	
	float result;
	const int nBinsEta = 2 ;
	const int nBinsPt = 10 ;
	
	float uncert[nBinsEta][nBinsPt] =
	{
		{ 12.0,  7.5,  6.3,  4.7,  3.7,  2.7,  2.6,  2.5,  2.4,  2.3 },
		{ 20.0,  17.2, 14.5,  10.7,  7.9,  6.9,  6.0,  5.2,  4.7,  4.4 }
	};
	
		
	int ipt = int(pt/10.);
	float rest = pt/10. - ipt;
	int ieta = 0;
	if (eta < -2.5 || eta > 2.5) {ieta = 1;}
		
	if (pt > 10. && pt < 100.) { 
		result = uncert[ieta][ipt-1] + rest*(uncert[ieta][ipt] - uncert[ieta][ipt-1]);
	} else if (pt >= 100.){
		result = uncert[ieta][9];
	} else {
		result = uncert[ieta][0];
	}
	
	
	return result*0.01;
}


double tools::JER (double eta) {

  double feta = fabs(eta);
  double scf = 1;
  if (feta < 0.5)
    scf = 1.052+0.062;
  else if (feta < 1.1)
    scf = 1.057+0.056;
  else if (feta < 1.7)
    scf = 1.096+0.063;
  else if (feta < 2.3)
    scf = 1.134+0.087;
  else scf = 1.288+0.155;
    
  return scf;
}













