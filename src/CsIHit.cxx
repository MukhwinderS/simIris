#include "CsIHit.h"

CsIHit::CsIHit()
{
	Thickness =12000.*4.51*0.1; //um*g/cm^3*0.1;
	mul = 0;
	fX.clear();
	fY.clear();
	fZ.clear();
	fPhiCalc.clear();
	fPhiRand.clear();
	//hit.clear();
	Seg.clear();
	dE.clear();
	dE_ideal.clear();
}

CsIHit::CsIHit(Double_t th)
{
	Thickness = th*4.51*0.1;
	mul = 0;
	fX.clear();
	fY.clear();
	fZ.clear();
	fPhiCalc.clear();
	fPhiRand.clear();
	//hit.clear();
	Seg.clear();
	dE.clear();
	dE_ideal.clear();
}

void CsIHit::Clear()
{
	mul = 0;
	fX.clear();
	fY.clear();
	fZ.clear();
	fPhiCalc.clear();
	fPhiRand.clear();
	//hit.clear();
	Seg.clear();
	dE.clear();
	dE_ideal.clear();
}

Double_t CsIHit::ThetaMin(Double_t CsIDistance)
{	
	const Double_t RIn = 50.;	// Inner radius in mm
	Double_t theta_min = TMath::ATan2(RIn,CsIDistance)*TMath::RadToDeg();
	return theta_min;
}

Double_t CsIHit::ThetaMax(Double_t CsIDistance)
{	
	const Double_t ROut = 150.;	// Outer radius in mm
	Double_t theta_max = TMath::ATan2(ROut,CsIDistance)*TMath::RadToDeg();
	return theta_max;
}

Bool_t CsIHit::Hit(Double_t theta, Double_t phi, Double_t CsIDistance, TVector3 targetPos)
{
	const Double_t RIn = 50.;	// Inner radius in mm
	const Double_t ROut = 150.;	// Outer radius in mm

	TRandom3 fRandom(0);
	Double_t fX0, fY0, fZ0;
	Double_t fPhiCalc0, fPhiRand0;
	Int_t Seg0;
	
	Bool_t hitTheta = 0;
	Bool_t hitPhi = 0;
	Bool_t hitBool = 0; // return value
	Double_t phiGap = 4.25*TMath::DegToRad(); //phi gap between CsI1s in rad 
	Double_t phiShift = -1.75*TMath::DegToRad();//shift from the vertical direction for the first CsI1
	Double_t phiRel;//Relative phi after phishift
	Double_t phiRange = TMath::Pi()/8-phiGap; //phi range for each CsI 

	phiRel = phi+phiShift;
	if (phiRel>TMath::Pi())  phiRel = phiRel - 2*TMath::Pi();
	
	// geometric efficiency
	fX0 = CsIDistance*tan(theta)*cos(phi);
	fY0 = CsIDistance*tan(theta)*sin(phi);
	
	TVector3 partVec(fX0,fY0,CsIDistance);
	
	partVec = partVec + targetPos; //taking into account the beam position at the target
	
	fX0 = partVec.X();
	fY0 = partVec.Y();
	fZ0 = partVec.Z();
	theta = partVec.Theta();
	phi = partVec.Phi();
	
	hitTheta = ((CsIDistance*tan(theta)>RIn) && (CsIDistance*tan(theta)<ROut));
	//hitPhi = kTRUE;
	//Seg0 = int((phiRel+TMath::Pi())/(TMath::Pi()/4));
  	Seg0 = int((TMath::Pi()-phiRel)*8./TMath::Pi());
	//hitPhi = fabs(phiRel +TMath::Pi() - Seg0*TMath::Pi()/4-TMath::Pi()/8) < phiRange/2;
	hitPhi = fabs(TMath::Pi() - phiRel - Seg0*TMath::Pi()/8.) < phiRange;
  	hitBool = (hitPhi && hitTheta); 
	
	if (hitBool){
		mul++;
		// hit0 = 1;
	  	fX0 = fX0;
		fY0 = fY0;	
 		fPhiCalc0 = 180.-((Seg0+0.5)*22.5+11.25+phiShift*TMath::RadToDeg());
 		if (fPhiCalc0<-180.) fPhiCalc0 = fPhiCalc0+360.;
      	Double_t random = 0.99*fRandom.Rndm();
 		fPhiRand0 = 180.-((Seg0+random)*22.5+11.25+phiShift*TMath::RadToDeg());
 		if (fPhiRand0<-180.) fPhiRand0 = fPhiRand0+360.;
		Seg0 = Seg0 + 12;
		if(Seg0>15) Seg0 = Seg0 - 16;
    	fX.push_back(fX0);
		fY.push_back(fY0);
		fZ.push_back(fZ0);
		fPhiCalc.push_back(fPhiCalc0);
		fPhiRand.push_back(fPhiRand0);
		//hit.push_back;
		Seg.push_back(Seg0);
	}
  	// else{
	// 	hit[mul] = 0;
	// 	fX[mul] = sqrt(-1);
	// 	fY[mul] = sqrt(-1);
	// 	fZ[mul] = sqrt(-1);
	// 	fPhiCalc[mul] = sqrt(-1);
	// 	Seg[mul] = sqrt(-1);
	// }
 	
	return hitBool;
}

Double_t CsIHit::ELoss(nucleus ncl, Double_t E, Double_t T)
{
  	Double_t dE0, dE_ideal0;
	//if(mul>0 && hit0==1){
		TRandom3 *rndm = new TRandom3(0);
		E -= eloss(ncl,15./31.,E,0.1*1.8219*0.1/cos(T),ncl.EL.eP, ncl.EL.dedxP);
		E -= eloss(ncl,13./27.,E,0.3*2.702*0.1/cos(T),ncl.EL.eAl, ncl.EL.dedxAl);
		E -= eloss(ncl,100./192.,E,6.*1.4*0.1/cos(T),ncl.EL.eMy, ncl.EL.dedxMy);
		dE0 = eloss(ncl,108./260.,E,Thickness/cos(T),ncl.EL.eCsI, ncl.EL.dedxCsI);
		E -= dE0;
		if(dE0<0.) dE0 = -dE0;
		dE_ideal0 = dE0;
		if(dE0!=0.) dE0 = rndm->Gaus(dE0,0.031*dE0*sqrt(14.1/dE0));
		if(dE0<0.) dE0 = 0.;
		dE.push_back(dE0);
		dE_ideal.push_back(dE_ideal0);
		//}
	return E;
}

void CsIHit::SortByEnergy()
{
	Bool_t have_swapped = true;
	while(have_swapped == true){
		for (Int_t x=0; x<mul; x++){
			have_swapped = false;
			for(Int_t y=0; y<mul-1; y++){
				if(dE[y]<dE[y+1]){
					std::swap(dE[y],dE[y+1]);
					std::swap(dE_ideal[y],dE_ideal[y+1]);
					std::swap(fX[y],fX[y+1]);
					std::swap(fY[y],fY[y+1]);
					std::swap(fZ[y],fZ[y+1]);
					std::swap(fPhiCalc[y],fPhiCalc[y+1]);
					std::swap(fPhiRand[y],fPhiRand[y+1]);
					std::swap(Seg[y],Seg[y+1]);
					have_swapped = true;
				}
			}
		}
	}
}
