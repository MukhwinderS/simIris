 #include <stdlib.h>
#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TStopwatch.h"

#include "header.h"
#include "detHits.h"

reacParams reacPrm;
geoParams geoPrm;
YYHit yd, yu;
CsIHit csi;
S3Hit sd1, sd2, su;
PTrack blP, tlP;
PTrack blDecP, tlDecP1, tlDecP2;
PTrack buP1, buP2, buP3, buP4;
IDet det;

TGenPhaseSpace PS0, PS1;
TLorentzVector LorVb;
TLorentzVector LorVB;
TLorentzVector LorVBdec;
TLorentzVector LorVcdec;
TLorentzVector LorVddec;
TLorentzVector LorVcbu;
TLorentzVector LorVdbu;
TLorentzVector LorVebu;
TLorentzVector LorVfbu;

Double_t mA=0.;	
Double_t ma=0.;	
Double_t mB=0.;	
Double_t mBR=0.;	
Double_t mbR=0.;	
Double_t mb=0.;	
Double_t mc=0.;
Double_t md=0.;
Double_t me=0.;
Double_t mf=0.;

Double_t Qgen=sqrt(-1.);
Double_t Qdet=sqrt(-1.);
Double_t Qdet_sd=sqrt(-1.);
Double_t EB_det = sqrt(-1.);
Double_t EB_det_sd = sqrt(-1.);
Double_t PB_det = sqrt(-1.);
Double_t PB_det_sd = sqrt(-1.);
Double_t beamE=0.;
Double_t beamBeta=0.;
Double_t beamGamma=0.;
Double_t beamEcm=0.;
TVector3 reacPos;
Double_t convertMM=0;
Double_t adjustZpos=0;
Double_t SSBdE=0.;

void clearEvt()
{
	//TargdE[0]=0.; TargdE[1]=0.;
	Qgen=sqrt(-1.); Qdet=sqrt(-1.); Qdet_sd=sqrt(-1.);
	EB_det=sqrt(-1.); PB_det=sqrt(-1.);
	EB_det_sd=sqrt(-1.); PB_det_sd=sqrt(-1.);
	mBR=0.;
	mbR=0.;
	tlP.Clear();
	blP.Clear();
	blDecP.Clear();
	tlDecP1.Clear();
	tlDecP2.Clear();
	buP1.Clear();
	buP2.Clear();
	buP3.Clear();
	buP4.Clear();
	SSBdE=0.;
	yd.Clear();
	yu.Clear();
	csi.Clear();
	sd1.Clear();
	sd2.Clear();
	su.Clear();
	det.Clear();
	LorVb.Clear();
	LorVB.Clear();
	LorVBdec.Clear();
	LorVcdec.Clear();
	LorVddec.Clear();
	LorVcbu.Clear();
	LorVdbu.Clear();
	LorVebu.Clear();
	LorVfbu.Clear();
	PS0.Clear();
	PS1.Clear();
	
	return;
}

int main(int argc, char *argv[])
{
	Bool_t have_reaction=kFALSE;
	Bool_t have_geometry=kFALSE;
	Bool_t have_dwba_xsec=kFALSE;
	Bool_t be_silent=kFALSE;
	char *endptr;
	Int_t nsim = 1E3;
	Int_t Run = 0;
	char *reacParamsname =NULL;
	char *geoParamsname =NULL;
	char *dedxpath =NULL;
	char *outputname =NULL;
	char *dwbaname =NULL;
	Bool_t isSHTReac = kFALSE;

	std::string binpath(argv[0]);
	printf("%s\n",binpath.data());

	if (argc > 1){
		for(int i=0; i<argc; i++){
			if(strncmp(argv[i],"--output=",9)==0){
				outputname = argv[i]+9;
			}
			else if(strncmp(argv[i],"--dedx_dir=",11)==0){
				dedxpath = argv[i]+11;
			}
			else if(strncmp(argv[i],"--reaction=",11)==0){
				have_reaction=kTRUE;
				reacParamsname = argv[i]+11;
			}
			else if(strncmp(argv[i],"--geo=",6)==0){
				have_geometry=kTRUE;
				geoParamsname = argv[i]+6;
			}
			else if(strncmp(argv[i],"--dwba=",7)==0){
				have_dwba_xsec=kTRUE;
				dwbaname = argv[i]+7;
			}
			else if(strncmp(argv[i],"--events=",9)==0){
	  			nsim = strtol(argv[i]+9,&endptr,10);//converting string to number
			}
			else if(strncmp(argv[i],"--run=",6)==0){
	  			Run = strtol(argv[i]+6,&endptr,10);//converting string to number
			}
			else if(strncmp(argv[i],"--silent",8)==0){
				be_silent=kTRUE;
			}
		}
	}
	TStopwatch timer;
	timer.Start();
		
	if(have_reaction==kTRUE){
		printf("Loading parameters from file %s.\n",reacParamsname);
		reacPrm.Load(reacParamsname);
	}
	reacPrm.Print();
	
	isSHTReac=reacPrm.SHT;
	
	TRandom3 *rndm = new TRandom3(0);
	Int_t Evnt=0; 
	Double_t chck=0.;
	Double_t chck2=0.;
	Double_t wght=0.;
	Double_t wght2=0.;

	Double_t ICdE;	
	YYHit *ipyd = &yd; 
	YYHit *ipyu = &yu; 
	CsIHit *ipcsi = &csi; 
	S3Hit *ipsd1 = &sd1; 
	S3Hit *ipsd2 = &sd2; 
	S3Hit *ipsu = &su; 
	PTrack *iptlP = &tlP;
	PTrack *ipblP = &blP;
	PTrack *ipblDecP = &blDecP;
	PTrack *iptlDecP1 = &tlDecP1;
	PTrack *iptlDecP2 = &tlDecP2;
	PTrack *ipbuP1 = &buP1;
	PTrack *ipbuP2 = &buP2;
	PTrack *ipbuP3 = &buP3;
	PTrack *ipbuP4 = &buP4;
	IDet *ipdet = &det;
	TLorentzVector *LVb = &LorVb;	
	TLorentzVector *LVB = &LorVB;	
	TLorentzVector *LVBdec = &LorVBdec;	
	TLorentzVector *LVcdec = &LorVcdec;	
	TLorentzVector *LVddec = &LorVddec;	
	TLorentzVector *LVcbu = &LorVcbu;	
	TLorentzVector *LVdbu = &LorVdbu;	
	TLorentzVector *LVebu = &LorVebu;	
	TLorentzVector *LVfbu = &LorVfbu;	
	
	nucleus A, a, B, b, c, d, e, f, decB,decc,decd;
	Double_t reacX, reacY, reacZ;

    if (isSHTReac)
    {
        A.getInfo(binpath, reacPrm.A);
        A.Print();
        a.getInfo(binpath, reacPrm.a);
        a.Print();
        B.getInfo(binpath, reacPrm.B);
        B.Print();
        b.getInfo(binpath, reacPrm.b);
        b.Print();
    }
    else
    {
        A.getInfo(binpath, reacPrm.A);
        A.Print();
        a.getInfo(binpath, reacPrm.foil);
        a.Print();
        B.getInfo(binpath, reacPrm.B);
        B.Print();
        b.getInfo(binpath, reacPrm.foil);
        b.Print();
    }
    
	if(reacPrm.N>2){
		c.getInfo(binpath, reacPrm.c);
		c.Print();
	}
	if(reacPrm.N>3){
		d.getInfo(binpath, reacPrm.d);
		d.Print();
	}
	if(reacPrm.N>4){
		e.getInfo(binpath, reacPrm.e);
		e.Print();
	}
	if(reacPrm.N>5){
		f.getInfo(binpath, reacPrm.f);
		f.Print();
	}

    mA = A.mass/1000.;
	ma = a.mass/1000.;	
	mB = B.mass/1000.+reacPrm.R1/1000.;	
	mBR = mB;	
	mb = b.mass/1000.+reacPrm.R2/1000.;	
	mbR = mb;	
	mc = c.mass/1000.;
	md = d.mass/1000.;
	me = e.mass/1000.;
	mf = f.mass/1000.;
   // Printf("print A mass %lf and a mass %lf \n",A.mass,a.mass);
    
// Check for sequential decays ****************************************
	Bool_t seqdec=kFALSE;
	Double_t S_low=50.;
	Int_t seqdecN=0;
	Int_t pick=0;
	Double_t mBdec=0., mcdec=0., mddec=0.;
	Double_t masses2[3];

	if(B.Sn!=0.&&B.Sn<S_low){
		S_low=B.Sn;
		pick=1;
	}
	if(B.Sp!=0.&&B.Sp<S_low){
		S_low=B.Sp;
		pick=2;
	}
	if(B.S2n!=0.&&B.S2n<S_low){
		S_low=B.S2n;
		pick=3;
	}
	if(B.S2p!=0.&&B.S2p<S_low){
		S_low=B.S2p;
		pick=4;
	}
	printf("\nResonance Energy: %.2lf\tlowest threshold: %.2lf\n",reacPrm.R1,S_low);

	if(S_low<reacPrm.R1){
		switch(pick){
			case 1 : 
				seqdec = kTRUE;
				printf("\nSequential 1n- decay!\n\n"); 
				seqdecN =2;
				decB.getInfo(binpath, B.N-1,B.Z);
				decB.Print();
				mBdec=decB.mass/1000.;
				decc.getInfo(binpath, "n");
				decc.Print();
				mcdec=decc.mass/1000.;
				break;
			case 2 : 
				seqdec = kTRUE;
				printf("\nSequential 1p- decay!\n\n"); 
				seqdecN =2;
				decB.getInfo(binpath, B.N,B.Z-1);
				decB.Print();
				mBdec=decB.mass/1000.;
				decc.getInfo(binpath, "p");
				decc.Print();
				mcdec=decc.mass/1000.;
				break;
			case 3 : 
				seqdec = kTRUE;
				printf("\nSequential 2n- decay!\n\n"); 
				seqdecN =3;
				decB.getInfo(binpath, B.N-2,B.Z);
				decB.Print();
				mBdec=decB.mass/1000.;
				decc.getInfo(binpath, "n");
				decc.Print();
				mcdec=decc.mass/1000.;
				decd.getInfo(binpath, "n");
				decd.Print();
				mddec=decd.mass/1000.;
				break;
			case 4 : 
				seqdec = kTRUE;
				printf("\nSequential 2p- decay!\n\n"); 
				seqdecN =3;
				decB.getInfo(binpath, B.N,B.Z-2);
				decB.Print();
				mBdec=decB.mass/1000.;
				decc.getInfo(binpath, "p");
				decc.Print();
				mcdec=decc.mass/1000.;
				decd.getInfo(binpath, "p");
				decd.Print();
				mddec=decd.mass/1000.;
				break;
			default : 
				seqdec = kFALSE;
				break;
		}
		masses2[0] = mBdec;
		masses2[1] = mcdec;
		masses2[2] = mddec;
	}	
	
	
	if(have_geometry==kTRUE){
		printf("Loading parameters from file %s.\n",geoParamsname);
		geoPrm.Load(geoParamsname);
	}
	geoPrm.AoZTgt = a.A/a.Z;
	geoPrm.Print();


// *******************************************************************
	// Set up output file and tree
	TFile *file = new TFile(outputname,"RECREATE");
	TTree *Iris = new TTree("Iris","Iris simulation");
	
	Iris->Branch("Evnt",&Evnt,"Evnt/I"); 
	Iris->Branch("Run",&Run,"Run/I"); 
	Iris->Branch("beamE",&beamE,"beamE/D"); 
	Iris->Branch("beamBeta",&beamBeta,"beamBeta/D"); 
	Iris->Branch("beamGamma",&beamGamma,"beamGamma/D"); 
	Iris->Branch("beamEcm",&beamEcm,"beamEcm/D"); 
	Iris->Branch("reacX",&reacX,"reacX/D"); 
	Iris->Branch("reacY",&reacY,"reacY/D"); 
	Iris->Branch("reacZ",&reacZ,"reacZ/D"); 
	Iris->Branch("tlP.",&iptlP,32000,99); 
	Iris->Branch("blP.",&ipblP,32000,99); 
	Iris->Branch("blDecP.",&ipblDecP,32000,99); 
	Iris->Branch("tlDecP1.",&iptlDecP1,32000,99); 
	Iris->Branch("tlDecP2.",&iptlDecP2,32000,99); 
	Iris->Branch("buP1.",&ipbuP1,32000,99); 
	Iris->Branch("buP2.",&ipbuP2,32000,99); 
	Iris->Branch("buP3.",&ipbuP3,32000,99); 
	Iris->Branch("buP4.",&ipbuP4,32000,99); 
	Iris->Branch("wght",&wght,"wght/D"); 
	Iris->Branch("Qgen",&Qgen,"Qgen/D"); 
	Iris->Branch("Qdet",&Qdet,"Qdet/D"); 
	Iris->Branch("Qdet_sd",&Qdet_sd,"Qdet_sd/D"); 
	Iris->Branch("ICdE",&ICdE,"ICdE/D"); 
	Iris->Branch("SSBdE",&SSBdE,"SSBdE/D"); 
	Iris->Branch("yd.",&ipyd,32000,99); 
	Iris->Branch("yu.",&ipyu,32000,99); 
	Iris->Branch("csi.",&ipcsi,32000,99); 
	Iris->Branch("sd1.",&ipsd1,32000,99); 
	Iris->Branch("sd2.",&ipsd2,32000,99); 
	Iris->Branch("su.",&ipsu,32000,99); 
	Iris->Branch("det",&ipdet,32000,99); 

	std::string dedxstr = dedxpath;
	A.EL.loadIncomingELoss(dedxstr,A.name.data(),geoPrm.MFoil,geoPrm.MTgt,A.mass);
	b.EL.loadOutgoingELoss(dedxstr,b.name.data(),geoPrm.MFoil,geoPrm.MTgt,b.mass);
	if(!seqdec) B.EL.loadOutgoingELoss(dedxstr,B.name.data(),geoPrm.MFoil,geoPrm.MTgt,B.mass);
	else{
	   	decB.EL.loadOutgoingELoss(dedxstr,decB.name.data(),geoPrm.MFoil,geoPrm.MTgt,decB.mass);
	   	if(decc.Z>0) decc.EL.loadOutgoingELoss(dedxstr,decc.name.data(),geoPrm.MFoil,geoPrm.MTgt,decc.mass);
	   	if(seqdecN>2&&decd.Z>0) decd.EL.loadOutgoingELoss(dedxstr,decd.name.data(),geoPrm.MFoil,geoPrm.MTgt,decd.mass);
	}
	if(reacPrm.N>2&&c.Z>0) c.EL.loadOutgoingELoss(dedxstr,c.name.data(),geoPrm.MFoil,geoPrm.MTgt,c.mass);
	if(reacPrm.N>3&&d.Z>0) c.EL.loadOutgoingELoss(dedxstr,d.name.data(),geoPrm.MFoil,geoPrm.MTgt,d.mass);
	if(reacPrm.N>4&&e.Z>0) c.EL.loadOutgoingELoss(dedxstr,e.name.data(),geoPrm.MFoil,geoPrm.MTgt,e.mass);
	if(reacPrm.N>5&&f.Z>0) c.EL.loadOutgoingELoss(dedxstr,f.name.data(),geoPrm.MFoil,geoPrm.MTgt,f.mass);

	Double_t BeamSpot=geoPrm.Bs/2.355; // FWHM->sigma 
	Double_t ICLength=22.9*0.00318*geoPrm.ICPressure; //cm*mg/cm^3 
	const Double_t ICWindow1=0.03*3.44*0.1; //mu*g/cm^3*0.1
	const Double_t ICWindow2=0.05*3.44*0.1; //mu*g/cm^3*0.1
	
	//tlP.nuc = b;
	//blP.nuc = B;

	yd.Init(geoPrm.TYY);
	yu.Init(geoPrm.TYYU);
	sd1.Init(0,geoPrm.TS3[0]);
	sd2.Init(1,geoPrm.TS3[1]);
	su.Init(0,geoPrm.TS3U);

	Bool_t LEHit, HEHit;

	Int_t LEHitcntr=0;
	Int_t HEHitcntr=0;

	Double_t LEeff, HEeff;
	Double_t E_after_IC=0.;
	Double_t E_before_SSB=0.;
	Double_t E_before_Tgt=0.;
	Double_t E_center_Tgt=0.;
	Double_t E_after_Tgt=0.;
	Double_t E_before_Foil=0.;
	Double_t E_center_Foil=0.;
	Double_t E_after_Foil=0.;

	TLorentzVector target, beam, Sys;
	TVector3 boostvect;

	Double_t wght_max,width1,width2,sh;
	Bool_t allowed;

	// Calculate energy loss up to center of the target
	Double_t EA = reacPrm.E;	
   	EA -= eloss(A,0.5,EA,ICWindow1,A.EL.eSi3N4, A.EL.dedxSi3N4);
   	ICdE = eloss(A,0.586,EA,ICLength,A.EL.eC4H10, A.EL.dedxC4H10);
   	EA -= ICdE;
   	EA -= eloss(A,0.5,EA,ICWindow2,A.EL.eSi3N4, A.EL.dedxSi3N4);
	E_after_IC = EA;
   // printf("energy loss in IC full%lf \n Energy after IC%lf",ICdE,EA);
	if(!isSHTReac){
		if(geoPrm.Orientation==1){
			E_before_Tgt = EA;
   			EA -= eloss(A,1.,EA,geoPrm.TTgt/2.,A.EL.eTgt, A.EL.dedxTgt);
			E_center_Tgt = EA;
   			EA -= eloss(A,1.,EA,geoPrm.TTgt/2.,A.EL.eTgt, A.EL.dedxTgt);
			E_after_Tgt = EA;
		}
		
		E_before_Foil = EA;
		E_center_Foil = EA - eloss(A,1./geoPrm.AoZFoil,EA,geoPrm.TFoil/2.,A.EL.eFoil, A.EL.dedxFoil);
		EA -= eloss(A,1./geoPrm.AoZFoil,EA,geoPrm.TFoil,A.EL.eFoil, A.EL.dedxFoil);
   		E_before_Tgt = EA;
   		E_after_Foil = EA;
		if(geoPrm.Orientation==1) E_before_SSB = E_after_Foil;
		
		if(geoPrm.Orientation==0){
			E_before_Tgt = E_after_Foil;
   			E_center_Tgt = E_after_Foil - eloss(A,1.,E_after_Foil,geoPrm.TTgt/2.,A.EL.eTgt, A.EL.dedxTgt);
   			E_after_Tgt = E_after_Foil-eloss(A,1.,E_after_Foil,geoPrm.TTgt,A.EL.eTgt, A.EL.dedxTgt);
			E_before_SSB = E_after_Tgt;
		}
	}
	else{
		if(geoPrm.Orientation==0){
			E_before_Foil = EA;
			EA -= eloss(A,1./geoPrm.AoZFoil,EA,geoPrm.TFoil/2.,A.EL.eFoil, A.EL.dedxFoil);
   			E_center_Foil = EA;
			//EA -= eloss(A,1./geoPrm.AoZFoil,EA,geoPrm.TFoil/2.,A.EL.eFoil, A.EL.dedxFoil);
   			E_after_Foil = EA;
            printf("Energy loss in foil%lf\n Energy after the foil%lf \n ",2*E_center_Foil,EA);
		}
   		
		E_before_Tgt = EA;
   		E_center_Tgt = EA - eloss(A,1.,EA,geoPrm.TTgt/2.,A.EL.eTgt, A.EL.dedxTgt);
   		E_after_Tgt = EA-eloss(A,1.,EA,geoPrm.TTgt,A.EL.eTgt, A.EL.dedxTgt);
        //cout<<"\nEnergy centre of target"<<" "<<E_center_Tgt;
        //cout<<"\nEnergy before tgt"<<" "<<E_before_Tgt<<" "<<"energy after tgt"<<E_after_Tgt;
		if(geoPrm.Orientation==0) E_before_SSB = E_after_Tgt;
		
		reacZ = geoPrm.TTgt/2.;
        if (isSHTReac)
        {
            EA -= eloss(A,1.,EA,reacZ,A.EL.eTgt, A.EL.dedxTgt);
        }
        else
        {
            printf("reaction on foil with energy %lf",EA);
        }
	
		if(geoPrm.Orientation==1){
			E_before_Foil = E_after_Tgt;
			E_center_Foil = E_after_Tgt - eloss(A,1./geoPrm.AoZFoil,E_after_Tgt,geoPrm.TFoil/2.,A.EL.eFoil, A.EL.dedxFoil);
			E_after_Foil = E_after_Tgt - eloss(A,1./geoPrm.AoZFoil,E_after_Tgt,geoPrm.TFoil,A.EL.eFoil, A.EL.dedxFoil);
			E_before_SSB = E_after_Foil;
		}
	}
	
	EA = EA/1000.; // convert to GeV for TGenPhaseSpace
	Double_t PA = sqrt(EA*EA+2*EA*mA);
	target.SetXYZT(0.0, 0.0, 0.0, ma);
	beam.SetXYZT(0.0, 0.0, PA, mA+EA);
	Sys = beam + target;
	beamE = EA*1000.;
	beamBeta = Sys.Beta();
	beamGamma = Sys.Gamma();
	beamEcm = EA*ma*1000./(mA+ma);

	printf("\nEnergy after IC window: %.2lf MeV\n", E_after_IC);
	if(geoPrm.Orientation==1){
		printf("Energy before target: %.2lf MeV\n", E_before_Tgt);
		printf("Energy at center of target: %.2lf MeV\n", E_center_Tgt);
		printf("Energy at behind target: %.2lf MeV\n", E_after_Tgt);
	}
	printf("Energy before foil: %.2lf MeV\n", E_before_Foil);
	printf("Energy at center of foil: %.2lf MeV\n", E_center_Foil);
	printf("Energy after foil: %.2lf MeV\n", E_after_Foil);
	if(geoPrm.Orientation==0){
		printf("Energy before target: %.2lf MeV\n", E_before_Tgt);
		printf("Energy at center of target: %.2lf MeV\n", E_center_Tgt);
		printf("Energy at behind target: %.2lf MeV\n", E_after_Tgt);
	}

	printf("\nBeta at center of target: %.3lf \n", beamBeta);
	printf("Gamma at center of target: %.3lf \n", beamGamma);
	printf("CM Energy at center of target: %.2lf MeV\n\n", beamEcm);

	printf("YY1 detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",geoPrm.DYY,yd.ThetaMin(geoPrm.DYY),yd.ThetaMax(geoPrm.DYY)); 
	printf("CsI detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",geoPrm.DYY+11.6,csi.ThetaMin(geoPrm.DYY+11.6),csi.ThetaMax(geoPrm.DYY+11.6)); 
	printf("First S3 detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",geoPrm.DS3,sd1.ThetaMin(geoPrm.DS3),sd1.ThetaMax(geoPrm.DS3)); 
	printf("Second S3 detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",geoPrm.DS3+14.8,sd2.ThetaMin(geoPrm.DS3+14.8),sd2.ThetaMax(geoPrm.DS3+14.8)); 
	printf("Upstream YY1 detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",geoPrm.DYYU,yu.ThetaMin(geoPrm.DYYU),yu.ThetaMax(geoPrm.DYYU)); 
	printf("Upstream S3 detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n\n",geoPrm.DS3U,su.ThetaMin(geoPrm.DS3U),su.ThetaMax(geoPrm.DS3U)); 
	
	Double_t masses[6] = { mb, mB, mc, md, me, mf};
	
	Double_t tht=0.; 
	Double_t xsec=0.; 
	Double_t xsec_chck=0.; 
	Double_t xsec_max=0.; 
	Double_t dwba_th[181]={0.}; 
	Double_t dwba_xsec[181]={0.}; 
	
	if(have_dwba_xsec==kTRUE){ 
		printf("Using DWBA cross section from %s!\n",dwbaname);
		xsec_max = load_dwba(dwbaname,dwba_th,dwba_xsec); 
	}
	
	allowed = PS0.SetDecay(Sys, reacPrm.N, masses);
	
	if(!allowed){
		printf("Impossible decay!\n");
		printf("Exiting...\n");
		exit(0);
	}
	else{
		printf("Starting...\n");
	}

	wght_max=PS0.GetWtMax();
	printf("%lf\t%lf\n",wght_max,xsec_max);
	width1 = reacPrm.W1/1000.;
	width2 = reacPrm.W2/1000.;
	printf("%lf\t%lf\t%lf\t%lf\n",mB,width1,mb,width2);

	Int_t whilecount;
	// Start of event loop
	while(Evnt<nsim) 
	{
        if(!isSHTReac){
            //reacZ = geoPrm.TFoil/2.;
            reacZ = rndm->Uniform(0,geoPrm.TFoil);
               EA = E_before_Foil - eloss(A,1./geoPrm.AoZFoil,E_before_Foil,reacZ/2,A.EL.eFoil, A.EL.dedxFoil);
        }
        else{
            reacZ = rndm->Uniform(0,geoPrm.TTgt);
            //reacZ = geoPrm.TTgt/2.;
               EA = E_before_Tgt - eloss(A,b.Z/b.A,E_before_Tgt,reacZ,A.EL.eTgt, A.EL.dedxTgt);
         }
		EA = EA/1000.; // convert to GeV for TGenPhaseSpace
		PA = sqrt(EA*EA+2*EA*mA);

		target.SetXYZT(0.0, 0.0, 0.0, ma);
		beam.SetXYZT(0.0, 0.0, PA, mA+EA);
		Sys = beam + target;
		boostvect = Sys.BoostVector();

		beamE = EA*1000.;
		beamBeta = Sys.Beta();
		beamGamma = Sys.Gamma();
		beamEcm = EA*ma*1000./(mA+ma);

		//width = reacPrm.W/1000.;

		wght = 0.;
		clearEvt();
		if(reacPrm.SHAPE==0) mBR = rndm->BreitWigner(mB,width1);
		else if(reacPrm.SHAPE==1) mBR = rndm->Gaus(mB,width1);
		else if(reacPrm.SHAPE==2) mBR = rndm->Uniform(mB-width1,mB+width1);
		else mBR = rndm->BreitWigner(mB,width1);
        sh=reacPrm.SHAPE;
       // printf("value%lf",sh);
		masses[1] =mBR;
		if(reacPrm.SHAPE==0) mbR = rndm->BreitWigner(mb,width2);
		else if(reacPrm.SHAPE==1) mbR = rndm->Gaus(mb,width2);
        else if(reacPrm.SHAPE==2) mbR = rndm->Uniform(mb-width2,mb+width2);
		else mbR = rndm->BreitWigner(mb,width2);
		masses[0] =mbR;
		PS0.SetDecay(Sys, reacPrm.N, masses); //recalculate with resonance energy
		//wght_max=PS0.GetWtMax();
        sh=reacPrm.SHAPE;
      //  printf("value%lf\n",sh);
		TLorentzVector *LTmp;
		whilecount=0;
		do{	
			wght = PS0.Generate();
			if(wght!=wght) continue; // catch NaNs
			chck = rndm->Uniform(0,1);
			if(have_dwba_xsec==kTRUE){
				LTmp = PS0.GetDecay(0);
				LTmp->Boost(-boostvect);
				tht = TMath::RadToDeg()*LTmp->Theta();
				xsec=eval_theta(tht,dwba_th,dwba_xsec);
				xsec_chck = rndm->Uniform(0,xsec_max);
				LTmp->Boost(boostvect);
			}
			else{
				xsec=1.;
				xsec_chck=0.;
			}
			whilecount++;
			//printf("%d\t%f\t%f\t%f\n",whilecount,tht,xsec,xsec_chck);
		}while(wght<chck||xsec<xsec_chck);

		LVb  = PS0.GetDecay(0);
		LVB  = PS0.GetDecay(1);

		TLorentzVector Frag  = Sys-*LVb;

		Qgen= (mA+ma-mb-Frag.M())*1000.;	
	
		tlP.T=LVb->Theta();	
		blP.T=LVB->Theta();
		tlP.E=(LVb->E()-mb)*1000.; 	
		blP.E=(LVB->E()-mB)*1000.;
		tlP.P=LVb->Phi();	
		blP.P=LVB->Phi();	
		
		// Convert angles to degrees for root file
		tlP.Tdeg=TMath::RadToDeg()*tlP.T;
		blP.Tdeg=TMath::RadToDeg()*blP.T;
		tlP.Pdeg=TMath::RadToDeg()*tlP.P;
		blP.Pdeg=TMath::RadToDeg()*blP.P;

		if(seqdec)
		{
			PS1.SetDecay(*LVB, seqdecN, masses2);
			do{
				wght2 = PS1.Generate();			
				chck2 = rndm->Uniform(0,1);
				LVBdec  = PS1.GetDecay(0);
				LVcdec  = PS1.GetDecay(1);
				if(seqdecN>2) LVddec  = PS1.GetDecay(2);
			}while(wght2<chck2);
			blDecP.T=LVBdec->Theta();
			blDecP.E=(LVBdec->E()-mBdec)*1000.;
			blDecP.P=LVBdec->Phi();	
			blDecP.Tdeg=TMath::RadToDeg()*blDecP.T;
			blDecP.Pdeg=TMath::RadToDeg()*blDecP.P;
			tlDecP1.T=LVcdec->Theta();
			tlDecP1.E=(LVcdec->E()-mcdec)*1000.;
			tlDecP1.P=LVcdec->Phi();	
			tlDecP1.Tdeg=TMath::RadToDeg()*tlDecP1.T;
			tlDecP1.Pdeg=TMath::RadToDeg()*tlDecP1.P;
			if(seqdecN>2){
				tlDecP2.T=LVcdec->Theta();
				tlDecP2.E=(LVcdec->E()-mddec)*1000.;
				tlDecP2.P=LVcdec->Phi();	
				tlDecP2.Tdeg=TMath::RadToDeg()*tlDecP2.T;
				tlDecP2.Pdeg=TMath::RadToDeg()*tlDecP2.P;
			}
		}
		if(reacPrm.N>2) // 3body
		{
			LVcbu = PS0.GetDecay(2);
		
			buP1.T=LVcbu->Theta();	
			buP1.E=(LVcbu->E()-mc)*1000.; 	
			buP1.P=LVcbu->Phi();	
			
			// Convert angles to degrees for root file
			buP1.Tdeg=TMath::RadToDeg()*buP1.T;
			buP1.Pdeg=TMath::RadToDeg()*buP1.P;
		}
		if(reacPrm.N>3) // 4body
		{
			LVdbu = PS0.GetDecay(3);
		
			buP2.T=LVdbu->Theta();
			buP2.E=(LVdbu->E()-md)*1000.;
			buP2.P=LVdbu->Phi();	
			
			// Convert angles to degrees for root file
			buP2.Tdeg=TMath::RadToDeg()*buP2.T;
			buP2.Pdeg=TMath::RadToDeg()*buP2.P;
		}
		if(reacPrm.N>4) // 5body
		{
			LVdbu = PS0.GetDecay(4);
		
			buP3.T=LVebu->Theta();
			buP3.E=(LVebu->E()-me)*1000.;
			buP3.P=LVebu->Phi();	
			
			// Convert angles to degrees for root file
			buP3.Tdeg=TMath::RadToDeg()*buP3.T;
			buP3.Pdeg=TMath::RadToDeg()*buP3.P;
		}	
		if(reacPrm.N>5) // 6body
		{
			LVdbu = PS0.GetDecay(5);
		
			buP4.T=LVfbu->Theta();
			buP4.E=(LVfbu->E()-mf)*1000.;
			buP4.P=LVfbu->Phi();	
			
			// Convert angles to degrees for root file
			buP4.Tdeg=TMath::RadToDeg()*buP4.T;
			buP4.Pdeg=TMath::RadToDeg()*buP4.P;
		}		

		// Position on target	
		reacX = BeamSpot*rndm->Gaus();
		reacY = BeamSpot*rndm->Gaus();
		reacPos.SetXYZ(reacX,reacY,reacZ);
		
		tlP = TgtELoss(tlP, b, geoPrm, reacZ, isSHTReac);// Calculate the energy loss of the particle in Foil and SHT
		LEHit = detHits(tlP, b, reacPos,geoPrm.Mask,geoPrm.Shield);
		
		if(!seqdec){
			blP = TgtELoss(blP, B, geoPrm, reacZ, isSHTReac);// Calculate the energy loss of the particle in Foil and SHT
			HEHit = detHits(blP, B, reacPos,geoPrm.Mask,geoPrm.Shield);
		}
		else{ 
			blDecP = TgtELoss(blDecP, decB, geoPrm, reacZ, isSHTReac);// Calculate the energy loss of the particle in Foil and SHT
			HEHit = detHits(blDecP, decB, reacPos,geoPrm.Mask,geoPrm.Shield);	
			if(decc.Z>0){
				tlDecP1 = TgtELoss(tlDecP1, decc, geoPrm, reacZ, isSHTReac);// Calculate the energy loss of the particle in Foil and SHT
				detHits(tlDecP1, decc, reacPos,geoPrm.Mask,geoPrm.Shield);
			}	
			if(seqdecN>2&&decd.Z>0){
				tlDecP2 = TgtELoss(tlDecP2, decd, geoPrm, reacZ, isSHTReac);// Calculate the energy loss of the particle in Foil and SHT
				detHits(tlDecP2, decd, reacPos,geoPrm.Mask,geoPrm.Shield);
			}
		}
		if(reacPrm.N>2&&c.Z>0){
			buP1 = TgtELoss(buP1, c, geoPrm, reacZ, isSHTReac);// Calculate the energy loss of the particle in Foil and SHT
			detHits(buP1, c, reacPos,geoPrm.Mask,geoPrm.Shield);
		}	
		if(reacPrm.N>3&&d.Z>0){
			buP2 = TgtELoss(buP2, d, geoPrm, reacZ, isSHTReac);// Calculate the energy loss of the particle in Foil and SHT
			detHits(buP2, d, reacPos,geoPrm.Mask,geoPrm.Shield);
		}	
		if(reacPrm.N>4&&e.Z>0){
			buP3 = TgtELoss(buP3, e, geoPrm, reacZ, isSHTReac);// Calculate the energy loss of the particle in Foil and SHT
			detHits(buP3, e, reacPos,geoPrm.Mask,geoPrm.Shield);
		}
		if(reacPrm.N>5&&f.Z>0){
			buP4 = TgtELoss(buP4, f, geoPrm, reacZ, isSHTReac);// Calculate the energy loss of the particle in Foil and SHT
			detHits(buP4, f, reacPos,geoPrm.Mask,geoPrm.Shield);
		}
		// Calculate energy loss in SSB
		SSBdE = eloss(A,14./28.,E_before_SSB,500.*2.3212*0.1,B.EL.eSi,B.EL.dedxSi);
		SSBdE =rndm->Gaus(SSBdE,0.05*SSBdE);

		sortEnergies(); // sort detector hits by energy
	
		//Calculating "measured" Q-Value
		if(LEHit && yd.dE.size()>0 && csi.dE.size()>0){
			if(csi.dE[0]>0. && yd.dE[0]>0.){
			   	LEHitcntr++;
				Double_t Eb = csi.dE[0];
				Double_t cosTheta = TMath::Cos(yd.fThetaRand[0]*TMath::DegToRad());
				Eb= Eb+elossFi(Eb,0.1*1.4*6./cosTheta,b.EL.eMy,b.EL.dedxMy); //Mylar
	      		Eb= Eb+elossFi(Eb,0.1*2.702*0.3/cosTheta,b.EL.eAl,b.EL.dedxAl); //0.3 u Al
	      		Eb= Eb+elossFi(Eb,0.1*1.822*0.1/cosTheta,b.EL.eP,b.EL.dedxP); // 0.1Phosphorus
				Eb+= yd.dE[0]; //use measured Yd // change june28
	      		Eb= Eb+elossFi(Eb,0.1*2.32*0.35/cosTheta,b.EL.eSi,b.EL.dedxSi); //0.3 u Al + 1 um B equivalent in 0.35 um Si
	    		Eb= Eb+elossFi(Eb,geoPrm.TTgt/2./cosTheta,b.EL.eTgt,b.EL.dedxTgt); //deuteron energy  in mid target midtarget
		
				Eb= Eb/1000.;
				Double_t E_center = E_center_Tgt/1000.;
				Double_t Pb = sqrt(Eb*Eb+2.*Eb*mb);	
				Double_t P_center = sqrt(E_center*E_center+2*E_center*mA);
				//EB_det = EA+mA+ma-Eb-mb;
				EB_det = E_center+mA+ma-Eb-mb;
				PB_det = sqrt(P_center*P_center+Pb*Pb-2.*P_center*Pb*cosTheta);
				Qdet = mA+ma-mb-sqrt(EB_det*EB_det-PB_det*PB_det);
				Qdet =Qdet*1000.;
			}
		}
		
		//Calculating "measured" Q-Value
		if(sd1.dE.size()>0 && sd2.dE.size()>0){
			if(sd1.dE[0]>0. && sd2.dE[0]>0.){
				Double_t Eb = sd2.dE[0];
				Double_t cosTheta = TMath::Cos(sd1.fThetaCalc[0]*TMath::DegToRad());
				//Sd2 ring side
				Eb = Eb+elossFi(Eb,0.1*2.35*0.5/cosTheta,b.EL.eB,b.EL.dedxB); //boron junction implant
				Eb = Eb+elossFi(Eb,0.1*2.7*0.3/cosTheta,b.EL.eAl,b.EL.dedxAl); //first metal
				Eb = Eb+elossFi(Eb,0.1*2.65*2.5/cosTheta,b.EL.eSiO2,b.EL.dedxSiO2); //SiO2
				Eb = Eb+elossFi(Eb,0.1*2.7*1.5/cosTheta,b.EL.eAl,b.EL.dedxAl); //second metal
				//Sd1 ring side
				Eb = Eb+elossFi(Eb,0.1*2.7*1.5/cosTheta,b.EL.eAl,b.EL.dedxAl); //second metal
				Eb = Eb+elossFi(Eb,0.1*2.65*2.5/cosTheta,b.EL.eSiO2,b.EL.dedxSiO2); //SiO2
				Eb = Eb+elossFi(Eb,0.1*2.7*0.3/cosTheta,b.EL.eAl,b.EL.dedxAl); //first metal
				Eb+= sd1.dE[0]; //use measured Sd // change june28
				//Sd1 sector side
				Eb = Eb+elossFi(Eb,0.1*1.822*0.5/cosTheta,b.EL.eP,b.EL.dedxP); //phosphorus implant                                                                                   
				Eb = Eb+elossFi(Eb,0.1*2.7*0.3/cosTheta,b.EL.eAl,b.EL.dedxAl); //metal

				Eb= Eb+elossFi(Eb,geoPrm.TTgt/2./cosTheta,b.EL.eTgt,b.EL.dedxTgt); //deuteron energy  in mid target midtarget
		
				Eb= Eb/1000.;

				Double_t Pb = sqrt(Eb*Eb+2.*Eb*mb);	
				EB_det_sd = EA+mA+ma-Eb-mb;
				PB_det_sd = sqrt(PA*PA+Pb*Pb-2.*PA*Pb*cosTheta);
				Qdet_sd = mA+ma-mb-sqrt(EB_det_sd*EB_det_sd-PB_det_sd*PB_det_sd);
				Qdet_sd =Qdet_sd*1000.;
			}
		}

		if(HEHit && sd1.dE.size()>0 && sd2.dE.size()>0.){
			if(sd1.dE[0]>0. && sd2.dE[0]>0.) HEHitcntr++;
		}
		tlP.Ecm = (LVb->E()-mb)*ma*1000./(mA+ma);
		blP.Ecm = (LVB->E()-mB)*ma*1000./(mA+ma);
		LVb->Boost(-boostvect);
		LVB->Boost(-boostvect);
		tlP.Tcm = TMath::RadToDeg()*(TMath::Pi()-LVb->Theta());
		blP.Tcm = TMath::RadToDeg()*LVB->Theta();
		
		if(seqdec)
		{
			blDecP.Ecm = (LVBdec->E()-mBdec)*ma*1000./(mA+ma);
			LVBdec->Boost(-boostvect);
			blDecP.Tcm = TMath::RadToDeg()*LVBdec->Theta();
			tlDecP1.Ecm = (LVcdec->E()-mcdec)*ma*1000./(mA+ma);
			LVcdec->Boost(-boostvect);
			tlDecP1.Tcm = TMath::RadToDeg()*(TMath::Pi()-LVcdec->Theta());
			if(seqdecN>2){
				tlDecP2.Ecm = (LVddec->E()-mddec)*ma*1000./(mA+ma);
				LVddec->Boost(-boostvect);
				tlDecP2.Tcm = TMath::RadToDeg()*LVddec->Theta();
			}
		}
		else if(reacPrm.N>2) // 3body
		{
			buP1.Ecm = (LVcbu->E()-mc)*ma*1000./(mA+ma);
			LVcbu->Boost(-boostvect);
			buP1.Tcm = TMath::RadToDeg()*(TMath::Pi()-LVcbu->Theta());
		}
		else if(reacPrm.N>3) // 4body
		{
			buP2.Ecm = (LVdbu->E()-md)*ma*1000./(mA+ma);
			LVdbu->Boost(-boostvect);
			buP2.Tcm = TMath::RadToDeg()*LVdbu->Theta();
		}
		else if(reacPrm.N>4) // 3body
		{
			buP3.Ecm = (LVcbu->E()-me)*ma*1000./(mA+ma);
			LVebu->Boost(-boostvect);
			buP3.Tcm = TMath::RadToDeg()*(TMath::Pi()-LVebu->Theta());
		}
		else if(reacPrm.N>5) // 3body
		{
			buP4.Ecm = (LVcbu->E()-mf)*ma*1000./(mA+ma);
			LVfbu->Boost(-boostvect);
			buP4.Tcm = TMath::RadToDeg()*LVfbu->Theta();
		}
		setIDet(ICdE,SSBdE);
		if(!be_silent) printf("Writing %s: %.6d of %.6d events processed. Last event: %d tries.\r",outputname,Evnt,nsim,whilecount);
		Evnt++;
		Iris->Fill();
	}

	Iris->AutoSave();
	file->Close();
	rndm->Delete();
	//Iris->Delete();
	//f->Delete();
	LEeff=Double_t(LEHitcntr)/Double_t(nsim)*100.;
	HEeff=Double_t(HEHitcntr)/Double_t(nsim)*100.;
	printf("\n\n********************\n");
	printf("Acceptance for target-like particles: %.1f\n",LEeff);
	printf("Acceptance for beam-like particles: %.1f\n",HEeff);
	Double_t time=timer.RealTime();
	printf("\nDone. %lf s\n",time);
	printf("\nOutput written to %s \n", outputname);
	printf("\n\n********************\n");
	
	return 0;
}
