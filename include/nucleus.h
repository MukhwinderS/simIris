// nucleus.h
#ifndef nucleus_H
#define nucleus_H
#include <TObject.h>
#include <TClass.h>
#include <string>
#include "dedx.h"

class nucleus : public TObject {
	public:
		nucleus(); 
		virtual ~nucleus() {} //! 

		std::string name; 
		Int_t N;
		Int_t Z;
		Int_t A; 
		Double_t mass; 
		Double_t Sp; 
		Double_t Sn; 
		Double_t S2p; 
		Double_t S2n; 
		dedx EL; 

		//virtual void ReadCalibPar(char* line);
		virtual void getInfo(std::string,std::string);
		virtual void getInfo(std::string, Int_t, Int_t);
		virtual void ReadFile(std::string, Int_t, std::string);
		//virtual void SetFilenames(std::string);
		virtual void Print();
		virtual void Clear();
//		ClassDef(nucleus,1)
};

#endif
// end
