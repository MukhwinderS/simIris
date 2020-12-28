// reacParams.h

#ifndef reacParams_H
#define reacParams_H
#include <TObject.h>
#include <TClass.h>
#include <string>

//Extern
//extern int gMesytecnitems;
class reacParams : public TObject {
	public:
		reacParams(); 
		virtual ~reacParams() {} //! 

		Int_t N;

		// Reaction A(a,b)B, 
		std::string A; 
		std::string a;
		std::string B;
		std::string b;
		std::string c;
		std::string d;
		std::string e;
		std::string f;
		std::string foil;

		Double_t E; // Beam energy
		Double_t R1; // Resonance energy
		Double_t R2; // Resonance energy
		Double_t W1; // Resonance width
		Double_t W2; // Resonance width
		Double_t SHAPE; // Resonance shape 0=Breit-Wigner, 1=Gaussian, 2=Square
		
		Bool_t SHT; // Reaction in solid Hydrogen/Deuterium target?
		
		
		//virtual void ReadCalibPar(char* line);
		virtual void ReadParams(char* line);
		virtual void Load(std::string filename);
		virtual void Print();
		virtual void Clear();
//		ClassDef(reacParams,1)
};

#endif
// end
