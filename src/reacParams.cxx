#include "reacParams.h"
#include <stdio.h>
#include <stdlib.h>

reacParams::reacParams(){
  //reacParams::Class()->IgnoreTObjectStreamer;
  reacParams::Clear();
}

void reacParams::ReadParams(char* line)
{
	bool expect_val=true;
	char *from=line;
	char *to=line;
	while (*from) {
		if (*from>32) {*to=*from;to++;}
		from++;
	}
	*to=0;
	if (*line==0) return; // line is empty
	
	char* val=strchr(line,'=');
	if (!val){ 
		val=strchr(line, '!');
		expect_val=false;
	}
	if (!val) printf("Missing = or ! in input file, line: '%s'\n",line);
	*val=0;
	
	// trim param name
	char* trm=val-1;
	while (*trm<=32) *(trm--)=0;
	
	val++;
	if (*val==0 && expect_val) printf("Value missing for parameter %s",line);

	char cval[256];	
	double fval;	
	std::string strval;
	sscanf(val,"%lf",&fval);
	sscanf(val,"%s",cval);
	strval=cval;
	
	//	parameter of type string:
	if (strcmp(line,"N")==0){
	   N=int(fval);
	}
	if (strcmp(line,"A")==0){
	   A=strval;	
	}
	if (strcmp(line,"a")==0){
	   a=strval;
	}
	if (strcmp(line,"B")==0){
	   B=strval;
	}
	if (strcmp(line,"b")==0){
	   b=strval;
	}
	if (strcmp(line,"c")==0){
	   c=strval;
	}
	if (strcmp(line,"d")==0){
	   d=strval;
	}
	if (strcmp(line,"e")==0){
	   e=strval;
	}
	if (strcmp(line,"f")==0){
	   f=strval;
	}
	if (strcmp(line,"R1")==0){
	   R1=fval;
	}
	if (strcmp(line,"R2")==0){
	   R2=fval;
	}
	if (strcmp(line,"W1")==0){
	   W1=fval;
	}
	if (strcmp(line,"W2")==0){
	   W2=fval;
	}
	if (strcmp(line,"SHAPE")==0){
	   SHAPE=fval;
	}
	if (strcmp(line,"E")==0){
	   E=fval;
	}
	if (strcmp(line,"SHT")==0){
	   SHT=(fval==1.);
	}
	if (strcmp(line,"foil")==0){
      	   foil=strval;
    }

}

void reacParams::Load(std::string filename){	

	char line[256];
	FILE* file=fopen(filename.data(),"rb");
	if (!file)
	{
		printf("Cannot open reacParams file '%s' for reading. Stop.\n",filename.data());
		exit(0);
	}
	
	printf("Reading reacParams file '%s'\n",filename.data());
	
	while (!feof(file))
	{
		if (!fgets(line,256,file)) break;
		printf("%s",line);
		// skip leading white spaces
		char* ptr=line;
		while ((*ptr>0) && (*ptr<32)) ptr++;
		//printf("%s\n",ptr[0]);
		switch (ptr[0])
		{
			case 0   :
			case '#' :
			case '/' :  continue;
			default  :  ReadParams(ptr);
		}
	}
	fclose(file);
	file=NULL;
}

void reacParams::Print(){
	printf("\n************************************************\n");
	switch(N){
		case 2:
			printf("%s(%s,%s(%.1lf))%s(%.1lf) @ %.1lf MeV, Width=%.2lf\n",A.data(),a.data(),b.data(),R2,B.data(),R1,E,W1);
			break;
		case 3:
			printf("%s(%s,%s)%s(%.1lf)+%s @ %.1lf MeV, Width=%.2lf\n",A.data(),a.data(),b.data(),B.data(),R1,c.data(),E,W1);
			break;
			case 4:
			printf("%s(%s,%s)%s(%.1lf)+%s+%s @ %.1lf MeV, Width=%.2lf\n",A.data(),a.data(),b.data(),B.data(),R1,c.data(),d.data(),E,W1);
			break;
		case 5:
			printf("%s(%s,%s)%s(%.1lf)+%s+%s+%s @ %.1lf MeV, Width=%.2lf\n",A.data(),a.data(),b.data(),B.data(),R1,c.data(),d.data(),e.data(),E,W1);
			break;
		case 6:
			printf("%s(%s,%s)%s(%.1lf)+%s+%s+%s+%s @ %.1lf MeV, Width=%.2lf\n",A.data(),a.data(),b.data(),B.data(),R1,c.data(),d.data(),e.data(),f.data(),E,W1);
			break;
		default:
			printf("%s(%s,%s)%s(%.1lf)+%s+%s @ %.1lf MeV, Width=%.2lf\n",A.data(),a.data(),b.data(),B.data(),R1,c.data(),d.data(),E,W1);
			break;
	}
	if(SHT) printf("Reaction in solid Hydrogen/Deuterium Target.\n");
	else printf("Reaction in . %s foil.\n",foil.data());
	printf("\n************************************************\n\n");
}

void reacParams::Clear(){

	A.clear(); 
	a.clear();
	B.clear();
	b.clear();
	c.clear();
	d.clear();
	e.clear();
	f.clear();

	E=0.;
	R1=0.;
	R2=0.;
	W1=0.;
	W2=0.;
	SHAPE=0;

	N=0;
	SHT=0;
}
