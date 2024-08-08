
#ifndef TopVar_H_
#define TopVar_H_



#include<vector>
#include <array>
#include<cmath>
namespace TopAnalysis {

struct ParticleType {
	
	std::string type;
   	int 	PdgID;
   	float	Charge;
   	float	Mass;
   	float	Eta;
   	float	Phi;
   	TLorentzVector LVec;
   		
   	};
		   	
		
		
}
#endif
