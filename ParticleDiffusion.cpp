#include "stdafx.h"
#include "GeneralPurpose.h"

#include "ParticleDiffusion.h"


void BufferOrOtherMolecule::GetPosition(BufferOrOtherMolecule & x){
	for (int i=0; i<3; i++){
		position[i] = x.position[i];
	}
}

int BufferOrOtherMolecule::CountNotDeleted(vector<BufferOrOtherMolecule> & MolVect){
	int result = 0;
	for (int i=0; i<MolVect.size(); i++){
		if (!MolVect[i].deleted){
			result++;
		}
	}
	return result;
}


ParticleDiffusion::ParticleDiffusion(){
	Randomize();
	HeightBoxInnm = 100;
	SizeXYBoxInnm = 200;
	TotalBufferInmM = 20;
	pKBuffer = 7.0;
	kAssociationBuffer = 1.0e10;
	pKGFP = 7.0;
	DBuffer = .5e-10;
	DProt = 9.3e-9;
	DOH = 5.3e-9;

	BulkpH = 7.0;

	dtInnSec = 50;

	FactorCaptureRadiusProtonOH = .64;
	FactorCaptureRadiusBuff = 0.245;
	FactorCaptureRadiusGFP = 0.019;
	FactorCaptureRadiusTransport = 0.15;

	NA = 6.022e23;

	kDissWater = 1.0e-14*1.3e11;  // 1.3e11 is the associtaion rate H+ + OH-1


}
ParticleDiffusion::~ParticleDiffusion(){
}

double ParticleDiffusion::GetProtonationProbability(double pH, double pK){
	double r = pow(10.0, (pH-pK));
	return 1.0/(1.0+r);
}

int ParticleDiffusion::GetIntegerNumberWithProb(double N){
	int result = N;
	double diff = N-int(result);
	double p=Rand();
	if (p<diff){
		result++;
	}
	return result;
}

void ParticleDiffusion::GenerateRandomPosition(BufferOrOtherMolecule & p){
	p.position[0]=Rand()*SizeXYBoxInnm;
	p.position[1]=Rand()*SizeXYBoxInnm;
	p.position[2]=Rand()*HeightBoxInnm;
}

void ParticleDiffusion::Init(){

	NTransport = 0;
	Pore.position[0] = Pore.position[1] = 0.5*SizeXYBoxInnm;
	Pore.position[2] = 0;
	GFP.position[0] = Pore.position[0] + 0.5;
	GFP.position[1] = Pore.position[1] + 0.5;
	GFP.position[2] = .5;

	Vol = HeightBoxInnm*SizeXYBoxInnm*SizeXYBoxInnm*1e-24; // in liters
	NBuff = Vol*TotalBufferInmM*1e-3*NA;

	pGFP = GetProtonationProbability(BulkpH, pKGFP);

	double p=Rand();
	if (p<pGFP){
		GFP.protonated = true;

	}
	else{
		GFP.protonated = false;

	}


	pBuffer = GetProtonationProbability(BulkpH, pKBuffer);
	BufferMolecules.resize(0);
	BufferMolecules.reserve(NBuff);
	int NBuffProt = 0;
	int NBuffNotProt = 0;
	for (int i=0; i<NBuff; i++){
		BufferOrOtherMolecule bm;
		GenerateRandomPosition(bm);
		double p=Rand();
		if (p<pBuffer){
			bm.protonated = true;
			NBuffProt++;
		}
		else{
			bm.protonated = false;
			NBuffNotProt++;
		}
		BufferMolecules.push_back(bm);
	}


	double ConcProt = pow(10.0, -BulkpH);
	double ConcOH = 1e-14/ConcProt;
	NProt = GetIntegerNumberWithProb(Vol*NA*ConcProt);
	NOH = GetIntegerNumberWithProb(Vol*NA*ConcOH);

	Protons.resize(0);
	Protons.reserve(NProt);
	for (int i=0;i<NProt; i++){
		BufferOrOtherMolecule prot;
		GenerateRandomPosition(prot);
		Protons.push_back(prot);
	}
	OH.resize(0);
	OH.reserve(NOH);
	for (int i=0;i<NOH; i++){
		BufferOrOtherMolecule oh;
		GenerateRandomPosition(oh);
		OH.push_back(oh);
	}
	kProt = sqrt(DProt*3.0*dtInnSec*1e-9)/1e-9;
	kOH = sqrt(DOH*3.0*dtInnSec*1e-9)/1e-9;
	kBuffer = sqrt(DBuffer*3.0*dtInnSec*1e-9)/1e-9;

	NWater = Vol*NA*55.5;

	NDissWater = NWater*dtInnSec*kDissWater*1e-9;

	kDissBuffer = kAssociationBuffer*pow(10.0, -pKBuffer);

	NDissBuffer = double(BufferMolecules.size())*dtInnSec*1e-9*kDissBuffer; // per time step

	pDissGFP = 10.0*kAssociationBuffer*pow(10.0, -pKGFP)*dtInnSec*1e-9;
	pDissGFP = kAssociationBuffer*pow(10.0, -pKGFP)*dtInnSec*1e-9;

	CaptureHOHInnm = FactorCaptureRadiusProtonOH*0.5*(kProt+kOH);

	CaptureBufferProtOHInnm = FactorCaptureRadiusBuff*kBuffer;

	CaptureGFPInnm = FactorCaptureRadiusGFP*kProt;
	CaptureTransportInnm = FactorCaptureRadiusTransport*kProt;



}
bool ParticleDiffusion::CheckPosition(double pos[3]){
	for (int i=0; i<3; i++){
		if (pos[i]<0){
			return false;
		}
	}
	for (int i=0; i<2; i++){
		if (pos[i]>SizeXYBoxInnm){
			return false;
		}
	}
	if (pos[2]>HeightBoxInnm){
		return false;
	}
	return true;
}

void ParticleDiffusion::Diffuse(){
	for (int i=0; i<BufferMolecules.size(); i++){
		BufferOrOtherMolecule & bm = BufferMolecules[i];
		BufferOrOtherMolecule Old = bm;
		for (int j=0; j<3; j++){
			bm.position[j] += kBuffer*(Rand()-0.5);
		}
		if (!CheckPosition(bm.position)){
			bool NotMovedBeyondMembrane = (bm.position[2]>=0);
			bm = Old; 
			if (NotMovedBeyondMembrane){// replace by bulk buffer if z>0
				double p=Rand();
				if (p<pBuffer){
					bm.protonated = true;
				}
				else{
					bm.protonated = false;
				}
			}
		}
	}
	for (int i=0; i<Protons.size(); i++){
		BufferOrOtherMolecule & bm = Protons[i];
		BufferOrOtherMolecule Old = bm;
		for (int j=0; j<3; j++){
			bm.position[j] += kProt*(Rand()-0.5);
		}
		if (!CheckPosition(bm.position)){
			bm = Old;
		}
	}
	for (int i=0; i<OH.size(); i++){
		BufferOrOtherMolecule & bm = OH[i];
		BufferOrOtherMolecule Old = bm;
		for (int j=0; j<3; j++){
			bm.position[j] += kOH*(Rand()-0.5);
		}
		if (!CheckPosition(bm.position)){
			bm = Old;
		}
	}
}

bool ParticleDiffusion::WithinCapture(BufferOrOtherMolecule & p, BufferOrOtherMolecule & q, double SquareRadius){
	double r=0;
	for (int i=0; i<2; i++){
		double dd = p.position[i]-q.position[i];
		r += dd*dd;
		if (r>=SquareRadius){
			return false;
		}
	}
	return true;
}

void ParticleDiffusion::EliminateDeleted(vector<BufferOrOtherMolecule> & Mols){
	vector<BufferOrOtherMolecule> NewMols;
	NewMols.resize(0);
	for (int i=0; i<Mols.size(); i++){
		if (!Mols[i].deleted){
			NewMols.push_back(Mols[i]);
		}
	}
	Mols.resize(0);
	for (int i=0; i<NewMols.size(); i++){
		Mols.push_back(NewMols[i]);
	}
}

void ParticleDiffusion::React(){
	// test H+ transport
	// test GFP reaction
	// test recombination of H+ with OH-
	// test raction of H+ with negative buffer
	// test reacttion of OH with neutral buffer

	//test recombination of H+ with OH- 

#define DoShowParticle false

	for (int j=0; j<BufferMolecules.size(); j++){
		BufferMolecules[j].JustReacted = false;
	}
	if (DoShowParticle){
	 	ShowFloat(BufferOrOtherMolecule::CountNotDeleted(Protons), "NProt before");
 		ShowFloat(BufferOrOtherMolecule::CountNotDeleted(OH), "NOH before");
	}

	int NRealDiss = GetIntegerNumberWithProb(NDissWater);
	if (NRealDiss>0){
		for (int i=0; i<NRealDiss; i++){
			BufferOrOtherMolecule q;
			GenerateRandomPosition(q);
			Protons.push_back(q);
			OH.push_back(q);
		}
	}
	if (DoShowParticle){
		ShowFloat(BufferOrOtherMolecule::CountNotDeleted(Protons), "p water diss");
		ShowFloat(BufferOrOtherMolecule::CountNotDeleted(OH), "o after water diss");
	}

	// Dissociate buffer 
	int NRealBufferToProt = GetIntegerNumberWithProb(this->NDissBuffer/2.0);
	int NB = BufferMolecules.size();
	double NBdouble = double(NB);
	for (int i=0; i<NRealBufferToProt; i++){
		bool notyetdiss = true;
		while (notyetdiss){
			int testi = int (Rand() * NBdouble);
			if (testi<NB){
				BufferOrOtherMolecule & b = BufferMolecules[testi];
				if (b.protonated){
					b.protonated = false;
					b.JustReacted = true;
					BufferOrOtherMolecule p;
					p.GetPosition(b);
					p.deleted = false;
					Protons.push_back(p);
					notyetdiss=false;
				}
			}
		}
	}
	int NRealBufferToOH = GetIntegerNumberWithProb(this->NDissBuffer/2.0);
	for (int i=0; i<NRealBufferToOH; i++){
		bool notyetass = true;
		while (notyetass){
			int testi = int (Rand() * NBdouble);
			if (testi<NB){
				BufferOrOtherMolecule & b = BufferMolecules[testi];
				if (!b.JustReacted){
					if (!b.protonated){
						b.protonated = true;
						b.JustReacted = true;
						BufferOrOtherMolecule p;
						p.GetPosition(b);
						p.deleted = false;
						OH.push_back(p);
						notyetass=false;
					}
				}
			}
		}
	}

	if (DoShowParticle){
		ShowFloat(BufferOrOtherMolecule::CountNotDeleted(Protons), "P buf creation");
		ShowFloat(BufferOrOtherMolecule::CountNotDeleted(OH), "O buf creation");
	}

	bool protonsdeleted = false; 
	for (int j=0; j<OH.size(); j++){
		OH[j].deleted = false;
		OH[j].JustReacted = false;
	}
	for (int j=0; j<Protons.size(); j++){
		Protons[j].deleted = false;
		Protons[j].JustReacted = false;
	}

	double SquareRadius = CaptureHOHInnm*CaptureHOHInnm;
	for (int i=0; i<Protons.size(); i++){
		BufferOrOtherMolecule & p = Protons[i];
		for (int j=0; j<OH.size(); j++){
			BufferOrOtherMolecule & q = OH[j];
			if (!q.deleted){
				if (WithinCapture(p, q, SquareRadius)){
					p.deleted = true;
					q.deleted = true;
					protonsdeleted = true;
					break; // end for j loop
				}
			}
		}
	}
	if (DoShowParticle){
		ShowFloat(BufferOrOtherMolecule::CountNotDeleted(Protons), "P water ass");
		ShowFloat(BufferOrOtherMolecule::CountNotDeleted(OH), "O water ass");
	}
	//test recombination of H+ with Buffer- and OH- with buffer
	SquareRadius = CaptureBufferProtOHInnm*CaptureBufferProtOHInnm;
	for (int i=0; i<Protons.size(); i++){
		BufferOrOtherMolecule & p = Protons[i];
		if (!p.deleted){
			for (int j=0; j<BufferMolecules.size(); j++){
				bool reacted = false;
				BufferOrOtherMolecule & q = BufferMolecules[j];
				if (!q.protonated){
					if (!q.JustReacted){
						if (WithinCapture(p, q, SquareRadius)){
							q.protonated = true;
							q.JustReacted = true;
							p.deleted = true;
							reacted = true;
							protonsdeleted = true;
						}
					}
				}
				if (reacted){
					break; // for loop
				}
			}
		}
	}
	if (DoShowParticle){
		ShowFloat(BufferOrOtherMolecule::CountNotDeleted(Protons), "P H+buffer");
		ShowFloat(BufferOrOtherMolecule::CountNotDeleted(OH), "O h+buffer");
	}

	SquareRadius = CaptureBufferProtOHInnm*CaptureBufferProtOHInnm;
	for (int i=0; i<OH.size(); i++){
		BufferOrOtherMolecule & p = OH[i];
		if (!p.deleted){
			for (int j=0; j<BufferMolecules.size(); j++){
				bool reacted = false;
				BufferOrOtherMolecule & q = BufferMolecules[j];
				if (q.protonated){
					if (!q.JustReacted){
						if (WithinCapture(p, q, SquareRadius)){
							q.protonated = false;
							q.JustReacted = true;
							p.deleted = true;
							reacted = true;
							protonsdeleted = true;
						}
					}
				}
				if (reacted){
					break; // for loop
				}
			}
		}
	}

		
	if (DoShowParticle){
		ShowFloat(BufferOrOtherMolecule::CountNotDeleted(Protons), "P OH+buff");
		ShowFloat(BufferOrOtherMolecule::CountNotDeleted(OH), "O OH+buff");
	}


	// Check Transport
	if (TransportOn){
		SquareRadius = CaptureTransportInnm*CaptureTransportInnm;

		for (int i=0; i<Protons.size(); i++){
			BufferOrOtherMolecule & p = Protons[i];
			if (!p.deleted){
				if (WithinCapture(p, Pore, SquareRadius)){
					p.deleted = true;
					protonsdeleted = true;
					NTransport++;
					//ShowFloat(NTransport, "NT");
				}
			}
		}
	}

	// Check GFP
	if (GFP.protonated){
		double p=Rand();
		if (p<this->pDissGFP){
			GFP.protonated = false;
			BufferOrOtherMolecule  p;
			p.GetPosition(GFP);
			p.deleted = false;
			Protons.push_back(p);
			//Beep(500,100);
		}
	}
	else{
		SquareRadius = CaptureGFPInnm*CaptureGFPInnm;
		for (int i=0; i<Protons.size(); i++){
			BufferOrOtherMolecule & p = Protons[i];
			if (!p.deleted){
				if (WithinCapture(p, GFP, SquareRadius)){
					p.deleted = true;
					protonsdeleted = true;
					GFP.protonated = true;
					//Beep(3000,100);
					//ShowFloat(NTransport, "NT");
				}
			}
		}
	}


	if (protonsdeleted){
		EliminateDeleted(Protons);
		EliminateDeleted(OH);
	}


}


