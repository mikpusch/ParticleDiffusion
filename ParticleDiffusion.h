#pragma once

class BufferOrOtherMolecule{
public:
	bool protonated;
	bool deleted;
	bool JustReacted;
	double position[3];

	void GetPosition(BufferOrOtherMolecule & x);
	static int CountNotDeleted(vector<BufferOrOtherMolecule> & MolVect);
};



class ParticleDiffusion{
public:
	ParticleDiffusion();
	~ParticleDiffusion();

	double HeightBoxInnm;
	double SizeXYBoxInnm;
	double TotalBufferInmM;
	double pKBuffer;
	double kAssociationBuffer;
	double pKGFP;
	double DBuffer;
	double DProt;
	double DOH;

	double BulkpH;


	double dtInnSec;

	double FactorCaptureRadiusProtonOH;
	double FactorCaptureRadiusBuff;
	double FactorCaptureRadiusGFP;
	double FactorCaptureRadiusTransport;

	double CaptureBufferProtOHInnm;
	double CaptureHOHInnm;
	double CaptureGFPInnm;
	double CaptureTransportInnm;

	vector<BufferOrOtherMolecule> BufferMolecules;
    vector<BufferOrOtherMolecule> Protons;
    vector<BufferOrOtherMolecule> OH;
	BufferOrOtherMolecule Pore;
	BufferOrOtherMolecule GFP;

	double NA;
	int NBuff;
	int NProt;
	int NOH;
	double kProt, kOH, kBuffer;
	double pBuffer;
	double pGFP; // Bulk
	double kDissWater;
	double Vol;
	double NWater;
	double NDissWater; // per time step

	double kDissBuffer;
	double NDissBuffer; // per time step

	double pDissGFP;

	int NTransport;

	bool TransportOn;

	void Init();

	void Diffuse();
	void React();
	bool CheckPosition(double pos[3]);
	void GenerateRandomPosition(BufferOrOtherMolecule & p);
	static bool WithinCapture(BufferOrOtherMolecule & p, BufferOrOtherMolecule & q, double SquareRadius);
	static void EliminateDeleted(vector<BufferOrOtherMolecule> & Mols);

	static double GetProtonationProbability(double pH, double pK);
	static int GetIntegerNumberWithProb(double N);
};

