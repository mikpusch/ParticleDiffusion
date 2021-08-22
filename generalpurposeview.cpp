GeneralPurposeView.cpp : implementation of the CGeneralPurposeView class
//

#include "stdafx.h"


void CGeneralPurposeView::OnViewParticlediffusion()
{
	// TODO: Add your command handler code here
	ViewType = 13;
	Invalidate();
}
int ScaleParticleX(double x, double SizeRectX, ParticleDiffusion & pd){
	return int (double(SizeRectX)*x/pd.SizeXYBoxInnm);
}

int ScaleParticleZ(double z, double SizeRectY, ParticleDiffusion & pd){
	return int (double(SizeRectY)*z/pd.HeightBoxInnm);
}

int ScaleParticleNProtX(int i, int Nsteps, double SizeRectX){
	return int (double(i)/double(Nsteps)*SizeRectX);
}
int ScaleParticleNProtY(double NProt, double SizeRectY){
	return int (SizeRectY- NProt*SizeRectY/10.0);
}

void CGeneralPurposeView::DrawParticleDiffusion(CDC* pDC){


	vector<int> indexVect;
	vector<double> NProtVect;
	vector<double> CurrentVect;
	vector<double> TransportVect;
	vector<double> PGFPVect;


	io myio;

	int NTransportOffOn = 1000000; // 50 ms
	if (myio.InInt(NTransportOffOn, "Numberof steps for each period (off/on)") != IDOK){
		return;
	}

	int NumberOnOff = 4;
	if (myio.InInt(NumberOnOff , "number of off/on events") != IDOK){
		return;
	}

	int NAv =             200000; // 10 us
	if (myio.InInt(NAv, "NAv") != IDOK){
		return;
	}

	if (myio.InDouble(PDiff.dtInnSec, "dtinnsec") != IDOK){
		return;
	}

	if (myio.InDouble(PDiff.SizeXYBoxInnm, "SizeXYBoxInnm") != IDOK){
		return;
	}
	if (myio.InDouble(PDiff.HeightBoxInnm, "HeightBoxInnm") != IDOK){
		return;
	}
	if (myio.InDouble(PDiff.TotalBufferInmM, "PDiff.TotalBufferInmM") != IDOK){
		return;
	}
	if (myio.InDouble(PDiff.FactorCaptureRadiusTransport, "FactorCaptureRadiusTransport") != IDOK){
		return;
	}

	if (myio.InDouble(PDiff.FactorCaptureRadiusGFP, "FactorCaptureRadiusGFP") != IDOK){
		return;
	}

	CString OutFile;
	CFileDialog diag( FALSE, NULL, OutFile, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, NULL, this);
	if (diag.DoModal() != IDOK){
	   return;
     }
    CString FileName = diag.GetPathName();

	CWaitCursor dummy;

	int NSteps = 2*NumberOnOff*NTransportOffOn;

	CRect Rect;
	GetClientRect(Rect);

	PDiff.TransportOn = false;
	PDiff.Init();
	int IAvg = 0;
	double NProtAvg = 0.0;

	double NTotalAvg = 0.0;
	int iTotalAvg = 0;

	int NGFPProtonated = 0;
	int NTransportOld = 0;

	for (int i=0; i<NSteps ; i++){
		int Episode = i/2/NTransportOffOn;
		int j = i - Episode*2*NTransportOffOn;

		if (j<NTransportOffOn){
			PDiff.TransportOn = false;
		}
		else{
			PDiff.TransportOn = true;
		}
		if (PDiff.GFP.protonated){
			NGFPProtonated++;
		}


		char s[20];

		if (i>NSteps/2){
			NTotalAvg += PDiff.Protons.size();
			iTotalAvg++;
		}

		NProtAvg += PDiff.Protons.size();
		IAvg++;

		if (IAvg>=NAv){
			NProtAvg /= double(NAv);
			double pGFP = double(NGFPProtonated)/double(NAv);

			int x = ScaleParticleNProtX( i, NSteps, Rect.Width());
			int y = ScaleParticleNProtY(NProtAvg, Rect.Height());
			if (i<NAv){
				pDC->MoveTo(x, y);
			}
			else{
				pDC->LineTo(x,y);
			}
			_itoa(i, s, 10);
			CString c = CString(s) + ": ";
			_gcvt(NProtAvg, 8, s);
			c += "NProts " + CString(s)+ CString(";  NTransport: ");

			_itoa(PDiff.NTransport, s, 10);
			c += CString(s)+ CString("; Current = ");

			double Current = double(PDiff.NTransport-NTransportOld)*1.6e-19/double(NAv)/(PDiff.dtInnSec*1e-9)/1e-15;
			NTransportOld = PDiff.NTransport;
			_gcvt(Current, 8, s);
			c += CString(s)+ CString(" fA; pGFP = ");
			
			_gcvt(pGFP, 8, s);
			c += CString(s);
			if (PDiff.TransportOn){
				c += "; Transport On               ";
			}
			else{
				c += "; Transport Off              ";
			}


			pDC->TextOutA(100, 100, c);

			indexVect.push_back(i);
			NProtVect.push_back(NProtAvg);
			CurrentVect.push_back(Current);
			TransportVect.push_back(PDiff.NTransport);
			PGFPVect.push_back(pGFP);


			IAvg = 0;
			NProtAvg = 0;
			NGFPProtonated = 0;

		}
		PDiff.Diffuse();
		PDiff.React();
	}


	NTotalAvg /= double(iTotalAvg);
	ShowFloat(NTotalAvg, "NTotalAvg");

	CString c;

	CWaitCursor dummy2;
	for (int i=0; i<indexVect.size(); i++){
		char s[20];
		_itoa(indexVect[i], s, 10);
		c+=CString(s)+"\t";
		_gcvt(NProtVect[i], 6, s);
		c+=CString(s)+"\t";
		_gcvt(CurrentVect[i], 6, s);
		c+=CString(s)+"\t";
		_gcvt(TransportVect[i], 6, s);
		c+=CString(s)+"\t";
		_gcvt(PGFPVect[i], 6, s);
		c+=CString(s)+"\r\n";
	}

    CFile file ( FileName, CFile::modeCreate | CFile::modeWrite);
	WriteStringOnFile(c, file);
	CopyTextToClipboard(c);
}
