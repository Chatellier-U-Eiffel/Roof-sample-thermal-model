/*======================================================================================================================
* FILE : Define.h							VERSION : beta 1.0		DATE : 18 december 2009		INITIAL DEVELOPER : LCPC
* DESCRIPTION : Défintion des variables de contrôle de la conpilation
* A detailed list of autors and modifications is given the end of the file
*
* This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
* Public License. This software is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
* merchandability or fitness for a particular purpose.
======================================================================================================================*/

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <string>
#include <math.h>

using namespace std;

#include "Utils.hh"
#include "ThermalConduction.hh"
//======================================================================================================================

cDiscretisationToolForThermalConductionEquation::cDiscretisationToolForThermalConductionEquation()
:
	dxMin(0.01),
	Ratio(1.0),
	nLayersMin(4)
{}

void cDiscretisationToolForThermalConductionEquation::Compute(
	const double& dxLayer,
	unsigned int& nbSubLayers,
	double& FirstThicknessSubLayers) const
{
	double RefLayerWidth, WidthCheck;

	if (Ratio != 1.0) {cout << "cDiscretisationToolForThermalConductionEquation::Compute == Error 01" << endl; exit(1);}

	nbSubLayers = dxLayer / dxMin;
	FirstThicknessSubLayers = dxMin;
	return;

	if (nLayersMin % 2 == 1) RefLayerWidth = dxMin; else RefLayerWidth = 0.0;
	for (unsigned int j = 0 ; j < IntegerPart(nLayersMin / 2.0) ; j++)
		RefLayerWidth = 1.5 * RefLayerWidth + 2.0 * dxMin;
	if (dxLayer > RefLayerWidth) {
		unsigned int SubLayerID = 0;
		double TmpLayerWidth = 0.0;
		while (TmpLayerWidth < dxLayer) {
			TmpLayerWidth += 2.0 * dxMin * pow(Ratio, SubLayerID);
			SubLayerID++;
		}
		nbSubLayers = 2 * SubLayerID;
		FirstThicknessSubLayers = dxMin * dxLayer / TmpLayerWidth;
	} else {
		nbSubLayers = nLayersMin;
		FirstThicknessSubLayers = dxMin * dxLayer / RefLayerWidth;
	}
}

void cDiscretisationToolForThermalConductionEquation::FillWidth(
	const double dxMin,
	const double Ratio,
	vector<double>& dx) const
{
	unsigned int EndUp, EndDown;
	if (dx.size() % 2 == 1) {
		EndUp = IntegerPart(dx.size() / 2.0);
		EndDown = IntegerPart(dx.size() / 2.0) + 1;
	}
	else {
		EndUp = IntegerPart(dx.size() / 2.0);
		EndDown = IntegerPart(dx.size() / 2.0);
	}
	dx[0] = dxMin;
	dx[dx.size() - 1] = dxMin;
	for (unsigned int i = 1 ; i < EndUp ; i++) {
		dx[i] = dx[i-1] * Ratio;
		dx[dx.size()-1-i] = dx[i];
	}
	if (dx.size() % 2 == 1) dx[EndUp] = dx[EndUp-1] * Ratio;
}

//======================================================================================================================

/*! \class cThermalConductionEquation1D
Une description détaillée du solveur<br>
Possible sur plusierus lignes<br>
\f$i=j\f$
*/

cThermalConductionEquation1D::cThermalConductionEquation1D()
:
	Area(1.0),
	FrontBoundaryCondition(UndifinedBOundaryCondition),
	BackBoundaryCondition(UndifinedBOundaryCondition)
{
	SetStabilityFactor(0.0, 0.0, 1000.0);
}

/*! \fn cThermalConductionEquation1D::Initialize
Description détaillée de la fonction<br>
Sur plusieurs lignes
*/

void cThermalConductionEquation1D::Initialize(
	const vector<double>& _dxLayer,
	const vector<double>& _TCLayer,
	const vector<double>& _TCInterface,
	const vector<double>& _VHCLayer,
	const vector<double>& _TLayer)
{
	unsigned int LastSublayer;
	//vector<double> FirstThicknessSubLayers;
	//Log("cThermalConductionEquation1D::Initialize == Entree");
	if (_dxLayer.size() != _TCLayer.size()) FatalError("cThermalConductionEquation1D::Setdx == Initialize 01");
	else if (_dxLayer.size() != _VHCLayer.size()) FatalError("cThermalConductionEquation1D::Setdx == Initialize 02");

	dxOfLayers = _dxLayer;

	Initializedx2(dxOfLayers, dxOfSublayers, dx);

	TC.resize(GetNumberOfSublayers(), 0.0);
	TCInterface.resize(GetNumberOfSublayers(), 0.0);
	VHC.resize(GetNumberOfSublayers(), 0.0);
	T.resize(GetNumberOfSublayers(), 0.0);
	Q.resize(GetNumberOfSublayers(), 0.0);
	LastSublayer = 0;
	for (unsigned int i = 0 ; i < _dxLayer.size(); i++)  {
		//vector<double> Localdx(nbSubLayers[i], 0.0);
		//DiscretisationTool.FillWidth(FirstThicknessSubLayers[i], DiscretisationTool.Ratio, Localdx);
		for (unsigned int j = 0 ; j < dxOfSublayers[i].size() ; j++) {
			TC[j+LastSublayer] = _TCLayer[i];
			VHC[j+LastSublayer] = _VHCLayer[i];
			T[j+LastSublayer] = _TLayer[i];
		}
		if (_TCInterface[i] != 0.0) TCInterface[LastSublayer+dxOfSublayers[i].size()-1] = _TCInterface[i];
		LastSublayer += dxOfSublayers[i].size();
	}
	//Log("cThermalConductionEquation1D::Initialize == Sortie");
}

void cThermalConductionEquation1D::Initializedx1(const vector<double>& dxLayer, vector<vector<double> >& dxOfSublayers, vector<double>& dx) {
	const unsigned int DefaultNbOfSublayers = 5;
	unsigned int TotalNumberOfSublayers, SublayerID;

	//nbSubLayers2.resize(GetNumberOfLayers(), 5);
	dxOfSublayers.resize(GetNumberOfLayers());

	for (unsigned int i = 0 ; i < GetNumberOfLayers(); i++) {
		dxOfSublayers[i].resize(DefaultNbOfSublayers);
		for (unsigned int j = 0 ; j < dxOfSublayers[i].size() ; j++) dxOfSublayers[i][j] = dxLayer[i] / dxOfSublayers[i].size();
	}

	TotalNumberOfSublayers = 0;
	for (unsigned int i = 0 ; i < GetNumberOfLayers(); i++) TotalNumberOfSublayers += dxOfSublayers[i].size();
	dx.resize(TotalNumberOfSublayers);

	SublayerID = 0;
	for (unsigned int i = 0 ; i < GetNumberOfLayers(); i++)
		for (unsigned int j = 0 ; j < dxOfSublayers[i].size() ; j++) {
			dx[SublayerID] = dxOfSublayers[i][j];
			SublayerID++;
		}
}

void cThermalConductionEquation1D::Initializedx2(const vector<double>& dxLayer, vector<vector<double> >& dxOfSublayers, vector<double>& dx) {
	const unsigned int DefaultNbOfSublayers = 10;
	bool Refine = true;
	unsigned int TotalNumberOfSublayers, SublayerID;

	//Log("cThermalConductionEquation1D::Initializedx2 == Entree");
	dxOfSublayers.resize(GetNumberOfLayers());
	//nbSubLayers2.resize(GetNumberOfLayers(), 5);

	for (unsigned int i = 0 ; i < GetNumberOfLayers(); i++) {
		dxOfSublayers[i].resize(DefaultNbOfSublayers);
		for (unsigned int j = 0 ; j < dxOfSublayers[i].size() ; j++) dxOfSublayers[i][j] = dxLayer[i] / dxOfSublayers[i].size();
	}
	while (Refine == true) {
		Refine = false;
		for (unsigned int i = 0 ; i < GetNumberOfLayers(); i++) {
			if ((i != (GetNumberOfLayers()-1)) && (dxOfSublayers[i].back() > (4.0 * dxOfSublayers[i+1].front()))) {
				//Log(dxOfSublayers[i], "Avant");
				dxOfSublayers[i].back() /= 2.0;
				dxOfSublayers[i].push_back(dxOfSublayers[i].back());
				//Log(dxOfSublayers[i], "Après");
				Refine = true;
			}
			if ((i != 0) && (dxOfSublayers[i].front() > (4.0 * dxOfSublayers[i-1].back()))) {
				//Log(dxOfSublayers[i], "Avant");
				dxOfSublayers[i].front() /= 2.0;
				dxOfSublayers[i].insert(dxOfSublayers[i].begin(), dxOfSublayers[i].front());
				//Log(dxOfSublayers[i], "Après");
				Refine = true;
			}
		}
	}
	TotalNumberOfSublayers = 0;
	for (unsigned int i = 0 ; i < GetNumberOfLayers(); i++) TotalNumberOfSublayers += dxOfSublayers[i].size();
	dx.resize(TotalNumberOfSublayers);

	SublayerID = 0;
	for (unsigned int i = 0 ; i < GetNumberOfLayers(); i++)  {
		for (unsigned int j = 0 ; j < dxOfSublayers[i].size() ; j++) {
			//dxOfSublayers[i][j] = dxLayer[j] / nbSubLayers[i];
			dx[SublayerID] = dxOfSublayers[i][j];
			SublayerID++;
		}
		//LastSublayer += nbSubLayers[i];
	}
	//Log("cThermalConductionEquation1D::Initializedx2 == Sortie");
}

void cThermalConductionEquation1D::SetStabilityFactor(
	const double& _ImplicitFactor,
	const double& _RelaxationFactor,
	const double& _MaxDelta)
{
	ImplicitFactor = _ImplicitFactor;
	ExplicitFactor = 1.0 - _ImplicitFactor;
	RelaxationFactor = _RelaxationFactor;
	MaxDelta = _MaxDelta;
}

void cThermalConductionEquation1D::SetFrontBoundaryConditionDataSet(const BoundaryConditionType& _BCFrontType, const vector<double>& t, vector<double>& Value) {
	if (t[0] != 0.0) FatalError("void cThermalConductionEquation1D::SetFrontBoundaryConditionValueDataSet");
	BCFrontType = _BCFrontType;
	BCFrontt = t;
	BCFrontValue = Value;
	LastBCFrontt = 0;
}

void cThermalConductionEquation1D::SetFrontBoundaryCondition(BoundaryConditionType _FrontBoundaryCondition, const double Value) {
	FrontBoundaryCondition = _FrontBoundaryCondition;
	switch (FrontBoundaryCondition) {
		case ValueImposed : FrontValue = Value; break;
		case FluxImposed : FrontFlux = Value; break;
		default : cout << "cThermalConductionEquation1D::SetFrontBoundaryCondition == Error 01" << endl; exit(1);
	}
}

void cThermalConductionEquation1D::SetBackBoundaryConditionDataSet(const BoundaryConditionType& _BCBackType, const vector<double>& t, vector<double>& Value) {
	if (t[0] != 0.0) FatalError("void cThermalConductionEquation1D::SetBackBoundaryConditionValueDataSet");
	BCBackType = _BCBackType;
	BCBackt = t;
	BCBackValue = Value;
	LastBCBackt = 0;
}

void cThermalConductionEquation1D::SetBackBoundaryCondition(BoundaryConditionType _BackBoundaryCondition, const double Value) {
	BackBoundaryCondition = _BackBoundaryCondition;
	switch (BackBoundaryCondition) {
		case ValueImposed : BackValue = Value; break;
		case FluxImposed : BackFlux = Value; break;
		default : cout << "cThermalConductionEquation1D::SetBackBoundaryCondition == Error 01" << endl; exit(1);
	}
}

void cThermalConductionEquation1D::Solve(const double& dt) {SolveFiniteVolume(dt);}

unsigned int cThermalConductionEquation1D::GetNumberOfLayers() const {return dxOfLayers.size();}

unsigned int cThermalConductionEquation1D::GetNumberOfSublayers() const {return dx.size();}

double cThermalConductionEquation1D::GetBCFront(const double& t) {
	if (t < BCFrontt[LastBCFrontt]) FatalError("double cThermalConductionEquation1D::GetBCFront == Error 00");
	while (t > BCFrontt[LastBCFrontt+1]) {
		LastBCFrontt++;
		if (LastBCFrontt == BCFrontt.size()-1) FatalError("double cThermalConductionEquation1D::GetBCFront == Error 01");
	}
	return
		BCFrontValue[LastBCFrontt] +
		(BCFrontValue[LastBCFrontt+1] - BCFrontValue[LastBCFrontt]) / (BCFrontt[LastBCFrontt+1] - BCFrontt[LastBCFrontt]) *
		(t - BCFrontt[LastBCFrontt]);
}


void cThermalConductionEquation1D::UpdateFrontBoundaryCondition(const double& t) {
	SetFrontBoundaryCondition(BCFrontType, GetBCFront(t));}

void cThermalConductionEquation1D::UpdateBackBoundaryCondition(const double& t) {
	SetBackBoundaryCondition(BCBackType, GetBCBack(t));
}


double cThermalConductionEquation1D::GetBCBack(const double& t) {
	if (t < BCBackt[LastBCBackt]) FatalError("double cThermalConductionEquation1D::GetBCBack == Error 00");
	while (t > BCBackt[LastBCBackt+1]) {
		LastBCBackt++;
		if (LastBCBackt == BCBackt.size()-1) FatalError("double cThermalConductionEquation1D::GetBCBack == Error 01");
	}
	return
		BCBackValue[LastBCBackt] +
		(BCBackValue[LastBCBackt+1] - BCBackValue[LastBCBackt]) / (BCBackt[LastBCBackt+1] - BCBackt[LastBCBackt]) *
		(t - BCBackt[LastBCBackt]);
}

double cThermalConductionEquation1D::GetTFront() const {return T.front();}

double cThermalConductionEquation1D::GetTBack() const {return T.back();}

void cThermalConductionEquation1D::GetT(vector<double>& _T) const {
	unsigned int LastLayer = 0;
	if (dxOfLayers.size() != _T.size())
		FatalError("cThermalConductionEquation1D::GetT == Error 01 ("+tsi(dxOfLayers.size())+" vs "+tsi(_T.size())+")");

	for (unsigned int i = 0 ; i < _T.size(); i++)  {
		double MeanTemperature = 0.0;
		double dxOfLayer = 0;
		for (unsigned int j = 0 ; j < dxOfSublayers[i].size() ; j++) {
			MeanTemperature += T[LastLayer+j] * dx[LastLayer+j];
			dxOfLayer += dx[LastLayer+j];
		}
		MeanTemperature /= dxOfLayer;
		_T[i] = MeanTemperature;
		LastLayer += dxOfSublayers[i].size();
	}
}

double cThermalConductionEquation1D::GetEnergy() const {
	double Energy = 0.0;
	for (unsigned int i = 0 ; i < GetNumberOfSublayers() ; i++) Energy += VHC[i] * dx[i] * Area * T[i];
	return Energy;
}

void cThermalConductionEquation1D::SolveFiniteVolume(const double& dt) {
	const unsigned int nMax = GetNumberOfSublayers()-1;
	vector<double> a0(GetNumberOfSublayers(), 0.0), a1(GetNumberOfSublayers(), 0.0), a2(GetNumberOfSublayers(), 0.0);
	vector<double> c(GetNumberOfSublayers(), 0.0);
	vector<double> LX(GetNumberOfSublayers(), 0.0);
	//Log("cThermalConductionEquation1D::Solve == Entree");
	if (T.size() != GetNumberOfSublayers()) {
		cout << "cThermalConductionEquation1D::Solve == Error 01" << endl;
		cout << T.size() << " vs " << GetNumberOfSublayers() << endl;
		exit(1);
	} else if (FrontBoundaryCondition == UndifinedBOundaryCondition) {cout << "cThermalConductionEquation1D::Solve == Error 02" << endl; exit(1);}
	else if (BackBoundaryCondition == UndifinedBOundaryCondition) {cout << "cThermalConductionEquation1D::Solve == Error 03" << endl; exit(1);}

	for (unsigned int i = 0 ; i < (GetNumberOfSublayers() - 1) ; i++) {
		if (TCInterface[i] == 0.0) LX[i] = TC[i] * 2.0 / (dx[i+1] + dx[i]);
		else LX[i] = TCInterface[i];
	}
	switch (FrontBoundaryCondition) {
		case ValueImposed :
			a1[0] = 1.0;
			c[0] = FrontValue;
			break;
		case FluxImposed :
			a1[0] = Area * dx[0] * VHC[0] / dt + ImplicitFactor * Area * LX[0];
			a2[0] = -ImplicitFactor * Area * LX[0];
			c[0] = T[0] * (Area * dx[0] * VHC[0] / dt - ExplicitFactor * Area * LX[0]) +
				ExplicitFactor * Area * LX[0] * T[1] + FrontFlux * Area +
				Q[0];
			break;
	}
	switch (BackBoundaryCondition) {
		case ValueImposed :
			a1[nMax] = 1.0;
			c[nMax] = BackValue;
			break;
		case FluxImposed :
			a0[nMax] = -ImplicitFactor * Area * LX[nMax-1];
			a1[nMax] =  Area * dx[nMax] * VHC[nMax] / dt + ImplicitFactor * Area * LX[nMax-1];
			c[nMax] = T[nMax] * (Area * dx[nMax] * VHC[nMax] / dt - ExplicitFactor * Area * LX[nMax-1]) +
				ExplicitFactor * Area * LX[nMax-1] * T[nMax-1] + BackFlux * Area +
				Q[nMax];
			//Log("ExplicitFactor = " + ts(ExplicitFactor) + " / " + ts(BackFlux));
			break;
	}

	for (unsigned int i = 1 ; i < nMax ; i++) {
		a0[i] = - ImplicitFactor * Area * LX[i-1];
		a1[i] = Area * dx[i] * VHC[i] / dt + ImplicitFactor * Area * (LX[i] + LX[i-1]);
		a2[i] = - ImplicitFactor * Area * LX[i];
		c[i] = T[i] * (Area * dx[i] * VHC[i] / dt - ExplicitFactor * Area * (LX[i] + LX[i-1])) +
			ExplicitFactor * Area * (LX[i] * T[i+1] + LX[i-1] * T[i-1]) +
			Q[i];
	}
	//Log(a0, "a0_"); Log(a1, "a1_"); Log(a2, "a2_"); Log(c, "c_");//*/
	for (int i = c.size()-2 ; i >= 0 ; i--) {
		c[i] = c[i] - a2[i] * c[i+1] / a1[i+1];
		a1[i] = a1[i] - a2[i] * a0[i+1] / a1[i+1];
	}
	for (unsigned int i = 1 ; i < c.size() ; i++) c[i] = c[i] - a0[i] * c[i-1] / a1[i-1];
	for (unsigned int i = 0 ; i < c.size() ; i++) {
		const double NewT = c[i] / a1[i];
		double NewTRelax;
		NewTRelax = RelaxationFactor * T[i] + (1.0 - RelaxationFactor) * NewT;
		if (NewTRelax > T[i] + MaxDelta) T[i] +=  MaxDelta;
		else if (NewTRelax < T[i] - MaxDelta) T[i] -= MaxDelta;
		else T[i] = NewTRelax;
	}
	//GetT(Tout);
	//DisplayVector(T, "T");
	//Log("cThermalConductionEquation1D::Solve == Sortie");
}

void cThermalConductionEquation1D::SolveFiniteDifference(const double& dt, vector<double>& T) {
	const unsigned int NL1 = GetNumberOfSublayers()-1;
	vector<double> a0(GetNumberOfSublayers()+1, 0.0), a1(GetNumberOfSublayers()+1, 0.0), a2(GetNumberOfSublayers()+1, 0.0);
	vector<double> c(GetNumberOfSublayers()+1, 0.0);
	//Log("cThermalConductionEquation1D::Solve == Entree");
	if (T.size() !=(GetNumberOfSublayers()+1)) {
		cout << "cThermalConductionEquation1D::Solve == Error 01" << endl;
		cout << T.size() << " vs " << GetNumberOfSublayers()+1 << endl;
		exit(1);
	} else if (FrontBoundaryCondition == UndifinedBOundaryCondition) {cout << "cThermalConductionEquation1D::Solve == Error 02" << endl; exit(1);}
	else if (BackBoundaryCondition == UndifinedBOundaryCondition) {cout << "cThermalConductionEquation1D::Solve == Error 03" << endl; exit(1);}

	switch (FrontBoundaryCondition) {
		case ValueImposed :
			a1[0] = 1.0;
			c[0] = FrontValue;
			break;
		case FluxImposed :
			a1[0] = 1.0;
			a2[0] = -ImplicitFactor;
			c[0] = ExplicitFactor * T[1] + FrontFlux * dx[0] / TC[0];/*
			a1[0] = 1.0;
			a2[0] = 1.0;
			c[0] = T[0] + T[1] + 2.0 / dx[0] / VHC[0] * dt * (FrontFlux + 0.5 * TC[0] / dx[0] * (T[1] - T[0]));//*/
			break;
	}
	switch (BackBoundaryCondition) {
		case ValueImposed :
			a1[GetNumberOfSublayers()] = 1.0;
			c[GetNumberOfSublayers()] = BackValue;
			break;
		case FluxImposed :
			a1[GetNumberOfSublayers()] = 1.0;
			a0[GetNumberOfSublayers()] = -ImplicitFactor;
			c[GetNumberOfSublayers()] = ExplicitFactor * T[GetNumberOfSublayers()-1] +
				BackFlux * dx[GetNumberOfSublayers()-1] / TC[GetNumberOfSublayers()-1];/*
			a1[NumberOfLayers] = 1.0;
			a0[NumberOfLayers] = 1.0;
			c[NumberOfLayers] = T[NL1] + T[NumberOfLayers] + 2.0 / dx[NL1] / VHC[NL1] * dt * (BackFlux + 0.5 * TC[NL1] * (T[NL1] - T[NumberOfLayers]) / dx[NL1]);//*/
			//Log("ExplicitFactor = " + ts(ExplicitFactor) + " / " + ts(BackFlux));
			break;
	}
	for (unsigned int i = 1 ; i < GetNumberOfSublayers() ; i++) {
		const double Coeffi = 2.0 * TC[i] * ImplicitFactor / (dx[i] + dx[i-1]);
		const double Coeffe = 2.0 * TC[i] * ExplicitFactor / (dx[i] + dx[i-1]);
		const double LocalVHC = (VHC[i-1] + VHC[i]) / 2.0;
		a0[i] = - Coeffi / dx[i-1];
		a1[i] = LocalVHC / dt + Coeffi / dx[i] + Coeffi / dx[i-1];
		a2[i] = - Coeffi / dx[i];
		c[i] = LocalVHC / dt * T[i] + Coeffe * ((T[i+1] - T[i]) / dx[i] - (T[i] - T[i-1]) / dx[i-1]);
	}
	//if (FrontBoundaryCondition == ValueImposed) Log("Front ValueImposed " + ts(FrontValue)); else Log("Front FluxImposed " + ts(FrontFlux));
	//if (BackBoundaryCondition == ValueImposed) Log("Back ValueImposed " + ts(BackValue)); else Log("Back FluxImposed " + ts(BackFlux));
	//Log(TC, "TC_"); Log(a0, "a0_"); Log(a1, "a1_"); Log(a2, "a2_"); Log(c, "c_");//*/
	for (int i = c.size()-2 ; i >= 0 ; i--) {
		c[i] = c[i] - a2[i] * c[i+1] / a1[i+1];
		a1[i] = a1[i] - a2[i] * a0[i+1] / a1[i+1];
	}
	for (unsigned int i = 1 ; i < c.size() ; i++) c[i] = c[i] - a0[i] * c[i-1] / a1[i-1];
	for (unsigned int i = 0 ; i < c.size() ; i++) {
		const double NewT = c[i] / a1[i];
		double NewTRelax;
		NewTRelax = RelaxationFactor * T[i] + (1.0 - RelaxationFactor) * NewT;
		if (NewTRelax > T[i] + MaxDelta) T[i] +=  MaxDelta;
		else if (NewTRelax < T[i] - MaxDelta) T[i] -= MaxDelta;
		else T[i] = NewTRelax;
	}
	//DisplayVector(T, "T");
	//Log("cThermalConductionEquation1D::Solve == Sortie");
}

//======================================================================================================================

cThermalConductionEquation1D_Test::cThermalConductionEquation1D_Test()
:
	dxLayer(0.05),
	FilePrefix("Test_cThermalConductionEquation1D")
{}

void cThermalConductionEquation1D_Test::Test01(
	const unsigned int NumberOfLayers,
	const double ValueFront,
	const double ValueBack
) const
{
	ofstream File;
	File.open(FilePrefix+"_Sortie/Test01.csv");
	Test0x(10, ValueImposed, ValueFront, ValueImposed, ValueBack, File);
	File.close();
}

void cThermalConductionEquation1D_Test::Test02(
	const unsigned int NumberOfLayers,
	const double FluxFront, // J / m2
	const double FluxBack	// J / m2
) const
{
	ofstream File;
	File.open(FilePrefix+"_Sortie/Test02.csv");
	Test0x(10, FluxImposed, FluxFront, FluxImposed, FluxBack, File);
	File << "Variation d'Energie imposée = " << ts(FluxFront+FluxBack) << "J/s" << endl;

	File.close();
}

void cThermalConductionEquation1D_Test::Test03(
	const unsigned int NumberOfLayers,
	const double ValueFront,
	const double FluxBack
) const
{
	ofstream File;
	File.open(FilePrefix+"_Sortie/Test03.csv");
	Test0x(10, ValueImposed, ValueFront, FluxImposed, FluxBack, File);
	File.close();
}

void cThermalConductionEquation1D_Test::Test04(
	const unsigned int NumberOfLayers,
	const double FluxFront,
	const double ValueBack)
const
{
	ofstream File;
	File.open(FilePrefix+"_Sortie/Test04.csv");
	Test0x(10, FluxImposed, FluxFront, ValueImposed, ValueBack, File);
	File.close();
}


void cThermalConductionEquation1D_Test::Test0x(const unsigned int NumberOfLayers,
	const BoundaryConditionType BCFront, const double XFront,
	const BoundaryConditionType BCBack, const double XBack,
	ofstream& File) const
{
	const double dt = 1.0;
	vector<double> Thickness(NumberOfLayers, dxLayer);					// m
	vector<double> ThermalConductivity(NumberOfLayers, 1.0);		// J / m / K
	vector<double> InterfaceThermalConductivity(NumberOfLayers, 0.0);		// J / m / K
	vector<double> VolumetricHeatCapacity(NumberOfLayers, 2.0e6);	// J / m3 / K
	vector<double> Temperature(NumberOfLayers, 0.0);
	cThermalConductionEquation1D Equation;
	double Energy;

	Equation.Initialize(Thickness, ThermalConductivity, InterfaceThermalConductivity, VolumetricHeatCapacity, Temperature);
	Equation.SetStabilityFactor(0.0, 0.0, 100.0);
	Equation.SetFrontBoundaryCondition(BCFront, XFront);
	Equation.SetBackBoundaryCondition(BCBack, XBack);
	File << "t";
	for (unsigned int i = 0 ; i < Temperature.size() ; i++) File << sep << "T_"+to_string(i);
	File << endl;
	for (unsigned int t = 0 ; t < 30 * 1000; t++) {
		Equation.Solve(dt);
		Equation.GetT(Temperature);
		if ((t % 1000) == 0) {
			File << t * dt;
			for (unsigned int i = 0 ; i < Temperature.size() ; i++) File << sep << Temperature[i];
			File <<  endl;
		}
	}
	File << endl;
	Energy = - Equation.GetEnergy();
	Equation.Solve(dt);
	Equation.GetT(Temperature);
	Energy += Equation.GetEnergy();
	File << "Variation d'Energie calculee = " << ts(Energy / dt) << "J/s" << endl;
}

