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

cThermalConductionEquationCylinder::cThermalConductionEquationCylinder () {}

void cThermalConductionEquationCylinder::Initialize(
			const vector<double>& _dRadius,
			const vector<double>& _dxLayer, //!< dx couche
			const vector<double>& TCLayer,
			const vector<double>& TCInterface,
			const vector<double>& VHCLayer,
			const vector<double>& TLayer)
{
	//double FirstRadius;
	Log("cThermalConductionEquationCylinder::Initialize == Entree");
	RadiusDelta = _dRadius;
	NumberOfRadius = RadiusDelta.size();
	Radius.resize(NumberOfRadius, 0.0);
	Annulus.resize(NumberOfRadius);

	Radius[0] = RadiusDelta[0];
	for (unsigned int j = 1 ; j < NumberOfRadius ; j++) Radius[j] = Radius[j-1] + RadiusDelta[j];
	for (unsigned int j = 0 ; j < NumberOfRadius ; j++) {
		Annulus[j].Initialize(_dxLayer, TCLayer, TCInterface, VHCLayer, TLayer);
		if (j == 0) Annulus[j].Area = M_PI * Radius[0] * Radius[0];
		else Annulus[j].Area = M_PI * (Radius[j] * Radius[j] - Radius[j-1] * Radius[j-1]);
	}
	for (unsigned int j = 1 ; j < NumberOfRadius ; j++) for (unsigned int i = 0 ; i < _dxLayer.size() ; i++)
		if (Annulus[j].dxOfSublayers[i].size() != Annulus[j-1].dxOfSublayers[i].size())
			FatalError("cThermalConductionEquationCylinder::Initialize == Error 01");
	Log("cThermalConductionEquationCylinder::Initialize == Sortie");
}

void cThermalConductionEquationCylinder::ReadFile(const string& FileName , const unsigned int col, vector<double>& Result) {
	ifstream File;
	string Line;

	File.open(FileName);
	Result.resize(0);
	while (getline(File, Line)) Result.push_back(Parse(Line, col));
	File.close();
}

void cThermalConductionEquationCylinder::Read(
	const string FileName,
	double& FinalTime, double& dt, double& Dt,
	vector<double>& _dRadius, vector<double>& _dxOfLayers,
	vector<double>& ThermalConductivity, vector<double>& InterfaceThermalConductivity,
	vector<double>& VolumetricHeatCapacity,
	vector<double>& InitialTemperature,
	BoundaryConditionType& FrontBoundaryConditionType, vector<double>& FrontBoundaryConditiont, vector<vector<double>>& FrontBoundaryConditionValue,
	BoundaryConditionType& BackBoundaryConditionType, vector<double>& BackBoundaryConditiont, vector<vector<double>>& BackBoundaryConditionValue,
	vector<double>& RadialBoundaryConditiont, vector<vector<double>>& RadialBoundaryConditionValueOfFlux
) {
	double BCType;
	string DataFileName;
	vector<double> Result;
	//Log("cThermalConductionEquationCylinder::Read == Entree");
	ParseXML(FileName, "TimeControl", Result);
	FinalTime = 3600.0 * Result[0];
	dt = Result[1];
	Dt = 3600.0 * Result[2];

	ParseXML(FileName, "Radius", Result);
	if (Result.size() != (Result[0]+1)) FatalError("cThermalConductionEquationCylinder::Read == 01 == Radius");
	_dRadius.resize(Result[0]);
	for (unsigned int i = 0 ; i < _dRadius.size() ; i++) _dRadius[i] = Result[i+1];
	//DisplayVector(_dRadius, "dRadius");

	ParseXML(FileName, "Height", Result);
	if (Result.size() != (Result[0]+1)) FatalError("cThermalConductionEquationCylinder::Read == 01 == Height");
	_dxOfLayers.resize(Result[0]);
	for (unsigned int i = 0 ; i < _dxOfLayers.size() ; i++) _dxOfLayers[i] = Result[i+1];
	//DisplayVector(_dxOfLayers, "dxOfLayers");

	ParseXML(FileName, "ThermalConductivity", ThermalConductivity);
	if (_dxOfLayers.size() != ThermalConductivity.size()) FatalError("cThermalConductionEquationCylinder::Read == 01 == ThermalConductivity");
	//DisplayVector(ThermalConductivity, "ThermalConductivity");

	ParseXML(FileName, "InterfaceThermalConductivity", InterfaceThermalConductivity);
	if (_dxOfLayers.size() != InterfaceThermalConductivity.size()) FatalError("cThermalConductionEquationCylinder::Read == 01 == InterfaceThermalConductivity");
	//DisplayVector(InterfaceThermalConductivity, "InterfaceThermalConductivity");

	ParseXML(FileName, "VolumetricHeatCapacity", VolumetricHeatCapacity);
	if (_dxOfLayers.size() != VolumetricHeatCapacity.size()) FatalError("cThermalConductionEquationCylinder::Read == 01 == VolumetricHeatCapacity");
	//DisplayVector(VolumetricHeatCapacity, "VolumetricHeatCapacity");

	ParseXML(FileName, "InitialTemperature", InitialTemperature);
	if (_dxOfLayers.size() != InitialTemperature.size()) FatalError("cThermalConductionEquationCylinder::Read == 01 == InitialTemperature");
	//DisplayVector(InitialTemperature, "InitialTemperature");

	ParseXML(FileName, "FrontBoundaryCondition", BCType, DataFileName);
	if (BCType == 1.0) FrontBoundaryConditionType = ValueImposed;
	else if (BCType == 2.0) FrontBoundaryConditionType = FluxImposed;
	else FatalError("cThermalConductionEquationCylinder::Read == 03");
	ReadFile(DataFileName, 0, FrontBoundaryConditiont);
	FrontBoundaryConditionValue.resize(_dRadius.size());
	for (unsigned int i = 0 ; i < _dRadius.size() ; i++) {
		ReadFile(DataFileName, i+1, Result);
		FrontBoundaryConditionValue[i] = Result;
	}



	ParseXML(FileName, "BackBoundaryCondition", BCType, DataFileName);
	if (BCType == 1.0) BackBoundaryConditionType = ValueImposed;
	else if (BCType == 2.0) BackBoundaryConditionType = FluxImposed;
	else FatalError("cThermalConductionEquationCylinder::Read == 03");
	ReadFile(DataFileName, 0, BackBoundaryConditiont);
	BackBoundaryConditionValue.resize(_dRadius.size());
	for (unsigned int i = 0 ; i < _dRadius.size() ; i++) {
		ReadFile(DataFileName, i+1, Result);
		BackBoundaryConditionValue[i] = Result;
	}

	ParseXML(FileName, "RadialBoundaryConditionValueOfFlux", BCType, DataFileName);
	ReadFile(DataFileName, 0, RadialBoundaryConditiont);
	RadialBoundaryConditionValueOfFlux.resize(_dxOfLayers.size());
	for (unsigned int i = 0 ; i < _dxOfLayers.size() ; i++) {
		ReadFile(DataFileName, i+1, Result);
		RadialBoundaryConditionValueOfFlux[i] = Result;
	}
}

void cThermalConductionEquationCylinder::SetStabilityFactor(const double& ImplicitFactor, const double& RelaxationFactor, const double& MaxDelta) {
	for (unsigned int j = 0 ; j < NumberOfRadius ; j++) Annulus[j].SetStabilityFactor(ImplicitFactor, RelaxationFactor, MaxDelta);
}

void cThermalConductionEquationCylinder::SetFrontBoundaryCondition(const BoundaryConditionType FrontBoundaryCondition, const double Value){
	for (unsigned int j = 0 ; j < NumberOfRadius ; j++) Annulus[j].SetFrontBoundaryCondition(FrontBoundaryCondition, Value);
}

void cThermalConductionEquationCylinder::SetFrontBoundaryCondition(const vector<BoundaryConditionType>& FrontBoundaryCondition, const vector<double>& Value) {
	for (unsigned int j = 0 ; j < NumberOfRadius ; j++) Annulus[j].SetFrontBoundaryCondition(FrontBoundaryCondition[j], Value[j]);
}

void cThermalConductionEquationCylinder::SetFrontBoundaryConditionDataSet(const BoundaryConditionType& BCFrontType, const vector<double>& t, vector<vector<double>>& Value){
	for (unsigned int j = 0 ; j < NumberOfRadius ; j++) Annulus[j].SetFrontBoundaryConditionDataSet(BCFrontType, t, Value[j]);
}


void cThermalConductionEquationCylinder::SetBackBoundaryCondition(const BoundaryConditionType BackBoundaryCondition, const double Value) {
	for (unsigned int j = 0 ; j < NumberOfRadius ; j++) Annulus[j].SetBackBoundaryCondition(BackBoundaryCondition, Value);
}

void cThermalConductionEquationCylinder::SetBackBoundaryConditionDataSet(const BoundaryConditionType& BCBackType, const vector<double>& t, vector<vector<double>>& Value){
	for (unsigned int j = 0 ; j < NumberOfRadius ; j++) Annulus[j].SetBackBoundaryConditionDataSet(BCBackType, t, Value[j]);
}

void cThermalConductionEquationCylinder::UpdateFrontBoundaryCondition(const double& t) {
	for (unsigned int j = 0 ; j < NumberOfRadius ; j++) Annulus[j].UpdateFrontBoundaryCondition(t);
}

void cThermalConductionEquationCylinder::UpdateBackBoundaryCondition(const double& t) {
	for (unsigned int j = 0 ; j < NumberOfRadius ; j++) Annulus[j].UpdateBackBoundaryCondition(t);
}

void cThermalConductionEquationCylinder::SetBackBoundaryCondition(const vector<BoundaryConditionType>& BackBoundaryCondition, const vector<double>& Value){
	for (unsigned int j = 0 ; j < NumberOfRadius ; j++) Annulus[j].SetBackBoundaryCondition(BackBoundaryCondition[j], Value[j]);
}

void cThermalConductionEquationCylinder::SetRadialFluxBoundaryConditionDataSet(const vector<double>& t, const vector<vector<double>>& Value) {	if (t[0] != 0.0) FatalError("void cThermalConductionEquationCylinder::SetRadialFluxBoundaryConditionDataSet");
	BCRadialt = t;
	BCRadialValue = Value;
	LastBCRadialt = 0;
}

void cThermalConductionEquationCylinder::SetRadialFluxBoundaryCondition(const vector<double>& RadiusFluxLayer) {
	unsigned int LastLayer;
	if (RadiusFluxLayer.size() != Annulus[NumberOfRadius-1].GetNumberOfLayers()) {
		FatalError("cThermalConductionEquationCylinder::SetBackBoundaryConditionFluxR == Error 01");
		exit(1);
	}
	Qbc.resize(Annulus[NumberOfRadius-1].GetNumberOfSublayers(), 0.0);
	LastLayer = 0;
	for (unsigned int i = 0 ; i < Annulus[NumberOfRadius-1].GetNumberOfLayers(); i++)  {
		for (unsigned int j = 0 ; j < Annulus[NumberOfRadius-1].dxOfSublayers[i].size() ; j++)
			Qbc[j+LastLayer] = RadiusFluxLayer[i] *
				2.0 * M_PI * Radius[NumberOfRadius-1] * Annulus[NumberOfRadius-1].dx[j+LastLayer];
		LastLayer += Annulus[NumberOfRadius-1].dxOfSublayers[i].size();
	}
}

void cThermalConductionEquationCylinder::UpdateRadialBoundaryCondition(const double& t) {
	vector<double> RadiusFluxLayer(Annulus[NumberOfRadius-1].GetNumberOfLayers());
	//Log("cThermalConductionEquationCylinder::UpdateRadialBoundaryCondition == Entree");
	if (t < BCRadialt[LastBCRadialt]) FatalError("cThermalConductionEquationCylinder::UpdateFrontBoundaryCondition == Error 00");

	while (t > BCRadialt[LastBCRadialt+1]) {
		LastBCRadialt++;
		if (LastBCRadialt == BCRadialt.size()-1) FatalError("cThermalConductionEquationCylinder::UpdateFrontBoundaryCondition == Error 01");
	}
	for (unsigned int i = 0 ; i < RadiusFluxLayer.size() ; i++) {
		RadiusFluxLayer[i] = BCRadialValue[i][LastBCRadialt] +
			(BCRadialValue[i][LastBCRadialt+1] - BCRadialValue[i][LastBCRadialt]) / (BCRadialt[LastBCRadialt+1] - BCRadialt[LastBCRadialt]) *
			(t - BCRadialt[LastBCRadialt]);
	}
	SetRadialFluxBoundaryCondition(RadiusFluxLayer);
	//Log("cThermalConductionEquationCylinder::UpdateRadialBoundaryCondition == Sortie");
}

void cThermalConductionEquationCylinder::Solve(const double& dt) {
	//Log("cThermalConductionEquationCylinder::Solve == Entree");
	for (unsigned int j = 0 ; j < NumberOfRadius ; j++) for (unsigned int i = 0 ; i < Annulus[j].dx.size() ; i++) {
		double FluxInterior, FluxExterior;
		if (j == 0)  {
			const double TCRj = (Annulus[j+1].TC[i] * RadiusDelta[j+1] + Annulus[j].TC[i] * RadiusDelta[j]) / (RadiusDelta[j+1] + RadiusDelta[j]);
			FluxInterior = 0.0;
			FluxExterior = Annulus[j].dx[i] * 2.0 * M_PI * Radius[j]   * TCRj  * (Annulus[j+1].T[i] - Annulus[j].T[i]) / ((RadiusDelta[j+1] + RadiusDelta[j]) / 2.0);
		} else if (j == (NumberOfRadius-1)) {
			const double TCRj1 = (Annulus[j].TC[i] * RadiusDelta[j] + Annulus[j-1].TC[i] * RadiusDelta[j-1]) / (RadiusDelta[j] + RadiusDelta[j-1]);
			FluxInterior = -Annulus[j].dx[i] * 2.0 * M_PI * Radius[j-1] * TCRj1 * (Annulus[j].T[i] - Annulus[j-1].T[i]) / ((RadiusDelta[j] + RadiusDelta[j-1]) / 2.0);;
			FluxExterior = Qbc[i];
		}
		else {
			const double TCRj = (Annulus[j+1].TC[i] * RadiusDelta[j+1] + Annulus[j].TC[i] * RadiusDelta[j]) / (RadiusDelta[j+1] + RadiusDelta[j]);
			const double TCRj1 = (Annulus[j].TC[i] * RadiusDelta[j] + Annulus[j-1].TC[i] * RadiusDelta[j-1]) / (RadiusDelta[j] + RadiusDelta[j-1]);
			FluxInterior = Annulus[j].dx[i] * 2.0 * M_PI * Radius[j]   * TCRj  * (Annulus[j+1].T[i] - Annulus[j].T[i]) / ((RadiusDelta[j+1] + RadiusDelta[j]) / 2.0);
			FluxExterior = -Annulus[j].dx[i] * 2.0 * M_PI * Radius[j-1] * TCRj1 * (Annulus[j].T[i] - Annulus[j-1].T[i]) / ((RadiusDelta[j] + RadiusDelta[j-1]) / 2.0);
		}
		Annulus[j].Q[i] = FluxInterior +  FluxExterior;
	}
	for (unsigned int j = 0 ; j < NumberOfRadius ; j++) Annulus[j].Solve(dt);
	//Log("cThermalConductionEquationCylinder::Solve == Sortie");
}

void cThermalConductionEquationCylinder::GetTx(const unsigned int iR, vector<double>& T) const {
	Annulus[iR].GetT(T);
}

void cThermalConductionEquationCylinder::GetTR(const unsigned int ix, vector<double>& T) const {
	vector<double> Temperature(T.size());
	for (unsigned int i = 0 ; i < NumberOfRadius ; i++) {
		Annulus[i].GetT(Temperature);
		T[i] = Temperature[i];
	}
}

double cThermalConductionEquationCylinder::GetFrontBoundaryConditionFlux() const {
	double Flux = 0;
	for (unsigned int j = 0 ; j < NumberOfRadius; j++) Flux += Annulus[j].FrontFlux * Annulus[j].Area;
	return Flux;
}

double cThermalConductionEquationCylinder::GetBackBoundaryConditionFlux() const {
	double Flux = 0;
	for (unsigned int j = 0 ; j < NumberOfRadius; j++) Flux += Annulus[j].BackFlux * Annulus[j].Area;
	return Flux;
}

double cThermalConductionEquationCylinder::GetEnergy() const {
	double Energy = 0;
	for (unsigned int j = 0 ; j < NumberOfRadius; j++) Energy += Annulus[j].GetEnergy();
	return Energy;
}

//======================================================================================================================

cThermalConductionEquationCylinder_Test::cThermalConductionEquationCylinder_Test()
:
	Radius(0.2),
	dxLayer(0.05),
	//FilePrefix("Test_cThermalConductionEquationCylinder")
	FilePrefix(".")
{}

void cThermalConductionEquationCylinder_Test::Test01(
	const unsigned int NumberOfLayers,
	const double ValueFront,
	const double ValueBack
) const
{
	ofstream File;
	File.open(FilePrefix+"_Sortie/Test01.csv");
	Test0x(NumberOfLayers, Radius, 0.0, ValueImposed, ValueFront, ValueImposed, ValueBack, File);
	File.close();
}

void cThermalConductionEquationCylinder_Test::Test02(
	const unsigned int NumberOfLayers,
	const double FluxFront, // W / m2
	const double FluxBack	// W / m2
) const
{
	const double RadiusFlux = 1.0;
	double CheckFlux;
	ofstream File;
	File.open(FilePrefix+"_Sortie/Test02.csv");
	Test0x(NumberOfLayers, Radius, RadiusFlux, FluxImposed, FluxFront, FluxImposed, FluxBack, File);
	File << "Variation d'Energie imposée = " << sep <<
		ts(
			(FluxFront+FluxBack) * M_PI * Radius * Radius +
			((NumberOfLayers * dxLayer) / 1.0) * (2.0 * M_PI * Radius) * RadiusFlux) << sep << "J/s" << endl;
	File.close();
}

void cThermalConductionEquationCylinder_Test::Test03(
	const unsigned int NumberOfLayers,
	const double ValueFront,
	const double FluxBack
) const
{
	ofstream File;
	File.open(FilePrefix+"_Sortie/Test03.csv");
	Test0x(NumberOfLayers, Radius, 1.0, ValueImposed, ValueFront, FluxImposed, FluxBack, File);
	File.close();
}

void cThermalConductionEquationCylinder_Test::Test04(
	const unsigned int NumberOfLayers,
	const double FluxFront,
	const double ValueBack)
const
{
	ofstream File;
	File.open(FilePrefix+"_Sortie/Test04.csv");
	Test0x(NumberOfLayers, Radius, 1.0, FluxImposed, FluxFront, ValueImposed, ValueBack, File);
	File.close();
}

void cThermalConductionEquationCylinder_Test::Test05() const
{
	double FinalTime, dt, Dt, LastTimeSaved;
	double Energy, RadialFluxImposed, FrontFluxImposed, BackFluxImposed;
	vector<double> LocaldRadius;					// m
	vector<double> LocaldxLayer;					// m
	vector<double> ThermalConductivity;				// W m-1 K-1
	vector<double> InterfaceThermalConductivity;	// W K-1
	vector<double> VolumetricHeatCapacity;			// W m-3 K-1
	vector<double> Temperature;
	BoundaryConditionType FrontBoundaryConditionType;
	vector<double> FrontBoundaryConditiont;
	vector<vector<double>>FrontBoundaryConditionValue;
	BoundaryConditionType BackBoundaryConditionType;
	vector<double> BackBoundaryConditiont;
	vector<vector<double>> BackBoundaryConditionValue;
	vector<vector<double>> RadialBoundaryConditionValueOfFlux;
	vector<double> RadialBoundaryConditiont;
	ofstream File;

	//File.open(FilePrefix+"_Sortie/Result.csv");
	File.open(FilePrefix+"/Result.csv");

	cThermalConductionEquationCylinder Equation;
	Log("cThermalConductionEquationCylinder_Test::Test05 == Entree");

	Equation.Read(FilePrefix+"/Input.xml",
		FinalTime,
		dt,
		Dt,
		LocaldRadius, LocaldxLayer,
		ThermalConductivity, InterfaceThermalConductivity,
		VolumetricHeatCapacity,
		Temperature,
		FrontBoundaryConditionType, FrontBoundaryConditiont, FrontBoundaryConditionValue,
		BackBoundaryConditionType, BackBoundaryConditiont, BackBoundaryConditionValue,
		RadialBoundaryConditiont, RadialBoundaryConditionValueOfFlux
		);
	Equation.Initialize(LocaldRadius, LocaldxLayer, ThermalConductivity, InterfaceThermalConductivity, VolumetricHeatCapacity, Temperature);
	Equation.SetStabilityFactor(0.0, 0.0, 100.0);
	Equation.SetFrontBoundaryConditionDataSet(FrontBoundaryConditionType, FrontBoundaryConditiont, FrontBoundaryConditionValue);
	Equation.SetBackBoundaryConditionDataSet(BackBoundaryConditionType, BackBoundaryConditiont, BackBoundaryConditionValue);
	Equation.SetRadialFluxBoundaryConditionDataSet(RadialBoundaryConditiont, RadialBoundaryConditionValueOfFlux);

	File << "t" << sep;
	for (unsigned j = 0 ; j < Equation.NumberOfRadius ; j++) {
		for (unsigned int i = 0 ; i < Temperature.size() ; i++) File << sep << "T_"+to_string(j)+"_"+to_string(i);
		if (j != (Equation.NumberOfRadius-1)) File << sep;
	}
	/*
	File << sep;
	for (unsigned int i = 0 ; i < Temperature.size() ; i++) File << sep << "Ti_"+to_string(i);
	*/
	File << endl;
	LastTimeSaved = -2.0 * Dt;
	for (unsigned int t = 0 ; t <= FinalTime; t += dt) {
		Equation.UpdateFrontBoundaryCondition(t/3600.0);
		Equation.UpdateBackBoundaryCondition(t/3600.0);
		Equation.UpdateRadialBoundaryCondition(t/3600.0);



		Equation.Solve(dt);
		if ((t - LastTimeSaved) >= Dt) {
			File << t << sep;
			for (unsigned j = 0 ; j < Equation.NumberOfRadius ; j++) {
				Equation.GetTx(j, Temperature);
				for (unsigned int i = 0 ; i < Temperature.size() ; i++) File << sep << Temperature[i];
				File <<  sep;
			}

			File <<  endl;
			LastTimeSaved = t;
		}
	}
	File << endl;
	Energy = - Equation.GetEnergy();
	Equation.Solve(dt);
	Energy += Equation.GetEnergy();

	RadialFluxImposed = 0.0;
	for (unsigned int i = 0 ; i < Equation.Annulus.back().GetNumberOfSublayers() ; i++) RadialFluxImposed += Equation.Qbc[i];

	FrontFluxImposed = Equation.GetFrontBoundaryConditionFlux();
	BackFluxImposed = Equation.GetBackBoundaryConditionFlux();
	Log("Flux radial imposé = " + ts(RadialFluxImposed) + "W");
	Log("Flux axial imposé = " + ts(FrontFluxImposed+BackFluxImposed) + "W");
	Log("Flux axial total = " + ts(RadialFluxImposed + FrontFluxImposed + BackFluxImposed) + "W");
	Log("Variation d'Energie calculee = " + ts(Energy / dt) + "W");
	File.close();

	Log("cThermalConductionEquationCylinder_Test::Test05 == Sortie");
}

void cThermalConductionEquationCylinder_Test::Test11(
	const unsigned int NumberOfLayers,
	const double FluxFront,
	const double ValueBack)
const
{
	ifstream FileInput;
	ofstream FileOutput;
	const double dt = 1.0;
	vector<double> LocaldRadius(10);					// m
	vector<double> LocaldxLayer(NumberOfLayers, dxLayer);					// m
	vector<double> ThermalConductivity(NumberOfLayers, 1.0);		// J / m / K
	vector<double> InterfaceThermalConductivity(NumberOfLayers, 0.0);		// J / m / K
	vector<double> VolumetricHeatCapacity(NumberOfLayers, 2.0e6);	// J / m3 / K
	vector<double> Temperature(NumberOfLayers, 0.0);
	vector<double> RadiusFluxLocal(NumberOfLayers, 0.0);
	cThermalConductionEquationCylinder Equation;
	double Energy, CheckFlux;
	Log("cThermalConductionEquationCylinder_Test::Test0x == Entree");


	for (unsigned int i = 0 ; i < LocaldRadius.size() ; i++) LocaldRadius[i] = Radius / LocaldRadius.size();
	Equation.Initialize(LocaldRadius, LocaldxLayer, ThermalConductivity, InterfaceThermalConductivity, VolumetricHeatCapacity, Temperature);
	Equation.SetStabilityFactor(0.0, 0.0, 100.0);
	Equation.SetFrontBoundaryCondition(FluxImposed, 0.0);
	Equation.SetBackBoundaryCondition(FluxImposed, 0.0);
	for (unsigned int i = 0 ; i < NumberOfLayers ; i++)
		if((i > (NumberOfLayers / 3.0)) && (i < (2.0 * NumberOfLayers / 3.0)))
			RadiusFluxLocal[i] = 0.0;
	//for (unsigned int i = 0 ; i < NumberOfLayers ; i++) cout << RadiusFluxLocal[i] << endl; exit(1);
	Equation.SetRadialFluxBoundaryCondition(RadiusFluxLocal);
	/*
	cout << Equation.Annulus.back().TotalNumberOfSublayers << endl;
	for (unsigned int i = 0 ; i < Equation.Annulus.back().TotalNumberOfSublayers ; i++) cout << "Fin == " << Equation.Annulus.back().dx[i] << endl;
	exit(0);//*/
	FileOutput << "t";
	for (unsigned int i = 0 ; i < Temperature.size() ; i++) FileOutput << sep << "Te_"+to_string(i);
	FileOutput << sep;
	for (unsigned int i = 0 ; i < Temperature.size() ; i++) FileOutput << sep << "Ti_"+to_string(i);
	FileOutput << endl;
	for (unsigned int t = 0 ; t < 30 * 1000; t++) {
		Equation.Solve(dt);
		if ((t % 1000) == 0) {
			FileOutput << t * dt;
			Equation.GetTx(Equation.NumberOfRadius-1, Temperature);
			for (unsigned int i = 0 ; i < Temperature.size() ; i++) FileOutput << sep << Temperature[i];
			FileOutput <<  sep;
			Equation.GetTx(0, Temperature);
			for (unsigned int i = 0 ; i < Temperature.size() ; i++) FileOutput << sep << Temperature[i];
			FileOutput <<  endl;
		}
	}
	FileOutput << endl;
	Energy = - Equation.GetEnergy();
	Equation.Solve(dt);
	Energy += Equation.GetEnergy();
	FileOutput << "Variation d'Energie calculee = " << sep << ts(Energy / dt) << sep << "J/s" << endl;
	CheckFlux = 0.0;
	for (unsigned int i = 0 ; i < Equation.Annulus.back().GetNumberOfSublayers() ; i++) {
		CheckFlux += Equation.Qbc[i];
		//cout << Equation.Qbc[i] << endl;
	}
	FileOutput << "Flux radial imposé = " << sep << ts(CheckFlux) << sep << "W" << endl;

	FileOutput.close();
	Log("cThermalConductionEquationCylinder_Test::Test0x == Sortie");
}

void cThermalConductionEquationCylinder_Test::Test0x(
	const unsigned int NumberOfLayers,
	const double Radius, const double RadiusFlux,
	const BoundaryConditionType BCFront, const double XFront,
	const BoundaryConditionType BCBack, const double XBack,
	ofstream& File) const
{
	const double dt = 1.0;
	vector<double> LocaldRadius(10);					// m
	vector<double> LocaldxLayer(NumberOfLayers, dxLayer);					// m
	vector<double> ThermalConductivity(NumberOfLayers, 1.0);		// J / m / K
	vector<double> InterfaceThermalConductivity(NumberOfLayers, 0.0);		// J / m / K
	vector<double> VolumetricHeatCapacity(NumberOfLayers, 2.0e6);	// J / m3 / K
	vector<double> Temperature(NumberOfLayers, 0.0);
	vector<double> RadiusFluxLocal(NumberOfLayers, 0.0);
	cThermalConductionEquationCylinder Equation;
	double Energy, CheckFlux;
	Log("cThermalConductionEquationCylinder_Test::Test0x == Entree");

	for (unsigned int i = 0 ; i < LocaldRadius.size() ; i++) LocaldRadius[i] = Radius / LocaldRadius.size();

	Equation.Initialize(LocaldRadius, LocaldxLayer, ThermalConductivity, InterfaceThermalConductivity, VolumetricHeatCapacity, Temperature);

	Equation.SetStabilityFactor(0.0, 0.0, 100.0);
	Equation.SetFrontBoundaryCondition(BCFront, XFront);
	Equation.SetBackBoundaryCondition(BCBack, XBack);
	for (unsigned int i = 0 ; i < NumberOfLayers ; i++)
		if((i > (NumberOfLayers / 3.0)) && (i < (2.0 * NumberOfLayers / 3.0)))
			RadiusFluxLocal[i] = RadiusFlux;
	//for (unsigned int i = 0 ; i < NumberOfLayers ; i++) cout << RadiusFluxLocal[i] << endl; exit(1);
	Equation.SetRadialFluxBoundaryCondition(RadiusFluxLocal);
	/*
	cout << Equation.Annulus.back().TotalNumberOfSublayers << endl;
	for (unsigned int i = 0 ; i < Equation.Annulus.back().TotalNumberOfSublayers ; i++) cout << "Fin == " << Equation.Annulus.back().dx[i] << endl;
	exit(0);//*/
	File << "t";
	for (unsigned int i = 0 ; i < Temperature.size() ; i++) File << sep << "Te_"+to_string(i);
	File << sep;
	for (unsigned int i = 0 ; i < Temperature.size() ; i++) File << sep << "Ti_"+to_string(i);
	File << endl;
	for (unsigned int t = 0 ; t < 30 * 1000; t++) {
		Equation.Solve(dt);
		if ((t % 1000) == 0) {
			File << t * dt;
			Equation.GetTx(Equation.NumberOfRadius-1, Temperature);
			for (unsigned int i = 0 ; i < Temperature.size() ; i++) File << sep << Temperature[i];
			File <<  sep;
			Equation.GetTx(0, Temperature);
			for (unsigned int i = 0 ; i < Temperature.size() ; i++) File << sep << Temperature[i];
			File <<  endl;
		}
	}
	File << endl;
	Energy = - Equation.GetEnergy();
	Equation.Solve(dt);
	Energy += Equation.GetEnergy();
	File << "Variation d'Energie calculee = " << sep << ts(Energy / dt) << sep << "J/s" << endl;
	CheckFlux = 0.0;
	for (unsigned int i = 0 ; i < Equation.Annulus.back().GetNumberOfSublayers() ; i++) {
		CheckFlux += Equation.Qbc[i];
		//cout << "Fin == " << Equation.Annulus.back().dx[i] << " / " << Equation.Qbc[i] << " / " << CheckFlux << endl;
	}
	File << "Flux radial imposé = " << sep << ts(CheckFlux) << "W" << endl;

	cout << "Flux radial imposé = " << ts(CheckFlux) << "W" << endl;
	cout << "Flux axial imposé = " << ts((XFront+XBack) * (M_PI * Radius * Radius)) << "W" << endl;
	cout << "Flux axial total = " << ts(CheckFlux + (XFront+XBack) * (M_PI * Radius * Radius)) << "W" << endl;
	cout << "Variation d'Energie calculee = " << ts(Energy / dt) << "W" << endl;
	Log("cThermalConductionEquationCylinder_Test::Test0x == Sortie");
}
