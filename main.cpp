#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <string>
#include <math.h>

using namespace std;

#include "Utils.hh"
#include "ThermalConduction.hh"

//ofstream FileErrorMessage;
//ofstream FileLogMessage;

void Test_ThermalConduction()
{
	cThermalConductionEquation1D Equation;
	cThermalConductionEquation1D_Test TestThermalConduction;

	//Equation.Initialize(Thickness, ThermalConductivity, VolumetricHeatCapacity, T); return 0;

	TestThermalConduction.Test01(10, 20.0, 10.0);
	TestThermalConduction.Test02(10, 10.0, 20.0);
	TestThermalConduction.Test03(10, 10.0, 20.0);
	TestThermalConduction.Test04(10, 10.0, 20.0);
}

void Test_ThermalConductionCylinder()
{
	cThermalConductionEquationCylinder_Test TestThermalConductionCylinder;

	//TestThermalConductionCylinder.Test01(10, 20.0, 10.0);
	//TestThermalConductionCylinder.Test02(10, 10.0, 10.0);
	//TestThermalConductionCylinder.Test03(10, 10.0, 20.0);
	//TestThermalConductionCylinder.Test04(10, 10.0, 20.0);
	TestThermalConductionCylinder.Test05();
}

int main()
{
	Test_ThermalConductionCylinder();
	return 0;
}


