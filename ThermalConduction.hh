/**=====================================================================================================================
* \file cRecherche.hh
* \version 0.1
* \date 6 septembre 2007
* \author Chatellier P. - Univertité Gustave Eiffel (France)
* \brief Programme de modélisation du climat urbain
*
* This library an open source software distributed under the terms of the GNU Lesser General Public License.
* This software is distributed WITHOUT ANY WARRANTY; without even the implied warranty of merchandability or fitness
* for a particular purpose.
*/

//=====================================================================================================================

//! \enum BoundaryConditionType \brief Boundary conditions for the resolution of teh Fourrier equation
typedef enum {ValueImposed, FluxImposed, UndifinedBOundaryCondition} BoundaryConditionType;

//=====================================================================================================================

class cDiscretisationToolForThermalConductionEquation {
	public :
		cDiscretisationToolForThermalConductionEquation();
	public :
		double
			dxMin, // m
			Ratio;
		unsigned int nLayersMin;
	public :
		void Compute(
			const double& dxLayer,
			unsigned int& nbSubLayers,
			double& FirstThicknessSubLayers) const;
		void FillWidth(const double dxMin, const double Ratio, vector<double>& dx) const;
	private :
		void Compute(
			const vector<double>& _dxLayer,
			vector<unsigned int>& nbSubLayers,
			vector<double>& FirstThicknessSubLayers) const;
		void Compute_Old(
			const vector<double>& _dxLayer,
			vector<unsigned int>& nbSubLayers,
			vector<double>& FirstThicknessSubLayers) const;
};

//! Résolution de léquation de Fourrier
class cThermalConductionEquation1D {
	public :
		cThermalConductionEquation1D();
	public :
		//! Initialisation du solveur
		void Initialize(
			const vector<double>& dxLayer, //!< dx couche
			const vector<double>& TCLayer,
			const vector<double>& TCInterface,
			const vector<double>& VHCLayer,
			const vector<double>& TLayer);
		void SetStabilityFactor(const double& ImplicitFactor, const double& RelaxationFactor, const double& MaxDelta);

		void SetFrontBoundaryConditionDataSet(const BoundaryConditionType& BCFrontType, const vector<double>& t, vector<double>& Value);
		void SetFrontBoundaryCondition(const BoundaryConditionType FrontBoundaryCondition, const double Value);
		void UpdateFrontBoundaryCondition(const double& t);

		void SetBackBoundaryCondition(const BoundaryConditionType BackBoundaryCondition, const double Value);
		void SetBackBoundaryConditionDataSet(const BoundaryConditionType& BCBackType, const vector<double>& t, vector<double>& Value);
		void UpdateBackBoundaryCondition(const double& t);
	public :
		void Solve(const double& dt);//, vector<double>& T);
		unsigned int GetNumberOfLayers() const;
		unsigned int GetNumberOfSublayers() const;
		double GetBCFront(const double& t);
		double GetBCBack(const double& t);
		double GetTFront() const;
		double GetTBack() const;
		void GetT(vector<double>& T) const;
	public :
		double GetEnergy() const;
	public :
		double Area;
		vector<vector<double> > dxOfSublayers;
		double ExplicitFactor, ImplicitFactor,RelaxationFactor, MaxDelta;
		BoundaryConditionType FrontBoundaryCondition, BackBoundaryCondition;
		double FrontValue, FrontFlux, BackValue, BackFlux;
		vector<double>
			dx,		//!< Thickness of sublayer (m)
			TC,		//!< thermal conductivity (W m-1 K-1)
			TCInterface,//!< thermal conductivity (W m-2 K-1)
			VHC,	//!< volumetric heat capacity (J m-3 K-1)
			Q,		//!< Added flux (W m-3)
			T;		//!< Temperature in sublayers
	private :
		void Initializedx1(const vector<double>& dxLayer, vector<vector<double> >& dxOfSublayers, vector<double>& dx);
		void Initializedx2(const vector<double>& dxLayer, vector<vector<double> >& dxOfSublayers, vector<double>& dx);
		void SolveFiniteDifference(const double& dt, vector<double>& T);
		void SolveFiniteVolume(const double& dt);//, vector<double>& T);
	private :
		vector<double> dxOfLayers;
		unsigned int LastBCFrontt, LastBCBackt;
		BoundaryConditionType BCFrontType, BCBackType;
		vector<double> BCFrontt, BCFrontValue, BCBackt, BCBackValue;
		//cDiscretisationToolForThermalConductionEquation DiscretisationTool2;
};

class cThermalConductionEquationCylinder {
	public :
		cThermalConductionEquationCylinder();
	public :
		void Initialize(
			const vector<double>& dRadius,
			const vector<double>& dxLayer, //!< dx couche
			const vector<double>& TCLayer,
			const vector<double>& TCInterface,
			const vector<double>& VHCLayer,
			const vector<double>& TLayer);
		void ReadFile(const string& FileName , const unsigned int col, vector<double>& Result);
		void Read(
			const string FileName,
			double& FinalTime, double& dt, double& Dt,
			vector<double>& dRadius, vector<double>& dxOfLayers,
			vector<double>& ThermalConductivity, vector<double>& InterfaceThermalConductivity,
			vector<double>& VolumetricHeatCapacity,
			vector<double>& InitialTemperature,
			BoundaryConditionType& FrontBoundaryConditionType, vector<double>& FrontBoundaryConditiont, vector<vector<double>>& FrontBoundaryConditionValue,
			BoundaryConditionType& BackBoundaryConditionType, vector<double>& BackBoundaryConditiont, vector<vector<double>>& BackBoundaryConditionValue,
			vector<double>& RadialBoundaryConditiont, vector<vector<double>>& RadialBoundaryConditionValueOfFlux);
		void SetStabilityFactor(const double& ImplicitFactor, const double& RelaxationFactor, const double& MaxDelta);

		void SetFrontBoundaryCondition(const BoundaryConditionType FrontBoundaryCondition, const double Value);
		void SetFrontBoundaryCondition(const vector<BoundaryConditionType>& FrontBoundaryCondition, const vector<double>& Value);
		void SetFrontBoundaryConditionDataSet(const BoundaryConditionType& BCFrontType, const vector<double>& t, vector<vector<double>>& Value);

		void SetBackBoundaryCondition(const BoundaryConditionType BackBoundaryCondition, const double Value);
		void SetBackBoundaryCondition(const vector<BoundaryConditionType>& BackBoundaryCondition, const vector<double>& Value);
		void SetBackBoundaryConditionDataSet(const BoundaryConditionType& BCFrontType, const vector<double>& t, vector<vector<double>>& Value);

		void SetRadialFluxBoundaryCondition(const vector<double>& RadiusFlux);
		void SetRadialFluxBoundaryConditionDataSet(const vector<double>& t, const vector<vector<double>>& BCradialValue);

		void UpdateFrontBoundaryCondition(const double& t);
		void UpdateBackBoundaryCondition(const double& t);
		void UpdateRadialBoundaryCondition(const double& t);


		void Solve(const double& dt);
	public :
		//void GetT(vector<double>& T) const;
		void GetTx(const unsigned int iR, vector<double>& T) const;
		void GetTR(const unsigned int ix, vector<double>& T) const;
		double GetFrontBoundaryConditionFlux() const; // w  intégrés sur la surface
		double GetBackBoundaryConditionFlux() const; // w intégrés sur la surface
		double GetEnergy() const;
	public :
		unsigned int NumberOfRadius;
		vector<double> Radius, RadiusDelta, Qbc;// RadiusFlux;
		vector<cThermalConductionEquation1D> Annulus;
	private :
		unsigned int LastBCRadialt;
		vector<double> BCRadialt;
		vector<vector<double>> BCRadialValue;
//cDiscretisationToolForThermalConductionEquation DiscretisationTool2;
};

//=====================================================================================================================

class cThermalConductionEquation1D_Test {
	public :
		cThermalConductionEquation1D_Test();
	public :
		void Test01(const unsigned int NumberOfLayers, const double ValueFront, const double ValueBack) const;
		void Test02(const unsigned int NumberOfLayers, const double FluxFront, const double FluxBack) const;
		void Test03(const unsigned int NumberOfLayers, const double ValueFront, const double FluxBack) const;
		void Test04(const unsigned int NumberOfLayers, const double FluxFront, const double ValueBack) const;

	private :
		void Test0x(const unsigned int NumberOfLayers,
			const BoundaryConditionType BCFront, const double XFront,
			const BoundaryConditionType BCBack, const double XBack,
			ofstream& File) const;
	private :
		double dxLayer;
		const string FilePrefix;
};

class cThermalConductionEquationCylinder_Test {
	public :
		cThermalConductionEquationCylinder_Test();
	public :
		void Test01(const unsigned int NumberOfLayers, const double ValueFront, const double ValueBack) const;
		void Test02(const unsigned int NumberOfLayers, const double FluxFront, const double FluxBack) const;
		void Test03(const unsigned int NumberOfLayers, const double ValueFront, const double FluxBack) const;
		void Test04(const unsigned int NumberOfLayers, const double FluxFront, const double ValueBack) const;
		void Test05() const;
	public :
		void Test11(const unsigned int NumberOfLayers, const double FluxFront, const double ValueBack) const;
	private :
		void Test0x(const unsigned int NumberOfLayers,
			const double Radius, const double RadiusFlux,
			const BoundaryConditionType BCFront, const double XFront,
			const BoundaryConditionType BCBack, const double XBack,
			ofstream& File) const;
	private :
		double Radius, dxLayer;
		const string FilePrefix;
};

