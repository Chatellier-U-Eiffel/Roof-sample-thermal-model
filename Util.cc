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
#include "cTemplate.hh"

//======================================================================================================================

string ts(const double X) {return std::to_string(X);}
string tsi(const int X) {return std::to_string(X);}

/*! \class cLog
* \brief A test class.
*
* A more detailed class description.
*/

/*! \fn cLog::cLog(const bool Reset)
* \brief Le constructeur
*
* Description détaillée du constructeur<br>
* Sur plusieurs lignes
*/
cLog::cLog(const bool Reset)
:
	LogFileName("FichierLog.txt"),
	Active(true)
{
	if (Reset == true) LogFile.open(LogFileName); else LogFile.open(LogFileName, ios::app);
}

cLog::~cLog() {LogFile.close();}

void cLog::SetActive(const bool _Active) {Active = _Active;}

void cLog::Reset() {LogFile.close(); LogFile.open(LogFileName);}

void cLog::operator()(const string LogMessage, const bool EndLine) {
	if (Active == false) return;
	else {
		LogFile << LogMessage;
		if (EndLine == true) LogFile << endl;
	}
}

void cLog::operator()(const vector<double>& V, const string Caption) {
	if (Active == false) return;
	else {
		WriteVectorHeader(V, Caption, LogFile); LogFile << endl;
		WriteVector(V, LogFile); LogFile << endl;
	}
}

ofstream& cLog::GetFile() {return LogFile;}

cLog Log(true);

//======================================================================================================================

/*! \fn int IntegerPart(const double X)
* \brief Integer part of a double
*
* \return Integer part of X
*/
int IntegerPart(const double X) {if (X >= 0) return floor(X); else return floor(X) + 1;}

/*! \fn int Sign(const double X)
* \brief Analysis of X sign
*
* \return 1.0 if X < 0.0<BR>
* -1.0 if X < 0.0<BR>
* 0 if X = 0.0
*/
double Sign(const double X) {if (X > 0.0) return 1.0; else if (X < 0.0) return -1.0; else return 0.0;}

/*! \fn void DisplayVector(const vector<double>& V, const string Caption)
* \brief Writes a vector to the standard output
*
* Writes first the string "Caption" to the stahdard output<BR>
* then writes the content of vector X as column
*/
void DisplayVector(const vector<double>& V, const string Caption) {
	cout << Caption << endl;
	for (unsigned int i = 0; i < V.size(); i++) cout << i << " => " << V[i] << endl;
}

/*! \fn void WriteVectorHeader(const vector<double>& V, const string Caption, ofstream& File)
* \brief Writes the meaning of vector velues to File outstream
*
* Writes first the string "Caption" to the stahdard output<BR>
* then writes the content of vector X as column
*/
void WriteVectorHeader(const vector<double>& V, const string Caption, ofstream& File) {
	for (unsigned int i = 0 ; i < V.size()-1 ; i++) File << Caption + "_" + tsi(i) << sep;
	File << Caption + tsi(V.size()-1);
}

/*! \fn void WriteVector(const vector<double>& V, ofstream& File, const bool Column)
* \brief Writes the values of a vector to File outstream
*
* Writes first the string "Caption" to the stahdard output<BR>
* then writes the content of vector X as column
* \param Column : Writes de values of the vector as a colum if true and as row if false
*/
void WriteVector(const vector<double>& V, ofstream& File, const bool Column) {
	if (Column == false) {
		for (unsigned int i = 0 ; i < V.size()-1 ; i++) File << V[i] << sep;
		File << V.back();
	} else for (unsigned int i = 0 ; i < V.size() ; i++) File << V[i] << endl;;
}

/*! \fn void FatalError(const string Message)
* \brief Stops the program execiution in case of fatal error
*
* Writes fatal error message to the stahdard output and stops the program with an exit code
* The fatal error message indicates the name of the funciton causing the fatal eror and an error code that locates
* the error in thefunction
* \param Exit code
*/
void FatalError(const string Message) {
	cout << Message << endl;
	exit(1);
}

/*! \fn string ParseS(const string S, unsigned int Col)
* \brief Anslyses a string containing information seperated by comas
*
* \param String to be analysed
* \return String extracted by the analysis of the arbument
*/
string ParseS(const string S, unsigned int Col) {
	unsigned int iCol, I, Colstart, Colend;
	bool StartValid, EndValid;
	string S1;
	//Log("ParseS == Entree");
	iCol = 0;
	I = 0;
	StartValid = false;
	EndValid = false;

	while (I < S.size()) {
		//cout << I << endl;
		if (S[I] == ',') {
			//cout << "Virgule = " << I << endl;
			iCol++;
			if (Col == 0) {Colstart = 0; Colend = I-1; StartValid = true; EndValid = true;}
			else if (iCol == Col) {Colstart = I+1; StartValid = true;}
			else if (iCol == (Col+1)) {Colend = I-1; EndValid = true;}
		}
		I++;
	}
	//cout << "Sortie While" << endl;
	if (Colend == 0) Colend = S.size() -1;
	if (StartValid == false) {cout << "Parse == Erreur 01" << endl; exit(1);}
	if (EndValid == false) {Colend = S.size()-1; EndValid = true;}
	//cout << Colstart << endl;
	//cout << Colend << endl;
	S1.resize(Colend-Colstart+1);
	for (iCol = 0; iCol < S1.size() ; iCol++) S1[iCol] = S[Colstart+iCol];
	//Log("ParseS == Sortie");
	return(S1);

}

double Parse(const string S, unsigned int Col) {return(::atof(ParseS(S, Col).c_str()));}

//======================================================================================================================

void ParseXML(const string FileName, const string Balise, std::vector<double>& Result) {
	const string StartBalise = string("<") + Balise + string(">");
	const string StopBalise = string("</") + Balise + string(">");
	ifstream File;
	string Line;

	File.open(FileName);
	if (File.is_open()) {
		while (getline(File, Line) && (Line != StartBalise)) {}
		Result.resize(0);
		while (getline(File, Line) && (Line != StopBalise)) {
			Result.push_back(::atof(Line.c_str()));
		}
		File.close();
	}
	else FatalError("ParseXML == Error 01");
}

void ParseXML(const string FileName, const string Balise, double& Parameter, string& String) {
	const string StartBalise = string("<") + Balise + string(">");
	const string StopBalise = string("</") + Balise + string(">");
	ifstream File;
	string Line;

	File.open(FileName);
	if (File.is_open()) {
		while (getline(File, Line) && (Line != StartBalise)) {}
		getline(File, Line);
		Parameter = ::atof(Line.c_str());
		getline(File, String);
		File.close();
	}
	else FatalError("ParseXML == Error 01");
}

//======================================================================================================================

template class cMatrix<double>;
template class cMatrix<string>;

