
/*======================================================================================================================
* FILE : Define.h							VERSION : beta 1.0		DATE : 18 december 2009		INITIAL DEVELOPER : LCPC
* DESCRIPTION : Défintion des variables de contrôle de la conpilation
* A detailed list of autors and modifications is given the end of the file
*
* This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
* Public License. This software is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
* merchandability or fitness for a particular purpose.
======================================================================================================================*/

/*
template <class tData>
cVector<tData>::cVector() {}

template <class tData>
cVector<tData>::cVector(const tData& D1, const tData& D2)
:
	vector<tData>(2)
{
	this->at(0) = D1;
	this->at(1) = D2;
}

template <class tData>
cVector<tData>::cVector(const tData& D1, const tData& D2, const tData& D3)
:
	vector<tData>(3)
{
	this->at(0) = D1;
	this->at(1) = D2;
	this->at(2) = D3;
}

template <class tData>
cVector<tData>::cVector(const tData& D1, const tData& D2, const tData& D3, const tData& D4)
:
	vector<tData>(4)
{
	this->at(0) = D1;
	this->at(1) = D2;
	this->at(2) = D3;
	this->at(3) = D4;
}

template <class tData>
cVector<tData>::cVector(const tData& D1, const tData& D2, const tData& D3, const tData& D4, const tData& D5)
:
	vector<tData>(5)
{
	this->at(0) = D1;
	this->at(1) = D2;
	this->at(2) = D3;
	this->at(3) = D4;
	this->at(4) = D5;
}

template <class tData>
cVector<tData>::cVector(const tData& D1, const tData& D2, const tData& D3, const tData& D4, const tData& D5, const tData& D6)
:
	vector<tData>(6)
{
	this->at(0) = D1;
	this->at(1) = D2;
	this->at(2) = D3;
	this->at(3) = D4;
	this->at(4) = D5;
	this->at(5) = D6;
}

template <class tData>
cVector<tData>::cVector(const tData& D1, const tData& D2, const tData& D3, const tData& D4, const tData& D5, const tData& D6, const tData& D7)
:
	vector<tData>(7)
{
	this->at(0) = D1;
	this->at(1) = D2;
	this->at(2) = D3;
	this->at(3) = D4;
	this->at(4) = D5;
	this->at(5) = D6;
	this->at(6) = D7;
}
*/
//======================================================================================================================

template <class tData>
cMatrix<tData>::cMatrix(size_t _nI, size_t _nJ)
:
	nI(_nI),
	nJ(_nJ),
	Data(_nI * _nJ)
{}

template <class tData>
void cMatrix<tData>::resize(size_t _nI, size_t _nJ) {
	nI = _nI;
	nJ = _nJ;
	Data.resize(_nI * _nJ);
}

template <class tData>
tData& cMatrix<tData>::operator()(size_t i, size_t j) {return Data[i * nJ + j];}

template <class tData>
tData cMatrix<tData>::operator()(size_t i, size_t j) const {return Data[i * nJ + j];}

//======================================================================================================================

template <class tData>
cMatrix3<tData>::cMatrix3(size_t _nI, size_t _nJ, size_t _nK)
:
	nI(_nI),
	nJ(_nJ),
	nK(_nK),
	mData(nI*nJ*nK)
{}

template <class tData>
void cMatrix3<tData>::resize(size_t _nI, size_t _nJ, size_t _nK) {
	nI = _nI;
	nJ = _nJ;
	nK = _nK;
	mData.resize(nI*nJ*nK);
}

template <class tData>
tData& cMatrix3<tData>::operator()(size_t i, size_t j, size_t k) {return mData[i*nJ*nK + j*nK + k];}

template <class tData>
tData cMatrix3<tData>::operator()(size_t i, size_t j, size_t k) const {
	const unsigned int Pos = i*nJ*nK + j*nK + k;
	//cout << "tData cMatrix3<tData>::operator() == Entree" << endl;
	//cout << "tData cMatrix3<tData>::operator() == Sortie" << endl;
	//return Data[i*nJ*nK + j*nK + k];
	return mData[Pos];
}

//======================================================================================================================

template <class tData>
cMatrix3ToPointer<tData>::cMatrix3ToPointer(size_t _nI, size_t _nJ, size_t _nK)
:
	cMatrix3<tData>(_nI, _nJ, _nK)
{}

template <class tData>
cMatrix3ToPointer<tData>::~cMatrix3ToPointer()
{
	for (unsigned int i = 0 ; i < cMatrix3<tData>::mData.size(); i++) delete cMatrix3<tData>::mData[i];
}

