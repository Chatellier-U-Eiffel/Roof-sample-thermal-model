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

//======================================================================================================================*/

//extern ofstream FileErrorMessage;
//extern ofstream FileLogMessage;

class cLog {
	public :
		cLog(const bool Reset = true);
		~cLog();
	public :
		void SetActive(const bool Active = true);
		void Reset();
		void operator()(const string LogMessage, const bool EndLine = true);
		void operator()(const vector<double>& V, const string Caption);
		ofstream& GetFile();
	private :
		const string LogFileName;
		ofstream LogFile;
		bool Active;
};

const char tab = '\t';
const string newline = " / ";
const char sep = ',';
extern cLog Log;

int IntegerPart(const double X);
double Sign(const double X);
string ts(const double X);
string tsi(const int X);
void DisplayVector(const vector<double>& V, const string Caption);
void WriteVectorHeader(const vector<double>& V, const string Caption, ofstream& File);
void WriteVector(const vector<double>& V, ofstream& File, const bool Column = false);
void FatalError(const string Message);

string ParseS(const string S, unsigned int Col);
double Parse(const string S, unsigned int Col);
void ParseXML(const string FileName, const string Balise, std::vector<double>& Result);
void ParseXML(const string FileName, const string Balise, double& Parameter, string& String);

template <class tData>
class cMatrix {
	public:
		cMatrix(size_t nI, size_t nJ);
		tData& operator()(size_t i, size_t j);
		tData operator()(size_t i, size_t j) const;
		void resize(size_t nI, size_t nJ);
	public:
		size_t nI, nJ;
		std::vector<tData> Data;
};

template <class tData>
class cMatrix3 {
	public:
		cMatrix3(size_t nI, size_t nJ, size_t nK);
		tData& operator()(size_t i, size_t j, size_t k);
		tData operator()(size_t i, size_t j, size_t k) const;
		void resize(size_t nI, size_t nJ, size_t nK);
	public:
		size_t nI, nJ, nK;
		std::vector<tData> mData;
};

template <class tData>
class cMatrix3ToPointer : public  cMatrix3<tData> {
	public:
		cMatrix3ToPointer(size_t nI, size_t nJ, size_t nK);
		~cMatrix3ToPointer();
};
