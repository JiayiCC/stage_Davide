// TRUST_NO_INDENT
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <Process.h>
#include <PrePostNN.h>

PrePostNN::PrePostNN(string filename)
{
  string buffer, tag;
  size_t npos;
  ifstream f(filename,ios::in);

  if(f){
    while(getline(f, buffer)){
		buffer = trim(buffer);
		tag = "OLD_ALPHA:";
		npos = buffer.find(tag);
		if(npos != string::npos) alpha = ReadDataFromLine(buffer,tag,npos);
		tag = "LAMBDA_MEAN:";
		npos = buffer.find(tag);
		if(npos != string::npos) lmean = ReadDataFromLine(buffer,tag,npos);
		tag = "LAMBDA_MAX:";
		npos = buffer.find(tag);
		if(npos != string::npos) lmax = ReadDataFromLine(buffer,tag,npos);
		tag = "T_FN:";
		npos = buffer.find(tag);
		if(npos != string::npos) tfn = ReadDataFromLine(buffer,tag,npos);
		tag = "B_SIGMA:";
		npos = buffer.find(tag);
		if(npos != string::npos) bsigma = ReadOneDataFromLine(buffer,tag,npos);
		tag = "LAMBDA_AU:";
		npos = buffer.find(tag);
		if(npos != string::npos) lambda_au = ReadDataFromSeveralLines(f,5);
		tag = "LAMBDA_AS:";
		npos = buffer.find(tag);
		if(npos != string::npos) lambda_as = ReadDataFromSeveralLines(f,5);
		tag = "T_THRESH:";
		npos = buffer.find(tag);
		if(npos != string::npos) t_thresh = ReadOneDataFromLine(buffer,tag,npos);
		tag = "PP_LAMBDA:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppl = ReadPPLFromLine(buffer,tag,npos);
		tag = "PP_T:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppt = ReadPPTFromLine(buffer,tag,npos);
		tag = "T0:";
		npos = buffer.find(tag);
		if(npos != string::npos) t0 = ReadDataFromLine(buffer,tag,npos); //from line belowmine start
		tag = "ALPHA_MAX:";
		npos = buffer.find(tag);
		if(npos != string::npos) alpha_max = ReadOneDataFromLine(buffer,tag,npos);
		tag = "ALPHA_KEPS_MAX:";
		npos = buffer.find(tag);
		if(npos != string::npos) alpha_keps_max = ReadOneDataFromLine(buffer,tag,npos);
		tag = "PP_ALPHA:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppalpha = ReadPPAlphaFromLine(buffer,tag,npos);
		tag = "PP_ALPHA_KEPS:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppalpha_keps = ReadPPAlphaKEPSFromLine(buffer,tag,npos);
		tag = "Y_PLUS_MAX_LOG:";
		npos = buffer.find(tag);
		if(npos != string::npos) y_plus_max_log = ReadOneDataFromLine(buffer,tag,npos);
		tag = "PP_Y_PLUS:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppy_plus = ReadPPYPlusFromLine(buffer,tag,npos);
		tag = "Y_PLUS_KEPS_MAX_LOG:";
		npos = buffer.find(tag);
		if(npos != string::npos) y_plus_keps_max_log = ReadOneDataFromLine(buffer,tag,npos);
		tag = "PP_Y_PLUS_KEPS:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppy_plus_keps = ReadPPYPlusKEPSFromLine(buffer,tag,npos);
		tag = "RE_TAU_MAX:";
		npos = buffer.find(tag);
		if(npos != string::npos) Re_t_max = ReadOneDataFromLine(buffer,tag,npos);
		tag = "RE_TAU_KEPS_MAX:";
		npos = buffer.find(tag);
		if(npos != string::npos) Re_t_keps_max = ReadOneDataFromLine(buffer,tag,npos);
		tag = "PP_RE_TAU_RST:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppre_tau = ReadPPReTauFromLine(buffer,tag,npos);
		tag = "PP_RE_TAU_KEPS:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppre_tau_keps = ReadPPReTauKEPSFromLine(buffer,tag,npos);
		tag = "PP_B:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppb = ReadPPBFromLine(buffer,tag,npos);
		tag = "NUM_T:";
		npos = buffer.find(tag);
		if(npos != string::npos) num_t = ReadOneDataFromLine(buffer,tag,npos);
		tag = "K_LS_MAX:";
		npos = buffer.find(tag);
		if(npos != string::npos) k_input_max = ReadOneDataFromLine(buffer,tag,npos);
		tag = "PP_K_LS:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppk_input = ReadPPKinputFromLine(buffer,tag,npos);
		tag = "EPS_LS_MAX:";
		npos = buffer.find(tag);
		if(npos != string::npos) eps_input_max = ReadOneDataFromLine(buffer,tag,npos);
		tag = "PP_EPS_LS:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppeps_input = ReadPPEPSinputFromLine(buffer,tag,npos);
		tag = "K_DNS_MAX:";
		npos = buffer.find(tag);
		if(npos != string::npos) k_output_max = ReadOneDataFromLine(buffer,tag,npos);
		tag = "PP_K_DNS:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppk_output = ReadPPKoutputFromLine(buffer,tag,npos);
		tag = "EPS_DNS_MAX:";
		npos = buffer.find(tag);
		if(npos != string::npos) eps_output_max = ReadOneDataFromLine(buffer,tag,npos);
		tag = "PP_EPS_DNS:";
		npos = buffer.find(tag);
		if(npos != string::npos) ppeps_output = ReadPPEPSoutputFromLine(buffer,tag,npos);
		tag = "DATADIR:";
		npos = buffer.find(tag);
		if(npos != string::npos) datadir = ReadPPLFromLine(buffer,tag,npos);

    }
    f.close();
  }
  else{
    cerr << "fichier: " << filename << " inexistant!" << endl;
    Process::exit();
  }
}

PrePostNN::~PrePostNN()
{
}

void PrePostNN::display(string tag,vector<double> vec)
{
  cout << tag << " ";
  for(unsigned int i=0;i<vec.size();i++)
    cout << vec[i] << " ";
  cout << endl;
}

void PrePostNN::display(string tag,vector<int> vec)
{
  cout << tag << " ";
  for(unsigned int i=0;i<vec.size();i++)
    cout << vec[i] << " ";
  cout << endl;
}

void PrePostNN::display(string tag,vector<vector<double>> mat)
{
  cout << tag << endl;
  for(unsigned int i=0;i<mat.size();i++){
    for(unsigned int j=0;j<mat[i].size();j++)
      cout << mat[i][j] << " ";
    cout << endl;
  }
}

string PrePostNN::trim(const std::string& str, const std::string& whitespace)
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

vector<double> PrePostNN::ReadDataFromLine(string buffer,string tag,size_t npos)
{
  vector<double> ret;
  double val;
  string sval;
  size_t ltag;

  ltag = tag.length();
  sval = buffer.substr(npos+ltag,buffer.length()-ltag);
  istringstream iss(sval);
  while(!iss.eof()){
    iss >> val;
    ret.push_back(val);
  }

  return(ret);
}

double PrePostNN::ReadOneDataFromLine(string buffer,string tag,size_t npos)
{
  double ret;
  string sval;
  size_t ltag;

  ltag = tag.length();
  sval = buffer.substr(npos+ltag,buffer.length()-ltag);
  istringstream iss(sval);
  iss >> ret;

  return(ret);
}

enum pp_T PrePostNN::ReadPPTFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_T ret = INDEFT;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("TR") == 0 ) ret = TR;
  else if( tmp.compare("TF") == 0 ) ret = TF;
  else if ( tmp.compare("FROT") == 0 ) ret = FROT;

  return(ret);
}

enum pp_lambda PrePostNN::ReadPPLFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_lambda ret = INDEFL;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("LNORM") == 0 ) ret = LNORM;
  else if( tmp.compare("LU") == 0 ) ret = LU;
  else if( tmp.compare("LUS") == 0 ) ret = LUS;

  return(ret);
}

enum pp_alpha PrePostNN::ReadPPAlphaFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_alpha ret = INDEFALPHA;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("MAXN") == 0 ) ret = MAXN_A;

  return(ret);
}

enum pp_alpha_keps PrePostNN::ReadPPAlphaKEPSFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_alpha_keps ret = INDEFALPHA_KEPS;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("MAXN") == 0 ) ret = MAXN_A_KEPS;

  return(ret);
}

enum pp_k_input PrePostNN::ReadPPKinputFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_k_input ret = INDEFK_INPUT;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("MAXN") == 0 ) ret = MAXN_K_INPUT;

  return(ret);
}

enum pp_eps_input PrePostNN::ReadPPEPSinputFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_eps_input ret = INDEFEPS_INPUT;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("MAXN") == 0 ) ret = MAXN_EPS_INPUT;

  return(ret);
}

enum pp_k_output PrePostNN::ReadPPKoutputFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_k_output ret = INDEFK_OUTPUT;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("MAXN") == 0 ) ret = MAXN_K_OUTPUT;

  return(ret);
}

enum pp_eps_output PrePostNN::ReadPPEPSoutputFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_eps_output ret = INDEFEPS_OUTPUT;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("MAXN") == 0 ) ret = MAXN_EPS_OUTPUT;

  return(ret);
}

enum pp_y_plus PrePostNN::ReadPPYPlusFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_y_plus ret = INDEFY_PLUS;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("MAXLOG") == 0 ) ret = MAXLOG;

  return(ret);
}

enum pp_y_plus_keps PrePostNN::ReadPPYPlusKEPSFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_y_plus_keps ret = INDEFY_PLUS_KEPS;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("MAXLOG") == 0 ) ret = MAXLOG_KEPS;

  return(ret);
}

enum pp_re_tau PrePostNN::ReadPPReTauFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_re_tau ret = INDEFRE_TAU;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("MAXN") == 0 ) ret = MAXN_RET;

  return(ret);
}

enum pp_re_tau_keps PrePostNN::ReadPPReTauKEPSFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_re_tau_keps ret = INDEFRE_TAU_KEPS;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("MAXN") == 0 ) ret = MAXN_RET_KEPS;

  return(ret);
}

enum pp_b PrePostNN::ReadPPBFromLine(string buffer,string tag,size_t npos)
{
  string tmp;
  size_t ltag;
  enum pp_b ret = INDEFB;

  ltag = tag.length();
  tmp = buffer.substr(npos+ltag,buffer.length()-ltag);
  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());

  if( tmp.compare("FROB") == 0 ) ret = FROB;

  return(ret);
}

string PrePostNN::ReadStringFromLine(string buffer,string tag,size_t npos)
{
  string ret;
  size_t ltag;

  ltag = tag.length();
  ret = buffer.substr(npos+ltag,buffer.length()-ltag);
  ret = trim(ret);

  return(ret);
}


vector<int> PrePostNN::ReadIndexFromLine(string buffer,string tag,size_t npos)
{
  vector<int> ret;
  string sval;
  int val;
  size_t ltag;

  ltag = tag.length();
  sval = buffer.substr(npos+ltag,buffer.length()-ltag);
  istringstream iss(sval);
  while(!iss.eof()){
    iss >> val;
    ret.push_back(val);
  }

  return(ret);
}

vector<vector<double>> PrePostNN::ReadDataFromSeveralLines(ifstream &f,int nblines)
{
  vector<vector<double>> ret;
  vector<double> tmp;
  string buffer;

  for(int i=0;i<nblines;i++){
    getline(f, buffer);
    tmp = ReadDataFromLine(buffer,"",0);
    ret.push_back(tmp);
  }

  return(ret);
}

void PrePostNN::AllDisplay()
{
  cout << "ALPHA_MAX: " << alpha_max << endl;
  cout << "Y_PLUS_MAX_LOG: " << y_plus_max_log << endl;
  cout << "RE_TAU_MAX: " << Re_t_max << endl;
  cout << "B_SIGMA: " << bsigma << endl;

}

void PrePostNN::AllDisplay_carre()
{
  cout << "Y_PLUS_MAX_LOG: " << y_plus_max_log << endl;
  cout << "RE_TAU_MAX: " << Re_t_max << endl;
  cout << "NUM_T: " << num_t << endl;
  display("T_FN:",tfn);
  cout << "B_SIGMA: " << bsigma << endl;

}
