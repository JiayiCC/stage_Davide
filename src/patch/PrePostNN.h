// TRUST_NO_INDENT
#include <string>
#include <vector>

#ifndef PrePostNN_included
#define PrePostNN_included

using namespace std;

enum pp_lambda {INDEFL,LNORM,LU,LUS};
enum pp_T {INDEFT,TF,TR,FROT};
enum pp_alpha {INDEFALPHA,MAXN_A};
enum pp_alpha_keps {INDEFALPHA_KEPS,MAXN_A_KEPS};
enum pp_y_plus {INDEFY_PLUS,MAXLOG};
enum pp_y_plus_keps {INDEFY_PLUS_KEPS,MAXLOG_KEPS};
enum pp_re_tau {INDEFRE_TAU,MAXN_RET};
enum pp_re_tau_keps {INDEFRE_TAU_KEPS,MAXN_RET_KEPS};
enum pp_b {INDEFB,FROB};
enum pp_k_input {INDEFK_INPUT,MAXN_K_INPUT};
enum pp_eps_input {INDEFEPS_INPUT,MAXN_EPS_INPUT};
enum pp_k_output {INDEFK_OUTPUT,MAXN_K_OUTPUT};
enum pp_eps_output {INDEFEPS_OUTPUT,MAXN_EPS_OUTPUT};

class PrePostNN
{
public:
  PrePostNN(string filename);
  ~PrePostNN();

  void AllDisplay();
  void AllDisplay_carre();

  vector<double> get_alpha() {return alpha;}
  vector<double> get_lmean() {return lmean;}
  vector<double> get_lmax() {return lmax;}
  vector<double> get_tfn() {return tfn;}
  vector<double> get_t0() {return t0;}
  double get_bsigma() {return bsigma;}
  double get_t_thresh() {return t_thresh;}
  vector<vector<double>> get_lambda_au() {return lambda_au;}
  vector<vector<double>> get_lambda_as() {return lambda_as;}
  enum pp_lambda get_ppl() {return ppl;}
  enum pp_T get_ppt() {return ppt;}
  vector<int> get_ilambda() {return ilambda;}
  vector<int> get_iT() {return iT;}

  double get_alpha_max() {return alpha_max;}
  double get_alpha_keps_max() {return alpha_keps_max;}
  double get_y_plus_max_log() {return y_plus_max_log;}
  double get_y_plus_keps_max_log() {return y_plus_keps_max_log;}
  double get_re_tau_max() {return Re_t_max;}
  double get_re_tau_keps_max() {return Re_t_keps_max;}
  double get_k_input_max() {return k_input_max;}
  double get_eps_input_max() {return eps_input_max;}
  double get_k_output_max() {return k_output_max;}
  double get_eps_output_max() {return eps_output_max;}
  int get_numT() {return num_t;}
  string get_datadir() {return datadir;}
  enum pp_alpha get_ppalpha() {return ppalpha;}
  enum pp_alpha_keps get_ppalpha_keps() {return ppalpha_keps;}
  enum pp_y_plus get_ppy_plus() {return ppy_plus;}
  enum pp_y_plus_keps get_ppy_plus_keps() {return ppy_plus_keps;}
  enum pp_re_tau get_ppre_tau() {return ppre_tau;}
  enum pp_re_tau_keps get_ppre_tau_keps() {return ppre_tau_keps;}
  enum pp_b get_ppb() {return ppb;}
  enum pp_k_input get_ppk_input() {return ppk_input;}
  enum pp_eps_input get_ppeps_input() {return ppeps_input;}
  enum pp_k_output get_ppk_output() {return ppk_output;}
  enum pp_eps_output get_ppeps_output() {return ppeps_output;}

private:
  void display(string tag,vector<double> vec);
  void display(string tag,vector<int> vec);
  void display(string tag,vector<vector<double>> mat);
  string trim(const std::string& str, const std::string& whitespace = " \t");
  vector<double> ReadDataFromLine(string buffer,string tag,size_t npos);
  double ReadOneDataFromLine(string buffer,string tag,size_t npos);
  enum pp_T ReadPPTFromLine(string buffer,string tag,size_t npos);
  enum pp_lambda ReadPPLFromLine(string buffer,string tag,size_t npos);
  vector<int> ReadIndexFromLine(string buffer,string tag,size_t npos);
  vector<vector<double>> ReadDataFromSeveralLines(ifstream &f,int nblines);
  enum pp_alpha ReadPPAlphaFromLine(string buffer,string tag,size_t npos);
  enum pp_alpha_keps ReadPPAlphaKEPSFromLine(string buffer,string tag,size_t npos);
  enum pp_y_plus ReadPPYPlusFromLine(string buffer,string tag,size_t npos);
  enum pp_y_plus_keps ReadPPYPlusKEPSFromLine(string buffer,string tag,size_t npos);
  enum pp_re_tau ReadPPReTauFromLine(string buffer,string tag,size_t npos);
  enum pp_re_tau_keps ReadPPReTauKEPSFromLine(string buffer,string tag,size_t npos);
  enum pp_b ReadPPBFromLine(string buffer,string tag,size_t npos);
  enum pp_k_input ReadPPKinputFromLine(string buffer,string tag,size_t npos);
  enum pp_eps_input ReadPPEPSinputFromLine(string buffer,string tag,size_t npos);
  enum pp_k_output ReadPPKoutputFromLine(string buffer,string tag,size_t npos);
  enum pp_eps_output ReadPPEPSoutputFromLine(string buffer,string tag,size_t npos);
  string ReadStringFromLine(string buffer,string tag,size_t npos);

  vector<double> alpha;
  vector<double> lmean;
  vector<double> lmax;
  vector<double> tfn;
  vector<double> t0;
  double bsigma;
  double t_thresh;
  vector<vector<double>> lambda_au;
  vector<vector<double>> lambda_as;
  enum pp_lambda ppl;
  enum pp_T ppt;
  vector<int> ilambda;
  vector<int> iT;

  double alpha_max;
  double alpha_keps_max;
  double y_plus_max_log;
  double y_plus_keps_max_log;
  double Re_t_max;
  double Re_t_keps_max;
  double k_input_max;
  double eps_input_max;
  double k_output_max;
  double eps_output_max;
  string datadir;
  enum pp_alpha ppalpha;
  enum pp_alpha_keps ppalpha_keps;
  enum pp_y_plus ppy_plus;
  enum pp_y_plus_keps ppy_plus_keps;
  enum pp_re_tau ppre_tau;
  enum pp_re_tau_keps ppre_tau_keps;
  enum pp_b ppb;
  enum pp_k_input ppk_input;
  enum pp_eps_input ppeps_input;
  enum pp_k_output ppk_output;
  enum pp_eps_output ppeps_output;
  int num_t;

};

#endif
