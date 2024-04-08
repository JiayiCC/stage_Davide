// TRUST_NO_INDENT
#include <string>
#include <vector>

using namespace std;

enum pp_lambda {INDEFL,LNORM,LU,LUS};
enum pp_T {INDEFT,TF,TR,FROT};
enum pp_alpha {INDEFALPHA,MAXN_A};
enum pp_y_plus {INDEFY_PLUS,MAXLOG};
enum pp_re_tau {INDEFRE_TAU,MAXN_RET};
enum pp_b {INDEFB,FROB};

class PrePostNN
{
public:
  PrePostNN(string filename);
  ~PrePostNN();

  void AllDisplay();

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
  double get_y_plus_max_log() {return y_plus_max_log;}
  double get_re_tau_max() {return Re_t_max;}
  int get_numT() {return num_t;}
  string get_datadir() {return datadir;}
  enum pp_alpha get_ppalpha() {return ppalpha;}
  enum pp_y_plus get_ppy_plus() {return ppy_plus;}
  enum pp_re_tau get_ppre_tau() {return ppre_tau;}
  enum pp_b get_ppb() {return ppb;}

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
  enum pp_y_plus ReadPPYPlusFromLine(string buffer,string tag,size_t npos);
  enum pp_re_tau ReadPPReTauFromLine(string buffer,string tag,size_t npos);
  enum pp_b ReadPPBFromLine(string buffer,string tag,size_t npos);
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
  double y_plus_max_log;
  double Re_t_max;
  string datadir;
  enum pp_alpha ppalpha;
  enum pp_y_plus ppy_plus;
  enum pp_re_tau ppre_tau;
  enum pp_b ppb;
  
  int num_t;

};
