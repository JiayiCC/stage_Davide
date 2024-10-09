// TRUST_NO_INDENT
#include <keras_model.h>
#include <PrePostNN.h>
#include <memory>


#ifndef TBNN_included
#define TBNN_included

//#include <fdeep/fdeep.hpp>

namespace fdeep
{
	class model;
}

class TBNN
{
public:
  TBNN(string keras_model_file,string preproc_file);
  ~TBNN();

  vector<double> predict(double alpha, double y_plus, double Re_t);
  vector<double> predict_carre(vector<double> lambda, vector<vector<double>> T, double y_plus, double z_plus, double Re_t);
  vector<double> predict_k(double alpha, double y_plus, double Re_t, double tke_input);
  vector<double> predict_k_carre(vector<double> lambda, double y_plus, double z_plus, double Re_t, double tke_input);
  vector<double> predict_eps(double alpha, double y_plus, double Re_t, double eps_input);
  vector<double> predict_eps_carre(vector<double> lambda, double y_plus, double z_plus, double Re_t, double eps_input);
  double get_g1(double b1, double alpha);
  vector<double> get_g_carre();
  vector<double> get_t0() { return _ppNN->get_t0(); };
  inline void canal_plan( bool val );
  inline void canal_carre( bool val);
  void output_processed_data();
  bool is_canal_plan_;
  bool is_canal_carre_;

  double compute_alpha( double k, double eps, double dudy);
  double compute_y_plus( double y_plus_wall, double h_maille_paroi, double h_elem);
  double compute_Re_t(double y_plus_wall, double h_maille_paroi);

private:

// Neural network
  string _model_file;              // Neural network
  std::unique_ptr<fdeep::model> _model_uploaded; // = fdeep::load_model(_model_file);
  vector<double> _g;              // output of the Neural network
  PrePostNN *_ppNN;               // objet pre et post processing
  
// Lambda pre-processing
  vector<double> _plambda;        // pre-processed lambda

  double _pp_alpha;
  double _pp_y_plus;
  double _pp_z_plus;
  double _pp_Re_t;
  double _pp_Re_t_keps;
  double _pp_k;
  double _pp_eps;

// T pre-processing
  double _normf1 = 0;                 // norm de froebenius pour le tenseur T1
  vector<vector<double>> _pT;     // pre-processed T

// b post-processing
  vector<double> _pb;             // result of the prediction not post-processed
  vector<double> _b;              // result of the prediction post-processed

  vector<double> _pk;             // result of the prediction not post-processed
  vector<double> _k;              // result of the prediction post-processed

  vector<double> _peps;             // result of the prediction not post-processed
  vector<double> _eps;              // result of the prediction post-processed

// process methods
  void process_lambda(vector<double> lambda); // calculate the pre-processed lambda
  void process_T(vector<vector<double>> T);   // calculate the pre-processed T
  void applyNN();                             // prédiction du réseau de neurones
  void applyNN_k();
  void applyNN_eps();
  void process_b();				  // calculate the post-processed b
  void process_k_output();
  void process_eps_output();
  void process_alpha(double alpha);
  void process_alpha_keps(double alpha);
  void process_Re_t(double Re_t);
  void process_Re_t_keps(double Re_t);
  void process_y_plus(double y_plus);
  void process_y_plus_keps(double y_plus);
  void process_z_plus(double z_plus);
  void process_k_input(double tke_input);
  void process_eps_input(double eps_input);
  //fdeep::model*  upload_model ();

};

inline void TBNN::canal_plan( bool val )
{
  is_canal_plan_ = val;
}

inline void TBNN::canal_carre( bool val )
{
  is_canal_carre_ = val;
}

#endif


