// TRUST_NO_INDENT
#include <TBNN.h>
#include <iostream>


#include <fdeep/fdeep.hpp>
#include <stdio.h>
#include <chrono>
using namespace std::chrono;

int sgn(double v) {
  return (v > 0) - (v < 0);
}

TBNN::TBNN(string keras_model_file,string preproc_file)
{
  _ppNN = new PrePostNN(preproc_file); //dovrebbe essere ok
  //_model_uploaded = std::unique_ptr(new fdeep::model());
  _model_file = keras_model_file;
  _model_uploaded = std::make_unique<fdeep::model>(fdeep::load_model( _model_file, false ));
  _pp_alpha = 0.;
  _pp_y_plus = 0.;
  _pp_Re_t = 0.;
  canal_plan( false );
  canal_carre( false );
  //_normf1 = 0.;
}

TBNN::~TBNN()
{
  delete(_ppNN);
}

vector<double> TBNN::predict(double alpha, double y_plus, double Re_t)
{
  //cerr << "Pre-processing in plane channel flow way" << endl;
  process_alpha(alpha);
  process_y_plus(y_plus);
  process_Re_t(Re_t);
  applyNN();
  process_b();
  return(_b);
}

vector<double> TBNN::predict_carre(vector<double> lambda, vector<vector<double>> T, double y_plus, double z_plus, double Re_t)
{
//  cout << "calling NN " << endl;
  process_lambda(lambda);
  process_T(T);
  process_y_plus(y_plus);
  process_z_plus(z_plus);
  process_Re_t(Re_t);
  applyNN();
  process_b();
  return(_b);
}

vector<double> TBNN::predict_k(double alpha, double y_plus, double Re_t, double tke_input)
{
  //cerr << "Pre-processing in plane channel flow way" << endl;
  process_alpha_keps(alpha);
  process_y_plus_keps(y_plus);
  process_Re_t_keps(Re_t);
  process_k_input(tke_input);
  applyNN_k();
  process_k_output();
  return(_k);
}

vector<double> TBNN::predict_k_carre(vector<double> lambda, double y_plus, double z_plus, double Re_t, double tke_input)
{

  process_lambda(lambda);
  process_y_plus(y_plus);
  process_z_plus(z_plus);
  process_Re_t_keps(Re_t);
  process_k_input(tke_input);
  applyNN_k();
  process_k_output();
  return(_k);
}

vector<double> TBNN::predict_eps(double alpha, double y_plus, double Re_t, double eps_input)
{
  //cerr << "Pre-processing in plane channel flow way" << endl;
  process_alpha_keps(alpha);
  process_y_plus_keps(y_plus);
  process_Re_t_keps(Re_t);
  process_eps_input(eps_input);
  applyNN_eps();
  process_eps_output();
  return(_eps);
}

vector<double> TBNN::predict_eps_carre(vector<double> lambda, double y_plus, double z_plus, double Re_t, double eps_input)
{
//  cout << "calling NN " << endl;
  process_lambda(lambda);
  process_y_plus(y_plus);
  process_z_plus(z_plus);
  process_Re_t_keps(Re_t);
  process_eps_input(eps_input);
  applyNN_eps();
  process_eps_output();
  return(_eps);
}

//void TBNN::output_processed_data()
//{
//    for (int j=0; j<5; j++)
//      {
//        ofstream fileOUT6("pl_" + std::to_string(j+1) + ".dat", ios::app); // open filename.txt in append mod
//        fileOUT6 << _plambda[j] << endl; // append "some stuff" to the end of the file
//        fileOUT6.close(); // close the file
//      }
//}

void TBNN::process_k_input(double tke_input)
{
	  switch(_ppNN->get_ppk_input())
	  {
	  case MAXN_K_INPUT:

	    if( _ppNN->get_k_input_max() > 0 ){
	    	_pp_k = tke_input / _ppNN->get_k_input_max();
//	    	cerr << "k_input_max=" <<  _ppNN->get_k_input_max() << endl;
	    }
	    else
	    	_pp_k = tke_input;
	    break;

	  default:
	    // on lance une exception
	    cerr << "Mauvaise methode de pre traitement des k_input" << endl;
	    break;
	  }
}

void TBNN::process_eps_input(double eps_input)
{

	  switch(_ppNN->get_ppeps_input())
	  {
	  case MAXN_EPS_INPUT:

	    if( _ppNN->get_eps_input_max() > 0 ){
	    	_pp_eps = eps_input / _ppNN->get_eps_input_max();
//	    	cerr << "eps_input_max=" <<  _ppNN->get_eps_input_max() << endl;
	    }
	    else
	    	_pp_eps = eps_input;
	    break;

	  default:
	    // on lance une exception
	    cerr << "Mauvaise methode de pre traitement des eps_input" << endl;
	    break;
	  }
}

void TBNN::process_k_output()
{
	size_t nk = 0;
	nk = _pk.size();
	_k.resize(nk);

	for(unsigned int i=0;i<nk;i++){
		_k[i] = _ppNN->get_k_output_max() * _pk[i];
//		cerr << "k_output_max=" <<  _ppNN->get_k_output_max() << endl;
	}
}

void TBNN::process_eps_output()
{
	size_t neps = 0;
	neps = _peps.size();
	_eps.resize(neps);
	for(unsigned int i=0;i<neps;i++){
		_eps[i] = _ppNN->get_eps_output_max() * _peps[i];
//		cerr << "eps_output_max=" <<  _ppNN->get_eps_output_max() << endl;
	}
}

void TBNN::process_alpha(double alpha)
{
  switch(_ppNN->get_ppalpha())
  {
  case MAXN_A:

    if( _ppNN->get_alpha_max() > 0 )
    	_pp_alpha = alpha / _ppNN->get_alpha_max();
    else
    	_pp_alpha = alpha;
    break;

  default:
    // on lance une exception
    cerr << "Mauvaise methode de pre traitement des alpha" << endl;
    break;
  }
}

void TBNN::process_alpha_keps(double alpha)
{
  switch(_ppNN->get_ppalpha_keps())
  {
  case MAXN_A_KEPS:

    if( _ppNN->get_alpha_keps_max() > 0 )
    	_pp_alpha = alpha / _ppNN->get_alpha_keps_max();
    else
    	_pp_alpha = alpha;
    break;

  default:
    // on lance une exception
    cerr << "Mauvaise methode de pre traitement des alpha_keps" << endl;
    break;
  }
}

void TBNN::process_y_plus(double y_plus)
{
	if (is_canal_plan_){
	    if( _ppNN->get_y_plus_max_log() > 0 )
	    	_pp_y_plus = log10(y_plus) / _ppNN->get_y_plus_max_log();
	    else
	    	_pp_y_plus = log10(y_plus);
	}

	else if (is_canal_carre_){
		_pp_y_plus = sgn(y_plus) * log10( 1 + abs(y_plus) );
	}
	else {
		cerr << "Mauvaise methode de pre traitement des y_plus" << endl;
	}

}

void TBNN::process_y_plus_keps(double y_plus)
{
	if (is_canal_plan_){
	    if( _ppNN->get_y_plus_keps_max_log() > 0 ){
	    	_pp_y_plus = log10(y_plus) / _ppNN->get_y_plus_keps_max_log();}
	    else
	    	_pp_y_plus = log10(y_plus);
	}

	else if (is_canal_carre_){
//		cerr << "y_plus = " << y_plus << endl;
		_pp_y_plus = sgn(y_plus) * log10( 1 + abs(y_plus) );
//		cerr << "_pp_y_plus = " << _pp_y_plus << endl;
	}
	else {
		cerr << "Mauvaise methode de pre traitement des y_plus" << endl;
	}

}


void TBNN::process_z_plus(double z_plus)
{
	if (is_canal_carre_){
//		cerr << "z_plus = " << z_plus << endl;
		_pp_z_plus = sgn(z_plus) * log10( 1 + abs(z_plus) );
//		cerr << "_pp_z_plus = " << _pp_z_plus << endl;
	}
	else {
		cerr << "Mauvaise methode de pre traitement des z_plus" << endl;
	}

}

void TBNN::process_Re_t(double Re_t)
{
  //cerr << "Re_t = " << Re_t << endl;
  switch(_ppNN->get_ppre_tau())
  {
  case MAXN_RET:

    if( _ppNN->get_re_tau_max() > 0 )
    	_pp_Re_t = Re_t / _ppNN->get_re_tau_max();
    else
    	_pp_Re_t = Re_t;
    break;

  default:
    // on lance une exception
    cerr << "Mauvaise methode de pre traitement des Re_t" << endl;
    break;
  }
  //cerr << "_pp_Re_t = " << _pp_Re_t << endl;
}

void TBNN::process_Re_t_keps(double Re_t)
{
  switch(_ppNN->get_ppre_tau_keps())
  {
  case MAXN_RET_KEPS:

    if( _ppNN->get_re_tau_keps_max() > 0 ){
//    	cout << "Re_t= " << Re_t << endl;
    	_pp_Re_t_keps = Re_t / _ppNN->get_re_tau_keps_max();
//    	cout << "re_tau_max " << _ppNN->get_re_tau_keps_max() << endl;
//    	cout << "_pp_Re_t= " << _pp_Re_t_keps << endl;
    }
    else
    	_pp_Re_t_keps = Re_t;
    break;

  default:
    // on lance une exception
    cerr << "Mauvaise methode de pre traitement des Re_t_keps" << endl;
    break;
  }
}

//
void TBNN::process_lambda(vector<double> lambda)
{
  size_t nbl = lambda.size();  // nombre d'invariants lambda
  _plambda.resize(nbl);

  for(unsigned int i=0;i<nbl;i++)
  	  _plambda[i] = sgn(lambda[i]) * log10( 1 + abs(lambda[i]) );
}

void TBNN::process_T(vector<vector<double>> T)
{
  size_t nbt = _ppNN->get_numT(); // nombre de tenseurs T
  size_t nbb = T[0].size(); // taille de chacun des tenseurs T
//
  _pT.resize(nbt);
  for(unsigned i=0;i<nbt;i++)
    _pT[i].resize(nbb);

// pre process de T
  switch(_ppNN->get_ppt())
  {
  	case FROT:
		// on divise le tenseur Ti par la norme globale
		for(unsigned int i=0;i<nbt;i++)
			for(unsigned int j=0;j<nbb;j++)
				_pT[i][j] = T[i][j] / _ppNN->get_tfn()[i];
		break;
	default:
		// on lance une exception
		cerr << "Mauvaise methode de pre traitement des tenseurs T" << endl;
		break;
  }
}

void TBNN::process_b()
{
	size_t nbb = 0;
	size_t nbt = 0;

	if (is_canal_plan_){
		nbb = 6;
		_pb.resize(nbb);

		_pb[1] = 0.5 * _pp_alpha * _g[1];
		_pb[2] = 0.0;
		_pb[4] = 0.0;
		if (_model_file.find("Cas5") != string::npos) {

			_pb[0] = -1.0 / 3.0 *_g[0] - 0.5* _pp_alpha*_pp_alpha*_g[2];
			_pb[3] = 1.0 / 6.0 *_g[0] + 0.5* _pp_alpha*_pp_alpha*_g[2];
			_pb[5] = 1.0 / 6.0 *_g[0];
			//cerr << "Entrato dentro Cas5" << endl;
		}
		else if (_model_file.find("Cas6") != string::npos) {

			_pb[0] = 1.0 / 6.0 *_g[0] - 0.5* _pp_alpha*_pp_alpha*_g[2];
			_pb[3] = -1.0 / 3.0 *_g[0] + 0.5* _pp_alpha*_pp_alpha*_g[2];
			_pb[5] = 1.0 / 6.0 *_g[0];
			//cerr << "Entrato dentro Cas6" << endl;
		}
		else if (_model_file.find("Cas7") != string::npos) {

			_pb[0] = 1.0 / 6.0 *_g[0] - 0.5* _pp_alpha*_pp_alpha*_g[2];
			_pb[3] = 1.0 / 6.0 *_g[0] + 0.5* _pp_alpha*_pp_alpha*_g[2];
			_pb[5] = -1.0 / 3.0 *_g[0];
			//cerr << "Entrato dentro Cas7" << endl;
		}
		else if (_model_file.find("Cas8") != string::npos) {

			_pb[0] = _g[0];
			_pb[3] = _g[2];
			_pb[5] = -(_g[0]+_g[2]);
			//cerr << "Entrato dentro Cas8" << endl;
		}
		else if (_model_file.find("CasLS8") != string::npos) {

			_pb[0] = _g[0];
			_pb[3] = _g[2];
			_pb[5] = -(_g[0]+_g[2]);
			//cerr << "Entrato dentro CasLS8" << endl;
		}
		else if (_model_file.find("Cas9") != string::npos) {

			_pb[0] = _g[0];
			_pb[3] = _g[2];
			_pb[5] = -(_g[0]+_g[2]);
			//cerr << "Entrato dentro Cas9" << endl;
		}
		else if (_model_file.find("Cas10") != string::npos) {

			_pb[0] = _g[0];
			_pb[3] = _g[2];
			_pb[5] = -(_g[0]+_g[2]);

			//cerr << "Entrato dentro Cas10" << endl;
		}
		else if (_model_file.find("Cas11") != string::npos) {

			_pb[0] = -1.0 / 3.0 *_g[0] - 0.5* _pp_alpha*_pp_alpha*_g[2];
			_pb[3] = 1.0 / 6.0 *_g[0] + 0.5* _pp_alpha*_pp_alpha*_g[2];
			_pb[5] = 1.0 / 6.0 *_g[0];
			//cerr << "Entrato dentro Cas11" << endl;
		}
		else if (_model_file.find("Cas12") != string::npos) {

			_pb[0] = _g[0];
			_pb[3] = _g[2];
			_pb[5] = -(_g[0]+_g[2]);

			//cerr << "Entrato dentro Cas12" << endl;
		}
		else if (_model_file.find("Cas13") != string::npos) {

			_pb[0] = -1.0 / 3.0 *_g[0] - 0.5* _pp_alpha*_pp_alpha*_g[2];
			_pb[3] = 1.0 / 6.0 *_g[0] + 0.5* _pp_alpha*_pp_alpha*_g[2];
			_pb[5] = 1.0 / 6.0 *_g[0];
			//cerr << "Entrato dentro Cas13" << endl;
		}
		else
			cerr << "Bad name of NN file .json (not 'Cas#')" << endl;
	}

	if (is_canal_carre_){
		//vector<int> iT = _ppNN->get_iT();
		nbt = _ppNN->get_numT();
		nbb = _pT[0].size();

		// calcul de _pb a partir de _g et de _pT
		_pb.resize(nbb);
		for(unsigned int i=0;i<nbb;i++){
			_pb[i] = 0;
			for(unsigned int j=0;j<nbt;j++)
				_pb[i] += _g[j] * _pT[j][i];
		}
		_pb[5] = -(_pb[0] + _pb[3]);
	}

  // post process de b
  _b.resize(nbb);
  for(unsigned int i=0;i<nbb;i++)
	_b[i] = _ppNN->get_bsigma() * _pb[i];
}

void TBNN::applyNN()
{
	//fdeep::model _model_uploaded = fdeep::load_model(_model_file);

	if (is_canal_plan_){
		//_ppNN->AllDisplay();
		const auto result_plan =  _model_uploaded->predict({ fdeep::tensor(fdeep::tensor_shape(static_cast<std::size_t>(3)), vector<float>{static_cast<float>(_pp_alpha), static_cast<float>(_pp_y_plus), static_cast<float>(_pp_Re_t)}) });

		_g.resize(result_plan[0].to_vector().size());
		for (unsigned int i =0; i < result_plan[0].to_vector().size(); i++)
			_g[i] = result_plan[0].to_vector()[i];
		_g[1] *= -1;

		if (result_plan[0].to_vector().size() != 3)
			cerr << "Bad treatment of model.predict" << endl;
	}

	if (is_canal_carre_){

		//_ppNN->AllDisplay_carre();


		// Sequential predictions
		// construct input vector
		vector<float> input_vector(_plambda.begin(),_plambda.end());
		input_vector.push_back(static_cast<float>(_pp_y_plus));
		input_vector.push_back(static_cast<float>(_pp_z_plus));
		input_vector.push_back(static_cast<float>(_pp_Re_t));
		const auto result_carre = _model_uploaded->predict({ fdeep::tensor(fdeep::tensor_shape(static_cast<std::size_t>(8)), input_vector) });

		_g.resize(result_carre[0].to_vector().size());
		for (unsigned int i =0; i < result_carre[0].to_vector().size(); i++)
			_g[i] = result_carre[0].to_vector()[i];
		_g[1] *= -1;

//	    // Check the size of the results
//	    if (result_carre[0].to_vector().size() != 6)
//	    {
//	        cerr << "Bad treatment of model.predict" << endl;
//	    }

	}

}

void TBNN::applyNN_k()
{
	//fdeep::model _model_uploaded = fdeep::load_model(_model_file);

	if (is_canal_plan_){
		//_ppNN->AllDisplay();

		const auto result_plan =  _model_uploaded->predict({ fdeep::tensor(fdeep::tensor_shape(static_cast<std::size_t>(4)), vector<float>{static_cast<float>(_pp_Re_t), static_cast<float>(_pp_y_plus), static_cast<float>(_pp_alpha), static_cast<float>(_pp_k)}) });

		_pk.resize(result_plan[0].to_vector().size());
		for (unsigned int i =0; i < result_plan[0].to_vector().size(); i++)
			_pk[i] = result_plan[0].to_vector()[i] ;
	}

	if (is_canal_carre_){
		vector<float> input_vector(_plambda.begin(),_plambda.end());
		input_vector.push_back(static_cast<float>(_pp_y_plus));
		input_vector.push_back(static_cast<float>(_pp_z_plus));
		input_vector.push_back(static_cast<float>(_pp_Re_t));
		input_vector.push_back(static_cast<float>(_pp_k));
		const auto result_carre = _model_uploaded->predict({ fdeep::tensor(fdeep::tensor_shape(static_cast<std::size_t>(9)), input_vector) });

		_pk.resize(result_carre[0].to_vector().size());
		for (unsigned int i =0; i < result_carre[0].to_vector().size(); i++)
			_pk[i] = result_carre[0].to_vector()[i] ;

	}

}

void TBNN::applyNN_eps()
{
	//fdeep::model _model_uploaded = fdeep::load_model(_model_file);

	if (is_canal_plan_){
		//_ppNN->AllDisplay();

		const auto result_plan =  _model_uploaded->predict({ fdeep::tensor(fdeep::tensor_shape(static_cast<std::size_t>(4)), vector<float>{static_cast<float>(_pp_Re_t), static_cast<float>(_pp_y_plus), static_cast<float>(_pp_alpha), static_cast<float>(_pp_eps)}) });

		_peps.resize(result_plan[0].to_vector().size());
		for (unsigned int i =0; i < result_plan[0].to_vector().size(); i++)
			_peps[i] = result_plan[0].to_vector()[i] ;
	}

	if (is_canal_carre_){

		vector<float> input_vector(_plambda.begin(),_plambda.end());
		input_vector.push_back(static_cast<float>(_pp_y_plus));
		input_vector.push_back(static_cast<float>(_pp_z_plus));
		input_vector.push_back(static_cast<float>(_pp_Re_t));
		input_vector.push_back(static_cast<float>(_pp_eps));
		const auto result_carre = _model_uploaded->predict({ fdeep::tensor(fdeep::tensor_shape(static_cast<std::size_t>(9)), input_vector) });

		_peps.resize(result_carre[0].to_vector().size());
		for (unsigned int i =0; i < result_carre[0].to_vector().size(); i++)
			_peps[i] = result_carre[0].to_vector()[i] ;

	}

}

double TBNN::get_g1(double b1, double alpha) //double y_elem, std::vector<double> C, std::vector<double> c_mu_DNS
{

	double ret = _g[1];
	ret = 2 * b1/ alpha;
	return ret; //-c_mu
}

vector<double> TBNN::get_g_carre()
{
	//cerr << "We get g1 carre" << endl;
	vector<double> ret = _g;
	for (unsigned int i=0; i<_g.size(); i++)
		ret[i] *= _ppNN->get_bsigma() / _ppNN->get_tfn()[i];
  return ret; //-c_mu
}

double TBNN::compute_alpha( double k, double eps, double dudy)
{
  double alpha;

  alpha = k / (eps + 1e-15) * dudy;

  return alpha;
}


double TBNN::compute_y_plus( double y_plus_wall, double h_maille_paroi, double h_elem)
{
  double y_plus;
  y_plus = abs(y_plus_wall / h_maille_paroi * h_elem);

  return y_plus;

}

double TBNN::compute_Re_t(double y_plus_wall, double h_maille_paroi)
{
  double Re_t;

  Re_t = abs(y_plus_wall / h_maille_paroi * 1.0);

  return Re_t;

}


