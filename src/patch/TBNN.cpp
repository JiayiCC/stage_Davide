// TRUST_NO_INDENT
#include <TBNN.h>
#include <iostream>


#include <fdeep/fdeep.hpp>
#include <stdio.h>


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
  //cerr << "Processing in square duct flow way" << endl;
  //cerr << "lambda is " << lambda[0] << lambda[1] << lambda[2] << lambda[3] << lambda[4] << endl;
  //cerr << "y_plus is " << y_plus << endl;
  //cerr << "z_plus is " << z_plus << endl;
  //cerr << "Re_t is " << Re_t << endl;
  process_lambda(lambda);
  process_T(T);
  process_y_plus(y_plus);
  process_z_plus(z_plus);
  process_Re_t(Re_t);
  applyNN();
  process_b();
  return(_b);
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

void TBNN::process_y_plus(double y_plus)
{
	if (is_canal_plan_){
	    if( _ppNN->get_y_plus_max_log() > 0 )
	    	_pp_y_plus = log10(y_plus) / _ppNN->get_y_plus_max_log();
	    else
	    	_pp_y_plus = log10(y_plus);
	}

	else if (is_canal_carre_){
		_pp_y_plus = sgn(y_plus) * log( 1 + abs(y_plus) );
	}
	else {
		cerr << "Mauvaise methode de pre traitement des y_plus" << endl;
	}

}


void TBNN::process_z_plus(double z_plus)
{
	if (is_canal_carre_){
		_pp_z_plus = sgn(z_plus) * log( 1 + abs(z_plus) );
	}
	else {
		cerr << "Mauvaise methode de pre traitement des z_plus" << endl;
	}

}

void TBNN::process_Re_t(double Re_t)
{
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
}
//
void TBNN::process_lambda(vector<double> lambda)
{
  vector<double> lc;                 // lambda centre
  vector<double> lcr;                // lambda centre reduit
  size_t nbl = lambda.size();  // nombre d'invariants lambda

  lc.resize(nbl);
  lcr.resize(nbl);
  _plambda.resize(nbl);


  for(unsigned int i=0;i<nbl;i++)
  	  _plambda[i] = sgn(lambda[i]) * log( 1 + abs(lambda[i]) );
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

		// construct input vector
		vector<float> input_vector(_plambda.begin(),_plambda.end());
	    input_vector.push_back(static_cast<float>(_pp_y_plus));
	    input_vector.push_back(static_cast<float>(_pp_z_plus));
	    input_vector.push_back(static_cast<float>(_pp_Re_t));

	    // resize output vector
		size_t nbt = _ppNN->get_numT();
		_g.resize(nbt);

		// on fait la prediction a l'aide du reseau de neurones
		const auto result_carre = _model_uploaded->predict({ fdeep::tensor(fdeep::tensor_shape(static_cast<std::size_t>(8)), input_vector) });
		//std::cout << fdeep::show_tensors(result_carre) << std::endl;
		// on stocke les sorties dans _g
		_g.resize(result_carre[0].to_vector().size());
		for (unsigned int i =0; i < result_carre[0].to_vector().size(); i++)
			_g[i] = result_carre[0].to_vector()[i];

	    // Check the size of the results
	    if (result_carre[0].to_vector().size() != 6)
	    {
	        cerr << "Bad treatment of model.predict" << endl;
	    }
	}

}

double TBNN::get_g1(double b1, double alpha) //double y_elem, std::vector<double> C, std::vector<double> c_mu_DNS
{

	double ret = _g[1];
	ret = 2 * b1/ alpha;
	return ret; //-c_mu
}

double TBNN::get_g1_carre()
{
	//cerr << "We get g1 carre" << endl;
	double ret = _g[1];
	ret *= _ppNN->get_bsigma() / _ppNN->get_tfn()[1];
  return ret; //-c_mu
}
