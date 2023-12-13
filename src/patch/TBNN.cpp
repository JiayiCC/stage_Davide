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
//  _model_uploaded = std::unique_ptr(new fdeep::model());
  _model_file = keras_model_file;
  _model_uploaded = std::make_unique<fdeep::model>(fdeep::load_model( _model_file ));
  _pp_alpha = 0.;
  _pp_y_plus = 0.;
  _pp_Re_t = 0.;
  //_normf1 = 0.;
}

TBNN::~TBNN()
{
  delete(_ppNN);
}

//fdeep::model*  TBNN::upload_model ()
//{
//	return &fdeep::load_model( _model_file );
//}


vector<double> TBNN::predict(double alpha, double y_plus, double Re_t)
{
  process_alpha(alpha);
  process_y_plus(y_plus);
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
  switch(_ppNN->get_ppy_plus())
  {
  case MAXLOG:

    if( _ppNN->get_y_plus_max_log() > 0 )
    	_pp_y_plus = log10(y_plus) / _ppNN->get_y_plus_max_log();
    else
    	_pp_y_plus = log10(y_plus);
    break;

  default:
    // on lance une exception
    cerr << "Mauvaise methode de pre traitement des y_plus" << endl;
    break;
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
//void TBNN::process_lambda(vector<double> lambda)
//{
//  vector<double> lc;                 // lambda centre
//  vector<double> lcr;                // lambda centre reduit
//  size_t nbl = lambda.size();  // nombre d'invariants lambda
//
//  lc.resize(nbl);
//  lcr.resize(nbl);
//  _plambda.resize(nbl);
//
//  switch(_ppNN->get_ppl())
//  {
//  case LNORM:
//    // centrage des lambda
//    if( _ppNN->get_lmean().size() == nbl )
//      for(unsigned int i=0;i<nbl;i++)
//        lc[i] = lambda[i] - _ppNN->get_lmean()[i];
//    else
//      for(unsigned int i=0;i<nbl;i++)
//        lc[i] = lambda[i];
//    // reduction des lambda
//    for(unsigned int i=0;i<nbl;i++)
//      _plambda[i] = lc[i] / _ppNN->get_lmax()[i];
//    break;
//  case LU:
//    // puissance alpha
//    if( _ppNN->get_alpha().size() == nbl )
//      for(unsigned int i=0;i<nbl;i++)
//        lc[i] = sgn(lambda[i]) * pow(std::fabs(lambda[i]),_ppNN->get_alpha()[i]);
//    else
//      for(unsigned int i=0;i<nbl;i++)
//        lc[i] = lambda[i];
//    // centrage des lambda
//    if( _ppNN->get_lmean().size() == nbl )
//      for(unsigned int i=0;i<nbl;i++)
//        lc[i] = lc[i] - _ppNN->get_lmean()[i];
//    // reduction des lambda
//    for(unsigned int i=0;i<nbl;i++)
//      lcr[i] = lc[i] / _ppNN->get_lmax()[i];
//    // multiplication par transpose(au)
//    for(unsigned int i=0;i<nbl;i++){
//      _plambda[i] = 0.;
//      for(unsigned int j=0;j<nbl;j++)
//	_plambda[i] += _ppNN->get_lambda_au()[j][i] * lcr[j];
//    }
//    break;
//  case LUS:
//    // puissance alpha
//    if( _ppNN->get_alpha().size() == nbl )
//      for(unsigned int i=0;i<nbl;i++)
//        lc[i] = sgn(lambda[i]) * pow(std::fabs(lambda[i]),_ppNN->get_alpha()[i]);
//    else
//      for(unsigned int i=0;i<nbl;i++)
//        lc[i] = lambda[i];
//    // centrage des lambda
//    if( _ppNN->get_lmean().size() == nbl )
//      for(unsigned int i=0;i<nbl;i++)
//        lc[i] = lc[i] - _ppNN->get_lmean()[i];
//    // reduction des lambda
//    for(unsigned int i=0;i<nbl;i++)
//      lcr[i] = lc[i] / _ppNN->get_lmax()[i];
//    // multiplication par transpose(as)
//    for(unsigned int i=0;i<nbl;i++){
//      _plambda[i] = 0.;
//      for(unsigned int j=0;j<nbl;j++)
//	_plambda[i] += _ppNN->get_lambda_as()[j][i] * lcr[j];
//    }
//    break;
//  default:
//    // on lance une exception
//    cerr << "Mauvaise methode de pre traitement des lambda" << endl;
//    break;
//  }
//}

//void TBNN::process_T(vector<vector<double>> T)
//{
//  size_t nbt = T.size();    // nombre de tenseurs T
//  size_t nbb = T[0].size(); // taille de chacun des tenseurs T
//
//  _pT.resize(nbt);
//  for(unsigned i=0;i<nbt;i++)
//    _pT[i].resize(nbb);
//
//  // le premier tenseur T0 reste inchange
//  for(unsigned int j=0;j<nbb;j++)
//    _pT[0][j] = T[0][j];
//
//  // pre process de T
//  switch(_ppNN->get_ppt())
//  {
//  case TF:
//    // on calcule la norme de Frobenius de chaque tenseur
//    for(unsigned int i=1;i<nbt;i++){
//      double normf = 0.;
//      for(unsigned int j=0;j<nbb;j++)
//	normf += T[i][j] * T[i][j];
//      normf = sqrt(normf);
//      if(i==1) _normf1 = normf;
//      for(unsigned int j=0;j<nbb;j++)
//	_pT[i][j] = T[i][j] / (normf + _ppNN->get_t_thresh());
//    }
//    break;
//  case TR:
//    // on divise le tenseur Ti par la norme globale
//    for(unsigned int i=1;i<nbt;i++)
//      for(unsigned int j=0;j<nbb;j++)
//	_pT[i][j] = T[i][j] / _ppNN->get_tfn()[i-1];
//    break;
//  default:
//    // on lance une exception
//    cerr << "Mauvaise methode de pre traitement des tenseurs T" << endl;
//    break;
//  }
//}

void TBNN::process_b()
{
//  vector<int> iT = _ppNN->get_iT();
//  size_t nbt = iT.size();
//  size_t nbb = _pT[0].size(); //IL FAUT CHANGER T0 GEN!!!!!!!!!!!!!!
//
//  // calcul de _pb a partir de _g et de _pT
//  _pb.resize(nbb);
//  for(unsigned int i=0;i<nbb;i++){
//    _pb[i] = 0;
//    for(unsigned int j=0;j<nbt;j++)
//      _pb[i] += _g[j] * _pT[iT[j]][i];
//  }

	size_t nbb = 6;
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
	else
		cerr << "Bad name of NN file .json (not 'Cas#')" << endl;


  // post process de b
  _b.resize(nbb);
  for(unsigned int i=0;i<nbb;i++)
    _b[i] = _ppNN->get_bsigma() * _pb[i];
}

void TBNN::applyNN()
{
//  vector<int> il = _ppNN->get_ilambda();
//  vector<int> iT = _ppNN->get_iT();
//  Tensor in((int)il.size());
//  Tensor out;
//  size_t nbt = iT.size();
//
//  _g.resize(nbt);
//
//  // l'entree du reseau est determinee par les indices des lambdas definis dans le vecteur il
//  for(unsigned int i=0;i<il.size();i++)
//    in.data_[i] = (float)_plambda[il[i]];
//
//  // on fait la prediction a l'aide du reseau de neurones
//  _model.Apply(&in,&out);
//
//  // on stocke les sorties dans _g
//  for(unsigned int i=0;i<nbt;i++) _g[i] = out(i);

	//fdeep::model _model_uploaded = fdeep::load_model(_model_file);
	const auto result = _model_uploaded->predict({ fdeep::tensor(fdeep::tensor_shape(static_cast<std::size_t>(3)), vector<float>{static_cast<float>(_pp_alpha), static_cast<float>(_pp_y_plus), static_cast<float>(_pp_Re_t)}) });

	_g.resize(result[0].to_vector().size());
	for (unsigned int i =0; i < result[0].to_vector().size(); i++)
		_g[i] = result[0].to_vector()[i];
	_g[1] *= -1;

	if (result[0].to_vector().size() != 3)
		cerr << "Bad treatment of model.predict" << endl;


}

double TBNN::get_g1(double b1, double alpha) //double y_elem, std::vector<double> C, std::vector<double> c_mu_DNS
{

	//double ret = _g[1];

//   //post process de g1
//  switch(_ppNN->get_ppt())
//  {
//  case TF:
//    // on multiplie par la norme de T1
//    ret *= _normf1 + _ppNN->get_t_thresh();
//    break;
//  case TR:
//    // on multiplie par la norme globale
//    ret *= _ppNN->get_tfn()[0];
//    break;
//  default:
//    // on lance une exception
//    cerr << "Mauvaise methode de pre traitement des tenseurs T" << endl;
//    break;
//  }
	//ret*= _ppNN->get_bsigma()/ _ppNN->get_alpha_max(); //circa 1/114


//	double c_mu= 0.0;
//
//	int line_inf = -1;
//
//	      for (unsigned int i = 0; i < C.size(); i++)
//	        {
//	          if (y_elem < C[i])
//	            {
//	              line_inf = i;
//	              break;
//	            }
//	        }
//
//	      // Altrimenti, imposta alpha_true al valore corrispondente da A.txt
//	      if (line_inf == 0)
//	        c_mu =  c_mu_DNS[line_inf];
//
//	      else if (line_inf == -1) // it means it is the highest cell
//	        c_mu = c_mu_DNS[C.size()-1];
//
//	      else
//	        c_mu = 0.5 * (c_mu_DNS[line_inf-1] + c_mu_DNS[line_inf]);

	// ret = 2 * b1/ alpha; GIUSTOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  return 2 * b1/ alpha; //-c_mu
}
