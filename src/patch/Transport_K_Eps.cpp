/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
//////////////////////////////////////////////////////////////////////////////
//
// File:        Transport_K_Eps.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////

#include <Transport_K_Eps.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Les_Pb_Turb.h>
#include <Param.h>
#include <EChaine.h>
#include <Fluide_Quasi_Compressible.h>
#include <TRUSTTrav.h>
#include <Debog.h>
#include <Discretisation_base.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Domaine_VDF.h>
#include <Champ_Face_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Nom.h>
#include <Dirichlet.h>
#include <Dirichlet_homogene.h>

Implemente_instanciable(Transport_K_Eps,"Transport_K_Eps",Transport_K_Eps_base);

/*! @brief Imprime le type de l'equation sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Transport_K_Eps::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

void Transport_K_Eps::set_param(Param& param)
{
  Transport_K_Eps_base::set_param(param);
  param.ajouter("with_nu",&with_nu_);
  param.dictionnaire("no",0);
  param.dictionnaire("yes",1);
  param.ajouter_non_std("recompute_from_NN",(this));
  param.ajouter_non_std("recompute_from_DNS",(this));
//  param.ajouter("recompute_from_DNS",&KEps_champ);
}

/*! @brief Lit les specifications d'une equation de transport K-epsilon a partir d'un flot d'entree.
 *
 *     Controle dynamique du type du terme source.
 *
 * @param (Entree& s) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Transport_K_Eps::readOn(Entree& s )
{
  with_nu_=0;
  // Lecture des attributs de l'equation
  Transport_K_Eps_base::readOn(s);

  // Ajout automatique du terme source si pas instancie dans le jeu de donnees
  if (les_sources.est_vide())
    {
      Source t;
      Source& so=les_sources.add(t);
      const Probleme_base& pb = probleme();
      Cerr << "Construction and typing for the source term of the Transport_K_Eps equation." << finl;
      Nom pbb = probleme().que_suis_je();
      if (sub_type(Pb_Hydraulique_Turbulent,pb) || milieu().que_suis_je()=="Fluide_Quasi_Compressible")
        {
          Nom typ = "Source_Transport_K_Eps";
          Cerr << "TYPAGE DES SOURCES : this " << *this << finl;
          so.typer(typ,*this);
        }
      else if (pbb.contient("ALE"))
        {
          Nom typ = "Source_Transport_K_Eps";
          Cerr << "TYPAGE DES SOURCES : this " << *this << finl;
          so.typer(typ,*this);
        }
      else if (sub_type(Pb_Thermohydraulique_Turbulent,pb))
        {
          Nom typ = "Source_Transport_K_Eps_anisotherme";
          so.typer(typ,*this);
        }
      else if (sub_type(Pb_Hydraulique_Concentration_Turbulent,pb))
        {
          Nom typ = "Source_Transport_K_Eps_aniso_concen";
          so.typer(typ,*this);
        }
      else if ( (sub_type(Pb_Thermohydraulique_Concentration_Turbulent,pb)) ) //|| (sub_type(Probleme_Combustion,pb)) )
        {
          Nom typ = "Source_Transport_K_Eps_aniso_therm_concen";
          so.typer(typ,*this);
        }
      so->associer_eqn(*this);
    }
  return s;
}


int Transport_K_Eps::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  if (mot=="diffusion")
    {
      Cerr << "Reading and typing of the diffusion operator : " << finl;

      if (with_nu_==0)
        {
          Cerr<<" On associe le champ de diffusion nul afin de faire comme avant !!!!!! " <<finl;
          EChaine tt("Champ_Uniforme 1 0");
          tt>> Champ_don_nul_;
          terme_diffusif.associer_diffusivite(Champ_don_nul_);
        }
      else
        {
          const Fluide_base& fluide_inc = ref_cast(Fluide_base,le_fluide.valeur());
          if (sub_type(Fluide_Quasi_Compressible,fluide_inc))
            terme_diffusif.associer_diffusivite(fluide_inc.viscosite_dynamique());
          else
            terme_diffusif.associer_diffusivite(fluide_inc.viscosite_cinematique());
        }
      is >> terme_diffusif;
      return 1;
    }
  else if (mot=="convection")
    {
      Cerr << "Reading and typing of the convection operator : " << finl;
      const Champ_Inc& vitesse_transportante = probleme().equation(0).inconnue();
      associer_vitesse(vitesse_transportante);
      terme_convectif.associer_vitesse(vitesse_transportante);
      is >> terme_convectif;
      return 1;
    }
  else if ((mot=="ecrire_fichier_xyz_valeur") || (mot=="ecrire_fichier_xyz_valeur_bin"))
    {
      Cerr << mot << " is not understood by " << que_suis_je() << finl;
      Cerr << "Use this keyword in the Navier Stokes equation, not in KEps equation, please." << finl;
      exit();
    }
  else if (mot=="recompute_from_NN")
    {
      Motcle mot2;
      is_recompute_from_NN_ = true;
      is >> mot2;
      if (mot2 != accolade_ouverte)
        {
          Cerr << "On attendait { pour commencer a lire les constantes de recompute_from_NN" << finl;
          exit();
        }

      Motcles les_mots(3);
      {
        les_mots[0] = "nom_fichiers";
        les_mots[1] = "canal_plan";
        les_mots[2] = "canal_carre";
      }
      is >> mot2;

      while (mot2 != accolade_fermee)
        {
          int rang=les_mots.search(mot2);
          switch(rang)
            {
            case 0 : // "nom_fichier"
              {
                is >> nn_casename_k;
                cout << "The read case name for k is " << nn_casename_k << endl;
                is >> nn_casename_eps;
                cout << "The read case name for eps is " << nn_casename_eps << endl;
                readNN();
                break;
              }
            case 1 : // "canal_plan"
              {
                tbnn_k->canal_plan(true);
                tbnn_eps->canal_plan(true);
                break;
              }
            case 2 : // "canal_carre"
              {
                tbnn_k->canal_carre(true);
                tbnn_eps->canal_carre(true);
                break;
              }
            default :
              {
                Cerr << "On ne comprend pas le mot cle : " << mot2 << " dans recompute_from_NN" << finl;
                exit();
              }
            }
          is >> mot2;
        }
    }
  else if (mot=="recompute_from_DNS")
    {
      Motcle mot3;
      is_recompute_from_DNS_ = true;
      is >> mot3;
      if (mot3 != accolade_ouverte)
        {
          Cerr << "On attendait { pour commencer a lire les constantes de recompute_from_DNS" << finl;
          exit();
        }
      Motcles les_mots(1);
      {
        les_mots[0] = "k_eps";
      }

      is >> mot3;

      while (mot3 != accolade_fermee)
        {
          int rang=les_mots.search(mot3);
          switch(rang)
            {
            case 0 : // "nom_fichier"
              {
                is >> KEps_champ;
                break;
              }
            default :
              {
                Cerr << "On ne comprend pas le mot cle : " << mot3 << " dans recompute_from_DNS" << finl;
                exit();
              }
            }
          is >> mot3;
        }

    }
  else
    return Transport_K_Eps_base::lire_motcle_non_standard(mot,is);
  return 1;
}

void Transport_K_Eps::readNN()
{
  string path_NN = string(getenv("project_directory")) + "/share/reseaux_neurones/";
  string model_NN_file = path_NN + string(nn_casename_k) + ".json";
  string ppp_NN_file = path_NN + string(nn_casename_k) + ".ppp";

  cout << "Chargement du reseau de neurones: " + model_NN_file << endl;
  tbnn_k = new TBNN(model_NN_file,ppp_NN_file);  //capire perchè è una giungla

  model_NN_file = path_NN + string(nn_casename_eps) + ".json";
  ppp_NN_file = path_NN + string(nn_casename_eps) + ".ppp";

  cout << "Chargement du reseau de neurones: " + model_NN_file << endl;
  tbnn_eps = new TBNN(model_NN_file,ppp_NN_file);  //capire perchè è una giungla

}

/*! @brief Associe un modele de turbulence K-epsilon a l'equation.
 *
 * @param (Modele_turbulence_hyd_K_Eps& modele) le modele de turbulence K-epsilon a asoocier a l'equation
 */
void Transport_K_Eps::associer_modele_turbulence(const Mod_turb_hyd_RANS& modele)
{
  const Equation_base& eqn_hydr = modele.equation();
  associer(eqn_hydr);
  associer_milieu_base(eqn_hydr.milieu());
  associer_vitesse(eqn_hydr.inconnue());
  mon_modele = ref_cast(Modele_turbulence_hyd_K_Eps,modele);
  discretiser();
}
//
//int Transport_K_Eps::reprendre(Entree& fich)
//{
//  Equation_base::reprendre(fich);
//  if (is_recompute_from_NN_)
//    {
//      cout << "The value of is_recompute_from_NN is " << is_recompute_from_NN_ << endl;
//      if (tbnn_k->is_canal_plan_)
//        {
//          Calcul_keps_NL_TBNN();
//        }
//      if (tbnn_k->is_canal_carre_)
//        {
//          Modele_turbulence_hyd_K_Eps modele_K_Eps = ref_cast(Modele_turbulence_hyd_K_Eps,mon_modele.valeur());
//          modele_K_Eps.Calcul_RSLambda();
//          Calcul_keps_NL_TBNN_carre();
//        }
//    }
//  return 1;
//}

int Transport_K_Eps::preparer_calcul()
{
  Equation_base::preparer_calcul();
  if (is_recompute_from_NN_)
    {
      cout << "The value of is_recompute_from_NN is " << is_recompute_from_NN_ << endl;
      if (tbnn_k->is_canal_plan_)
        {
          cout << "Preparing calculation with TKENN and EPSNN" << endl;
          Calcul_keps_NN();
          cout << "Calculation prepared with TKENN and EPSNN" << endl;
        }
      if (tbnn_k->is_canal_carre_)
        {
          cout << "Preparing calculation with TKENN and EPSNN" << endl;
          Calcul_keps_NN_carre();
          cout << "Calculation prepared with TKENN and EPSNN" << endl;
        }
    }
  if (KEps_champ.non_nul())
    {
      cout << "Preparing calculation with k and eps from DNS" << endl;
      const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF,le_dom_dis.valeur().valeur());
      const IntTab& les_elems=domaine_VDF.domaine().les_elems();
      int nb_som_elem=les_elems.dimension(1);
      DoubleTab& K_eps = le_champ_K_Eps.valeurs();
      int nb_elem_tot=domaine_VDF.nb_elem_tot();
      for (int ele=0; ele<nb_elem_tot; ele++)
        {
          for (int s=0; s<nb_som_elem; s++)
            {
              int sglob=les_elems(ele,s);
              K_eps(ele,0)+=KEps_champ->valeurs()(sglob,0);
              K_eps(ele,1)+=KEps_champ->valeurs()(sglob,1);
            }
        }
      double inv_nb_som_elem=1./(nb_som_elem);
      K_eps*=inv_nb_som_elem;
      cout << "Calculation prepared with k and eps from DNS" << endl;
    }
  return 1;
}


int Transport_K_Eps::Calcul_keps_NN()
{

  const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF,le_dom_dis.valeur().valeur());
  const Domaine_Cl_VDF& domaine_Cl_VDF = ref_cast(Domaine_Cl_VDF,le_dom_Cl_dis.valeur());

  vector<double> k_ec, eps;
  double alpha, y_plus, Re_t;
  double y_maille_paroi =-1, y_elem;
  //  vector<vector<double>> T;
  DoubleTab& K_eps = le_champ_K_Eps.valeurs();
  double u_t;

  double dudy_paroi ;

  vector<double> y_plus_wall, xpx, xpz, Re_true;
  double x_elem, z_elem;
  int position_base = -1;

  int nb_elem_tot=domaine_VDF.nb_elem_tot();
  Modele_turbulence_hyd_K_Eps modele_K_Eps = ref_cast(Modele_turbulence_hyd_K_Eps,mon_modele.valeur());
  const DoubleTab& vitesse = modele_K_Eps.equation().inconnue().valeurs();
  Champ_Face_VDF& ch = ref_cast(Champ_Face_VDF,modele_K_Eps.equation().inconnue().valeur());
  assert (vitesse.line_size() == 1);
  DoubleTab gij(nb_elem_tot,dimension,dimension, vitesse.line_size());
  ch.calcul_duidxj( vitesse,gij,domaine_Cl_VDF );

  const Fluide_base& fluide_inc = ref_cast(Fluide_base,le_fluide.valeur());
  double nu = fluide_inc.viscosite_cinematique().valeurs()[0];

  y_plus_wall.resize(9);
  Re_true.resize(9);
  xpx.resize(9);
  xpz.resize(9);

  // Boucle sur les bords pour traiter les conditions aux limites
  for (int n_bord=0; n_bord<domaine_VDF.nb_front_Cl(); n_bord++)
    {

      // pour chaque Condition Limite on regarde son type
      // Si face de Dirichlet ou de Symetrie on ne fait rien
      // Si face de Neumann on calcule la contribution au terme source

      const Cond_lim& la_cl = domaine_Cl_VDF.les_conditions_limites(n_bord);

      int elem_paroi;

      if ( (sub_type(Dirichlet,la_cl.valeur()))   ||  (sub_type(Dirichlet_homogene,la_cl.valeur())))
        {

          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();

          for ( int num_face=ndeb; num_face<nfin; num_face++)
            {
              elem_paroi = domaine_VDF.face_voisins(num_face,1);//1 means that I want the elem above the face (which is the bottom paroi so ok)

              xpx[num_face-ndeb] =domaine_VDF.xp(elem_paroi, 0) ; // x=0
              xpz[num_face-ndeb] =domaine_VDF.xp(elem_paroi, 2) ; // z=2

              y_maille_paroi =domaine_VDF.xp(elem_paroi, 1) ;

              dudy_paroi = gij(elem_paroi,0,1,0);

              //u_t = sqrt(fabs(nu* dudy_paroi)+1e-20);
              //TODO: calculate average u_t
              u_t = sqrt(nu* dudy_paroi)+1e-40;
              y_plus_wall[num_face-ndeb] = y_maille_paroi * u_t / nu;
              Re_true [num_face-ndeb] =  u_t / nu;
            }
        }
    }



  for (int elem=0; elem<domaine_VDF.nb_elem(); elem++)
    {

      x_elem =domaine_VDF.xp(elem, 0) ; // x=0
      y_elem =domaine_VDF.xp(elem, 1) ; // y=1
      z_elem =domaine_VDF.xp(elem, 2) ; // z=2

      for (unsigned int i = 0; i < y_plus_wall.size(); i++)
        if ( fabs (x_elem-xpx[i]) <1e-3 and fabs (z_elem-xpz[i]) <1e-3)
          {
            position_base = i;
            break;
          }

      alpha = tbnn_k->compute_alpha(K_eps(elem, 0), K_eps(elem,1), gij(elem,0,1,0) ) + 1e-40;
      y_plus = tbnn_k->compute_y_plus(y_plus_wall[position_base],y_maille_paroi, y_elem);
      Re_t = tbnn_k->compute_Re_t(y_plus_wall[position_base],y_maille_paroi);

      // prediction by neural network
      //Cerr << " Testing PCF TRE_VDF_Face.cpp line 1213 " << finl;
      k_ec = tbnn_k->predict_k(alpha, y_plus, Re_t, K_eps(elem, 0)); //;
      eps = tbnn_eps->predict_eps(alpha, y_plus, Re_t, K_eps(elem, 1)); //, K_eps(elem, 1));

      K_eps(elem, 0) = k_ec[0];
      K_eps(elem, 1) = eps[0];
    }


  return 1;
}

int Transport_K_Eps::Calcul_keps_NN_carre()
{
//  Cout << "calling Calcul_keps_NN_carre " << endl;
  const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF,le_dom_dis.valeur().valeur());
  const Domaine_Cl_VDF& domaine_Cl_VDF = ref_cast(Domaine_Cl_VDF,le_dom_Cl_dis.valeur());

  vector<double> lambda;
  vector<double> k_ec, eps;

  double y_plus, z_plus, Re_t;
//  double y_maille_paroi =-1;
//  double z_maille_paroi=-1;
  vector<vector<double>> T;
  double u_t = 0., sum_surface=0.;
  double dudy_paroi, dudz_paroi ;
  double y_elem, z_elem;
//  int position_base = -1;
//  int position_gauche = -1;
  int ori;

  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const IntVect& orientation = domaine_VDF.orientation();
//  const DoubleTab& xv = domaine_VDF.xv();
  const DoubleTab& xp = domaine_VDF.xp();
  const DoubleVect& face_surfaces = domaine_VDF.face_surfaces();
//  DoubleTab y_plus_wall, xp_paroi;
  DoubleTab& K_eps = le_champ_K_Eps.valeurs();

//  y_plus_wall.copy(xv,Array_base::NOCOPY_NOINIT);
//  y_plus_wall = 0.;

//  xp_paroi.copy(xv,Array_base::NOCOPY_NOINIT);
//  xp_paroi = 0.;

  int nb_elem_tot=domaine_VDF.nb_elem_tot();

  Modele_turbulence_hyd_K_Eps modele_K_Eps = ref_cast(Modele_turbulence_hyd_K_Eps,mon_modele.valeur());
  const DoubleTab& vitesse = modele_K_Eps.equation().inconnue().valeurs();
  Champ_Face_VDF& ch = ref_cast(Champ_Face_VDF,modele_K_Eps.equation().inconnue().valeur());
  assert (vitesse.line_size() == 1);
  DoubleTab gij(nb_elem_tot,dimension,dimension, vitesse.line_size());
  ch.calcul_duidxj( vitesse,gij,domaine_Cl_VDF );

  const Fluide_base& fluide_inc = ref_cast(Fluide_base,le_fluide.valeur());
  double nu = fluide_inc.viscosite_cinematique().valeurs()[0];

  // Boucle sur les bords pour calculer u_tau
  for (int n_bord=0; n_bord<domaine_VDF.domaine().nb_bords(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VDF.les_conditions_limites(n_bord);

      int elem_paroi;
      if ( (sub_type(Dirichlet,la_cl.valeur()))   ||  (sub_type(Dirichlet_homogene,la_cl.valeur())))
        {

          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());


          int ndeb = le_bord.num_premiere_face();
          int nfaces = le_bord.nb_faces_tot();
          int nfin = ndeb + nfaces;

          for (int num_face=ndeb; num_face<nfin; num_face++)
            {

              if (face_voisins(num_face, 0) != -1)
                elem_paroi = face_voisins(num_face, 0);
              else
                elem_paroi = face_voisins(num_face, 1);

              ori = orientation(num_face);

              if (ori == 1)
                {
                  dudy_paroi = gij(elem_paroi,0,1,0);
                  u_t += sqrt( abs(nu* dudy_paroi))*face_surfaces(num_face);
                  sum_surface += face_surfaces(num_face);
                }

              if ( ori == 2 )
                {
                  dudz_paroi = gij(elem_paroi,0,2,0);
                  u_t += sqrt( abs(nu* dudz_paroi))*face_surfaces(num_face);
                  sum_surface += face_surfaces(num_face);
                }
            }
        }
    }

  u_t = mp_sum(u_t);
  u_t = u_t / mp_sum(sum_surface);
//  cout << "utau=" << u_t << endl;
  Re_t = u_t/nu;
//  cout << "Re_t=" << Re_t << endl;
  //  init of lambda and T arrays
  lambda.resize(5);
  modele_K_Eps.Calcul_RSLambda();
  DoubleTab& lambda_1_etoile = modele_K_Eps.lambda_1_etoile();
  DoubleTab& lambda_2_etoile = modele_K_Eps.lambda_2_etoile();
  DoubleTab& lambda_3_etoile = modele_K_Eps.lambda_3_etoile();
  DoubleTab& lambda_4_etoile = modele_K_Eps.lambda_4_etoile();
  DoubleTab& lambda_5_etoile = modele_K_Eps.lambda_5_etoile();

  //  static bool executed = false;  // Static variable to track execution
  for (int elem=0; elem<domaine_VDF.nb_elem(); elem++)
    {

      lambda[0] = lambda_1_etoile(elem);
      lambda[1] = lambda_2_etoile(elem);
      lambda[2] = lambda_3_etoile(elem);
      lambda[3] = lambda_4_etoile(elem);
      lambda[4] = lambda_5_etoile(elem);

//      x_elem = xp(elem, 0) ; // x=0
      y_elem = xp(elem, 1) ; // y=1
      z_elem = xp(elem, 2) ; // z=2

      y_plus = (y_elem + 1)*u_t/nu;  //xp_paroi(position_base,1);
      z_plus =  (z_elem + 1)*u_t/nu;

      eps = tbnn_eps->predict_eps_carre(lambda, y_plus, z_plus, Re_t, K_eps(elem, 1));
      k_ec = tbnn_k->predict_k_carre(lambda,  y_plus, z_plus, Re_t, K_eps(elem, 0));

      K_eps(elem, 0) = k_ec[0];
      K_eps(elem, 1) = eps[0];
    }
//  cout << "Re_t=" << Re_t << endl;
  return 1;
}

//int Transport_K_Eps::controler_K_Eps()
//{
//  DoubleTab& K_Eps = le_champ_K_Eps.valeurs();
//  int size=K_Eps.dimension(0);
//  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF,domaine_dis().valeur());
//  double LeEPS_MIN = modele_turbulence().get_LeEPS_MIN();
//  double LeEPS_MAX = modele_turbulence().get_LeEPS_MAX();
//  double LeK_MIN = modele_turbulence().get_LeK_MIN();
//  const IntTab& face_voisins = domaine_vf.face_voisins();
//  const IntTab& elem_faces = domaine_vf.elem_faces();
//  Nom position;
//  // PL on ne fixe au seuil minimum que si negatifs
//  // car la loi de paroi peut fixer a des valeurs tres petites
//  // et le rapport K*K/eps est coherent
//  // Changement: 13/12/07: en cas de valeurs negatives pour k OU eps
//  // on fixe k ET eps a une valeur moyenne des 2 elements voisins
//
//  for (int n=0; n<size; n++)
//    {
//      double& k   = K_Eps(n,0);
//      double& eps = K_Eps(n,1);
//      if (k < 0 || eps < 0)
//        {
//          k = 0;
//          eps = 0;
//          int nk = 0;
//          int neps = 0;
//          int nb_faces_elem = elem_faces.line_size();
//          if (size==face_voisins.dimension(0))
//            {
//              // K-Eps on faces (eg:VEF)
//              for (int i=0; i<2; i++)
//                {
//                  int elem = face_voisins(n,i);
//                  if (elem!=-1)
//                    for (int j=0; j<nb_faces_elem; j++)
//                      if (j != n)
//                        {
//                          double& k_face = K_Eps(elem_faces(elem,j),0);
//                          if (k_face > LeK_MIN)
//                            {
//                              k += k_face;
//                              nk++;
//                            }
//                          double& e_face = K_Eps(elem_faces(elem,j),1);
//                          if (e_face > LeEPS_MIN)
//                            {
//                              eps += e_face;
//                              neps++;
//                            }
//                        }
//                }
//            }
//          if (nk!=0) k /= nk;
//          else k = LeK_MIN;
//          if (neps!=0) eps /= neps;
//          else eps = LeEPS_MIN;
//        }
//      else if (eps > LeEPS_MAX)
//        {
//          eps = LeEPS_MAX;
//        }
//      else if (eps < LeEPS_MIN || k < LeK_MIN )
//        {
//          Cerr << "eps and k < MIN, but we do nothing." << endl;
//        }
//    }
//
//
//  return 1;
//}

/*! @brief Renvoie le nombre d'operateurs de l'equation.
 *
 * Ici 2.
 *
 * @return (int) le nombre d'operateurs de l'equation
 */
int Transport_K_Eps::nombre_d_operateurs() const
{
  return 2;
}

/*! @brief Renvoie l'operateur specifie par son index: renvoie terme_diffusif si i = 0
 *
 *      renvoie terme_convectif si i = 1
 *      exit si i>1
 *     (version const)
 *
 * @param (int i) l'index de l'operateur a renvoyer
 * @return (Operateur&) l'operateur specifie
 * @throws l'equation n'a pas plus de 2 operateurs
 */
const Operateur& Transport_K_Eps::operateur(int i) const
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    case 1:
      return terme_convectif;
    default :
      Cerr << "Error for Transport_K_Eps::operateur("<<i<<") !! " << finl;
      Cerr << "Transport_K_Eps has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // Pour les compilos!!
  return terme_diffusif;
}

/*! @brief Renvoie l'operateur specifie par son index: renvoie terme_diffusif si i = 0
 *
 *      renvoie terme_convectif si i = 1
 *      exit si i>1
 *
 * @param (int i) l'index de l'operateur a renvoyer
 * @return (Operateur&) l'operateur specifie
 * @throws l'equation n'a pas plus de 2 operateurs
 */
Operateur& Transport_K_Eps::operateur(int i)
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    case 1:
      return terme_convectif;
    default :
      Cerr << "Error for Transport_K_Eps::operateur("<<i<<") !! " << finl;
      Cerr << "Transport_K_Eps has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // Pour les compilos!!
  return terme_diffusif;
}


/*! @brief Renvoie le nom du domaine d'application de l'equation.
 *
 * Ici "Transport_Keps".
 *
 * @return (Motcle&) le nom du domaine d'application de l'equation
 */
const Motcle& Transport_K_Eps::domaine_application() const
{
  static Motcle domaine = "Transport_Keps";
  return domaine;
}

DoubleTab& Transport_K_Eps::corriger_derivee_impl(DoubleTab& d)
{
  Nom pbb = probleme().que_suis_je();
  if (pbb.contient("ALE")) corriger_derivee_impl_ALE(d);

  const Turbulence_paroi_base& loi_paroi=modele_turbulence().loi_paroi().valeur();
  loi_paroi.corriger_derivee_impl(d);
  return Transport_K_Eps_base::corriger_derivee_impl(d);
}

