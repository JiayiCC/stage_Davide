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
// File:        Tenseur_Reynolds_Externe_VDF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Kernel/IA
//
//////////////////////////////////////////////////////////////////////////////

#include <Tenseur_Reynolds_Externe_VDF_Face.h>
#include <Champ_Uniforme.h>
#include <Domaine_Cl_dis.h>
#include <Domaine_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Neumann_sortie_libre.h>
#include <Dirichlet.h>
#include <Dirichlet_homogene.h>
#include <Symetrie.h>
#include <Periodique.h>
#include <Navier_Stokes_Turbulent.h>
#include <Probleme_base.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Modele_Fonc_Bas_Reynolds.h>
#include <TRUSTTrav.h>
#include <Dirichlet_paroi_defilante.h>
#include <Echange_externe_impose.h>
#include <Neumann.h>
#include <Neumann_homogene.h>
#include <Champ_Face_VDF.h>
#include <Transport_K_Eps.h>
#include <TRUST_Ref.h>
#include <Process.h>
#include <Modifier_nut_pour_fluide_dilatable.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Fluide_base.h>

#include <Navier_Stokes_std.h>
#include <iostream>
#include <algorithm>
#include <fstream>
using namespace std;
#include <chrono>
using namespace std::chrono;

Implemente_instanciable_sans_constructeur_ni_destructeur(Tenseur_Reynolds_Externe_VDF_Face,"Tenseur_Reynolds_Externe_VDF_Face",Source_base);

// XD tenseur_Reynolds_externe source_base tenseur_Reynolds_externe 1 Use a neural network to estimate the values of the Reynolds tensor. The structure of the neural networks is stored in a file located in the share/reseaux_neurones directory.
// XD  attr nom_fichier chaine nom_fichier 0 The base name of the file.

Tenseur_Reynolds_Externe_VDF_Face::Tenseur_Reynolds_Externe_VDF_Face()
{

}

Tenseur_Reynolds_Externe_VDF_Face::~Tenseur_Reynolds_Externe_VDF_Face()
{
  delete tbnn;
}


//// printOn
//

Sortie& Tenseur_Reynolds_Externe_VDF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

Entree& Tenseur_Reynolds_Externe_VDF_Face::readOn(Entree& is )
{

  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;

  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "On attendait { pour commencer a lire les constantes de Tenseur_Reynolds_Externe" << finl;
      exit();
    }
  Motcles les_mots(22);
  {
    les_mots[0] = "nom_fichier";
    les_mots[1] = "canal_plan";
    les_mots[2] = "canal_carre";
    les_mots[3] = "0";
    les_mots[4] = "1";
    les_mots[5] = "2";
    les_mots[6] = "3";
    les_mots[7] = "4";
    les_mots[8] = "5";
    les_mots[9] = "6";
    les_mots[10] = "7";
    les_mots[11] = "8";
    les_mots[12] = "9";
    les_mots[13] = "10";
    les_mots[14] = "mettre_a_jour";
    les_mots[15] = "T_stab";
    les_mots[16] = "sous_relaxation_RST";
    les_mots[17] = "sous_relaxation_nut";
    les_mots[18] = "only_g1";
    les_mots[19] = "K_EPS_fichiers";
    les_mots[20] = "k_eps";
    les_mots[21] = "Nisizima";
  }
  is >> motlu;
  //Cerr << "We are here" << finl;
// Cerr << "The input stream is " << is << finl;
  while (motlu != accolade_fermee)
    {
      int rang=les_mots.search(motlu);
      switch(rang)
        {
        case 0 : // "nom_fichier"
          {
            is >> nn_casename;
            cout << "The read case name is " << nn_casename << endl;
            readNN();
            break;
          }
        case 1 : // "canal_plan"
          {
            tbnn->canal_plan(true);
            break;
          }
        case 2 : // "canal_carre"
          {
            tbnn->canal_carre(true);
            break;
          }
        case 3 :
          {
            T_list_.push_back(0);
//            modele_K_Eps_.valeur().T_list().push_back(0);
            break;
          }
        case 4 :
          {
//            cout << "we got 1 for T_list " << endl;
            T_list_.push_back(1);
            break;
          }
        case 5 :
          {
            T_list_.push_back(2);
            break;
          }
        case 6 :
          {
            T_list_.push_back(3);
            break;
          }
        case 7 :
          {
            T_list_.push_back(4);
            break;
          }
        case 8 :
          {
            T_list_.push_back(5);
            break;
          }
        case 9 :
          {
            T_list_.push_back(6);
            break;
          }
        case 10 :
          {
            T_list_.push_back(7);
            break;
          }
        case 11 :
          {
            T_list_.push_back(8);
            break;
          }
        case 12 :
          {
            T_list_.push_back(9);
            break;
          }
        case 13 :
          {
            T_list_.push_back(10);
            break;
          }
        case 14 :
          {
            is >> idt_mettre_a_jour;
            cout << "The external RST source will be updated every " << idt_mettre_a_jour << " time step" << endl;
            break;
          }
        case 15 :
          {
            is >> T_stab;
            cout << "The stability function is 1-exp(-t/" << T_stab << ")" << endl;
            break;
          }
        case 16 :
          {
            is >> relax_RST;
            cout << "The RST relaxation factor is " << relax_RST << endl;
            break;
          }
        case 17 :
          {
            is >> relax_g1;
            cout << "The g1 relaxation factor is " << relax_g1 << endl;
            break;
          }
        case 18 : // "canal_plan"
          {
            only_g1(true);
            break;
          }
        case 19 :
          {
            apply_keps_NN(true);
            is >> nn_casename_k;
            cout << "The read case name for k is " << nn_casename_k << endl;
            is >> nn_casename_eps;
            cout << "The read case name for eps is " << nn_casename_eps << endl;
            readNN_keps();
            break;
          }
        case 20 :
          {
            apply_keps_DNS(true);
            is >> KEps_champ;
            break;
          }
        case 21 :
          {
            apply_Nisizima(true);
            break;
          }
        default :
          {
            Cerr << "On ne comprend pas le mot cle : " << motlu << " dans Tenseur_Reynolds_Externe" << finl;
            exit();
          }
        }

      is >> motlu;
    }
  return is;

}

void Tenseur_Reynolds_Externe_VDF_Face::readNN()
{
  string path_NN = string(getenv("project_directory")) + "/share/reseaux_neurones/";
  string model_NN_file = path_NN + string(nn_casename) + ".json";
  string ppp_NN_file = path_NN + string(nn_casename) + ".ppp";

  cout << "Chargement du reseau de neurones: " + model_NN_file << endl;
  tbnn = new TBNN(model_NN_file,ppp_NN_file);  //capire perchè è una giungla

}

void Tenseur_Reynolds_Externe_VDF_Face::readNN_keps()
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

void Tenseur_Reynolds_Externe_VDF_Face::associer_pb(const Probleme_base& pb)
{
  const Equation_base& eqn = pb.equation(0);
  if  ( !sub_type(Navier_Stokes_Turbulent,eqn) )
    {
      Cerr << "Error TRUST in " << que_suis_je()  << finl;
      Cerr << "Hydraulic equation not found" << finl;
      exit();
    }
  else
    {
      probleme_ = pb;

      eqn_NS_ = ref_cast(Navier_Stokes_Turbulent,eqn);

      const Mod_turb_hyd& modele_turbulence = eqn_NS_.valeur().modele_turbulence();

      const Modele_turbulence_hyd_K_Eps& modele_turbulence_keps = ref_cast(Modele_turbulence_hyd_K_Eps,modele_turbulence.valeur());

      modele_K_Eps_ = modele_turbulence_keps;

      eqn_transport_K_Eps_ = ref_cast(Transport_K_Eps,modele_K_Eps_.valeur().eqn_transp_K_Eps());
    }
}

void Tenseur_Reynolds_Externe_VDF_Face::associer_domaines(const Domaine_dis& domaine_dis,
                                                          const Domaine_Cl_dis& domaine_Cl_dis)
{
  le_dom_VDF = ref_cast(Domaine_VDF, domaine_dis.valeur());
  le_dom_Cl_VDF = ref_cast(Domaine_Cl_VDF, domaine_Cl_dis.valeur());

  nelem_ = le_dom_VDF.valeur().nb_elem();
}

void Tenseur_Reynolds_Externe_VDF_Face::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const Domaine_Cl_VDF& domaine_Cl_VDF = le_dom_Cl_VDF.valeur();
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const IntVect& orientation = domaine_VDF.orientation();
  const DoubleVect& porosite_surf = domaine_Cl_VDF.equation().milieu().porosite_face();
  const DoubleVect& volumes_entrelaces = domaine_VDF.volumes_entrelaces();

  int ndeb,nfin,ncomp,num_face,elem1,elem2;
  double vol;

  if (sub_type(Champ_Uniforme,la_source.valeur()))
    {
      const DoubleVect& s = la_source->valeurs();

      // Boucle sur les conditions limites pour traiter les faces de bord

      for (int n_bord=0; n_bord<domaine_VDF.nb_front_Cl(); n_bord++)
        {

          // pour chaque Condition Limite on regarde son type
          // Si face de Dirichlet ou de Symetrie on ne fait rien
          // Si face de Neumann on calcule la contribution au terme source

          const Cond_lim& la_cl = domaine_Cl_VDF.les_conditions_limites(n_bord);

          if (sub_type(Periodique,la_cl.valeur()))
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              ndeb = le_bord.num_premiere_face();
              nfin = ndeb + le_bord.nb_faces();

              for (num_face=ndeb; num_face<nfin; num_face++)
                {
                  vol = volumes_entrelaces(num_face);
                  //  for (ncomp=0; ncomp<dimension ; ncomp++)
                  ncomp = orientation(num_face);
                  secmem(num_face)+= s(ncomp)*vol;
                }
            }
          else if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
            {
              //              const Neumann_sortie_libre& la_cl_neumann = ref_cast(Neumann_sortie_libre,la_cl.valeur());
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              ndeb = le_bord.num_premiere_face();
              nfin = ndeb + le_bord.nb_faces();

              for (num_face=ndeb; num_face<nfin; num_face++)
                {
                  vol = volumes_entrelaces(num_face)*porosite_surf(num_face);
                  ncomp = orientation(num_face);
                  secmem(num_face)+= s(ncomp)*vol;
                }

            }
          else if (sub_type(Symetrie,la_cl.valeur()))
            ;
          else if ( (sub_type(Dirichlet,la_cl.valeur()))
                    ||
                    (sub_type(Dirichlet_homogene,la_cl.valeur()))
                  )
            {
              // do nothing
              ;
            }
        }


      // Boucle sur les faces internes

      ndeb = domaine_VDF.premiere_face_int();
      for (num_face =domaine_VDF.premiere_face_int(); num_face<domaine_VDF.nb_faces(); num_face++)
        {
          vol = volumes_entrelaces(num_face)*porosite_surf(num_face);
          ncomp = orientation(num_face);
          secmem(num_face) += s(ncomp)*vol;

        }
    }
  else // le champ source n'est plus uniforme
    {
      const DoubleTab& s = la_source->valeurs();

      // Boucle sur les conditions limites pour traiter les faces de bord

      for (int n_bord=0; n_bord<domaine_VDF.nb_front_Cl(); n_bord++)
        {

          // pour chaque Condition Limite on regarde son type
          // Si face de Dirichlet ou de Symetrie on ne fait rien
          // Si face de Neumann on calcule la contribution au terme source

          const Cond_lim& la_cl = domaine_Cl_VDF.les_conditions_limites(n_bord);

          if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
            {

              //            const Neumann_sortie_libre& la_cl_neumann = ref_cast(Neumann_sortie_libre,la_cl.valeur());
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              ndeb = le_bord.num_premiere_face();
              nfin = ndeb + le_bord.nb_faces();

              for (num_face=ndeb; num_face<nfin; num_face++)
                {
                  vol = volumes_entrelaces(num_face)*porosite_surf(num_face);
                  ncomp = orientation(num_face);
                  elem1 = face_voisins(num_face,0);

                  if (elem1 != -1)
                    secmem(num_face)+= s(elem1,ncomp)*vol;
                  else
                    {
                      elem2 = face_voisins(num_face,1);
                      secmem(num_face)+= s(elem2,ncomp)*vol;
                    }
                }

            }
          else if (sub_type(Symetrie,la_cl.valeur()))
            ;
          else if ( (sub_type(Dirichlet,la_cl.valeur()))
                    ||
                    (sub_type(Dirichlet_homogene,la_cl.valeur()))
                  )
            ;
          else if (sub_type(Periodique,la_cl.valeur()))
            {
              double s_face;
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              ndeb = le_bord.num_premiere_face();
              nfin = ndeb + le_bord.nb_faces();
              for (num_face=ndeb; num_face<nfin; num_face++)
                {
                  vol = volumes_entrelaces(num_face)*porosite_surf(num_face);
                  ncomp = orientation(num_face);
                  s_face = 0.5*( s(face_voisins(num_face,0),ncomp) + s(face_voisins(num_face,1),ncomp) );
                  secmem(num_face) += s_face*vol;
                }

            }
        }
      // Boucle sur les faces internes

      double s_face;
      ndeb = domaine_VDF.premiere_face_int();
      for (num_face =domaine_VDF.premiere_face_int(); num_face<domaine_VDF.nb_faces(); num_face++)
        {

          vol = volumes_entrelaces(num_face)*porosite_surf(num_face);
          ncomp = orientation(num_face);
          s_face = 0.5*( s(face_voisins(num_face,0),ncomp) + s(face_voisins(num_face,1),ncomp) );
          secmem(num_face) += s_face*vol;

        }
    }
}

DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}

const DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::get_bij() const
{
  return bij_;
}

const DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::get_bij_NL() const
{
  return bij_NL_;
}

const DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::get_g1() const
{
  return g1_;
}

void Tenseur_Reynolds_Externe_VDF_Face::Calcul_RSLambdaT()
{
//  Cout << "Calling Calcul_RSLambdaT " << endl;
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const Domaine_Cl_VDF& domaine_Cl_VDF = le_dom_Cl_VDF.valeur();

  int nb_elem_tot=domaine_VDF.nb_elem_tot();
  const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF,eqn_NS_->inconnue().valeur() );
  assert (vitesse.valeurs().line_size() == 1);
  DoubleTab gij(nb_elem_tot,dimension,dimension, vitesse.valeurs().line_size());
  ref_cast_non_const(Champ_Face_VDF,vitesse).calcul_duidxj( vitesse.valeurs(),gij,domaine_Cl_VDF );


  DoubleTab& lambda_1 = modele_K_Eps_.valeur().get_l1();
  DoubleTab& lambda_2 = modele_K_Eps_.valeur().get_l2();
  DoubleTab& lambda_3 = modele_K_Eps_.valeur().get_l3();
  DoubleTab& lambda_4 = modele_K_Eps_.valeur().get_l4();
  DoubleTab& lambda_5 = modele_K_Eps_.valeur().get_l5();

  const DoubleTab& K_eps = eqn_transport_K_Eps_->inconnue().valeurs();

  DoubleTab S_etoile(nb_elem_tot,dimension,dimension);
  DoubleTab R_etoile(nb_elem_tot,dimension,dimension);

  DoubleTab S2(nb_elem_tot,dimension,dimension);
  DoubleTab R2(nb_elem_tot,dimension,dimension);
  DoubleTab RS(nb_elem_tot,dimension,dimension);
  DoubleTab SR(nb_elem_tot,dimension,dimension);

  DoubleTab S3(nb_elem_tot,dimension,dimension);
  DoubleTab R2S(nb_elem_tot,dimension,dimension);
  DoubleTab RS2(nb_elem_tot,dimension,dimension);
  DoubleTab S2R(nb_elem_tot,dimension,dimension);
  DoubleTab SR2(nb_elem_tot,dimension,dimension);

  DoubleTab R2S2(nb_elem_tot,dimension,dimension);
  DoubleTab S2R2(nb_elem_tot,dimension,dimension);
  DoubleTab SRS2(nb_elem_tot,dimension,dimension);
  DoubleTab R2SR(nb_elem_tot,dimension,dimension);
  DoubleTab RSR2(nb_elem_tot,dimension,dimension);
  DoubleTab S2RS(nb_elem_tot,dimension,dimension);

  DoubleTab RS2R2(nb_elem_tot,dimension,dimension);
  DoubleTab R2S2R(nb_elem_tot,dimension,dimension);

  DoubleTab L1Id(nb_elem_tot,dimension,dimension);
  DoubleTab L2Id(nb_elem_tot,dimension,dimension);
  DoubleTab L4Id(nb_elem_tot,dimension,dimension);
  DoubleTab L5Id(nb_elem_tot,dimension,dimension);

  DoubleTab tab_zero(nb_elem_tot,dimension,dimension);

  for (int elem=0; elem<domaine_VDF.nb_elem(); elem++)
    {

      double k_sur_eps = K_eps(elem,0) / ( K_eps(elem,1) + 1.e-15 );

      lambda_1(elem) = 0.;
      lambda_2(elem) = 0.;
      lambda_3(elem) = 0.;
      lambda_4(elem) = 0.;
      lambda_5(elem) = 0.;

      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            tab_zero(elem,i,j) = 0.;
            S_etoile(elem,i,j) = 0.5*( gij(elem,i,j,0) + gij(elem,j,i,0) ) * k_sur_eps;
            R_etoile(elem,i,j) = 0.5*( gij(elem,i,j,0) - gij(elem,j,i,0) ) * k_sur_eps;
          }

      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            S2(elem,i,j)   = 0.;
            R2(elem,i,j)   = 0.;
            SR(elem,i,j)   = 0.;
            RS(elem,i,j)   = 0.;
            L1Id(elem,i,j) = 0.;
            L2Id(elem,i,j) = 0.;

            for (int k=0; k<Objet_U::dimension; k++)
              {
                S2(elem,i,j) +=  S_etoile(elem,i,k)*S_etoile(elem,k,j) ;
                R2(elem,i,j) +=  R_etoile(elem,i,k)*R_etoile(elem,k,j) ;
                SR(elem,i,j) +=  S_etoile(elem,i,k)*R_etoile(elem,k,j) ;
                RS(elem,i,j) +=  R_etoile(elem,i,k)*S_etoile(elem,k,j) ;
              }

            if (i==j)
              {
                lambda_1(elem) += S2(elem,i,j);
                lambda_2(elem) += R2(elem,i,j);
              }
          }

      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            if (i==j)
              {
                L1Id(elem,i,j) = lambda_1(elem)/3.;
                L2Id(elem,i,j) = lambda_2(elem)/3.;
              }
          }

      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            S3(elem,i,j)   = 0.;
            R2S(elem,i,j)  = 0.;
            RS2(elem,i,j)  = 0.;
            S2R(elem,i,j)  = 0.;
            SR2(elem,i,j)  = 0.;
            L4Id(elem,i,j) = 0.;

            for (int k=0; k<Objet_U::dimension; k++)
              {
                S3(elem,i,j)  +=  S_etoile(elem,i,k)*S2(elem,k,j) ;
                R2S(elem,i,j) +=  R2(elem,i,k)*S_etoile(elem,k,j) ;
                RS2(elem,i,j) +=  R_etoile(elem,i,k)*S2(elem,k,j) ;
                S2R(elem,i,j) +=  S2(elem,i,k)*R_etoile(elem,k,j) ;
                SR2(elem,i,j) +=  S_etoile(elem,i,k)*R2(elem,k,j) ;
              }

            if (i==j)
              {
                lambda_3(elem) += S3(elem,i,j);
                lambda_4(elem) += R2S(elem,i,j);
              }
          }

      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            if (i==j)
              {
                L4Id(elem,i,j) = lambda_4(elem)*2./3.;
              }
          }

      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            S2R2(elem,i,j) = 0.;
            R2S2(elem,i,j) = 0.;
            SRS2(elem,i,j) = 0.;
            R2SR(elem,i,j) = 0.;
            RSR2(elem,i,j) = 0.;
            S2RS(elem,i,j) = 0.;
            L5Id(elem,i,j) = 0.;

            for (int k=0; k<Objet_U::dimension; k++)
              {
                S2R2(elem,i,j) +=  S2(elem,i,k)*R2(elem,k,j) ;
                R2S2(elem,i,j) +=  R2(elem,i,k)*S2(elem,k,j) ;
                SRS2(elem,i,j) +=  SR(elem,i,k)*S2(elem,k,j) ;
                R2SR(elem,i,j) +=  R2(elem,i,k)*SR(elem,k,j) ;
                RSR2(elem,i,j) +=  RS(elem,i,k)*R2(elem,k,j) ;
                S2RS(elem,i,j) +=  S2(elem,i,k)*RS(elem,k,j) ;
              }

            if (i==j)
              {
                lambda_5(elem) += R2S2(elem,i,j);
              }
          }

      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            if (i==j)
              {
                L5Id(elem,i,j) = lambda_5(elem)*2./3.;
              }
          }

      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            RS2R2(elem,i,j) = 0.;
            R2S2R(elem,i,j) = 0.;
            for (int k=0; k<Objet_U::dimension; k++)
              {
                RS2R2(elem,i,j) +=  RS2(elem,i,k)*R2(elem,k,j) ;
                R2S2R(elem,i,j) +=  R2(elem,i,k)*S2R(elem,k,j) ;
              }
          }
    }

  lambda_1_etoile_ = lambda_1;
  lambda_2_etoile_ = lambda_2;
  lambda_3_etoile_ = lambda_3;
  lambda_4_etoile_ = lambda_4;
  lambda_5_etoile_ = lambda_5;
  T0_etoile_ = tab_zero;
  for (int elem=0; elem<domaine_VDF.nb_elem(); elem++)
    {
      T0_etoile_(elem,0,0) = 1./3.;
      T0_etoile_(elem,0,1) = 0.;
      T0_etoile_(elem,0,2) = 0.;
      T0_etoile_(elem,1,0) = 0.;
      T0_etoile_(elem,1,1) = -1./6.;
      T0_etoile_(elem,1,2) = 0.;
      T0_etoile_(elem,2,0) = 0.;
      T0_etoile_(elem,2,1) = 0.;
      T0_etoile_(elem,2,2) = -1./6.;
    }


  T1_etoile_  = S_etoile;

  if (apply_Nisizima_)
    {
      T2_etoile_  = SR;
      T2_etoile_ -= RS;
      T3_etoile_  = S2;
      T3_etoile_ -= L1Id;
      T4_etoile_  = R2;
      T4_etoile_ -= L2Id;
      T5_etoile_  = tab_zero;
      T6_etoile_  = tab_zero;
      T7_etoile_ = tab_zero;
      T8_etoile_  = tab_zero;
      T9_etoile_ = tab_zero;
      T10_etoile_ = tab_zero;
    }
  else
    {
      if (!T_list_.empty())
        {
          if ( std::find(T_list_.begin(), T_list_.end(), 2) != T_list_.end() )
            {
              //          cout << "using T2" << endl;
              T2_etoile_  = SR;
              T2_etoile_ -= RS;
            }
          else
            T2_etoile_  = tab_zero;
          if ( std::find(T_list_.begin(),T_list_.end(), 3) != T_list_.end() )
            {
              //          cout << "using T3" << endl;
              T3_etoile_  = S2;
              T3_etoile_ -= L1Id;
            }
          else
            T3_etoile_  = tab_zero;
          if ( std::find(T_list_.begin(), T_list_.end(), 4) != T_list_.end() )
            {
              //          cout << "using T4" << endl;
              T4_etoile_  = R2;
              T4_etoile_ -= L2Id;
            }
          else
            T4_etoile_  = tab_zero;
          if ( std::find(T_list_.begin(), T_list_.end(), 5) != T_list_.end() )
            {
              //          cout << "using T5" << endl;
              T5_etoile_  = RS2;
              T5_etoile_ -= S2R;
            }
          else
            T5_etoile_  = tab_zero;
          if ( std::find(T_list_.begin(), T_list_.end(), 6) != T_list_.end() )
            {
              //          cout << "using T6" << endl;
              T6_etoile_  = R2S;
              T6_etoile_ += SR2;
              T6_etoile_ -= L4Id;
            }
          else
            T6_etoile_  = tab_zero;
          if ( std::find(T_list_.begin(), T_list_.end(), 7) != T_list_.end() )
            {
              //          cout << "using T7" << endl;
              T7_etoile_  = RSR2;
              T7_etoile_ -= R2SR;
            }
          else
            T7_etoile_  = tab_zero;
          if ( std::find(T_list_.begin(), T_list_.end(), 8) != T_list_.end() )
            {
              //          cout << "using T8" << endl;
              T8_etoile_  = SRS2;
              T8_etoile_ -= S2RS;
            }
          else
            T8_etoile_  = tab_zero;
          if ( std::find(T_list_.begin(), T_list_.end(), 9) != T_list_.end() )
            {
              //          cout << "using T9" << endl;
              T9_etoile_  = R2S2;
              T9_etoile_ += S2R2;
              T9_etoile_ -= L5Id;
            }
          else
            T9_etoile_  = tab_zero;
          if ( std::find(T_list_.begin(), T_list_.end(), 10) != T_list_.end() )
            {
              //          cout << "using T10" << endl;
              T10_etoile_  = RS2R2;
              T10_etoile_ -= R2S2R;
            }
          else
            T10_etoile_  = tab_zero;
        }
      else
        {
          cerr << "T_list is empty."  << endl;
          T2_etoile_  = tab_zero;
          T3_etoile_  = tab_zero;
          T4_etoile_  = tab_zero;
          T5_etoile_  = tab_zero;
          T6_etoile_  = tab_zero;
          T7_etoile_ = tab_zero;
          T8_etoile_  = tab_zero;
          T9_etoile_ = tab_zero;
          T10_etoile_ = tab_zero;
        }
    }

}

void Tenseur_Reynolds_Externe_VDF_Face::mettre_a_jour(double temps)
{
  int nb_elem_tot = le_dom_VDF.valeur().nb_elem_tot();
  int nb_faces_tot = le_dom_VDF.valeur().nb_faces_tot();
  const IntTab& face_vois = le_dom_VDF.valeur().face_voisins();
  const DoubleVect& volumes = le_dom_VDF.valeur().volumes();

  DoubleTab valeurs_source(nb_elem_tot,dimension);
  DoubleTab tenseur_reynolds_elements(nb_elem_tot,dimension,dimension);
  DoubleTab tenseur_reynolds_faces(nb_faces_tot,dimension,dimension);

  valeurs_source = 0;
  tenseur_reynolds_elements = 0;
  tenseur_reynolds_faces = 0.;

  if (apply_Nisizima_)
    {
      Calcul_RSLambdaT();
      Calcul_Tenseur_Reynolds_NL( tenseur_reynolds_elements );

      for (int fac=0; fac<nb_faces_tot; fac++)
        {
          int elem1 = face_vois(fac,0);
          int elem2 = face_vois(fac,1);
          double vol = 0;
          if (elem1!=-1)
            {
              for (int i=0; i<Objet_U::dimension; i++)
                for (int j=0; j<Objet_U::dimension; j++)
                  {
                    tenseur_reynolds_faces(fac,i,j) +=  tenseur_reynolds_elements(elem1,i,j)*volumes(elem1);
                  }
              vol += volumes(elem1);
            }
          if (elem2!=-1)
            {
              for (int i=0; i<Objet_U::dimension; i++)
                for (int j=0; j<Objet_U::dimension; j++)
                  {
                    tenseur_reynolds_faces(fac,i,j) +=  tenseur_reynolds_elements(elem2,i,j)*volumes(elem2);
                  }
              vol += volumes(elem2);
            }
          for (int i=0; i<Objet_U::dimension; i++)
            for (int j=0; j<Objet_U::dimension; j++)
              {
                tenseur_reynolds_faces(fac,i,j) /= vol;
              }
        }

      const IntTab& elem_faces   = le_dom_VDF.valeur().elem_faces();
      const IntTab& face_voisins = le_dom_VDF.valeur().face_voisins();
      const DoubleVect& inverse_vol  = le_dom_VDF.valeur().inverse_volumes();
      const int nb_faces_elem=elem_faces.line_size();


      int facei;
      double signe=0.;
      double div=0.;

      for (int elem=0; elem<nb_elem_tot; elem++)
        {
          //Calcul de la divergence par element
          for (int i=0; i<Objet_U::dimension; i++)
            {
              div=0.;
              for (int facei_loc=0; facei_loc<nb_faces_elem; facei_loc++)
                {
                  facei=elem_faces(elem,facei_loc);
                  signe=(face_voisins(facei,0)==elem)? 1.:-1.;

                  for (int j=0; j<Objet_U::dimension; j++)
                    {
                      div+=signe*le_dom_VDF.valeur().face_normales(facei,j)*tenseur_reynolds_faces(facei,i,j);
                    }
                }
              div*= inverse_vol(elem);
              valeurs_source(elem,i) = div;
            }
        }
      valeurs_source.echange_espace_virtuel();
      la_source->valeurs( ) = valeurs_source;

    }
  else
    {
      if ((idt_2 ==0) || (idt_2 % idt_mettre_a_jour == 0))
        {

          if (tbnn->is_canal_carre_)
            {
              Calcul_RSLambdaT();
            }
          Calcul_Tenseur_Reynolds_NL( tenseur_reynolds_elements );
          for (int fac=0; fac<nb_faces_tot; fac++)
            {
              int elem1 = face_vois(fac,0);
              int elem2 = face_vois(fac,1);
              double vol = 0;
              if (elem1!=-1)
                {
                  for (int i=0; i<Objet_U::dimension; i++)
                    for (int j=0; j<Objet_U::dimension; j++)
                      {
                        tenseur_reynolds_faces(fac,i,j) +=  tenseur_reynolds_elements(elem1,i,j)*volumes(elem1);
                      }
                  vol += volumes(elem1);
                }
              if (elem2!=-1)
                {
                  for (int i=0; i<Objet_U::dimension; i++)
                    for (int j=0; j<Objet_U::dimension; j++)
                      {
                        tenseur_reynolds_faces(fac,i,j) +=  tenseur_reynolds_elements(elem2,i,j)*volumes(elem2);
                      }
                  vol += volumes(elem2);
                }
              for (int i=0; i<Objet_U::dimension; i++)
                for (int j=0; j<Objet_U::dimension; j++)
                  {
                    tenseur_reynolds_faces(fac,i,j) /= vol;
                  }
            }

          const IntTab& elem_faces   = le_dom_VDF.valeur().elem_faces();
          const IntTab& face_voisins = le_dom_VDF.valeur().face_voisins();
          const DoubleVect& inverse_vol  = le_dom_VDF.valeur().inverse_volumes();
          const int nb_faces_elem=elem_faces.line_size();


          int facei;
          double signe=0.;
          double div=0.;
          double factor=0.;

          if (idt_2 <= 0)
            factor = 1;
          else if (T_stab > 0)
            {
              factor = 1 - exp(-t_stab/T_stab);
            }
          else
            {
              factor = 1;
            }


          for (int elem=0; elem<nb_elem_tot; elem++)
            {
              //Calcul de la divergence par element
              for (int i=0; i<Objet_U::dimension; i++)
                {
                  div=0.;
                  for (int facei_loc=0; facei_loc<nb_faces_elem; facei_loc++)
                    {
                      facei=elem_faces(elem,facei_loc);
                      signe=(face_voisins(facei,0)==elem)? 1.:-1.;

                      for (int j=0; j<Objet_U::dimension; j++)
                        {
                          div+=signe*le_dom_VDF.valeur().face_normales(facei,j)*tenseur_reynolds_faces(facei,i,j);
                        }
                    }
                  div*= inverse_vol(elem);
                  valeurs_source(elem,i) = div;
                  valeurs_source(elem,i) *= factor;
                }
            }
          valeurs_source.echange_espace_virtuel();
          la_source->valeurs( ) = valeurs_source;
          if ((idt_2 ==0) & (!is_only_g1_))
            cout << "Source term updated by NN at time step " << idt_2  << " S^{n} = " << " S^{n} * " << factor << endl;

          if ((idt_2>0) & (!is_only_g1_))
            {
              for (int elem=0; elem<nb_elem_tot; elem++)
                {
                  for (int i=0; i<Objet_U::dimension; i++)
                    {
                      la_source->valeurs()(elem,i) = relax_RST * valeurs_source(elem,i) + (1-relax_RST) * s_avant_(elem,i) ;
                    }
                }
              cout << "Source term updated by NN at time step " << idt_2  << " S^{n} = " << relax_RST << " * S^{n} * " << factor << " + (1-" << relax_RST << ") * S^{n-1} " << endl;
            }

          if (is_only_g1_)
            {
              for (int elem=0; elem<nb_elem_tot; elem++)
                {
                  for (int i=0; i<Objet_U::dimension; i++)
                    {
                      la_source->valeurs()(elem,i) = 0. ;
                    }
                }
              cout << "Source term updated by NN at time step " << idt_2  << " S^{n} = 0"  << endl;
            }
        }
    }

  s_avant_ = la_source->valeurs();
  idt_2 = idt_2 + 1;
  t_stab = t_stab + 1;
}

void Tenseur_Reynolds_Externe_VDF_Face::completer()
{
  Source_base::completer();

  la_source.DERIV(Champ_Don_base)::typer( "Champ_Fonc_xyz" );
}


DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::Calcul_Tenseur_Reynolds_NL(DoubleTab& resu)
{
  DoubleTab bij_NL;
  if (apply_Nisizima_)
    {
      bij_NL = Calcul_bij_NL_Nisizima(resu);
    }
  else
    {

      if (tbnn->is_canal_plan_)
        {
          bij_NL = Calcul_bij_NL_TBNN(resu);
        }
      else if (tbnn->is_canal_carre_)
        {
          bij_NL = Calcul_bij_NL_TBNN_carre(resu);
        }

    }

  const DoubleTab& K_eps = eqn_transport_K_Eps_->inconnue().valeurs();

  for (int elem=0; elem<nelem_; elem++)
    {
      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            resu(elem,i,j) = bij_NL(elem,i,j) ;

            if (i==j)
              {
                resu(elem,i,j) += 1./3.;
              }

            if ((apply_keps_NN_) || (apply_keps_DNS_))
              {
                resu(elem,i,j) *= -2. * KEps_true_(elem,0); // -Rij = -ui'uj' = -2k*(bij + 1/3 delta_ij)
              }

            if ((!apply_keps_NN_) && (!apply_keps_DNS_))
              {
                resu(elem,i,j) *= -2. * K_eps(elem,0);
              }
            if (T_list_.size() == 1 && T_list_[0] == 1)
              {
                resu(elem,i,j) = 0.;
              }

          }
    }


  return resu;
}


DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::Calcul_bij_NL_TBNN(DoubleTab& resu)
{

  DoubleTab& resu_NL = Calcul_bij_TBNN(resu);

  double alpha;
  const DoubleTab& K_eps = eqn_transport_K_Eps_->inconnue().valeurs();

  const Domaine_Cl_VDF& domaine_Cl_VDF = le_dom_Cl_VDF.valeur();
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();

  int nb_elem_tot=domaine_VDF.nb_elem_tot();
  const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF,eqn_NS_->inconnue().valeur() );
  assert (vitesse.valeurs().line_size() == 1);
  DoubleTab gij(nb_elem_tot,dimension,dimension, vitesse.valeurs().line_size());
  ref_cast_non_const(Champ_Face_VDF,vitesse).calcul_duidxj( vitesse.valeurs(),gij,domaine_Cl_VDF );

  for (int elem=0; elem<nelem_; elem++)
    {
      alpha = tbnn->compute_alpha(K_eps(elem, 0), K_eps(elem,1), gij(elem,0,1,0) );
      resu_NL(elem,0,1) -= 0.5 * g1_(elem) * alpha;
      resu_NL(elem,1,0) -= 0.5 * g1_(elem) * alpha;
    }

  return resu_NL;
}

DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::Calcul_bij_NL_TBNN_carre(DoubleTab& resu)
{

  Calcul_bij_TBNN_carre();
//  DoubleTab resu_NL(nelem_,dimension,dimension);
  resu = 0;

  vector<DoubleTab*> g_map =
  {
    &g0_, &g1_, &g2_, &g3_, &g4_, &g5_, &g6_, &g7_, &g8_, &g9_, &g10_
  };

  vector<DoubleTab*> T_etoile_map =
  {
    &T0_etoile_, &T1_etoile_, &T2_etoile_, &T3_etoile_, &T4_etoile_,
    &T5_etoile_, &T6_etoile_, &T7_etoile_, &T8_etoile_, &T9_etoile_, &T10_etoile_,
  };

  for (int elem=0; elem<nelem_; elem++)
    {
      for (int i=0; i<Objet_U::dimension; i++)
        {
          for (int j=0; j<Objet_U::dimension; j++)
            {
              for (size_t k = 0; k < T_list_.size(); k++)
                {
                  int t_index = T_list_[k];
                  if (t_index != 1)
                    {
                      resu(elem,i,j) += (*g_map[t_index])(elem) *  (*T_etoile_map[t_index])(elem,i,j);
                    }
                }
            }
        }
    }
  DoubleTab& bij_NL_tab = modele_K_Eps_.valeur().get_bij_NL();
  for (int elem=0; elem<nelem_; elem++)
    {
      bij_NL_tab(elem,0) = resu(elem,0,0);
      bij_NL_tab(elem,1) = resu(elem,0,1);
      bij_NL_tab(elem,2) = resu(elem,0,2);
      bij_NL_tab(elem,3) = resu(elem,1,1);
      bij_NL_tab(elem,4) = resu(elem,1,2);
      bij_NL_tab(elem,5) = resu(elem,2,2);
    }
  bij_NL_ = resu;


  return resu;
}

DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::Calcul_bij_TBNN(DoubleTab& resu)
{
  vector<double> b;
  double alpha, y_plus, Re_t;
  double y_maille_paroi =-1, y_elem;
  //  vector<vector<double>> T;
  DoubleTab g1(nelem_);
  const DoubleTab& K_eps = eqn_transport_K_Eps_->inconnue().valeurs();
  double u_t = 0.;

  double dudy_paroi ;

  vector<double> y_plus_wall, xpx, xpz, Re_true;
  double x_elem, z_elem;
  int position_base = -1;

  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const Domaine_Cl_VDF& domaine_Cl_VDF = le_dom_Cl_VDF.valeur();
  int nb_elem_tot=domaine_VDF.nb_elem_tot();
  const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF,eqn_NS_->inconnue().valeur() );
  assert (vitesse.valeurs().line_size() == 1);

  DoubleTab gij(nb_elem_tot,dimension,dimension, vitesse.valeurs().line_size());
  ref_cast_non_const(Champ_Face_VDF,vitesse).calcul_duidxj( vitesse.valeurs(),gij,domaine_Cl_VDF );

  double nu = eqn_NS_->fluide().viscosite_cinematique().valeurs()[0];

  y_plus_wall.resize(9);
  Re_true.resize(9);
  xpx.resize(9);
  xpz.resize(9);
  double g1_max=0;

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
              elem_paroi = le_dom_VDF->face_voisins(num_face,1);//1 means that I want the elem above the face (which is the bottom paroi so ok)

              xpx[num_face-ndeb] =le_dom_VDF->xp(elem_paroi, 0) ; // x=0
              xpz[num_face-ndeb] =le_dom_VDF->xp(elem_paroi, 2) ; // z=2

              y_maille_paroi =le_dom_VDF->xp(elem_paroi, 1) ;

              dudy_paroi = gij(elem_paroi,0,1,0);
              u_t = sqrt( abs(nu* dudy_paroi));
              y_plus_wall[num_face-ndeb] = y_maille_paroi * u_t / nu;
              Re_true [num_face-ndeb] =  u_t / nu;
            }
        }
    }



  for (int elem=0; elem<nelem_; elem++)
    {

      x_elem =le_dom_VDF->xp(elem, 0) ; // x=0
      y_elem =le_dom_VDF->xp(elem, 1) ; // y=1
      z_elem =le_dom_VDF->xp(elem, 2) ; // z=2

      for (unsigned int i = 0; i < y_plus_wall.size(); i++)
        if ( fabs (x_elem-xpx[i]) <1e-3 and fabs (z_elem-xpz[i]) <1e-3)
          {
            position_base = i;
            break;
          }

      alpha = tbnn->compute_alpha(K_eps(elem, 0), K_eps(elem,1), gij(elem,0,1,0) ) + 1e-40;
      y_plus = tbnn->compute_y_plus(y_plus_wall[position_base],y_maille_paroi, y_elem);
      Re_t = tbnn->compute_Re_t(y_plus_wall[position_base],y_maille_paroi);


      // prediction by neural network
      //Cerr << " Testing PCF TRE_VDF_Face.cpp line 1213 " << finl;
      b = tbnn->predict(alpha, y_plus, Re_t);

      int idx_diag[] = {0,3,5};
      int idx_nondiag[] = {1,2,4};

      for (int i=0; i<3; i++ )
        {
          if (b[idx_diag[i]] < -1./3.)
            {
              b[idx_diag[i]] = -1./3.;
            }
          else if (b[idx_diag[i]] > 2./3.)
            {
              b[idx_diag[i]] = 2./3.;
            }
        }

      for (int i=0; i<3; i++ )
        {
          if (b[idx_nondiag[i]] < -1./2.)
            {
              b[idx_nondiag[i]] = -1./2.;
            }
          else if (b[idx_nondiag[i]] > 1./2.)
            {
              b[idx_nondiag[i]] = 1./2.;
            }
        }
      // copy of the predictions in resu array
      resu(elem,0,0) = b[0];
      resu(elem,0,1) = b[1];
      resu(elem,0,2) = b[2];
      resu(elem,1,0) = b[1];
      resu(elem,1,1) = b[3];
      resu(elem,1,2) = b[4];
      resu(elem,2,0) = b[2];
      resu(elem,2,1) = b[4];
      resu(elem,2,2) = b[5];

      g1(elem) = tbnn->get_g1(resu(elem,0,1), alpha);

      if (g1(elem) > g1_max)
        {
          g1(elem) = g1_max;
        }
    }

  // save predictions and g1 values
  g1_ = g1;

  DoubleTab& bij_tab = modele_K_Eps_.valeur().get_bij();
  for (int elem=0; elem<nelem_; elem++)
    {
      bij_tab(elem,0) = resu(elem,0,0);
      bij_tab(elem,1) = resu(elem,0,1);
      bij_tab(elem,2) = resu(elem,0,2);
      bij_tab(elem,3) = resu(elem,1,1);
      bij_tab(elem,4) = resu(elem,1,2);
      bij_tab(elem,5) = resu(elem,2,2);
    }
  bij_ = resu;

  return resu;
}

void Tenseur_Reynolds_Externe_VDF_Face::Calcul_bij_TBNN_carre()
{
  vector<double> k_ec, eps;
  DoubleTab& K_eps = eqn_transport_K_Eps_->inconnue().valeurs();
  vector<double> b, lambda;
  double y_plus, z_plus, Re_t;
//  double y_maille_paroi =-1;
//  double z_maille_paroi=-1;
  vector<vector<double>> T;
  vector<vector<double>> T_used;
  DoubleTab g0(nelem_);
  DoubleTab g1(nelem_);
  DoubleTab g2(nelem_);
  DoubleTab g3(nelem_);
  DoubleTab g4(nelem_);
  DoubleTab g5(nelem_);
  DoubleTab g6(nelem_);
  DoubleTab g7(nelem_);
  DoubleTab g8(nelem_);
  DoubleTab g9(nelem_);
  DoubleTab g10(nelem_);

  DoubleTab resu(nelem_,dimension,dimension);
  DoubleTab resu_keps(nelem_, 2);
  double u_t = 0., sum_surface=0.;
  double dudy_paroi, dudz_paroi ;
  double y_elem, z_elem;
//  int position_base = -1;
//  int position_gauche = -1;
  int ori;

  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const Domaine_Cl_VDF& domaine_Cl_VDF = le_dom_Cl_VDF.valeur();
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const IntVect& orientation = domaine_VDF.orientation();
  const DoubleTab& xp = domaine_VDF.xp();
  const DoubleVect& face_surfaces = domaine_VDF.face_surfaces();


  int nb_elem_tot=domaine_VDF.nb_elem_tot();
  const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF,eqn_NS_->inconnue().valeur() );
  assert (vitesse.valeurs().line_size() == 1);

  DoubleTab gij(nb_elem_tot,dimension,dimension, vitesse.valeurs().line_size());
  ref_cast_non_const(Champ_Face_VDF,vitesse).calcul_duidxj( vitesse.valeurs(),gij,domaine_Cl_VDF );

  double nu = eqn_NS_->fluide().viscosite_cinematique().valeurs()[0];

  // Boucle sur les bords pour calculer u_tau

//  auto start = high_resolution_clock::now();
//  cout << "computing u_tau..." << endl;

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
//  cout << "u_t _RST = " << u_t << endl;
//  end = high_resolution_clock::now();
//  duration = duration_cast<microseconds>(end - start).count();
//  cout << "y_plus_wall computed..." << endl;
//  cout << "Time for y_plus_wall " << duration << " microseconds" << endl;

//  auto start3 = high_resolution_clock::now();
//  cout << "applying NN..." << endl;
  //  init of lambda and T arrays
  lambda.resize(5);
  T.resize(11);
  for (int i=0; i<11; i++)
    T[i].resize(6);
  T_used.resize(T_list_.size());

  if ( std::find(T_list_.begin(), T_list_.end(), 0) != T_list_.end() )
    {
      T[0][0] = 1./3.;
      T[0][1] = 0.;
      T[0][2] = 0.;
      T[0][3] = -1./6.;
      T[0][4] = 0.;
      T[0][5] = -1./6.;
    }
  else
    {
      T[0][0] = 0.;
      T[0][1] = 0.;
      T[0][2] = 0.;
      T[0][3] = 0.;
      T[0][4] = 0.;
      T[0][5] = 0.;
    }

  for (int elem=0; elem<nelem_; elem++)
    {

      lambda[0] = lambda_1_etoile_(elem);
      lambda[1] = lambda_2_etoile_(elem);
      lambda[2] = lambda_3_etoile_(elem);
      lambda[3] = lambda_4_etoile_(elem);
      lambda[4] = lambda_5_etoile_(elem);
      unsigned int k = 0;
      for (int i=0; i<Objet_U::dimension; i++)
        {
          for (int j=i; j<Objet_U::dimension; j++)
            {
              T[1][k] = T1_etoile_(elem,i,j);
              T[2][k] = T2_etoile_(elem,i,j);
              T[3][k] = T3_etoile_(elem,i,j);
              T[4][k] = T4_etoile_(elem,i,j);
              T[5][k] = T5_etoile_(elem,i,j);
              T[6][k] = T6_etoile_(elem,i,j);
              T[7][k] = T7_etoile_(elem,i,j);
              T[8][k] = T8_etoile_(elem,i,j);
              T[9][k] = T9_etoile_(elem,i,j);
              T[10][k] = T10_etoile_(elem,i,j);
              k++;
            }
        }


      for (size_t idx = 0; idx < T_list_.size(); ++idx)
        {
          int T_index = T_list_[idx];
          if (static_cast<size_t>(T_index) < T.size())
            {
              T_used[idx] = T[T_index];
//              cout << "T_used[" << idx << "] (from T[" << T_index << "]): "<< endl;
            }
          else
            {
              cerr << "Index " << T_index << " out of bounds for T vector." << endl;
            }
        }

//      x_elem = xp(elem, 0) ; // x=0
      y_elem = xp(elem, 1) ; // y=1
      z_elem = xp(elem, 2) ; // z=2
      y_plus = (y_elem + 1)*u_t/nu;  //xp_paroi(position_base,1);
      z_plus =  (z_elem + 1)*u_t/nu;
      Re_t = u_t/nu;
//      cout << "Re_t_RST=" << Re_t << endl;
      // prediction by neural network
//      start = high_resolution_clock::now();
//      cout << "Making predictions on elem " << elem << "..." << endl;

      b = tbnn->predict_carre(lambda, T_used, y_plus, z_plus, Re_t);

      if (!keps_NN_applied)  // Check if the block has already been executed
        {
          if (apply_keps_NN_)
            {
              tbnn_eps->canal_carre(true);
              tbnn_k->canal_carre(true);
              eps = tbnn_eps->predict_eps_carre(lambda, y_plus, z_plus, Re_t, K_eps(elem, 1));
              k_ec = tbnn_k->predict_k_carre(lambda,  y_plus, z_plus, Re_t, K_eps(elem, 0));
              resu_keps(elem, 0) = k_ec[0];
              resu_keps(elem, 1) = eps[0];
            }
        }

//      end = high_resolution_clock::now();
//      duration = duration_cast<microseconds>(end - start).count();
//      cout << "Time for making predictions " << duration << " microseconds" << endl;
      // add realizability constraints
      int idx_diag[] = {0,3,5};
      int idx_nondiag[] = {1,2,4};

      for (int i=0; i<3; i++ )
        {
          if (b[idx_diag[i]] < -1./3.)
            {
              b[idx_diag[i]] = -1./3.;
            }
          else if (b[idx_diag[i]] > 2./3.)
            {
              b[idx_diag[i]] = 2./3.;
            }
        }

      b[5] = -(b[0] + b[3]);

      for (int i=0; i<3; i++ )
        {
          if (b[idx_nondiag[i]] < -1./2.)
            {
              b[idx_nondiag[i]] = -1./2.;
            }
          else if (b[idx_nondiag[i]] > 1./2.)
            {
              b[idx_nondiag[i]] = 1./2.;
            }
        }

      resu(elem,0,0) = b[0];
      resu(elem,0,1) = b[1];
      resu(elem,0,2) = b[2];
      resu(elem,1,0) = b[1];
      resu(elem,1,1) = b[3];
      resu(elem,1,2) = b[4];
      resu(elem,2,0) = b[2];
      resu(elem,2,1) = b[4];
      resu(elem,2,2) = b[5];

      g0(elem) = 0;
      g1(elem) = 0;
      g2(elem) = 0;
      g3(elem) = 0;
      g4(elem) = 0;
      g5(elem) = 0;
      g6(elem) = 0;
      g7(elem) = 0;
      g8(elem) = 0;
      g9(elem) = 0;
      g10(elem) = 0;

      vector<DoubleTab*> g_map =
      {
        &g0, &g1, &g2, &g3, &g4, &g5, &g6, &g7, &g8, &g9, &g10
      };


      for (size_t idx_g = 0; idx_g < T_list_.size(); idx_g++)
        {
          int t_index = T_list_[idx_g];
          (*g_map[t_index])(elem) = tbnn->get_g_carre()[idx_g];
        }
    }

  // save predictions and g1 values
//  cout << "idt_1=" << idt_1 << endl;
  if ((apply_keps_NN_) && (!keps_NN_applied))
    {
      KEps_true_ = resu_keps;
      keps_NN_applied = true;
    }

  if ((apply_keps_DNS_) && (!keps_DNS_applied))
    {
//      const Domaine_dis_base& domaine_dis_base=eqn_transport_K_Eps_->inconnue().domaine_dis_base();
      const IntTab& les_elems=domaine_VDF.domaine().les_elems();
      int nb_som_elem=les_elems.dimension(1);

      const int nn = eqn_transport_K_Eps_->inconnue().valeurs().dimension(0);
      KEps_true_.resize(nn,2);
      KEps_true_ = 0.;
      assert(KEps_true_.dimension(0) == nelem_);
      for (int ele=0; ele<nelem_; ele++)
        {
          for (int s=0; s<nb_som_elem; s++)
            {
              int sglob=les_elems(ele,s);
              KEps_true_(ele,0)+=KEps_champ->valeurs()(sglob,0);
              KEps_true_(ele,1)+=KEps_champ->valeurs()(sglob,1);
            }
        }
      double inv_nb_som_elem=1./(nb_som_elem);
      KEps_true_*=inv_nb_som_elem;
      keps_DNS_applied = true;
    }

  if ((idt_1 == 0) || (idt_1 % idt_mettre_a_jour == 0))
    {
//      cout << "idt_1 in the loop=" << idt_1 << endl;
      g0_ = g0;
      g1_ = g1;
      g2_ = g2;
      g3_ = g3;
      g4_ = g4;
      g5_ = g5;
      g6_ = g6;
      g7_ = g7;
      g8_ = g8;
      g9_ = g9;
      g10_ = g10;
      if (idt_1==0)
        {
          Cout <<  "g1 updated by NN at time step " << idt_1 << " g1^{n} = g1^{n} "  << endl;
        }

      if (idt_1>0)
        {
          for (int elem=0; elem<nelem_; elem++)
            {g1_(elem) = relax_g1 * g1(elem) + (1-relax_g1)* g1_avant_(elem);}
          Cout <<  "g1 updated by NN at time step " << idt_1 << " g1^{n} = " << relax_g1 << " * g1^{n} " << " + (1-" << relax_g1 << ") * g1^{n-1} " << endl;
        }

      DoubleTab& bij_tab = modele_K_Eps_.valeur().get_bij();
      for (int elem=0; elem<nelem_; elem++)
        {
          bij_tab(elem,0) = resu(elem,0,0); //b11
          bij_tab(elem,1) = resu(elem,0,1); //b12
          bij_tab(elem,2) = resu(elem,0,2); //b13
          bij_tab(elem,3) = resu(elem,1,1); //b22
          bij_tab(elem,4) = resu(elem,1,2); //b23
          bij_tab(elem,5) = resu(elem,2,2); //b33

        }
      bij_ = resu;

    }
  idt_1 = idt_1 + 1;
  g1_avant_ = g1_;

}

void Tenseur_Reynolds_Externe_VDF_Face::Calcul_bij_Nisizima()
{
  DoubleTab resu(nelem_,dimension,dimension);
  DoubleTab g1(nelem_);
  DoubleTab g2(nelem_);
  DoubleTab g3(nelem_);
  DoubleTab g4(nelem_);

  const Domaine_dis& le_dom_dis = eqn_transport_K_Eps_.valeur().domaine_dis();
  const Domaine_Cl_dis& le_dom_Cl_dis = eqn_transport_K_Eps_.valeur().domaine_Cl_dis();
  const Champ_Don ch_visco=ref_cast(Fluide_base,eqn_transport_K_Eps_.valeur().milieu()).viscosite_cinematique();
  const DoubleTab& tab_K_Eps = eqn_transport_K_Eps_->inconnue().valeurs();
  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc=ref_cast(Modele_turbulence_hyd_K_Eps,eqn_transport_K_Eps_->modele_turbulence()).associe_modele_fonction();

  DoubleTrav Fmu;
  Fmu.resize(nelem_);
  mon_modele_fonc.Calcul_Fmu(Fmu,le_dom_dis,le_dom_Cl_dis,tab_K_Eps,ch_visco);
  for (int elem=0; elem<nelem_; elem++)
    {

      g1(elem) = -0.09;
      g2(elem) = -0.0324 * Fmu(elem);
      g3(elem) = -0.1368 * Fmu(elem);
      g4(elem) = 0.1872 * Fmu(elem);

      for (int i=0; i<Objet_U::dimension; i++)
        {
          for (int j=0; j<Objet_U::dimension; j++)
            {
              resu(elem,i,j) = g1(elem)  * Fmu(elem) *  T1_etoile_(elem,i,j)
                               + g2(elem) *  T2_etoile_(elem,i,j)
                               + g3(elem) *  T3_etoile_(elem,i,j) + g4(elem) *  T4_etoile_(elem,i,j);
            }
        }
    }
  g1_ = g1;
  g2_ = g2;
  g3_ = g3;
  g4_ = g4;

  DoubleTab& bij_tab = modele_K_Eps_.valeur().get_bij();
  for (int elem=0; elem<nelem_; elem++)
    {
      bij_tab(elem,0) = resu(elem,0,0);
      bij_tab(elem,1) = resu(elem,0,1);
      bij_tab(elem,2) = resu(elem,0,2);
      bij_tab(elem,3) = resu(elem,1,1);
      bij_tab(elem,4) = resu(elem,1,2);
      bij_tab(elem,5) = resu(elem,2,2);
    }
  bij_ = resu;
}

DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::Calcul_bij_NL_Nisizima(DoubleTab& resu)
{
  Calcul_bij_Nisizima();
  resu = 0;

  for (int elem=0; elem<nelem_; elem++)
    {

      for (int i=0; i<Objet_U::dimension; i++)
        {
          for (int j=0; j<Objet_U::dimension; j++)
            {
              resu(elem,i,j) = g2_(elem) *  T2_etoile_(elem,i,j)
                               + g3_(elem) *  T3_etoile_(elem,i,j) + g4_(elem) *  T4_etoile_(elem,i,j);
            }
        }
    }
  return resu;
}

std::vector<double> Tenseur_Reynolds_Externe_VDF_Face::read_file(const std::string& filename)
{
  std::ifstream file(filename);
  std::vector<double> vec;
  double num;
  while (file >> num)
    {
      vec.push_back(num);
    }
  return vec;
}
