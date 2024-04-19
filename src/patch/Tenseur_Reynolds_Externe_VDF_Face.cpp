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

//#define NON_LINEAR_BIJ

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


//// readOn
//

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
  Motcles les_mots(3);
  {
    les_mots[0] = "nom_fichier";
    les_mots[1] = "canal_plan";
    les_mots[2] = "canal_carre";
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
            cout << "Using neural network trained on plane channel flow" << endl;
            //Cerr << "The value of is_canal_plan_ bool is: " << static_cast<unsigned int>(f) << finl; // Cast CANAL_PLAN_ to unsigned int
            //Cerr << "The value of is_canal_carre_ bool is: " << static_cast<unsigned int>(tbnn->is_canal_carre_) << finl; // Cast CANAL_PLAN_ to unsigned int
            break;
          }
        case 2 : // "canal_carre"
          {
            tbnn->canal_carre(true);
            cout << "Using neural network trained on square duct flow" << endl;
            //Cerr << "The value of is_canal_plan_ bool is: " << static_cast<unsigned int>(tbnn->is_canal_plan_) << finl; // Cast CANAL_PLAN_ to unsigned int
            //Cerr << "The value of is_canal_carre_ bool is: " << static_cast<unsigned int>(tbnn->is_canal_carre_) << finl; // Cast CANAL_PLAN_ to unsigned int
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

const DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::get_g1() const
{
  return g1_;
}

void Tenseur_Reynolds_Externe_VDF_Face::mettre_a_jour(double temps)
{
  int nb_elem_tot = le_dom_VDF.valeur().nb_elem_tot();
  int nb_faces_tot = le_dom_VDF.valeur().nb_faces_tot();
  const IntTab& face_vois = le_dom_VDF.valeur().face_voisins();
  const DoubleVect& volumes = le_dom_VDF.valeur().volumes();

  DoubleTab valeurs_source(nb_elem_tot,dimension);
  valeurs_source = 0;

  if (tbnn->is_canal_carre_)
    {
      Calcul_RSLambda();
    }

  DoubleTab tenseur_reynolds_elements(nb_elem_tot,dimension,dimension);
  tenseur_reynolds_elements = 0;

  Calcul_Tenseur_Reynolds_NL( tenseur_reynolds_elements );

  DoubleTab tenseur_reynolds_faces(nb_faces_tot,dimension,dimension);
  tenseur_reynolds_faces = 0.;

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
          valeurs_source(elem,i) += div;
        }
    }

  valeurs_source.echange_espace_virtuel();
  la_source->valeurs( ) = valeurs_source;

}

void Tenseur_Reynolds_Externe_VDF_Face::completer()
{
  Source_base::completer();

  la_source.DERIV(Champ_Don_base)::typer( "Champ_Fonc_xyz" );
}

void Tenseur_Reynolds_Externe_VDF_Face::Calcul_RSLambda()
{
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const Domaine_Cl_VDF& domaine_Cl_VDF = le_dom_Cl_VDF.valeur();

  int nb_elem_tot=domaine_VDF.nb_elem_tot();
  const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF,eqn_NS_->inconnue().valeur() );
  assert (vitesse.valeurs().line_size() == 1);
  DoubleTab gij(nb_elem_tot,dimension,dimension, vitesse.valeurs().line_size());
  ref_cast_non_const(Champ_Face_VDF,vitesse).calcul_duidxj( vitesse.valeurs(),gij,domaine_Cl_VDF );

  DoubleTab lambda_1(nb_elem_tot);
  DoubleTab lambda_2(nb_elem_tot);
  DoubleTab lambda_3(nb_elem_tot);
  DoubleTab lambda_4(nb_elem_tot);
  DoubleTab lambda_5(nb_elem_tot);

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

  for (int elem=0; elem<nelem_; elem++)
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

                L1Id(elem,i,j) += lambda_1(elem)/3.;
                L2Id(elem,i,j) += lambda_2(elem)/3.;
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

                L4Id(elem,i,j) += lambda_4(elem)*2./3.;
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

                L5Id(elem,i,j) += lambda_5(elem)*2./3.;
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

  T1_etoile_  = S_etoile;

  T2_etoile_  = SR;
  T2_etoile_ -= RS;

  T3_etoile_  = S2;
  T3_etoile_ -= L1Id;

  T4_etoile_  = R2;
  T4_etoile_ -= L2Id;

  T5_etoile_  = RS2;
  T5_etoile_ -= S2R;

  T6_etoile_  = R2S;
  T6_etoile_ += SR2;
  T6_etoile_ -= L4Id;

  T7_etoile_  = RSR2;
  T7_etoile_ -= R2SR;

  T8_etoile_  = SRS2;
  T8_etoile_ -= S2RS;

  T9_etoile_  = R2S2;
  T9_etoile_ += S2R2;
  T9_etoile_ -= L5Id;

  T10_etoile_  = RS2R2;
  T10_etoile_ -= R2S2R;
}

DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::Calcul_Tenseur_Reynolds(DoubleTab& resu)
{
  DoubleTab bij;

  if (tbnn->is_canal_plan_)
    {
      bij = Calcul_bij_TBNN(resu);
    }
  if (tbnn->is_canal_carre_)
    {
      bij = Calcul_bij_TBNN_carre(resu);
    }

  const DoubleTab& K_eps = eqn_transport_K_Eps_->inconnue().valeurs();

  for (int elem=0; elem<nelem_; elem++)
    {
      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            resu(elem,i,j) = bij(elem,i,j) ;

            if (i==j)
              {
                resu(elem,i,j) += 1./3.;
              }

            resu(elem,i,j) *= 2. * K_eps(elem,0);
          }
    }


  return resu;
}

DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::Calcul_Tenseur_Reynolds_NL(DoubleTab& resu)
{
  if (tbnn->is_canal_plan_)
    {
      Calcul_bij_NL_TBNN(resu);
    }
  if (tbnn->is_canal_carre_)
    {
      Calcul_bij_NL_TBNN_carre(resu);
    }

  const DoubleTab& K_eps = eqn_transport_K_Eps_->inconnue().valeurs();

  for (int elem=0; elem<nelem_; elem++)
    {
      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            resu(elem,i,j) *= 2. * K_eps(elem,0);

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
      alpha = compute_alpha(K_eps(elem, 0), K_eps(elem,1), gij(elem,0,1,0) );
      resu_NL(elem,0,1) -= 0.5 * g1_(elem) * alpha;
      resu_NL(elem,1,0) -= 0.5 * g1_(elem) * alpha;
    }

  return resu_NL;
}

DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::Calcul_bij_NL_TBNN_carre(DoubleTab& resu)
{

  DoubleTab& resu_NL = Calcul_bij_TBNN_carre(resu);

  for (int elem=0; elem<nelem_; elem++)
    {
      for (int i=0; i<Objet_U::dimension; i++)
        {
          for (int j=0; j<Objet_U::dimension; j++)
            {
# ifdef NON_LINEAR_BIJ
              resu_NL(elem,i,j) -= abs(g1_(elem)) * T1_etoile_(elem,i,j);
# else
              resu_NL(elem,i,j) = 0;
# endif
            }
        }
    }

  return resu_NL;
}
DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::Calcul_bij_TBNN(DoubleTab& resu)
{
  vector<double> b;
  double alpha, y_plus, Re_t;
  double y_maille_paroi =-1, y_elem;
  //  vector<vector<double>> T;
  DoubleTab g1(nelem_);
  const DoubleTab& K_eps = eqn_transport_K_Eps_->inconnue().valeurs();
  double u_t;

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

              //u_t = sqrt(fabs(nu* dudy_paroi)+1e-20);
              u_t = sqrt(nu* dudy_paroi)+1e-40;
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

      alpha = compute_alpha(K_eps(elem, 0), K_eps(elem,1), gij(elem,0,1,0) ) + 1e-40;
      y_plus = compute_y_plus(y_plus_wall[position_base],y_maille_paroi, y_elem);
      Re_t = compute_Re_t(y_plus_wall[position_base],y_maille_paroi);


      // prediction by neural network
      //Cerr << " Testing PCF TRE_VDF_Face.cpp line 1213 " << finl;
      b = tbnn->predict(alpha, y_plus, Re_t);

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

DoubleTab& Tenseur_Reynolds_Externe_VDF_Face::Calcul_bij_TBNN_carre(DoubleTab& resu)
{
  vector<double> b, lambda;
  double y_plus, z_plus, Re_t;
  double y_maille_paroi =-1;
  double z_maille_paroi=-1;
  vector<vector<double>> T;
  DoubleTab g1(nelem_);
  double u_t;

  double dudy_paroi, dudz_paroi ;

  //  std::vector<double> A = read_file("/volatile/catB/dr274584/stage_sanae/dr_essais/dr_dossier/TBNN/NN_DNS/DNS_1000/A_1000.dat");
  //  std::vector<double> C = read_file("/volatile/catB/dr274584/stage_sanae/dr_essais/dr_dossier/TBNN/NN_DNS/DNS_1000/C_1000.dat");
  //  std::vector<double> B0 = read_file("/volatile/catB/dr274584/stage_sanae/dr_essais/dr_dossier/TBNN/NN_DNS/DNS_1000/B0_1000.dat");
  //  std::vector<double> B1 = read_file("/volatile/catB/dr274584/stage_sanae/dr_essais/dr_dossier/TBNN/NN_DNS/DNS_1000/B1_1000.dat");
  //  std::vector<double> B2 = read_file("/volatile/catB/dr274584/stage_sanae/dr_essais/dr_dossier/TBNN/NN_DNS/DNS_1000/B2_1000.dat");
  //  std::vector<double> B3 = read_file("/volatile/catB/dr274584/stage_sanae/dr_essais/dr_dossier/TBNN/NN_DNS/DNS_1000/B3_1000.dat");
  //  std::vector<double> B4 = read_file("/volatile/catB/dr274584/stage_sanae/dr_essais/dr_dossier/TBNN/NN_DNS/DNS_1000/B4_1000.dat");
  //  std::vector<double> B5 = read_file("/volatile/catB/dr274584/stage_sanae/dr_essais/dr_dossier/TBNN/NN_DNS/DNS_1000/B5_1000.dat");
  //  std::vector<double> c_mu_DNS = read_file("/volatile/catB/dr274584/stage_sanae/dr_essais/dr_dossier/TBNN/Cas_1000/c_mu_DNS_1000.dat");

  //const Champ_base& y_plus_champ = eqn_NS_.valeur().get_champ("y_plus");
  //DoubleTab y_plus_elem;
  //const DoubleTab& positions=le_dom_VDF.valeur().xp();
  //IntVect les_polys(positions.dimension(0));
  //le_dom_VDF.valeur().domaine().chercher_elements(positions, les_polys);
  //y_plus_champ.valeur_aux_elems(positions, les_polys, y_plus_elem);

  //const DoubleTab& y_plus_elem = eqn_NS_->get_champ("y_plus").valeurs();

  double x_elem, y_elem, z_elem;
  int position_base = -1;
  int position_gauche = -1;

  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const Domaine_Cl_VDF& domaine_Cl_VDF = le_dom_Cl_VDF.valeur();

  int nb_elem_tot=domaine_VDF.nb_elem_tot();
  const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF,eqn_NS_->inconnue().valeur() );
  assert (vitesse.valeurs().line_size() == 1);
  vector<double> y_plus_wall, z_plus_wall, x1, z1, x2, y2, Re_true;

  //const DoubleTab& vitesse_faces = eqn_NS_->inconnue().valeurs() ;
  //double vitesse_sx, vitesse_dx;

  DoubleTab gij(nb_elem_tot,dimension,dimension, vitesse.valeurs().line_size());
  ref_cast_non_const(Champ_Face_VDF,vitesse).calcul_duidxj( vitesse.valeurs(),gij,domaine_Cl_VDF );

  double nu = eqn_NS_->fluide().viscosite_cinematique().valeurs()[0];

  int direction = -1;

  // Boucle sur les bords pour traiter les conditions aux limites
  for (int n_bord=0; n_bord<domaine_VDF.domaine().nb_bords(); n_bord++)
    {
      // Si face de Dirichlet (les parois) on calcule y+ ou z+

      const Cond_lim& la_cl = domaine_Cl_VDF.les_conditions_limites(n_bord);

      int elem_paroi;

      if ( (sub_type(Dirichlet,la_cl.valeur()))   ||  (sub_type(Dirichlet_homogene,la_cl.valeur())))
        {

          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());


          int ndeb = le_bord.num_premiere_face();
          int nfaces = le_bord.nb_faces();
          int nfin = ndeb + nfaces;
          //Cout << "ndeb " << ndeb << finl;

          elem_paroi = le_dom_VDF->face_voisins(ndeb, 1);
          double coordy =le_dom_VDF->xp(elem_paroi, 1);
          double coordz =le_dom_VDF->xp(elem_paroi, 2);

          if (fabs(coordy - le_dom_VDF->xv(ndeb, 1)) < 1e-8 )
            {
              direction = 2;
              z_plus_wall.resize(nfaces);
              x2.resize(nfaces);
              y2.resize(nfaces);
              //Cout << "nfaces direction 2 " << nfaces << finl;
            }
          if (fabs(coordz - le_dom_VDF->xv(ndeb, 2)) < 1e-8 )
            {
              direction = 1;
              y_plus_wall.resize(nfaces);
              Re_true.resize(nfaces);
              x1.resize(nfaces);
              z1.resize(nfaces);
              //Cout << "nfaces direction 1 " << nfaces << finl;
            }

          for ( int num_face=ndeb; num_face<nfin; num_face++)
            {
              elem_paroi = le_dom_VDF->face_voisins(num_face, 1);//1 means that I want the elem above the face (which is the bottom paroi so ok)

              if (direction == 1)
                {
                  // compute y+ from the y wall
                  y_maille_paroi =le_dom_VDF->xp(elem_paroi, 1) ;
                  //Cerr << " y_maille_paroi= "<< y_maille_paroi << finl; VERIFICATO
                  dudy_paroi = gij(elem_paroi,0,1,0);
                  u_t = sqrt( abs(nu* dudy_paroi))+1e-40;

                  y_plus_wall[num_face-ndeb] = abs(y_maille_paroi) * u_t / nu;
                  Re_true [num_face-ndeb] =  u_t / nu;
                  x1[num_face-ndeb] =le_dom_VDF->xv(num_face, 0) ;
                  z1[num_face-ndeb] =le_dom_VDF->xv(num_face, 2) ;
                }

              if ( direction == 2 )
                {
                  // compute y+ from the z wall
                  z_maille_paroi =le_dom_VDF->xp(elem_paroi, 2) ;
                  dudz_paroi = gij(elem_paroi,0,2,0);
                  u_t = sqrt( abs(nu* dudz_paroi))+1e-40;

                  z_plus_wall[num_face-ndeb] = abs(z_maille_paroi) * u_t / nu;
                  x2[num_face-ndeb] =le_dom_VDF->xv(num_face, 0) ;
                  y2[num_face-ndeb] =le_dom_VDF->xv(num_face, 1) ;
                }

            }
        }
    }

  //
  //  init of lambda and T arrays
  lambda.resize(5);
  T.resize(11);
  for (int i=0; i<11; i++)
    T[i].resize(6);
  T[0][0] = 1./6.;
  T[0][1] = 0.;
  T[0][2] = 0.;
  T[0][3] = 1./6.;
  T[0][4] = 0.;
  T[0][5] = -1./3.;


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

      x_elem =le_dom_VDF->xp(elem, 0) ; // x=0
      y_elem =le_dom_VDF->xp(elem, 1) ; // y=1
      z_elem =le_dom_VDF->xp(elem, 2) ; // z=2

      for (unsigned int i = 0; i < y_plus_wall.size(); i++)
        if ( fabs (x_elem-x1[i]) <1e-8 and fabs (z_elem-z1[i]) <1e-8)
          {
            position_base = i;
            //Cout << "position_base " << position_base << finl;
            break;
          }

      for (unsigned int i = 0; i < z_plus_wall.size(); i++)
        if ( fabs (x_elem-x2[i]) <1e-8 and fabs (y_elem-y2[i]) <1e-8)
          {
            position_gauche = i;
            //Cout << "position_gauche " << position_gauche << finl;
            break;
          }

      y_plus = compute_y_plus(y_plus_wall[position_base],y_maille_paroi, y_elem);
      z_plus = compute_y_plus(z_plus_wall[position_gauche],z_maille_paroi, z_elem);
      Re_t = compute_Re_t(y_plus_wall[position_base],y_maille_paroi);
      //Cerr << "Re_t is: "<< Re_t << finl;

      // prediction by neural network
      b = tbnn->predict_carre(lambda, T, y_plus, z_plus, Re_t);
      resu(elem,0,0) = b[0];
      resu(elem,0,1) = b[1];
      resu(elem,0,2) = b[2];
      resu(elem,1,0) = b[1];
      resu(elem,1,1) = b[3];
      resu(elem,1,2) = b[4];
      resu(elem,2,0) = b[2];
      resu(elem,2,1) = b[4];
      resu(elem,2,2) = b[5];

      g1(elem) = tbnn->get_g1_carre();
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

double Tenseur_Reynolds_Externe_VDF_Face::compute_alpha( double k, double eps, double dudy)
{
  double alpha;

  alpha = k / (eps + 1e-15) * dudy;

  return alpha;
}


double Tenseur_Reynolds_Externe_VDF_Face::compute_y_plus( double y_plus_wall, double h_maille_paroi, double h_elem)
{
  double y_plus;

  y_plus = abs(y_plus_wall / h_maille_paroi * h_elem);

  return y_plus;

}

double Tenseur_Reynolds_Externe_VDF_Face::compute_Re_t(double y_plus_wall, double h_maille_paroi)
{
  double Re_t;

  Re_t = abs(y_plus_wall / h_maille_paroi * 1.0);

  return Re_t;

}
