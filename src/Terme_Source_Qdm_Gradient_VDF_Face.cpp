/****************************************************************************
* Copyright (c) 2022, CEA
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

#include <Terme_Source_Qdm_Gradient_VDF_Face.h>
#include <Champ_Uniforme.h>
#include <Domaine_Cl_dis.h>
#include <Domaine_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Neumann_sortie_libre.h>
#include <Milieu_base.h>
#include <Dirichlet.h>
#include <Dirichlet_homogene.h>
#include <Symetrie.h>
#include <Periodique.h>
#include <Navier_Stokes_std.h>


Implemente_instanciable(Terme_Source_Qdm_Gradient_VDF_Face,"Source_Qdm_Gradient_VDF_Face",Source_base);

//// printOn

Sortie& Terme_Source_Qdm_Gradient_VDF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

//// readOn
//

Entree& Terme_Source_Qdm_Gradient_VDF_Face::readOn(Entree& s )
{
  s >> la_source;
  if (la_source->nb_comp() != dimension)
    {
      Cerr << "Erreur a la lecture du terme source de type " << que_suis_je() << finl;
      Cerr << "le champ source doit avoir " << dimension << " composantes" << finl;
      exit();
    }
  return s;
}

void Terme_Source_Qdm_Gradient_VDF_Face::associer_domaines(const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis)
{
  le_dom_VDF = ref_cast(Domaine_VDF, domaine_dis.valeur());
  le_dom_Cl_VDF = ref_cast(Domaine_Cl_VDF, domaine_Cl_dis.valeur());
}

void Terme_Source_Qdm_Gradient_VDF_Face::associer_pb(const Probleme_base& )
{
  ;
}

void Terme_Source_Qdm_Gradient_VDF_Face::ajouter_blocs(matrices_t matrices, DoubleTab& resu, const tabs_t& semi_impl) const
{
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  //const Domaine_Cl_VDF& domaine_Cl_VDF = le_dom_Cl_VDF.valeur();
  //const IntTab& face_voisins = domaine_VDF.face_voisins();
//  const IntVect& orientation = domaine_VDF.orientation();
//  const DoubleVect& porosite_surf = equation().milieu().porosite_face();
  // const DoubleVect& face_surfaces = domaine_VDF.face_surfaces();
  //const DoubleVect& volumes_entrelaces = domaine_VDF.volumes_entrelaces();

  // pour interpoler la valeur dans l'élément
  //const DoubleTab& xp = domaine_VDF.xp();
  //const IntTab& les_elems = le_dom_VDF->domaine().les_elems();
  //const int nb_som_elem=les_elems.dimension(1);
  //const int nb_elem = le_dom_VDF->nb_elem();
  //DoubleTab valeurs_source(nb_elem,dimension);


  //Boucle que sur les faces internes
  for (int num_face =domaine_VDF.premiere_face_int(); num_face<domaine_VDF.nb_faces(); num_face++)
    {



//  valeurs_source = 0.;
//  double surf_check;
//
//  for (int ele=0; ele<nb_elem; ele++)
//    {
//      int face_glob_1_x = domaine_VDF.elem_faces(ele,0);
//      int face_glob_2_x = domaine_VDF.elem_faces(ele,dimension);
//
//      for (int i = 0; i < dimension; i++)
//        {
//
//          //we must verify that the opposite face is face+dim
//          int face_glob_1 = domaine_VDF.elem_faces(ele,i);
//          int face_glob_2 = domaine_VDF.elem_faces(ele,i + dimension);
////          double vol = domaine_VDF.volumes(ele);
//
//
//          //this check is not necessary anymore bc the code works as expected
//          surf_check = fabs( domaine_VDF.face_surfaces(face_glob_2) - domaine_VDF.face_surfaces(face_glob_1))/ domaine_VDF.face_surfaces(face_glob_1);
//          if ( surf_check >0.05)
//            cerr<< "Not matching surfaces, error in percentage:"<< surf_check<< endl;
//
//
//          double valeurs_source_face_1 = 0;
//          double valeurs_source_face_2 = 0;
//
//          //compute value at one face and the opposite
//          for (int somloc=0; somloc<4; somloc++)
//            {
//              int sglob1 = domaine_VDF.face_sommets( face_glob_1, somloc);
//              int sglob2 = domaine_VDF.face_sommets( face_glob_2, somloc);         //face_loc_2 = face_glob_1 + dim
//
//              valeurs_source_face_1 += 0.25 * la_source->valeurs()(sglob1,i);  //i is the dimension taken into account
//              valeurs_source_face_2 += 0.25 * la_source->valeurs()(sglob2,i);
//            }
//
//          //compute gradient
//          // check if 2-1 or viceversa
//          // instead of doing ./disance_btw_faces*volume_ele we just *surface to optimise the code
//          double grad_times_vol = (valeurs_source_face_2 - valeurs_source_face_1) * domaine_VDF.face_surfaces(face_glob_1);
//
////          if (resu(face_glob_1)!= 0)
////            Cout << "Grad_vol impact on resu : " << fabs(grad_times_vol/resu(face_glob_1))*100  << "%" << finl;
////
////          else
////        	  Cout << "resu (face_glob) ==0" << finl;
//
//          // 0.5 assigned to face 1 and 0.5 to face 2
//          // - because Sanae told me that she checked it and it was - on her code
//          resu(face_glob_1_x) -= 0.5 * grad_times_vol;
//          resu(face_glob_2_x) -= 0.5 * grad_times_vol;
//


    }


}

DoubleTab& Terme_Source_Qdm_Gradient_VDF_Face::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}

void Terme_Source_Qdm_Gradient_VDF_Face::mettre_a_jour(double temps)
{
  la_source->mettre_a_jour(temps);
}




////      for (int som=0; som < nb_som_elem; som++)
////        {
////          int sglob=les_elems(ele,som);
////          valeurs_source(ele,i) += la_source->valeurs()(sglob,i);
////        }
////	  double inv_nb_som_elem = 1./(nb_som_elem);
////	  valeurs_source*=inv_nb_som_elem;
//
//  // Boucle sur les bords pour traiter les conditions aux limites
//  for (int n_bord = 0; n_bord < domaine_VDF.nb_front_Cl(); n_bord++)
//    {
//      const Cond_lim& la_cl = domaine_Cl_VDF.les_conditions_limites(n_bord);
//      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
//      const int ndeb = le_bord.num_premiere_face(), nfin = ndeb + le_bord.nb_faces();
//
//      if ( sub_type(Neumann_sortie_libre,la_cl.valeur()) )
//        {
////          for (int num_face=ndeb; num_face<nfin; num_face++)
////            {
////              const double coef = face_surfaces(num_face)*porosite_surf(num_face);
////              int ncomp = orientation(num_face);
////              int elem1 = face_voisins(num_face,0);
////
////              if (elem1 != -1)
////                resu(num_face) -= valeurs_source(elem1,ncomp)*coef;
////              else
////                {
////                  int elem2 = face_voisins(num_face,1);
////                  resu(num_face) += valeurs_source(elem2,ncomp)*coef;
////                }
////            }
//          ;
//        }
//      else if (sub_type(Symetrie,la_cl.valeur()))
//        {
//          ;
//        }
//      else if ( (sub_type(Dirichlet,la_cl.valeur()))
//                ||
//                (sub_type(Dirichlet_homogene,la_cl.valeur()))
//              )
//        ;
//      else if (sub_type(Periodique,la_cl.valeur()))
//        {
//
//          for (int num_face=ndeb; num_face<nfin; num_face++)
//            {
//
//              double vol = volumes_entrelaces(num_face)*porosite_surf(num_face);
//              int ncomp = orientation(num_face);
//              const int n0 = face_voisins(num_face, 0), n1 = face_voisins(num_face, 1);
//              // the distANCE BTW THE TWO CENTER OF GRAVITY CHANGES
//              // Il faut changer 1 si on change le domaine (géneraliser)
//              double dist = 1 - fabs(xp(n1,ncomp)- xp(n0,ncomp));
//              // here 0.5* because the for cycle counts the same face at right end and at left end
//              resu(num_face) += 0.5 * (valeurs_source(n1, ncomp) - valeurs_source(n0, ncomp)) / dist * vol;
//
//
////              const int n0 = face_voisins(num_face, 0), n1 = face_voisins(num_face, 1);
////              const double coef = face_surfaces(num_face) * porosite_surf(num_face) ;
////              int ncomp = orientation(num_face);
////              resu(num_face) += coef * (valeurs_source(n1, ncomp) - valeurs_source(n0, ncomp));
//            }
//        }
//    }
//
//  // Boucle sur les faces internes
//
//  for (int num_face =domaine_VDF.premiere_face_int(); num_face<domaine_VDF.nb_faces(); num_face++)
//    {
//      double vol = volumes_entrelaces(num_face)*porosite_surf(num_face);
//      int ncomp = orientation(num_face);
//      const int n0 = face_voisins(num_face, 0), n1 = face_voisins(num_face, 1);
//      double dist = fabs(xp(n1,ncomp)- xp(n0,ncomp));
//      resu(num_face) += (valeurs_source(n1, ncomp) - valeurs_source(n0, ncomp)) / dist * vol;
////      Cout << (valeurs_source(n1, ncomp) - valeurs_source(n0, ncomp)) / dist * vol << finl;
////      Cout << "resu(num_face) = " << resu(num_face) << finl;
//    }

