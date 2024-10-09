/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Modele_turbulence_hyd_K_Eps.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Eps.h>
#include <Probleme_base.h>
#include <Debog.h>
#include <Modifier_nut_pour_fluide_dilatable.h>
#include <Schema_Temps_base.h>
#include <Schema_Temps.h>
#include <stat_counters.h>
#include <Modele_turbulence_scal_base.h>
#include <Param.h>
#include <communications.h>
#include <Fluide_base.h>
#include <Champ_Face_VDF.h>
#include <TRUSTTrav.h>
#include <Champ_Uniforme.h>
#include <TRUSTTab_parts.h>
#include <Champ_Inc_P0_base.h>
#include <Tenseur_Reynolds_Externe_VDF_Face.h>
#include <Domaine_VDF.h>
#include <Domaine_Cl_VDF.h>

Implemente_instanciable(Modele_turbulence_hyd_K_Eps,"Modele_turbulence_hyd_K_Epsilon",Mod_turb_hyd_RANS);
// XD k_epsilon mod_turb_hyd_rans k_epsilon -1 Turbulence model (k-eps).

/*! @brief Ecrit le type de l'objet sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Modele_turbulence_hyd_K_Eps::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}


/*! @brief Simple appel a Mod_turb_hyd_RANS::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Modele_turbulence_hyd_K_Eps::readOn(Entree& s )
{
  return Mod_turb_hyd_RANS::readOn(s);
}

void Modele_turbulence_hyd_K_Eps::set_param(Param& param)
{
  Mod_turb_hyd_RANS::set_param(param);
  param.ajouter_non_std("Transport_K_Epsilon",(this),Param::REQUIRED); // XD_ADD_P transport_k_epsilon Keyword to define the (k-eps) transportation equation.
  param.ajouter_non_std("Modele_Fonc_Bas_Reynolds",(this)); // XD_ADD_P modele_fonction_bas_reynolds_base This keyword is used to set the bas Reynolds model used.
  param.ajouter("CMU",&LeCmu); // XD_ADD_P double Keyword to modify the Cmu constant of k-eps model : Nut=Cmu*k*k/eps Default value is 0.09
  param.ajouter("PRANDTL_K",&Prandtl_K); // XD_ADD_P double Keyword to change the Prk value (default 1.0).
  param.ajouter("PRANDTL_EPS",&Prandtl_Eps); // XD_ADD_P double Keyword to change the Pre value (default 1.3).
  //added
  param.ajouter("CMU_champ",&LeCmu_champ); // XD_ADD_P field_base todo
}

int Modele_turbulence_hyd_K_Eps::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="Transport_K_Epsilon")
    {
      eqn_transp_K_Eps().associer_modele_turbulence(*this);
      is >> eqn_transp_K_Eps();
      return 1;
    }
  else if (mot=="Modele_Fonc_Bas_Reynolds")
    {
      Cerr << "Lecture du modele bas reynolds associe " << finl;
      mon_modele_fonc.associer_eqn(eqn_transp_K_Eps());
      is >> mon_modele_fonc;
      Cerr << "mon_modele_fonc.que_suis_je() avant discretisation " << mon_modele_fonc.que_suis_je() << finl;
      mon_modele_fonc.valeur().discretiser();
      Cerr << "mon_modele_fonc.que_suis_je() " << mon_modele_fonc.valeur().que_suis_je() << finl;
      mon_modele_fonc.valeur().lire_distance_paroi();
      return 1;
    }
  else
    return Mod_turb_hyd_RANS::lire_motcle_non_standard(mot,is);
}

/*! @brief Discretise le modele de turbulence.
 *
 */
void Modele_turbulence_hyd_K_Eps::discretiser()
{

  Mod_turb_hyd_RANS::discretiser();

  const Discretisation_base& dis=ref_cast(Discretisation_base, mon_equation->discretisation());
  dis.discretiser_champ("champ_elem",mon_equation->domaine_dis().valeur(),"Cmu","",1,mon_equation->schema_temps().temps_courant(),Cmu_);
  champs_compris_.ajoute_champ(Cmu_);

  dis.discretiser_champ("champ_elem",mon_equation->domaine_dis().valeur(),"Lambda1","",1,mon_equation->schema_temps().temps_courant(),lambda1_);
  champs_compris_.ajoute_champ(lambda1_);

  dis.discretiser_champ("champ_elem",mon_equation->domaine_dis().valeur(),"Lambda2","",1,mon_equation->schema_temps().temps_courant(),lambda2_);
  champs_compris_.ajoute_champ(lambda2_);

  dis.discretiser_champ("champ_elem",mon_equation->domaine_dis().valeur(),"Lambda3","",1,mon_equation->schema_temps().temps_courant(),lambda3_);
  champs_compris_.ajoute_champ(lambda3_);

  dis.discretiser_champ("champ_elem",mon_equation->domaine_dis().valeur(),"Lambda4","",1,mon_equation->schema_temps().temps_courant(),lambda4_);
  champs_compris_.ajoute_champ(lambda4_);

  dis.discretiser_champ("champ_elem",mon_equation->domaine_dis().valeur(),"Lambda5","",1,mon_equation->schema_temps().temps_courant(),lambda5_);
  champs_compris_.ajoute_champ(lambda5_);


  Noms noms(6);
  Noms unit(6);
  noms[0]="bij_00";
  noms[1]="bij_01";
  noms[2]="bij_02";
  noms[3]="bij_11";
  noms[4]="bij_12";
  noms[5]="bij_22";

  unit[0]="";
  unit[1]="";
  unit[2]="";
  unit[3]="";
  unit[4]="";
  unit[5]="";

  dis.discretiser_champ("champ_elem",mon_equation->domaine_dis().valeur(),multi_scalaire,noms,unit,6,mon_equation->schema_temps().temps_courant(),bij_);
  bij_.valeur().nommer("bij");
  champs_compris_.ajoute_champ(bij_);

  Noms noms2(6);
  noms2[0]="bij_NL_00";
  noms2[1]="bij_NL_01";
  noms2[2]="bij_NL_02";
  noms2[3]="bij_NL_11";
  noms2[4]="bij_NL_12";
  noms2[5]="bij_NL_22";

  dis.discretiser_champ("champ_elem",mon_equation->domaine_dis().valeur(),multi_scalaire,noms2,unit,6,mon_equation->schema_temps().temps_courant(),bij_NL_);
  bij_NL_.valeur().nommer("bij_NL");
  champs_compris_.ajoute_champ(bij_NL_);
}

/*! @brief Calcule la viscosite turbulente au temps demande.
 *
 * @param (double temps) le temps auquel il faut calculer la viscosite
 * @return (Champ_Fonc&) la viscosite turbulente au temps demande
 * @throws erreur de taille de visco_turb_K_eps
 */
Champ_Fonc& Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente(double temps)
{
  if (LeCmu_champ.non_nul() && !is_initialized)
    {

      const Domaine_dis_base& domaine_dis_base=eqn_transp_K_Eps().inconnue()->domaine_dis_base();
      const IntTab& les_elems=domaine_dis_base.domaine().les_elems();
      int nb_som_elem=les_elems.dimension(1);
      int nb_elem=domaine_dis_base.nb_elem();

      const int nn = eqn_transp_K_Eps().inconnue()->valeurs().dimension(0);
      LeCmu_tab.resize(nn,1);
      LeCmu_tab = 0.;
      assert(LeCmu_tab.dimension(0) == nb_elem);
      for (int ele=0; ele<nb_elem; ele++)
        {
          for (int s=0; s<nb_som_elem; s++)
            {
              int sglob=les_elems(ele,s);
              LeCmu_tab(ele,0)+=LeCmu_champ->valeurs()(sglob,0);
            }
        }
      double inv_nb_som_elem=1./(nb_som_elem);
      LeCmu_tab*=inv_nb_som_elem;

      is_initialized = 1;
    }

  const Champ_base& chK_Eps=eqn_transp_K_Eps().inconnue().valeur();
  Nom type=chK_Eps.que_suis_je();
  const DoubleTab& tab_K_Eps = chK_Eps.valeurs();
  DoubleTab& visco_turb =  la_viscosite_turbulente.valeurs();
  DoubleTrav Fmu;

//  K_Eps(i,0) = K au noeud i
//  K_Eps(i,1) = Epsilon au noeud i
//  int n = tab_K_Eps.dimension(0);
  int n = tab_K_Eps.dimension(0);
  if (n<0)
    {
      if (sub_type(Champ_Inc_P0_base, chK_Eps))
        n = eqn_transp_K_Eps().domaine_dis().domaine().nb_elem();
      else
        {
          Cerr << "Unsupported K_Eps field in Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente" << finl;
          Process::exit(-1);
        }
    }

  if (tenseur_de_Reynolds_externe_)
    {
      // Cerr<< " On utilise un Cmu non constant donne par le reseau de neurones (qui supplante eventuellement le Cmu non constant du modele fonctionnel) "<< finl;
      const DoubleTab& g1 = get_source_tenseur_de_Reynolds_NN( ).get_g1( );
      for (int i=0; i<n; i++)
        {
          Cmu_.valeurs()[i] = abs( g1( i ));
        }
    }
  if (LeCmu_champ.non_nul() && is_initialized)
    {
      for (int i=0; i<n; i++)
        {
          Cmu_.valeurs()[i] = LeCmu_tab(i,0);
        }
    }

  int is_modele_fonc=(mon_modele_fonc.non_nul());
  if (is_modele_fonc)
    {
      const Domaine_dis& le_dom_dis = eqn_transp_K_Eps().domaine_dis();
      const Domaine_Cl_dis& le_dom_Cl_dis = eqn_transp_K_Eps().domaine_Cl_dis();
      const Champ_Don ch_visco=ref_cast(Fluide_base,eqn_transp_K_Eps().milieu()).viscosite_cinematique();
      Fmu.resize(tab_K_Eps.dimension_tot(0));

      mon_modele_fonc.Calcul_Fmu( Fmu,le_dom_dis,le_dom_Cl_dis,tab_K_Eps,ch_visco);
    }

  // dans le cas d'un domaine nul on doit effectuer le dimensionnement
  double non_prepare=1;
  Debog::verifier("Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente la_viscosite_turbulente before",la_viscosite_turbulente.valeurs());
  Debog::verifier("Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente tab_K_Eps",tab_K_Eps);
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bas_Reynolds::calculer_viscosite_turbulente Fmu",Fmu);
  if (visco_turb.size() == n)
    non_prepare=0.;
  non_prepare=mp_max(non_prepare);

  if (non_prepare==1)
    {
      //Cerr << "non_prepare=1" << finl;

      Champ_Inc visco_turb_au_format_K_eps;

      visco_turb_au_format_K_eps.typer(type);
      Champ_Inc_base& ch_visco_turb_K_eps=visco_turb_au_format_K_eps.valeur();
      ch_visco_turb_K_eps.associer_domaine_dis_base(eqn_transp_K_Eps().domaine_dis().valeur());
      ch_visco_turb_K_eps.nommer("diffusivite_turbulente");
      ch_visco_turb_K_eps.fixer_nb_comp(1);
      ch_visco_turb_K_eps.fixer_nb_valeurs_nodales(n);
      ch_visco_turb_K_eps.fixer_unite("inconnue");
      ch_visco_turb_K_eps.changer_temps(0.);
      DoubleTab& visco_turb_K_eps =  ch_visco_turb_K_eps.valeurs();

      if(visco_turb_K_eps.size() != n)
        {
          Cerr << "visco_turb_K_eps size is " << visco_turb_K_eps.size()
               << " instead of " << n << finl;
          exit();
        }
      // A la fin de cette boucle, le tableau visco_turb_K_eps
      // contient les valeurs de la viscosite turbulente
      // au centre des faces du maillage.
      // Debog::verifier("Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente visco_turb_K_eps before",visco_turb_K_eps);
      for (int i=0; i<n; i++)
        {
          if (tenseur_de_Reynolds_externe_)
            {
              if (get_source_tenseur_de_Reynolds_NN( ).get_Nisizima())
                {
//                  cout << "nut computed with Nisizima" << endl;
                  visco_turb_K_eps[i] =Cmu_.valeurs()[i]*Fmu(i)*tab_K_Eps(i,0)*tab_K_Eps(i,0)/(tab_K_Eps(i,1)+1e-20);
                }
              else
                {
                  visco_turb_K_eps[i] =Cmu_.valeurs()[i]*tab_K_Eps(i,0)*tab_K_Eps(i,0)/(tab_K_Eps(i,1)+1e-20);
                }
            }

          else if (LeCmu_champ.non_nul() && is_initialized)
            {
              visco_turb_K_eps[i] =Cmu_.valeurs()[i]*tab_K_Eps(i,0)*tab_K_Eps(i,0)/(tab_K_Eps(i,1)+1e-20);
            }
          else if (is_modele_fonc)
            {
              visco_turb_K_eps[i] = Fmu(i)*LeCmu*tab_K_Eps(i,0)*tab_K_Eps(i,0)/(tab_K_Eps(i,1)+1e-20);
            }
          else
            visco_turb_K_eps[i] =LeCmu*tab_K_Eps(i,0)*tab_K_Eps(i,0)/(tab_K_Eps(i,1)+1e-20);
        }
      // Debog::verifier("Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente visco_turb_K_eps after",visco_turb_K_eps);

      // On connait donc la viscosite turbulente au centre des faces de chaque element
      // On cherche maintenant a interpoler cette viscosite turbulente au centre des
      // elements.
      la_viscosite_turbulente->affecter(visco_turb_au_format_K_eps.valeur());
      // Debog::verifier("Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente visco_turb_au_format_K_eps",visco_turb_au_format_K_eps.valeur());
    }
  else
    {
      //Cerr << "non_prepare=0" << finl;
      for (int i=0; i<n; i++)
        {
          if (tenseur_de_Reynolds_externe_)
            {
              if (get_source_tenseur_de_Reynolds_NN( ).get_Nisizima())
                {
//                  cout << "nut computed with Nisizima" << endl;
                  visco_turb[i] =Cmu_.valeurs()[i]*Fmu(i)*tab_K_Eps(i,0)*tab_K_Eps(i,0)/(tab_K_Eps(i,1)+1e-20);
                }
              else
                {
                  visco_turb[i] =Cmu_.valeurs()[i]*tab_K_Eps(i,0)*tab_K_Eps(i,0)/(tab_K_Eps(i,1)+1e-20);
                }
            }
          else if (LeCmu_champ.non_nul() && is_initialized)
            {
              visco_turb[i] =Cmu_.valeurs()[i]*tab_K_Eps(i,0)*tab_K_Eps(i,0)/(tab_K_Eps(i,1)+1e-20);
            }
          else if (is_modele_fonc)
            {
              visco_turb[i] = Fmu(i)*LeCmu*tab_K_Eps(i,0)*tab_K_Eps(i,0)/(tab_K_Eps(i,1)+1e-20);
            }
          else
            visco_turb[i] =LeCmu*tab_K_Eps(i,0)*tab_K_Eps(i,0)/(tab_K_Eps(i,1)+1e-20);
        }
    }

  la_viscosite_turbulente.changer_temps(temps);
  Debog::verifier("Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente la_viscosite_turbulente after",la_viscosite_turbulente.valeurs());
  return la_viscosite_turbulente;
}

void imprimer_evolution_keps(const Champ_Inc& le_champ_K_Eps, const Schema_Temps_base& sch, double LeCmu, int avant)
{
  if (sch.nb_pas_dt()==0 || sch.limpr())
    {
      const DoubleTab& K_Eps = le_champ_K_Eps.valeurs();
      double k_min=DMAXFLOAT;
      double eps_min=DMAXFLOAT;
      double nut_min=DMAXFLOAT;
      double k_max=0;
      double eps_max=0;
      double nut_max=0;
      int loc_k_min=-1;
      int loc_eps_min=-1;
      int loc_nut_min=-1;
      int loc_k_max=-1;
      int loc_eps_max=-1;
      int loc_nut_max=-1;
      int size=K_Eps.dimension(0);
      if (size<0)
        {
          if (sub_type(Champ_Inc_P0_base, le_champ_K_Eps.valeur()))
            size = le_champ_K_Eps.valeur().equation().domaine_dis().domaine().nb_elem();
          else
            {
              Cerr << "Unsupported K_Eps field in Modele_turbulence_hyd_K_Eps::imprimer_evolution_keps()" << finl;
              Process::exit(-1);
            }
        }
      ConstDoubleTab_parts parts(le_champ_K_Eps.valeurs());
      for (int n=0; n<size; n++)
        {
          const double k = K_Eps(n,0);
          const double eps = K_Eps(n,1);
          double nut = 0;
          if (eps > 0) nut = LeCmu*k*k/eps;
          if (k < k_min)
            {
              k_min = k;
              loc_k_min = n;
            }
          else if (k > k_max)
            {
              k_max = k;
              loc_k_max = n;
            }
          if (eps < eps_min)
            {
              eps_min = eps;
              loc_eps_min = n;
            }
          else if (eps > eps_max)
            {
              eps_max = eps;
              loc_eps_max = n;
            }
          if (nut < nut_min)
            {
              nut_min = nut;
              loc_nut_min = n;
            }
          else if (nut > nut_max)
            {
              nut_max = nut;
              loc_nut_max = n;
            }
        }
      /*
      k_min = Process::mp_min(k_min);
      eps_min = Process::mp_min(eps_min);
      nut_min = Process::mp_min(nut_min);
      k_max = Process::mp_max(k_max);
      eps_max = Process::mp_max(eps_max);
      nut_max = Process::mp_max(nut_max);
      */
      ArrOfDouble values(3);

      values[0]=k_min;
      values[1]=eps_min;
      values[2]=nut_min;
      mp_min_for_each_item(values);
      k_min=values[0];
      eps_min=values[1];
      nut_min=values[2];

      values[0]=k_max;
      values[1]=eps_max;
      values[2]=nut_max;
      mp_max_for_each_item(values);
      k_max=values[0];
      eps_max=values[1];
      nut_max=values[2];
      if (Process::je_suis_maitre())
        {
          Cout << finl << "K_Eps evolution (" << (avant?"before":"after") << " law of the wall applies) at time " << le_champ_K_Eps.temps() << ":" << finl;
          Cout << "std::min(k)=" << k_min;
          if (Process::nproc()==1) Cout << " located at node " << loc_k_min;
          Cout << finl;
          Cout << "std::min(eps)=" << eps_min;
          if (Process::nproc()==1) Cout << " located at node " << loc_eps_min;
          Cout << finl;
          Cout << "std::min(nut)=" << nut_min;
          if (Process::nproc()==1) Cout << " located at node " << loc_nut_min;
          Cout << finl;
          Cout << "std::max(k)=" << k_max;
          if (Process::nproc()==1) Cout << " located at node " << loc_k_max;
          Cout << finl;
          Cout << "std::max(eps)=" << eps_max;
          if (Process::nproc()==1) Cout << " located at node " << loc_eps_max;
          Cout << finl;
          Cout << "std::max(nut)=" << nut_max;
          if (Process::nproc()==1) Cout << " located at node " << loc_nut_max;
          Cout << finl;
        }
    }
}

int Modele_turbulence_hyd_K_Eps::preparer_calcul()
{
  eqn_transp_K_Eps().preparer_calcul();
  Mod_turb_hyd_base::preparer_calcul();
  // GF pour initialiser la loi de paroi thermique en TBLE
  if (equation().probleme().nombre_d_equations()>1)
    {
      const RefObjU& modele_turbulence = equation().probleme().equation(1).get_modele(TURBULENCE);
      if (sub_type(Modele_turbulence_scal_base,modele_turbulence.valeur()))
        {
          Turbulence_paroi_scal& loi_paroi_T = ref_cast_non_const(Modele_turbulence_scal_base,modele_turbulence.valeur()).loi_paroi();
          loi_paroi_T.init_lois_paroi();
        }
    }
  // GF quand on demarre un calcul il est bon d'utliser la ldp
  // encore plus quand on fait une reprise !!!!!!!!
  Champ_Inc& ch_K_Eps = K_Eps();

  const Milieu_base& mil=equation().probleme().milieu();
  if (equation().probleme().is_dilatable()) diviser_par_rho_si_dilatable(ch_K_Eps.valeurs(),mil);
  imprimer_evolution_keps(ch_K_Eps,eqn_transp_K_Eps().schema_temps(),LeCmu,1);
  loipar.calculer_hyd(ch_K_Eps);
  eqn_transp_K_Eps().controler_K_Eps();
  calculer_viscosite_turbulente(ch_K_Eps.temps());
  limiter_viscosite_turbulente();
  // on remultiplie K_eps par rho
  if (equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K_Eps.valeurs(),mil);
      correction_nut_et_cisaillement_paroi_si_qc(*this);
    }
  la_viscosite_turbulente.valeurs().echange_espace_virtuel();
  Debog::verifier("Modele_turbulence_hyd_K_Eps::preparer_calcul la_viscosite_turbulente",la_viscosite_turbulente.valeurs());
  imprimer_evolution_keps(ch_K_Eps,eqn_transp_K_Eps().schema_temps(),LeCmu,0);
  return 1;
}

bool Modele_turbulence_hyd_K_Eps::initTimeStep(double dt)
{
  return eqn_transport_K_Eps.initTimeStep(dt);
}

/*! @brief Effectue une mise a jour en temps du modele de turbulence.
*
* Met a jour l'equation de transport K-epsilon,
*     calcule la loi de paroi et la viscosite turbulente
*     au nouveau temps.
*
* @param (double temps) le temps de mise a jour
*/
void Modele_turbulence_hyd_K_Eps::mettre_a_jour(double temps)
{
  Champ_Inc& ch_K_Eps = K_Eps();
  Schema_Temps_base& sch =eqn_transp_K_Eps().schema_temps();
  // Voir Schema_Temps_base::faire_un_pas_de_temps_pb_base
  eqn_transp_K_Eps().domaine_Cl_dis().mettre_a_jour(temps);
  if (!eqn_transp_K_Eps().equation_non_resolue())
    sch.faire_un_pas_de_temps_eqn_base(eqn_transp_K_Eps());
  eqn_transp_K_Eps().mettre_a_jour(temps);

  statistiques().begin_count(nut_counter_);
  const Milieu_base& mil=equation().probleme().milieu();
  Debog::verifier("Modele_turbulence_hyd_K_Eps::mettre_a_jour la_viscosite_turbulente before",la_viscosite_turbulente.valeurs());
  // on divise K_eps par rho en QC pour revenir a K et Eps
  if (equation().probleme().is_dilatable()) diviser_par_rho_si_dilatable(ch_K_Eps.valeurs(),mil);
  imprimer_evolution_keps(ch_K_Eps,eqn_transp_K_Eps().schema_temps(),LeCmu,1);
  loipar.calculer_hyd(ch_K_Eps);
  eqn_transp_K_Eps().controler_K_Eps();
  calculer_viscosite_turbulente(ch_K_Eps.temps());
  limiter_viscosite_turbulente();
  // on remultiplie K_eps par rho
  if (equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K_Eps.valeurs(),mil);
      correction_nut_et_cisaillement_paroi_si_qc(*this);
    }
  la_viscosite_turbulente.valeurs().echange_espace_virtuel();
  Debog::verifier("Modele_turbulence_hyd_K_Eps::mettre_a_jour la_viscosite_turbulente after",la_viscosite_turbulente.valeurs());
  imprimer_evolution_keps(ch_K_Eps,eqn_transp_K_Eps().schema_temps(),LeCmu,0);
  statistiques().end_count(nut_counter_);
}

const Equation_base& Modele_turbulence_hyd_K_Eps::equation_k_eps(int i) const
{
  assert ((i==0));
  return eqn_transport_K_Eps;

}
const Champ_base&  Modele_turbulence_hyd_K_Eps::get_champ(const Motcle& nom) const
{

  try
    {
      return Mod_turb_hyd_RANS::get_champ(nom);
    }
  catch (Champs_compris_erreur)
    {
    }
  if (mon_modele_fonc.non_nul())
    {
      try
        {
          return  mon_modele_fonc.valeur().get_champ(nom);
        }
      catch (Champs_compris_erreur)
        {
        }
    }
  throw Champs_compris_erreur();

}
void Modele_turbulence_hyd_K_Eps::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  Mod_turb_hyd_RANS::get_noms_champs_postraitables(nom,opt);
  if (mon_modele_fonc.non_nul())
    mon_modele_fonc.valeur().get_noms_champs_postraitables(nom,opt);

}
void Modele_turbulence_hyd_K_Eps::verifie_loi_paroi()
{
  Nom lp=loipar.valeur().que_suis_je();
  if (lp=="negligeable_VEF" || lp=="negligeable_VDF")
    if (!associe_modele_fonction().non_nul())
      {
        Cerr<<"The turbulence model of type "<<que_suis_je()<<finl;
        Cerr<<"must not be used with a wall law of type negligeable or with a modele_function."<<finl;
        Cerr<<"Another wall law must be selected with this kind of turbulence model."<<finl;
      }
}

void Modele_turbulence_hyd_K_Eps::Calcul_RSLambda()
{
//  Cout << "calling Calcul_RSLambda " << endl;
  const Domaine_Cl_dis& le_dom_Cl_dis = eqn_transport_K_Eps.domaine_Cl_dis();
  const Domaine_dis& domaine_dis = le_dom_Cl_dis.valeur().domaine_dis();

  const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF,domaine_dis.valeur());
  const Domaine_Cl_VDF& dom_Cl_VDF = ref_cast(Domaine_Cl_VDF,le_dom_Cl_dis.valeur());

  int nb_elem_tot=domaine_VDF.nb_elem_tot();
  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();
  Champ_Face_VDF& ch = ref_cast(Champ_Face_VDF,mon_equation->inconnue().valeur());

  assert (vitesse.line_size() == 1);

  DoubleTab gij(nb_elem_tot,dimension,dimension, vitesse.line_size());
  ch.calcul_duidxj( vitesse,gij,dom_Cl_VDF );


  DoubleTab& lambda_1 = get_l1();
  DoubleTab& lambda_2 = get_l2();
  DoubleTab& lambda_3 = get_l3();
  DoubleTab& lambda_4 = get_l4();
  DoubleTab& lambda_5 = get_l5();

  const DoubleTab& K_eps = eqn_transport_K_Eps.inconnue().valeurs();

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
            tab_zero(elem, i, j) = 0.;
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

}
