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
// File:        Modele_turbulence_hyd_K_Eps.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_turbulence_hyd_K_Eps_included
#define Modele_turbulence_hyd_K_Eps_included

#include <Transport_K_Eps.h>
#include <Modele_Fonc_Bas_Reynolds.h>
#include <Tenseur_Reynolds_Externe_VDF_Face.h>
#include <Champ_Don.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Modele_turbulence_hyd_K_Eps
//    Cette classe represente le modele de turbulence (k,eps) pour les
//    equations de Navier-Stokes.
// .SECTION voir aussi
//    Mod_turb_hyd_base Mod_turb_hyd_ss_maille
//////////////////////////////////////////////////////////////////////////////

class Modele_turbulence_hyd_K_Eps : public Mod_turb_hyd_RANS
{

  Declare_instanciable(Modele_turbulence_hyd_K_Eps);

public:

  void set_param(Param& param) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int preparer_calcul() override;
  void verifie_loi_paroi() override;
  bool initTimeStep(double dt) override;
  void mettre_a_jour(double ) override;
  virtual inline Champ_Inc& K_Eps();
  virtual inline const Champ_Inc& K_Eps() const;
  void discretiser() override;

  inline int nombre_d_equations() const override;
  inline Transport_K_Eps_base& eqn_transp_K_Eps() override;
  inline const Transport_K_Eps_base& eqn_transp_K_Eps() const override;
  const Equation_base& equation_k_eps(int) const override ;

  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;

  inline bool get_if_Cmu_DNS() const
  {
    return LeCmu_champ.non_nul() && !is_initialized;
  }
//  inline void set_bij(const DoubleTab& bij);
  inline DoubleTab& get_bij();
  inline DoubleTab& get_bij_NL();
  inline DoubleTab& get_l1();
  inline DoubleTab& get_l2();
  inline DoubleTab& get_l3();
  inline DoubleTab& get_l4();
  inline DoubleTab& get_l5();

  void Calcul_RSLambda();

  inline DoubleTab& lambda_1_etoile();
  inline DoubleTab& lambda_2_etoile();
  inline DoubleTab& lambda_3_etoile();
  inline DoubleTab& lambda_4_etoile();
  inline DoubleTab& lambda_5_etoile();

protected:
  Transport_K_Eps  eqn_transport_K_Eps;
  virtual Champ_Fonc& calculer_viscosite_turbulente(double temps);

  Champ_Don LeCmu_champ;
  DoubleTab LeCmu_tab;
  Champ_Fonc Cmu_, bij_, bij_NL_;
  Champ_Fonc lambda1_, lambda2_, lambda3_, lambda4_, lambda5_;
  int is_initialized = 0;

  DoubleTab lambda_1_etoile_;
  DoubleTab lambda_2_etoile_;
  DoubleTab lambda_3_etoile_;
  DoubleTab lambda_4_etoile_;
  DoubleTab lambda_5_etoile_;
};

//inline void Modele_turbulence_hyd_K_Eps::set_bij(const DoubleTab& bij)
//{
//  bij_.valeurs() = bij;
//}

inline DoubleTab& Modele_turbulence_hyd_K_Eps::get_bij()
{
  return bij_.valeurs();
}

inline DoubleTab& Modele_turbulence_hyd_K_Eps::get_bij_NL()
{
  return bij_NL_.valeurs();
}

inline DoubleTab& Modele_turbulence_hyd_K_Eps::get_l1()
{
  return lambda1_.valeurs();
}

inline DoubleTab& Modele_turbulence_hyd_K_Eps::get_l2()
{
  return lambda2_.valeurs();
}

inline DoubleTab& Modele_turbulence_hyd_K_Eps::get_l3()
{
  return lambda3_.valeurs();
}

inline DoubleTab& Modele_turbulence_hyd_K_Eps::get_l4()
{
  return lambda4_.valeurs();
}

inline DoubleTab& Modele_turbulence_hyd_K_Eps::get_l5()
{
  return lambda5_.valeurs();
}

inline DoubleTab& Modele_turbulence_hyd_K_Eps::lambda_1_etoile()
{
  return lambda_1_etoile_;
}
inline DoubleTab& Modele_turbulence_hyd_K_Eps::lambda_2_etoile()
{
  return lambda_2_etoile_;
}

inline DoubleTab& Modele_turbulence_hyd_K_Eps::lambda_3_etoile()
{
  return lambda_3_etoile_;
}

inline DoubleTab& Modele_turbulence_hyd_K_Eps::lambda_4_etoile()
{
  return lambda_4_etoile_;
}

inline DoubleTab& Modele_turbulence_hyd_K_Eps::lambda_5_etoile()
{
  return lambda_5_etoile_;
}

// Description:
//    Renvoie le champ inconnue du modele de turbulence
//    i.e. : (K,Epsilon). Cette inconnue est portee
//    par l'equation de transport K-eps porte par le modele.
//    (version const)
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Champ_Inc&
//    Signification: le champ inconnue (K,epsilon)
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
inline const Champ_Inc& Modele_turbulence_hyd_K_Eps::K_Eps() const
{
  return eqn_transport_K_Eps.inconnue();
}


// Description:
//    Renvoie le champ inconnue du modele de turbulence
//    i.e. : (K,Epsilon). Cette inconnue est portee
//    par l'equation de transport K-eps porte par le modele.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Champ_Inc&
//    Signification: le champ inconnue (K,epsilon)
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
inline Champ_Inc& Modele_turbulence_hyd_K_Eps::K_Eps()
{
  return eqn_transport_K_Eps.inconnue();
}

// Description:
//    Renvoie l'equation du modele de turbulence
//    i.e. : (K,Epsilon).
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Transport_K_Eps&
//    Signification: equation (K,epsilon)
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
inline Transport_K_Eps_base& Modele_turbulence_hyd_K_Eps::eqn_transp_K_Eps()
{
  return eqn_transport_K_Eps;
}

// Description:
//    Renvoie l'equation du modele de turbulence
//    i.e. : (K,Epsilon).
//    (version const)
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Transport_K_Eps&
//    Signification: equation (K,epsilon)
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
inline const Transport_K_Eps_base& Modele_turbulence_hyd_K_Eps::eqn_transp_K_Eps() const
{
  return eqn_transport_K_Eps;
}
inline int Modele_turbulence_hyd_K_Eps::nombre_d_equations() const
{
  return 1;
}
#endif
