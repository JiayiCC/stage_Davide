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
// File:        Tenseur_Reynolds_Externe_VDF_Face.h
// Directory:   $TURBULENCE_ROOT/src/Kernel/IA
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Tenseur_Reynolds_Externe_VDF_Face_included
#define Tenseur_Reynolds_Externe_VDF_Face_included

#include <Source_base.h>
#include <Terme_Source_Qdm.h>
#include <Equation.h>
#include <TBNN.h>
#include <TRUST_Ref.h>
#include <TRUSTTab.h>

class Modele_turbulence_hyd_K_Eps;
class Navier_Stokes_Turbulent;
class Transport_K_Eps;
class Domaine_Cl_VDF;
class Domaine_Cl_dis;
class Domaine_VDF;
class Probleme_base;

/*! @brief class Tenseur_Reynolds_Externe_VDF_Face
 *
 *
 *
 * @sa Source_base
 */
class Tenseur_Reynolds_Externe_VDF_Face : public Source_base, public Terme_Source_Qdm
{

  Declare_instanciable_sans_constructeur_ni_destructeur(Tenseur_Reynolds_Externe_VDF_Face);

public:
  Tenseur_Reynolds_Externe_VDF_Face();
  ~Tenseur_Reynolds_Externe_VDF_Face() override;

  void associer_pb(const Probleme_base& ) override;
  inline int has_interface_blocs() const override { return 1; }
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const override;
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const override {}
  DoubleTab& calculer(DoubleTab& ) const override;
  void mettre_a_jour(double ) override;
  void completer() override;
  void Calcul_RSLambdaT();
  const DoubleTab& get_bij() const;
  const DoubleTab& get_bij_NL() const;
  const DoubleTab& get_g1() const ;
  inline bool get_Nisizima() const ;
  inline void only_g1( bool val );
  inline void apply_keps_NN( bool val );
  inline void apply_keps_DNS( bool val );
  inline void apply_Nisizima( bool val );
  void set_param(Param& param);

protected:
  void readNN();
  void readNN_keps();

  REF(Navier_Stokes_Turbulent)           eqn_NS_;
  REF(Modele_turbulence_hyd_K_Eps)       modele_K_Eps_;
  REF(Probleme_base)                     probleme_;
  REF(Transport_K_Eps)                   eqn_transport_K_Eps_;

  REF(Domaine_VDF) le_dom_VDF;
  REF(Domaine_Cl_VDF) le_dom_Cl_VDF;
  Champ_Don KEps_champ;

  void associer_domaines(const Domaine_dis& ,const Domaine_Cl_dis& ) override;
  DoubleTab& Calcul_bij_TBNN(DoubleTab& resu);
  void Calcul_bij_TBNN_carre();
  DoubleTab& Calcul_bij_NL_TBNN(DoubleTab& );
  DoubleTab& Calcul_bij_NL_TBNN_carre(DoubleTab& );
//  DoubleTab& Calcul_Tenseur_Reynolds( DoubleTab& );
  DoubleTab& Calcul_Tenseur_Reynolds_NL( DoubleTab& );
  DoubleTab& Calcul_bij_NL_Nisizima(DoubleTab& resu);
  void Calcul_bij_Nisizima();
  DoubleTab g0_;
  DoubleTab g1_;
  DoubleTab g2_;
  DoubleTab g3_;
  DoubleTab g4_;
  DoubleTab g5_;
  DoubleTab g6_;
  DoubleTab g7_;
  DoubleTab g8_;
  DoubleTab g9_;
  DoubleTab g10_;
  DoubleTab bij_;
  DoubleTab bij_NL_;
  DoubleTab KEps_true_;

  DoubleTab lambda_1_etoile_;
  DoubleTab lambda_2_etoile_;
  DoubleTab lambda_3_etoile_;
  DoubleTab lambda_4_etoile_;
  DoubleTab lambda_5_etoile_;
  DoubleTab T0_etoile_;
  DoubleTab T1_etoile_;
  DoubleTab T2_etoile_;
  DoubleTab T3_etoile_;
  DoubleTab T4_etoile_;
  DoubleTab T5_etoile_;
  DoubleTab T6_etoile_;
  DoubleTab T7_etoile_;
  DoubleTab T8_etoile_;
  DoubleTab T9_etoile_;
  DoubleTab T10_etoile_;
  std::vector<int> T_list_;

  int nelem_;

  int idt_mettre_a_jour=2147483647;
  int idt_1=0;
  int idt_2=0;
  double t_stab=-4;
  double T_stab=-1;
  double relax_g1=1;
  double relax_RST=1;
  DoubleTab g1_avant_;
  DoubleTab s_avant_;
  Nom nn_casename, nn_casename_k, nn_casename_eps;                   // nom du reseau de neurones a charger
  TBNN *tbnn, *tbnn_k, *tbnn_eps;                        // objet reseau de neurones
  std::vector<double> read_file(const std::string& filename);
  bool is_only_g1_=false;
  bool apply_keps_DNS_=false;
  bool apply_keps_NN_=false;
  bool keps_NN_applied=false;
  bool keps_DNS_applied=false;
  bool apply_Nisizima_=false;
};

inline void Tenseur_Reynolds_Externe_VDF_Face::only_g1( bool val )
{
  is_only_g1_ = val;
}

inline void Tenseur_Reynolds_Externe_VDF_Face::apply_keps_NN( bool val )
{
  apply_keps_NN_ = val;
}

inline void Tenseur_Reynolds_Externe_VDF_Face::apply_keps_DNS( bool val )
{
  apply_keps_DNS_ = val;
}

inline void Tenseur_Reynolds_Externe_VDF_Face::apply_Nisizima( bool val )
{
  apply_Nisizima_ = val;
}

inline bool Tenseur_Reynolds_Externe_VDF_Face::get_Nisizima() const
{
  return apply_Nisizima_;
}

#endif
