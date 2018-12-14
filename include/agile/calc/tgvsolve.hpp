
//==================================================================================
//
//                    TGVSOLVE.HPP
//
//==================================================================================

#ifndef AGILE_IRGN_TGVSOLVE_HPP
#define AGILE_IRGN_TGVSOLVE_HPP


#include "agile/calc/irgn.hpp"

namespace agile
{

  template <typename TType, typename TType2>
  class TGVSolve : public IRGN<TType, TType2>
  {
    public:

      //! \brief Default constructor.
      //!
      //! The default constructor creates an empty TGVSolve Class.
      TGVSolve()
        : IRGN<TType,TType2>()
      {
      }


      //! \brief Constructor.
      //!
      //! The Constructor creates an TGVSolve.
      TGVSolve(GPUMatrix<TType>* coil, unsigned int num_coils, IRGN_Params param)
        : IRGN<TType,TType2>(coil, num_coils, param)
      {
        Init();
      }

      //! \brief Constructor.
      //!
      //! The Constructor creates an initialized TGVSolve.
      TGVSolve(GPUMatrix<TType>* coil, unsigned int num_coils)
          : IRGN<TType,TType2>(coil, num_coils)
      {
        Init();
      }

      //! \brief Destructor.
      virtual ~TGVSolve()
      {
        delete[] cb_mat_;
        cb_mat_ = 0;
        delete[] Mc_mat_;
        Mc_mat_ = 0;
        delete[] eta4_mat_;
        eta4_mat_ = 0;
        delete[] dc_old_mat_;
        dc_old_mat_ = 0;
      }

      void Solve(const GPUMatrix<TType>* u, const GPUMatrix<TType>* c,
                                     const GPUMatrix<TType>* rhs, const GPUMatrix<TType>* u0,
                                     unsigned maxits, typename to_real_type<TType>::type alpha,
                                     typename to_real_type<TType>::type beta,
                                     GPUMatrix<TType>* du, GPUMatrix<TType>* dc);


    private:

      void Init();


      GPUMatrix<TType> ub_mat_;
      GPUMatrix<TType> ukp_mat_;
      GPUMatrix<TType> Mu_mat_;
      GPUMatrix<TType> eta3_mat_;
      GPUMatrix<TType> eta1_mat_;
      GPUMatrix<TType> eta2_mat_;
      GPUMatrix<TType> y1p_mat_;
      GPUMatrix<TType> y2p_mat_;
      GPUMatrix<TType> p1_mat_;
      GPUMatrix<TType> p2_mat_;
      GPUMatrix<TType> du_old_mat_;
      GPUMatrix<TType> safe_;
      GPUMatrix<TType> safe2_;
      GPUMatrix<TType> safe3_;

      GPUMatrix<TType> v1_mat_;
      GPUMatrix<TType> v2_mat_;
      GPUMatrix<TType> vb1_mat_;
      GPUMatrix<TType> vb2_mat_;
      GPUMatrix<TType> zeta1_mat_;
      GPUMatrix<TType> zeta2_mat_;
      GPUMatrix<TType> zeta3_mat_;
      GPUMatrix<TType> zeta4_mat_;
      GPUMatrix<TType> zeta5_mat_;
      GPUMatrix<TType> q1_mat_;
      GPUMatrix<TType> q2_mat_;
      GPUMatrix<TType> q3_mat_;
      GPUMatrix<TType> v1_old_mat_;
      GPUMatrix<TType> v2_old_mat_;


      GPUMatrix<TType>* cb_mat_;
      GPUMatrix<TType>* Mc_mat_;
      GPUMatrix<TType>* eta4_mat_;
      GPUMatrix<TType>* dc_old_mat_;

      typename agile::to_real_type<TType>::type norm_;
      typename agile::to_real_type<TType>::type one_;
      typename agile::to_real_type<TType>::type two_;

  };


}
#endif //AGILE_IRGN_TGVSOLVE_HPP
