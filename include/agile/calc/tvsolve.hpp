//==================================================================================
//
//                    TVSOLVE.HPP
//
//==================================================================================

#ifndef AGILE_IRGN_TVSOLVE_HPP
#define AGILE_IRGN_TVSOLVE_HPP


#include "agile/calc/irgn.hpp"

namespace agile
{
  using agile::GPUMatrix;
  using agile::GPUVector;


  template <typename TType, typename TType2>
  class TVSolve : public IRGN<TType, TType2>
  {
    public:

      //! \brief Default constructor.
      //!
      //! The default constructor creates an empty TVSolve Class.
      TVSolve()
        : IRGN<TType,TType2>()
      {
      }


      //! \brief Constructor.
      //!
      //! The Constructor creates an TVSolve.
      TVSolve(GPUMatrix<TType>* coil, unsigned int num_coils, IRGN_Params param)
        : IRGN<TType,TType2>(coil, num_coils, param)
      {
        Init();
      }

      //! \brief Constructor.
      //!
      //! The Constructor creates an initialized TVSolve.
      TVSolve(GPUMatrix<TType>* coil, unsigned int num_coils)
          : IRGN<TType,TType2>(coil, num_coils)
      {
        Init();
      }

      //! \brief Destructor.
      virtual ~TVSolve()
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
                                     unsigned maxits, typename TType::value_type alpha, typename TType::value_type beta,
                                     GPUMatrix<TType>* du, GPUMatrix<TType>* dc);


    private:

      void Init();

      GPUMatrix<TType> x1_mat_;
      GPUMatrix<TType> y1_mat_;
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

      GPUMatrix<TType>* x2_mat_;
      GPUMatrix<TType>* y2_mat_;
      GPUMatrix<TType>* cb_mat_;
      GPUMatrix<TType>* Mc_mat_;
      GPUMatrix<TType>* eta4_mat_;
      GPUMatrix<TType>* dc_old_mat_;

      typename agile::to_real_type<TType>::type norm_;
      typename agile::to_real_type<TType>::type one_;
      typename agile::to_real_type<TType>::type two_;

      std::complex<typename agile::to_real_type<TType>::type> l1_;
      std::complex<typename agile::to_real_type<TType>::type> l2_;
      std::complex<typename agile::to_real_type<TType>::type> l2_safe_;

  };

}
#endif // AGILE_IRGN_TVSOLVE_HPP
