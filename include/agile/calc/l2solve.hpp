
//==================================================================================
//
//                    L2SOLVE.HPP
//
//==================================================================================

#ifndef AGILE_IRGN_L2SOLVE_HPP
#define AGILE_IRGN_L2SOLVE_HPP


#include "agile/calc/irgn.hpp"


namespace agile
{
  using agile::GPUMatrix;
  using agile::GPUVector;


  template <typename TType, typename TType2>
  class L2Solve : public IRGN<TType, TType2>
  {
    public:

      //! \brief Default constructor.
      //!
      //! The default constructor creates an empty L2Solve Class.
      L2Solve()
        : IRGN<TType,TType2>()
      {
      }


      //! \brief Constructor.
      //!
      //! The Constructor creates an L2Solve.
      L2Solve(GPUMatrix<TType>* coil, unsigned int num_coils, IRGN_Params param)
        : IRGN<TType,TType2>(coil, num_coils, param)
      {
        Init();
      }

      //! \brief Constructor.
      //!
      //! The Constructor creates an initialized L2Solve.
      L2Solve(GPUMatrix<TType>* coil, unsigned int num_coils)
          : IRGN<TType,TType2>(coil, num_coils)
      {
        Init();
      }

      //! \brief Destructor.
      virtual ~L2Solve()
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

      GPUMatrix<TType> ub_mat_;
      GPUMatrix<TType> Mu_mat_;
      GPUMatrix<TType> eta3_mat_;
      GPUMatrix<TType> du_old_mat_;
      GPUMatrix<TType> safe_;

      GPUMatrix<TType>* cb_mat_;
      GPUMatrix<TType>* Mc_mat_;
      GPUMatrix<TType>* eta4_mat_;
      GPUMatrix<TType>* dc_old_mat_;

      typename agile::to_real_type<TType>::type norm_;
      typename agile::to_real_type<TType>::type one_;
      typename agile::to_real_type<TType>::type two_;
  };
}
#endif //AGILE_IRGN_L2SOLVE_HPP
