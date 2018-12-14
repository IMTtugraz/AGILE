#include "agile/calc/genkspacefov.hpp"


namespace agile
{

  // ---------------------------------------------------------------
  //! \brief genkspace_fov(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat)
  //!  crops the fov from kspace
  // ---------------------------------------------------------------
  template <typename TType>
  int KSpaceFOV<TType>::genkspace_fov(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat)
  {
    unsigned num_rows = in_mat.getNumRows();
    unsigned num_cols = in_mat.getNumColumns();

    out_mat.resize(num_rows, num_rows);                     //size of num_rows x num_rows
    // cut the picture
    if(num_rows > num_cols)
    {
      std::cerr<<"\n num_rows > num_cols: "<<num_rows<<" > "<<num_cols;
      return -1;
    }
    else if(num_rows == num_cols)
    {
      agile::copy(in_mat, out_mat);
      //std::cerr<<"\n num_rows == num_cols: "<<num_rows<<" == "<<num_cols;
      return 0;
    }

    GPUMatrix<TType> pic_mat;
    GPUMatrix<TType> patt_mat;
    GPUMatrix<TType> patt_mat_big;


    pic_mat.resize(num_rows, num_cols);
    patt_mat.resize(num_rows, num_rows);
    patt_mat_big.resize(num_rows, num_cols);

    fftobj_->setfftplan(num_rows,num_cols);
    fftobj_->calc_pattern(in_mat);
    agile::copy(*fftobj_->get_pattern(),patt_mat_big);

    // Generate Picture-Mat
    fftobj_->CenterdIFFT(in_mat,pic_mat);

    int coloffset_left = (num_cols-num_rows)/2;

    agile::get_content(pic_mat,0,coloffset_left,out_mat);
    agile::get_content(patt_mat_big,0,coloffset_left,patt_mat);
    fftobj_->set_pattern(patt_mat);

    //generate kspace from picture
    fftobj_->setfftplan(num_rows,num_rows);
    fftobj_->CenterdFFTpattern(out_mat,out_mat);
  }
}
