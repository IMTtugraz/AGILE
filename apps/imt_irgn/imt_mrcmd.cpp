#include "imt_mrcmd.h"
// -- Agile - includes
#include "agile/io/file.hpp"
#include "agile/io/dicom.hpp"
//#include "agile/matrixhelper.h"
#include "agile/calc/irgn.hpp"


IMT_MRCMD::IMT_MRCMD()
{
  // Initialize the first GPU
  GPU0_.allocateGPU(0);
  GPU0_.printInformation(std::cout);

  size_t free;
  size_t total;
  cuMemGetInfo(&free, &total);
  std::cout << "free memory: " << free / 1024 / 1024 << "mb, total memory: " << total / 1024 / 1024 << "mb" << std::endl;


  fftobj_ = new agile::FFT<TType>();

  _postprocess = new agile::PostProcess<TType, TType_real>();
}

IMT_MRCMD::~IMT_MRCMD()
{
  delete fftobj_;
  fftobj_ = 0;

  delete _postprocess;
}

void IMT_MRCMD::setIRGNParams(unsigned int maxit,
                         unsigned char tvtype,
                         unsigned int tvits,
                         unsigned int tvmax,
                         float alpha0,
                         float alpha_min,
                         float alpha_q,
                         float beta0,
                         float beta_min,
                         float beta_q)
{
    irgnParams.maxit = maxit;               // maximum number of IRGN iterations
    irgnParams.tvtype = tvtype;             // regularization term: 0: L2, 1: TV, 2: TGV
    irgnParams.tvits = tvits;               // initial number of gradient steps
    irgnParams.tvmax = tvmax;               // upper bound on number of gradient steps
    irgnParams.alpha_min = alpha_min;       // final value of alpha
    irgnParams.beta_min = beta_min;         // final value of beta: 0: no T(G)V effect, >0 effect
    irgnParams.alpha0 = alpha0;             // initial penalty alpha_0 (L2, sensitivites)
    irgnParams.beta0 = beta0;               // initial penalty beta_0 (image)
    irgnParams.alpha_q = alpha_q;           // reduction factor for alpha
    irgnParams.beta_q = beta_q;             // reduction factor for beta
}

void IMT_MRCMD::setFormat(unsigned format)
{
    this->_format = format;
}

//--------------------
// readMatlab
//--------------------
void IMT_MRCMD::readMatlab(std::vector<TType_data>& matrix_read,
                           unsigned& num_rows_data,
                           unsigned& num_columns_data,
                           unsigned& num_coils_data,
                           const std::string& path)
{


  bool ok = agile::readMatrixFile3D(path.c_str(),
                                    num_rows_data,
                                    num_columns_data,
                                    num_coils_data,
                                    matrix_read);

  if(!ok)
  {
    std::cerr<<"\n Error: reading Matlab-file - error-number: "<<ok;
  }
}


//--------------------
// start calculation
//--------------------
int IMT_MRCMD::startCalculation(unsigned num_rows,
                                unsigned num_columns,
                                unsigned num_coils,
                                const std::vector<TType_data>& matrix_read,
                                std::vector<TType_image>& image_data)
{
  std::vector<TType_real> image_cpu_vec;
  std::vector<TType> matrix_calc(matrix_read.begin(), matrix_read.end());

  if(matrix_calc.empty())
    return -1;               // return 1 - if matrix is empty

  agile::GPUMatrix<TType>* Coil;
  Coil = new agile::GPUMatrix<TType> [num_coils];
  agile::GPUMatrix<TType>* Coil_erg;
  Coil_erg = new agile::GPUMatrix<TType> [num_coils];

  _postprocess->set_size(num_rows,num_columns,num_coils);

  for(unsigned i=0; i<num_coils;++i)
  {
    Coil[i].assignFromHost(num_rows, num_columns, &matrix_calc[num_rows*num_columns*i]);
    Coil_erg[i].resize(Coil[i].getNumRows(),Coil[i].getNumColumns());
  }

  for(int i=0; i<num_columns*num_rows;  ++i)
    zero_vec_cpu_.push_back(TType_real(0));

  //--------------------------- Start IRGN, IFFT, RAW Calculation  -------------------------
  if(irgnParams.tvtype == 3)
  {
    IRGNobj = new agile::TGVSolve<TType,TType_real>(Coil,num_coils, irgnParams);
    irgncalc(IRGNobj, image_data);

    delete IRGNobj;
  }
  else if(irgnParams.tvtype == 2)
  {
    IRGNobj = new agile::TVSolve<TType,TType_real>(Coil,num_coils, irgnParams);
    irgncalc(IRGNobj, image_data);

    delete IRGNobj;
  }
  else if(irgnParams.tvtype == 1)
  {
    IRGNobj = new agile::L2Solve<TType,TType_real>(Coil,num_coils, irgnParams);
    irgncalc(IRGNobj, image_data);

    delete IRGNobj;
  }
  else if(irgnParams.tvtype == 0)
  {
    agile::GPUMatrix<TType_real> erg_val(num_rows,num_columns,&zero_vec_cpu_[0]);

    int error = 0;

    fftobj_->setfftplan(num_rows,num_columns);

    for(unsigned i=0;i<num_coils;++i)
    {
      fftobj_->calc_pattern(Coil[0]);
      error = fftobj_->CenterdIFFT(Coil[i],Coil_erg[i]);       //ifft each Coil Matrix

      if(error != 0)
      {
        std::cerr<<"\n !err: "<<error;
        return error;
      }
    }

    calc_postprocess(&Coil_erg[0],erg_val);
    erg_val.copyToHost(image_cpu_vec);                  //copy solution back to host
    std::vector<TType_image> img(image_cpu_vec.begin(), image_cpu_vec.end());

    image_data.assign(img.begin(), img.end());
  }

  delete[] Coil;
  delete[] Coil_erg;

  Coil=0;
  Coil_erg=0;

  return 0;
}


//--------------------
// irgn calculations
//--------------------
void IMT_MRCMD::irgncalc(agile::IRGN<TType,TType_real>* IRGNobj, std::vector<TType_image>& image_data)
{
  agile::GPUMatrix<TType>* us_matrix;
  agile::GPUMatrix<TType_real> erg_val;
  std::vector<TType_real> image_cpu_vec;

  IRGNobj->HighFreqPenalty();
  IRGNobj->Normalize();
  IRGNobj->Iteration();
  IRGNobj->Postprocess();

  us_matrix = IRGNobj->get_us_mat();

  _postprocess->set_size(IRGNobj->get_coil()->getNumRows(),IRGNobj->get_coil()->getNumColumns());
  calc_postprocess(&us_matrix[0],erg_val);
  erg_val.copyToHost(image_cpu_vec);
  image_data.assign(image_cpu_vec.begin(), image_cpu_vec.end());
}

//--------------------
// calc postprocess
//--------------------
void IMT_MRCMD::calc_postprocess(const agile::GPUMatrix<TType>* in_mat, agile::GPUMatrix<TType_real>& out_mat)
{
  switch(_format)
  {
  //real:
   case 0: _postprocess->calc_real(in_mat,out_mat);
           break;
  //imag:
   case 1: _postprocess->calc_imag(in_mat,out_mat);
           break;
  //abs:
   case 2: _postprocess->calc_abs(in_mat,out_mat);
           break;
  //phase:
   case 3: _postprocess->calc_phase(in_mat,out_mat);
           break;
  }
}


//--------------------
// write out file
//--------------------
int IMT_MRCMD::writeoutfile(unsigned num_rows,
                            unsigned num_columns,
                            const std::vector<TType_real_data>& image_data,
                            const std::string &path)
{
  agile::DICOM dicomfile;
  dicomfile.gendicom(path.c_str(), image_data, num_rows, num_columns);   //generate Dicom file
}

