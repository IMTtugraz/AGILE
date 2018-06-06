#ifndef IMT_MRCMD_H
#define IMT_MRCMD_H

#include <complex>

// -- std - C/C++ - includes
#include <iostream>
#include <vector>
#include <time.h>

#include "agile/calc/postprocess.hpp"
#include "agile/calc/fft.hpp"

#include "agile/calc/l2solve.hpp"
#include "agile/calc/tvsolve.hpp"
#include "agile/calc/tgvsolve.hpp"

//typedef for calculation - float or double cuda calculation
typedef std::complex<double> TType;

typedef std::complex<float> TType_data;//typedef for data - float or double data format
typedef float TType_image;  //don't change image type - no need for double-type image quality
typedef agile::to_real_type<TType_data>::type TType_real_data;
typedef agile::to_real_type<TType>::type TType_real;


class IMT_MRCMD
{

public:
    explicit IMT_MRCMD();
    ~IMT_MRCMD();

    void setIRGNParams(unsigned int maxit,
                       unsigned char tvtype,
                       unsigned int tvits,
                       unsigned int tvmax,
                       float alpha0,
                       float alpha_min,
                       float alpha_q,
                       float beta0,
                       float beta_min,
                       float beta_q);

    void setFormat(unsigned format);

    //------------
    void readMatlab(std::vector<TType_data>& matrix_read,
                    unsigned& num_rows_data,
                    unsigned& num_columns_data,
                    unsigned& num_coils_data,
                    const std::string &path);

    int startCalculation(unsigned num_rows,
                         unsigned num_columns,
                         unsigned num_coils,
                         const std::vector<TType_data>& matrix_read,
                         std::vector<TType_image>& image_data);

    void calc_postprocess(const agile::GPUMatrix<TType>* in_mat, agile::GPUMatrix<TType_real>& out_mat);

    int writeoutfile(unsigned num_rows,
                     unsigned num_columns,
                     const std::vector<TType_real_data>& image_data,
                     const std::string &path);


private:
    agile::GPUEnvironment GPU0_;



    void irgncalc(agile::IRGN<TType,TType_real>* IRGNobj, std::vector<TType_image>& image_data);

    //---
    unsigned _format;

    agile::IRGN<TType,TType_real>* IRGNobj;
    agile::IRGN_Params irgnParams;

    agile::FFT<TType>* fftobj_;

    std::vector<TType_real> zero_vec_cpu_;

    unsigned _num_rows_calc;
    unsigned _num_columns_calc;

    agile::PostProcess<TType, TType_real>* _postprocess;

    bool _readmatlab;
    std::vector<TType_data> _matlab_read;

    unsigned _matlab_num_rows;
    unsigned _matlab_num_columns;
    unsigned _matlab_num_coils;

    time_t _seconds;

};


#endif // IMT_MRCMD

