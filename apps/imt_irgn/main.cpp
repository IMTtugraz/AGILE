// -- QtGui - includes
#include <QtGui/QApplication>
#include "gui/imt_mrgui.h"
#include "imt_mrcmd.h"

//#include <getopt.h>
#include <stdio.h>     /* for printf */
#include <stdlib.h>    /* for exit */
#include <string.h>
#include "optionparser.h"



struct Arg: public option::Arg
{
  static void printError(const char* msg1, const option::Option& opt, const char* msg2)
  {
    fprintf(stderr, "%s", msg1);
    fwrite(opt.name, opt.namelen, 1, stderr);
    fprintf(stderr, "%s", msg2);
  }

  static option::ArgStatus Unknown(const option::Option& option, bool msg)
  {
    if (msg) printError("Unknown option '", option, "'\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Required(const option::Option& option, bool msg)
  {
    if (option.arg != 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires an argument\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
  {
    if (option.arg != 0 && option.arg[0] != 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires a non-empty argument\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Numeric(const option::Option& option, bool msg)
  {
    char* endptr = 0;
    if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
    if (endptr != option.arg && *endptr == 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires a numeric argument\n");
    return option::ARG_ILLEGAL;
  }
};

enum  optionIndex { UNKNOWN,
                    GUI,
                    HELP,
                    INPUT,
                    OUTPUT,
                    FORMAT,
                    CALC,
                    A0,
                    AMIN,
                    AQ,
                    B0,
                    BMIN,
                    BQ,
                    TVITS,
                    TVMAX,
                    MAXIT};

const option::Descriptor usage[] = {
{ UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: example_arg [options]\n\n"
                                          "Options:" },
{ HELP,    0,"", "help",    Arg::None,    "  \t--help  \tPrint usage and exit." },
{ GUI,    0,"", "GUI",    Arg::None,    "  \t--GUI  \tstart the GUI" },
{ INPUT,  0,"i","input", Arg::Required,"  -i[<arg>], \t--input[=<arg>]"
                                          "  \tinput matlab.bin filename" },
{ OUTPUT,  0,"o","output",Arg::Required,"  -o[<arg>], \t--output[=<arg>]"
                                              "  \touput filename" },
{ FORMAT    ,0,"f","format",Arg::Required,"  -f[<arg>], \t--format[=<arg>]"
                                                "  \toutput format: '0..real', '1..imag', '2..abs', '3..phase'" },
{ CALC  ,0,"c","calc",Arg::Required,"  -c[<arg>], \t--calc[=<arg>]"
                                    "  \tcalculation type: '0..inv', '1..l2', '2..tv', '3..tgv'" },
{ A0, 0,"a0","alpha0", Arg::Numeric, "  -a0 <num>, \t--alpha0=<num>  \talpha0 value" },
{ AMIN, 0,"amin","alphamin", Arg::Numeric, "  -amin <num>, \t--alphamin=<num>  \talphamin value" },
{ AQ, 0,"aq","alphaq", Arg::Numeric, "  -aq <num>, \t--alphaq=<num>  \talphaq value" },
{ B0, 0,"b0","beta0", Arg::Numeric, "  -b0 <num>, \t--beta0=<num>  \tbeta0 value" },
{ BMIN, 0,"bmin","betamin", Arg::Numeric, "  -bmin <num>, \t--betamin=<num>  \tbetamin value" },
{ BQ, 0,"bq","betaq", Arg::Numeric, "  -bq <num>, \t--betaq=<num>  \tbetaq value" },
{ TVITS, 0,"tvi","tvits", Arg::Numeric, "  -ti <num>, \t--tvits=<num>  \tstart tv iterations value" },
{ TVMAX, 0,"tvm","tvmax", Arg::Numeric, "  -tvm <num>, \t--tvmax=<num>  \tmax tv iterations value" },
{ MAXIT, 0,"mi","maxit", Arg::Numeric, "  -mi <num>, \t--maxit=<num>  \tmax iterations value" },
{ UNKNOWN, 0,"", "",        Arg::None,
 "\nExamples:\n"
 "  ./IMT_IRGN -fmeas.dat \n"
 "  ./IMT_IRGN --format=meas.dat \n"
 "  ./IMT_IRGN -ibraindata.bin -orecon.dcm -fabs -cinv \n"
 "  ./IMT_IRGN --input=braindata.bin --output=recon.dcm --format=2 --calc=0 \n"

},
{ 0, 0, 0, 0, 0, 0 } };

int main(int argc, char *argv[])
{
    //---
    //Initial IRGN values
    unsigned int maxit = 6;
    unsigned char tvtype = 0;       //inv
    unsigned int tvits = 20;
    unsigned int tvmax = 1000;
    float alpha0 = 1;
    float alpha_min = 0;
    float alpha_q = 0.1;
    float beta0 = 1;
    float beta_min = 0;
    float beta_q = 0.2;

    //declarations
    std::string matlabPath;
    std::string outPath;

    std::vector<TType_data> matrix_read;
    std::vector<TType_image> image_data;
    unsigned num_rows_data;
    unsigned num_columns_data;
    unsigned num_coils_data;

    matrix_read.reserve(2*num_columns_data*num_rows_data*num_coils_data);
    image_data.reserve(num_columns_data*num_rows_data);
    //------------------------------
    unsigned format = 2;        //abs
    bool startGUI = false;


    argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
    option::Stats stats(usage, argc, argv);

    #ifdef __GNUC__
        // GCC supports C99 VLAs for C++ with proper constructor calls.
      option::Option options[stats.options_max], buffer[stats.buffer_max];
    #else
        // use calloc() to allocate 0-initialized memory. It's not the same
        // as properly constructed elements, but good enough. Obviously in an
        // ordinary C++ program you'd use new[], but this file demonstrates that
        // TLMC++OP can be used without any dependency on the C++ standard library.
      option::Option* options = (option::Option*)calloc(stats.options_max, sizeof(option::Option));
      option::Option* buffer  = (option::Option*)calloc(stats.buffer_max,  sizeof(option::Option));
    #endif

    option::Parser parse(usage, argc, argv, options, buffer);

    if (parse.error())
        return 1;

    if (options[HELP] || argc == 0)
    {
        int columns = getenv("COLUMNS")? atoi(getenv("COLUMNS")) : 80;
        option::printUsage(fwrite, stdout, usage, columns);

        startGUI = true;
    }

    for (int i = 0; i < parse.optionsCount(); ++i)
    {
        char* endptr = 0;

        option::Option& opt = buffer[i];
//        fprintf(stdout, "Argument #%d is ", i);
        switch (opt.index())
        {
            case GUI:
                startGUI = true;
            break;
            case INPUT:
                matlabPath = opt.arg;
            break;
            case OUTPUT:
                outPath = opt.arg;
            break;
            case FORMAT:
                format = (unsigned int)strtoul(opt.arg, &endptr, 10);
            break;
            case CALC:
                tvtype = (unsigned char)strtoul(opt.arg, &endptr, 10);
            break;
            case A0:
                alpha0 = (float)strtod(opt.arg, &endptr);
            break;
            case AMIN:
                alpha_min = (float)strtod(opt.arg, &endptr);
            break;
            case AQ:
                alpha_q = (float)strtod(opt.arg, &endptr);
            break;
            case B0:
                beta0 = (float)strtod(opt.arg, &endptr);
            break;
            case BMIN:
                beta_min = (float)strtod(opt.arg, &endptr);
            break;
            case BQ:
                beta_q = (float)strtod(opt.arg, &endptr);
            break;
            case TVMAX:
                tvmax = (unsigned int)strtoul(opt.arg, &endptr, 10);
            break;
            case TVITS:
                tvits = (unsigned int)strtoul(opt.arg, &endptr, 10);
            break;
            case MAXIT:
                maxit = (unsigned int)strtoul(opt.arg, &endptr, 10);
            break;
            case HELP:
            // not possible, because handled further above and exits the program
            case UNKNOWN:
                std::cout<<std::endl<<"UNKNOWN----";
            // not possible because Arg::Unknown returns ARG_ILLEGAL
            // which aborts the parse with an error
            break;
        }
    }

    for (int i = 0; i < parse.nonOptionsCount(); ++i)
    {
        fprintf(stdout, "Non-option argument #%d is %s\n", i, parse.nonOption(i));
        int columns = getenv("COLUMNS")? atoi(getenv("COLUMNS")) : 80;
        option::printUsage(fwrite, stdout, usage, columns);
        return 0;
    }

    if(startGUI)
    {
        Q_INIT_RESOURCE(images);

        // start GUI Application
        QApplication a(argc, argv);
        IMT_MRGui w;
        //w.resize(400,640);
        w.setFixedSize(705,570);
        w.show();

        QIcon icon(":/images/icon24x24.png");
        w.setWindowIcon(icon);

        a.exec();
    }
    else
    {

        std::cout<<std::endl<<"Configuration:";
        std::cout<<std::endl<<"input: "<<matlabPath;
        std::cout<<std::endl<<"output: "<<outPath;
        std::cout<<std::endl<<"format: "<<format;

        std::cout<<std::endl<<"-";
        std::cout<<std::endl<<"calc: "<<(unsigned int)tvtype;
        std::cout<<std::endl<<"maxit: "<<maxit;
        std::cout<<std::endl<<"tvits: "<<tvits;
        std::cout<<std::endl<<"tvmax: "<<tvmax;
        std::cout<<std::endl<<"alpha0: "<<alpha0;
        std::cout<<std::endl<<"alphaq: "<<alpha_q;
        std::cout<<std::endl<<"alphamin: "<<alpha_min;
        std::cout<<std::endl<<"beta0: "<<beta0;
        std::cout<<std::endl<<"betaq: "<<beta_q;
        std::cout<<std::endl<<"betamin: "<<beta_min<<std::endl<<std::endl;


        //------------- set imt_mrcmd -------------
        IMT_MRCMD imt_mrcmd;

        imt_mrcmd.setIRGNParams(maxit,
                                tvtype,
                                tvits,
                                tvmax,
                                alpha0,
                                alpha_min,
                                alpha_q,
                                beta0,
                                beta_min,
                                beta_q);

        imt_mrcmd.setFormat(format);

        //---
        imt_mrcmd.readMatlab(matrix_read,
                             num_rows_data,
                             num_columns_data,
                             num_coils_data,
                             matlabPath);

        std::cout<<std::endl<<"Rawdata info:";
        std::cout<<"num_rows_data: "<<num_rows_data;
        std::cout<<std::endl;
        std::cout<<"num_columns_data: "<<num_columns_data;
        std::cout<<std::endl;
        std::cout<<"num_coils_data: "<<num_coils_data;
        std::cout<<std::endl;

        imt_mrcmd.startCalculation(num_rows_data,
                                   num_columns_data,
                                   num_coils_data,
                                   matrix_read,
                                   image_data);

        imt_mrcmd.writeoutfile(num_rows_data,
                               num_columns_data,
                               image_data,
                               outPath);
    }

    exit(EXIT_SUCCESS);

  // ===============================================================================
/*
  unsigned int num_rows = 0;
  unsigned int num_columns = 0;
  unsigned int num_coils = 0;
  std::vector<TType> matrix_read;

  Uint8* pData;
  unsigned long length = 0;

  std::ofstream dcmmeas_file("dcm-meas.dat", std::ofstream::binary);
  if (!dcmmeas_file.is_open()) {
    std::cerr << "not able to write to file: " << "dcm-meas.dat" << std::endl;
    return false;
  }

  readdicom("001_000004_000001.dcm",pData, length, dcmmeas_file);
  readdicom("001_000004_000002.dcm",pData, length, dcmmeas_file);
  readdicom("001_000005_000001.dcm",pData, length, dcmmeas_file);
  readdicom("001_000005_000002.dcm",pData, length, dcmmeas_file);
  readdicom("001_000005_000003.dcm",pData, length, dcmmeas_file);
  readdicom("001_000005_000004.dcm",pData, length, dcmmeas_file);


  dcmmeas_file.close();


  // ========= ReadSiemens =========

  unsigned short int acq, sli, par, echo, pha, rep, set, seg;

  ReadSiemensVD11* read_siemens;
  //read_siemens = new ReadSiemens("meas_MID253_t2_tse_384_singleSlice_triple_FID36961.dat");
  //read_siemens = new ReadSiemensVD11("meas_MID131_localizer_FID7889.dat");

  //read_siemens = new ReadSiemensVD11("meas_MID93_t1_fl2d_tra_FID11734.dat");
  //read_siemens = new ReadSiemensVD11("meas_MID94_t2_tse_tra_FID11735.dat");

  //read_siemens = new ReadSiemensVD11("meas_MID95_t1_fl2d_tra_2averages_FID11736.dat");
  //read_siemens = new ReadSiemensVD11("meas_MID131_localizer_FID7889.dat");
  read_siemens = new ReadSiemensVD11("dcm-meas.dat");

  int error = read_siemens->readfile(true);   //true for dicom-read
  if (error == 0)
  {
    read_siemens->getRawdata_info(acq, sli, par, echo, pha, rep, set, seg);

    std::cout<<"\n acq: "<<acq<<"  sli: "<<sli<<"  par: "<<par<<"  echo: "<<echo<<"  pha: "<<pha<<"  rep: "<<rep<<"  set: "<<set<<"  seg: "<<seg;

    read_siemens->getRawdata(num_rows, num_columns, num_coils, matrix_read,0,1,0,0,0,0,0,0);

    std::cout<<"\n num: "<<num_rows<<"  "<<num_columns<<"  "<<num_coils;

    agile::writeMatrixFile3D("measrawimage.dat",num_rows,num_columns,num_coils,matrix_read);


  }
  //output("\n rawdata: ",10,10,matrix_read);

  //agile::readMatrixFile3D("braindata.bin", num_rows, num_columns, num_coils, matrix_read);

  delete read_siemens;
  read_siemens = 0;


*/
  return 0;
}

