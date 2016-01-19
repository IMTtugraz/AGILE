// Copyright (C) 2010-2011 Institute of Medical Engineering,
// Graz University of Technology
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses>.



#include "agile/gpu_matrix.hpp"
#include "agile/gpu_complex.hpp"
#include "agile/gpu_vector.hpp"

#include <fstream>
#include <limits>
#include <iostream>
#include <vector>
#include <list>
#include <math.h>
#include <stdexcept>



#ifndef VD11_IMAGE
#define VD11_IMAGE

class VD11_Image
{
  friend std::ostream &operator<<(std::ostream &, const VD11_Image &);

  public:

    struct sChanneldata
    {
      std::vector<std::complex<float> > rawdata;    //rawdata to corresponding properties
    };


    VD11_Image();
    //VD11_Image(const VD11_Image &);
    ~VD11_Image(){};
    //VD11_Image &operator=(const VD11_Image &rhs);
    int operator==(const VD11_Image &rhs) const;
    int operator<(const VD11_Image &rhs) const;
    //% properties:
    unsigned short int   ushSamplesInScan;                     //# of samples acquired in scan
    unsigned short int   ushUsedChannels;                      //# of channels used in scan
    unsigned short int   ushKSpaceCentreColumn;                //centre of echo
    float                fReadOutOffcentre;                    //ReadOut offcenter value
    unsigned             ulTimeSinceLastRF;                    //Sequence time stamp since last RF pulse
    unsigned short int   ushKSpaceCentreLineNo;                //number of K-space centre line
    unsigned short int   ushKSpaceCentrePartitionNo;           //number of K-space centre partition    };

    // -----------------------------
    // SET - METHODS:
    // -----------------------------
    void set_channeldata(unsigned short int channelid, const std::vector<std::complex<float> >& rawdata)
    {
     sChanneldata chdata;
      std::vector<sChanneldata>::iterator it;

      chdata.rawdata = rawdata;

      it = this->channeldata.begin();

      this->channeldata.insert(it+channelid,chdata);
    }

    void set_ushLine(unsigned short int   ushLine)
    {
      this->ushLine = ushLine;
      ushLine_length = (ushLine+1 >= ushLine_length) ? ushLine+1 : ushLine_length;
    }

    void set_ushAcquisition(unsigned short int   ushAcquisition)
    {
      this->ushAcquisition = ushAcquisition;
      ushAcquisition_length = (ushAcquisition >= ushAcquisition_length) ? ushAcquisition : ushAcquisition_length;
    }

    void set_ushSlice(unsigned short int   ushSlice)
    {
      this->ushSlice = ushSlice;
      ushSlice_length = (ushSlice >= ushSlice_length) ? ushSlice : ushSlice_length;
    }

    void set_ushPartition(unsigned short int   ushPartition)
    {
      this->ushPartition = ushPartition;
      ushPartition_length = (ushPartition >= ushPartition_length) ? ushPartition : ushPartition_length;
    }

    void set_ushEcho(unsigned short int   ushEcho)
    {
      this->ushEcho = ushEcho;
      ushEcho_length = (ushEcho >= ushEcho_length) ? ushEcho : ushEcho_length;
    }

    void set_ushPhase(unsigned short int   ushPhase)
    {
      this->ushPhase = ushPhase;
      ushPhase_length = (ushPhase >= ushPhase_length) ? ushPhase : ushPhase_length;
    }

    void set_ushRepetition(unsigned short int   ushRepetition)
    {
      this->ushRepetition = ushRepetition;
      ushRepetition_length = (ushRepetition >= ushRepetition_length) ? ushRepetition : ushRepetition_length;
    }

    void set_ushSet(unsigned short int   ushSet)
    {
      this->ushSet = ushSet;
      ushSet_length = (ushSet >= ushSet_length) ? ushSet : ushSet_length;
    }

    void set_ushSeg(unsigned short int   ushSeg)
    {
      this->ushSeg = ushSeg;
      ushSeg_length = (ushSeg >= ushSeg_length) ? ushSeg : ushSeg_length;
    }

  // -----------------------------
  // GET - METHODS:
  // -----------------------------
    std::vector<std::complex<float> > get_channeldata(unsigned short int channelid)
    {

      if(this->channeldata.size() <= channelid)
        std::cerr<<" wrong channelid ";

      return this->channeldata[channelid].rawdata;
    }

    unsigned short int get_ushLine_length(void)
    {
      return this->ushLine_length;
    }

    unsigned short int get_ushAcquisition_length(void)
    {
      return this->ushAcquisition_length;
    }

    unsigned short int get_ushSlice_length(void)
    {
      return this->ushSlice_length;
    }

    unsigned short int get_ushPartition_length(void)
    {
      return this->ushPartition_length;
    }

    unsigned short int get_ushEcho_length(void)
    {
      return this->ushEcho_length;
    }

    unsigned short int get_ushPhase_length(void)
    {
      return this->ushPhase_length;
    }

    unsigned short int get_ushRepetition_length(void)
    {
      return this->ushRepetition_length;
    }

    unsigned short int get_ushSet_length(void)
    {
      return this->ushSet_length;
    }

    unsigned short int get_ushSeg_length(void)
    {
      return this->ushSeg_length;
    }

  private:
    std::vector<sChanneldata > channeldata; //data % holding the raw data

    unsigned short int   ushLine;                             //line index
    unsigned short int   ushAcquisition;                      //acquisition index
    unsigned short int   ushSlice;                            //slice index
    unsigned short int   ushPartition;                        //partition index
    unsigned short int   ushEcho;                             //echo index
    unsigned short int   ushPhase;                            //phase index
    unsigned short int   ushRepetition;                       //measurement repeat index
    unsigned short int   ushSet;                              //set index
    unsigned short int   ushSeg;                              //segment index  (for TSE)

    static unsigned short int   ushLine_length;                 //line length
    static unsigned short int   ushAcquisition_length;          //acquisition length
    static unsigned short int   ushSlice_length;                //slice length
    static unsigned short int   ushPartition_length;            //partition length
    static unsigned short int   ushEcho_length;                 //echo length
    static unsigned short int   ushPhase_length;                //phase length
    static unsigned short int   ushRepetition_length;           //measurement repeat length
    static unsigned short int   ushSet_length;                  //set length
    static unsigned short int   ushSeg_length;                  //segment length  (for TSE)

};

unsigned short int VD11_Image::ushLine_length = 0;
unsigned short int VD11_Image::ushAcquisition_length = 0;
unsigned short int VD11_Image::ushSlice_length = 0;
unsigned short int VD11_Image::ushPartition_length = 0;
unsigned short int VD11_Image::ushEcho_length = 0;
unsigned short int VD11_Image::ushPhase_length = 0;
unsigned short int VD11_Image::ushRepetition_length = 0;
unsigned short int VD11_Image::ushSet_length = 0;
unsigned short int VD11_Image::ushSeg_length = 0;

VD11_Image::VD11_Image()   // Constructor
{
  this->ushSamplesInScan = 0;
  this->ushUsedChannels = 0;
  this->ushKSpaceCentreColumn = 0;
  this->fReadOutOffcentre = 0;
  this->ulTimeSinceLastRF = 0;
  this->ushKSpaceCentreLineNo = 0;
  this->ushKSpaceCentrePartitionNo = 0;

  this->ushLine         = 0;
  this->ushAcquisition  = 0;
  this->ushSlice        = 0;
  this->ushPartition    = 0;
  this->ushEcho         = 0;
  this->ushPhase        = 0;
  this->ushRepetition   = 0;
  this->ushSet          = 0;
  this->ushSeg          = 0;
}
/*
VD11_Image::VD11_Image(const VD11_Image &copyin)   // Copy constructor to handle pass by value.
{
  this->ushSamplesInScan      = copyin.ushSamplesInScan;
  this->ushUsedChannels       = copyin.ushUsedChannels;
  this->ushKSpaceCentreColumn = copyin.ushKSpaceCentreColumn;
  this->fReadOutOffcentre     = copyin.fReadOutOffcentre;
  this->ulTimeSinceLastRF     = copyin.ulTimeSinceLastRF;
  this->ushKSpaceCentreLineNo = copyin.ushKSpaceCentreLineNo;
  this->ushKSpaceCentrePartitionNo = copyin.ushKSpaceCentrePartitionNo;

  set_ushLine(copyin.ushLine);
  set_ushAcquisition(copyin.ushAcquisition);
  set_ushSlice(copyin.ushSlice);
  set_ushPartition(copyin.ushPartition);
  set_ushEcho(copyin.ushEcho);
  set_ushPhase(copyin.ushPhase);
  set_ushRepetition(copyin.ushRepetition);
  set_ushSet(copyin.ushSet);
  set_ushSeg(copyin.ushSeg);
}
*/

std::ostream &operator<<(std::ostream &output, const VD11_Image &image)
{
  output << image.ushLine << ' ' << image.ushAcquisition << ' ' << image.ushSlice << ' '
         << image.ushPartition << ' ' << image.ushEcho << ' ' << image.ushPhase << ' ' << image.ushRepetition << ' '
         << image.ushSet << ' ' << image.ushSeg;
  return output;
}
/*
VD11_Image& VD11_Image::operator=(const VD11_Image &rhs)
{
  this->ushSamplesInScan      = rhs.ushSamplesInScan;
  this->ushUsedChannels       = rhs.ushUsedChannels;
  this->ushKSpaceCentreColumn = rhs.ushKSpaceCentreColumn;
  this->fReadOutOffcentre     = rhs.fReadOutOffcentre;
  this->ulTimeSinceLastRF     = rhs.ulTimeSinceLastRF;
  this->ushKSpaceCentreLineNo = rhs.ushKSpaceCentreLineNo;
  this->ushKSpaceCentrePartitionNo = rhs.ushKSpaceCentrePartitionNo;

  set_ushLine(rhs.ushLine);
  set_ushAcquisition(rhs.ushAcquisition);
  set_ushSlice(rhs.ushSlice);
  set_ushPartition(rhs.ushPartition);
  set_ushEcho(rhs.ushEcho);
  set_ushPhase(rhs.ushPhase);
  set_ushRepetition(rhs.ushRepetition);
  set_ushSet(rhs.ushSet);
  set_ushSeg(rhs.ushSeg);

  return *this;
}
*/

int VD11_Image::operator==(const VD11_Image &rhs) const
{
  //if( this->ushLine         != rhs.ushLine)         return 0;
  if( this->ushAcquisition  != rhs.ushAcquisition)  return 0;
  if( this->ushSlice        != rhs.ushSlice)        return 0;
  if( this->ushPartition    != rhs.ushPartition)    return 0;
  if( this->ushEcho         != rhs.ushEcho)         return 0;
  if( this->ushPhase        != rhs.ushPhase)        return 0;
  if( this->ushRepetition   != rhs.ushRepetition)   return 0;
  if( this->ushSet          != rhs.ushSet)          return 0;
  if( this->ushSeg          != rhs.ushSeg)          return 0;

  return 1;
}

// This function is required for built-in STL list functions like sort
int VD11_Image::operator<(const VD11_Image &rhs) const
{
  if( this->ushLine == rhs.ushLine && this->ushAcquisition == rhs.ushAcquisition && this->ushSlice == rhs.ushSlice
      && this->ushPartition == rhs.ushPartition && this->ushEcho == rhs.ushEcho && this->ushPhase == rhs.ushPhase
      && this->ushRepetition == rhs.ushRepetition && this->ushSet == rhs.ushSet && this->ushSeg < rhs.ushSeg) return 1;
  if( this->ushLine == rhs.ushLine && this->ushAcquisition == rhs.ushAcquisition && this->ushSlice == rhs.ushSlice
      && this->ushPartition == rhs.ushPartition && this->ushEcho == rhs.ushEcho && this->ushPhase == rhs.ushPhase
      && this->ushRepetition == rhs.ushRepetition && this->ushSet < rhs.ushSet) return 1;
  if( this->ushLine == rhs.ushLine && this->ushAcquisition == rhs.ushAcquisition && this->ushSlice == rhs.ushSlice
      && this->ushPartition == rhs.ushPartition && this->ushEcho == rhs.ushEcho && this->ushPhase == rhs.ushPhase
      && this->ushRepetition < rhs.ushRepetition) return 1;
  if( this->ushLine == rhs.ushLine && this->ushAcquisition == rhs.ushAcquisition && this->ushSlice == rhs.ushSlice
      && this->ushPartition == rhs.ushPartition && this->ushEcho == rhs.ushEcho && this->ushPhase < rhs.ushPhase) return 1;
  if( this->ushLine == rhs.ushLine && this->ushAcquisition == rhs.ushAcquisition && this->ushSlice == rhs.ushSlice
      && this->ushPartition == rhs.ushPartition && this->ushEcho < rhs.ushEcho) return 1;
  if( this->ushLine == rhs.ushLine && this->ushAcquisition == rhs.ushAcquisition && this->ushSlice == rhs.ushSlice
      && this->ushPartition < rhs.ushPartition) return 1;
  if( this->ushLine == rhs.ushLine && this->ushAcquisition == rhs.ushAcquisition && this->ushSlice < rhs.ushSlice) return 1;
  if( this->ushLine == rhs.ushLine && this->ushAcquisition < rhs.ushAcquisition) return 1;
  if( this->ushLine < rhs.ushLine ) return 1;

  return 0;
}

#endif  //VD11_IMAGE

//==================================================================================
//
//                    ReadSiemensVD11.HPP
//
//==================================================================================

#ifndef READSIEMENSVD11_HPP
#define READSIEMENSVD11_HPP

class ReadSiemensVD11
{
  public:
    /*
     % inlining of evalInfoMask
    */
    struct sMask_info
    {
        double MDH_ACQEND;
        double MDH_RTFEEDBACK;
        double MDH_HPFEEDBACK;
        double MDH_REFPHASESTABSCAN;
        double MDH_PHASESTABSCAN;
        double MDH_PHASCOR;
        double MDH_PATREFSCAN;
        double MDH_PATREFANDIMASCAN;
        double MDH_REFLECT;
        double MDH_RAWDATACORRECTION;
        double MDH_SIGNREV;
        double MDH_NOISEADJSCAN;
        double MDH_IMASCAN;
    };

    //%%--------------------------------------------------------------------------%%
    //%% Definition of loop counter structure                                     %%
    //%% Note: any changes of this structure affect the corresponding swapping    %%
    //%%       method of the measurement data header proxy class (MdhProxy)       %%
    //%%--------------------------------------------------------------------------%%
    struct sLoopCounter
    {
      unsigned short int ushLine;                             //line index
      unsigned short int ushAcquisition;                      //acquisition index
      unsigned short int ushSlice;                            //slice index
      unsigned short int ushPartition;                        //partition index
      unsigned short int ushEcho;                             //echo index
      unsigned short int ushPhase;                            //phase index
      unsigned short int ushRepetition;                       //measurement repeat index
      unsigned short int ushSet;                              //set index
      unsigned short int ushSeg;                              //segment index  (for TSE)
      unsigned short int ushIda;                              //IceDimension a index
      unsigned short int ushIdb;                              //IceDimension b index
      unsigned short int ushIdc;                              //IceDimension c index
      unsigned short int ushIdd;                              //IceDimension d index
      unsigned short int ushIde;                              //IceDimension e index
    };

    //%%--------------------------------------------------------------------------%%
    //%%  Definition of measurement data header                                   %%
    //%%--------------------------------------------------------------------------%%
    struct sMDH
    {
        unsigned int          aulEvalInfoMask[2];                   // evaluation info mask field
        unsigned short int    ushSamplesInScan;                     //# of samples acquired in scan
        unsigned short int    ushUsedChannels;                      //# of channels used in scan
        sLoopCounter          sLC;                                  //loop counters
        char                  ushDummy1[4];                            //for swapping
        unsigned short int    ushKSpaceCentreColumn;                //centre of echo
        unsigned short int    ushDummy;                             //for swapping
        float                 fReadOutOffcentre;                    //ReadOut offcenter value
        unsigned              ulTimeSinceLastRF;                    //Sequence time stamp since last RF pulse
        unsigned short int    ushKSpaceCentreLineNo;                //number of K-space centre line
        unsigned short int    ushKSpaceCentrePartitionNo;           //number of K-space centre partition    };
    };


// =====================================================================================


    //! \brief Default constructor.
    //!
    //! The default constructor creates an empty ReadSiemens.
    ReadSiemensVD11(const char* filename)
    {
      setFileName(filename);

      Init();
    }

    //! \brief Destructor.
    virtual ~ReadSiemensVD11()
    {
      myfile_.close();
    }

    //! \brief sets the filename
    void setFileName(const char* filename)
    {
      if(filename != "")
        file_name_ = filename;
      else
        std::cerr<<"wrong filename!";
    }

    int readfile(bool dicom=false);

    //! \brief returns the rawdata for given image - index
    void getRawdata(unsigned& num_rows, unsigned& num_columns, unsigned& num_coils,std::vector<std::complex<float> >& data,
                    unsigned short int acq, unsigned short int sli, unsigned short int par, unsigned short int echo,
                    unsigned short int pha, unsigned short int rep, unsigned short int set, unsigned short int seg );

    void getRawdata_info(unsigned short int& acq, unsigned short int& sli, unsigned short int& par, unsigned short int& echo,
                    unsigned short int& pha, unsigned short int& rep, unsigned short int& set, unsigned short int& seg );

    int checkRawdata_info(unsigned short int acq, unsigned short int sli, unsigned short int par, unsigned short int echo,
                    unsigned short int pha, unsigned short int rep, unsigned short int set, unsigned short int seg);



  private:

    void Init();
    int openfile(std::ifstream &myfile, const char * file_name);
    unsigned int getMDHinfo(unsigned int cPos, sMDH& mdh_data);
    void readrawdata_line(unsigned int cPos, const sMDH& mdh_data);

    std::ifstream myfile_;

    const char *file_name_;

    std::list<VD11_Image> image_list_;
    sMask_info mask_;

    int fileerror_;

};

#endif // READSIEMENSVD11_HPP

//==================================================================================
//
//                    ReadSiemensVD11.CPP
//
//==================================================================================


//#include "ReadSiemensVD11.hpp"

//===========================================================================
//
//              M E T H O D E S
//
//===========================================================================

// ---------------------------------------------------------------
//! \brief Init()
//!  special initialize for class ReadSiemens.
// ---------------------------------------------------------------

void ReadSiemensVD11::Init()
{
  fileerror_ = openfile(myfile_,file_name_);

  mask_.MDH_ACQEND = 0;
}


//-------------------------------------------------------------
//  Open the File
//
//  call example:
//             std::fstream myfile;
//             init_matrixlog(myfile,"meas.dat");
//-------------------------------------------------------------
int ReadSiemensVD11::openfile(std::ifstream &myfile, const char * file_name)
{
  myfile.open(file_name, std::fstream::in | std::fstream::binary);

  if (!myfile.is_open())
  {
      std::cerr << "File not found: " << file_name << std::endl;
      return -1;
  }
  return 0;
}

// ---------------------------------------------------------------
//! \brief Init()
//!  special initialize for class ReadSiemens.
// ---------------------------------------------------------------
int ReadSiemensVD11::readfile(bool dicom)
{
  sMDH mdh_data;
  unsigned long cPos;

 if (fileerror_ != 0)
    return -1;

  myfile_.seekg (0, std::ios::end);
  long file_length = myfile_.tellg();
  myfile_.seekg (0, std::ios::beg);  // back to the start of the file

  if(dicom == false)
  {
    unsigned int firstInt;
    unsigned int secondInt;
    myfile_.read((char*)&firstInt, sizeof(firstInt));
    myfile_.read((char*)&secondInt, sizeof(secondInt));

    myfile_.seekg (16, std::ios::beg);            //move 16bytes from beginning

    unsigned long meas_offset;
    unsigned long meas_length;
    myfile_.read((char*)&meas_offset, sizeof(meas_offset));
    myfile_.read((char*)&meas_length, sizeof(meas_length));

    myfile_.seekg(meas_offset, std::ios::beg);

    unsigned int hdr_len;
    myfile_.read((char*)&hdr_len, sizeof(hdr_len));

    unsigned long dat_start = meas_offset + hdr_len;
    cPos = dat_start;
  }
  else
    cPos = 0;

  myfile_.seekg(cPos, std::ios::beg);   //set to end of Measurement-Header = start of first Scan-Header (MDH)

  while (long(myfile_.tellg())+128 < file_length) // % fail-safe; in case we miss MDH_ACQEND
  {

    unsigned int nBytes = getMDHinfo(cPos, mdh_data);

    if(mask_.MDH_IMASCAN == 1)
    {
      readrawdata_line(cPos, mdh_data);
    }

    if (mask_.MDH_ACQEND)
    {
        break;
    }

    cPos = cPos + nBytes;

    if(cPos > file_length)
        break;
  }

  /*VD11_Image image_obj;
  std::cout<<"\n get_ushLine_length: "<<image_obj.get_ushLine_length();
  std::cout<<"\n get_ushSlice_length: "<<image_obj.get_ushSlice_length();
  std::cout<<"\n get_ushAcquisition_length: "<<image_obj.get_ushAcquisition_length();
  std::cout<<"\n get_ushEcho_length: "<<image_obj.get_ushEcho_length();
  std::cout<<"\n get_ushPartition_length: "<<image_obj.get_ushPartition_length();
  std::cout<<"\n get_ushPhase_length: "<<image_obj.get_ushPhase_length();
  std::cout<<"\n get_ushRepetition_length: "<<image_obj.get_ushRepetition_length();
  std::cout<<"\n get_ushSeg_length: "<<image_obj.get_ushSeg_length();
  std::cout<<"\n get_ushSet_length: "<<image_obj.get_ushSet_length();
*/

  return 0;


}


// ---------------------------------------------------------------
//! \brief unsigned int getMDHinfo(unsigned int cPos)
//!
//! %% Sort RawData according to number of used channels and PE Lines
//! [mdh mask nBytes] = getMDHinfo(fid,cPos);
// ---------------------------------------------------------------
unsigned int ReadSiemensVD11::getMDHinfo(unsigned int cPos, sMDH& mdh_data)
{
    unsigned szScanHeader    = 192; //% [bytes]
    unsigned szChannelHeader = 32;  //% [bytes]

    //% inlining of readScanHeader
    myfile_.seekg(cPos+40, std::ios::beg);

    myfile_.read((char*)&mdh_data, sizeof(mdh_data));

    mask_.MDH_ACQEND             = std::min( (int(mdh_data.aulEvalInfoMask[0]) & int(pow(double(2),int(0)))) , int(1));
    mask_.MDH_RTFEEDBACK         = std::min( (int(mdh_data.aulEvalInfoMask[0]) & int(pow(double(2),int(1)))) , int(1));
    mask_.MDH_HPFEEDBACK         = std::min( (int(mdh_data.aulEvalInfoMask[0]) & int(pow(double(2),int(2)))) , int(1));
    mask_.MDH_REFPHASESTABSCAN   = std::min( (int(mdh_data.aulEvalInfoMask[0]) & int(pow(double(2),int(14)))) , int(1));
    mask_.MDH_PHASESTABSCAN      = std::min( (int(mdh_data.aulEvalInfoMask[0]) & int(pow(double(2),int(15)))) , int(1));
    mask_.MDH_PHASCOR            = std::min( (int(mdh_data.aulEvalInfoMask[0]) & int(pow(double(2),int(21)))) , int(1));
    mask_.MDH_PATREFSCAN         = std::min( (int(mdh_data.aulEvalInfoMask[0]) & int(pow(double(2),int(22)))) , int(1));
    mask_.MDH_PATREFANDIMASCAN   = std::min( (int(mdh_data.aulEvalInfoMask[0]) & int(pow(double(2),int(23)))) , int(1));
    mask_.MDH_REFLECT            = std::min( (int(mdh_data.aulEvalInfoMask[0]) & int(pow(double(2),int(24)))) , int(1));
    mask_.MDH_RAWDATACORRECTION  = std::min( (int(mdh_data.aulEvalInfoMask[0]) & int(pow(double(2),int(10)))) , int(1));
    mask_.MDH_SIGNREV            = std::min( (int(mdh_data.aulEvalInfoMask[0]) & int(pow(double(2),int(17)))) , int(1));
    mask_.MDH_NOISEADJSCAN       = std::min( (int(mdh_data.aulEvalInfoMask[0]) & int(pow(double(2),int(25)))) , int(1));
    mask_.MDH_IMASCAN            = 1;

    if (mask_.MDH_ACQEND || mask_.MDH_RTFEEDBACK || mask_.MDH_HPFEEDBACK || mask_.MDH_REFPHASESTABSCAN
            || mask_.MDH_PHASESTABSCAN || mask_.MDH_PHASCOR || mask_.MDH_NOISEADJSCAN)
        mask_.MDH_IMASCAN = 0;

    unsigned int nBytes = szScanHeader + mdh_data.ushUsedChannels * (szChannelHeader + 2*4* mdh_data.ushSamplesInScan);

    return nBytes;
}


// ---------------------------------------------------------------
//! \brief unsigned int readrawdata_line(unsigned int cPos, sMDH mdh_data)
//!
//! read rawdata with corresponding mdh-info from startposition cPos
//!
// ---------------------------------------------------------------
void ReadSiemensVD11::readrawdata_line(unsigned int cPos, const sMDH& mdh_data)
{
  unsigned szScanHeader    = 192; //% [bytes]
  unsigned szChannelHeader = 32;  //% [bytes]

  //VB-version
  //szScanHeader    = 128; //% [bytes]
  //szChannelHeader = 0;  //% [bytes]

  std::vector<std::complex<float> > rawdata;
  VD11_Image image_obj;
  std::list<VD11_Image> image_list;
  long datasize = sizeof(std::complex<float>) * mdh_data.ushSamplesInScan;

  rawdata.resize(mdh_data.ushSamplesInScan,std::complex<float>(0,0));

  image_obj.ushSamplesInScan      = mdh_data.ushSamplesInScan;
  image_obj.ushUsedChannels       = mdh_data.ushUsedChannels;
  image_obj.ushKSpaceCentreColumn = mdh_data.ushKSpaceCentreColumn;
  image_obj.fReadOutOffcentre     = mdh_data.fReadOutOffcentre;
  image_obj.ulTimeSinceLastRF     = mdh_data.ulTimeSinceLastRF;
  image_obj.ushKSpaceCentreLineNo = mdh_data.ushKSpaceCentreLineNo;
  image_obj.ushKSpaceCentrePartitionNo = mdh_data.ushKSpaceCentrePartitionNo;

  image_obj.set_ushLine(mdh_data.sLC.ushLine);
  image_obj.set_ushAcquisition(mdh_data.sLC.ushAcquisition);
  image_obj.set_ushSlice(mdh_data.sLC.ushSlice);
  image_obj.set_ushPartition(mdh_data.sLC.ushPartition);
  image_obj.set_ushEcho(mdh_data.sLC.ushEcho);
  image_obj.set_ushPhase(mdh_data.sLC.ushPhase);
  image_obj.set_ushRepetition(mdh_data.sLC.ushRepetition);
  image_obj.set_ushSet(mdh_data.sLC.ushSet);
  image_obj.set_ushSeg(mdh_data.sLC.ushSeg);

  //set startposition
  myfile_.seekg(cPos+szScanHeader, std::ios::beg);

  for(int i=0; i < mdh_data.ushUsedChannels; ++i)
  {
    myfile_.seekg(szChannelHeader, std::ios::cur);

    myfile_.read((char*)&rawdata[0], datasize);

    image_obj.set_channeldata(i,rawdata);
  }

/*
  unsigned long* rawint=(unsigned long*)&rawdata[0];
  std::cout<<"\n             "<<cPos+szScanHeader+szChannelHeader<<"    "<< image_obj;
  printf ("   %#x  ",*rawint);
*/




  image_list.push_back(image_obj);
  this->image_list_.merge(image_list);

}



// ---------------------------------------------------------------
//! \brief unsigned int getMDHinfo(unsigned int cPos)
//!
//! \brief returns the rawdata for given image - index
//! in:
//! acq           ...  wanted Acquisition index
//! sli           ...  wanted Slice index
//! par           ...  wanted Partition index
//! echo          ...  wanted Echo index
//! pha           ...  wanted Phase index
//! rep           ...  wanted Repedition index
//! set           ...  wanted Set index
//! set           ...  wanted Seg index
//!
//! out:
//! num_rows      ...  return number of rows
//! num_columns   ...  return number of columns
//! num_coils     ...  return number of coils
//! data          ...  returns the rawdata
// ---------------------------------------------------------------
void ReadSiemensVD11::getRawdata(unsigned& num_rows, unsigned& num_columns, unsigned& num_coils,std::vector<std::complex<float> >& data,
                unsigned short int acq, unsigned short int sli, unsigned short int par, unsigned short int echo,
                unsigned short int pha, unsigned short int rep, unsigned short int set, unsigned short int seg )
{
  std::list<VD11_Image>::iterator i;
  i=image_list_.begin();
  VD11_Image image_obj;
  VD11_Image comp_obj;

  if(checkRawdata_info(acq,sli,par,echo,pha,rep,set,seg) < 0)
  {
    std::cerr<<"\n check check Rawdata_info index";
  }

  comp_obj.set_ushAcquisition(acq);
  comp_obj.set_ushEcho(echo);
  comp_obj.set_ushPartition(par);
  comp_obj.set_ushPhase(pha);
  comp_obj.set_ushRepetition(rep);
  comp_obj.set_ushSeg(seg);
  comp_obj.set_ushSet(set);
  comp_obj.set_ushSlice(sli);

  image_obj = *i;
  num_rows = image_obj.get_ushLine_length();
  num_columns = image_obj.ushSamplesInScan;
  num_coils = image_obj.ushUsedChannels;

  data.clear(); //force the data-vector to be emty

  for (int chx=0; chx < num_coils; ++chx )
  {
    for(i=image_list_.begin(); i != image_list_.end(); ++i)
    {
      image_obj = *i;

      if (image_obj == comp_obj)
      {
        std::vector<std::complex<float> >::iterator it_chxdata_begin;
        std::vector<std::complex<float> >::iterator it_chxdata_end;
        it_chxdata_begin = image_obj.get_channeldata(chx).begin();
        it_chxdata_end = image_obj.get_channeldata(chx).end();

        data.insert(data.end(), it_chxdata_begin, it_chxdata_end);
      }
    }
  }

}


// ---------------------------------------------------------------
//! \brief void getRawdata_info
//!
//! \brief returns the readed rawdata-info index
//! out:
//! acq           ...  Acquisition index
//! sli           ...  Slice index
//! par           ...  Partition index
//! echo          ...  Echo index
//! pha           ...  Phase index
//! rep           ...  Repedition index
//! set           ...  Set index
//! set           ...  Seg index
// ---------------------------------------------------------------
void ReadSiemensVD11::getRawdata_info(unsigned short int& acq, unsigned short int& sli, unsigned short int& par, unsigned short int& echo,
                unsigned short int& pha, unsigned short int& rep, unsigned short int& set, unsigned short int& seg)
{
  VD11_Image image_obj;
  acq = image_obj.get_ushAcquisition_length();
  sli = image_obj.get_ushSlice_length();
  par = image_obj.get_ushPartition_length();
  echo = image_obj.get_ushEcho_length();
  pha = image_obj.get_ushPhase_length();
  rep = image_obj.get_ushRepetition_length();
  set = image_obj.get_ushSet_length();
  seg = image_obj.get_ushSeg_length();
}

// ---------------------------------------------------------------
//! \brief void checkRawdata_info
//!
//! \brief checks for correct index
//! in:
//! acq           ...  Acquisition index
//! sli           ...  Slice index
//! par           ...  Partition index
//! echo          ...  Echo index
//! pha           ...  Phase index
//! rep           ...  Repedition index
//! set           ...  Set index
//! set           ...  Seg index
//! out:
//!   returns -1 if index is out of possible range
//!   returns 0  if index is ok
// ---------------------------------------------------------------
int ReadSiemensVD11::checkRawdata_info(unsigned short int acq, unsigned short int sli, unsigned short int par, unsigned short int echo,
                unsigned short int pha, unsigned short int rep, unsigned short int set, unsigned short int seg)
{
  VD11_Image image_obj;
  if (acq < 0 || acq > image_obj.get_ushAcquisition_length())
    return -1;
  if (sli < 0 || sli > image_obj.get_ushSlice_length())
    return -1;
  if (par < 0 || par > image_obj.get_ushPartition_length())
    return -1;
  if (echo < 0 || echo > image_obj.get_ushEcho_length())
    return -1;
  if (pha < 0 || pha > image_obj.get_ushPhase_length())
    return -1;
  if (rep < 0 || rep > image_obj.get_ushRepetition_length())
    return -1;
  if (set < 0 || set > image_obj.get_ushSet_length())
    return -1;
  if (seg < 0 || seg > image_obj.get_ushSeg_length())
    return -1;

  return 0;
}
