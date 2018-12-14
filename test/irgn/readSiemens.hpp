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


#ifndef READSIEMENS_HPP
#define READSIEMENS_HPP

#include "agile/gpu_matrix.hpp"
#include "agile/gpu_complex.hpp"
#include "agile/gpu_vector.hpp"

#include <fstream>
#include <limits>
#include <iostream>
#include <vector>

class ReadSiemens
{
  public:

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
    //%%  Definition of slice vectors                                             %%
    //%%--------------------------------------------------------------------------%%
    struct sVector
    {
      float flSag;
      float flCor;
      float flTra;
    };
    struct sSliceData
    {
      sVector sSlicePosVec;                                       //slice position vector
      float aflQuaternion[4];                         //rotation matrix as quaternion
    };

    //%%--------------------------------------------------------------------------%%
    //%%  Definition of cut-off data                                              %%
    //%%--------------------------------------------------------------------------%%
    struct sCutOffData
    {
      unsigned short int ushPre;                              //write ushPre zeros at line start
      unsigned short int ushPost;                             //write ushPost zeros at line end
    };

    //%%--------------------------------------------------------------------------%%
    //%%  Definition of measurement data header                                   %%
    //%%--------------------------------------------------------------------------%%
    struct sMDH
    {
      unsigned              ulDMALength;                          //DMA length [bytes] must be
      int                   lMeasUID;                             //measurement user ID
      unsigned              ulScanCounter;                        //scan counter [1...]
      unsigned              ulTimeStamp;                          //time stamp [2.5 ms ticks since 00:00]
      unsigned              ulPMUTimeStamp;                       //PMU time stamp [2.5 ms ticks since last trigger]
      unsigned              aulEvalInfoMask[2];                   // evaluation info mask field        !!! zeros(1,MDH_NUMBEROFEVALINFOMASK),.
      unsigned short int    ushSamplesInScan;                     //# of samples acquired in scan
      unsigned short int    ushUsedChannels;                      //# of channels used in scan
      sLoopCounter          sLC;                                  //loop counters
      sCutOffData           sCutOff;                              //cut-off values
      unsigned short int    ushKSpaceCentreColumn;                //centre of echo
      unsigned short int    ushDummy;                             //for swapping
      float                 fReadOutOffcentre;                    //ReadOut offcenter value
      unsigned              ulTimeSinceLastRF;                    //Sequence time stamp since last RF pulse
      unsigned short int    ushKSpaceCentreLineNo;                //number of K-space centre line
      unsigned short int    ushKSpaceCentrePartitionNo;           //number of K-space centre partition
      unsigned short int    aushIceProgramPara[4];                //free parameter for IceProgram      // !!! zeros(1,MDH_NUMBEROFICEPROGRAMPARA)
      //unsigned short int aushIceProgramPara[MDH_NUMBEROFICEPROGRAMPARA];       //free parameter for IceProgram      // !!! zeros(1,MDH_NUMBEROFICEPROGRAMPARA)
      unsigned short int    aushFreePara[4];                      //free parameter      // !!! zeros(1,MDH_FREEHDRPARA)
      sSliceData            sSD;                                  //Slice Data
      unsigned              ulChannelId;                          //channel Id must be the last parameter
    };


    struct ADC_dimension
    {
      sLoopCounter aADC;
      unsigned  ulChannelId;
    };

    struct ADC_info
    {
      ADC_dimension dimensions;

      std::vector<std::complex<float> > adc_data_;      //NOT NEEDED ???
    };

    struct ADC_rawdata
    {
      unsigned  ulChannelId;
      std::vector<std::complex<float> > rawdata_;
    };



// =====================================================================================


    //! \brief Default constructor.
    //!
    //! The default constructor creates an empty ReadSiemens.
    ReadSiemens(const char* filename)
    {
      unsigned IDX_DIM = 20;

      setFileName(filename);

      Init();
    }

    //! \brief Destructor.
    virtual ~ReadSiemens()
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

    //! \brief returns the ADC_Info - Vector
    std::vector<ADC_info> getADC_Info()
    {
      return this->adc_info_;
    }

    //! \brief returns the sMDH - Data
    sMDH getsMDH()
    {
      return this->mdh_data_;
    }

    //! \brief returns the ADC_rawdata - Data
    std::vector<ADC_rawdata> getsADC_rawdata()
    {
      return this->adc_rawdata_;
    }


    void readfile();

  private:

    unsigned IDX_DIM;

    void Init();
    void openfile(std::ifstream &myfile, const char * file_name);
    void max_loopcounters(std::vector<ADC_info> adc_info);
    void sort_rawdata(std::vector<ADC_info> adc_info);


    std::ifstream myfile_;

    const char *file_name_;
    sMDH mdh_data_;
    std::vector<ADC_info> adc_info_;
    std::vector<ADC_rawdata> adc_rawdata_;

    //std::vector<std::complex<float> > rawdata_;
};

//===========================================================================
//
//              M E T H O D E S
//
//===========================================================================

// ---------------------------------------------------------------
//! \brief Init()
//!  special initialize for class ReadSiemens.
// ---------------------------------------------------------------

void ReadSiemens::Init()
{
  openfile(myfile_,file_name_);

}

//-------------------------------------------------------------
//  Open the File
//
//  call example:
//             std::fstream myfile;
//             init_matrixlog(myfile,"meas.dat");
//-------------------------------------------------------------
void ReadSiemens::openfile(std::ifstream &myfile, const char * file_name)
{
  myfile.open(file_name, std::fstream::in);

  if (!myfile.is_open())
  {
      std::cerr << "File not found: " << file_name << std::endl;
  }
}

// ---------------------------------------------------------------
//! \brief Init()
//!  special initialize for class ReadSiemens.
// ---------------------------------------------------------------
void ReadSiemens::readfile()
{

  bool lastADCdata = false;
  ADC_info adc_info;
  std::vector<std::complex<float> > adc_data;

  myfile_.seekg (0, std::ios::end);
  unsigned file_length = myfile_.tellg();
  myfile_.seekg (16, std::ios::beg);            //move 16bytes from beginning

  /*
  fseek(fid,16,'bof');
  % measOffset: points to beginning of header, usually at 10240 bytes
  measOffset = fread(fid,1,'uint64');
  measLength = fread(fid,1,'uint64');
  fseek(fid,measOffset,'bof');
  hdrLength  = fread(fid,1,'uint32');
  datStart   = measOffset + hdrLength;
  */
  unsigned double meas_offset;
  unsigned double meas_length;
  myfile_.read((char*)&meas_offset, sizeof(meas_offset));
  myfile_.read((char*)&meas_length, sizeof(meas_length));

  myfile_.seekg(meas_offset, std::ios::beg);

  unsigned int hdr_len;
  myfile_.read((char*)&hdr_len, sizeof(hdr_len));

  unsigned double dat_start = meas_offset + hdr_len;

  myfile_.seekg(dat_start, std::ios::beg);

  unsigned size = 0;
  unsigned oldsize = 0;
  unsigned lcounter = 0;

  while(!lastADCdata && !myfile_.eof())
  {
    myfile_.read((char*)&mdh_data_, sizeof(mdh_data_));

    adc_info.dimensions.aADC.ushLine        = mdh_data_.sLC.ushLine;
    adc_info.dimensions.aADC.ushAcquisition = mdh_data_.sLC.ushAcquisition;
    adc_info.dimensions.aADC.ushSlice       = mdh_data_.sLC.ushSlice;
    adc_info.dimensions.aADC.ushPartition   = mdh_data_.sLC.ushPartition;
    adc_info.dimensions.aADC.ushEcho        = mdh_data_.sLC.ushEcho;
    adc_info.dimensions.aADC.ushPhase       = mdh_data_.sLC.ushPhase;
    adc_info.dimensions.aADC.ushRepetition  = mdh_data_.sLC.ushRepetition;
    adc_info.dimensions.aADC.ushSet         = mdh_data_.sLC.ushSet;
    adc_info.dimensions.aADC.ushSeg         = mdh_data_.sLC.ushSeg;
    adc_info.dimensions.aADC.ushIda         = mdh_data_.sLC.ushIda;
    adc_info.dimensions.aADC.ushIdb         = mdh_data_.sLC.ushIdb;
    adc_info.dimensions.aADC.ushIdc         = mdh_data_.sLC.ushIdc;
    adc_info.dimensions.aADC.ushIdd         = mdh_data_.sLC.ushIdd;
    adc_info.dimensions.aADC.ushIde         = mdh_data_.sLC.ushIde;
    adc_info.dimensions.ulChannelId         = mdh_data_.ulChannelId;

    size = mdh_data_.ushSamplesInScan;
    adc_data.resize(size, std::complex<float>(0));
    myfile_.read((char*)&adc_data[0], size*sizeof(std::complex<float>));
    adc_info.adc_data_ = adc_data;    // NOT NEEDED ????

    ADC_rawdata adcrawdata;
    if (lcounter < mdh_data_.ushUsedChannels)
    {
      adcrawdata.rawdata_ = adc_data;
      adcrawdata.ulChannelId = mdh_data_.ulChannelId;
      this->adc_rawdata_.push_back(adcrawdata);
    }
    else
    {
      for(int i=0; i < mdh_data_.ushUsedChannels; ++i)
      {
        if(mdh_data_.ulChannelId == adc_rawdata_[i].ulChannelId)
        {
          adc_rawdata_[i].rawdata_.insert(
                                adc_rawdata_[i].rawdata_.end(),
                                adc_data.begin(), adc_data.end());
        }
      }
    }



    if (lcounter == 0)
      oldsize = size;

    if (oldsize == size)
      adc_info_.push_back(adc_info);


    //if curr_pos >= Nbytes - 384, LastADCData = 1; end
    unsigned curr_pos = myfile_.tellg();
    if(curr_pos >= file_length-384)
      lastADCdata = true;

    lcounter++;
  }
  std::cout<<"\n lcounter++; "<<lcounter;
  max_loopcounters(adc_info_);
}


// ---------------------------------------------------------------
//! \brief max_loopcounters(std::vector<ADC_info> adc_info)
//!
//! Let sMDH (mdh_data_) have the maximum values of all loop counters
// ---------------------------------------------------------------
void ReadSiemens::max_loopcounters(std::vector<ADC_info> adc_info)
{
  ADC_dimension dimensions;
  dimensions = adc_info[0].dimensions;

  for(int i=1; i<adc_info_.size(); ++i)
  {
    //(a<b)?b:a
    dimensions.aADC.ushLine = std::max(adc_info_[i].dimensions.aADC.ushLine, dimensions.aADC.ushLine);

    dimensions.aADC.ushAcquisition = std::max(adc_info_[i].dimensions.aADC.ushAcquisition, dimensions.aADC.ushAcquisition);

    dimensions.aADC.ushSlice = std::max(adc_info_[i].dimensions.aADC.ushSlice, dimensions.aADC.ushSlice);

    dimensions.aADC.ushPartition = std::max(adc_info_[i].dimensions.aADC.ushPartition, dimensions.aADC.ushPartition);

    dimensions.aADC.ushEcho = std::max(adc_info_[i].dimensions.aADC.ushEcho, dimensions.aADC.ushEcho);

    dimensions.aADC.ushPhase = std::max(adc_info_[i].dimensions.aADC.ushPhase, dimensions.aADC.ushPhase);

    dimensions.aADC.ushRepetition = std::max(adc_info_[i].dimensions.aADC.ushRepetition, dimensions.aADC.ushRepetition);

    dimensions.aADC.ushSet = std::max(adc_info_[i].dimensions.aADC.ushSet, dimensions.aADC.ushSet);

    dimensions.aADC.ushSeg = std::max(adc_info_[i].dimensions.aADC.ushSeg, dimensions.aADC.ushSeg);

    dimensions.aADC.ushIda = std::max(adc_info_[i].dimensions.aADC.ushIda, dimensions.aADC.ushIda);

    dimensions.aADC.ushIdb = std::max(adc_info_[i].dimensions.aADC.ushIdb, dimensions.aADC.ushIdb);

    dimensions.aADC.ushIdc = std::max(adc_info_[i].dimensions.aADC.ushIdc, dimensions.aADC.ushIdc);

    dimensions.aADC.ushIdd = std::max(adc_info_[i].dimensions.aADC.ushIdd, dimensions.aADC.ushIdd);

    dimensions.aADC.ushIde = std::max(adc_info_[i].dimensions.aADC.ushIde, dimensions.aADC.ushIde);

    dimensions.ulChannelId = std::max(adc_info_[i].dimensions.ulChannelId, dimensions.ulChannelId);
  }

  this->mdh_data_.ulChannelId = dimensions.ulChannelId;
  this->mdh_data_.sLC = dimensions.aADC;
}

// ---------------------------------------------------------------
//! \brief max_loopcounters(std::vector<ADC_info> adc_info)
//!
//! %% Sort RawData according to number of used channels and PE Lines
// ---------------------------------------------------------------
void ReadSiemens::sort_rawdata(std::vector<ADC_info> adc_info)
{
  ADC_info adc_info_safe;


  std::vector<std::complex<float> > rawdata;
  /*
    nCh = sMDH.ushUsedChannels;
    for ii = 1:nCh
        rawdata_temp(:,:,ii) = rawdata(ii:nCh:size(rawdata,1),:);
        loopcounters_temp1 = loopcounters(ii:nCh:size(loopcounters,1),1);
        loopcounters_temp2 = loopcounters(ii:nCh:size(loopcounters,1),3);
        loopcounters_temp3 = loopcounters(ii:nCh:size(loopcounters,1),9);
    end
  */
  unsigned short int nCh = this->mdh_data_.ushUsedChannels;

  for(int i=0; i<nCh; ++i)
  {
   // if (adc_info_safe.dimensions.ulChannelId == adc_info_[i]






  }




}



#endif // READSIEMENS_HPP
