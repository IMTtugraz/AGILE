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


#ifndef READSIEMENSVD11_HPP
#define READSIEMENSVD11_HPP

#include "agile/gpu_matrix.hpp"
#include "agile/gpu_complex.hpp"
#include "agile/gpu_vector.hpp"

#include <fstream>
#include <limits>
#include <iostream>
#include <ostream>
#include <vector>
#include <list>
#include <math.h>
#include <stdexcept>
#include <stdint.h>

namespace agile
{
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
      ~VD11_Image(){}
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
      unsigned short int   ushKSpaceCentrePartitionNo;           //number of K-space centre partition    

      // -----------------------------
      // SET - METHODS:
      // -----------------------------
      void set_channeldata(unsigned short int channelid, const std::vector<std::complex<float> >& rawdata)
      {
       sChanneldata chdata;
        std::vector<sChanneldata>::iterator it;

        chdata.rawdata = rawdata;

        it = this->channeldata_.begin();

        this->channeldata_.insert(it+channelid,chdata);
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

      unsigned short int getUshPhase() 
      {
        return this->ushPhase;
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

      void set_ushKSpaceCentreColumn(unsigned short int ushKSpaceCentreColumn)
      {
        this->ushKSpaceCentreColumn = ushKSpaceCentreColumn;
        ushConsistentKSpaceCentreColumn = ushKSpaceCentreColumn;
      }

      void set_ushKSpaceCentreLineNo(unsigned short int ushKSpaceCentreLineNo)
      {
        this->ushKSpaceCentreLineNo = ushKSpaceCentreLineNo;
        ushConsistentKSpaceCentreLineNo = ushKSpaceCentreLineNo;
      }

      void set_num_columns(unsigned  num_columns)
      {
        this->_num_columns = num_columns;
      }

      void set_num_channels(unsigned  num_channels)
      {
        this->_num_channels = num_channels;
      }

      unsigned short int get_ushLine()
      {
        return this->ushLine;
      }

      void set_all_length_zero()
      {
        ushLine_length = 0;                 //line length
        ushAcquisition_length = 0;          //acquisition length
        ushSlice_length = 0;                //slice length
        ushPartition_length = 0;            //partition length
        ushEcho_length = 0;                 //echo length
        ushPhase_length = 0;                //phase length
        ushRepetition_length = 0;           //measurement repeat length
        ushSet_length = 0;                  //set length
        ushSeg_length = 0;                  //segment length  (for TSE)
        ushConsistentKSpaceCentreColumn = 0;
        ushConsistentKSpaceCentreLineNo = 0;
      }



    // -----------------------------
    // GET - METHODS:
    // -----------------------------
      std::vector<std::complex<float> >& get_channeldata(unsigned short int channelid)
      {

        if(this->channeldata_.size() <= channelid)
          std::cerr<<" wrong channelid ";

        return this->channeldata_[channelid].rawdata;
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

      unsigned short int get_num_columns(void)
      {
        return this->_num_columns;
      }

      unsigned short int get_num_channels(void)
      {
        return this->_num_channels;
      }

      unsigned short int getUshConsistentKSpaceCentreLineNo()
      {
        return this->ushConsistentKSpaceCentreLineNo;
      }

      unsigned short int getUshConsistentKSpaceCentreColumn()
      {
        return this->ushConsistentKSpaceCentreColumn;
      }

    private:
      std::vector<sChanneldata > channeldata_; //data % holding the raw data

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
      static unsigned short int   ushConsistentKSpaceCentreColumn;
      static unsigned short int   ushConsistentKSpaceCentreLineNo;

      static unsigned _num_columns;                           //ushSamplesInScan;
      static unsigned _num_channels;                          //ushUsedChannels;
  };


//==================================================================================
//
//                    ReadSiemensVD11
//                    Work by H. Heigl, extended by A. Schwarzl
//                    see: https://github.com/cjohnevans/Gannet2.0/blob/master/mapVBVD.m
//
//==================================================================================

  class ReadSiemensVD11
  {
    public:
      /*
       % inlining of evalInfoMask
      */
      struct EvalInfoMask
      {
          bool ACQEND;
          bool RTFEEDBACK;
          bool HPFEEDBACK;
          bool ONLINE;
          bool OFFLINE;
          bool SYNCDATA;
          bool six;
          bool seven;
          bool LASTSCANINCONCAT;
          bool nine;
          bool RAWDATACORRECTION;
          bool LASTSCANINMEAS;
          bool SCANSSCALEFACTOR;
          bool SECONDHADAMARPULSE;
          bool REFPHASESTABSCAN;
          bool PHASESTABSCAN;
          bool D3FFT;
          bool SIGNREV;
          bool PHASEFFT;
          bool SWAPPED;
          bool POSTSHAREDLINE;
          bool PHASCOR;
          bool PATREFSCAN;
          bool PATREFANDIMASCAN;
          bool REFLECT;
          bool NOISEADJSCAN;
          bool SHARENOW;
          bool LASTMEASUREDLINE;
          bool FIRSTSCANINSLICE;
          bool LASTSCANINSLICE;
          bool TREFFECTIVEBEGIN;
          bool TREFFECTIVEEND;
          bool thirty_two;
          bool thirty_three;
          bool thirty_four;
          bool thirty_five;
          bool thirty_six;
          bool thirty_seven;
          bool thirty_eight;
          bool thirty_nine;
          bool FIRST_SCAN_IN_BLADE;
          bool LAST_SCAN_IN_BLADE;
          bool LAST_BLADE_IN_TR;
          bool RETRO_LASTPHASE;
          bool RETRO_ENDOFMEAS;
          bool RETRO_REPEATTHISHEARTBEAT;
          bool RETRO_REPEATPREVHEARTBEAT;
          bool RETRO_ABORTSCANNOW;
          bool RETRO_LASTHEARTBEAT;
          bool RETRO_DUMMYSCAN;
          bool RETRO_ARRDETDISABLED;

          // TODO check existence
          bool IMASCAN;

          EvalInfoMask() {
          };

          EvalInfoMask(const unsigned int data[2])
          {
            ACQEND             = this->getMaskBitValue(int(data[0]),0);
            RTFEEDBACK         = this->getMaskBitValue(int(data[0]),1);
            HPFEEDBACK         = this->getMaskBitValue(int(data[0]),2);
            ONLINE             = this->getMaskBitValue(int(data[0]),3);
            OFFLINE            = this->getMaskBitValue(int(data[0]),4);
            SYNCDATA           = this->getMaskBitValue(int(data[0]),5);
            six                = this->getMaskBitValue(int(data[0]),6);
            seven              = this->getMaskBitValue(int(data[0]),7);
            LASTSCANINCONCAT   = this->getMaskBitValue(int(data[0]),8);
            nine               = this->getMaskBitValue(int(data[0]),9);
            RAWDATACORRECTION  = this->getMaskBitValue(int(data[0]),10);
            LASTSCANINMEAS     = this->getMaskBitValue(int(data[0]),11);
            SCANSSCALEFACTOR   = this->getMaskBitValue(int(data[0]),12);
            SECONDHADAMARPULSE = this->getMaskBitValue(int(data[0]),13);
            REFPHASESTABSCAN   = this->getMaskBitValue(int(data[0]),14);
            PHASESTABSCAN      = this->getMaskBitValue(int(data[0]),15);
            D3FFT              = this->getMaskBitValue(int(data[0]),16);
            SIGNREV            = this->getMaskBitValue(int(data[0]),17);
            PHASEFFT           = this->getMaskBitValue(int(data[0]),18);
            SWAPPED            = this->getMaskBitValue(int(data[0]),19);
            POSTSHAREDLINE     = this->getMaskBitValue(int(data[0]),20);
            PHASCOR            = this->getMaskBitValue(int(data[0]),21);
            PATREFSCAN         = this->getMaskBitValue(int(data[0]),22);
            PATREFANDIMASCAN   = this->getMaskBitValue(int(data[0]),23);
            REFLECT            = this->getMaskBitValue(int(data[0]),24);
            NOISEADJSCAN       = this->getMaskBitValue(int(data[0]),25);
            SHARENOW           = this->getMaskBitValue(int(data[0]),26);
            LASTMEASUREDLINE   = this->getMaskBitValue(int(data[0]),27);
            FIRSTSCANINSLICE   = this->getMaskBitValue(int(data[0]),28);
            LASTSCANINSLICE    = this->getMaskBitValue(int(data[0]),29);
            TREFFECTIVEBEGIN   = this->getMaskBitValue(int(data[0]),30);
            TREFFECTIVEEND     = this->getMaskBitValue(int(data[0]),31);
            thirty_two         = this->getMaskBitValue(int(data[1]),0);
            thirty_three       = this->getMaskBitValue(int(data[1]),1);
            thirty_four        = this->getMaskBitValue(int(data[1]),2);
            thirty_five        = this->getMaskBitValue(int(data[1]),3);
            thirty_six         = this->getMaskBitValue(int(data[1]),4);
            thirty_seven       = this->getMaskBitValue(int(data[1]),5);
            thirty_eight       = this->getMaskBitValue(int(data[1]),6);
            thirty_nine        = this->getMaskBitValue(int(data[1]),7);
            FIRST_SCAN_IN_BLADE= this->getMaskBitValue(int(data[1]),8);
            LAST_SCAN_IN_BLADE = this->getMaskBitValue(int(data[1]),9);
            LAST_BLADE_IN_TR   = this->getMaskBitValue(int(data[1]),10);
            RETRO_LASTPHASE    = this->getMaskBitValue(int(data[1]),13);
            RETRO_ENDOFMEAS    = this->getMaskBitValue(int(data[1]),14);
            RETRO_REPEATTHISHEARTBEAT = this->getMaskBitValue(int(data[1]),15);
            RETRO_REPEATPREVHEARTBEAT = this->getMaskBitValue(int(data[1]),16);
            RETRO_ABORTSCANNOW  = this->getMaskBitValue(int(data[1]),17);
            RETRO_LASTHEARTBEAT = this->getMaskBitValue(int(data[1]),18);
            RETRO_DUMMYSCAN     = this->getMaskBitValue(int(data[1]),19);
            RETRO_ARRDETDISABLED= this->getMaskBitValue(int(data[1]),20);

            IMASCAN            = 1;
          }

        private: 
          int getMaskBitValue(unsigned int mask, int bit)
          {
            return std::min( (uint32_t(mask) & uint32_t(std::pow(double(2),int(bit)))) , 1u);
          }

      };

      //%%--------------------------------------------------------------------------%%
      //%% Definition of loop counter structure                                     %%
      //%% Note: any changes of this structure affect the corresponding swapping    %%
      //%%       method of the measurement data header proxy class (MdhProxy)       %%
      //%%--------------------------------------------------------------------------%%
      struct LoopCounter
      {
        unsigned short int Line;                             //line index
        unsigned short int Acquisition;                      //acquisition index
        unsigned short int Slice;                            //slice index
        unsigned short int Partition;                        //partition index
        unsigned short int Echo;                             //echo index
        unsigned short int Phase;                            //phase index
        unsigned short int Repetition;                       //measurement repeat index
        unsigned short int Set;                              //set index
        unsigned short int Seg;                              //segment index  (for TSE)
        unsigned short int Ida;                              //IceDimension a index
        unsigned short int Idb;                              //IceDimension b index
        unsigned short int Idc;                              //IceDimension c index
        unsigned short int Idd;                              //IceDimension d index
        unsigned short int Ide;                              //IceDimension e index
      };

      //%%--------------------------------------------------------------------------%%
      //%%  Definition of measurement data header                                   %%
      //%%--------------------------------------------------------------------------%%
      struct sMDH
      {
          unsigned int          DMALengthAndFlags;
          unsigned int          MeasUID;
          unsigned int          ScanCounter;
          unsigned int          TimeStamp;
          unsigned int          PMUTimeStamp;
          unsigned short        SystemType;
          unsigned short        PTABPosDelay;
          unsigned int          PTABPosX;
          unsigned int          PTABPosY;
          unsigned int          PTABPosZ;
          unsigned int          Reserved1;
          unsigned int          sEvalInfoMask[2];                   // evaluation info mask field
          unsigned short int    SamplesInScan;                     //# of samples acquired in scan
          unsigned short int    UsedChannels;                      //# of channels used in scan
          LoopCounter           sLoopCounter;                                  //loop counters
          char                  sCutOffData[4];                            //for swapping
          unsigned short int    KSpaceCentreColumn;                //centre of echo
          unsigned short int    CoilSelect;                //centre of echo
          float                 ReadOutOffcentre;                    //ReadOut offcenter value
          unsigned int          TimeSinceLastRF;                    //Sequence time stamp since last RF pulse
          unsigned short int    KSpaceCentreLineNo;                //number of K-space centre line
          unsigned short int    KSpaceCentrePartitionNo;           //number of K-space centre partition 
          float                 SliceData[7];                             //for swapping
          unsigned int          IceProgramPara[12];
          unsigned int          ReservedPara[2];
          unsigned short int    ApplicationCounter;
          unsigned short int    ApplicationMask;
          unsigned int          CRC;

          /** 
           * \brief DMALength is contained in first 3 byte of field
           */
          unsigned int getDMALength() const
          {
            return ((DMALengthAndFlags << 8) >> 8);
          }

          EvalInfoMask getEvalInfoMask() const
          {
            return EvalInfoMask(sEvalInfoMask);
          }

      };

  // =====================================================================================

      //! \brief default Constructor.
      //!
      //! The default constructor creates an empty ReadSiemens.
      ReadSiemensVD11()
      {
      }

      //! \brief Constructor.
      //!
      //! The constructor creates an ReadSiemens Obj.
      explicit ReadSiemensVD11(std::string filename)
      {
        if(filename != "")
          file_name_ = filename;
        else
          std::cerr<<"wrong filename!";

        Init();
      }

      ReadSiemensVD11(const ReadSiemensVD11& obj)
      {
        this->file_name_ = obj.file_name_;
        this->Init();
      }

      ReadSiemensVD11& operator=(const ReadSiemensVD11& obj)
      {
        this->file_name_ = obj.file_name_;
        this->Init();
        return *this;
      }

      //! \brief Destructor.
      virtual ~ReadSiemensVD11()
      {
        VD11_Image image;

        image.set_all_length_zero();

        myfile_.close();
      }

      //! \brief sets the filename
      void setFile(std::string filename)
      {
        if(filename != "")
          file_name_ = filename;
        else
          std::cerr<<"wrong filename!";

        //close old myfile_ and init image-vals to zero
        VD11_Image image;
        image.set_all_length_zero();
        myfile_.close();

        Init();
      }

      int readfile(bool dicom=false, bool readdata=true);

      int readfile_sizeinfo(unsigned& num_coils, unsigned& nBytes, bool dicom=false);

      //! \brief returns the rawdata for given image - index
      void getRawdata(std::vector<std::complex<float> >& data,
                      unsigned short int acq, unsigned short int sli, unsigned short int par, unsigned short int echo,
                      unsigned short int pha, unsigned short int rep, unsigned short int set, unsigned short int seg );

      void getRawdata_info(unsigned int& num_rows, unsigned int& num_columns, unsigned int& num_coils,
                           unsigned short int& acq, unsigned short int& sli, unsigned short int& par, unsigned short int& echo,
                           unsigned short int& pha, unsigned short int& rep, unsigned short int& set, unsigned short int& seg);

      int checkRawdata_info(unsigned short int acq, unsigned short int sli, unsigned short int par, unsigned short int echo,
                      unsigned short int pha, unsigned short int rep, unsigned short int set, unsigned short int seg);


      unsigned int getMDHinfo(unsigned int cPos, sMDH& mdh_data);

      //! \brief return list of acquired raw data lines
      std::list<VD11_Image> getRawdataList() 
      {
        return image_list_;
      }

    private:

      void Init();
      int openfile(std::ifstream &myfile, std::string file_name);
      void readrawdata_line(unsigned int cPos, const sMDH& mdh_data);

      std::ifstream myfile_;

      std::string file_name_;

      std::list<VD11_Image> image_list_;
      EvalInfoMask mask_;

      int fileerror_;

  };
  
  /**
   * <<operator declarations.
   */
  std::ostream& operator<<(std::ostream& os, const ReadSiemensVD11::sMDH& obj);
  std::ostream& operator<<(std::ostream& os, const ReadSiemensVD11::EvalInfoMask& obj);
  std::ostream& operator<<(std::ostream& os, const ReadSiemensVD11::LoopCounter& obj);

}
#endif // READSIEMENSVD11_HPP
