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


#include "agile/io/readSiemensVD11.hpp"
  

namespace agile
{
  std::ostream& operator<<(std::ostream& os, const ReadSiemensVD11::LoopCounter& obj)
  {
    os << "\tLine:        " << std::setw(4) << obj.Line;
    os << "\t\tAcquisition: " << std::setw(4)<< obj.Acquisition << std::endl;
    os << "\tSlice:       " << std::setw(4)<< obj.Slice;
    os << "\t\tPartition:   " << std::setw(4) << obj.Partition << std::endl;
    os << "\tEcho:        " << std::setw(4)<< obj.Echo;
    os << "\t\tPhase:       " << std::setw(4)<< obj.Phase << std::endl;
    os << "\tRepetition:  " << std::setw(4)<< obj.Repetition;
    os << "\t\tSet:         " << std::setw(4)<< obj.Set << std::endl;
    os << "\tSeg:         " << std::setw(4)<< obj.Seg;
    os << "\t\tIda:         " << std::setw(4)<< obj.Ida << std::endl;
    os << "\tIdb:         " << std::setw(4)<< obj.Idb;
    os << "\t\tIdc:         " << std::setw(4)<< obj.Idc << std::endl;
    os << "\tIdd:         " << std::setw(4)<< obj.Idd;
    os << "\t\tIde:         " << std::setw(4)<< obj.Ide;
    return os;
  }
  
  std::ostream& operator<<(std::ostream& os, const ReadSiemensVD11::EvalInfoMask& obj)
  {
    if (obj.ACQEND             )        os << "\tACQEND                   :true " << std::endl;
    if (obj.ONLINE             )        os << "\tONLINE                   :true " << std::endl;
    if (obj.OFFLINE            )        os << "\tOFFLINE                  :true " << std::endl;
    if (obj.SYNCDATA           )        os << "\tSYNCDATA                 :true " << std::endl;
    if (obj.six                )        os << "\tsix                      :true " << std::endl;
    if (obj.seven              )        os << "\tseven                    :true " << std::endl;
    if (obj.LASTSCANINCONCAT   )        os << "\tLASTSCANINCONCAT         :true " << std::endl;
    if (obj.nine               )        os << "\tnine                     :true " << std::endl;
    if (obj.RAWDATACORRECTION  )        os << "\tRAWDATACORRECTION        :true " << std::endl;
    if (obj.LASTSCANINMEAS     )        os << "\tLASTSCANINMEAS           :true " << std::endl;
    if (obj.SCANSSCALEFACTOR   )        os << "\tSCANSSCALEFACTOR         :true " << std::endl;
    if (obj.SECONDHADAMARPULSE )        os << "\tSECONDHADAMARPULSE       :true " << std::endl;
    if (obj.REFPHASESTABSCAN   )        os << "\tREFPHASESTABSCAN         :true " << std::endl;
    if (obj.PHASESTABSCAN      )        os << "\tPHASESTABSCAN            :true " << std::endl;
    if (obj.D3FFT              )        os << "\tD3FFT                    :true " << std::endl;
    if (obj.SIGNREV            )        os << "\tSIGNREV                  :true " << std::endl;
    if (obj.PHASEFFT           )        os << "\tPHASEFFT                 :true " << std::endl;
    if (obj.SWAPPED            )        os << "\tSWAPPED                  :true " << std::endl;
    if (obj.POSTSHAREDLINE     )        os << "\tPOSTSHAREDLINE           :true " << std::endl;
    if (obj.PHASCOR            )        os << "\tPHASCOR                  :true " << std::endl;
    if (obj.PATREFSCAN         )        os << "\tPATREFSCAN               :true " << std::endl;
    if (obj.PATREFANDIMASCAN   )        os << "\tPATREFANDIMASCAN         :true " << std::endl;
    if (obj.REFLECT            )        os << "\tREFLECT                  :true " << std::endl;
    if (obj.NOISEADJSCAN       )        os << "\tNOISEADJSCAN             :true " << std::endl;
    if (obj.SHARENOW           )        os << "\tSHARENOW                 :true " << std::endl;
    if (obj.LASTMEASUREDLINE   )        os << "\tLASTMEASUREDLINE         :true " << std::endl;
    if (obj.FIRSTSCANINSLICE   )        os << "\tFIRSTSCANINSLICE         :true " << std::endl;
    if (obj.LASTSCANINSLICE    )        os << "\tLASTSCANINSLICE          :true " << std::endl;
    if (obj.TREFFECTIVEBEGIN   )        os << "\tTREFFECTIVEBEGIN         :true " << std::endl;
    if (obj.TREFFECTIVEEND     )        os << "\tTREFFECTIVEEND           :true " << std::endl;
    if (obj.thirty_two         )        os << "\tthirty_two               :true " << std::endl;
    if (obj.thirty_three       )        os << "\tthirty_three             :true " << std::endl;
    if (obj.thirty_four        )        os << "\tthirty_four              :true " << std::endl;
    if (obj.thirty_five        )        os << "\tthirty_five              :true " << std::endl;
    if (obj.thirty_six         )        os << "\tthirty_six               :true " << std::endl;
    if (obj.thirty_seven       )        os << "\tthirty_seven             :true " << std::endl;
    if (obj.thirty_eight       )        os << "\tthirty_eight             :true " << std::endl;
    if (obj.thirty_nine        )        os << "\tthirty_nine              :true " << std::endl;
    if (obj.FIRST_SCAN_IN_BLADE)        os << "\tFIRST_SCAN_IN_BLADE      :true " << std::endl;
    if (obj.LAST_SCAN_IN_BLADE )        os << "\tLAST_SCAN_IN_BLADE       :true " << std::endl;
    if (obj.LAST_BLADE_IN_TR   )        os << "\tLAST_BLADE_IN_TR         :true " << std::endl;
    if (obj.RETRO_LASTPHASE    )        os << "\tRETRO_LASTPHASE          :true " << std::endl;
    if (obj.RETRO_ENDOFMEAS    )        os << "\tRETRO_ENDOFMEAS          :true " << std::endl;
    if (obj.RETRO_REPEATTHISHEARTBEAT ) os << "\tRETRO_REPEATTHISHEARTBEAT:true " << std::endl;
    if (obj.RETRO_REPEATPREVHEARTBEAT ) os << "\tRETRO_REPEATPREVHEARTBEAT:true " << std::endl;
    if (obj.RETRO_ABORTSCANNOW  )       os << "\tRETRO_ABORTSCANNOW       :true " << std::endl;
    if (obj.RETRO_LASTHEARTBEAT )       os << "\tRETRO_LASTHEARTBEAT      :true " << std::endl;
    if (obj.RETRO_DUMMYSCAN     )       os << "\tRETRO_DUMMYSCAN          :true " << std::endl;
    if (obj.RETRO_ARRDETDISABLED)       os << "\tRETRO_ARRDETDISABLED     :true" << std::endl;
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const ReadSiemensVD11::sMDH& obj) 
  {
    os << "DMALength:      " << std::setw(10) << obj.getDMALength();
    os << "\t\tMeasUID:       " <<std::setw(10) <<  obj.MeasUID << std::endl;
    os << "ScanCounter:    " <<std::setw(10) <<  obj.ScanCounter;
    os << "\t\tTimeStamp:     " <<std::setw(10) <<  obj.TimeStamp << std::endl;
    os << "PMUTimeStamp:   " <<std::setw(10) <<  obj.PMUTimeStamp;
    os << "\t\tSystemType:    " << std::setw(10) <<  obj.SystemType << std::endl;
    os << "PTABPosDelay:   " << std::setw(10) <<  obj.PTABPosDelay;
    os << "\t\tPTABPosX:      " << std::setw(10) <<  obj.PTABPosX << std::endl;
    os << "PTABPosY:       " << std::setw(10) <<  obj.PTABPosY;
    os << "\t\tPTABPosZ:      " << std::setw(10) <<  obj.PTABPosZ << std::endl;

    os << "EvalInfoMask:   " << std::endl << obj.getEvalInfoMask() << std::endl;

    os << "SamplesInScan:  " << std::setw(5) <<  obj.SamplesInScan;
    os << "\t\tUsedChannels:     " << std::setw(5) <<  obj.UsedChannels << std::endl;
    os << "SLoopCounter:   " << std::endl <<  obj.sLoopCounter << std::endl;

    os << "SCutOffData:    " <<  obj.sCutOffData << std::endl;
    os << "KSpaceCentreColumn:    " << std::setw(10) <<  obj.KSpaceCentreColumn;
    os << "\t\tCoilSelect:            " << std::setw(10) <<  obj.CoilSelect << std::endl;
    os << "ReadOutOffcentre:      " << std::setw(10) <<  obj.ReadOutOffcentre;
    os << "\t\tTimeSinceLastRF:       " << std::setw(10) <<  obj.TimeSinceLastRF << std::endl;
    os << "KSpaceCentreLineNo:    " << std::setw(10) <<  obj.KSpaceCentreLineNo;
    os << "\t\tKSpaceCentrePartNo:    " << std::setw(10) <<  obj.KSpaceCentrePartitionNo << std::endl;
    os << "SliceData:            [" << obj.SliceData[0] 
      << "," << obj.SliceData[1] << "," << obj.SliceData[2] << "," << obj.SliceData[3] 
      << "," << obj.SliceData[4] << "," << obj.SliceData[5] << "," << obj.SliceData[6] << "]" << std::endl;

    os << "IceProgramPara:       [" << obj.IceProgramPara[0] 
      << "," << obj.IceProgramPara[1] << "," << obj.IceProgramPara[2] << "," << obj.IceProgramPara[3] 
      << "," << obj.IceProgramPara[4] << "," << obj.IceProgramPara[5] << "," << obj.IceProgramPara[6] 
      << "," << obj.IceProgramPara[7] << "," << obj.IceProgramPara[8] << "," << obj.IceProgramPara[9] 
      << "," << obj.IceProgramPara[10] << "," << obj.IceProgramPara[11]<< "]" << std::endl;

    os << "ApplicationCounter:    " << std::setw(10) <<  obj.ApplicationCounter;
    os << "\t\tApplicationMask:       " << std::setw(10) <<  obj.ApplicationMask << std::endl;

    os << "CRC:                   " << std::setw(10) <<  obj.CRC << std::endl;

    return os;
  }


  unsigned short int VD11_Image::ushLine_length = 0;
  unsigned short int VD11_Image::ushAcquisition_length = 0;
  unsigned short int VD11_Image::ushSlice_length = 0;
  unsigned short int VD11_Image::ushPartition_length = 0;
  unsigned short int VD11_Image::ushEcho_length = 0;
  unsigned short int VD11_Image::ushPhase_length = 0;
  unsigned short int VD11_Image::ushRepetition_length = 0;
  unsigned short int VD11_Image::ushSet_length = 0;
  unsigned short int VD11_Image::ushSeg_length = 0;
  unsigned  VD11_Image::_num_columns = 0;
  unsigned  VD11_Image::_num_channels = 0;
  unsigned short int VD11_Image::ushConsistentKSpaceCentreColumn = 0;
  unsigned short int VD11_Image::ushConsistentKSpaceCentreLineNo = 0;

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

  // This function is required for built-in STL list functions like sort, merge

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


//==================================================================================
//
//                    ReadSiemensVD11
//
//==================================================================================

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

    mask_.ACQEND = 0;
  }


  //-------------------------------------------------------------
  //  Open the File
  //
  //  call example:
  //             std::fstream myfile;
  //             init_matrixlog(myfile,"meas.dat");
  //-------------------------------------------------------------

  int ReadSiemensVD11::openfile(std::ifstream &myfile, std::string file_name)
  {
    myfile.open(file_name.c_str(), std::fstream::in | std::fstream::binary);

    if (!myfile.is_open())
    {
        std::cerr << "File not found: " << file_name << std::endl;
        return -1;
    }
    return 0;
  }


  // ---------------------------------------------------------------
  //! \brief int ReadSiemensVD11::readfile(bool dicom)
  //!  function: reads the file-header information and the corresponding rawdata to a vector
  //!
  //!  in: bool dicom:   true for a dicom-file
  //!                    false for a non dicom-file header (meas.dat)
  //!      bool readdata true for loading data
  //!                    false for not
  //!  return: -1 for error  ,  0 for ok  ,  1 for acqend no reached
  // ---------------------------------------------------------------

  int ReadSiemensVD11::readfile(bool dicom, bool readdata)
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

    while (long(myfile_.tellg())+128 < file_length) // % fail-safe; in case we miss ACQEND
    {
      unsigned int nBytes = getMDHinfo(cPos, mdh_data);

      if((mask_.IMASCAN == 1) && (readdata == true))
      {
        readrawdata_line(cPos, mdh_data);
      }

      if (mask_.ACQEND)
      {
        //std::cout<<"\n mdh_acqend: "<<mask_.ACQEND;
          break;
      }

      cPos = cPos + nBytes;

      if(cPos > file_length)
          break;
    }

    if (mask_.ACQEND)
      return 0;             //fileend=true;
    else
      return 1 ;            //fileend=false;

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
      myfile_.read((char*)&mdh_data, sizeof(mdh_data));

      mask_ = mdh_data.getEvalInfoMask();

      if (mask_.ACQEND || mask_.RTFEEDBACK || mask_.HPFEEDBACK || mask_.REFPHASESTABSCAN
              || mask_.PHASESTABSCAN || mask_.PHASCOR || mask_.NOISEADJSCAN)
          mask_.IMASCAN = 0;

      unsigned int nBytes = szScanHeader + mdh_data.UsedChannels * (szChannelHeader + 2*4* mdh_data.SamplesInScan);

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

    std::vector<std::complex<float> > rawdata;
    VD11_Image image_obj;
    std::list<VD11_Image> image_list;
    long datasize = sizeof(std::complex<float>) * mdh_data.SamplesInScan;

    rawdata.resize(mdh_data.SamplesInScan,std::complex<float>(0,0));

    image_obj.ushSamplesInScan      = mdh_data.SamplesInScan;
    image_obj.ushUsedChannels       = mdh_data.UsedChannels;
    image_obj.ushKSpaceCentrePartitionNo= mdh_data.KSpaceCentrePartitionNo;
    image_obj.fReadOutOffcentre     = mdh_data.ReadOutOffcentre;
    image_obj.ulTimeSinceLastRF     = mdh_data.TimeSinceLastRF;
    image_obj.set_ushKSpaceCentreLineNo(mdh_data.KSpaceCentreLineNo);
    image_obj.set_ushKSpaceCentreColumn(mdh_data.KSpaceCentreColumn);

    image_obj.set_ushLine(mdh_data.sLoopCounter.Line);
    image_obj.set_ushAcquisition(mdh_data.sLoopCounter.Acquisition);
    image_obj.set_ushSlice(mdh_data.sLoopCounter.Slice);
    image_obj.set_ushPartition(mdh_data.sLoopCounter.Partition);
    image_obj.set_ushEcho(mdh_data.sLoopCounter.Echo);
    image_obj.set_ushPhase(mdh_data.sLoopCounter.Phase);
    image_obj.set_ushRepetition(mdh_data.sLoopCounter.Repetition);
    image_obj.set_ushSet(mdh_data.sLoopCounter.Set);
    image_obj.set_ushSeg(mdh_data.sLoopCounter.Seg);
    image_obj.set_num_columns(mdh_data.SamplesInScan);
    image_obj.set_num_channels(mdh_data.UsedChannels);

    //set startposition
    myfile_.seekg(cPos+szScanHeader, std::ios::beg);

    for(int i=0; i < mdh_data.UsedChannels; ++i)
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
  //! data          ...  returns the rawdata
  // ---------------------------------------------------------------

  void ReadSiemensVD11::getRawdata(std::vector<std::complex<float> >& data,
                  unsigned short int acq, unsigned short int sli, unsigned short int par, unsigned short int echo,
                  unsigned short int pha, unsigned short int rep, unsigned short int set, unsigned short int seg )
  {
    std::list<VD11_Image>::iterator i;
    i=image_list_.begin();
    VD11_Image image_obj;
    VD11_Image comp_obj;

    unsigned act_line=0;

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
    unsigned num_coils = image_obj.ushUsedChannels;

//    data.clear(); //force the data-vector to be emty
//    std::cout<<"\n data.capacity(): "<<data.capacity()<<" data.size: "<<data.size();

    int test=0;


    for (unsigned chx=0; chx < num_coils; ++chx )
    {
      i=image_list_.begin();
      act_line=0;
      while(i != image_list_.end())
      {
        image_obj = *i;

//        std::cout<<"\n i_obj: "<<image_obj;
//        std::cout<<"\n comp_obj: "<<comp_obj;

        if (image_obj == comp_obj)
        {
//          std::cout<<"\n actline                : "<<act_line;
//          std::cout<<"\n image_obj.get_ushLine(): "<<image_obj.get_ushLine();

          if(act_line < image_obj.get_ushLine())    //for accelerated images put complex zeros in between
          {
            data.insert(data.end(),image_obj.get_channeldata(chx).size() , std::complex<float>(0,0));
            ++act_line;
          }
          else
          {
//            std::vector<std::complex<float> >::iterator it_chxdata_begin;
//            std::vector<std::complex<float> >::iterator it_chxdata_end;
//            it_chxdata_begin = image_obj.get_channeldata(chx).begin();
//            it_chxdata_end = image_obj.get_channeldata(chx).end();
//            data.insert(data.end(), it_chxdata_begin, it_chxdata_end);

            for(unsigned lauf=0;lauf<image_obj.get_channeldata(chx).size();++lauf)
            {
              data.push_back(image_obj.get_channeldata(chx)[lauf]);
            }
            ++i;
            ++act_line;
          }
        }
        else
          ++i;
      }
      //for accelerated data - add zero-lines to fill missing lines (num_columns/2)
      for(int ii=0; ii<image_obj.get_num_columns()/2-act_line; ++ii)
        data.insert(data.end(),image_obj.get_channeldata(chx).size() , std::complex<float>(0,0));
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

  void ReadSiemensVD11::getRawdata_info(unsigned int& num_rows, unsigned int& num_columns, unsigned int& num_coils,
                                        unsigned short int& acq, unsigned short int& sli, unsigned short int& par, unsigned short int& echo,
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

    //num_rows = image_obj.get_num_columns()/2;   //image_obj.get_ushLine_length();
    num_rows = image_obj.get_ushLine_length();
    num_columns = image_obj.get_num_columns();
    num_coils = image_obj.get_num_channels();
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


  // ---------------------------------------------------------------
  //! int readfile_sizeinfo(unsigned& num_coils, unsigned& nBytes, bool dicom)
  //!
  //! function: reads the num_coils info of an dicom or meas.dat VD11 file
  //! in:
  //! bool dicom      - if file is a dicom-file (different start for header read)
  //! out:
  //!   num_coils     - number of coils
  //!   nBytes        - byte-size per line
  //! return: -1 for error  ,  0 vor ok
  // ---------------------------------------------------------------
  int ReadSiemensVD11::readfile_sizeinfo(unsigned& num_coils, unsigned& nBytes, bool dicom)
  {
    sMDH mdh_data;
    unsigned long cPos;

    if (fileerror_ != 0)
      return -1;

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

    nBytes = getMDHinfo(cPos, mdh_data);

    num_coils = mdh_data.UsedChannels;

    if((num_coils != 0) && (nBytes != 0))
      return 0;
    else
    {
      std::cerr<<"\n reading header: num_coils || nBytes == 0";
      return -1;
    }
  }

}   //namespace agile


