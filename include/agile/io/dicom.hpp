#ifndef DICOM_HPP
#define DICOM_HPP

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <vector>
#include <string>
#include <algorithm>

#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "dcmtk/oflog/logger.h"


namespace agile
{
  #define PRIVATE_CREATOR_NAME_RAWData "MR Siemens RAW Data"

  #define PRIVATE_CREATOR_TAG_RAWData  0x7fe1, 0x0010
  #define PRIVATE_ELEMENT1_TAG_RAWData 0x7fe1, 0x1010

  #define PRV_PrivateCreator_RAWData   DcmTag(PRIVATE_CREATOR_TAG_RAWData)
  #define PRV_PrivateElement1_RAWData  DcmTag(PRIVATE_ELEMENT1_TAG_RAWData, PRIVATE_CREATOR_NAME_RAWData)


  #define PRIVATE_CREATOR_NAME_FOV "Siemens MR HEADER"

  #define PRIVATE_CREATOR_TAG_FOV  0x0051, 0x0010
  #define PRIVATE_ELEMENT1_TAG_FOV 0x0051, 0x100C

  #define PRV_PrivateCreator_FOV   DcmTag(PRIVATE_CREATOR_TAG_FOV)
  #define PRV_PrivateElement1_FOV  DcmTag(PRIVATE_ELEMENT1_TAG_FOV, PRIVATE_CREATOR_NAME_FOV)


  struct DicomInfo{
    std::string str_studyinstuid;
    std::string str_seriesinstuid;
    std::string str_seriesdescription;
    std::string str_patientname;
    std::string str_repetitiontime;
    std::string str_echotime;
    std::string str_flipangle;
    std::string str_nbofphaseencodingsteps;
    std::string str_slicethickness;
    std::string str_fieldofview;
    unsigned short ush_rows;
    unsigned short ush_columns;
  };


  class DICOM
  {
    public:

      //! \brief default Constructor.
      //!
      //! The default constructor creates an empty DICOM.
      DICOM(){
        _str_studyinstuid = " ";
        _str_seriesinstuid = " ";
        _str_seriesdescription = " ";
        _str_patientname = " ";
        _str_repetitiontime = " ";
        _str_echotime = " ";
        _str_flipangle = " ";
        _str_nbofphaseencodingsteps = " ";
        _str_slicethickness = " ";
        _str_fieldofview = " ";
      }

      //! \brief Destructor.
      //!
      //! The Destructor
      ~DICOM(){
      }

      int readdicom(std::string file_name, std::ofstream& dcmmeas_file, unsigned long& length, const Uint8* pdata = NULL);

      int readdicom_info(std::string file_name);


      int gendicom(const char* file_name, const std::vector<float> &data,
                   unsigned int rows, unsigned int cols);

      int gendicom_info(const char* file_name, DicomInfo dcminfo);

      void set_filename(std::string filename)
      {
        this->_file_name=filename;
      }

      std::vector<float> get_data()
      {
        return _data;
      }

      //returns some header information
      DicomInfo get_dicominfo()
      {
        DicomInfo dcminfo;
        dcminfo.str_studyinstuid          = _str_studyinstuid;
        dcminfo.str_seriesinstuid         = _str_seriesinstuid;
        dcminfo.str_seriesdescription     = _str_seriesdescription;
        dcminfo.str_patientname           = _str_patientname;
        dcminfo.str_repetitiontime        = _str_repetitiontime;
        dcminfo.str_echotime              = _str_echotime;
        dcminfo.str_flipangle             = _str_flipangle;
        dcminfo.str_nbofphaseencodingsteps= _str_nbofphaseencodingsteps;
        dcminfo.str_slicethickness        = _str_slicethickness;
        dcminfo.str_fieldofview           = _str_fieldofview;
        dcminfo.ush_rows = _ush_rows;
        dcminfo.ush_columns = _ush_columns;

        return dcminfo;
      }

      void set_dicominfo(DicomInfo dcminfo)
      {
        _str_studyinstuid = dcminfo.str_studyinstuid;
        _str_seriesinstuid = dcminfo.str_seriesinstuid;
        _str_seriesdescription = dcminfo.str_seriesdescription;
        _str_patientname = dcminfo.str_patientname;
        _str_repetitiontime = dcminfo.str_repetitiontime;
        _str_echotime = dcminfo.str_echotime;
        _str_flipangle = dcminfo.str_flipangle;
        _str_nbofphaseencodingsteps = dcminfo.str_nbofphaseencodingsteps;
        _str_slicethickness = dcminfo.str_slicethickness;
        _str_fieldofview = dcminfo.str_fieldofview;

      }

      DcmDataset* getDicomDataSet() 
      {
        return _dataset;
      }


    private:
      std::string _file_name;

      std::ofstream _dcmmeas_file;

      std::string _str_studyinstuid;
      std::string _str_seriesinstuid;
      std::string _str_seriesdescription;
      std::string _str_patientname;
      std::string _str_repetitiontime;
      std::string _str_echotime;
      std::string _str_flipangle;
      std::string _str_nbofphaseencodingsteps;
      std::string _str_slicethickness;
      std::string _str_fieldofview;
      unsigned short _ush_rows;
      unsigned short _ush_columns;

      std::vector<float> _data;
      unsigned int rows;
      unsigned int cols;

      DcmFileFormat _dcmfile;
      DcmDataset * _dataset;
  };

}   //namespace agile
#endif // DICOM_HPP
