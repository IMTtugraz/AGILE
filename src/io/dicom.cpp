#include "agile/io/dicom.hpp"

//==============================================================================================
// int readdicom(const char* file_name, const Uint8* pdata, unsigned long& length, std::ofstream& dcmmeas_file)
// function: reads a dicomfile and writes the private MR Raw Data to an ofstream-dcmmeas_file
//
// in:      file_name   - dicomfilename
// in/out:  dcmeas_file - ofstream to which the Raw-data is written
// out:     length      - length of Rawdata (in bytes)
// out:     pdata       - const Uint8 pointer to first data-element (const= value is const)
//==============================================================================================

namespace agile
{
  int DICOM::readdicom(std::string file_name, std::ofstream& dcmmeas_file, unsigned long& length, const Uint8* pdata)
  {
    int error = 0;
    DcmDataset *dataset = _dcmfile.getDataset();
    OFCondition result = EC_Normal;

    const Uint8* rawData = NULL;

    //const char* strModality = NULL;
    //Uint16 rows, columns;

    length = 0;

    OFCondition status = _dcmfile.loadFile(file_name.c_str());
    if (status.good())
    {
      if (result.good()) result = dataset->findAndGetUint8Array(PRV_PrivateElement1_RAWData,rawData,&length);

      //if (result.good()) result = dataset->findAndGetString(DCM_Modality,strModality);
      //if (result.good()) result = dataset->findAndGetUint16(DCM_Rows,rows);
      //if (result.good()) result = dataset->findAndGetUint16(DCM_Columns,columns);

      //size_t printFlags = DCMTypes::PF_shortenLongTagValues;
      //_dcmfile.print(std::cout, printFlags, 0 /*level*/);
    }

    pdata=rawData;

    dcmmeas_file.write((char*) pdata, length);

    if(dcmmeas_file.bad())
    {
      std::cerr << "Error: cannot write data raw file - readdicom"<< std::endl;
      error = -1;
    }

    if (result.bad())
    {
      std::cerr << "Error: cannot read DICOM file Tag: (" << result.text() << ")" << std::endl;
      error = -1;
    }

    if (status.bad())
    {
      std::cerr << "Error: cannot read DICOM file (" << status.text() << ")" << std::endl;
      error = -1;
    }

  return error;
  }

  //==============================================================================================
  // int readdicom_info(const char* file_name, std::string& str_studyinstuid, std::string& str_seriesinstuid)
  // function: reads a dicomfile and writes the private MR Raw Data to an ofstream-dcmmeas_file
  //
  // in:      file_name   - dicomfilename
  // out:     str_studyinstuid  - DCM_StudyInstanceUID
  //          str_seriesinstuid - DCM_SeriesInstanceUID
  //          0010,0010  Patient's Name: TUG_Rawdata_Headcoil
  //          0018,0080  Repetition Time: 250
  //          0018,0081  Echo Time: 2.48
  //          0008,103E  Series Description: t1_fl2d_tra
  //          0018,1314  Flip Angle: 70
  //          0018,0089  Number of Phase Encoding Steps: 256
  //          0018,0050  Slice Thickness: 4
  //          private 0051,100C  ---: FoV 200*200
  //==============================================================================================
  int DICOM::readdicom_info(std::string file_name)
  {
    int error = 0;
    _dataset = _dcmfile.getDataset();
    OFCondition result = EC_Normal;
    OFString studyinstuid="", seriesinstuid="";
    OFString seriesdescription="", patientname="", repetitiontime="", echotime="", flipangle="", nbofphaseencodingsteps="",
        slicethickness="", fieldofview="";

    OFCondition status = _dcmfile.loadFile(file_name.c_str());
    if (status.good())
    {
      _dataset->findAndGetOFString(DCM_StudyInstanceUID,studyinstuid);
      _dataset->findAndGetOFString(DCM_SeriesInstanceUID,seriesinstuid);
      _dataset->findAndGetOFString(DCM_PatientName,patientname);
      _dataset->findAndGetOFString(DCM_RepetitionTime,repetitiontime);
      _dataset->findAndGetOFString(DCM_SeriesDescription,seriesdescription);
      _dataset->findAndGetOFString(DCM_EchoTime,echotime);
      _dataset->findAndGetOFString(DCM_FlipAngle,flipangle);
      _dataset->findAndGetOFString(DCM_NumberOfPhaseEncodingSteps,nbofphaseencodingsteps);
      _dataset->findAndGetOFString(DCM_SliceThickness,slicethickness);
      _dataset->findAndGetOFString(PRV_PrivateElement1_FOV,fieldofview);
      _dataset->findAndGetUint16(DCM_Rows, _ush_rows);
      _dataset->findAndGetUint16(DCM_Columns, _ush_columns);
    }

    _str_studyinstuid = studyinstuid.c_str();
    _str_seriesinstuid = seriesinstuid.c_str();
    _str_seriesdescription = seriesdescription.c_str();
    _str_patientname = patientname.c_str();
    _str_repetitiontime = repetitiontime.c_str();
    _str_echotime = echotime.c_str();
    _str_flipangle = flipangle.c_str();
    _str_nbofphaseencodingsteps = nbofphaseencodingsteps.c_str();
    _str_slicethickness = slicethickness.c_str();
    _str_fieldofview = fieldofview.c_str();


    if (result.bad())
    {
      std::cerr << "Error: cannot read DICOM file Tag: (" << result.text() << ")" << std::endl;
      error = -1;
    }

    if (status.bad())
    {
      std::cerr << "Error: cannot load DICOM file (" << status.text() << ")" << std::endl;
      error = -1;
    }

    return error;
  }

  //==============================================================================================
  // int gendicom(const char* file_name, const std::vector<float> &data, unsigned int rows, unsigned int cols)
  // function: writes the std::vector<float> data with Header Information to a generated Dicom-file
  //
  // in:      file_name   - dicomfilename
  // in:      data        - data to be written
  // in:      rows        - number of rows
  // in:      cols        - number of columns
  //==============================================================================================

  int DICOM::gendicom(const char* file_name, const std::vector<float> &data,
               unsigned int rows, unsigned int cols)
  {

    char uid[100];
    DcmFileFormat fileformat;
    DcmDataset *dataset = fileformat.getDataset();
    OFCondition result = EC_Normal;

    if (result.good()) result = dataset->putAndInsertString(DCM_SOPClassUID, UID_SecondaryCaptureImageStorage);

    //if (result.good()) result = dataset->putAndInsertString(DCM_MediaStorageSOPClassUID, UID_SecondaryCaptureImageStorage); //same as SOP
    if (result.good()) result = dataset->putAndInsertString(DCM_SOPInstanceUID, dcmGenerateUniqueIdentifier(uid, SITE_INSTANCE_UID_ROOT));
    //if (result.good()) result = dataset->putAndInsertString(DCM_MediaStorageSOPInstanceUID, dcmGenerateUniqueIdentifier(uid, SITE_INSTANCE_UID_ROOT));
    if (result.good()) result = dataset->putAndInsertString(DCM_SeriesInstanceUID, dcmGenerateUniqueIdentifier(uid, SITE_INSTANCE_UID_ROOT));
    if (result.good()) result = dataset->putAndInsertString(DCM_StudyInstanceUID, dcmGenerateUniqueIdentifier(uid, SITE_INSTANCE_UID_ROOT));
    //if (result.good()) result = dataset->putAndInsertString(DCM_TransferSyntaxUID, UID_LittleEndianExplicitTransferSyntax);

    //if (result.good()) result = dataset->putAndInsertUint16(DCM_FileMetaInformationVersion, 1);

    if (result.good()) result = dataset->putAndInsertString(DCM_ImageType, "DERIVED");
    if (result.good()) result = dataset->putAndInsertString(DCM_Modality, "MR");
    if (result.good()) result = dataset->putAndInsertString(DCM_ConversionType, "WSD");
    if (result.good()) result = dataset->putAndInsertString(DCM_DerivationDescription, "IMT Processed MR Reconstruction");
    if (result.good()) result = dataset->putAndInsertString(DCM_SecondaryCaptureDeviceManufacturer, "IMT TUGRAZ");
    if (result.good()) result = dataset->putAndInsertString(DCM_SecondaryCaptureDeviceManufacturerModelName, "IMT Workstation");
    if (result.good()) result = dataset->putAndInsertString(DCM_PatientName, "Doe^John");

    // set instance creation date and time
    OFString s;
    if (result.good()) result = DcmDate::getCurrentDate(s);
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_InstanceCreationDate, s);
    if (result.good()) result = DcmTime::getCurrentTime(s);
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_InstanceCreationTime, s);


    //--- Write image-data ---
    std::vector<Uint16> uint16_data;
    uint16_data.clear();

    float val=0;
    float max_val_float = *std::max_element(data.begin(),data.end());
    Uint16 min_val;
    Uint16 max_val;


    //std::cout<<"\n max: "<<*std::max_element(data.begin(),data.end());
    //std::cout<<"\n min: "<<*std::min_element(data.begin(),data.end());

    //for(int x=data.size()/2;x<data.size()/2+10;++x)
    //  std::cout<<"\n image data-gendicom: "<<data[x];

    for(unsigned int i=0; i<data.size(); ++i)
    {
      val = (data[i]/max_val_float)*float(65535);
      uint16_data.push_back(Uint16(val));
    }

    max_val = *std::max_element(uint16_data.begin(),uint16_data.end());
    min_val = *std::min_element(uint16_data.begin(),uint16_data.end());

    //std::cout<<"\n uint16_data.size(): "<<uint16_data.size();
    //std::cout<<"\n max_val: "<<max_val;
    //std::cout<<"\n min_val: "<<min_val;

    unsigned bits=16;
    Uint16 bitsAllocated=((bits-1)/8+1)*8;
    Uint16 bitsStored=bits;
    Uint16 highBit=bits-1;
    if (result.good()) result = dataset->putAndInsertUint16(DCM_BitsAllocated, bitsAllocated);
    if (result.good()) result = dataset->putAndInsertUint16(DCM_BitsStored, bitsStored);
    if (result.good()) result = dataset->putAndInsertUint16(DCM_HighBit, highBit);

    if (result.good()) result = dataset->putAndInsertUint16(DCM_Rows, rows);
    if (result.good()) result = dataset->putAndInsertUint16(DCM_Columns, cols);

    if (result.good()) result = dataset->putAndInsertUint16(DCM_PixelRepresentation, 0);  // 1 signed, 0 unsigned
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_PhotometricInterpretation, "MONOCHROME2");
    if (result.good()) result = dataset->putAndInsertUint16(DCM_SamplesPerPixel, 1);

    if (result.good()) result = dataset->putAndInsertUint16(DCM_SmallestImagePixelValue, min_val);
    if (result.good()) result = dataset->putAndInsertUint16(DCM_LargestImagePixelValue, max_val);

    std::string patientname = _str_patientname;
    std::string repetitiontime = _str_repetitiontime;
    std::string echotime = _str_echotime;
    std::string seriesdescription = _str_seriesdescription;
    std::string flipangle = _str_flipangle;
    std::string nbofphaseencodingsteps = _str_nbofphaseencodingsteps;
    std::string slicethickness = _str_slicethickness;
    std::string fieldofview = _str_fieldofview;

    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_FieldOfViewDimensions, OFString(&fieldofview[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_SliceThickness, OFString(&slicethickness[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_NumberOfPhaseEncodingSteps, OFString(&nbofphaseencodingsteps[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_FlipAngle, OFString(&flipangle[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_SeriesDescription, OFString(&seriesdescription[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_EchoTime, OFString(&echotime[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_RepetitionTime, OFString(&repetitiontime[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_PatientName, OFString(&patientname[0]));
    //if (result.good()) result = dataset->putAndInsertOFStringArray(PRV_PrivateElement1_FOV, OFString(&fieldofview[0]));


    Uint8* pixelData = (Uint8*)&uint16_data[0];
    Uint32 pixelLength;

    pixelLength = uint16_data.size()*2;   //number of elements in vector * 2bytes (for Uint16)

    dataset->putAndInsertUint8Array(DCM_PixelData, pixelData, pixelLength);
    OFCondition status = fileformat.saveFile(file_name, EXS_LittleEndianExplicit);

    if (result.bad())
      std::cerr << "Error: cannot write DICOM file (" << result.text() << ")" << std::endl;

    if (status.bad())
      std::cerr << "Error: cannot save DICOM file (" << status.text() << ")" << std::endl;

    return 0;

  }

  //==============================================================================================
  // int gendicom_info(const char* file_name, std::string patientname, std::string repetitiontime,
  // std::string echotime, std::string seriesdescription, std::string flipangle,
  // std::string nbofphaseencodingsteps, std::string slicethickness, std::string str_fieldofview)
  //
  // function: writes  Header Information to an existing Dicom-file
  //
  // in:      file_name   - dicomfilename
  // out:     str_studyinstuid  - DCM_StudyInstanceUID
  //          str_seriesinstuid - DCM_SeriesInstanceUID
  //          0010,0010  Patient's Name: TUG_Rawdata_Headcoil
  //          0018,0080  Repetition Time: 250
  //          0018,0081  Echo Time: 2.48
  //          0008,103E  Series Description: t1_fl2d_tra
  //          0018,1314  Flip Angle: 70
  //          0018,0089  Number of Phase Encoding Steps: 256
  //          0018,0050  Slice Thickness: 4
  //          private 0051,100C  ---: FoV 200*200
  //==============================================================================================
  int DICOM::gendicom_info(const char* file_name, DicomInfo dcminfo)
  {

    char uid[100];
    DcmFileFormat fileformat;
    OFCondition result = EC_Normal;
    OFCondition readstat = fileformat.loadFile(file_name);

    DcmDataset *dataset = fileformat.getDataset();

    std::string patientname = dcminfo.str_patientname;
    std::string repetitiontime = dcminfo.str_repetitiontime;
    std::string echotime = dcminfo.str_echotime;
    std::string seriesdescription = dcminfo.str_seriesdescription;
    std::string flipangle = dcminfo.str_flipangle;
    std::string nbofphaseencodingsteps = dcminfo.str_nbofphaseencodingsteps;
    std::string slicethickness = dcminfo.str_slicethickness;
    std::string fieldofview = dcminfo.str_fieldofview;

    //std::cout<<"\n patientname: "<<patientname;

    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_SliceThickness, OFString(&slicethickness[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_NumberOfPhaseEncodingSteps, OFString(&nbofphaseencodingsteps[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_FlipAngle, OFString(&flipangle[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_SeriesDescription, OFString(&seriesdescription[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_EchoTime, OFString(&echotime[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_RepetitionTime, OFString(&repetitiontime[0]));
    if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_PatientName, OFString(&patientname[0]));
    //if (result.good()) result = dataset->putAndInsertOFStringArray(PRV_PrivateElement1_FOV, OFString(&fieldofview[0]));

    OFCondition status = fileformat.saveFile(file_name, EXS_LittleEndianExplicit);

    if (readstat.bad())
      std::cerr << "Error: cannot read DICOM file (" << readstat.text() << ")" << std::endl;

    if (result.bad())
      std::cerr << "Error: cannot write DICOM file (" << result.text() << ")" << std::endl;

    if (status.bad())
      std::cerr << "Error: cannot save DICOM file (" << status.text() << ")" << std::endl;

    return 0;

    }

}   //namespace agile
