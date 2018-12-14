#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <vector>

#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/oflog/logger.h"

//==============================================================================================
// int gendicom(const char* file_name, const std::vector<float> &data, unsigned int rows, unsigned int cols)
// function: writes the std::vector<float> data with Header Information to a generated Dicom-file
//
// in:      file_name   - dicomfilename
// in:      data        - data to be written
// in:      rows        - number of rows
// in:      cols        - number of columns
//==============================================================================================

int gendicom(const char* file_name, const std::vector<float> &data,
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
if (result.good()) result = dataset->putAndInsertString(DCM_DerivationDescription, "IRGN Processed MR Reconstruction");
if (result.good()) result = dataset->putAndInsertString(DCM_SecondaryCaptureDeviceManufacturer, "IMT TUGRAZ");
if (result.good()) result = dataset->putAndInsertString(DCM_SecondaryCaptureDeviceManufacturerModelName, "IMT Cuda Workstation");
if (result.good()) result = dataset->putAndInsertString(DCM_PatientName, "Doe^John");

// set instance creation date and time
OFString s;
if (result.good()) result = DcmDate::getCurrentDate(s);
if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_InstanceCreationDate, s);
if (result.good()) result = DcmTime::getCurrentTime(s);
if (result.good()) result = dataset->putAndInsertOFStringArray(DCM_InstanceCreationTime, s);


//--- Write image-data ---
std::vector<Uint16> uint16_data;

float val=0;
float min_val;
float max_val = *std::max_element(data.begin(),data.end());
for(unsigned int i=0; i<data.size(); ++i)
{
  val = (data[i]/max_val)*65535;
  uint16_data.push_back(Uint16(val));
}

max_val = *std::max_element(uint16_data.begin(),uint16_data.end());
min_val = *std::min_element(uint16_data.begin(),uint16_data.end());
std::cout<<"\n max-val: "<<max_val;
std::cout<<"\n min-val: "<<min_val;

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

Uint8* pixelData = (Uint8*)&uint16_data[0];
Uint32 pixelLength;

pixelLength = uint16_data.size()*2;   //number of elements in vector * 2bytes (for Uint16)

dataset->putAndInsertUint8Array(DCM_PixelData, pixelData, pixelLength);
OFCondition status = fileformat.saveFile(file_name, EXS_LittleEndianExplicit);

if (result.bad())
  std::cerr << "Error: cannot write DICOM file (" << result.text() << ")" << std::endl;

if (status.bad())
  std::cerr << "Error: cannot write DICOM file (" << status.text() << ")" << std::endl;

return 0;

}
