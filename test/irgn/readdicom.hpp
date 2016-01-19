#ifndef READDICOM_HPP
#define READDICOM_HPP

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <vector>
#include <string>

#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "dcmtk/oflog/logger.h"

//==============================================================================================
// int readdicom(const char* file_name, const Uint8* pdata, unsigned long& length, std::ofstream& dcmmeas_file)
// function: reads a dicomfile and writes the private MR Raw Data to an ofstream-dcmmeas_file
//
// in:      file_name   - dicomfilename
// in/out:  dcmeas_file - ofstream to which the Raw-data is written
// out:     length      - length of Rawdata (in bytes)
// out:     pdata       - const Uint8 pointer to first data-element (const= value is const)
//==============================================================================================

#define PRIVATE_CREATOR_NAME "MR Siemens RAW Data"

#define PRIVATE_CREATOR_TAG  0x7fe1, 0x0010
#define PRIVATE_ELEMENT1_TAG 0x7fe1, 0x1010

#define PRV_PrivateCreator   DcmTag(PRIVATE_CREATOR_TAG)
#define PRV_PrivateElement1  DcmTag(PRIVATE_ELEMENT1_TAG, PRIVATE_CREATOR_NAME)

int readdicom(const char* file_name, const Uint8* pdata, unsigned long& length, std::ofstream& dcmmeas_file)
{
  int error = 0;
  DcmFileFormat fileformat;
  DcmDataset *dataset = fileformat.getDataset();
  OFCondition result = EC_Normal;

  const Uint8* rawData = NULL;

  //const char* strModality = NULL;
  //Uint16 rows, columns;

  length = 0;

  OFCondition status = fileformat.loadFile(file_name);
  if (status.good())
  {
    if (result.good()) result = dataset->findAndGetUint8Array(PRV_PrivateElement1,rawData,&length);

    //if (result.good()) result = dataset->findAndGetString(DCM_Modality,strModality);
    //if (result.good()) result = dataset->findAndGetUint16(DCM_Rows,rows);
    //if (result.good()) result = dataset->findAndGetUint16(DCM_Columns,columns);

    //size_t printFlags = DCMTypes::PF_shortenLongTagValues;
    //fileformat.print(std::cout, printFlags, 0 /*level*/);
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

#endif // READDICOM_HPP
