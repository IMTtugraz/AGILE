#ifndef AGILE_MATRIXHELPER_HPP
#define AGILE_MATRIXHELPER_HPP

#include <fstream>
#include <limits>
#include <iostream>
#include <vector>

// Function to print a vector to \p std::cout.
template <typename TType>
void output(const char* string, const std::vector<TType>& x)
{
  std::cout << string;
  for (unsigned counter = 0; counter < x.size(); ++counter)
    std::cout << x[counter] <<std::setprecision(20) << " ";
  std::cout << std::endl;
};


//Function to print a matrix to \p std::cout.
template <typename TType>
void output(const char* string, unsigned num_rows, unsigned num_columns,
            const std::vector<TType>& data)
{
  typename std::vector<TType>::const_iterator iter = data.begin();
  std::cout << string << std::endl;
  for (unsigned row = 0; row < num_rows; ++row)
  {
    std::cout << "  ";
    for (unsigned column = 0; column < num_columns; ++column)
      std::cout <<std::setprecision(3) << std::setw(10) << *iter++;
    std::cout << std::endl;
  }
  std::cout << std::endl;
}


//-------------------------------------------------------------
//  INIT Matrix-Logger
//
//  call example:
//             std::fstream myfile;
//             init_matrixlog(myfile,"matrixlog.txt");
//-------------------------------------------------------------
void init_matrixlog(std::fstream &myfile, const char * file_name)
{

  //file_name = "matrixlog.txt";
  myfile.open (file_name, std::fstream::out);

  if (!myfile.is_open())
  {
      std::cerr << "File not found: " << file_name << std::endl;
  }

}

//-------------------------------------------------------------
//  Close Matrix-Logger
//
//  call example:
//            close_matrixlog(myfile);
//-------------------------------------------------------------
void close_matrixlog(std::fstream &myfile)
{

  myfile.close();

}

//-------------------------------------------------------------
//  Matrix-Logger for Vector/Matrix
//
//  call example:
//            matrixlog(myfile,"Matrix 1 ",num_rows,num_columns, data_vector);
//-------------------------------------------------------------
template <typename TType>
void matrixlog(std::fstream &myfile, const char* string, unsigned num_rows, unsigned num_columns,
            const std::vector<TType>& data)
{
  //write to file:
  typename std::vector<TType>::const_iterator iter = data.begin();

  myfile << std::endl<< std::endl<<string;
  myfile << std::endl << "num_rows: "<<num_rows << "\t" << "num_columns: "<<num_columns<<std::endl;
  myfile << std::endl;
  for (unsigned row = 0; row < num_rows; ++row)
  {
    myfile << "  ";
    for (unsigned column = 0; column < num_columns; ++column)
      myfile <<std::setprecision(10) << std::setw(17) << *iter++;
    myfile << std::endl;
  }
  myfile << std::endl;

}

//-------------------------------------------------------------
//  Matrix-Logger for single value
//
//  call example:
//            matrixlog(myfile,"Matrix 1 ", value);
//-------------------------------------------------------------
template <typename TType>
void matrixlog(std::fstream &myfile, const char* string, const TType data)
{

  //write to file:

  myfile << std::endl<< std::endl<<string;
  myfile << std::endl;
  myfile << "  ";
  myfile <<std::setprecision(3) << std::setw(10) << data;
  myfile << std::endl;
  myfile << std::endl;

}

//-------------------------------------------------------------
//  Matrix-Logger for string
//
//  call example:
//            matrixlog(myfile,"Matrix 1 ");
//-------------------------------------------------------------
template <typename TType>
void matrixlog(std::fstream &myfile, const char* string)
{

  //write to file:

  myfile <<string;
  myfile << std::endl;

}


#endif


/*

 std::vector<TType> data_host;
 Coil[0].copyToHost(data_host);
 std::cout<<"\n data: ";
 output("  ",10,10, data_host);


 std::vector<TType> data_host;
 Coil[0].copyToHost(data_host);
 matrixlog(myfile,"Matrix 1 ",num_rows,num_columns, data_host);

 matrixlog(myfile,"Value 1 ", value);

*/
