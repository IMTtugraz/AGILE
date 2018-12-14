
#include "agile/gpu_environment.hpp"
#include "agile/network/gpu_communicator.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/operator/binary_measure_operator.hpp"
#include <iostream>
#include <time.h>

#define TEST_VECTOR_SIZE 5



inline void seperator()
{
  std::cout << "------------------------------------------------------------"
            << std::endl;
}

float getRandomNumber()
{
  float num = rand() % 10;
  num += (rand() % 1000) / 1000.0f;
  return num;
}

template <typename TType>
void initWithDemoData(std::vector<TType> &vec)
{
  for (unsigned i=0; i<vec.size(); ++i)
    vec[i] = TType(getRandomNumber());
}
template <typename TType>
void initWithDemoData(std::vector<std::complex<TType> > &vec)
{
  for (unsigned i=0; i<vec.size(); ++i)
    vec[i] = std::complex<TType>(getRandomNumber(), getRandomNumber());
}

template <typename TType>
TType getCPUScalarProduct(std::vector<TType> &vec_x, std::vector<TType> &vec_y)
{
  TType result = 0;
  for (unsigned i=0; i<vec_x.size(); ++i)
    result += vec_x[i] * vec_y[i];
  return result;
}

template <typename TType>
std::complex<TType> getCPUScalarProduct(
    std::vector<std::complex<TType> > &vec_x,
    std::vector<std::complex<TType> > &vec_y)
{
  std::complex<TType> result(0,0);
  for (unsigned i=0; i<vec_x.size(); ++i)
    result += conj(vec_x[i]) * vec_y[i];
  return result;
}

template <typename TCommunicatorType, typename TDataType>
void testScalarProduct(TCommunicatorType &com, const unsigned testVectorSize)
{
  seperator();
    
  // build testdata
  std::vector<TDataType> vec_x(testVectorSize), vec_y(testVectorSize);
  
  // init test vectors
  initWithDemoData(vec_x);
  initWithDemoData(vec_y);

  agile::GPUVector<TDataType> gpu_vec_x, gpu_vec_y;
  gpu_vec_x.assignFromHost(vec_x.begin(), vec_x.end());
  gpu_vec_y.assignFromHost(vec_y.begin(), vec_y.end());
  
    // generate a binary measure
  typedef agile::ScalarProductMeasure<TCommunicatorType> measure_type;
  measure_type scalar_product(com);
  
  // calculate dot product
  TDataType dot_product = scalar_product(gpu_vec_x, gpu_vec_y);
  
  std::cout << std::endl << "x: ";
  for (unsigned i=0; i<vec_x.size(); ++i)
  {
    std::cout << vec_x[i] << " ";
  }
  std::cout << std::endl << "y: ";
  for (unsigned i=0; i<vec_y.size(); ++i)
  {
    std::cout << vec_y[i] << " ";
  }
  
  // Test solution
  std::cout << std::endl
            << "dot product GPU: " << dot_product
            << std::endl;

  // Reference solution
  std::cout << "dot product CPU: " << getCPUScalarProduct(vec_x, vec_y)
            << std::endl << std::endl;            
}


int main(int argc, char* argv[])
{
  // init random number generator
  srand ( time(NULL) );
  
  // init the network
  agile::NetworkEnvironment environment(argc, argv);
  
  // allocate a GPU
  typedef agile::GPUCommunicator<unsigned, float, float> communicator_type;
  communicator_type com;
  com.allocateGPU();

  // GPU Information
  agile::GPUEnvironment::printInformation(std::cout);
  std::cout << std::endl;
  
  // Tests
  testScalarProduct<communicator_type, float>(com, TEST_VECTOR_SIZE);
  testScalarProduct<communicator_type, std::complex<float> >(com, TEST_VECTOR_SIZE);

  return 0;
}
