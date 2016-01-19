#include "agile/gpu_config.hpp"
#include "agile/gpu_complex.hpp"
#include "agile/gpu_environment.hpp"
#include "agile/gpu_type_traits.hpp"
#include "agile/gpu_vector.hpp"
#include <cuda.h>

// this construct is needed to get the texture in an automated fashion
// for every extract/insert routine you instantiate, it is also needed to
// add such a function here
template <typename TType>
texture<TType>& getGPUCommunicatorTexture();

//template <typename TType>
//static __inline__ __device__ texture<TType> getGPUCommunicatorTextureDevice();

template <typename TType>
struct GPUCommunicatorTextureDevice 
{
  inline __device__ texture<TType> getTexture()
  {
    return NULL;
  };
};

texture<float> extract_insert_texture_float;
template <>
texture<float>& getGPUCommunicatorTexture<float>()
{
  return extract_insert_texture_float;
}
/*template <>
static __inline__ __device__ texture<float>
  getGPUCommunicatorTextureDevice<float>()
{
  return extract_insert_texture_float;
}
*/


texture<float2> extract_insert_texture_float2;
template <>
texture<float2>& getGPUCommunicatorTexture<float2>()
{
  return extract_insert_texture_float2;
}
/*template <>
static __inline__ __device__ texture<float2>
  getGPUCommunicatorTextureDevice<float2>()
{
  return extract_insert_texture_float2;
}
*/

//BEGIN Schwarzl

template <> struct GPUCommunicatorTextureDevice<float> {
  __inline__ __device__ texture<float> getTexture()
  {
    return extract_insert_texture_float;
  };
};

template <> struct GPUCommunicatorTextureDevice<float2> {
  __inline__ __device__ texture<float2> getTexture()
  {
    return extract_insert_texture_float2;
  };
};

//END Schwarzl

//! \brief Extract shared nodes from a vector.
//!
//! This function extracts a set of vector elements and stores them in a
//! new vector \p shared, i.e.
//! \f$ shared[i] \leftarrow x[{map_shared_to_vector}[i]] \f$, where
//! \p x is a texture. The length of \p x has to be larger than
//! \f$ {max}_{0 \le i < N} {map_shared_to_vector}[i] \f$, with
//! \f$ N \f$ being the length of \p map_shared_to_vector.
//! \param[in] map_shared_to_vector Vector defining a map from the shared
//! elements to the big vector.
//! \param[out] shared A vector holding the shared vector elements. This vector
//! has to have the same length as \p map_shared_to_vector.
//! \param[in] num_shared_elements The length of \p map_shared_to_vector and
//! \p shared.
template <typename TType1, typename TType2>
__global__ void extractSharedNodes_GPU(const TType1* map_shared_to_vector,
                                       TType2* shared,
                                       unsigned num_shared_elements)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  if (thread_id < num_shared_elements)
  {
    TType2 vector_value = tex1Dfetch(
      GPUCommunicatorTextureDevice<typename agile::to_tuple_type<TType2>::type>().getTexture(),
      map_shared_to_vector[thread_id]);
    shared[thread_id] = vector_value;
  }
}

//! \brief Insert shared nodes back in a large vector.
//!
//! This function takes the shared nodes from a texture and inserts them
//! back into a large vector, i.e.
//! \f$ x[i] \leftarrow y[{map_vector_to_shared}[i]] \f$.
//! \param[in] map_vector_to_shared A mapping from the vector elements to the
//! shared elements. If a vector element is not shared, use a destination index
//! which is larger than \p num_shared_indices.
//! \param[in] map_length The length of the vector \p map_vector_to_shared.
//! \param[in] num_shared_indices The amount of shared indices.
//! \param[out] vector The large vector holding all elements. This vector has to
//! be at least as long as \p map_vector_to_shared.
template <typename TType1, typename TType2>
__global__ void insertSharedNodes_GPU(const TType1* map_vector_to_shared,
                                      unsigned map_length,
                                      unsigned num_shared_indices,
                                      TType2* vector)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  if (thread_id < map_length)
  {
    // get the index of the shared node
    unsigned shared_index = map_vector_to_shared[thread_id];
    // not every vector element has a corresponding shared element so make sure
    // the index is valid
    if (shared_index < num_shared_indices)
    {
      TType2 shared_value =tex1Dfetch(
        GPUCommunicatorTextureDevice<typename agile::to_tuple_type<TType2>::type>().getTexture(),
        shared_index);
      vector[thread_id] = shared_value;
    }
  }
}

//! \brief Multiply x and y element-wise and store result in y (GPU function).
template <typename TType1, typename TType2>
__global__ void multiplySharedNodes_GPU(
  const TType1* x, TType2* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  if (thread_id < size)
    y[thread_id] *= x[thread_id];
}

namespace agile
{
  //! \brief Extract vector elements (host function).
  //!
  //! Extracts the shared nodes of vector \p v and stores them in \p shared.
  //! \param[in] map_shared_to_vector This is a vector containing the indices
  //! of the shared nodes.
  //! \param[in] v The vector to extract the shared elements from. The vector's
  //! length has to be larger than the largest value in \p m_shared_to_vector.
  //! \param[out] shared A pointer to the receive buffer in GPU memory. The
  //! byte-size of the pointer is the length of \p map_shared_to_vector times
  //! the byte-size of \p TType2.
  template <typename TType1, typename TType2>
  void extractSharedNodes(const GPUVector<TType1>& map_shared_to_vector,
                          const GPUVector<TType2>& v, void* shared)
  {
    // bind a texture to the vector as we need random access to it
    cudaBindTexture(
      0, getGPUCommunicatorTexture<
           typename to_tuple_type<TType2>::type>(), v.data(),
      v.size() * sizeof(typename substitute_gpu_complex<TType2>::type));

    unsigned grid_size = (map_shared_to_vector.size()
                          + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
                         / GPUEnvironment::getMaxNumThreadsPerBlock();

    extractSharedNodes_GPU<<<grid_size,
                             GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)
        map_shared_to_vector.data(),
      (typename substitute_gpu_complex<TType2>::type*)shared,
      map_shared_to_vector.size());

    // free the texture again
    cudaUnbindTexture(
      getGPUCommunicatorTexture<typename to_tuple_type<TType2>::type>());
  }

  //! \brief Insert vector elements (host function).
  //!
  //! Inserts the shared nodes stored in the buffer \p shared back into the
  //! vector \p v.
  //! \param[in] map_vector_to_shared This vector holds the mapping from vector
  //! indices to shared indices. If a vector element \p i is not to be shared,
  //! the \p i-th element in this vector has to be set to a value larger than
  //! \p num_shared_indices.
  //! \param[in] shared Pointer to GPU memory holding the shared elements.
  //! \param[in] num_shared_indices The amount of indices that are shared.
  //! \param[in,out] v The vector with the shared elements inserted.
  template <typename TType1, typename TType2>
  void insertSharedNodes(const GPUVector<TType1>& map_vector_to_shared,
                         const void* shared, unsigned num_shared_indices,
                         GPUVector<TType2>& v)
  {
    // bind a texture to the vector as we need random access to it
    cudaBindTexture(
      0, getGPUCommunicatorTexture<
           typename to_tuple_type<TType2>::type>(), shared,
      num_shared_indices
        * sizeof(typename substitute_gpu_complex<TType2>::type));

    unsigned grid_size = (map_vector_to_shared.size()
                          + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
                         / GPUEnvironment::getMaxNumThreadsPerBlock();

    insertSharedNodes_GPU<<<grid_size,
                            GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)
        map_vector_to_shared.data(),
      map_vector_to_shared.size(), num_shared_indices,
      (typename substitute_gpu_complex<TType2>::type*)v.data());

    // free the texture again
    cudaUnbindTexture(
      getGPUCommunicatorTexture<typename to_tuple_type<TType2>::type>());
  }

  //! \brief Multiply shared nodes with a factor (host function).
  template <typename TType1, typename TType2, typename TType3>
  void multiplySharedNodes(
    const GPUVector<TType1>& shared_node_indices,
    const GPUVector<TType1>& vector_to_shared_node_indices,
    const GPUVector<TType2>& shared_node_multiplicand,
    void* gpu_buffer, GPUVector<TType3>& v)
  {
    // first extract the shared nodes
    extractSharedNodes(shared_node_indices, v, gpu_buffer);
    // multiply
    unsigned grid_size = (shared_node_indices.size()
                          + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
                         / GPUEnvironment::getMaxNumThreadsPerBlock();
    multiplySharedNodes_GPU<<<grid_size,
                              GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType2>::type*)
        shared_node_multiplicand.data(),
      (typename substitute_gpu_complex<TType3>::type*)gpu_buffer,
      shared_node_multiplicand.size());

    // insert the elements again
    insertSharedNodes(vector_to_shared_node_indices, gpu_buffer,
      shared_node_indices.size(), v);
  }

  // explicit instantiation
  template void extractSharedNodes<unsigned, float>(
    const GPUVector<unsigned>& map_shared_to_vector,
    const GPUVector<float>& v, void* shared);
  template void extractSharedNodes<unsigned, std::complex<float> >(
    const GPUVector<unsigned>& map_shared_to_vector,
    const GPUVector<std::complex<float> >& v, void* shared);

  template void insertSharedNodes<unsigned, float>(
    const GPUVector<unsigned>& map_vector_to_shared,
    const void* shared, unsigned num_shared_indices,
    GPUVector<float>& v);
  template void insertSharedNodes<unsigned, std::complex<float> >(
    const GPUVector<unsigned>& map_vector_to_shared,
    const void* shared, unsigned num_shared_indices,
    GPUVector<std::complex<float> >& v);

  template void multiplySharedNodes<unsigned, float, float>(
    const GPUVector<unsigned>& shared_node_indices,
    const GPUVector<unsigned>& vector_to_shared_node_indices,
    const GPUVector<float>& shared_node_multiplicand,
    void* gpu_buffer, GPUVector<float>& v);
  template void multiplySharedNodes<unsigned, float, std::complex<float> >(
    const GPUVector<unsigned>& shared_node_indices,
    const GPUVector<unsigned>& vector_to_shared_node_indices,
    const GPUVector<float>& shared_node_multiplicand,
    void* gpu_buffer, GPUVector<std::complex<float> >& v);

} // namespace agile
