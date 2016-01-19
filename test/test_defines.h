// $Id: test_defines.h 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_TEST_DEFINES_H
#define AGILE_TEST_DEFINES_H

// print debug output for std::vector
//---------------------------------------------------------------
#define PRINT_VEC(name, vec) do {                               \
  std::cout << std::endl << name                                \
            << " (size: " << vec.size() << "):" << std::endl;   \
  for (unsigned i=0; i<vec.size(); ++i)                         \
    std::cout << vec[i] << " ";                                 \
  std::cout << std::endl;                                       \
} while(0);

#define PRINT_VEC_AS_MAT(name, vec, cols, rows, rowMajor) do {  \
  std::cout << std::endl << name                                \
            << " (size: " << vec.size() << "):" << std::endl;   \
  if (rowMajor) {                                               \
    for (unsigned row=0; row<rows; ++row) {                     \
      for (unsigned col=0; col<cols; ++col) {                   \
        std::cout << "\t" << vec[row*cols+col];                 \
      }                                                         \
      std::cout << std::endl;                                   \
    }                                                           \
  } else {                                                      \
    for (unsigned col=0; col<cols; ++col) {                     \
      for (unsigned row=0; row<rows; ++row) {                   \
        std::cout << "\t" << vec[row*cols+col];                 \
      }                                                         \
      std::cout << std::endl;                                   \
    }                                                           \
  }                                                             \
  std::cout << std::endl;                                       \
} while(0);

// print gpu vector (inline copy to host before)
//---------------------------------------------------------------
#define PRINT_GPU_VEC(name, gpu_vec) do {                       \
  std::vector<float> host_vec;                                  \
  gpu_vec.copyToHost(host_vec);                                 \
  PRINT_VEC(name, host_vec)                                     \
} while(0);

// print delimiter
//---------------------------------------------------------------
#define PRINT_DELIMITER() do {                                  \
 for (int i=5; i--;) std::cout << "----------------";           \
} while(0);

// print section heading (for better readability)
//---------------------------------------------------------------
#define PRINT_SECTION(title) do {                               \
  std::cout << std::endl << std::endl << std::endl;             \
  PRINT_DELIMITER();                                            \
  std::cout << std::endl << "| " << title << std::endl;         \
  PRINT_DELIMITER();                                            \
  std::cout << std::endl;                                       \
} while(0);

#endif // AGILE_TEST_DEFINES_H

// End of $Id: test_defines.h 476 2011-06-16 08:54:14Z freiberger $.
