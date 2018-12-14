// config.hpp

#define WITH_TIMER 1

#define LSQR_MAX_ITERATIONS 200
#define LSQR_ABS_TOLERANCE 1e-10
#define SHOW_LSQR_DETAILS 1

#if WITH_TIMER
  #include <iomanip>
  #include <sys/time.h>
#endif