#pragma once 
#include <math.h>

namespace mag_tools{
    namespace constants{
        static const double pi = std::acos(0.0)*2;
        static const double mu0 = 4e-7*pi;
        static const double nu0 = 1/mu0;
        static const double mu0inv = nu0;
    }
}