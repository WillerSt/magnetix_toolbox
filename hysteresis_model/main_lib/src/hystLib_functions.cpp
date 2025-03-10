// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#include "hystLib_functions.h"

    double hystLib::functions::acos_x(const double& alpha)
    {
        if (alpha > 1)
        {
            return 0.0;
        }
        else if(alpha < -1)
        {
            return hystLib::constants::pi;
        }
        else
        {
            return std::acos(alpha);
        }
    }
            
        