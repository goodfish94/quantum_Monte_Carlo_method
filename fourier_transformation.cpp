//
//  fourier_transformation.cpp
//  cthyb
//
//  Created by 胡昊昱 on 8/2/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#include "fourier_transformation.hpp"
#include <stdio.h>
#include <cmath>
#include <vector>
std::vector<double> boson_inv_fourier(int N_t, double beta, std::vector<double> &F_w)
{
    using namespace std;
    vector<double> ft(N_t+1,0.0);
    
    double t=0.0;
    double omega_n;
    
    for( int i=0; i< N_t + 1; i++ )
    {
        ft[ i ] += F_w[0];
    }
    
    for( int n=1; n<F_w.size(); n++)
    {
        omega_n = 2.0 *n * M_PI/beta;
        
        for( int i=0; i< N_t + 1; i++ )
        {
            t= i*beta/N_t;
            ft[ i ] += 2.0* F_w[n]*cos( omega_n * t);
            
        }
    }
    for( int i=0; i< N_t + 1; i++ )
    {
        
        ft[ i ]  = ft[i]/beta;
        
    }
    
    return ft;
}

