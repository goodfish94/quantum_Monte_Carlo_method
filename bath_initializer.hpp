//
//  bath_initializer.hpp
//  cthyb
//
//  Created by 胡昊昱 on 7/28/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#ifndef bath_initializer_hpp
#define bath_initializer_hpp

#include <stdio.h>
#include <cmath>
#include <vector>
using namespace std;

vector<double> fermion_bath_flat(int N_t, double D, double v, double beta)
{// rho(e) = 1/(2D) theta( 2D-|e| )
    int N_e = 10*N_t;
    vector< double > G_(N_t+1,0.0);
    
    double nf = 0.0;
    double de,t,ep;
    
    de = 2*D/double(N_e);
    for( int i=0;i<N_e;i++)
    {
        ep = -D + i*de;
        nf = 1.0/(exp(beta*ep)+1);
        
        for( int it =0 ; it<N_t+1;it++)
        {
            t= it/double(N_t) * beta;
            G_[it] += nf/(2.0*D) * exp((beta-t)*ep) * de;
        }
    }
    for( int it=0; it<N_t+1; it++)
    {
        G_[it] = -G_[it]*v*v;
    }
    return G_;
    
}

vector<double> boson_bath(int N_t, double s, double lambda, double coeff, double beta)
{// rho(w) = k_0 w^s e^{-w/lambda} theta(w) theta( 10*lambda- w)
    // bath = bath * coeff
    int N_e = 50*N_t;
    int N_w = 10*N_t;
    vector< double > J_(N_t+1,0.0);
    vector< double > J_w(N_w,0.0);
    
    double de,t,x,ep;
    double temp;

    
    double K0 =1.0/( pow(lambda,1.0+s) *  tgamma(1.0+s) );
    
    de = 1.0/double(N_e);
    for( int i=0;i<N_e;i++)
    {
        x = (i+0.5)*de;
        ep = (1.0-x)/x;
        
        temp =coeff*2.0 * de* K0 * pow(ep,1.0+s) * exp( -ep/lambda )/x/x;
        
        for( int iw =0 ; iw<N_w;iw++)
        {
            double nu_n = 2.0 * iw * M_PI/beta;
            J_w[iw] += temp / ( ep*ep + nu_n*nu_n );
        }
    }
    
    
    
    for( int it =0 ;it<N_t+1; it++)
    {
        t= it/double(N_t) * beta;
        
        for( int iw = 1; iw<N_w; iw++ )
        {
            double nu_n = 2.0 * iw * M_PI/beta;
            J_[it] += 2.0* J_w[iw] * ( cos(nu_n * t ) );
        }
        J_[it] += J_w[0] ;
        J_[it] /= beta;
    }
//
//    de = 1.0/double(N_e);
//    for( int i=0;i<N_e;i++)
//    {
//        x = (i+0.5)*de;
//        ep = (1.0-x)/x;
//
//        temp =coeff*de* K0* pow(ep,s) * exp( -ep/lambda )/x/x;
//
//
//        for( int it =0 ; it<N_t+1;it++)
//        {
//            t= it/double(N_t) * beta;
//
//
//            J_[it] += temp * ( exp(-t*ep) + exp( -(beta-t) * ep )    )/(1.0-1.0/(exp(beta*ep)));;
//        }
//    }
   
    

    return J_;
    
}


vector<double> boson_bath_w(int N_w, double s, double lambda, double coeff, double beta)
{// rho(w) = k_0 w^s e^{-w/lambda} theta(w) theta( 10*lambda- w)
    // bath = bath * coeff
    // matsu bara
    double N_e = 10*N_w;
    vector< double > J_w(N_w,0.0);
    
    double K0 =1.0/( pow(lambda,1.0+s) *  tgamma(1.0+s) );
    
    double x;
    double de,ep;
    double temp;
    double nu_n;
    
    de = 1.0/double(N_e);
    for( int i=0;i<N_e;i++)
    {
        x = (i+0.5)*de;
        ep = (1.0-x)/x;
        
        temp =coeff*2.0 * de* K0 * pow(ep,1.0+s) * exp( -ep/lambda )/x/x;
        
        for( int iw =0 ; iw<N_w;iw++)
        {
            nu_n = 2.0 * iw * M_PI/beta;
            J_w[iw] += temp / ( ep*ep + nu_n*nu_n );
        }
    }
    
    
    return J_w;
    
}


vector<double> boson_bath_diff(int N_t, double s, double lambda, double coeff, double beta)
{// rho(w) = k_0 w^s e^{-w/lambda} theta(w) theta( 10*lambda- w)
    // differnetial bath
    // bath = bath * coeff
    int N_e = 100*N_t;
    int N_w= 200*N_t;
    vector< double > J_w(N_w,0.0);
    vector< double > B_t(N_t+1,0.0);
   
    double K0 =1.0/( pow(lambda,1.0+s) *  tgamma(1.0+s) );
    
    double x;
    double de,ep;
    double temp;
    double nu_n;
    
    double t;
    
    de = 1.0/double(N_e);
    for( int i=0;i<N_e;i++)
    {
        x = (i+0.5)*de;
        ep = (1.0-x)/x;
        
        temp =coeff*2.0 * de* K0 * pow(ep,1.0+s) * exp( -ep/lambda )/x/x;
        
        for( int iw =0 ; iw<N_w;iw++)
        {
            nu_n = 2.0 * iw * M_PI/beta;
            J_w[iw] += temp / ( ep*ep + nu_n*nu_n );
        }
    }
    
    
    
    for( int it =0 ;it<N_t+1; it++)
    {
        t= it/double(N_t) * beta;
        
        for( int iw = 1; iw<N_w; iw++ )
        {
            nu_n = 2.0 * iw * M_PI/beta;
            B_t[it] += 2.0* J_w[iw] * ( 1.0 - cos(nu_n * t ) )/( nu_n * nu_n);
        }
        B_t[it] += J_w[0] * t* (t-beta)*0.5;
        B_t[it] /= beta;
    }
    return B_t;
    
}



#endif /* bath_initializer_hpp */
