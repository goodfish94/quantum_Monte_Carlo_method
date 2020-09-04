////
////  main.cpp
////  cthyb
////
////  Created by 胡昊昱 on 7/15/19.
////  Copyright © 2019 胡昊昱. All rights reserved.
////
//
#include <stdio.h>
#include <iostream>
#include <utility>
#include <vector>

#include "random.hpp"

#include "fermi_bath.hpp"
#include "local_config.hpp"
#include "trace.hpp"
#include "random.hpp"
#include "solver.hpp"
#include "diag_flip_bath.hpp"
#include "bath_initializer.hpp"


using namespace std;

void solver::test_update( int n_sweep)
{// for debug
    
    double rand_num;
    
    for( int i=0; i<n_sweep ;i++)
    {
//        diag_flip_bath_update(1);
        fermion_bath_update(1);
//        if( !check_sanity() )exit(100);
        continue;
        rand_num = rand.random_double(1.0);
        if( rand_num<0.3 )
        {
            fermion_bath_update(1);
            
            
        }
        else if(rand_num <0.6)
        {
            diag_flip_bath_update(1);
        }
        else
        {
            fermion_bath_diag_flip_updae(1);
        }
        
        
    }
//    cout<<"------------"<<endl;
//    for( int i=0; i<norb;i++)
//    {
//        config_vec[i].print_config();
//    }
//    cout<<"------------"<<endl;
    
   

    
}




int main()
{
    double U = 0.1;
    double ed = -U/2;
    double v;
    double beta_in = 20.0;
    int N_t = int( beta_in);
    int N_w = N_t;
    int num_of_block;
    int norb_in = 2;

    vector< vector<double> > U_in(norb_in,vector<double>(norb_in,0.0));

    double gamma = 0.3;
    double D=1.0;
    double lambda = 1.0, s= 0.6;
    double g = 0.5;
//    g=0.0;
    v = sqrt( gamma*2.0*D/M_PI );
    vector<double> delta_t = fermion_bath_flat(N_t, D, v, beta_in);
    vector<double> J_t = boson_bath(N_t, s, lambda, g*g/1.0, beta_in);
    vector<double> B_t = boson_bath_diff(N_t, s, lambda, g*g/4.0, beta_in);
   
    ed = ed - 0.25*g*g/(lambda * s );
    U = U +  0.5*g*g/(lambda * s );
  


    

    for( int i=0; i<norb_in/2; i++)
    {

        U_in[i*2][i*2] = ed;
        U_in[i*2+1][i*2+1] = ed;
        U_in[i*2+1][i*2] = U/2.0;
        U_in[i*2][i*2+1] = U/2.0;

    }




    std::vector< std::vector< fermi_type > >  delta_in(norb_in*norb_in,std::vector< fermi_type >(0) );

    vector< set<int> > vec_index_field;
    vec_index_field.resize(norb_in);
    for( int i=0; i<norb_in; i++)
    {
        vec_index_field[i] ={i};
        delta_in[i*norb_in+i] = delta_t;
    }
    num_of_block = norb_in;

    vector< vector<double> > B_t_in(norb_in + norb_in*(norb_in-1)/2, vector<double>(N_t +1 ,0.0));
    map<pair<int, int>, int> orb_to_index ={};

    int index = 0;
    for( int i=0; i<norb_in/2; i++)
    {
        for( int j=0;j<N_t+1; j++){ B_t_in[index][j] = +B_t[j];}
        orb_to_index.insert( make_pair(pair<int,int>(i*2,i*2) , index));
        index++;

        for( int j=0;j<N_t+1; j++){ B_t_in[index][j] = +B_t[j];}
        orb_to_index.insert( make_pair(pair<int,int>(i*2+1,i*2+1) , index));
        index++;

        for( int j=0;j<N_t+1; j++){ B_t_in[index][j] = -B_t[j];}
        orb_to_index.insert( make_pair(pair<int,int>(i*2,i*2+1) , index));
        index++;

    }

//    for( int i=0; i<norb_in/2-1; i++)
//    {
//        for( int j=i+1; j<norb_in/2-1;j++)
//        {
//            B_t_in[index] = B_t;
//            orb_to_index.insert( make_pair( pair<int,int>(i*2,j*2) , index));
//            index++;
//
//            B_t_in[index] = B_t;
//            orb_to_index.insert( make_pair( pair<int,int>(i*2+1,j*2+1) , index));
//            index++;
//
//            for( int k=0;k<N_t+1; k++){ B_t_in[index][k] = -B_t[k];}
//            orb_to_index.insert( make_pair( pair<int,int>(i*2,j*2+1) , index));
//            index++;
//
//            for( int k=0;k<N_t+1; k++){ B_t_in[index][k] = -B_t[k];}
//            orb_to_index.insert( make_pair( pair<int,int>(i*2+1,j*2) , index));
//            index++;
//
//        }
//    }




//
//    vector< vector< double> > J_t_in(norb_in/2,vector<double>(N_t+1, 0.0) );
//    map<pair<int, int>, int> orb_to_index_diag_flip ={};
//    for( int i=0; i<norb_in/2; i++)
//    {
//        J_t_in[i] = J_t;
//        orb_to_index_diag_flip.insert(make_pair( pair<int,int>(i*2,i*2+1), i) );
//    }

//    vector< vector< double> > J_t_off_in( (norb_in/2)*(norb_in/2-1),vector<double>(N_t+1, 0.0) );
//    vector< vector<int> > orb_to_index_off_diag_flip;
//    vector<int> orb_ind;
//    orb_to_index_off_diag_flip.clear();
//
//    index = 0;
//    for( int i=0; i<norb_in/2; i++)
//    {
//        for( int j=i+1; j<norb_in/2; j++)
//        {
//
//            J_t_off_in[index] = J_t;
//            orb_ind = { 2*i, 2*i+1, 2*j+1, 2*j  };
//            orb_to_index_off_diag_flip.push_back(orb_ind);
//            index++;
//
//            J_t_off_in[index] = J_t;
//            orb_ind = { 2*j, 2*j+1, 2*i+1, 2*i  };
//            orb_to_index_off_diag_flip.push_back(orb_ind);
//            index++;
//
//
//
//
//        }
//
//    }
//




    int seed = 40;


    solver test( beta_in, U_in ,norb_in, delta_in,  num_of_block, vec_index_field,seed,N_t,N_w);

    test.set_density_bath(B_t_in, orb_to_index, N_t);
//    test.set_diag_flip_bath(J_t_in, orb_to_index_diag_flip, N_t, true,true);
//    test.set_off_diag_flip_bath(J_t_off_in, orb_to_index_off_diag_flip, N_t, true,true);

    vector< vector< int> > g_two_orb_index;
    g_two_orb_index.clear();
    for( int i=0; i<norb_in/2; i++)
    {
        for( int j=i ; j<norb_in/2; j++)
        {
            g_two_orb_index.push_back({ 2*i, 2*i+1, 2*j+1,2*j});
        }
    }
    g_two_orb_index.push_back({0,1,1,0});

    int n_thermal = 300000;
    int i_meas = 100 ;
    int meas_total = 500000;

    test.set_measure_n_n();
    test.set_measure_g_one();
    test.set_meas_g_two(g_two_orb_index);

//    test.set_empty_orb_full(1);
//    test.set_empty_orb_full(3);

    test.test_update(n_thermal);
//    test.fermion_bath_update(n_thermal );
//    test.off_diag_flip_bath_update(n_thermal);
//    test.diag_flip_bath_update(n_thermal);
//    test.fermion_bath_diag_flip_updae(n_thermal);

    test.iter_numer = 0;
    for( int i=0;i<meas_total;i++)
    {
//        cout<<" iter="<<i<<endl;
        test.iter_numer ++;
        test.test_update(i_meas);
//        if( !test.check_sanity() )exit(100);
//        test.diag_flip_bath_update(i_meas);
//        test.off_diag_flip_bath_update(1);
//        test.fermion_bath_off_diag_flip_updae(1);
//        test.fermion_bath_update(i_meas);
//        test.fermion_bath_diag_flip_updae(i_meas);
        
        test.meas_static();
//        test.meas_g_two_diag_flip_bath();
//        test.meas_g_two_off_diag_flip_bath();
//        test.meas_g_two();
        test.meas_n_n();
        test.meas_g_one();
//
        
    }


    test.print();



    return 1;


}
