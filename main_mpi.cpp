////
////  main.cpp
////  cthyb
////
////  Created by 胡昊昱 on 7/15/19.
////  Copyright © 2019 胡昊昱. All rights reserved.
////
//
//#include <stdio.h>
//#include <iostream>
//#include <mpi.h>
//
//#include "fermi_bath.hpp"
//#include "local_config.hpp"
//#include "trace.hpp"
//#include "random.hpp"
//#include "solver.hpp"
//#include "bath_initializer.hpp"
//
//
//
//int main()
//{
//    double U =0.1;
//    double ed = -U/2;
//    double v;
//    vector< vector<double> > U_in(2,vector<double>(2,0.0));
//    double beta_in = 800.0;
//    int N_t = 800;
//    int N_w = 800;
//    int num_of_block = 2;
//    int norb_in = 2;
//
//    double gamma = 0.5;
//    double D=1.0;
//    double lambda = 1.0, s= 0.6;
//    double g = 0.5;
//
//    v = gamma*2.0*D/M_PI;
//    v= sqrt(v);
//    vector<double> delta_t = fermion_bath_flat(N_t, D, v, beta_in);
//    vector<double> J_t = boson_bath(N_t, s, lambda, g*g, beta_in);
//    vector<double> B_t = boson_bath_diff(N_t, s, lambda, g*g/4.0, beta_in);
//    vector<double> J_w_for_B = boson_bath_w(N_w, s, lambda, g*g/4.0, beta_in);
//
//
//    ed = ed - 0.25*g*g/(lambda * s );
//    U = U +  0.5*g*g/(lambda * s );
//
//
//
//    U_in[0][1] = U/2;
//    U_in[1][0] = U/2;
//    U_in[0][0] = ed;
//    U_in[1][1] = ed;
//
//
//    std::vector< std::vector< fermi_type > >  delta_in(4,std::vector< fermi_type >(0) );
//
//    vector< set<int> > vec_index_field;
//    vec_index_field.resize(2);
//    vec_index_field[0] = { 0 };
//    vec_index_field[1] = { 1 };
//    delta_in[0] = delta_t;
//    delta_in[3] = delta_t;
//
//    vector< vector<double> > B_t_in(3, vector<double>(N_t +1 ,0.0));
//    for( int i=0; i<2; i++){ B_t_in[i] = B_t;}
//    for( int i=0;i<N_t+1; i++){ B_t_in[2][i] = -B_t[i];}
//    map<pair<int, int>, int> orb_to_index ={};
//    orb_to_index.insert( make_pair(pair<int,int>(0,0), 0) );
//    orb_to_index.insert( make_pair(pair<int,int>(1,1), 1) );
//    orb_to_index.insert( make_pair(pair<int,int>(0,1), 2) );
//
//
//
//
//    vector< vector< double> > J_t_in(1,vector<double>(N_t+1, 0.0) );
//    J_t_in[0] = J_t;
//    map<pair<int, int>, int> orb_to_index_diag_flip ={};
//    orb_to_index_diag_flip.insert(make_pair( pair<int,int>(0,1), 0) );
//
//    vector< vector< int> > g_two_orb_index;
//    g_two_orb_index.clear();
//    g_two_orb_index.push_back({0,1,1,0});
//
//
//    int seed = 10;
//
//
//
//
//
//
//    int n_thermal = 10000;
//    int i_meas = 100;
//    int meas_total = 50000;
//
//
//
//    int rank = 0;
//    int size = 0;
//
//    MPI_Init(NULL,NULL);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    
//    seed = seed + rank * 1000;
//
//
//    solver main_solver(beta_in, U_in ,norb_in, delta_in,  num_of_block, vec_index_field,seed,N_t,N_w);
//    main_solver.set_density_bath(B_t_in, orb_to_index, N_t);
//    main_solver.set_diag_flip_bath(J_t_in, orb_to_index_diag_flip, N_t, true,true);
//
//    main_solver.set_measure_n_n();
//    main_solver.set_measure_g_one();
//    main_solver.set_meas_g_two(g_two_orb_index);
//
//    main_solver.mixed_update(n_thermal );
//    
//    i_meas = main_solver.obatin_i_meas();
//    for( int i=0;i<meas_total;i++)
//    {
//        main_solver.mixed_update(i_meas);
//        main_solver.meas_static();
//        main_solver.meas_g_two_diag_flip_bath();
//        main_solver.meas_g_two();
//        main_solver.meas_n_n();
//        main_solver.meas_g_one();
//    }
//
//    main_solver.average();
//    main_solver.average_mpi();
//    if( rank==0){ main_solver.print_all_mpi(); }
////    main_solver.print();
//
//
//    MPI_Finalize();
//
//    return 1;
//
//
//}
