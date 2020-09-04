//
//  solver.hpp
//  cthyb
//
//  Created by 胡昊昱 on 7/15/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#ifndef solver_hpp
#define solver_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <fstream>

#include "accumulator.hpp"
#include "random.hpp"
#include "local_config.hpp"
#include "fermi_bath.hpp"
#include "trace.hpp"
#include "density_bath.hpp"



using namespace std;
struct solver
{
    
    
public:
    solver( double beta_in, vector<vector<double> > &U_in , int norb_in, vector<vector<double> > &delta_in, int num_of_block, vector< set<int> > &vec_index_field,int seed, int N_t, int N_w)
    :beta(beta_in),norb(norb_in),
    config_vec(norb_in, local_config( 4,norb_in*norb_in,beta_in) ),fermi_bath(num_of_block,fermi_bath_block(beta_in, norb_in,N_t)),trace_(U_in, beta_in,norb_in),rand(seed), accu_(beta_in,norb_in, vec_index_field,N_t,N_w), B_bath(norb_in,beta_in), diag_J_bath(beta_in), off_diag_J_bath(beta_in)
    {
        fermi_bath[0].input_delta(delta_in);
        for( int i=0; i< num_of_block; i++)
        {
            fermi_bath[i].input_index_field( vec_index_field[i]);
        }
        trace_.init_log_trace(config_vec);
        orbital_to_block.clear();
        index_to_orb_diag_flip.clear();
        
        for( int i=0; i<num_of_block;i++)
        {
            for( set<int>::iterator  it=vec_index_field[i].begin(); it!=vec_index_field[i].end(); it++)
            {
                orbital_to_block.insert( make_pair(*(it), i ) );
            }
        }
        propose_fermi_bath.resize(3,0);
        accept_fermi_bath.resize(3,0);
        propose_diag_flip_bath.resize(4,0);
        accept_diag_flip_bath.resize(4,0);
        propose_fermion_diag_flip_bath.resize(3,0);
        accept_fermion_diag_flip_bath.resize(3,0);
        propose_off_diag_flip_bath.resize(4,0);
        accept_off_diag_flip_bath.resize(4,0);
        propose_fermion_off_diag_flip_bath.resize(3,0);
        accept_fermion_off_diag_flip_bath.resize(3,0);
        sign_ = 1;
    }
    
    void set_density_bath( vector< vector< double> > &B_ij_tau_in,  map< pair<int,int>, int>& orb_to_index_in, int &N_t_in )
    {
        B_bath.set_density_bath(B_ij_tau_in,  orb_to_index_in, N_t_in,beta ,norb );
    }
    void set_diag_flip_bath( vector< vector< double> > &J_tau_in,  map< pair<int,int>, int>& orb_to_index_in, int &N_t_in, bool if_meas, bool if_meas_tau )
    {
        diag_J_bath.set_diag_flip_bath(J_tau_in, orb_to_index_in, N_t_in, beta);
        index_to_orb_diag_flip.resize( static_cast<int>(orb_to_index_in.size() )  );
        for( map< pair<int,int> ,int> ::iterator it = orb_to_index_in.begin(); it!=orb_to_index_in.end(); it++ )
        {
            index_to_orb_diag_flip[it->second] = it->first;
        }
        if( if_meas )
        {
            if( if_meas_tau){accu_.set_meas_diag_flip_tau( diag_J_bath.get_orb_to_index() );}
            else{ accu_.set_meas_diag_flip_w(diag_J_bath.get_orb_to_index() );  }
        }
        
    }
    
    void set_off_diag_flip_bath( vector< vector< double> > &J_tau_in,  vector< vector<int> >& index_to_orb_off_diag_flip_in , int &N_t_in, bool if_meas, bool if_meas_tau )
    {
        off_diag_J_bath.set_off_diag_flip_bath(J_tau_in, index_to_orb_off_diag_flip_in, N_t_in, beta);
        index_to_orb_off_diag_flip = index_to_orb_off_diag_flip_in;

        if( if_meas )
        {
            if( if_meas_tau){accu_.set_meas_off_diag_flip_tau( index_to_orb_off_diag_flip_in  );}
            else{ accu_.set_meas_off_diag_flip_w( index_to_orb_off_diag_flip_in );  }
        }

    }


    
    
    void set_measure_g_one(){ accu_.set_meas_one_tau(); }
    void set_meas_g_two( vector< vector<int> > G_two_index_in  ){  accu_.set_meas_two_tau( G_two_index_in  );  }
    void set_measure_n_n(){ accu_.set_meas_n_n_tau(); }
   
    int obatin_i_meas();// get i_meas according to the current information
    
    bool fermion_bath_update_insert();
    bool fermion_bath_update_remove();
    bool fermion_bath_update_shift();
    
    bool diag_flip_bath_update_insert();
    bool diag_flip_bath_update_remove();
    bool diag_flip_bath_update_shift();
    bool diag_flip_bath_update_swap();
    
    bool off_diag_flip_bath_update_insert();
    bool off_diag_flip_bath_update_remove();
    bool off_diag_flip_bath_update_shift();
    bool off_diag_flip_bath_update_swap();


    bool fermion_bath_to_off_diag_flip_update();
    bool off_diag_flip_to_fermion_bath_update();
    bool swap_fermion_bath_off_diag_flip_update();
    
    bool fermion_bath_to_diag_flip_update();
    bool diag_flip_to_fermion_bath_update();
    bool swap_fermion_bath_diag_flip_update();

    
    void mixed_update( int n_sweep );
    void test_update( int n_sweep); // for debug 
    
    void fermion_bath_update( int n_sweep);// update for n_sweep times;
    void diag_flip_bath_update( int n_sweep);// update for n_sweep times;
    void fermion_bath_diag_flip_updae(int n_sweep);
    void off_diag_flip_bath_update(int n_sweep);
    void fermion_bath_off_diag_flip_updae(int n_sweep);
    
    int get_sign() const {return sign_;};
    
    void meas_g_two_diag_flip_bath(){ accu_.measure_diag_flip_bath(diag_J_bath, sign_); }
    void meas_g_two_off_diag_flip_bath(){ accu_.measure_off_diag_flip_bath(off_diag_J_bath, sign_); }
    void meas_g_one(){  accu_.measure_g_one(fermi_bath, sign_); }
    void meas_g_two(){  accu_.measure_g_two(fermi_bath, sign_ ); }
    void meas_static(){ accu_.measure_static(config_vec, trace_,diag_J_bath,off_diag_J_bath, sign_);   }
    void meas_n_n(){ accu_.measure_n_n(config_vec,sign_);}
    
    
    void print();
    
    void average_mpi(){ accu_.average_mpi();}
    void print_all_mpi(){accu_.print_all_mpi();}
    
    bool set_empty_orb_full( int i)
    {
        if( !config_vec[i].set_empty_full())
        {
            return false;
        }
        trace_.init_log_trace(config_vec);
        return true;
    }
    
    bool check_sanity( );
    
    int iter_numer ; 
    
    
    
    
    
    
private:
    double beta;
    int norb;
    
    int sign_;
    
    map<int, int> orbital_to_block;
    vector< pair<int,int> > index_to_orb_diag_flip;// index to orbital for diag flip bath
    vector< vector<int> > index_to_orb_off_diag_flip;
    
    trace_cal trace_;
    
    vector<local_config> config_vec;
    vector<fermi_bath_block> fermi_bath;
    density_bath B_bath;
    diag_flip_bath diag_J_bath;
    off_diag_flip_bath off_diag_J_bath;
    
    
    accumulator accu_;
    rand_gen rand;
    
    
    op op_one, op_two;
    op op_1,op_2,op_3,op_4;
    op op_1_dag,op_2_dag,op_3_dag,op_4_dag;
    pair<int,double> pair_one,pair_two;
    double lmax, lmax2;
    double lmax3,lmax4;
    
    vector<int> propose_fermi_bath;
    vector<int> accept_fermi_bath;
    vector<int> propose_diag_flip_bath;
    vector<int> accept_diag_flip_bath;
    
    vector<int> propose_fermion_diag_flip_bath;
    vector<int> accept_fermion_diag_flip_bath;
    
    vector<int> propose_off_diag_flip_bath;
    vector<int> accept_off_diag_flip_bath;
    
    vector<int> propose_fermion_off_diag_flip_bath;
    vector<int> accept_fermion_off_diag_flip_bath;
   
    
    fermi_type try_insert_multi_op_to_fermi_bath( vector<int> &orb,  vector<double> &ts, vector<double> &te, int &sign_change);
    fermi_type try_remove_multi_op_to_fermi_bath( vector<int> &orb, vector<double> &ts, vector<double> &te, int &sign_change);
    fermi_type try_replace_multi_op_to_fermi_bath( vector<int> &orb, vector<bool> &dag, vector<double> &told, vector<double> &tnew,int &sign_change);
    
    bool do_insert_multi_op_to_fermi_bath(int &sign_change, vector<int> &orb, bool if_sucess);
    bool do_remove_multi_op_to_fermi_bath(int &sign_change, vector<int> &orb, bool if_sucess);
    bool do_replace_multi_op_to_fermi_bath(int &sign_change, vector<int> &orb,bool if_sucess);
    
    // insert/remove/replace according to input
    // case one: 2 pairs of dddag, 2 pairs may belong to the same fermi block
    // case two: 4 pairs of ddag, 2 pairs spin up, and 2 pairs spin down
    // for replace, dag denotes whether replace d or ddag
    // notice that if tnew change the flag if_wrap, then the sign need to be updated. this function doesn't handle it
    // return the ratio of operation
    // set up the sign_change
    

    
    int calculate_perm_sign();
    
    
    vector<double> st_vec ;
    vector<double> en_vec ;
    vector<int> orb_vec ;
    vector<bool> dag_vec ;
    
    double log_trace_diff;
    double log_trace_diff_B_bath;
    
    const double zero_plus = 1e-7;// 0+, used to enforce time order for c^dag c like operator
    
    
};


#endif /* solver_hpp */
