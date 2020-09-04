//
//  accumulator.hpp
//  cthyb
//
//  Created by 胡昊昱 on 7/16/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#ifndef accumulator_hpp
#define accumulator_hpp

#include <cmath>
#include <utility>
#include <vector>
#include <set>
#include <map>
#include <fstream>

#include "fermi_bath.hpp"
#include "local_config.hpp"
#include "trace.hpp"
#include "diag_flip_bath.hpp"
#include "off_diag_flip_bath.hpp"
#include "fourier_transformation.hpp"

using namespace std;
class accumulator
{// measure green's function and etc.
    
public:
    accumulator( double beta_in,  int norb_in, vector< set<int> > &vec_index_field_in , int N_t_in, int N_w_in ):vec_index_field(vec_index_field_in),beta(beta_in),norb(norb_in)
    {
        N_t = N_t_in;
        N_w = N_w_in;
        dt = beta_in/N_t;
        
        bose_nu_n.resize( N_w,0.0);
        for( int i=0; i<N_w; i++){ bose_nu_n[i] = 2.0*M_PI*i/beta_in; }
        
        sign_sum_static=0;
        sign_sum_G_one =0;
        sign_sum_diag_flip = 0;
        n_meas_static = 0;
        n_meas_G_one = 0;
        n_meas_diag_flip = 0;
        
        n_i_n_j.resize(norb_in,vector<double>(norb_in,0.0));
        n_i_n_j_mpi.resize(norb_in,vector<double>(norb_in,0.0));
        n_i_n_j_mpi_error.resize(norb_in,vector<double>(norb_in,0.0));
        max_fermi = 1000;
        max_diag_flip = 1000;
        max_off_diag_flip = 1000;
        
        histogram_fermion_bath.resize(norb_in, vector<int>(max_fermi,0));
        histogram_diag_flip_bath.resize(norb_in*norb_in, vector<int>(max_fermi,0));
        
        if_meas_one_tau=false;
        if_meas_two_tau=false;
        if_meas_one_w=false;
        if_meas_two_w=false;
        if_meas_nn_tau = false;
        if_meas_diag_flip_tau = false;
//        if_average = false;
        
        if_avg_mpi = false;
        
        orb_to_block.clear();
        
        for( int i=0 ;i<static_cast<int>(vec_index_field_in.size());i++)
        {
            for( set<int>::iterator it = vec_index_field_in[i].begin(); it!=vec_index_field_in[i].end();it++ )
            {
                orb_to_block.insert(  pair<int,int> (*it, i) );
            }
            
        }
        
    }
    
   
    
    void set_meas_one_tau(  );
    void set_meas_two_tau( vector< vector<int> > G_two_index_in  );
    void set_meas_diag_flip_tau( map< pair<int,int>,int> orb_to_index);
    void set_meas_diag_flip_w( map< pair<int,int>,int> orb_to_index );
    void set_meas_n_n_tau(  );
    void set_meas_off_diag_flip_tau(vector< vector<int> > index_to_orb_off_diag_flip_in  );
    void set_meas_off_diag_flip_w(vector< vector<int> > index_to_orb_off_diag_flip_in  );
    
    void measure_static( vector<local_config> &config_vec  , trace_cal& trace_, diag_flip_bath& diag_J_bath,off_diag_flip_bath& off_diag_J_bath, int sign);
    void measure_g_one( vector<fermi_bath_block> &fermi_bath , int sign);
    void measure_diag_flip_bath( diag_flip_bath &diag_J_bath,int sign);
    void measure_off_diag_flip_bath( off_diag_flip_bath &off_diag_J_bath,int sign);
    void measure_n_n( vector<local_config> &config_vec , int sign);
    void measure_g_two( vector<fermi_bath_block> &fermi_bath, int sign);
//    void average_all();
    
    bool get_if_meas_one_tau() const{ return if_meas_one_tau; }
    bool get_if_meas_two_tau() const{ return if_meas_two_tau; }
    bool get_if_meas_one_w() const{ return if_meas_one_w; }
    bool get_if_meas_two_w() const{ return if_meas_two_w; }
    bool get_if_meas_diag_flip_tau() const{ return if_meas_diag_flip_tau;}
    
    void print_G_one();
    void print_G_two();
    void print_G_two_diag_flip();
    void print_G_two_off_diag_flip();
    void print_static();
    void print_n_n();
    
    void print_G_one_mpi();
    void print_G_two_mpi();
    void print_G_two_diag_flip_mpi();
    void print_static_mpi();
    void print_n_n_mpi();
    
     bool average_mpi();
    
    void print_all()
    {
        print_static();
        if( if_meas_one_tau){ print_G_one();}
        if( if_meas_diag_flip_tau || if_meas_diag_flip_w ){  print_G_two_diag_flip();}
        if( if_meas_off_diag_flip_tau || if_meas_off_diag_flip_w  ){print_G_two_off_diag_flip(); }
        if( if_meas_two_tau){ print_G_two();}
        if( if_meas_nn_tau){ print_n_n(); }
    }
    
    void print_all_mpi();
    

    double avg_fermion_bath_num();
    double avg_diag_flip_bath_num();

    
    
    
private:
    
    vector< double > bose_nu_n;// omega_n for boson
    
    vector< vector< fermi_type> > G_one;// one particle green's function
    vector< vector< fermi_type> > G_one_mpi;// one particle green's function
    vector< vector< fermi_type> > G_one_mpi_error;// one particle green's function
    vector< vector< fermi_type> > G_two;// two particle green's function
    vector< vector< fermi_type> > G_two_mpi;// two particle green's function
    vector< vector< fermi_type> > G_two_mpi_error;// two particle green's function
    vector< vector< double > > G_two_diag_flip;// two particle green's function corresponding to the diag flip bath;/
    vector< vector< double > > G_two_diag_flip_mpi;// two particle green's function corresponding to the diag flip bath;/
    vector< vector< double > > G_two_diag_flip_mpi_error;// two particle green's function corresponding to the diag flip bath;/
    
    vector< vector< double > > G_two_off_diag_flip;// two particle green's function corresponding to the diag flip bath;/
    vector< vector< double > > G_two_off_diag_flip_mpi;// two particle green's function corresponding to the diag flip bath;/
    vector< vector< double > > G_two_off_diag_flip_mpi_error;// two particle green's function corresponding to the diag flip bath;/
    
    vector< vector< double > > n_n_;// density density correlation
    vector< vector< double > > n_n_mpi;// density density correlation
    vector< vector< double > > n_n_mpi_error;// density density correlation
    vector< vector< int > > n_i_t;// n_i_t[i][t], fock space  at orbital i time t at this measurement. used to cal n_n_tau;
    map< pair<int,int> ,int> orb_to_index_nn;// orbital to n_n_t index
    
    
    
    
    vector< vector< double > > n_i_n_j;// <n_in_j>
    vector< vector< double > > n_i_n_j_mpi;// <n_in_j>
    vector< vector< double > > n_i_n_j_mpi_error;// <n_in_j>
    
    
    
    map< pair<int, int>, int > G_one_index;// store <i,j>,k  d_i^\dag d_j store in G_one[k]
    map< pair<int,int> , int > diag_flip_index;
    vector< vector<int> > G_two_orb_index;//c_1^\dag c_2(\tau) c_3^\dag c_4
    vector< vector<int> > G_two_block_index;// block index of the orb
    
    vector< vector<int> > histogram_fermion_bath; // for each orbital
    vector< vector<int> > histogram_diag_flip_bath;
    vector< vector<int> > histogram_off_diag_flip_bath;
    
    vector< set<int> > vec_index_field; //orb index for each block
    map< int, int> orb_to_block;// orbital to block index;
    
    vector< vector<int> >  index_to_orb_off_diag_flip;
    
    bool if_avg_mpi;
    
    bool if_meas_one_tau;
    bool if_meas_two_tau;
    bool if_meas_one_w;
    bool if_meas_two_w;
    bool if_meas_nn_tau;
    
    
    bool if_meas_diag_flip_tau;
    bool if_meas_diag_flip_w;
    bool if_meas_off_diag_flip_tau;
    bool if_meas_off_diag_flip_w;
//
//    bool if_average;
    
    int norb;
    
    int N_t;
    double dt;
    double beta;
    int N_w;
    
    int max_fermi; // max num of fermi bath segement
    int max_diag_flip; // max num of daig flip bath
    int max_off_diag_flip;
    
    double sign_sum_static_mpi;
    int n_meas_static_mpi;
    
    double sign_sum_static;//
    double sign_sum_G_one;//
    double sign_sum_G_two;//
    double sign_sum_diag_flip;
    double sign_sum_off_diag_flip;
    double sign_sum_n_n_;
    int n_meas_static;//
    int n_meas_G_one;
    int n_meas_G_two;
    int n_meas_diag_flip;
    int n_meas_off_diag_flip;
    int n_meas_n_n_;
    

    
    
    void measure_two_aux_tau(vector<fermi_bath_block> &fermi_bath, int sign, int type);// measure type-th two
    
    
    
    
};

#endif /* accumulator_hpp */
