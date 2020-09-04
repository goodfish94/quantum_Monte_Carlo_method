//
//  trace.hpp
//  cthyb
//
//  Created by 胡昊昱 on 7/15/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#ifndef trace_hpp
#define trace_hpp

#include <stdio.h>
#include <cmath>
#include <vector>

#include "local_config.hpp"
using namespace std;

struct trace_cal
{// trace calculator
    trace_cal( vector< vector< double> > &U_in, double beta_in, int norb_in):U_(U_in),norb(norb_in),fock_space(norb_in,0), beta(beta_in)
    {
        config_it_vec.resize(norb_in);
        overlap.resize( norb,vector<double>(norb,0.0) );
        delta_overlap.resize( norb, vector<double>(norb,0.0) );
        update_orb.resize(norb,0);
        update_dag.resize(norb,0);
    }
    
    double get_eve( ); // return <n|U|n>
    double get_update_eve( int flip_orb); // <n|U|n> after flip
    double get_update_eve( ); // <n|U|n> after flip update_orb
    double log_trace_diff(  vector<local_config> &config_vec, int orb_index, op &node_old, op &node_new );
    // calculate the trace diff after insert node_old node_new segment works only if node_old.dag != node_new.dag
    // if node_old< node_new, then insert [node_old, beta] + [0,node_new]
    // return log(trace diff)
    
    double log_trace_diff(  vector<local_config> &config_vec, vector<int> &orb_vec, vector<bool> &dag_vec, double ts, double te );
    // try to insert O^\dag(ts) O(te)
    // with O = \pi_i d_{orb_vec[i]}^{dag_vec[i]}
    // multi op version
    
    double log_trace_diff(  vector<local_config> &config_vec, vector<int> &orb_vec, vector<bool> &dag_at_ts_vec, vector<double> ts, vector<double> te );
    // try to insert multi sgements
    // works only if all the orb index are diff
    
    
    double cal_log_trace( vector<local_config> &config_vec);// directly calculate log trace
    
    double init_log_trace(  vector<local_config> &config_vec );// calculate log trace and let log_trace be it;, set overlap
    
    inline void update( double delta_log_trace){ log_trace+=delta_log_trace; update_overlap( );  }
   
    
    
    inline double get_log_trace(){ return log_trace;};
    //vector< vector< double> >& get_overlap(){ return overlap;}
    double get_overlap(int i, int j)
    {
        if( i==j ) return overlap[i][i];
        return overlap[i][j] + overlap[j][i];
        
    }
    
    bool check_sanity(vector<local_config> &config_vec );
    
    
private:
    vector< vector<double> > U_;
    double beta;
    int norb;
    double log_trace;
    vector< vector<double> > overlap;// overlap between i-th and j-th orbital = overlap[i][j] + overlap[j][i]
    
    double log_trace_aux( vector<local_config> &config_vec, int orb_index, op &node_old, op &node_new );
    // calculate the trace diff after insert node_old node_new segment works only if node_old.dag != node_new.dag and node_old.t<node_new.t
    // return log(trace diff)
    
    double log_trace_aux(  vector<local_config> &config_vec, vector<int> &orb_vec, vector<bool> &dag_vec, double ts, double te );

    vector< vector<double> > delta_overlap;// difference of overlapping while update certain orbital overlap between i,j orbital delta_overlap[i][j] +delta_overlap[j][i]
    
    vector<int> fock_space; // fock_space
    vector< config_type::iterator > config_it_vec;
    vector<int> update_orb;
    vector<int> update_dag;// =1 if insert dag at te, =-1 others
    
    void update_overlap( )
    {
        for( int i=0; i<norb ;i++)
        {
            for( int j=0;j<norb;j++)
            {
                overlap[i][j]+= delta_overlap[i][j];
            }
        }
        
    }
    
    void set_overlap( vector<local_config> &config_vec);// set up overlap from config_vecs
};





#endif /* trace_hpp */
