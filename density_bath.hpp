//
//  density_bath.hpp
//  cthyb
//
//  Created by 胡昊昱 on 7/25/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#ifndef density_bath_hpp
#define density_bath_hpp

#include <cmath>
#include <vector>
#include <map>
#include <utility>
#include <list>

#include "local_config.hpp"
#include "operator.hpp"

using namespace std;

class density_bath
{
    
public:
    
    
    density_bath( vector< vector< double> > &B_ij_tau_in,  map< pair<int,int>, int>& orb_to_index_in, int &N_t_in,double beta_in,int norb_in):B_ij_tau(B_ij_tau_in), orb_to_index(orb_to_index_in), N_t(N_t_in),beta(beta_in),norb(norb_in)
    {
        log_trace=0.0;
    }
    density_bath(int norb_in,double beta_in)
    {
        B_ij_tau.clear();
        orb_to_index.clear();
        N_t=1;
        beta=beta_in;
        norb=norb_in;
        log_trace = 0.0;
        
    
        
    }
    
    void set_density_bath( vector< vector< double> > &B_ij_tau_in,  map< pair<int,int>, int>& orb_to_index_in, int &N_t_in,double beta_in,int norb_in )
    {
        B_ij_tau = B_ij_tau_in;
        orb_to_index = orb_to_index_in;
        N_t = N_t_in;
        beta = beta_in;
        norb = norb_in;
        log_trace=0.0;
        
        
//        for( auto it = orb_to_index.begin(); it!= orb_to_index.end(); it++)
//        {
//            cout<<it->first.first<<" "<<it->first.second<<" "<< it->second<<endl;
//        }
//        for( auto it1 = B_ij_tau.begin(); it1!=B_ij_tau.end();it1++)
//        {
//            for( auto it2 = it1->begin(); it2 != it1->end(); it2++)
//            {
//                cout<<*it2<<endl;
//            }
//            cout<<"-----"<<endl;
//        }
//
        
    }
    
    double log_trace_diff_insert( vector<local_config> &config_vec, op &op_one, op &op_two , int& orb);
    double log_trace_diff_remove( vector<local_config> &config_vec, op &op_one, op &op_two , int& orb);
    double log_trace_diff_replace( vector<local_config> &config_vec, op &op_old, op &op_new, int& orb);
    // trace diff after insert, remove ,replace;
    
    double log_trace_diff_insert( vector<local_config> &config_vec, vector<int>& orb_vec, vector<bool>& dag_vec , vector<double> & t);
    // insert multi orb
    
    double log_trace_diff_remove( vector<local_config> &config_vec, vector<int>& orb_vec, vector<bool>& dag_vec , vector<double> & t);
    // remove
    
    double log_trace_diff_replace( vector<local_config> &config_vec, vector<int>& orb_vec, vector<bool>& dag_vec , vector<double> & t_old, vector<double> &t_new);
    
    void update_log_trace( double& log_trace_diff ){ log_trace += log_trace_diff; }
    double get_log_trace() const{return log_trace;}
    
    bool check_sanity(vector<local_config> &config_vec );
    
    
    
private:
    
    vector< vector< double> > B_ij_tau;// store B_ij_tau // first index is the orbital indx, second index is the tau index
    map< pair<int,int>, int> orb_to_index;// orbital index to bath index
    
    inline double interpolate_B_ij_t(double t, int index)
    {
        
        if(t<0.0){ t = -t;}
        double dt = beta/N_t;
        int k = floor(t/dt);
        return ( B_ij_tau[index][k] + (t-k*dt)*(B_ij_tau[index][k+1]-B_ij_tau[index][k])/dt );
    }
    int N_t;
    double beta;
    int norb;
    
    double log_trace;// current log trace;
};

#endif /* density_bath_hpp */
