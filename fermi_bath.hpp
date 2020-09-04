//
//  fermi_bath.hpp
//  cthyb
//
//  Created by 胡昊昱 on 2019/7/4.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#ifndef fermi_bath_hpp
#define fermi_bath_hpp
#include <set>
#include <map>
#include <utility>
#include <vector>
#include <boost/qvm/mat_operations.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "mat_vec.hpp"
#include "operator.hpp"

using namespace std;




class fermi_bath_block
{// fermion bath for each block
    
    
public:
    fermi_bath_block(int norb_in, double beta_in, set< int > index_field_in,int N_t_in):norb(norb_in),beta(beta_in),index_field(index_field_in), N_inv(0,0),N_(0,0)
    {
        nsize = 0;
        time_start.clear();
        time_end.clear();
        det = 1.0;
        N_t = N_t_in;
        dt = beta_in / N_t_in;

        insert_lock = true; remove_lock = true;replace_lock = true;
    }
    
    fermi_bath_block( double beta_in, int norb_in, int N_t_in):norb(norb_in),beta(beta_in), N_inv(0,0),N_(0,0)
    {
        nsize = 0;
        time_start.clear();
        time_end.clear();
        det = 1.0;
        N_t = N_t_in;
        dt = beta_in / N_t_in;
        insert_lock = true; remove_lock = true;replace_lock = true;
        index_field={};
        
    }
    
    void input_delta( vector<vector<fermi_type> > &delta_in)
    {
        delta_tau=delta_in;
    }
    void input_index_field( set<int> &index_field_in)
    {
        using namespace std;
        index_field.clear();
        for(set< int >::iterator it =index_field_in.begin(); it!=index_field_in.end(); it++)
        {
            index_field.insert(*it);
        }
    }
    
    fermi_type try_insert(int orb_one, int orb_second, op op_one,  op op_two)
    {
        pair<int,double> node_st, node_en;
        if( op_one.dag ){  node_st= { orb_one, op_one.t}; node_en = {orb_second, op_two.t};  }
        else{ node_st= {orb_second, op_two.t};  node_en = { orb_one, op_one.t};}

        return try_insert( node_st, node_en );
    }
    
    fermi_type try_insert(int orb_st, int orb_en, double t_st,  double t_en)
    {
        pair<int,double> node_st = { orb_st, t_st}, node_en = {orb_en,t_en};
        return try_insert( node_st, node_en );
    }
    
    fermi_type try_remove(int orb_one, int orb_second, op op_one,  op op_two)
    {
        pair<int,double> node_st, node_en;
        if( op_one.dag ){  node_st= { orb_one, op_one.t}; node_en = {orb_second, op_two.t};  }
        else{ node_st= {orb_second, op_two.t};  node_en = { orb_one, op_one.t};}
   
        return try_remove( node_st, node_en );
    }
    
    fermi_type try_remove(int orb_one, int orb_second, double t_st, double t_en)
    {
        pair<int,double> node_st={orb_one,t_st};
        pair<int,double> node_en={orb_second,t_en};        
        return try_remove( node_st, node_en );
    }
    
    fermi_type try_replace(int orb, double t_old, double t_new, bool if_st)
    {
        pair<int,double> node_old ={ orb,t_old };
        pair<int,double> node_new ={orb,t_new};
        return try_replace(node_old,node_new, if_st);
    }
    
    fermi_type try_insert(pair< int,double> &node_start, pair<int,double> &node_end);
    fermi_type try_remove(pair< int,double> &node_start, pair<int,double> &node_end);
    fermi_type try_remove(int k_st, int k_en);// remove the k-th element
    fermi_type try_replace( pair<int,double> &node_st_old, pair<int,double> &node_st_new,pair<int,double> &node_en_old, pair<int,double> &node_en_new);
    fermi_type try_replace( pair<int,double> &node_old, pair<int,double> &node_new, bool if_st);// if_st = true replace start;
    
    
    // try_** set the needed vector and var
    // return the ratio of weight
    bool do_insert( );
    bool do_remove( int &swap_sign );
    bool do_replace( );
    // do after try
    //do the update
    //swap_sign  = sign change from swap row/col and determinatn
    
    inline void lock_insert(){ insert_lock = true;}
    inline void lock_remove(){ remove_lock = true;}
    inline void lock_replace(){ replace_lock = true;}
    inline void lock_all(){insert_lock = true; remove_lock = true;replace_lock = true;}
    // lock to prevent do_* operation
    
    
    int get_size() const{
        return nsize;
    }
    
   
    
    int get_det_sign() const{  return (det>0 ? (1): (-1));  }
    double get_det() const{ return det;}
    bool check_sanity();// check wheter det, N, N_inv is correct
    //int get_permute_sign();// permute the current time_start, time_end to the time ordered time_start,time_end
    
    map< pair<int,double> ,int>  get_time_start() const{ return time_start; }
    map< pair<int,double> ,int>  get_time_end() const{ return time_end; }
    
    friend class accumulator;
    friend class solver;// debug only
private:
    
    double beta;
    int nsize;// num of segment
    int norb; // number of total orb
    set< int > index_field; // store the orb index in this block
    
    map< pair<int,double> ,int> time_start;
    map< pair<int,double> ,int> time_end;//store time array, map between op and index
    
    static std::vector< std::vector< fermi_type > > delta_tau;//bath function first index, type of bath, second index  tau_index
    
    int N_t;
    double dt;
    
    
    
    inline fermi_type interpolate_delta(const pair<int,double> &t1 ,const pair<int,double> &t2)
    {// interpolate and get delta(t2-t1);
        int index = (t1.first*norb ) + t2.first;
        double t = t1.second  - t2.second;
        int sign = 1.0;
        if(t<0.0)
        {
            t = t+beta;
            sign = -1.0;
        }
 
        int k = static_cast<int> (t/dt);
//        cout<<"dt="<<dt<<endl;
//        cout<<"index = "<<index<<endl;
//        cout<<t1.first<<" "<<t2.first<<endl;
        return sign * ( delta_tau[index][k] + (t-k*dt)*(delta_tau[index][k+1]-delta_tau[index][k])/dt );

    }


    fermi_mat N_inv;
    fermi_mat N_;
    fermi_vec R_insert,Q_insert,R_tilde_insert,Q_tilde_insert;
    fermi_type S_insert,S_tilde_insert;
    fermi_vec R_remove,Q_remove,R_tilde_remove,Q_tilde_remove;
    fermi_type S_remove,S_tilde_remove;
    
    fermi_vec aux_vec;
    fermi_mat aux_mat;
    

    
    

    fermi_type det;

    pair<int,double> node_start_insert;
    pair<int,double> node_end_insert;
    
    pair<int,double> node_start_remove;
    pair<int,double> node_end_remove;
    
    int row_remove, col_remove;// row and column to remove
    
    
    fermi_type try_remove_aux();// remove after set up
    
    bool insert_lock;// if true, then can't do_insert();
    bool remove_lock;// if true, then can't do_remove();
    bool replace_lock;// if true, then can't do_replace();
    
    
};


double cal_det(boost::numeric::ublas::matrix<double> m);
#endif /* fermi_bath_hpp */
