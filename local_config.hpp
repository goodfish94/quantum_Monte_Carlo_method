//
//  local_config.hpp
//  cthyb
//
//  Created by 胡昊昱 on 2019/7/6.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#ifndef local_config_hpp
#define local_config_hpp

#include <stdio.h>
#include <vector>
#include <set>
#include <utility>
#include <iostream>
#include <numeric>

#include "operator.hpp"

using namespace std;
typedef set<op,comp_op> config_type;

struct local_config
{// store segment pic for each orbital
    
public:
    local_config( int num_of_type, int num_of_subtype, double beta_in):config(),beta(beta_in)
    {
        n_type.resize(num_of_type,vector<int>(num_of_subtype,0));
        n_=0;
        remove_lock = true; insert_lock = true;shift_lock = true;find_lock = true; find_lock_2 = true;
        if_full=false;if_wrap=false;
        
    }
    
    double try_insert( op &new_node );// return l_max
    bool insert( op &node_st, op &node_en );// insert node_st node_en, if node_st>node_en, then wrap=!wrap;
    
    op find_k_th_op( int k, bool dag, int type, int subtype,double &down_lim, double &up_lim);// find k-th op with given type dag
    op find_k_th_op( int k, bool dag, int type, int subtype);// find k-th op with given type dag
    //bool find_op_after_t(double t, bool dag, int type, int subtype, op &op_find);// find the op after time t with given type dag and transform value to op_find
    bool find_op_before_t(double t, bool dag, int type, int subtype, op &op_find);// find the op before time t with given type dag and transform value to op_find
    bool find_op(op op_find,double &down_lim, double &up_lim);// try to find op
    bool find_op(op op_find);// try to find op
    bool reaplace_find( op new_op); // replace  the find op with new one
    
    op find_k_th_op_2( int k, bool dag, int type, int subtype,double &down_lim, double &up_lim);// find k-th op with given type dag
    op find_k_th_op_2( int k, bool dag, int type, int subtype);// find k-th op with given type dag
    //bool find_op_after_t(double t, bool dag, int type, int subtype, op &op_find);// find the op after time t with given type dag and transform value to op_find
    bool find_op_before_t_2(double t, bool dag, int type, int subtype, op &op_find);// find the op before time t with given type dag and transform value to op_find
    bool find_op_2(op op_find);// try to find op
    bool find_op_2(op op_find,double &down_lim, double &up_lim);// try to find op
    bool reaplace_find_2( op new_op); // replace  the find op with new one
    
    // if needed, obtain the limi of shift that op
    
    double try_remove(int type_remove, int subtype_remove, int k_remove, bool if_anti_seg,op &node_one, op&node_two  );    //try remove  k-th  given-type element, false if can't remove
    double try_remove( op node_one, op node_two  ); // try remove seg node_one, node_two, return lmax
    
    bool  enforce_remove(op &node_1, op &node_2);// remove node_1 and node after it w/o propose
    bool remove();// remove node according to remove_it_2;
    // all the insert remove action won't check sanity()
    
    bool try_shift(int type_shift,int subtype_shift, int k_shift,bool if_dag,double &down_lim, double &up_lim ,op &node);
    // shift type_shift, k-th dag node, return corresponding node
    bool try_shift(op shift_node,double &down_lim, double &up_lim);// return shift region t-> t+[down_lim,up_lim];
    bool shift(double delta_tau);// shift after try, delta_tau \in[ down_lim, up_lim]
    bool shift(op node_new);// shift to node_new,may also change the type of node
    
    
    inline void lock_insert(){ insert_lock = true;}
    inline void lock_remove(){ remove_lock = true;}
    inline void lock_shift(){ shift_lock = true;}
    inline void lock_find(){ find_lock = true;}
    inline void lock_find_2(){ find_lock_2 = true;}
    inline void lock_all(){remove_lock = true; insert_lock = true;shift_lock = true;find_lock = true;find_lock_2=true;}
    
    bool set_empty_full(){ if( n_==0 && if_full == false ){ if_full = true; return true; } else{ return false; }}
    
    
    double get_n_avg();// calculate len(segment)/beta;
    
    int get_seg_num() const{ return n_/2;}
    int get_op_num() const{  return n_;}
    int get_seg_num(int i) const{ return accumulate( n_type[i].begin(), n_type[i].end(), 0)/2;}
    int get_seg_num(int i, int j) const{ return n_type[i][j]/2; }
    // seg num = op num/2;
    int get_num_of_type() const{ return static_cast<int> ( n_type.size());}
    double get_beta() const{ return beta;}
    bool get_if_wrap() const{ if( n_==0) return if_full; return if_wrap;}
    bool get_if_full() const{ return if_full;}
    config_type& get_config() { return config;}
    
    
    bool check_sanity();
    void print_config();
    
    friend struct trace_cal; // trace calculator
    friend struct density_bath; // density bath trace calculator
    
    
    
private:
    config_type config;
    int n_;// total number of op;
    vector< vector< int > > n_type;// ( number of i-th type p ,j_th subtype);
    double beta;
    bool if_wrap;
    bool if_full;// if full filled
    
    op node_remove_1, node_remove_2;// save the node to be removed

    
    
    config_type::iterator insert_hint, remove_hint, shift_hint;
    config_type::iterator find_hint, find_hint_2; // iterator for find k_th function

    bool shift_lock;
    bool remove_lock;
    bool insert_lock;
    bool find_lock;
    bool find_lock_2;
    
    inline void update_wrap(){
        if(n_ ==0 ){if_wrap =false; return;}
        if( config.begin()->dag ){if_wrap =false;}
        else{ if_wrap = true;}
    }
    
};

#endif /* local_config_hpp */
