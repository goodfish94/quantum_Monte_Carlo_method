//
//  density_bath.cpp
//  cthyb
//
//  Created by 胡昊昱 on 7/25/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#include "density_bath.hpp"

double density_bath::log_trace_diff_insert( vector<local_config> &config_vec, op &op_one, op &op_two , int &orb)
{
    double diff_log_trace = 0.0;
    int s_1 = op_one.dag ? 1:(-1);
    int s_2 = op_two.dag ? 1:(-1);
    int index;
    


    
    
    for( int i=0;i<norb;i++ )
    {
        map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(   (i<orb)? make_pair(i,orb) : make_pair(orb,i)    );
        if( it_index ==  orb_to_index.end() ){ continue;}
        else{ index = it_index->second; }
        
        for( config_type::iterator it = config_vec[i].config.begin(); it!=config_vec[i].config.end(); it++)
        {
            if( it->dag )
            {
                diff_log_trace -= s_1*interpolate_B_ij_t( op_one.t - it->t, index ) + s_2*interpolate_B_ij_t( op_two.t - it->t, index );
            }
            else
            {
                diff_log_trace += s_1*interpolate_B_ij_t( op_one.t - it->t, index ) + s_2*interpolate_B_ij_t( op_two.t - it->t,index );
            }
        }
        
    }
    
    map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(  make_pair(orb, orb)    );
    if( it_index !=  orb_to_index.end() )
    {
        diff_log_trace -=  s_1*s_2*interpolate_B_ij_t(op_one.t-op_two.t, it_index->second) ;
    }
    diff_log_trace = 2.0*diff_log_trace;
    
    
    
//
//
//    double diff_loc_trace_test = 0.0 ;
//    double sigma_new, sigma;
//    double s;
//    if(orb == 0)
//    {
//        sigma_new = 1;
//    }
//    else
//    {
//        sigma_new = -1;
//    }
//
//
//    for( int i=0;i<norb;i++ )
//    {
//        if( i==0 ){ sigma = 1;}
//        else{ sigma = -1; }
//        for( config_type::iterator it = config_vec[i].config.begin(); it!=config_vec[i].config.end(); it++)
//        {
//            if( it->dag )
//            {
//                s=1;
//            }
//            else
//            {
//                s=-1;
//            }
//            diff_loc_trace_test +=(-2.0) * sigma*sigma_new*s*s_1*( interpolate_B_ij_t(it->t - op_one.t,0) - interpolate_B_ij_t(it->t - op_two.t, 0)  );
//        }
//
//    }
//    diff_loc_trace_test += (-2.0) * (-1.0)*interpolate_B_ij_t(op_two.t - op_one.t, 0);
//    if( abs(diff_loc_trace_test - diff_log_trace) > 1e-6)
//    {
//        cout<<"trace error "<<endl;
//        exit(2);
//    }
//
//
    return diff_log_trace;
    
}

double density_bath::log_trace_diff_remove( vector<local_config> &config_vec, op &op_one, op &op_two,int& orb)
{
    double diff_log_trace = 0.0;
    int s_1 = op_one.dag ? 1:(-1);
    int s_2 = op_two.dag ? 1:(-1);
    int index;
    
    
    for( int i=0;i<norb;i++ )
    {
        map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(   (i<orb)? make_pair(i,orb) : make_pair(orb,i)    );
        if( it_index ==  orb_to_index.end() ){ continue;}
        else{ index =  it_index->second; }
        
        for( config_type::iterator it = config_vec[i].config.begin(); it!=config_vec[i].config.end(); it++)
        {
            if( it->dag )
            {
                diff_log_trace -= s_1*interpolate_B_ij_t( op_one.t - it->t,index ) + s_2*interpolate_B_ij_t( op_two.t - it->t, index );
            }
            else
            {
                diff_log_trace += s_1*interpolate_B_ij_t( op_one.t - it->t, index ) + s_2*interpolate_B_ij_t( op_two.t - it->t, index );
            }
        }
        
    }
    map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(  make_pair(orb, orb) );
    if( it_index !=  orb_to_index.end() )
    {
     diff_log_trace -= -2.0* interpolate_B_ij_t( 0.0 ,  it_index->second ) - s_1*s_2* interpolate_B_ij_t(op_one.t-op_two.t, it_index->second) ;
    }
     // subtract double counting
    diff_log_trace = -2.0*diff_log_trace;// remove
    
    
    
    
    
    
    
    
    
    
    
    double diff_loc_trace_test = 0.0 ;
    double sigma_new, sigma;
    double s;
    if(orb == 0)
    {
        sigma_new = 1;
    }
    else
    {
        sigma_new = -1;
    }
    
    
    for( int i=0;i<norb;i++ )
    {
        if( i==0 ){ sigma = 1;}
        else{ sigma = -1; }
        for( config_type::iterator it = config_vec[i].config.begin(); it!=config_vec[i].config.end(); it++)
        {
            if( it->dag )
            {
                s=1;
            }
            else
            {
                s=-1;
            }
            if( *it == op_one|| *it == op_two ) continue;
            diff_loc_trace_test +=(-2.0) * sigma*sigma_new*s*s_1*( interpolate_B_ij_t(it->t - op_one.t,0) - interpolate_B_ij_t(it->t - op_two.t, 0)  );
        }
        
    }
    diff_loc_trace_test += (-2.0) * (-1.0)*interpolate_B_ij_t(op_two.t - op_one.t, 0);
    
    diff_log_trace = - diff_log_trace;
    if( abs(diff_loc_trace_test - diff_log_trace) > 1e-6)
    {
        cout<<"trace error "<<endl;
        exit(2);
    }
    
    

    
    
    
    
    return diff_log_trace;
}


double density_bath::log_trace_diff_replace( vector<local_config> &config_vec, op &op_old, op &op_new, int& orb)
{
    double diff_log_trace = 0.0;
    int s = op_old.dag ? 1:(-1);
    int index;
    
    
    for( int i=0;i<norb;i++ )
    {
        map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(   (i<orb)? make_pair(i,orb) : make_pair(orb,i)    );
        if( it_index ==  orb_to_index.end() ){ continue;}
        else{ index = it_index->second; }
        
        for(config_type::iterator it = config_vec[i].config.begin(); it!=config_vec[i].config.end(); it++)
        {
            if( it->dag )
            {
                diff_log_trace -= s*(interpolate_B_ij_t( op_new.t - it->t, index ) - interpolate_B_ij_t( op_old.t - it->t, index ) );
            }
            else
            {
                diff_log_trace += s*(interpolate_B_ij_t( op_new.t - it->t, index ) - interpolate_B_ij_t( op_old.t - it->t, index ) );
            }
        }
        
    }
    map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(  make_pair(orb, orb)    );
    if( it_index !=  orb_to_index.end() )
    {
        diff_log_trace -=( -interpolate_B_ij_t(op_new.t - op_old.t, it_index->second) + interpolate_B_ij_t(0.0, it_index->second) )  ;
    }
    
    diff_log_trace = 2.0*diff_log_trace;
    
    return diff_log_trace;
    
    
}




double density_bath::log_trace_diff_insert( vector<local_config> &config_vec, vector<int>& orb_vec, vector<bool>& dag_vec , vector<double> &t)
{
 // insert multi operator
    int index;
    double diff_log_trace = 0.0;
    vector<int> s(static_cast<int>(dag_vec.size()),0 );
    for( int i=0; i<static_cast<int>(dag_vec.size()) ; i++) { s[i] = dag_vec[i]?(1):(-1); }
    
    for( int j=0 ; j<static_cast<int>(dag_vec.size()) ; j++ )
    {
        for( int i=0;i<norb;i++ )
        {
            map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(   (i<orb_vec[j])? make_pair(i,orb_vec[j]) : make_pair(orb_vec[j],i)    );
            if( it_index ==  orb_to_index.end() ){ continue;}
            else{ index = it_index->second; }
        
            for( config_type::iterator it = config_vec[i].config.begin(); it!=config_vec[i].config.end(); it++)
            {
                if( it->dag ){ diff_log_trace -= s[j]*interpolate_B_ij_t( t[j] - it->t, index );}
                else{ diff_log_trace += s[j]*interpolate_B_ij_t( t[j] - it->t, index );}
            }
        
        }
    }
    for( int i=0 ; i<static_cast<int>(dag_vec.size()) ; i++ )
    {
        for( int j=i+1; j<static_cast<int>(dag_vec.size()) ; j++ )
        {
            map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(   (orb_vec[i]<orb_vec[j])? make_pair(orb_vec[i],orb_vec[j]) : make_pair(orb_vec[j],orb_vec[i])    );
            if( it_index ==  orb_to_index.end() ){ continue;}
            else{ index = it_index->second; }
            diff_log_trace -= s[i] * s[j] * interpolate_B_ij_t(t[j] - t[i], index);
        }
    }
    diff_log_trace = 2.0*diff_log_trace;
    return diff_log_trace;
}



double density_bath::log_trace_diff_remove( vector<local_config> &config_vec, vector<int>& orb_vec, vector<bool>& dag_vec , vector<double> & t)
{
    // insert multi operator
    int index;
    double diff_log_trace = 0.0;
    vector<int> s(static_cast<int>(dag_vec.size()),0 );
    for( int i=0; i<static_cast<int>(dag_vec.size()) ; i++) { s[i] = dag_vec[i]?(1):(-1); }
    
    for( int j=0 ; j<static_cast<int>(dag_vec.size()) ; j++ )
    {
        for( int i=0;i<norb;i++ )
        {
            map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(   (i<orb_vec[j])? make_pair(i,orb_vec[j]) : make_pair(orb_vec[j],i)    );
            if( it_index ==  orb_to_index.end() ){ continue;}
            else{ index = it_index->second; }
            
            for( config_type::iterator it = config_vec[i].config.begin(); it!=config_vec[i].config.end(); it++)
            {
                if( it->dag ){ diff_log_trace -= s[j]*interpolate_B_ij_t( t[j] - it->t, index );}
                else{ diff_log_trace += s[j]*interpolate_B_ij_t( t[j] - it->t, index );}
            }
            
        }
    }
    for( int i=0 ; i<static_cast<int>(dag_vec.size()) ; i++ )
    {
        map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(  make_pair(orb_vec[i], orb_vec[i]) );
        if( it_index !=  orb_to_index.end() ){ diff_log_trace -= -interpolate_B_ij_t(0.0, it_index->second) ; }
        for( int j=i+1; j<static_cast<int>(dag_vec.size()) ; j++ )
        {
            map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(   (orb_vec[i]<orb_vec[j])? make_pair(orb_vec[i],orb_vec[j]) : make_pair(orb_vec[j],orb_vec[i])    );
            if( it_index ==  orb_to_index.end() ){ continue;}
            else{ index = it_index->second; }
            diff_log_trace -= -s[i] * s[j] * interpolate_B_ij_t(t[j] - t[i], index);
        }
    }
    diff_log_trace = -2.0*diff_log_trace;
    
    return diff_log_trace;
}



double density_bath::log_trace_diff_replace( vector<local_config> &config_vec, vector<int>& orb_vec, vector<bool>& dag_vec , vector<double> & t_old, vector<double> &t_new)
{
    
    int index;
    double diff_log_trace = 0.0;
    vector<int> s(static_cast<int>(dag_vec.size()),0 );
    for( int i=0; i<static_cast<int>(dag_vec.size()) ; i++) { s[i] = dag_vec[i]?(1):(-1); }
    
    for( int j=0 ; j<static_cast<int>(dag_vec.size()) ; j++ )
    {
        for( int i=0;i<norb;i++ )
        {
            map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(   (i<orb_vec[j])? make_pair(i,orb_vec[j]) : make_pair(orb_vec[j],i)    );
            if( it_index ==  orb_to_index.end() ){ continue;}
            else{ index = it_index->second; }
            
            for( config_type::iterator it = config_vec[i].config.begin(); it!=config_vec[i].config.end(); it++)
            {
                if( it->dag ){ diff_log_trace += -s[j]*interpolate_B_ij_t( t_new[j] - it->t, index )   + s[j]*interpolate_B_ij_t( t_old[j] - it->t, index )  ;}
                else{ diff_log_trace += s[j]*interpolate_B_ij_t( t_new[j] - it->t, index )   - s[j]*interpolate_B_ij_t( t_old[j] - it->t, index ) ;}
            }
            
        }
    }
    for( int i=0 ; i<static_cast<int>(dag_vec.size()) ; i++ )
    {
        map< pair<int,int>, int>::iterator it_index =  orb_to_index.find( make_pair(orb_vec[i],orb_vec[i]) );
        if( it_index !=  orb_to_index.end() )
        {
            diff_log_trace += -interpolate_B_ij_t( 0.0 , it_index->second) + interpolate_B_ij_t(t_old[i]-t_new[i], it_index->second );
        }
        
        for( int j=i+1; j<static_cast<int>(dag_vec.size()) ; j++ )
        {
            map< pair<int,int>, int>::iterator it_index =  orb_to_index.find(   (orb_vec[i]<orb_vec[j])? make_pair(orb_vec[i],orb_vec[j]) : make_pair(orb_vec[j],orb_vec[i])    );
            if( it_index ==  orb_to_index.end() ){ continue;}
            else{ index = it_index->second; }
            diff_log_trace += s[i] * s[j] * interpolate_B_ij_t(t_new[j] - t_old[i], index) +s[i] * s[j] * interpolate_B_ij_t(t_new[i] - t_old[j], index);
            diff_log_trace += -s[i]*s[j] * interpolate_B_ij_t(t_new[i]-t_new[j], index);
            diff_log_trace += -s[i]*s[j]*interpolate_B_ij_t(t_old[i]-t_old[j], index);
        }
    }
    diff_log_trace = 2.0*diff_log_trace;
    return diff_log_trace;
    
}



bool density_bath::check_sanity(vector<local_config> &config_vec )
{
    
    double test_log_trace = 0.0;
    int sign;
    int index;
    for( int i=0; i<norb; i++)
    {
        for( int j=0; j<norb; j++)
        {
            for( config_type::iterator it1 = config_vec[i].config.begin(); it1 !=config_vec[i].config.end(); it1 ++)
            {
                for( config_type::iterator it2 = config_vec[j].config.begin(); it2 !=config_vec[j].config.end(); it2 ++)
                {
                    if( i==j && it1 == it2 ) continue;
                    
                    map< pair<int,int>, int>::iterator it_index =  orb_to_index.find( i<j ?  make_pair(i,j) :make_pair(j,i) );
                    if( it_index ==  orb_to_index.end() ){ continue;}
                    else{ index = it_index->second; }
                    
                    if( it1 ->dag == it2 ->dag ){ sign = -1.0;}
                    else{ sign = 1.0;}
                    test_log_trace +=  sign*interpolate_B_ij_t( it1->t - it2->t, index );
                    
                }
                
                
            }
            
            
        }
    }
    if( abs(test_log_trace - log_trace) >1e-8 )
    {
        cerr<<" density bosonic bath error "<<endl;
        cerr<<test_log_trace<<" "<<log_trace<<endl;
        return false;
    }
    return true;
    
}


