//
//  sovler_fermi.cpp
//  cthyb
//
//  Created by 胡昊昱 on 7/29/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#include "solver.hpp"


bool solver::fermion_bath_update_insert()
{
    int k_orb = rand.random_int(norb);
    
    
    op_one.t = rand.random_double(beta);
    op_one.subtype = 0; op_one.type = 0;
    if( rand.random_double(1.0)<0.5){op_one.dag=true;}
    else{op_one.dag=false;}
    
    lmax = config_vec[k_orb].try_insert(op_one);
    if( lmax <= 0.0 ){ return false; }// directly failed
    
    op_two.t = op_one.t+ rand.random_double(lmax);
    if( op_two.t>beta){ op_two.t -= beta;}
    op_two.subtype = 0; op_two.type = 0;
    op_two.dag = !(op_one.dag);
    
    pair_one.first = k_orb;
    pair_two.first = k_orb;
    
    if( op_one.dag)
    {
        pair_one.second = op_one.t;
        pair_two.second = op_two.t;
    }
    else
    {
        pair_one.second = op_two.t;
        pair_two.second = op_one.t;
    }
    int k_block = orbital_to_block[k_orb];
    
    int sign_change=1;
    if( pair_one.second >pair_two.second )
    {
        sign_change = -sign_change;
    }
    
   log_trace_diff = trace_.log_trace_diff(config_vec, k_orb, op_one, op_two);
    
    
   log_trace_diff_B_bath = B_bath.log_trace_diff_insert(config_vec, op_one, op_two, k_orb);
    
    
    double p = exp(log_trace_diff + log_trace_diff_B_bath) ;
    
    
    p = p* ( (fermi_bath[k_block]).try_insert(pair_one,pair_two) );
    p = p* beta*lmax/( config_vec[k_orb].get_seg_num(0,0) +1);
    
    if( p<0.0 )
    {
        p=-p;
        sign_change = -sign_change;
    }
    
    
    
    
    if( p>=1 ||  rand.random_double(1.0) < p )
    {
        config_vec[k_orb].insert(op_one, op_two);
        trace_.update(log_trace_diff);
        fermi_bath[k_block].do_insert();
        B_bath.update_log_trace(log_trace_diff_B_bath);
        //        cout<<" insert "<<k_orb<<endl;
        //        cout<<" current sign, sign change "<<sign_<<" "<<sign_change<<endl;
        //        cout<<op_one.dag<<" "<<op_one.t<<endl;
        //        cout<<op_two.dag<<" "<<op_two.t<<endl;
        
        sign_ *= sign_change;
        return true;
    }
    else
    {
        config_vec[k_orb].lock_insert();
        fermi_bath[k_block].lock_insert();
        return false;
        
    }
    

    return true;
    
}

bool solver::fermion_bath_update_remove()
{
    int k_orb = rand.random_int(norb);
    int n_type = config_vec[k_orb].get_seg_num(0,0);
    int k_th_seg = rand.random_int(n_type);
    
    if( 0.5< rand.random_double(1.0))
    {
        lmax = config_vec[k_orb].try_remove(0, 0, k_th_seg ,true, op_one, op_two  );
    }
    else
    {
        lmax = config_vec[k_orb].try_remove(0, 0, k_th_seg ,false, op_one, op_two  );
    }
    if( lmax<=0.0){ return false; }// directly failed
    
    int k_block = orbital_to_block[k_orb];
    int sign_change=1;
    
    op_one.dag = !op_one.dag;
    op_two.dag = !op_two.dag; // remove seg == insert anti seg
    log_trace_diff = trace_.log_trace_diff( config_vec, k_orb, op_one, op_two);
    op_one.dag = !op_one.dag;
    op_two.dag = !op_two.dag; // flip back
    
    log_trace_diff_B_bath = B_bath.log_trace_diff_remove(config_vec, op_one, op_two, k_orb);
    
    pair_one.first = k_orb;
    pair_two.first = k_orb;
    if( op_one.dag)
    {
        pair_one.second = op_one.t;
        pair_two.second = op_two.t;
    }
    else
    {
        pair_one.second = op_two.t;
        pair_two.second = op_one.t;
    }
    
    if( pair_one.second > pair_two.second )
    {
        sign_change = -sign_change;
    }
    
    double p = exp(log_trace_diff + log_trace_diff_B_bath);
    p = p*fermi_bath[k_block].try_remove(pair_one, pair_two);
    p = p*n_type/(beta*lmax);
    
    if( p<0.0)
    {
        sign_change = -sign_change;
        p=-p;
    }
    
    if( p>=1 ||  rand.random_double(1.0) < p )
    {
        
        //        cout<<" remove "<<k_orb<<" "<<fermi_bath[k_orb].get_size()<<endl;
        //        cout<<"row and col:"<<fermi_bath[k_orb].row_remove<<" "<<fermi_bath[k_orb].col_remove<<endl;
        //        cout<<" sign "<<sign_<<" "<<sign_change<<endl;
        config_vec[k_orb].remove();
        trace_.update(log_trace_diff);
        B_bath.update_log_trace(log_trace_diff_B_bath);
        sign_ *= sign_change;
        fermi_bath[k_block].do_remove(sign_change);// swap(rol/col) sign
        sign_ *= sign_change;
        
        
        //        cout<<op_one.dag<<" "<<op_one.t<<endl;
        //        cout<<op_two.dag<<" "<<op_two.t<<endl;
        
        return true;
    }
    else
    {
        config_vec[k_orb].lock_remove();
        fermi_bath[k_block].lock_remove();
        return false;
        
    }
    

    return true;
    
}

bool solver::fermion_bath_update_shift()
{
    int k_orb = rand.random_int(norb);
    int n_type = config_vec[k_orb].get_seg_num(0,0);
    int k_th_seg = rand.random_int(n_type);
    bool if_dag = true;
    if( rand.random_double(1.0) < 0.5){ if_dag =false;}
    double down_lim =0.0, up_lim = 0.0;
    
    if( !( config_vec[k_orb].try_shift(0,0, k_th_seg, if_dag, down_lim, up_lim, op_one) ) ){return false;}
    op_two = op_one;
    
    double dt = down_lim + rand.random_double(up_lim-down_lim);

    op_two.t += dt;
    int sign_change = 1;
    if( op_two.t<0){ op_two.t+=beta; sign_change = -1;}
    else if(op_two.t>beta){op_two.t-=beta;sign_change = -1;}
    op_one.dag = !(op_two.dag);
    
    if( dt>0.0 )
    {
        log_trace_diff = trace_.log_trace_diff(config_vec, k_orb, op_one, op_two);
    }
    else
    {
        log_trace_diff = trace_.log_trace_diff( config_vec, k_orb,  op_two,op_one);
    }
    
    op_one.dag = op_two.dag;
    log_trace_diff_B_bath = B_bath.log_trace_diff_replace(config_vec, op_one, op_two, k_orb);
    
    
    
    double p = exp(log_trace_diff + log_trace_diff_B_bath );
    int k_block = orbital_to_block[k_orb];
    pair_one.first = k_orb;
    pair_two.first = k_orb;
    pair_one.second = op_one.t;
    pair_two.second = op_two.t;
    
    p *= fermi_bath[k_block].try_replace(pair_one, pair_two, op_one.dag);
    
    if(p<0.0)
    {
        sign_change = -sign_change;
        p=-p;
    }
    
    if( p>=1 ||  rand.random_double(1.0) < p )
    {
        
        config_vec[k_orb].shift(op_two);
        trace_.update(log_trace_diff);
        B_bath.update_log_trace( log_trace_diff_B_bath );
        fermi_bath[k_block].do_replace();
        sign_ *= sign_change;
        //        cout<<"sign update "<<swap_sign<<" "<<sign_change<<endl;
        //        cout<<op_one.dag<<" "<<op_one.t<<endl;
        //        cout<<op_two.dag<<" "<<op_two.t<<endl;
        return true;
    }
    else
    {
        config_vec[k_orb].lock_shift();
        fermi_bath[k_block].lock_replace();
        return false;
        
    }

    return true;
}



void solver::fermion_bath_update( int n_sweep)
{
    double rand_num ;
    for(int i=0; i<n_sweep ;i++)
    {
//        cout<<i<<"-th--------------"<<endl;
        rand_num = rand.random_double(1.0);
        
        if( rand_num<0.3)
        {
//            cout<<" insert fm--------------"<<endl;
            propose_fermi_bath[0]++;
            if(fermion_bath_update_insert()){ accept_fermi_bath[0]++;}
        }
        else if( rand_num<0.6)
        {
//            cout<<" remove fm--------------"<<endl;
            propose_fermi_bath[1]++;
            if(fermion_bath_update_remove()){ accept_fermi_bath[1]++;}
            
        }
        else
        {
            
            propose_fermi_bath[2]++;
            if(fermion_bath_update_shift()){ accept_fermi_bath[2]++;}

           

        }

//        if( !( check_sanity()) )
//        {
//            exit(10);
//        }
        
        //        cout<<"current sign = "<<permute_sign<<endl;
        //        cout<<"config one "<<endl;
        //        config_vec[0].print_config();
        //        cout<<"config two "<<endl;
        //        config_vec[1].print_config();
        //
        //        for( int i=0; i<norb; i++)
        //        {
        //            for( int j=0; j<norb ; j++)
        //            {
        //                cout<<i<<" "<<j<<" "<<trace_.get_overlap()[i][j]<<endl;
        //            }
        //        }
        
        //
        //        if( !( check_sanity()) )
        //        {
        //            cout<<" not sanity "<<endl;
        //            cout<<"config one "<<endl;
        //            config_vec[0].print_config();
        //            cout<<"config two "<<endl;
        //            config_vec[1].print_config();
        //
        //            cout<<" det "<<endl;
        //            cout<< fermi_bath[0].get_det()<<" "<<fermi_bath[1].get_det()<<endl;
        //            cout<<"time start array"<<endl;
        //            for( auto it = fermi_bath[1].time_start.begin() ; it != fermi_bath[1].time_start.end() ; it++)
        //            {
        //                cout<<(it->first).second<<" "<<(it->second)<<" "<<endl;
        //            }
        //            cout<<endl;
        //            cout<<"time end array"<<endl;
        //            for( auto it = fermi_bath[1].time_end.begin() ; it != fermi_bath[1].time_end.end() ; it++)
        //            {
        //                cout<<(it->first).second<<" "<<(it->second)<<" "<<endl;
        //            }
        //            exit(10);
        //        }
        //        if( get_sign() != 1)
        //        {
        //
        //            cout<<"config one "<<endl;
        //            config_vec[0].print_config();
        //            cout<<"config two "<<endl;
        //            config_vec[1].print_config();
        //
        //            cout<<" det "<<endl;
        //            cout<< fermi_bath[0].get_det()<<" "<<fermi_bath[1].get_det()<<endl;
        //            cout<< fermi_bath[1].N_<<" "<<fermi_bath[1].N_<<endl;
        //            cout<<"time start array"<<endl;
        //            for( auto it = fermi_bath[1].time_start.begin() ; it != fermi_bath[1].time_start.end() ; it++)
        //            {
        //                cout<<(it->first).second<<" "<<(it->second)<<" "<<endl;
        //            }
        //            cout<<endl;
        //            cout<<"time end array"<<endl;
        //            for( auto it = fermi_bath[1].time_end.begin() ; it != fermi_bath[1].time_end.end() ; it++)
        //            {
        //                cout<<(it->first).second<<" "<<(it->second)<<" "<<endl;
        //            }
        //            cout<<endl;
        //
        //            exit(1);
        //
        //
        //        }
        
    }
    
}

