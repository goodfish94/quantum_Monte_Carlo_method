//
//  trace.cpp
//  cthyb
//
//  Created by 胡昊昱 on 7/16/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#include "trace.hpp"


double trace_cal::get_eve( )
{
    double eve = 0.0;
    for( int i=0; i<fock_space.size();i++ )
    {
        for( int j=i; j<fock_space.size(); j++)
        {
            if( i!=j)
                eve += (U_[i][j]+U_[j][i])*fock_space[i]*fock_space[j];
            else
                eve += U_[i][j]*fock_space[i]*fock_space[j];
        }
    }
    return eve;
}

double trace_cal::get_update_eve( int flip_orb)
{// trace_new = trace + return after flip flip_orb// doesn't flip the spin in this function
    double eve = 0.0;
    for( int i=0; i<fock_space.size();i++ )
    {
        if( i!=flip_orb)
            eve += (U_[i][flip_orb]+U_[flip_orb][i])*(1-2*fock_space[flip_orb])*fock_space[i];
        else
            eve += (U_[flip_orb][flip_orb])*(1-2*fock_space[flip_orb]);
        
    }
    return eve;
}
double trace_cal::log_trace_aux(  vector<local_config> &config_vec, int orb_index, op &node_old, op &node_new )
{
    // calculate the trace diff after insert node_old node_new segment works only if node_new->t > node_old_t and node_old.dag != node_new.dag (use log_trace to judge whether insert or remove) 
    // return log(trace diff)
    using namespace std;
    
    double delta_log_trace =0.0;
 
    for( int i=0;i<config_vec.size(); i++)
    {
        if( i!=orb_index)
        {
            config_it_vec[i] =config_vec[i].config.lower_bound(node_old);
            if(config_it_vec[i] == config_vec[i].config.end() )
            {
                if( config_vec[i].get_seg_num() == 0)
                {
                    if(config_vec[i].get_if_full() ){ fock_space[i] = 1;}
                    else { fock_space[i] = 0;}
                }
                else
                {
                    if( config_vec[i].get_if_wrap() ){ fock_space[i] = 1;}
                    else{ fock_space[i]=0; }
                }
            }
            else
            {
                if( config_it_vec[i]-> dag == true ){ fock_space[i] = 0;}
                else { fock_space[i] = 1;}
            }
        }
        else
        {
            fock_space[i] = 0;
        }
    }
  
    
    int sign_overlap = 1;
    if( node_old.dag  == false )
    {//remove
        sign_overlap = -1;
    }

    bool if_update = true;
    double time_ind_one = node_old.t;
    double time_ind_two =  node_old.t;
    
    double current_eve = get_update_eve(orb_index);
   
    
    int index_j=0;
    double dt;
    while(1)
    {
        if_update = false;
        time_ind_one = time_ind_two;
        time_ind_two =  node_new.t;
        for( int j=0;j<config_vec.size();j++)
        {
            
            if( (j!= orb_index)&& (config_it_vec[j] != config_vec[j].config.end() ) &&  (config_it_vec[j]->t < time_ind_two ) )
            {
                index_j = j ;
                time_ind_two =  config_it_vec[j]->t;
                if_update = true;
            }
            
        }
        dt =  time_ind_two - time_ind_one;
        delta_log_trace += current_eve*dt ;
        for( int j=0;j<config_vec.size();j++)
        {
            if( j==orb_index ){ delta_overlap[orb_index][j] +=  sign_overlap*dt;  }
            else { delta_overlap[orb_index][j] += sign_overlap*fock_space[j] * dt; }
        }
        
        if( !if_update) break;
        fock_space[index_j] = 1- fock_space[index_j];
        current_eve = get_update_eve(orb_index);

        if( config_it_vec[index_j] != config_vec[index_j].config.end() ) { config_it_vec[index_j] ++;}
        
    }
    
    if( node_old.dag  ==true )
    {//insert
        delta_log_trace = -delta_log_trace;
    }
    
    return delta_log_trace;
    
    
    
}


double trace_cal::log_trace_diff(   vector<local_config> &config_vec, int orb_index, op &node_old, op &node_new )
{
    // return log(trace diff) after insert [node_old, node_new segment] wrap if node_old.t>node_new.t
    using namespace std;
    double delta_log_trace=0.0;
    fill( update_orb.begin() , update_orb.end(), 0);
    update_orb[orb_index] = 1;
    for(int i=0;i<norb;i++)
        fill(delta_overlap[i].begin(), delta_overlap[i].end(),0.0);
    if( node_old.t<node_new.t)
    {
        delta_log_trace  = log_trace_aux( config_vec, orb_index,  node_old,  node_new );
    }
    else
    {
        
        op op_temp = node_new;
        op_temp.t=beta;
        delta_log_trace = log_trace_aux( config_vec, orb_index,  node_old, op_temp  );
  
        op_temp = node_old;
        op_temp.t = 0.0;
 
        delta_log_trace+=log_trace_aux( config_vec, orb_index,  op_temp,node_new );


    }
    
    return delta_log_trace;
    
    
    
}



double trace_cal::log_trace_diff(  vector<local_config> &config_vec, vector<int> &orb_vec, vector<bool> &dag_vec, double ts, double te )
{
    
    using namespace std;
    double delta_log_trace=0.0;
    
    
    fill( update_orb.begin() , update_orb.end(), 0);
    fill( update_dag.begin() , update_dag.end(), 0);
    for(int i=0;i<norb;i++)
        fill(delta_overlap[i].begin(), delta_overlap[i].end(),0.0);
    for( int i=0; i<static_cast<int>(orb_vec.size()); i++){ update_orb[orb_vec[i]]=1; update_dag[orb_vec[i]] = dag_vec[i]?(1):(-1); }
    
    
    if( ts<te)
    {
        delta_log_trace  = log_trace_aux( config_vec,orb_vec,dag_vec, ts, te);
    }
    else
    {
        delta_log_trace = log_trace_aux(  config_vec,orb_vec,dag_vec, ts, beta);
        delta_log_trace +=log_trace_aux( config_vec, orb_vec, dag_vec, 0,te);
    }
    
    return delta_log_trace;
    
    
    
}

double trace_cal::get_update_eve( )
{// trace_new = trace + return after flip flip_orb// doesn't flip the spin in this function
    double eve = 0.0;
    double old_eve = 0.0, new_eve = 0.0;
    for( int i=0; i<norb;i++ )
    {
        for( int j=0; j<norb ; j++)
        {
            old_eve= U_[i][j] *  (fock_space[i]*fock_space[j]);
            new_eve = U_[i][j] * ( (update_orb[i]==1) ? fock_space[i]:(1-fock_space[i]) ) * ( (update_orb[j]==1)? fock_space[j]:(1-fock_space[j]) );
            eve += new_eve-old_eve;
        }
        
    }
    return eve;
}

double trace_cal::log_trace_aux(  vector<local_config> &config_vec, vector<int> &orb_vec, vector<bool> &dag_vec, double ts, double te )
{
    
    
    double delta_log_trace =0.0;
    for( int i=0;i<config_vec.size(); i++)
    {
        if( update_orb[i] != 1 )
        {
            config_it_vec[i] =config_vec[i].config.lower_bound( op(ts, true,0,0));
            if(config_it_vec[i] == config_vec[i].config.end() )
            {
                if( config_vec[i].get_seg_num() == 0)
                {
                    if(config_vec[i].get_if_full() ){ fock_space[i] = 1;}
                    else { fock_space[i] = 0;}
                }
                else
                {
                    if( config_vec[i].get_if_wrap() ){ fock_space[i] = 1;}
                    else{ fock_space[i]=0; }
                }
            }
            else
            {
                if( config_it_vec[i]-> dag == true ){ fock_space[i] = 0;}
                else { fock_space[i] = 1;}
            }
        }
        else
        {
            fock_space[i] = (1+update_dag[i])/2; // =0 if insert ddag d, =1 is insert d ddag
        }
    }
    
    
    bool if_update = true;
    double time_ind_one = ts;
    double time_ind_two =  ts;
    
    double current_eve = get_update_eve();
    
    
    int index_j=0;
    double dt;
    while(1)
    {
        if_update = false;
        time_ind_one = time_ind_two;
        time_ind_two = te;
        for( int j=0;j<config_vec.size();j++)
        {
            
            if( (update_orb[j] != 1)&& (config_it_vec[j] != config_vec[j].config.end() ) &&  (config_it_vec[j]->t < time_ind_two ) )
            {
                index_j = j ;
                time_ind_two =  config_it_vec[j]->t;
                if_update = true;
            }
            
        }
        dt =  time_ind_two - time_ind_one;
        delta_log_trace += current_eve*dt ;
        for( int j=0;j<norb;j++)
        {
            if( update_orb[j] == 1 )
            {
                for( int i=0; i<norb ;i++)
                {
                    if( update_orb[i] != 1)
                    {
                        delta_overlap[i][j] += -update_dag[j] * fock_space[i] * dt;// update overlap btw non flip and flip orb
                    }
                }
            }
        }
        
        
        if( !if_update) break;
        fock_space[index_j] = 1- fock_space[index_j];
        current_eve = get_update_eve();
        if( config_it_vec[index_j] != config_vec[index_j].config.end() ) { config_it_vec[index_j] ++;}
        
    }
    for( int i=0; i<norb ; i++)
    {
        if( update_orb[i] != 1 ) continue;
        for( int j=i; j<norb;j++)
        {
            if( update_orb[j] != 1) continue;
            if( i==j ){ delta_overlap[i][i] += -update_dag[i] * (te-ts);}
            else{ delta_overlap[i][j] += ( -update_dag[i] - update_dag[j] ) * (te-ts);  }
                
        }
    }// update overlap btw flipped orb;
    return delta_log_trace;
    
}































void trace_cal::set_overlap( vector<local_config> &config_vec)
{
    
    for( int i=0; i<norb; i++){ fill(overlap[i].begin(), overlap[i].end(), 0.0);}
    
    for(int i=0; i <norb ;i++)
    {
        config_it_vec[i] = config_vec[i].config.begin();
        if( config_vec[i].n_ == 0)
        {
            if( config_vec[i].if_full == true)
            {
                fock_space[i]=1;
            }
            else
            {
                fock_space[i]=0;
            }
        }
        else
        {
            if( config_vec[i].if_wrap ){ fock_space[i] = 1;}
            else{ fock_space[i] = 0; }
        }
    }
    

    double t_st = 0.0;
    double t_en = 0.0;
    int index_j=0;
    double dt;
    bool if_update = true;
    while(1)
    {
        if_update = false;
        t_st = t_en;
        t_en = beta;
        
        for( int j=0;j<config_vec.size();j++)
        {
            
            if( (config_it_vec[j] != config_vec[j].config.end() ) &&  (config_it_vec[j]->t <= t_en ) )
            {
                index_j = j ;
                t_en =  config_it_vec[j]->t;
                if_update = true;
            }
        }
        dt = t_en-t_st;
        
        for( int i=0; i<norb; i++)
        {
            for( int j=i; j<norb; j++)
            {
                overlap[i][j] += dt*fock_space[i]*fock_space[j];
            }
        }
        if( !if_update) break;
        
        fock_space[index_j] = 1- fock_space[index_j];
        if( config_it_vec[index_j] != config_vec[index_j].config.end() ) { config_it_vec[index_j] ++;}
    }
    
}





double trace_cal::cal_log_trace( vector<local_config> &config_vec)
{
    using namespace std;

    double trace =0.0;
    for(int i=0; i <norb ;i++)
    {
        config_it_vec[i] = config_vec[i].config.begin();
        if( config_vec[i].n_ == 0)
        {
            if( config_vec[i].if_full == true)
            {
                fock_space[i]=1;
            }
            else
            {
                fock_space[i]=0;
            }
        }
        else
        {
            if( config_vec[i].if_wrap ){ fock_space[i] = 1;}
            else{ fock_space[i] = 0; }
        }
    }
    

    double eve=get_eve( );
    double t_st = 0.0;
    double t_en = 0.0;
    int index_j=0;
    
    bool if_update = true;
    while(1)
    {
        if_update = false;
        t_st = t_en;
        t_en = beta;
        
        for( int j=0;j<config_vec.size();j++)
        {
            
            if( (config_it_vec[j] != config_vec[j].config.end() ) &&  (config_it_vec[j]->t < t_en ) )
            {
                index_j = j ;
                t_en =  config_it_vec[j]->t;
                if_update = true;
            }
        }
        trace += eve* ( t_en - t_st ) ;
        if( !if_update) break;
        
        eve += get_update_eve( index_j);
        fock_space[index_j] = 1- fock_space[index_j];
        if( config_it_vec[index_j] != config_vec[index_j].config.end() ) { config_it_vec[index_j] ++;}
    }
    trace = -trace;// exp(-betaH)
    
    return trace;
}


double trace_cal::init_log_trace(  vector<local_config> &config_vec )
{
    double temp = cal_log_trace(config_vec);
    log_trace = temp;
    set_overlap(config_vec);
    return temp;
}// calculate log trace and let log_trace be it;



bool trace_cal::check_sanity(vector<local_config> &config_vec )
{
  
    double test_log_trace = 0.0;
    for( int i=0; i<norb;i++)
    {
        for( int j=0; j<norb;j++)
        {
            if( i==j ) { test_log_trace += overlap[i][i] * U_[i][i]; }
            else { test_log_trace += overlap[i][j] *( U_[i][j] + U_[j][i] );}
        }
    }
    test_log_trace = -test_log_trace;
    double test_log_trace_2 = cal_log_trace(config_vec);
    if( std::abs(test_log_trace - test_log_trace_2)>1e-8 || std::abs(test_log_trace_2 - log_trace) > 1e-8)
    {

        cerr<<" log trace error "<<endl;
        cerr<< test_log_trace<<" "<< test_log_trace_2<<" "<<log_trace<<endl;
        for( int i=0; i<norb;i++)
        {
            cout<<i<<"  "<<i<<"  "<<overlap[i][i]<<endl;
            for( int j=i+1; j<norb;j++)
            {
                cerr<<i<<"  "<<j<<"  "<<overlap[i][j] + overlap[j][i]<<endl;
            }
        }
        return false;
        
        
    }

    return true;
    
}




double trace_cal::log_trace_diff( vector<local_config> &config_vec, vector<int> &orb_vec, vector<bool> &dag_at_ts_vec, vector<double> ts, vector<double> te )
{// works only if all the orb index are diff
// try to insert multi sgements
    

    double delta_log_trace=0.0;
    
    for(int i=0;i<norb;i++)
        fill(delta_overlap[i].begin(), delta_overlap[i].end(),0.0);
    
    op node_old;
    op node_new;

    for( int i=0; i< static_cast<int>( orb_vec.size() ) ; i++  )
    {
        node_old.set_value( ts[i], dag_at_ts_vec[i],0,0);
        node_new.set_value( te[i], !dag_at_ts_vec[i], 0,0);
        fill( update_orb.begin() , update_orb.end(), 0);
        update_orb[orb_vec[i]] = 1;
        if( node_old.t<node_new.t)
        {
            delta_log_trace += log_trace_aux( config_vec, orb_vec[i],  node_old,  node_new );
        }
        else
        {
            
            op op_temp = node_new;
            op_temp.t=beta;
            delta_log_trace += log_trace_aux( config_vec, orb_vec[i],  node_old, op_temp  );
            
            op_temp = node_old;
            op_temp.t = 0.0;
            
            delta_log_trace +=log_trace_aux( config_vec, orb_vec[i],  op_temp,node_new );
            
        }
//        cout<<" !!!!!!!!!!!!!!!!!!!!!!! "<<endl;
//        cout<<"orb "<<orb_vec[i]<<" "<<ts[i]<<" "<<te[i]<<" "<<dag_at_ts_vec[i]<<endl;
//        for( int i=0; i< norb; i++)
//        {
//            cout<<i<<" "<<i<<" "<<delta_overlap[i][i]<<endl;
//            for( int j=i+1;j<norb;j++)
//            {
//                cout<<i<<" "<<j<<" "<<delta_overlap[i][j]+delta_overlap[j][i]<<endl;
//            }
//        }
//        cout<<" !!!!!!!!!!!!!!!!!!!!!!! "<<endl;
//        
        
    }
    double overlap_i_j;
    double ti_e_extend;
    double tj_e_extend;
    for( int i=0; i<static_cast<int>( orb_vec.size() ) ; i++)
    {
        
        if( te[i] < ts[i] ){ ti_e_extend = beta + te[i];}
        else{ ti_e_extend = te[i]; }
        
        for( int j=i+1 ; j<static_cast<int>(orb_vec.size()); j++)
        {
        
            
            if( te[j] < ts[j] ){ tj_e_extend = beta+ te[j];}
            else{ tj_e_extend = te[j] ;}
            overlap_i_j = max( 0.0, min( tj_e_extend , ti_e_extend) - max(ts[i], ts[j]) );
            if( dag_at_ts_vec[i] == dag_at_ts_vec[j] ){}
            else{ overlap_i_j = -overlap_i_j; }
            
            delta_overlap[orb_vec[i]][orb_vec[j]] += overlap_i_j;
            delta_log_trace += - ( U_[orb_vec[i]][orb_vec[j]] + U_[orb_vec[j]][orb_vec[i]] )* overlap_i_j;
           
            
        }
    }
    return delta_log_trace;
    


}
