//
//  solver.cpp
//  cthyb
//
//  Created by 胡昊昱 on 7/15/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#include "solver.hpp"



int solver::obatin_i_meas()
{// set i_meas according to the current information
    double n_seg=0.0;
    double accept = 0.3;
    
    n_seg += accu_.avg_fermion_bath_num() + accu_.avg_diag_flip_bath_num();

    return floor(n_seg/accept);

}


void solver::print()
{
    accu_.print_all();
    ofstream file;
    file.open("qmc.txt", ios::trunc);

    if( ! file.is_open() ){ cerr<<"can't open qmc.txt"<<endl; return ; }
    
    file<<"upate    prop    accep    ratio"<<endl;
    file<<" fermion bath "<<endl;
    file<<"insert    "<<propose_fermi_bath[0]<<"    "<<accept_fermi_bath[0]<<"    "<<double(accept_fermi_bath[0])/propose_fermi_bath[0]<<endl;
    file<<"remove    "<<propose_fermi_bath[1]<<"    "<<accept_fermi_bath[1]<<"    "<<double(accept_fermi_bath[1])/propose_fermi_bath[1]<<endl;
    file<<"shift     "<<propose_fermi_bath[2]<<"    "<<accept_fermi_bath[2]<<"    "<<double(accept_fermi_bath[2])/propose_fermi_bath[2]<<endl;
    
    file<<" diag flip bath "<<endl;
    file<<"insert    "<<propose_diag_flip_bath[0]<<"    "<<accept_diag_flip_bath[0]<<"    "<<double(accept_diag_flip_bath[0])/propose_diag_flip_bath[0]<<endl;
    file<<"remove    "<<propose_diag_flip_bath[1]<<"    "<<accept_diag_flip_bath[1]<<"    "<<double(accept_diag_flip_bath[1])/propose_diag_flip_bath[1]<<endl;
    file<<"shift      "<<propose_diag_flip_bath[2]<<"    "<<accept_diag_flip_bath[2]<<"    "<<double(accept_diag_flip_bath[2])/propose_diag_flip_bath[2]<<endl;
     file<<"swap      "<<propose_diag_flip_bath[3]<<"    "<<accept_diag_flip_bath[3]<<"    "<<double(accept_diag_flip_bath[3])/propose_diag_flip_bath[3]<<endl;
    
    file<<"update btw diag flip bath and fermion bath "<<endl;
    file<<"fermi to diag flip    "<<propose_fermion_diag_flip_bath[0]<<"    "<<accept_fermion_diag_flip_bath[0]<<"    "<<double(accept_fermion_diag_flip_bath[0])/propose_fermion_diag_flip_bath[0]<<endl;
    file<<"diag flip to fermi    "<<propose_fermion_diag_flip_bath[1]<<"    "<<accept_fermion_diag_flip_bath[1]<<"    "<<double(accept_fermion_diag_flip_bath[1])/propose_fermion_diag_flip_bath[1]<<endl;
    file<<"swap                  "<<propose_fermion_diag_flip_bath[2]<<"    "<<accept_fermion_diag_flip_bath[2]<<"    "<<double(accept_fermion_diag_flip_bath[2])/propose_fermion_diag_flip_bath[2]<<endl;
    
    file<<"off diag flip bath "<<endl;
    file<<"insert    "<<propose_off_diag_flip_bath[0]<<"    "<<accept_off_diag_flip_bath[0]<<"    "<<double(accept_off_diag_flip_bath[0])/propose_off_diag_flip_bath[0]<<endl;
    file<<"remove    "<<propose_off_diag_flip_bath[1]<<"    "<<accept_off_diag_flip_bath[1]<<"    "<<double(accept_off_diag_flip_bath[1])/propose_off_diag_flip_bath[1]<<endl;
    file<<"shift      "<<propose_off_diag_flip_bath[2]<<"    "<<accept_off_diag_flip_bath[2]<<"    "<<double(accept_off_diag_flip_bath[2])/propose_off_diag_flip_bath[2]<<endl;
    file<<"swap      "<<propose_off_diag_flip_bath[3]<<"    "<<accept_off_diag_flip_bath[3]<<"    "<<double(accept_off_diag_flip_bath[3])/propose_off_diag_flip_bath[3]<<endl;
    
    file<<"update btw off diag flip bath and fermion bath "<<endl;
    file<<"fermi to off diag flip    "<<propose_fermion_off_diag_flip_bath[0]<<"    "<<accept_fermion_off_diag_flip_bath[0]<<"    "<<double(accept_fermion_off_diag_flip_bath[0])/propose_fermion_off_diag_flip_bath[0]<<endl;
    file<<"off diag flip to fermi    "<<propose_fermion_off_diag_flip_bath[1]<<"    "<<accept_fermion_off_diag_flip_bath[1]<<"    "<<double(accept_fermion_off_diag_flip_bath[1])/propose_fermion_off_diag_flip_bath[1]<<endl;
    file<<"swap                  "<<propose_fermion_off_diag_flip_bath[2]<<"    "<<accept_fermion_off_diag_flip_bath[2]<<"    "<<double(accept_fermion_off_diag_flip_bath[2])/propose_fermion_off_diag_flip_bath[2]<<endl;
   
    file.close();

}




bool solver::check_sanity( )
{
    bool if_sanity = true;

    for( vector<fermi_bath_block>::iterator it = fermi_bath.begin(); it!= fermi_bath.end();it++)
    {
        bool fermi_bath_sanity = it->check_sanity();
        if( !fermi_bath_sanity)
        {
            if_sanity=false;
            cerr<<"fermi bath error"<<endl;

            cerr<<static_cast<int>(it-fermi_bath.begin()) <<"-th bath"<<endl;
        }

    }
    for( vector<local_config>::iterator it = config_vec.begin(); it != config_vec.end(); it++)
    {
        bool local_config_sanity = it->check_sanity();
        if(!local_config_sanity)
        {
            if_sanity=false;
            cerr<<"local config error"<<endl;
            cerr<<static_cast<int>(it-config_vec.begin()) <<"-th orbital"<<endl;
        }
    }

    if( !trace_.check_sanity(config_vec) )
    {
        if_sanity=false;
        cerr<<" trace error "<<endl;
    }
    
    for( int i=0; i<norb; i++)
    {
        if( std::abs(trace_.get_overlap(i,i)/beta - config_vec[i].get_n_avg() ) >1e-7 )
        {
            if_sanity=false;
            cerr<<" overlap error "<<i<<" "<<trace_.get_overlap(i,i)/beta<<" "<<config_vec[i].get_n_avg() <<endl;
        }
    }

    if( !if_sanity)
    {
         for(vector<local_config>::iterator it=config_vec.begin(); it!=config_vec.end() ; it++)
             it->print_config();    }


    int n_node =0;
    int n_fm=0;

    for( vector<fermi_bath_block>::iterator it = fermi_bath.begin(); it!= fermi_bath.end() ; it++)
    {
        n_fm += (it->get_size());
    }


    for( int i =0 ; i< static_cast<int>(config_vec.size()); i++)
    {
        config_type test_config = config_vec[i].get_config();
        int block_index = orbital_to_block[ i ];
        fermi_bath_block test_bath = fermi_bath[block_index];
        map< pair<int,double> ,int>  start_ = test_bath.get_time_start();
        map< pair<int,double> ,int>  end_ = test_bath.get_time_end();



        for( config_type::iterator it = test_config.begin(); it!= test_config.end(); it++ )
        {
            op test_op = *(it);

            if( test_op.type !=0 )
                continue;
            else
            {
                n_node ++;
                if( test_op.dag == true )
                {
                    map< pair<int,double> ,int>::iterator tmp_it = start_.find( make_pair(i,test_op.t) );
                    if( tmp_it == start_.end())
                    {
                        cerr<<" can't find op in fermion bath"<<endl;
                        if_sanity = false;
                    }
                }
                else
                {
                    map< pair<int,double> ,int>::iterator tmp_it = end_.find( make_pair(i,test_op.t) );
                    if( tmp_it == end_.end())
                    {
                        cerr<<" can't find op in fermion bath"<<endl;
                        if_sanity = false;
                    }

                }


            }
        }

    }

    if( n_fm*2 != n_node)
    {
        cerr<<" fermi bath seg number != config seg num "<<n_fm*2<<" "<<n_node<<endl;
        if_sanity=false;
    }

    int test_sign = calculate_perm_sign();
    for( vector<fermi_bath_block>::iterator it= fermi_bath.begin(); it != fermi_bath.end() ; it++)
    {
        if( it->get_det() <0.0)
        {
            test_sign = - test_sign;
        }
    }
    if( test_sign !=  sign_ )
    {
        cerr<<" perm sign error "<<test_sign<<" "<<sign_<<endl;
        if_sanity = false;

    }
    
    if( B_bath.check_sanity(config_vec) == false)
    {
        cerr<< " B bath errror "<<endl;
        if_sanity = false;
    }
    
    if( diag_J_bath.check_sanity(config_vec) == false)
    {
        cerr<< " daig J bath error "<<endl;
        if_sanity = false;
    }
    
    if( off_diag_J_bath.check_sanity(config_vec) == false)
    {
        cerr<< "off  daig J bath error "<<endl;
        if_sanity = false;
    }
//    
//    if( ! if_sanity )
//    {
//        for( int i=0; i<norb; i++)
//        {
//            cerr<<" config -"<<i<<endl;
//            config_vec[i].print_config();
//        }
//    }
    
    return if_sanity;


}


int solver::calculate_perm_sign()
{// calculate fermionic perm sign;
    
    int perm_sign = 1;
    
    
    for( vector<fermi_bath_block>::iterator it = fermi_bath.begin(); it != fermi_bath.end(); it++ )
    {
        
        vector<double> tmp(2*(it->nsize) , 0.0 );
        
        for( map< pair<int,double> ,int>::iterator it2 = (it->time_start).begin(); it2 != (it->time_start).end() ; it2++ )
        {
            tmp[ 2*it2->second] = (it2->first).second;
            
        }
        for( map< pair<int,double> ,int>::iterator it2 = (it->time_end).begin(); it2 != (it->time_end).end() ; it2++ )
        {
            tmp[ 2*it2->second +1] = (it2->first).second;
            
        }
        
       
        for( int i=0; i< 2*(it->nsize); i++)
        {
            for( int j=0; j<2*(it->nsize)-i-1; j++ )
            {
                if( tmp[j]>tmp[j+1])
                {
                    swap( tmp[j], tmp[j+1]);
                    perm_sign = - perm_sign;
                }
            }
        }
        
    }
    

    
    return perm_sign;
    
    

    
    
    
}



void solver::mixed_update( int n_sweep )
{
    double rand_num;
    
    for( int i=0; i<n_sweep ;i++)
    {
        rand_num = rand.random_double(1.0);
        if( rand_num<0.4 )
        {
            
            fermion_bath_update(1);
//            cout<<"fm"<<endl;
//            cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        }
        else if(rand_num <0.6)
        {
            
            diag_flip_bath_update(1);
//            cout<<"diag "<<endl;
//            cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        }
        else
        {
//            if( !check_sanity()){ cerr<<" not sanity "<<endl; exit(2);}
//            cout<<"sanity"<<endl;
//            cout<<"diag J config"<<endl;
//            diag_J_bath.print_config();
//            cout<<endl<<endl;
            fermion_bath_diag_flip_updae(1);
//            cout<<"fm diag"<<endl;
//            cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
//            cout<<endl<<endl;
        }
        
//        if( !check_sanity()){ cerr<<" not sanity "<<endl; exit(2);}
        
    }
    
    
    
}
