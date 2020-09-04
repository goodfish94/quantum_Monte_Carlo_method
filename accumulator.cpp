//
//  accumulator.cpp
//  cthyb
//
//  Created by 胡昊昱 on 7/16/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#include "accumulator.hpp"

//static double test_n = 0.0;
void accumulator::set_meas_off_diag_flip_tau(vector< vector<int> > index_to_orb_off_diag_flip_in  )
{
    if_meas_off_diag_flip_w = false;
    if_meas_off_diag_flip_tau = true;
    index_to_orb_off_diag_flip = index_to_orb_off_diag_flip_in;
    sign_sum_off_diag_flip = 0;
    n_meas_off_diag_flip = 0;
    
    G_two_off_diag_flip.resize(static_cast<int>(index_to_orb_off_diag_flip.size()), vector<double>(N_t+1,0.0));
    G_two_off_diag_flip_mpi.resize(static_cast<int>(index_to_orb_off_diag_flip.size()), vector<double>(N_t+1,0.0));
    G_two_off_diag_flip_mpi_error.resize(static_cast<int>(index_to_orb_off_diag_flip.size()), vector<double>(N_t+1,0.0));
    
    histogram_off_diag_flip_bath.resize(index_to_orb_off_diag_flip_in.size(), vector<int>( max_off_diag_flip,0));
    
    
}
void accumulator::set_meas_off_diag_flip_w(vector< vector<int> > index_to_orb_off_diag_flip_in  )
{
    
    if_meas_off_diag_flip_w = true;
    if_meas_off_diag_flip_tau = false;
    index_to_orb_off_diag_flip = index_to_orb_off_diag_flip_in;
    sign_sum_off_diag_flip = 0;
    n_meas_off_diag_flip = 0;
    
    G_two_off_diag_flip.resize(static_cast<int>(index_to_orb_off_diag_flip.size()), vector<double>(N_w,0.0));
    G_two_off_diag_flip_mpi.resize(static_cast<int>(index_to_orb_off_diag_flip.size()), vector<double>(N_w,0.0));
    G_two_off_diag_flip_mpi_error.resize(static_cast<int>(index_to_orb_off_diag_flip.size()), vector<double>(N_w,0.0));
    histogram_off_diag_flip_bath.resize(index_to_orb_off_diag_flip_in.size(), vector<int>( max_off_diag_flip,0));
    
}

void accumulator::set_meas_diag_flip_tau( map<pair<int,int>,int> orb_to_index )
{
    if_meas_diag_flip_w = false;
    if_meas_diag_flip_tau = true;
    diag_flip_index = orb_to_index;
    sign_sum_diag_flip = 0;
    n_meas_diag_flip = 0;
    
    G_two_diag_flip.resize(static_cast<int>(orb_to_index.size()), vector<double>(N_t+1,0.0));
    G_two_diag_flip_mpi.resize(static_cast<int>(orb_to_index.size()), vector<double>(N_t+1,0.0));
    G_two_diag_flip_mpi_error.resize(static_cast<int>(orb_to_index.size()), vector<double>(N_t+1,0.0));
    
}
void accumulator::set_meas_diag_flip_w( map<pair<int,int>,int> orb_to_index )
{
    if_meas_diag_flip_w = true;
    if_meas_diag_flip_tau = false;
    diag_flip_index = orb_to_index;
    sign_sum_diag_flip = 0;
    n_meas_diag_flip = 0;
    
    G_two_diag_flip.resize(static_cast<int>(orb_to_index.size()), vector<double>(N_w,0.0));
    G_two_diag_flip_mpi.resize(static_cast<int>(orb_to_index.size()), vector<double>(N_w,0.0));
    G_two_diag_flip_mpi_error.resize(static_cast<int>(orb_to_index.size()), vector<double>(N_w,0.0));
    
}
void accumulator::set_meas_one_tau( )
{
    if_meas_one_tau = true;
    sign_sum_G_one =0;
    n_meas_G_one = 0;
    
    int num_of_nonzero_g_one=0;
    
    
    for( vector< set<int> >::iterator it=vec_index_field.begin(); it!=vec_index_field.end(); it++)
    {
        int size_block = static_cast<int>( it->size() );
        num_of_nonzero_g_one += size_block*size_block;
        
    }
    G_one.resize( num_of_nonzero_g_one, vector<fermi_type>( N_t+1,0.0));
    G_one_mpi.resize( num_of_nonzero_g_one, vector<fermi_type>( N_t+1,0.0));
    G_one_mpi_error.resize( num_of_nonzero_g_one, vector<fermi_type>( N_t+1,0.0));
    
    G_one_index.clear();
    int index = 0 ;
    for( vector< set<int> >::iterator it=vec_index_field.begin(); it!=vec_index_field.end(); it++)
    {
        for( set<int>::iterator it1 = it->begin(); it1 != it->end() ; it1++)
        {
            for( set<int>::iterator it2 = it->begin(); it2 != it->end() ; it2++)
            {
                pair<int,int> ij_(*it1, *it2);
    
                G_one_index.insert( make_pair(ij_,index) );
                index++;
                
            }
        }
    }
  
    
}

void accumulator::set_meas_two_tau(  vector< vector<int> > G_two_index_in )
{// only works for the case with S^+ S^-
    
    if_meas_two_tau = true;
    
    G_two_orb_index.clear();
    G_two_block_index.clear();
    
    int size = 0;
    for( int i=0; i< static_cast<int>(G_two_index_in.size()); i++)
    {
        vector<int> block(4,-1);
        for( int j=0; j<4;j++)
        {
            block[j] = orb_to_block[G_two_index_in[i][j]];
        }
        
        if( block[1] != block[0] && G_two_index_in[i][0] == G_two_index_in[i][3] &&  G_two_index_in[i][1] == G_two_index_in[i][2]  )
        {
            size ++;
            G_two_orb_index.push_back(G_two_index_in[i]);
            G_two_block_index.push_back(block);
        }
    }
    G_two.resize(size,vector<fermi_type>(N_t+1,0.0));
    G_two_mpi.resize(size,vector<fermi_type>(N_t+1,0.0));
    G_two_mpi_error.resize(size,vector<fermi_type>(N_t+1,0.0));
    
    n_meas_G_two=0;
    sign_sum_G_two=0;
}


void accumulator::set_meas_n_n_tau(  )
{
    if_meas_nn_tau = true;
    sign_sum_n_n_ = 0;
    n_meas_n_n_ = 0;
    
    orb_to_index_nn.clear();
    
    int index = 0;
    for( int i=0;i<norb;i++)
    {
        for( int j=i; j<norb; j++)
        {
            orb_to_index_nn.insert(make_pair(pair<int,int>(i,j), index));
            index++;
        }
    }
    n_n_.resize(orb_to_index_nn.size(), vector<double>(N_t+1,0.0) );
    n_n_mpi.resize(orb_to_index_nn.size(), vector<double>(N_t+1,0.0) );
    n_n_mpi_error.resize(orb_to_index_nn.size(), vector<double>(N_t+1,0.0) );
    n_i_t.resize(norb, vector<int>(N_t+1,0) );
    
}

void accumulator::measure_n_n(vector<local_config> &config_vec, int sign)
{
    int n_pt = 0;
    int st = 0;
    int en = 0;
    int index;
    int t;
    double double_N_t = static_cast<double>(N_t);
    
    
    n_meas_n_n_ += 1;
    sign_sum_n_n_ += sign;
    for( int i=0 ;i<norb; i++)
    {
        if( config_vec[i].get_if_full()){ fill(n_i_t[i].begin(), n_i_t[i].end(),  1); continue; }
        if( config_vec[i].get_if_wrap() && config_vec[i].get_op_num()!=0 ){ n_pt = 1;}
        else{ n_pt = 0; }
        st = 0;
        
        for( config_type::iterator it = config_vec[i].get_config().begin(); it !=(config_vec[i].get_config()).end(); it++)
        {
            en = floor( it->t/dt );
            if( en < st )
            {
                n_pt = 1- n_pt;
                continue;
            }
            fill( n_i_t[i].begin()+st, n_i_t[i].begin()+en+1, n_pt );
            n_pt = 1- n_pt;
            st = en+1;
        }
        fill( n_i_t[i].begin()+st, n_i_t[i].end(), n_pt);
    
    }
   
//    test_n +=n_i_t[0][0];
//    if( n_meas_n_n_ %1000 ==0)
//    {
//        cout<<test_n/n_meas_n_n_<<endl;
//    }

    
    for( int i=0; i<norb; i++)
    {
        
        for( int t_i =0 ;t_i <N_t ; t_i++ )
        {
            if( n_i_t[i][t_i] == 0){ continue; }
            for( int j=i; j<norb; j++)
            {
                index = orb_to_index_nn[ pair<int,int>(i,j)  ];
                for( int t_j = 0; t_j<N_t; t_j++ )
                {
                    if( n_i_t[j][t_j]==0){continue;}
                    t = t_i-t_j;
                    if( t< 0){ t += N_t; }
                    n_n_[index][t] += sign/double_N_t;

                    
                }
            
            }
        }
            
    }
        
        
    
}


    
    



void accumulator::measure_static( vector<local_config> &config_vec , trace_cal& trace_, diag_flip_bath& diag_J_bath, off_diag_flip_bath& off_diag_J_bath,  int sign )
{
    
    n_meas_static += 1;
    sign_sum_static += sign;
    for( int i=0; i<norb ;i++)
    {
        histogram_fermion_bath[i][ (config_vec[i].get_seg_num(0,0)<max_fermi) ? (config_vec[i].get_seg_num(0,0)):(max_fermi-1) ]++;
    }
    for( int k=0; k< static_cast<int>(diag_flip_index.size() ); k++)
    {
        histogram_diag_flip_bath[k][ diag_J_bath.get_nsize(k) < max_diag_flip ?  ( diag_J_bath.get_nsize(k) ): ( max_diag_flip ) ]++;
    }
    
    for( int k=0; k< static_cast<int>( index_to_orb_off_diag_flip.size() ); k++)
    {
        histogram_off_diag_flip_bath[k][ off_diag_J_bath.get_nsize(k) < max_off_diag_flip ?  ( off_diag_J_bath.get_nsize(k) ): ( max_off_diag_flip ) ]++;
    }
    for( int i=0; i<norb; i++)
    {
        for( int j=0; j<norb; j++)
        {
            n_i_n_j[i][j] += sign*( trace_.get_overlap(i , j) )/beta;
            
        }
    }
}




void accumulator::measure_g_one( vector<fermi_bath_block>& fermi_bath , int sign)
{
    if( !if_meas_one_tau ) return;
    int alpha,gamma;
    double t_s,t_e;
    int i,j;
    double re;
    double t;
    sign_sum_G_one += sign;
    n_meas_G_one += 1;
    
    for( vector< fermi_bath_block>::iterator it = fermi_bath.begin(); it!= fermi_bath.end(); it++)
    {
        
        for( map< pair<int,double> ,int>::iterator it_st = (it->time_start).begin(); it_st != (it->time_start).end() ;it_st++)
        {
            alpha = (it_st->first).first;
            t_s = (it_st->first).second;
            i= it_st -> second;
            
            for(  map< pair<int,double> ,int>::iterator it_en = (it->time_end).begin() ; it_en != (it->time_end).end(); it_en++)
            {
                gamma = (it_en->first).first;
                t_e = (it_en->first).second;
                j = it_en -> second;
                re = -sign * it->N_(j,i)/(beta*dt);
                
                t = t_s-t_e;
                if( t<0 ){ t += beta; re = -re;}
                G_one[ G_one_index[ make_pair(alpha, gamma)] ][ floor( t/dt )] += re;
                
            }
            
        }
        
        
    }
    
    
}

void accumulator::measure_diag_flip_bath( diag_flip_bath& diag_J_bath, int sign)
{
    double t;
    sign_sum_diag_flip += sign;
    n_meas_diag_flip += 1;
    if( if_meas_diag_flip_tau )
    {// measure tau
        for( int i=0; i< static_cast<int>( diag_flip_index.size()); i++)
        {
            for( list<Ldag_L>::iterator it = diag_J_bath.Ldag_L_config[i].begin(); it != diag_J_bath.Ldag_L_config[i].end(); it++)
            {
                t = ( it->te-it->ts );
                if( t<0 ){ t += beta;}
                G_two_diag_flip[i][floor( t/dt )] += sign/(beta*dt*it->jt);
                
            }
        
        }
    }
    else if(if_meas_diag_flip_w)
    {// measure w
        
        for( int i=0; i< static_cast<int>( diag_flip_index.size()); i++)
        {
            for( int n=0; n<N_w; n++)
            {
                for( list<Ldag_L>::iterator it = diag_J_bath.Ldag_L_config[i].begin(); it != diag_J_bath.Ldag_L_config[i].end(); it++)
                {
                    t = ( it->te-it->ts );
                    
                    G_two_diag_flip[i][n] += sign/(beta*it->jt)*cos(bose_nu_n[n]*t);
                    
                }
            }
        }
        
    }
}

void accumulator::measure_off_diag_flip_bath( off_diag_flip_bath &off_diag_J_bath,int sign)
{
    
    double t;
    sign_sum_off_diag_flip += sign;
    n_meas_off_diag_flip += 1;
    if( if_meas_diag_flip_tau )
    {// measure tau
        for( int i=0; i< static_cast<int>( index_to_orb_off_diag_flip.size()); i++)
        {
            for( list<Ldag_L>::iterator it = off_diag_J_bath.Ldag_L_config[i].begin(); it != off_diag_J_bath.Ldag_L_config[i].end(); it++)
            {
                t = ( it->te-it->ts );
                if( t<0 ){ t += beta;}
                G_two_off_diag_flip[i][floor( t/dt )] += sign/(beta*dt*it->jt);
                
            }
            
        }
    }
    else if(if_meas_diag_flip_w)
    {// measure w
        
        for( int i=0; i< static_cast<int>( diag_flip_index.size()); i++)
        {
            for( int n=0; n<N_w; n++)
            {
                for( list<Ldag_L>::iterator it = off_diag_J_bath.Ldag_L_config[i].begin(); it != off_diag_J_bath.Ldag_L_config[i].end(); it++)
                {
                    t = ( it->te-it->ts );
                    
                    G_two_off_diag_flip[i][n] += sign/(beta*it->jt)*cos(bose_nu_n[n]*t);
                    
                }
            }
        }
        
    }
    
    
}


void accumulator::print_G_one()
{
    if( !if_meas_one_tau ) return;
    
    ofstream file;
    file.open("G_one_tau.txt",ios::trunc);
    if(!file.is_open())
    {
        cerr<<" can't open G_one_tau.txt"<<endl;
        return;
    }
    
    for( map< pair<int, int>, int >::iterator it=G_one_index.begin(); it!=G_one_index.end(); it++)
    {
        
        int i= (it->first).first;
        int j= (it->first).second;
        int k= it->second;
        
       
        
        file<<endl;
        file<<endl;
        file<<"d_"<<i<<"^dag(tau) "<<"d_"<<j<<endl;
        file<<"tau    "<<"G(tau)"<<endl;
    
        
        if( i== j)
            file<<0<<"    "<<n_i_n_j[i][i]/sign_sum_static<<endl;
        for(int i_t = 0; i_t<N_t; i_t++)
        {
            file<<(i_t+0.5)*dt<<"    "<<G_one[k][i_t]/sign_sum_G_one<<endl;
        }
        if( i== j)
            file<<beta<<"    "<<1.0-(n_i_n_j[i][i])/sign_sum_static<<endl;
    
    
        
        
        
    }
    file.close();
    
    
}

void accumulator::print_G_two_diag_flip()
{
    if( !if_meas_diag_flip_w && !if_meas_diag_flip_tau ) return;
    
    ofstream file;
    file.open("G_two_diag_flip.txt",ios::trunc);
    if(!file.is_open())
    {
        cerr<<" can't open G_two_diag_flip.txt"<<endl;
        return;
    }
    int n_ = sign_sum_diag_flip;
    int k;
    
    if( if_meas_diag_flip_tau)
    {
        file<<" L(tau) L^dag(0) "<<endl;
        file<<endl;
        
        for( map< pair<int,int>,int >::iterator it = diag_flip_index.begin(); it != diag_flip_index.end(); it++)
        {
            k=it->second;
            file<<"L = "<<"d_"<<(it->first).first<<"^dag d_"<<(it->first).second<<endl;
            for(int i_t = 0; i_t<N_t; i_t++)
            {
                file<<(i_t+0.5)*dt<<"    "<<G_two_diag_flip[k][i_t]/n_<<endl;
            }
            file<<endl;
            file<<endl;
        }
    }
    else if(if_meas_diag_flip_w)
    {
        file<<" L(tau) L^dag(0) "<<endl;
        file<<endl;

        
        for( map< pair<int,int>,int >::iterator it = diag_flip_index.begin(); it != diag_flip_index.end(); it++)
        {
            k=it->second;
            file<<"L = "<<"d_"<<(it->first).first<<"^dag d_"<<(it->first).second<<endl;
            file<<"omega_n"<<"   "<<"L Ldag(omega)"<<endl;
            for(int i_n = 0; i_n<N_w; i_n++)
            {
                file<<bose_nu_n[i_n]<<"    "<<G_two_diag_flip[k][i_n]/n_<<endl;
            }
            file<<endl;
            vector< double > ft = boson_inv_fourier(N_t, beta, G_two_diag_flip[k] );
            
            
            file<<"tau"<<"   "<<"L Ldag(tau)"<<endl;
            for( int i_t =0; i_t<N_t; i_t++)
            {
                file<<i_t*dt<<"    "<<ft[i_t]/n_<<endl;
            }
            
            file<<endl;
        }
        
        
    }
    
    
    file.close();
    
    
}


void accumulator::print_G_two_off_diag_flip()
{
    if( !if_meas_off_diag_flip_w && !if_meas_off_diag_flip_tau ) return;
    
    ofstream file;
    file.open("G_two_off_diag_flip.txt",ios::trunc);
    if(!file.is_open())
    {
        cerr<<" can't open G_two_off_diag_flip.txt"<<endl;
        return;
    }
    int n_ = sign_sum_off_diag_flip;

    
    if( if_meas_off_diag_flip_tau)
    {
        file<<" L_1(tau) L_2(0) "<<endl;
        file<<endl;
        
        for(  int i=0 ;i< static_cast<int>( index_to_orb_off_diag_flip.size());i++)
        {
            file<<"d_"<<index_to_orb_off_diag_flip[i][0]<<"^dag(tau) "<<"d_"<<index_to_orb_off_diag_flip[i][1]<<"(tau)  ";
            file<<"d_"<<index_to_orb_off_diag_flip[i][2]<<"^dag(0) "<<"d_"<<index_to_orb_off_diag_flip[i][3]<<"(0)  "<<endl;
            for(int i_t = 0; i_t<N_t; i_t++)
            {
                file<<(i_t+0.5)*dt<<"    "<<G_two_off_diag_flip[i][i_t]/n_<<endl;
            }
            file<<endl;
            file<<endl;
        }
    }
    else if(if_meas_diag_flip_w)
    {
        file<<" L_1(tau) L_2(0) "<<endl;
        file<<endl;
        
        
        for(  int i=0 ;i< static_cast<int>( index_to_orb_off_diag_flip.size());i++)
        {
           
            file<<"L_1 = "<<"d_"<<index_to_orb_off_diag_flip[i][0]<<"^dag"<<" "<<"d_"<<index_to_orb_off_diag_flip[i][1]<<endl;
            file<<"L_2 = "<<"d_"<<index_to_orb_off_diag_flip[i][2]<<"^dag"<<" "<<"d_"<<index_to_orb_off_diag_flip[i][3]<<endl;
            
            file<<"omega_n"<<"   "<<"L L (omega)"<<endl;
            for(int i_n = 0; i_n<N_w; i_n++)
            {
                file<<bose_nu_n[i_n]<<"    "<<G_two_off_diag_flip[i][i_n]/n_<<endl;
            }
            file<<endl;
            vector< double > ft = boson_inv_fourier(N_t, beta, G_two_off_diag_flip[i] );
            
            
            file<<"tau"<<"   "<<"L L (tau)"<<endl;
            for( int i_t =0; i_t<N_t; i_t++)
            {
                file<<i_t*dt<<"    "<<ft[i_t]/n_<<endl;
            }
            
            file<<endl;
        }
        
        
    }
    
    
    file.close();
    
    
}
void accumulator::print_static()
{
    ofstream file;
    file.open("static.txt",ios::trunc);
    
    double n_total =sign_sum_static;
   
    
    file<<"<d^dag d>"<<endl;
    
    for( int k=0; k< norb;k++)
    {
        file<<k<<"    "<<n_i_n_j[k][k]/n_total<<endl;
    }
    file<<endl;
    
    file<<"<n_i n_j>"<<endl;
    for( int i=0;i<norb; i++)
    {
        for( int j=0;j<norb; j++)
        {
             file<<i<<"   "<<j<<"    "<<n_i_n_j[i][j]/n_total<<endl;
        }
    }
    
    
    file<<endl;
    file<<endl;
    
    file<<"n_meas_static    average_sign"<<endl;
    file<<n_meas_static<<"    "<<sign_sum_static/n_meas_static<<endl;
    
    if( if_meas_one_tau)
    {
        file<<"n_meas_g_one_tau    average_sign"<<endl;
        file<<n_meas_G_one<<"    "<<sign_sum_G_one/n_meas_G_one<<endl;;
    }
    
    file<<endl;
    file<<endl;
    file<<"histogram for fermion bath"<<endl<<endl;
    for( int i=0; i< static_cast<int>(histogram_fermion_bath.size()); i++)
    {
        file<<i<<"-th fermion bath"<<endl;
        file<<"order     number"<<endl;
        for(int j=0; j< static_cast<int>(histogram_fermion_bath[i].size()); j++)
        {
            if( histogram_fermion_bath[i][j]!=0)
            {
                file<<j<<"    "<<histogram_fermion_bath[i][j]<<endl;
            }
        }
        file<<endl;
        
        
    }
    
    file<<"histogram for diag flip bath"<<endl<<endl;
    for( int k=0; k<static_cast<int>( diag_flip_index.size() ); k++)
    {
        file<<k<<"-th diag flip bath "<<endl;
        file<<"order     number"<<endl;
        for(int j=0; j< static_cast<int>(histogram_diag_flip_bath[k].size()); j++)
        {
            if(histogram_diag_flip_bath[k][j]!=0)
            {
                file<<j<<"    "<<histogram_diag_flip_bath[k][j]<<endl;
            }
        }
        file<<endl;
    }
    
    file<<"histogram for off diag flip bath"<<endl<<endl;
    for( int k=0; k<static_cast<int>( index_to_orb_off_diag_flip.size() ); k++)
    {
        file<<k<<"-th off diag flip bath "<<endl;
        file<<"order     number"<<endl;
        for(int j=0; j< static_cast<int>(histogram_off_diag_flip_bath[k].size()); j++)
        {
            if(histogram_off_diag_flip_bath[k][j]!=0)
            {
                file<<j<<"    "<<histogram_off_diag_flip_bath[k][j]<<endl;
            }
        }
        file<<endl;
    }
    

    file.close();
    
    
}


void accumulator::print_n_n()
{
    if( !if_meas_nn_tau ) return;
    
    ofstream file;
    file.open("n_n.txt",ios::trunc);
    if(!file.is_open())
    {
        cerr<<" can't open n_n.txt"<<endl;
        return;
    }
    int n_ = sign_sum_n_n_;
    int k;
    
    if( if_meas_nn_tau)
    {
        file<<" n_i(tau) n_j(0) "<<endl;
        file<<endl;
        
        for( map< pair<int,int>,int >::iterator it = orb_to_index_nn.begin(); it != orb_to_index_nn.end(); it++)
        {
            k=it->second;
            n_n_[k][N_t] = n_n_[k][0];
            file<<" n_"<<(it->first).first<<"(tau) n_"<<(it->first).second<<"(0)"<<endl;
            for(int i_t = 0; i_t<N_t; i_t++)
            {
                file<<(i_t)*dt<<"    "<<n_n_[k][i_t]/n_<<endl;
            }
            file<<endl;
            file<<endl;
        }
        
        file<<" S_Z(tau) S_Z(0) "<<endl;
        file<<endl;
        for( int i_t = 0; i_t<N_t;i_t++)
        {
            file<<(i_t)*dt<<"   "<<0.25*(n_n_[0][i_t]+n_n_[2][i_t]-n_n_[1][i_t]-n_n_[1][(N_t-i_t)%N_t])/n_<<endl;
            
        }
       
    }
    
    
    
    file.close();
    
    
}




void accumulator::measure_g_two( vector<fermi_bath_block> &fermi_bath, int sign)
{
    n_meas_G_two ++;
    sign_sum_G_two += sign;
    
    if( if_meas_two_tau )
    {
        for( int i=0; i< static_cast<int>(G_two_orb_index.size()); i++)
        {
            measure_two_aux_tau(fermi_bath,sign,i);
        }
    }
    return ;
}

void accumulator::measure_two_aux_tau( vector<fermi_bath_block> &fermi_bath, int sign, int type)
{
    int bk_1 = G_two_block_index[type][0];
    int bk_2 = G_two_block_index[type][2];
    double t1,t2;
    double tau;
    int k;
    int sign_order;
    // assume bk1 != bk2
    
    vector< vector< double > > t_st_vec( 2, vector<double>(0));
    vector< vector< double > > t_en_vec( 2, vector<double>(0));
    vector< vector< double > > index_st( 2, vector<double>(0));
    vector< vector< double > > index_en( 2, vector<double>(0)); // index of the N_ matrix
//    vector< vector< int > > sign_vec( 2, vector<int>(0));
    
    // find the c^dag c at the same time point

    for(  map< pair<int,double> ,int>::iterator it_st = fermi_bath[bk_1].time_start.begin() ;
        it_st != fermi_bath[bk_1].time_start.end(); it_st ++)
    {
        if( it_st->first.first != G_two_orb_index[type][0]) continue;
        t1 = it_st->first.second;
        for(  map< pair<int,double> ,int>::iterator it_en = fermi_bath[bk_2].time_end.begin() ;
            it_en != fermi_bath[bk_2].time_end.end(); it_en++)
        {
            if( it_en->first.first != G_two_orb_index[type][1]) continue;
            
            t2 = it_en->first.second;
            if( abs( t1-t2) > 0.5*dt && abs(t1-t2)<beta-0.5*dt) continue;
            
            t_st_vec[0].push_back(t1);
            t_en_vec[0].push_back(t2);
//            if( t1-t2<0.0){ sign_vec[0].push_back(-1);}
//            else{  sign_vec[0].push_back(1);}
            index_st[0].push_back(it_st->second);
            index_en[0].push_back(it_en->second);
            
            
            
        }
        
        
    }
    
    for(  map< pair<int,double> ,int>::iterator it_st = fermi_bath[bk_2].time_start.begin() ;
        it_st != fermi_bath[bk_2].time_start.end(); it_st ++)
    {
        if( it_st->first.first != G_two_orb_index[type][2]) continue;
        t1 = it_st->first.second;
        for(  map< pair<int,double> ,int>::iterator it_en = fermi_bath[bk_1].time_end.begin() ;
            it_en != fermi_bath[bk_1].time_end.end(); it_en++)
        {
            if( it_en->first.first != G_two_orb_index[type][3]) continue;
            
            t2 = it_en->first.second;
            if( abs( t1-t2) > 0.5*dt && abs(t1-t2)<beta-0.5*dt) continue;
            
            t_st_vec[1].push_back(t1);
            t_en_vec[1].push_back(t2);
//            if( t1-t2<0.0){ sign_vec[1].push_back(-1);}
//            else{  sign_vec[1].push_back(1);}
            index_st[1].push_back(it_st->second);
            index_en[1].push_back(it_en->second);
            
            
        }
        
        
    }
    for( int it_0 = 0; it_0 < static_cast<int>(t_st_vec[0].size()); it_0 ++ )
    {
        
        for( int it_1 = 0; it_1 < static_cast<int>(t_st_vec[1].size()); it_1 ++ )
        {
            tau = t_st_vec[0][it_0] - t_st_vec[1][it_1];
            if(tau<0) tau+= beta;
            k = floor( tau/dt );
            sign_order = 1;
            if( t_st_vec[0][it_0] < t_en_vec[1][it_1]){ sign_order = -sign_order; }
            if( t_en_vec[0][it_0] < t_st_vec[1][it_1]){sign_order = - sign_order;}
           
            G_two[type][k] += sign* sign_order * (-1.0) * fermi_bath[bk_1].N_( index_en[1][it_1] , index_st[0][it_0]   ) * fermi_bath[bk_2].N_( index_en[0][it_0] , index_st[1][it_1]   ) / (dt*dt*dt*beta);
            
        }
        
    }
    
    
    
    
}






void accumulator::print_G_two()
{
    if( !if_meas_one_tau ) return;
    
    ofstream file;
    file.open("G_two_tau.txt",ios::trunc);
    if(!file.is_open())
    {
        cerr<<" can't open G_two_tau.txt"<<endl;
        return;
    }
    
    int N_ = sign_sum_G_two;
    
    for( int i=0; i< static_cast<int>(G_two_orb_index.size()); i++ )
    {
        file<<"d^dag_"<<G_two_orb_index[i][0]<<"(tau) ";
        file<<"d_"<<G_two_orb_index[i][1]<<"(tau) ";
        file<<"d^dag_"<<G_two_orb_index[i][2]<<"(0) ";
        file<<"d_"<<G_two_orb_index[i][3]<<"(0) ";
        file<<endl;
        file<<"tau         G_two(tau) "<<endl;
        
        for( int t=0; t<N_t; t++)
        {
            file<<(t+0.5)*dt<<"   "<<G_two[i][t]/N_<<endl;
        }
        
        file<<endl;
        file<<endl;
        
        
        
    }
    file.close();
    
    
}


double accumulator::avg_fermion_bath_num()
{
    double seg_num=0.0;
    int N_=0;
    double avg_fm = 0.0;
    for( int i=0; i<norb ;i++)
    {
        seg_num = 0.0;
        N_= 0;
        for( int n =0 ; n<max_fermi; n++ )
        {
            if(histogram_fermion_bath[i][n]!=0)
            {
                N_ += histogram_fermion_bath[i][n];
                seg_num += histogram_fermion_bath[i][n]*n;
            }
        }
        avg_fm += seg_num/N_;
        
    }
    

    return avg_fm;
    
    
}
double accumulator::avg_diag_flip_bath_num()
{
    double seg_num=0.0;
    int N_=0;
    double avg_diag_flip = 0.0;
    
    for( int i=0; i< static_cast<int>(diag_flip_index.size() ); i++)
    {
        seg_num = 0.0;
        N_ = 0;
        for( int n =0 ; n<max_diag_flip; n++ )
        {
            if(histogram_diag_flip_bath[i][n]!=0)
            {
                N_ += histogram_diag_flip_bath[i][n];
                seg_num += histogram_diag_flip_bath[i][n]*n;
            }
        }
        avg_diag_flip += seg_num/N_;
        
    }
    return avg_diag_flip;
}
