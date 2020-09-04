////
////  accumulator_mpi.cpp
////  cthyb
////
////  Created by 胡昊昱 on 8/7/19.
////  Copyright © 2019 胡昊昱. All rights reserved.
////
//
//#include "accumulator.hpp"
//
//#include <mpi.h>
//
//bool mpi_reduce_double_vec(vector<double> &vec_loc_sum, int sum_sign, vector<double> &vec_global, vector<double> &vec_error)
//{// average each core and then mpi average;
//    int size;
//    int rank;
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//
//    vector<double> vec_loc(vec_loc_sum.size(),0.0);
//    for( int i=0; i<static_cast<int>(vec_loc.size()); i++) {  vec_loc[i] = vec_loc_sum[i]/sum_sign; }
//
//    vector<double> vec_square(vec_loc.size(),0.0);
//    for( int i=0; i<static_cast<int>(vec_loc.size()) ; i++)
//    {
//        vec_square[i] = vec_loc[i]*vec_loc[i];
//        MPI_Reduce(&vec_loc[i],&vec_global[i],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//        MPI_Reduce(&vec_square[i],&vec_error[i],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    }
//    if( rank==0 )
//    {
//        for( int i=0; i<static_cast<int>(vec_loc.size()) ; i++)
//        {
//            vec_global[i] = vec_global[i]/size;
//            vec_error[i] = vec_error[i] - size*vec_global[i]*vec_error[i];
//            vec_error[i] = sqrt( vec_error[i]/(size) * (size-1) );
//        }
//    }
//    return true;
//}
//
//bool mpi_reduce_double_vec_vec(vector< vector<double> > &vec_loc_sum,int sum_sign ,vector< vector<double> >&vec_global, vector< vector<double> > &vec_error)
//{
//    for( int i=0; i< static_cast<int>(vec_loc_sum.size());i++){   mpi_reduce_double_vec(vec_loc_sum[i], sum_sign, vec_global[i],vec_error[i]);          }
//    return true;
//}
//
//bool accumulator::average_mpi()
//{
//    if( if_average ){ return false;}
//    if_avg_mpi = true;
//
//    MPI_Reduce(&sign_sum_static,&sign_sum_static_mpi,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    MPI_Reduce(&n_meas_static,&n_meas_static_mpi,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
//
//
//
//    mpi_reduce_double_vec_vec(n_i_n_j,sign_sum_static, n_i_n_j_mpi, n_i_n_j_mpi_error);
//
//    if( if_meas_one_tau )
//    {
//        mpi_reduce_double_vec_vec(G_one, sign_sum_G_one, G_one_mpi, G_one_mpi_error);
//    }
//
//    if( if_meas_two_tau )
//    {
//        mpi_reduce_double_vec_vec(G_two,sign_sum_G_two, G_two_mpi, G_two_mpi_error);
//    }
//
//    if( if_meas_diag_flip_tau || if_meas_diag_flip_w )
//    {
//        mpi_reduce_double_vec_vec(G_two_diag_flip,sign_sum_diag_flip, G_two_diag_flip_mpi, G_two_diag_flip_mpi_error);
//
//    }
//
//    if( if_meas_nn_tau )
//    {
//        mpi_reduce_double_vec_vec(n_n_,sign_sum_n_n_, n_n_mpi, n_n_mpi_error);
//    }
//
//    return true;
//}
//
//
//
//void accumulator::print_G_one_mpi()
//{
//    if( !if_meas_one_tau ) return;
//
//    ofstream file;
//    file.open("G_one_tau_mpi.txt",ios::trunc);
//    if(!file.is_open())
//    {
//        cerr<<" can't open G_one_tau_mpi.txt"<<endl;
//        return;
//    }
//
//    for( map< pair<int, int>, int >::iterator it=G_one_index.begin(); it!=G_one_index.end(); it++)
//    {
//
//        int i= (it->first).first;
//        int j= (it->first).second;
//        int k= it->second;
//
//
//
//        file<<endl;
//        file<<endl;
//        file<<"d_"<<i<<"^dag(tau) "<<"d_"<<j<<endl;
//        file<<"tau    "<<"G(tau)    "<<"error "<<endl;
//
//
//        if( i== j)
//            file<<0<<"    "<<n_i_n_j[i][i]<<endl;
//        for(int i_t = 0; i_t<N_t; i_t++)
//        {
//            file<<(i_t+0.5)*dt<<"    "<<G_one_mpi[k][i_t]<<" "<<G_one_mpi_error[k][i_t]<<endl;
//        }
//        if( i== j)
//            file<<beta<<"    "<<1.0-n_i_n_j_mpi[i][i]<<"  "<<n_i_n_j_mpi_error[i][j]<<endl;
//
//
//
//
//    }
//    file.close();
//
//
//}
//
//void accumulator::print_G_two_diag_flip_mpi()
//{
//    if( !if_meas_diag_flip_w && !if_meas_diag_flip_tau ) return;
//
//    ofstream file;
//    file.open("G_two_diag_flip_mpi.txt",ios::trunc);
//    if(!file.is_open())
//    {
//        cerr<<" can't open G_two_diag_flip_mpi.txt"<<endl;
//        return;
//    }
//    int k;
//
//    if( if_meas_diag_flip_tau)
//    {
//        file<<" L(tau) L^dag(0) "<<endl;
//        file<<endl;
//
//        for( map< pair<int,int>,int >::iterator it = diag_flip_index.begin(); it != diag_flip_index.end(); it++)
//        {
//            k=it->second;
//            file<<"L = "<<"d_"<<(it->first).first<<"^dag d_"<<(it->first).second<<endl;
//            for(int i_t = 0; i_t<N_t; i_t++)
//            {
//                file<<(i_t+0.5)*dt<<"    "<<G_two_diag_flip_mpi[k][i_t]<<"   "<<G_two_diag_flip_mpi_error[k][i_t]<<endl;
//            }
//            file<<endl;
//            file<<endl;
//        }
//    }
//    else if(if_meas_diag_flip_w)
//    {
//        file<<" L(tau) L^dag(0) "<<endl;
//        file<<endl;
//
//
//        for( map< pair<int,int>,int >::iterator it = diag_flip_index.begin(); it != diag_flip_index.end(); it++)
//        {
//            k=it->second;
//            file<<"L = "<<"d_"<<(it->first).first<<"^dag d_"<<(it->first).second<<endl;
//            file<<"omega_n"<<"   "<<"L Ldag(omega)"<<endl;
//            for(int i_n = 0; i_n<N_w; i_n++)
//            {
//                file<<bose_nu_n[i_n]<<"    "<<G_two_diag_flip_mpi[k][i_n]<<"    "<<G_two_diag_flip_mpi_error[k][i_n]<<endl;
//            }
//            file<<endl;
//            vector< double > ft = boson_inv_fourier(N_t, beta, G_two_diag_flip_mpi[k] );
//
//
//            file<<"tau"<<"   "<<"L Ldag(tau)"<<endl;
//            for( int i_t =0; i_t<N_t; i_t++)
//            {
//                file<<i_t*dt<<"    "<<ft[i_t]<<endl;
//            }
//
//            file<<endl;
//        }
//
//
//    }
//
//
//    file.close();
//
//
//}
//void accumulator::print_static_mpi()
//{
//    ofstream file;
//    file.open("static_mpi.txt",ios::trunc);
//
//
//    file<<"<d^dag d>"<<endl;
//
//    for( int k=0; k< norb;k++)
//    {
//        file<<k<<"    "<<n_i_n_j_mpi[k][k]<<"   "<<n_i_n_j_mpi_error[k][k]<<endl;
//    }
//    file<<endl;
//
//    file<<"<n_i n_j>"<<endl;
//    for( int i=0;i<norb; i++)
//    {
//        for( int j=0;j<norb; j++)
//        {
//            file<<i<<"   "<<j<<"    "<<n_i_n_j_mpi[i][j]<<"   "<<n_i_n_j_mpi_error[i][j]<<endl;
//        }
//    }
//
//
//    file<<endl;
//    file<<endl;
//
//    file<<"n_meas_static    average_sign"<<endl;
//    file<<n_meas_static_mpi<<"    "<<sign_sum_static_mpi/n_meas_static_mpi<<endl;
//
//
//    file<<endl;
//    file<<endl;
//    file<<"histogram for fermion bath"<<endl<<endl;
//    for( int i=0; i< static_cast<int>(histogram_fermion_bath.size()); i++)
//    {
//        file<<i<<"-th fermion bath"<<endl;
//        file<<"order     number"<<endl;
//        for(int j=0; j< static_cast<int>(histogram_fermion_bath[i].size()); j++)
//        {
//            if( histogram_fermion_bath[i][j]!=0)
//            {
//                file<<j<<"    "<<histogram_fermion_bath[i][j]<<endl;
//            }
//        }
//        file<<endl;
//
//
//    }
//
//    file<<"histogram for diag flip bath"<<endl<<endl;
//    for( int k=0; k<static_cast<int>( diag_flip_index.size() ); k++)
//    {
//        file<<k<<"-th diag flip bath "<<endl;
//        file<<"order     number"<<endl;
//        for(int j=0; j< static_cast<int>(histogram_diag_flip_bath[k].size()); j++)
//        {
//            if(histogram_diag_flip_bath[k][j]!=0)
//            {
//                file<<j<<"    "<<histogram_diag_flip_bath[k][j]<<endl;
//            }
//        }
//        file<<endl;
//    }
//
//
//    file.close();
//
//
//}
//
//
//void accumulator::print_n_n_mpi()
//{
//    if( !if_meas_nn_tau ) return;
//
//    ofstream file;
//    file.open("n_n_mpi.txt",ios::trunc);
//    if(!file.is_open())
//    {
//        cerr<<" can't open n_n_mpi.txt"<<endl;
//        return;
//    }
//    int k;
//
//    if( if_meas_nn_tau)
//    {
//        file<<" n_i(tau) n_j(0) "<<endl;
//        file<<endl;
//
//        for( map< pair<int,int>,int >::iterator it = orb_to_index_nn.begin(); it != orb_to_index_nn.end(); it++)
//        {
//            k=it->second;
//            file<<" n_"<<(it->first).first<<"(tau) n_"<<(it->first).second<<"(0)"<<endl;
//            for(int i_t = 0; i_t<N_t; i_t++)
//            {
//                file<<(i_t)*dt<<"    "<<n_n_mpi[k][i_t]<<"    "<<n_n_mpi_error[k][i_t]<<endl;
//            }
//            file<<endl;
//            file<<endl;
//        }
//
//        file<<" S_Z(tau) S_Z(0) "<<endl;
//        file<<endl;
//        for( int i_t = 0; i_t<N_t;i_t++)
//        {
//            double data=0.0,  data_error=0.0;
//            int orb1_up,orb2_up, orb1_down,orb2_down;
//            for( int i=0 ;i<norb/2; i++)
//            {
//                orb1_up = 2*i;
//                orb1_down = 2*i+1;
//                for( int j=0; j<norb/2;j++)
//                {
//                    orb2_up = 2*j;
//                    orb2_down = 2*j+1;
//
//                    int index_upup = orb1_up<orb2_up ? ( orb_to_index_nn[ pair<int,int>(orb1_up,orb2_up)] ): ( orb_to_index_nn[ pair<int,int>(orb2_up,orb1_up)] );
//                    int index_dndn = orb1_down<orb2_down ? ( orb_to_index_nn[ pair<int,int>(orb1_down,orb2_down)] ): ( orb_to_index_nn[ pair<int,int>(orb2_down,orb1_down)] );
//                    int index_updown =  orb1_up<orb2_down ? ( orb_to_index_nn[ pair<int,int>(orb1_up,orb2_down)] ): ( orb_to_index_nn[ pair<int,int>(orb2_down,orb1_up)] );
//                    int index_downup =  orb1_down<orb2_up ? ( orb_to_index_nn[ pair<int,int>(orb1_down,orb2_up)] ): ( orb_to_index_nn[ pair<int,int>(orb2_up,orb1_down)] );
//
//
//
//                    data += 0.25*( n_n_mpi[ index_upup ][i_t] + n_n_mpi[index_dndn][i_t] - n_n_mpi[index_updown][i_t] - n_n_mpi[index_downup][i_t]   );
//                    data_error += 0.25*( n_n_mpi_error[ index_upup ][i_t] + n_n_mpi_error[index_dndn][i_t] - n_n_mpi_error[index_updown][i_t] - n_n_mpi_error[index_downup][i_t]   );
//                }
//
//            }
//            file<<(i_t)*dt<<"   "<< data<<"   "<<data_error<<endl;
//
//        }
//
//    }
//
//
//
//    file.close();
//
//
//}
//
//
//
//
//
//void accumulator::print_G_two_mpi()
//{
//    if( !if_meas_one_tau ) return;
//
//    ofstream file;
//    file.open("G_two_tau_mpi.txt",ios::trunc);
//    if(!file.is_open())
//    {
//        cerr<<" can't open G_two_tau_mpi.txt"<<endl;
//        return;
//    }
//
//
//    for( int i=0; i< static_cast<int>(G_two_orb_index.size()); i++ )
//    {
//        file<<"d^dag_"<<G_two_orb_index[i][0]<<"(tau) ";
//        file<<"d_"<<G_two_orb_index[i][1]<<"(tau) ";
//        file<<"d^dag_"<<G_two_orb_index[i][2]<<"(0) ";
//        file<<"d_"<<G_two_orb_index[i][3]<<"(0) ";
//        file<<endl;
//        file<<"tau         G_two(tau) "<<endl;
//
//        for( int t=0; t<N_t; t++)
//        {
//            file<<(t+0.5)*dt<<"   "<<G_two_mpi[i][t]<<"    "<<G_two_mpi_error[i][t]<<endl;
//        }
//
//        file<<endl;
////        file<<endl;
//
//
//
//    }
//    file.close();
//
//
//}
//
//
//void accumulator::print_all_mpi()
//{
//    if( !if_avg_mpi) return;
//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    if(rank !=0 )return;
//    print_static_mpi();
//    if( if_meas_one_tau){ print_G_one_mpi();}
//    if( if_meas_diag_flip_tau || if_meas_diag_flip_w ){  print_G_two_diag_flip_mpi();}
//    if( if_meas_two_tau){ print_G_two_mpi();}
//    if( if_meas_nn_tau){ print_n_n_mpi(); }
//}
