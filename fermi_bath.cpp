//
//  fermi_bath.cpp
//  cthyb
//
//  Created by 胡昊昱 on 2019/7/4.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#include "fermi_bath.hpp"

std::vector< std::vector< fermi_type > > fermi_bath_block::delta_tau={};

fermi_type fermi_bath_block::try_insert(pair< int,double> &node_start, pair<int,double> &node_end)
{

    
    insert_lock = false;
    
    
    node_start_insert = node_start;
    node_end_insert = node_end;
    S_insert= interpolate_delta(node_start,node_end);
    S_tilde_insert = 1.0/S_insert;
    if( nsize == 0)
    {// empty
        return S_insert;
    }
    
    R_insert.resize(nsize);
    Q_insert.resize(nsize);
    
    for( map< pair<int,double> ,int>::iterator it=time_start.begin(); it!=time_start.end(); it++)
    {
        Q_insert[it->second] = interpolate_delta( it->first, node_end );
    }
    for( map< pair<int,double> ,int>::iterator it=time_end.begin(); it!=time_end.end(); it++)
    {
        R_insert[it->second] = interpolate_delta( node_start, it->first);
    }
    aux_vec = prod(N_,Q_insert);// temp vecs
    S_tilde_insert = inner_prod(R_insert,aux_vec);
    S_tilde_insert = 1.0/( S_insert-S_tilde_insert );
    
    return 1.0/S_tilde_insert;
    
}
 

bool fermi_bath_block::do_insert(  )
{
    if( insert_lock ) return false;
    
    insert_lock = true;
    
    using namespace boost::numeric::ublas;
    
    int last_ind = nsize;
    nsize++;
    
    
    if(nsize == 1)
    {// insert from empty
        time_start.insert( make_pair(node_start_insert,last_ind) );
        time_end.insert( make_pair(node_end_insert,last_ind) );
        det = S_insert;
       
        N_inv.resize(nsize,nsize);
        N_.resize(nsize, nsize);
        N_inv(last_ind,last_ind) = S_insert;
        N_(last_ind,last_ind) = 1.0/S_insert;
        return true;
    }
    
    
    // insert time
    time_start.insert( make_pair(node_start_insert,last_ind) );
    time_end.insert( make_pair(node_end_insert,last_ind) );
    
    det = det/S_tilde_insert;
    
//    Q_tilde_insert = aux_vec*(-S_tilde_insert);
//    R_tilde_insert = prod(R_insert,N_)*(-S_tilde_insert);
//    aux_mat = outer_prod(Q_tilde_insert,R_tilde_insert)/S_tilde_insert;
    // s_tilde may be very large so this part may give large error bar.
    
    Q_tilde_insert = aux_vec;
    R_tilde_insert = prod(R_insert,N_);
    aux_mat = outer_prod(Q_tilde_insert,R_tilde_insert)*S_tilde_insert;
    
    Q_tilde_insert = Q_tilde_insert * (-S_tilde_insert);
    R_tilde_insert = R_tilde_insert * (-S_tilde_insert);;
    
    
    N_ = N_ + aux_mat;
    
    N_.resize(nsize, nsize,true);// resize while preserving ele;
    N_inv.resize(nsize, nsize,true);

    
    R_insert.resize( nsize, true );
    Q_insert.resize( nsize, true );
    R_tilde_insert.resize( nsize, true );
    Q_tilde_insert.resize( nsize, true );
    
    row(N_inv,last_ind) = R_insert;
    column(N_inv,last_ind)=  Q_insert;
    row(N_,last_ind) = R_tilde_insert;
    column(N_,last_ind) = Q_tilde_insert;
    N_(last_ind,last_ind) = S_tilde_insert;
    N_inv(last_ind,last_ind) = S_insert;
    
    
    
    return true;
    
    
}


fermi_type fermi_bath_block::try_remove(int k_st, int k_en)
{
 
    remove_lock = false;
    
    if( nsize ==0 )
    {
        return 0.0;
    }
    
    map< pair<int,double> ,int>::iterator it_st = time_start.begin();
    for( int i=0; i<k_st;i++){ it_st++;}
    it_st ++;
    map< pair<int,double> ,int>::iterator it_en = time_end.begin();
    for( int i=0; i<k_en;i++){ it_en++;}
    
    row_remove = it_st->second;
    col_remove = it_en->second;
    
    node_start_remove= it_st->first;
    node_end_remove = it_en->first;
    return try_remove_aux();
}

fermi_type fermi_bath_block::try_remove(pair< int,double> &node_start, pair<int,double> &node_end)
{
    remove_lock = false;
    
    if( nsize ==0 )
    {
        return 0.0;
    }
    
    row_remove = time_start[node_start];
    col_remove = time_end[node_end];
    node_start_remove = node_start;
    node_end_remove = node_end;
    return try_remove_aux();
    
    
}

fermi_type fermi_bath_block::try_remove_aux()
{
    
    fermi_type ratio = N_(col_remove, row_remove);
    if( col_remove != nsize -1 ){ ratio = -ratio;}
    if( row_remove != nsize -1 ){ ratio = -ratio;}
    
    return ratio;
    
}


bool fermi_bath_block::do_remove(int &swap_sign )
{
    if( remove_lock ) return false;
    remove_lock = true;
    
    using namespace boost::numeric::ublas;
    swap_sign = 1;
    
    S_tilde_remove =N_(col_remove, row_remove);
    
    if(nsize == 1)
    {
        time_start.erase(node_start_remove);
        time_end.erase(node_end_remove);// erase time
        nsize = 0;
        N_inv.resize(0,0);
        N_.resize(0, 0);
        
        det = 1.0;
        return true;
    }
    nsize--;
    det = det * S_tilde_remove;
    //det and size
    
   
    
    time_start.erase(node_start_remove);
    time_end.erase(node_end_remove);// erase time
    
    
    
    if( row_remove != nsize )
    {
        row(N_inv,row_remove).swap( row(N_inv, nsize) );
        column(N_,row_remove).swap( column(N_, nsize) ); // swap
        
        det = -det;
        swap_sign = -swap_sign;
        for(map< pair<int,double> ,int>::iterator it = time_start.begin(); it!=time_start.end(); it++)
        {
            if( it->second == nsize ){ it->second = row_remove;  break;}
        }
    }
    if(col_remove != nsize)
    {
        column(N_inv,col_remove).swap( column(N_inv,nsize) );
        row(N_,col_remove).swap( row(N_,nsize) ); // swap
        
        det = -det;
        swap_sign = -swap_sign;
        for( map< pair<int,double> ,int>::iterator it = time_end.begin(); it!=time_end.end(); it++)
        {
            if( it->second == nsize ){ it->second = col_remove; break;}
        }
    }// swap time array

    
    
    
    N_inv.resize(nsize,nsize,true);
    aux_mat = project(N_,range(0,nsize),range(0,nsize) );
    R_tilde_remove = row(N_,nsize);
    R_tilde_remove.resize(nsize);
    Q_tilde_remove = column(N_,nsize);
    Q_tilde_remove.resize(nsize);
    N_ = aux_mat - outer_prod(Q_tilde_remove, R_tilde_remove)/S_tilde_remove;
    // update N,N_inv;
    
    
    return true;
    
}




fermi_type fermi_bath_block::try_replace( pair<int,double> &node_st_old, pair<int,double> &node_st_new,pair<int,double> &node_en_old, pair<int,double> &node_en_new)
{
    replace_lock = false;

    using namespace boost::numeric::ublas;
    
    if(nsize == 0) return 0.0;
    
    node_start_insert = node_st_new;
    node_end_insert = node_en_new;
    node_start_remove = node_st_old;
    node_end_remove= node_en_old;
    
//    if(  time_start.find( node_st_old ) == time_start.end() )
//    {
//        cout<<node_st_old.second<<endl;
//        cout<<"can't find one"<<endl;
//    }
//    
//    if(  time_end.find( node_en_old ) == time_end.end() )
//    {
//        cout<<node_en_old.second<<endl;
//        cout<<"can't find two"<<endl;
//    }
//    
    
    row_remove = time_start[node_st_old];
    col_remove = time_end[node_en_old];

    
    S_insert=interpolate_delta(node_start_insert, node_end_insert);
    
    if(nsize == 1 ){ S_tilde_insert = 1.0/S_insert; return S_insert/N_inv(0,0);}
    
    R_insert.resize(nsize-1);
    Q_insert.resize(nsize-1);
    for( map< pair<int,double> ,int>::iterator it=time_start.begin(); it!=time_start.end(); it++)
    {
        if( it->second == row_remove )
            continue;
        else if( it->second == nsize-1)
            Q_insert[row_remove] = interpolate_delta( it->first, node_end_insert );
        else
            Q_insert[it->second] = interpolate_delta( it->first, node_end_insert );
            
    }
    for( map< pair<int,double> ,int>::iterator it=time_end.begin(); it!=time_end.end(); it++)
    {
        if( it->second == col_remove )
            continue;
        else if(it->second == nsize-1)
            R_insert[col_remove] = interpolate_delta( node_start_insert, it->first);
        else
            R_insert[it->second] = interpolate_delta( node_start_insert, it->first);
    }
    
    aux_mat = N_;
    row( aux_mat, col_remove).swap( row(aux_mat,nsize-1) );
    column( aux_mat, row_remove).swap( column(aux_mat,nsize-1) );
    Q_tilde_remove = column( aux_mat, nsize-1 );
    R_tilde_remove = row( aux_mat, nsize-1 );
    S_tilde_remove = aux_mat(nsize-1,nsize-1);
    
    aux_mat.resize(nsize-1, nsize-1,true);
    Q_tilde_remove.resize(nsize-1,true);
    R_tilde_remove.resize(nsize-1,true);
    
  
    aux_mat = aux_mat - outer_prod(Q_tilde_remove,R_tilde_remove)/S_tilde_remove;
    
    
    
    aux_vec = prod(aux_mat,Q_insert);// temp vecs
    S_tilde_insert = inner_prod(R_insert,aux_vec);
    S_tilde_insert = 1.0/(S_insert-S_tilde_insert );
    

    
    return S_tilde_remove/S_tilde_insert;
    
}

fermi_type fermi_bath_block::try_replace( pair<int,double> &node_old, pair<int,double> &node_new, bool if_st)
{// if_st = true replace start;

    if( if_st)
    {
        pair<int,double> node_old_two = time_end.begin()->first;
        pair<int,double> node_new_two = node_old_two;
        return try_replace( node_old, node_new,node_old_two, node_new_two);
    }
    else
    {
        pair<int,double> node_old_two = time_start.begin()->first;
        pair<int,double> node_new_two = node_old_two;
        return try_replace( node_old_two, node_new_two,node_old, node_new);
    }
}

bool fermi_bath_block::do_replace()
{
   
    if( replace_lock ) return false;
    replace_lock = true;
    
    using namespace boost::numeric::ublas;
    
    //cout<<"nsize ="<<nsize<<" time size="<<time_start.size()<<" "<<time_end.size()<<endl;
    time_start.erase(node_start_remove);
    time_end.erase(node_end_remove);
    
    //cout<<"nsize ="<<nsize<<" time size="<<time_start.size()<<" "<<time_end.size()<<endl;
    
    time_start.insert( make_pair(node_start_insert, row_remove));
    time_end.insert(make_pair(node_end_insert, col_remove));
    
    if(nsize == 1 )
    {
        N_inv(0,0)=S_insert;
        N_(0,0)= S_tilde_insert;
        
        det = S_insert;
        
        return true;
    }
    
    
    det = det * S_tilde_remove/S_tilde_insert;
    
    
    Q_tilde_insert = aux_vec*(-S_tilde_insert);
    R_tilde_insert = prod(R_insert,aux_mat)*(-S_tilde_insert);
    aux_mat = aux_mat + outer_prod(Q_tilde_insert, R_tilde_insert)/(S_tilde_insert);
    project(N_,range(0,nsize-1),range(0,nsize-1)) = aux_mat;
    Q_tilde_insert.resize(nsize,true);
    R_tilde_insert.resize(nsize,true);
    
    row(N_,nsize-1)=R_tilde_insert;
    column(N_, nsize-1)=Q_tilde_insert;
    N_(nsize-1,nsize-1) = S_tilde_insert;
    row( N_, col_remove ).swap( row(N_, nsize-1) );
    column( N_, row_remove ).swap( column(N_, nsize-1) );
    
    
    row( N_inv, row_remove ).swap( row(N_inv, nsize-1) );
    column( N_inv, col_remove ).swap( column(N_inv, nsize-1) );
    
    
    
    R_insert.resize(nsize,true);
    Q_insert.resize(nsize,true);
    row(N_inv,nsize-1)  = R_insert;
    column(N_inv,nsize-1) = Q_insert;
    N_inv(nsize-1,nsize-1) = S_insert;
    row( N_inv, row_remove ).swap( row(N_inv, nsize-1) );
    column( N_inv, col_remove ).swap( column(N_inv, nsize-1) );
    
    
    return true;
    
}









//------------------------------------------------------------debug---------------------------------------------------------------------------------------------------------
double cal_det(boost::numeric::ublas::matrix<double> m)
{
    assert(m.size1() == m.size2() && "Can only calculate the determinant of square matrices");
    boost::numeric::ublas::permutation_matrix<std::size_t> pivots(m.size1() );
    
    const int is_singular = boost::numeric::ublas::lu_factorize(m, pivots);
    
    if (is_singular) return 0.0;
    
    double d = 1.0;
    const std::size_t sz = pivots.size();
   
    for (std::size_t i=0; i != sz; ++i)
    {
        
        if (pivots(i) != i)
        {
            d *= -1.0;
        }
       
        d *= m(i,i);
    }
   
    return d;
}


bool fermi_bath_block::check_sanity()
{
    int i,j;
    pair<int,double> t_start,t_end;
    
    bool if_sanity=true;
    
    if( nsize != static_cast<int>(time_start.size()) || nsize != static_cast<int>(time_end.size()) )
    {
        cerr<<"size of fermi time array is wrong"<<endl;
        cerr<<" nsize, size(time_start), size(time_end): "<<endl;
        cerr<<nsize<<"  "<<time_start.size()<<" "<<time_end.size()<<endl;
        return false;
    }
    if( N_.size1() != nsize || N_.size2() != nsize || N_inv.size1() != nsize || N_inv.size2() != nsize )
    {
        cerr<<"size of fermi mat is wrong"<<endl;
        cerr<<" nsize, size(N_) size(N_inv): "<<nsize<<N_.size1()<<N_.size2()<<N_inv.size1()<<N_inv.size2()<<endl;
        return false;
    }
    
    if( nsize ==0 )
    {
        if(  det != 1.0)
        {
            cerr<<"nsize 0 wrong ";
            cerr<<" det="<<det;
            return false;
        }
        return true;
    }
    
    aux_mat = N_inv;
    if( std::abs( det - static_cast<fermi_type> ( cal_det(aux_mat) ) ) >0.000001)
    {
        
        cerr<<"determinant of fermi bath is wrong"<<endl;
        cerr<<det<<" "<<cal_det(N_inv)<<endl;
        if_sanity=false;
    }
    
    
    
    
    for(map< pair<int,double> ,int>::iterator it_st = time_start.begin(); it_st != time_start.end(); it_st ++ )
    {
        t_start = it_st->first;
        i = it_st->second;
        for(map< pair<int,double> ,int>::iterator it_en = time_end.begin(); it_en != time_end.end(); it_en ++ )
        {
            t_end = it_en -> first;
            j = it_en->second;
            
            if( std::abs( interpolate_delta(t_start, t_end)  - N_inv(i,j)) >0.000001)
            {
                cerr<<"error, N_inv element wrong"<<endl;
                cerr<<i<<" "<<j<<endl;
                cerr<< N_inv(i,j)<<" "<<interpolate_delta(t_start, t_end) <<endl;
                if_sanity=false;
            }
            
        }
    }
    
    aux_mat = prod( N_inv, N_);
    double err=0.0;
    for( i=0;i<aux_mat.size1();i++)
    {
        for(j=0;j<aux_mat.size2();j++)
        {
            if(i==j)
            {
                err += std::abs( aux_mat(i,j)-1.0);
            }
            else
            {
                err += std::abs(aux_mat(i,j));
            }
        }
    }
    
    if(err/max(1,nsize*nsize)>0.00001)
    {
        cerr<<"error inv"<<endl;
        cerr<<err<<endl;
        for(int i=0; i<aux_mat.size2();i++)
        {
            for( int j=0; j<aux_mat.size2();j++)
            {
    
                if( i==j)
                {
                    if( std::abs( aux_mat(i,j)-1.0) > 0.00001)
                        cerr<<i<<" "<<j<<" "<<std::abs( aux_mat(i,j)-1.0)<<endl;
                }
                else
                {
                    if( std::abs( aux_mat(i,j)) > 0.00001)
                        cerr<<i<<" "<<j<<" "<<std::abs( aux_mat(i,j))<<endl;
                }
            }
        }
        if_sanity=false;
    }
    
    for(map< pair<int,double> ,int>::iterator it=time_start.begin();it!=time_start.end();it++)
    {
        if( index_field.find( (it->first).first ) == index_field.end() )
        {
            cerr<<"index not in the fild, time st array"<<endl;
            cerr<<(it->first).first<<endl;
            if_sanity=false;
        }
    }
    for(map< pair<int,double> ,int>::iterator it=time_end.begin();it!=time_end.end();it++)
    {
        if( index_field.find( (it->first).first ) == index_field.end() )
        {
            cerr<<"index not in the fild, time en array"<<endl;
            cerr<<(it->first).first<<endl;
            if_sanity=false;
        }
    }
    
    return if_sanity;
    
}

