//
//  main.cpp
//  cthyb
//
//  Created by 胡昊昱 on 2019/7/4.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#include <iostream>
#include <set>
#include <utility>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "fermi_bath.hpp"
#include "mat_vec.hpp"
#include "random.hpp"


using namespace boost::numeric::ublas;
using namespace boost::qvm;


int test_main_fm_bath ()
{
    int norb_in=1;
    double beta_in = 100.0;
    std::set<int> index_field_in={0};
    std::vector< std::vector< fermi_type > >  delta_in(1,std::vector< fermi_type >(0) );
    int N_t = 100;
    
    delta_in[0].resize(N_t+1);
    
    double nf = 0.0;
    double de,t,ep;
    int N_e = 10;
    
    de = 2.0/double(N_e);
    for( int i=0;i<N_e;i++)
    {
        ep = -1.0 + i*de;
        nf = 1.0/(exp(beta_in*ep)+1);
        
        for( int it =0 ; it<N_t+1;it++)
        {
            t= it/double(N_t) * beta_in;
            delta_in[0][it] +=nf/2.0 * exp(t*ep) * de;
        }
    }
    for( int it =0 ; it<N_t+1;it++)
    {
        t= it/double(N_t) * beta_in;
        //cout<<t<<" "<<delta_in[0][it] <<endl;
    }
    
    fermi_bath_block test_fm(norb_in,  beta_in, index_field_in);
    test_fm.input_delta(delta_in);
    
    rand_gen rand(10);
    double beta = beta_in;
    
    std::vector< pair<int,double> > time_st(0);
    std::vector< pair<int,double> > time_en(0);
    
    
    for( int i=0;i<5;i++)
    {
        cout<<i<<endl;
        pair<int,double> node_st(0, (2.0*i)/20.0 * beta );
        pair<int,double> node_en(0, (2.0*i+0.1)/20.0 * beta );
//        pair<int,double> node_st(0, rand.random_double(beta) );
//        pair<int,double> node_en(0, rand.random_double(beta)  );
        
        time_st.push_back(node_st);
        time_en.push_back(node_en);
        double ratio = test_fm.try_insert(node_st, node_en);
        cout<<"weight "<<test_fm.get_det()<<" "<<ratio<<" "<<test_fm.get_det()*ratio<<" ";
        test_fm.do_insert();
        cout<<test_fm.get_det()<<endl;
        if( not test_fm.check_sanity() )
        {
            exit(100);
        }
        
    }
    
    for( int i=0;i<2;i++)
    {
        cout<<i<<endl;
        int n = static_cast<int>(time_st.size());
        int k1 = rand.random_int(n);
        int k2 = rand.random_int(n);
//        int k1 = (n-3<0)? 0: n-3;
//        int k2 = (n-1<0)? 0: n-1;
        pair<int,double> node_st = time_st[ k1 ];
        pair<int,double> node_en = time_en[ k2 ];
        time_st.erase(time_st.begin() + k1);
        time_en.erase(time_en.begin() + k2);
        
        //cout<<"weight ratio  "<<test_fm.try_remove(node_st, node_en)<<" "<<test_fm.get_det()<<endl;
        double ratio = test_fm.try_remove(node_st,node_en);
        cout<<"weight "<<test_fm.get_det()<<" "<<ratio<<" "<<test_fm.get_det()*ratio<<" ";
        
        test_fm.do_remove();
        if( not test_fm.check_sanity() )
        {
            exit(100);
        }
        cout<<test_fm.get_det()<<endl;
        
    }
    cout<<"replace"<<endl;
    for( int i=0;i<3;i++)
    {
        cout<<i<<endl;
        int n = static_cast<int>(time_st.size());
        int k1 = rand.random_int(n);
        int k2 = rand.random_int(n);
        //        int k1 = (n-3<0)? 0: n-3;
        //        int k2 = (n-1<0)? 0: n-1;
        pair<int,double> node_st = time_st[ k1 ];
        pair<int,double> node_en = time_en[ k2 ];
        pair<int,double> node_st_new = make_pair( 0, rand.random_double(beta));
        pair<int,double> node_en_new = make_pair( 0, rand.random_double(beta));
        
        time_st.erase(time_st.begin() + k1);
        time_en.erase(time_en.begin() + k2);
        
        time_st.push_back(node_st_new);
        time_en.push_back(node_en_new);
        
        double ratio = test_fm.try_replace(node_st, node_st_new, node_en, node_en_new);
        cout<<"weight "<<test_fm.get_det()<<" "<<ratio<<" "<<test_fm.get_det()*ratio<<" ";
        test_fm.do_replace();
        cout<<test_fm.get_det()<<endl;
        if( not test_fm.check_sanity() )
        {
            exit(100);
        }
        
        
    }
    
    return 0;
    

}
