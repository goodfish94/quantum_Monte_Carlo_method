//
//  main.cpp
//  cthyb
//
//  Created by 胡昊昱 on 7/11/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#include <stdio.h>
#include <iostream>

#include "local_config.hpp"
#include "random.hpp"

using namespace std;

void insert_test(rand_gen& rand, local_config& test )
{
    double beta = test.get_beta();
    int ntype =test.get_num_of_type();
    bool st_dag;
    int type = rand.random_int(ntype);
    if( rand.random_double(1.0)<0.5)
    {
        st_dag = true;
    }
    else
    {
        st_dag = false;
    }
    op st_node( rand.random_double(beta), st_dag, type,0);
    double lmax = test.try_insert( st_node );
    double t_new = st_node.t;
    if( lmax>0.0)
    {
        t_new+= rand.random_double(lmax);
    }
    else
    {
        cout<<"lmax = 0"<<endl;
        return;
    }
    
    if(t_new>beta){t_new-= beta;}
    
    op en_node(t_new,!st_dag,type,0);
    test.insert( st_node, en_node );
    
    bool sanity = test.check_sanity();
    cout<<"sanity="<<sanity<<endl;
    test.print_config();
    if( sanity == false)
    {
        cout<<"lmax = "<<lmax<<endl;
        cout<<"opst = "<<st_node.t<<" "<<st_node.dag<<" "<<st_node.type<<endl;
        cout<<"open = "<<en_node.t<<" "<<en_node.dag<<" "<<en_node.type<<endl;
        exit(0);
    }
}

void remove_test(rand_gen& rand, local_config& test )
{
    double beta = test.get_beta();
    int ntype =test.get_num_of_type();
    int type = rand.random_int(ntype);
    bool st_dag;
    if( rand.random_double(1.0)<0.5)
    {
        st_dag = true;
    }
    else
    {
        st_dag = false;
    }
    int k_remove = rand.random_int( test.get_seg_num(type) );
    op node_one,node_two;
    double lmax = test.try_remove(type, k_remove, st_dag, node_one, node_two );
    if( lmax<=0.0)
    {
        cout<<"lmax="<<lmax<<endl;
        return;
    }
    
    cout<<"remove="<<test.remove()<<endl;;
    bool sanity = test.check_sanity();
    cout<<"sanity="<<sanity<<endl;
    test.print_config();
    if( sanity == false)
    {
        cout<<"lmax = "<<lmax<<endl;
        cout<<"op one = "<<node_one.t<<" "<<node_one.dag<<" "<<node_one.type<<endl;
        cout<<"op two= "<<node_two.t<<" "<<node_two.dag<<" "<<node_two.type<<endl;
        exit(0);
    }
}


void shift_test(rand_gen& rand, local_config& test )
{
    double beta = test.get_beta();
    int ntype =test.get_num_of_type();
    op modified_node;
    double up_lim,down_lim;
    bool if_dag;
    int type = rand.random_int(ntype);
    int k_shift = rand.random_int( test.get_seg_num(type) );
    if( rand.random_double(1.0)<0.5)
    {
        if_dag = true;
    }
    else
    {
        if_dag = false;
    }
    
    test.try_shift(type, k_shift,if_dag,down_lim, up_lim, modified_node  );

    double delta_tau = down_lim + rand.random_double(up_lim-down_lim);
    
    test.shift( delta_tau);
    bool sanity = test.check_sanity();
    test.print_config();
    if( sanity == false)
    {
        cout<<"up, down dtau= "<<up_lim<<" "<<down_lim<<"  "<<delta_tau<<endl;
        cout<<"old node = "<<modified_node .t<<" "<<modified_node .dag<<" "<<modified_node .type<<endl;
        exit(0);
    }
}

int test_local_config()
{
    double beta = 100.0;
    int ntype = 4;
    local_config test( ntype, beta);
    rand_gen rand(10);
    
    
    for(int i=0;i<1000;i++)
    {
        double num = rand.random_double(1.0);
        if( num <0.25)
        {
            cout<<"insert------------------"<<i<<"------------------"<<endl;
            insert_test( rand, test );

        }
        else if(num<0.5)
        {
            cout<<"remove------------------"<<i<<"------------------"<<endl;
            remove_test( rand, test );
            
        }
        else
        {
            cout<<"shift------------------"<<i<<"------------------"<<endl;
            shift_test( rand, test );
            
            
        }
    }
    return 1;
}
