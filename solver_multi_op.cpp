//
//  solver_multi_op.cpp
//  cthyb
//
//  Created by 胡昊昱 on 8/13/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#include "solver.hpp"
fermi_type solver::try_insert_multi_op_to_fermi_bath( vector<int> &orb, vector<double> &ts, vector<double> &te, int &sign_change)
{// let sign_change_this_function be the sign change from this insertion, then sign_change = sign_change * sign_change_this_function
    int n_op = static_cast<int>(orb.size());
    double re= 1.0;
    if(n_op == 2)
    {// insert two op
        int bk1=orbital_to_block[orb[0]];
        int bk2=orbital_to_block[orb[1]];
        if( bk1!=bk2)
        {// different block
            
            re *= fermi_bath[bk1].try_insert(orb[0], orb[0], ts[0], te[0] );
            re *= fermi_bath[bk2].try_insert(orb[1], orb[1], ts[1], te[1] );
            if( ts[0]>te[0]){ sign_change = -sign_change;}
            if( ts[1]>te[1]){ sign_change = -sign_change;}
            
        }
        else
        {
            cerr<<"insert multi op!"<<endl;
            cerr<<"havn's finished "<<endl;
            exit(9);
        }
    }
    else if( n_op==4)
    {
       
        vector<int> temp_orb(2,0);
        vector<double> temp_ts(2,0);
        vector<double> temp_te(2,0);
        if( orb[0]%2 == orb[1]%2 )
        {// 1,2 same spin; 3,4 same spin
            
            temp_orb={orb[0],orb[1]};
            temp_ts={ts[0],ts[1]};
            temp_te={te[0],te[1]};
            re *= try_insert_multi_op_to_fermi_bath(temp_orb,  temp_ts, temp_te, sign_change);
            temp_orb={orb[2],orb[3]};
            temp_ts={ts[2],ts[3]};
            temp_te={te[2],te[3]};
            re *= try_insert_multi_op_to_fermi_bath(temp_orb,  temp_ts, temp_te, sign_change);
        }
        else if( orb[0]%2 == orb[2]%2 )
        {// 1,3 same spin; 2,4 same spin
            
            temp_orb={orb[0],orb[2]};
            temp_ts={ts[0],ts[2]};
            temp_te={te[0],te[2]};
            re *= try_insert_multi_op_to_fermi_bath(temp_orb,  temp_ts, temp_te, sign_change);
            temp_orb={orb[1],orb[3]};
            temp_ts={ts[1],ts[3]};
            temp_te={te[1],te[3]};
            re *= try_insert_multi_op_to_fermi_bath(temp_orb, temp_ts, temp_te, sign_change);
        }
        else
        {// 1,4 same spin; 2,3 same spin
            temp_orb={orb[0],orb[3]};
            temp_ts={ts[0],ts[3]};
            temp_te={te[0],te[3]};
            re *= try_insert_multi_op_to_fermi_bath(temp_orb,  temp_ts, temp_te, sign_change);
            temp_orb={orb[1],orb[2]};
            temp_ts={ts[1],ts[2]};
            temp_te={te[1],te[2]};
            re *= try_insert_multi_op_to_fermi_bath(temp_orb, temp_ts, temp_te, sign_change);
        }
        
        
    }
    else
    {
        cerr<<" try to insert "<<n_op<<" ops"<<endl;
        exit(1);
    }
    return re;
    
    
    
    
}



fermi_type solver::try_remove_multi_op_to_fermi_bath( vector<int> &orb, vector<double> &ts, vector<double> &te, int &sign_change)
{
    int n_op = static_cast<int>(orb.size());
    double re= 1.0;
    if(n_op == 2)
    {// two op
        int bk1=orbital_to_block[orb[0]];
        int bk2=orbital_to_block[orb[1]];
        if( bk1!=bk2)
        {// different block
    
            re *= fermi_bath[bk1].try_remove(orb[0], orb[0], ts[0], te[0] );
            re *= fermi_bath[bk2].try_remove(orb[1], orb[1], ts[1], te[1] );
            if( ts[0]>te[0]){ sign_change = -sign_change;}
            if( ts[1]>te[1]){ sign_change = -sign_change;}
            
        }
        else
        {
            cerr<<"remove multi op!"<<endl;
            cerr<<"havn't finished "<<endl;
            exit(9);
        }
    }
    else if( n_op==4)
    {
       
        vector<int> temp_orb(2,0);
        vector<double> temp_ts(2,0);
        vector<double> temp_te(2,0);
       
        if( orb[0]%2 == orb[1]%2 )
        {// 1,2 same spin; 3,4 same spin
            
            temp_orb={orb[0],orb[1]};
            temp_ts={ts[0],ts[1]};
            temp_te={te[0],te[1]};
            re *= try_remove_multi_op_to_fermi_bath(temp_orb,  temp_ts, temp_te, sign_change);
            temp_orb={orb[2],orb[3]};
            temp_ts={ts[2],ts[3]};
            temp_te={te[2],te[3]};
            re *= try_remove_multi_op_to_fermi_bath(temp_orb,  temp_ts, temp_te, sign_change);
        }
        else if( orb[0]%2 == orb[2]%2 )
        {// 1,3 same spin; 2,4 same spin
            
            temp_orb={orb[0],orb[2]};
            temp_ts={ts[0],ts[2]};
            temp_te={te[0],te[2]};
            re *= try_remove_multi_op_to_fermi_bath(temp_orb,  temp_ts, temp_te, sign_change);
            temp_orb={orb[1],orb[3]};
            temp_ts={ts[1],ts[3]};
            temp_te={te[1],te[3]};
            re *= try_remove_multi_op_to_fermi_bath(temp_orb, temp_ts, temp_te, sign_change);
        }
        else
        {// 1,4 same spin; 2,3 same spin
            temp_orb={orb[0],orb[3]};
            temp_ts={ts[0],ts[3]};
            temp_te={te[0],te[3]};
            re *= try_remove_multi_op_to_fermi_bath(temp_orb,  temp_ts, temp_te, sign_change);
            temp_orb={orb[1],orb[2]};
            temp_ts={ts[1],ts[2]};
            temp_te={te[1],te[2]};
            re *= try_remove_multi_op_to_fermi_bath(temp_orb, temp_ts, temp_te, sign_change);
        }
        
        
        
    }
    else
    {
        cerr<<" try to remove "<<n_op<<" ops"<<endl;
        exit(1);
    }
    return re;
    
}

fermi_type solver::try_replace_multi_op_to_fermi_bath( vector<int> &orb, vector<bool> &dag,vector<double> &told, vector<double> &tnew,int &sign_change)
{
    
    int n_op = static_cast<int>(orb.size());
    double re= 1.0;
    if(n_op == 2)
    {// insert two op
        int bk1=orbital_to_block[orb[0]];
        int bk2=orbital_to_block[orb[1]];
//        cout<<"orb1,2"<<orb[0]<<" "<<orb[1]<<endl;
//        cout<<"bk"<<bk1<<" "<<bk2<<endl;
        if( bk1!=bk2)
        {// different block
            
            re *= fermi_bath[bk1].try_replace(orb[0],told[0], tnew[0] ,dag[0]);
            re *= fermi_bath[bk2].try_replace(orb[1], told[1], tnew[1], dag[1]);
        }
        else
        {
            throw 1;
            cerr<<"replace multi op!"<<endl;
            cerr<<"havn't finished "<<endl;
            exit(9);
        }
    }
    else if( n_op==4)
    {
        
        vector<int> temp_orb(2,0);
        vector<double> temp_old(2,0);
        vector<double> temp_new(2,0);
        vector<bool> temp_dag(2,true);
        if( orb[0]%2 == orb[1]%2 )
        {// 1,2 same spin; 3,4 same spin
            
            temp_orb={orb[0],orb[1]};
            temp_old={told[0],told[1]};
            temp_new={tnew[0],tnew[1]};
            temp_dag={dag[0],dag[1]};
            re *= try_replace_multi_op_to_fermi_bath(temp_orb,temp_dag,  temp_old, temp_new,sign_change);
            temp_orb={orb[2],orb[3]};
            temp_old={told[2],told[3]};
            temp_new={tnew[2],tnew[3]};
            temp_dag={dag[2],dag[3]};
            re *= try_replace_multi_op_to_fermi_bath(temp_orb,temp_dag,  temp_old, temp_new,sign_change);
        }
        else if( orb[0]%2 == orb[2]%2 )
        {// 1,3 same spin; 2,4 same spin
            
            temp_orb={orb[0],orb[2]};
            temp_old={told[0],told[2]};
            temp_new={tnew[0],tnew[2]};
            temp_dag={dag[0],dag[2]};
            re *= try_replace_multi_op_to_fermi_bath(temp_orb,temp_dag,  temp_old, temp_new,sign_change);
            temp_orb={orb[1],orb[3]};
            temp_old={told[1],told[3]};
            temp_new={tnew[1],tnew[3]};
            temp_dag={dag[1],dag[3]};
            re *= try_replace_multi_op_to_fermi_bath(temp_orb,temp_dag,  temp_old, temp_new,sign_change);
        }
        else
        {// 1,4 same spin; 2,3 same spin
            temp_orb={orb[0],orb[3]};
            temp_old={told[0],told[3]};
            temp_new={tnew[0],tnew[3]};
            temp_dag={dag[0],dag[3]};
            re *= try_replace_multi_op_to_fermi_bath(temp_orb,temp_dag,  temp_old, temp_new,sign_change);
            temp_orb={orb[1],orb[2]};
            temp_old={told[1],told[2]};
            temp_new={tnew[1],tnew[2]};
            temp_dag={dag[1],dag[2]};
            re *= try_replace_multi_op_to_fermi_bath(temp_orb,temp_dag,  temp_old, temp_new,sign_change);
        }
        
        
        
    }
    else
    {
        cerr<<" try to remove "<<n_op<<" ops"<<endl;
        exit(1);
    }
    return re;
    
}
    

bool solver::do_insert_multi_op_to_fermi_bath(int &sign_change, vector<int> &orb ,bool if_sucess)
{
    
    int n_op = static_cast<int>(orb.size());
    if(n_op ==2 )
    {
        int bk1=orbital_to_block[orb[0]];
        int bk2=orbital_to_block[orb[1]];
        if( bk1 != bk2 )
        {
            if( if_sucess )
            {
                fermi_bath[bk1].do_insert();
                fermi_bath[bk2].do_insert();
            }
            else
            {
                fermi_bath[bk1].lock_insert();
                fermi_bath[bk2].lock_insert();
            }
        }
        else
        {
            cerr<<" do insert multi op "<<endl;
            cerr<<"havn't finished yet"<<endl;
        }
    }
    else if( n_op==4)
    {
       
        vector<int> temp_orb(2,0);
        if( orb[0]%2 == orb[1]%2 )
        {// 1,2 same spin; 3,4 same spin
            
            temp_orb={orb[0],orb[1]};
            do_insert_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
        
            temp_orb={orb[2],orb[3]};
            do_insert_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
        }
        else if( orb[0]%2 == orb[2]%2 )
        {// 1,3 same spin; 2,4 same spin
            
            temp_orb={orb[0],orb[2]};
            do_insert_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
            
            temp_orb={orb[1],orb[3]};
            do_insert_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
        }
        else
        {// 1,4 same spin; 2,3 same spin
            temp_orb={orb[0],orb[3]};
            do_insert_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
            
            temp_orb={orb[1],orb[2]};
            do_insert_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
        }
        
        
    }
    else
    {
        cerr<<" try to insert "<<n_op<<" ops"<<endl;
        exit(1);
    }
    return true;
}
bool solver::do_remove_multi_op_to_fermi_bath(int &sign_change, vector<int> &orb, bool if_sucess)
{
    
    int n_op = static_cast<int>(orb.size());
    int temp_sign;
    if(n_op ==2 )
    {
        int bk1=orbital_to_block[orb[0]];
        int bk2=orbital_to_block[orb[1]];
        if( bk1 != bk2 )
        {
            if( if_sucess )
            {
                temp_sign = 1;
                fermi_bath[bk1].do_remove(temp_sign);
                sign_change *= temp_sign;
                temp_sign = 1;
                fermi_bath[bk2].do_remove(temp_sign);
                sign_change *= temp_sign;
            }
            else
            {
                fermi_bath[bk1].lock_remove();
                fermi_bath[bk2].lock_remove();
            }
        }
        else
        {
            cerr<<" do remove multi op "<<endl;
            cerr<<"havn't finished yet"<<endl;
        }
    }
    else if( n_op==4)
    {
    
        vector<int> temp_orb(2,0);
        if( orb[0]%2 == orb[1]%2 )
        {// 1,2 same spin; 3,4 same spin
            
            temp_orb={orb[0],orb[1]};
            do_remove_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
            
            temp_orb={orb[2],orb[3]};
            do_remove_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
        }
        else if(orb[0]%2 == orb[2]%2)
        {// 1,3 same spin; 2,4 same spin
            
            temp_orb={orb[0],orb[2]};
            do_remove_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
            
            temp_orb={orb[1],orb[3]};
            do_remove_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
        }
        else
        {// 1,4 same spin; 2,3 same spin
            temp_orb={orb[0],orb[3]};
            do_remove_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
            
            temp_orb={orb[1],orb[2]};
            do_remove_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
        }
        
        
    }
    else
    {
        cerr<<" try to remove "<<n_op<<" ops"<<endl;
        exit(1);
    }
    return true;
}
bool solver::do_replace_multi_op_to_fermi_bath(int &sign_change, vector<int> &orb,bool if_sucess)
{
    
    int n_op = static_cast<int>(orb.size());
    if(n_op ==2 )
    {
        int bk1=orbital_to_block[orb[0]];
        int bk2=orbital_to_block[orb[1]];
        if( bk1 != bk2 )
        {
            if( if_sucess )
            {
                fermi_bath[bk1].do_replace();
                fermi_bath[bk2].do_replace();
            }
            else
            {
                fermi_bath[bk1].lock_replace();
                fermi_bath[bk2].lock_replace();
            }
        }
        else
        {
            cerr<<" do replace multi op "<<endl;
            cerr<<"havn't finished yet"<<endl;
            }
    }
    else if( n_op==4)
    {
       
        vector<int> temp_orb(2,0);
        if( orb[0]%2 == orb[1]%2 )
        {// 1,2 same spin; 3,4 same spin
            
            temp_orb={orb[0],orb[1]};
            do_replace_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
            
            temp_orb={orb[2],orb[3]};
            do_replace_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
        }
        else if( orb[0]%2 == orb[2]%2 )
        {// 1,3 same spin; 2,4 same spin
            
            temp_orb={orb[0],orb[2]};
            do_replace_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
            
            temp_orb={orb[1],orb[3]};
            do_replace_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
        }
        else
        {// 1,4 same spin; 2,3 same spin
            temp_orb={orb[0],orb[3]};
            do_replace_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
            
            temp_orb={orb[1],orb[2]};
            do_replace_multi_op_to_fermi_bath( sign_change,temp_orb,  if_sucess);
        }
        
        
    }
    else
    {
        cerr<<" try to replace "<<n_op<<" ops"<<endl;
        exit(1);
    }
    return true;
}
