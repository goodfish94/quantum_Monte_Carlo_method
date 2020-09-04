//
//  local_config.cpp
//  cthyb
//
//  Created by 胡昊昱 on 2019/7/6.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#include "local_config.hpp"

double local_config::try_insert( op &new_node )
{
    if( n_ == 0 )
    {
        if( new_node.dag != if_full )
        {
            insert_lock=false;
            insert_hint = config.begin();
            return beta;
        }
        else{ return 0.0;}
    }
    insert_hint = config.lower_bound(new_node);// find iterator to the element right after new_node;
    if( insert_hint != config.end() )
    {
        if( insert_hint->dag == new_node.dag )
        {
            insert_lock=false;
            return (insert_hint->t - new_node.t);
        }
        else{ return 0.0; }
    }
    else
    {
        if( if_wrap == (! new_node.dag) )
        {
            insert_lock =false;
            return beta-new_node.t+(config.begin()->t);
        }
        else{return 0.0;}
    }
    
}


bool local_config::insert(op &node_st , op &node_en)
{
    if(insert_lock){ return false;}
    else{ insert_lock =true;}
    if( n_==0 ){if_full = false;}
    insert_hint = config.insert(insert_hint, node_st);//! warning: use hint to insert, useful only for C++11;
    insert_hint++;
    if( node_st.t>node_en.t)
    {
        config.insert(config.begin(), node_en);//! warning: use hint to insert, useful only for C++11;
    }
    else
    {
        config.insert(insert_hint, node_en);//! warning: use hint to insert, useful only for C++11;
    }
    n_ = n_ +2;
    n_type[node_st.type][node_st.subtype]=n_type[node_st.type][node_st.subtype]+2;
    update_wrap();
    return true;
}




double local_config::try_remove(int type_remove, int subtype_remove, int k_remove, bool if_anti_seg, op &node_one, op&node_two  )
{// -1 if can't remove
    //try remove  k-th  given-type element, k\in 0, 1, 2,...
    // return corresponding lmax
    int k = k_remove;
    double lmax = 0.0;
    if( n_type[type_remove][subtype_remove]==0) return -1;
    
    for( remove_hint =config.begin(); remove_hint != config.end(); remove_hint++)
    {
        if( remove_hint->type == type_remove && remove_hint->subtype == subtype_remove && (if_anti_seg!=(remove_hint->dag)) )
        {
            if( k== 0){ break;}
            k--;
        }
    }
    config_type::iterator it =remove_hint;
    it++;
    if( it== config.end() ){ lmax = lmax+beta; it = config.begin();}
    if( it->type != remove_hint->type || it->subtype != remove_hint->subtype ) { return -1;}
    node_one = *remove_hint;
    node_two = *it;
    node_remove_1 = node_one;
    node_remove_2 = node_two;
    remove_lock=false;
    
    it++;
    if( it== config.end() ){ lmax = lmax+beta; it = config.begin();}
    lmax = lmax  + it->t - remove_hint->t;

    
    return lmax;
    
}


double local_config::try_remove( op node_one, op node_two  )
{// -1 if can't remove
    //try remove node_one, node_two ( keep order )
    // return corresponding lmax
    double lmax = 0.0;
    if( n_type[node_one.type][node_one.subtype]==0) return -1;
    remove_hint = config.find(node_one);
    if( remove_hint == config.end() ){ return -1;}
    config_type::iterator it =remove_hint;
    it++;
    if( it == config.end() ){ it = config.begin(); lmax = beta;}
    if( *it != node_two ) {return -1;}
    it ++;
    if( it == config.end() ){ it = config.begin(); lmax = beta;}
    lmax += it->t - node_one.t;

    
    remove_lock = false;
    return lmax;
    
    
}
bool local_config::enforce_remove(op &node_1, op &node_2)
{// remove node;
    remove_hint = config.find( node_1 );
    node_remove_1 =*remove_hint;
    config_type::iterator it = remove_hint;
    if( it== config.end() ){ it = config.begin();}
    node_remove_2 = *(it);
    if( node_remove_2!= node_2){ return false; }
    remove_lock = false;
    return remove();
    
}
bool local_config::remove()
{/// remove node according to remove hint
    if( remove_lock ){return false;}
    else{ remove_lock = true;}
    
    n_ = n_ -2;
    n_type[ remove_hint->type][remove_hint->subtype]= n_type[ remove_hint->type][remove_hint->subtype] - 2;
    
    if( n_==0 )
    {
        if( remove_hint->dag ){/* remove seg */ }
        else{ if_full = true; }
    }
    
    remove_hint = config.erase( remove_hint );// warning: valid operation for C++11 ;
    if( remove_hint == config.end())
    {
        remove_hint =config.begin();
    }

    config.erase(remove_hint);// warning: valid operation for C++11 ;
    
    update_wrap();
    return true;
}


bool local_config::try_shift(int type_shift, int subtype_shift, int k_shift,bool if_dag,double &down_lim, double &up_lim,op &node )
{
// shift type_shift, k-th node
    if(n_type[type_shift][subtype_shift] == 0 ) return false;
    
    config_type::iterator it=config.end();
    
    down_lim = (config.begin()->t) ; up_lim = (config.begin()->t);
    
    for( shift_hint = config.begin(); shift_hint!= config.end();shift_hint++)
    {
        if(shift_hint->dag == if_dag && shift_hint->type == type_shift  && shift_hint->subtype == subtype_shift )
        {
            if( k_shift ==0 )
            {
                break;
            }
            k_shift--;
        }
        down_lim = shift_hint->t;
    }
    if( shift_hint == config.begin() )
    {
        it--;
        down_lim = -( beta - it->t + shift_hint->t );
    }
    else
    {
        down_lim = down_lim - shift_hint->t;
    }
    it = shift_hint;
    it++;
    if( it== config.end() )
    {
        up_lim = config.begin()->t + beta - shift_hint->t;
    }
    else
    {
        up_lim = it->t - shift_hint->t;
    }
    node = *(shift_hint);
    shift_lock = false;
    return true;
}

 bool local_config::try_shift(op shift_node,double &down_lim, double &up_lim)
{// return shift region t-> t+[down_lim,up_lim];
    shift_hint = config.find(shift_node);
    if( shift_hint == config.end() ){return false;}
    
    config_type::iterator it = shift_hint;
    
    if( it==config.begin() )
    {
        it = config.end();
        it--;
        down_lim = -( beta - it->t + shift_hint->t );
    }
    else
    {
        it--;
        down_lim = -( shift_hint->t - it->t);
    }
    it = shift_hint;
    it++;
    if( it == config.end() )
    {
        up_lim = beta - shift_hint->t + config.begin()->t;
    }
    else
    {
        up_lim = it->t - shift_hint->t;
    }
    shift_lock = false;
    return true;

    
}
bool local_config::shift(double delta_tau)
{
    // shift after try, delta_tau \in[ down_lim, up_lim]
    if( shift_lock )return false;
    shift_lock = true;
    
    op new_op = *(shift_hint);
    double t_new = shift_hint->t +delta_tau;
    
    shift_hint = config.erase(shift_hint);
    
    if( t_new >beta )
    {
        new_op.t =  t_new - beta;
        shift_hint = config.begin();
    }
    else if(t_new<0)
    {
        new_op.t = t_new+beta;
        shift_hint = config.end();
    }
    else
    {
        new_op.t = t_new;
    }
    config.insert(shift_hint,new_op);
    update_wrap();

    return true;

}



bool local_config::shift( op new_op)
{
    // shift after try, delta_tau \in[ down_lim, up_lim]
    if( shift_lock )return false;
    shift_lock = true;
   
    
    n_type[ shift_hint->type][shift_hint->subtype]--;
    n_type[ new_op.type][new_op.subtype] ++;
    
    shift_hint = config.erase(shift_hint);
    config.insert(shift_hint,new_op);
    update_wrap();
    
    
    return true;
    
}


double local_config::get_n_avg()
{
    if(n_==0)
    {
        if(if_full) return 1.0;
        else return 0.0;
    }
    
    double re =0.0;
    
    double t_1 =0.0, t_2 =0.0;
    for( config_type::iterator it=config.begin(); it !=config.end(); it++)
    {
        t_1 = t_2;
        t_2 = it->t;
        if( !(it->dag) ) { re += t_2 -t_1;}
        
    }
    if(if_wrap){ re += beta-t_2;}
    return re/beta;
}


op local_config::find_k_th_op( int k, bool dag, int type, int subtype)
{// find k-th op with given type dag
    find_lock = false;
    int count = k;
    op re(0.0, dag, type,subtype);
    for( find_hint = config.begin(); find_hint != config.end(); find_hint++)
    {
        if(  find_hint->type == type && find_hint->subtype == subtype && find_hint->dag == dag  )
        {
            if( count == 0)
            {
                re.t = find_hint -> t;
                return re;
            }
            count -- ;
        }
        
    }
    cerr<<" local config can't find error "<<endl;
    exit(4);
    
}

op local_config::find_k_th_op( int k, bool dag, int type, int subtype,double &down_lim, double &up_lim)
{// find k-th op with given type dag
    find_lock = false;
    int count = k;
    op re(0.0, dag, type,subtype);
    for( find_hint = config.begin(); find_hint != config.end(); find_hint++)
    {
        if(  find_hint->type == type && find_hint->subtype == subtype && find_hint->dag == dag  )
        {
            if( count == 0)
            {
                re.t = find_hint -> t;
                break;
            }
            count -- ;
        }
        
    }
    config_type::iterator it= find_hint;
    it ++;
    if( it == config.end()){ it = config.begin();up_lim = beta;}
    else{ up_lim = 0.0;}
    up_lim += it->t - find_hint->t;
    
    it = find_hint;
    if( it == config.begin() ){ down_lim = beta; it = config.end(); }
    else{ down_lim =  0.0; }
    it--;
    down_lim += find_hint->t - it->t;
    down_lim = -down_lim;
    return re;
    
}

bool local_config::find_op(op op_find)
{// try to find op
    find_hint = config.find(op_find);
    if( find_hint == config.end()){ return false;}
    
    find_lock = false;
    return true;
}
bool local_config::find_op(op op_find,double &down_lim, double &up_lim)
{
    find_hint = config.find(op_find);
    if( find_hint == config.end()){ return false;}
    config_type::iterator it= find_hint;
    it ++;
    if( it == config.end()){ it = config.begin();up_lim = beta;}
    else{ up_lim = 0.0;}
    up_lim += it->t - find_hint->t;
    
    it = find_hint;
    if( it == config.begin() ){ down_lim = beta; it = config.end(); }
    else{ down_lim =  0.0; }
    it--;
    down_lim += find_hint->t - it->t;
    down_lim = -down_lim;
    find_lock = false;
    return true;
}
bool local_config::find_op_before_t(double t, bool dag, int type, int subtype, op &op_find)
{// find the op before time t with given type dag and transform value to op_find
    if( n_type[type][subtype] == 0 ){ return false;}
    
    op index_op(t,dag,type,subtype);
    find_hint = config.lower_bound(index_op);
    
    if( find_hint == config.begin() ){ find_hint = config.end();}
    find_hint --; 
    
    while( 1 )
    {
        if( find_hint->type == type && find_hint->dag == dag && find_hint->subtype == subtype )
        {
            find_lock = false;
            op_find = *find_hint;
            return true;
        }
        if( find_hint == config.begin() ){ find_hint = config.end();}
        find_hint --;
    }
}

bool local_config::reaplace_find( op new_op)
{// replace  the find op with new one
    if( find_lock == true ){ return false;}
    else{ find_lock = true;}
    
    n_type[ find_hint->type][find_hint->subtype]--;
    n_type[ new_op.type][new_op.subtype] ++;
    

    find_hint = config.erase(find_hint);
    config.insert(find_hint,new_op);
    update_wrap();
    return true;
}



op local_config::find_k_th_op_2( int k, bool dag, int type, int subtype)
{// find k-th op with given type dag
    find_lock_2 = false;
    int count = k;
    op re(0.0, dag, type,subtype);
    for( find_hint_2 = config.begin(); find_hint_2 != config.end(); find_hint_2++)
    {
        if(  find_hint_2->type == type && find_hint_2->subtype == subtype && find_hint_2->dag == dag  )
        {
            if( count == 0)
            {
                re.t = find_hint_2 -> t;
                return re;
            }
            count -- ;
        }
        
    }
    cerr<<" local config can't find error "<<endl;
    exit(4);
    
}


op local_config::find_k_th_op_2( int k, bool dag, int type, int subtype,double &down_lim, double &up_lim)
{// find k-th op with given type dag
    find_lock_2 = false;
    int count = k;
    op re(0.0, dag, type,subtype);
    for( find_hint_2 = config.begin(); find_hint_2 != config.end(); find_hint_2++)
    {
        if(  find_hint_2->type == type && find_hint_2->subtype == subtype && find_hint_2->dag == dag  )
        {
            if( count == 0)
            {
                re.t = find_hint_2 -> t;
                break;
            }
            count -- ;
        }
        
    }
    config_type::iterator it= find_hint_2;
    it ++;
    if( it == config.end()){ it = config.begin();up_lim = beta;}
    else{ up_lim = 0.0;}
    up_lim += it->t - find_hint_2->t;
    
    it = find_hint_2;
    if( it == config.begin() ){ down_lim = beta; it = config.end(); }
    else{ down_lim =  0.0; }
    it--;
    down_lim += find_hint_2->t - it->t;
    down_lim = -down_lim;
    return re;
    
}


bool local_config::find_op_2(op op_find)
{// try to find op
    find_hint_2 = config.find(op_find);
    if( find_hint_2 == config.end()){ return false;}
    
    find_lock_2 = false;
    return true;
}
bool local_config::find_op_2(op op_find,double &down_lim, double &up_lim)
{
    find_hint_2 = config.find(op_find);
    if( find_hint_2 == config.end()){ return false;}
    
    config_type::iterator it= find_hint_2;
    it ++;
    if( it == config.end()){ it = config.begin();up_lim = beta;}
    else{ up_lim = 0.0;}
    up_lim += it->t - find_hint_2->t;
    
    it = find_hint_2;
    if( it == config.begin() ){ down_lim = beta; it = config.end(); }
    else{ down_lim =  0.0; }
    it--;
    down_lim += find_hint_2->t - it->t;
    down_lim = -down_lim;
    
    find_lock_2 = false;
    return true;
    
}

bool local_config::find_op_before_t_2(double t, bool dag, int type, int subtype, op &op_find)
{// find the op before time t with given type dag and transform value to op_find
    if( n_type[type][subtype] == 0 ){ return false;}
    
    op index_op(t,dag,type,subtype);
    find_hint_2 = config.lower_bound(index_op);
    
    if( find_hint_2 == config.begin() ){ find_hint_2 = config.end();}
    find_hint_2 --;
    
    while( 1 )
    {
        if( find_hint_2->type == type && find_hint_2->dag == dag && find_hint_2->subtype == subtype )
        {
            find_lock_2 = false;
            op_find = *find_hint_2;
            return true;
        }
        if( find_hint_2 == config.begin() ){ find_hint_2 = config.end();}
        find_hint_2 --;
    }
}

bool local_config::reaplace_find_2( op new_op)
{// replace  the find op with new one
    if( find_lock_2 == true ){ return false;}
    else{ find_lock_2 = true;}
    
    n_type[ find_hint_2->type][find_hint_2->subtype]--;
    n_type[ new_op.type][new_op.subtype] ++;
    
    
    find_hint_2 = config.erase(find_hint_2);
    config.insert(find_hint_2,new_op);
    update_wrap();
    return true;
}


//--------------------------------------debug-----------------------------------

bool local_config::check_sanity()
{
    bool if_correct = true;
    if(n_==0)
    {
        for(vector< vector<int> >::iterator it = n_type.begin(); it!= n_type.end(); it++)
        {
            for( vector< int >::iterator it2 = it->begin(); it2 != it->end(); it2 ++ )
            {
                if( *(it2) != 0)
                {
                    cerr<<" n_type num error"<<endl;
                    cerr<<"n_, n_type[i] "<< n_<<" "<<*(it2)<<endl;
                    if_correct = false;
                    return if_correct;
                }
            }
        }
        return if_correct;
    }
    if( if_full == true && n_!=0)
    {
        cerr<<"config is full, but n_= "<<n_<<endl;
        if_correct = false;
    }
    if( config.begin()->dag == if_wrap )
    {
        cerr<<"warp flag error"<<endl;
        cerr<<"dag ="<<config.begin()->dag<<" "<<"if_wrap="<<if_wrap<<endl;
        if_correct = false;
    }
    int n_sum =0;
    for(vector< vector<int> >::iterator it = n_type.begin(); it!= n_type.end(); it++)
    {
        for( vector< int >::iterator it2 = it->begin(); it2 != it->end(); it2 ++ )
        {
            n_sum += *it2;
        }
    }
    
    if( n_sum != n_)
    {
        cerr<<"sum(n_type) != n_"<<endl;
        cerr<<n_<<" "<<n_sum<<endl;
        
        cerr<<endl;
        if_correct = false;
    }
    if( config.size() != n_ )
    {
        cerr<<"opt number != n_ "<<endl;
        cerr<<"cofing.size, n_= "<<config.size()<<" "<<n_<<endl;
        if_correct = false;
    }
    
    
    bool dag_flag = config.begin()->dag;
    bool segment_flag = true;
    config_type::iterator it=config.begin();
    vector<int> n_test( n_type.size(), 0);
    
    n_test[ it->type]++;
    
    it++;
    
    for(;it!=config.end();it++)
    {
        if( it->dag == dag_flag)
        {
            cerr<<" zero trace "<<endl;
            segment_flag = false;
        }
        n_test[it->type]++;
        dag_flag = it->dag;
    }
    for( int i=0;i<n_test.size();i++)
    {
        if( n_test[i] != accumulate( n_type[i].begin(), n_type[i].end(),0) )
        {
            segment_flag = false;
            cerr<<" n_type number != segment number "<<endl;
            
        }
    }
    
    if(segment_flag == false)
    {
        if_correct=false;
        for(config_type::iterator it=config.begin(); it!= config.end(); it++)
        {
            cout<<it->dag<<" ";
        }
        cout<<endl;
        for(config_type::iterator it=config.begin(); it!= config.end(); it++)
        {
            cout<<it->type<<" ";
        }
        cout<<endl;
    }
    
    
    return if_correct;
    
}
void local_config::print_config()
{
    cout<<"dag"<<endl;
    for(config_type::iterator it=config.begin(); it!= config.end(); it++)
    {
        cout<<it->dag<<" ";
    }
    cout<<endl;
    cout<<"type"<<endl;
    for(config_type::iterator it=config.begin(); it!= config.end(); it++)
    {
        cout<<it->type<<" ";
    }
    cout<<endl;
    cout<<"subtype"<<endl;
    for(config_type::iterator it=config.begin(); it!= config.end(); it++)
    {
        cout<<it->subtype<<" ";
    }
    cout<<endl;
    
    cout<<"time"<<endl;
    for(config_type::iterator it=config.begin(); it!= config.end(); it++)
    {
        cout<<it->t<<" ";
    }
    cout<<endl;
    
    
}
