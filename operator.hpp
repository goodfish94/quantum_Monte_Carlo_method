//
//  operator.hpp
//  cthyb
//
//  Created by 胡昊昱 on 7/15/19.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#ifndef operator_hpp
#define operator_hpp



struct op
{// operator
    op(double time_in, bool dag_in,int type_in, int subtype_in):t(time_in),type(type_in),subtype(subtype_in),dag(dag_in){}
    op()
    {
        t=0.0;
        dag=true;
        type=0;
        subtype = 0;
    }
    
    double t;
    //time
    int type;
    // op type
    // 1 , comes form c^\dag d or  d^\dag c
    // 2 , comes from L^\dag \phi or L^\dag \phi
    // 3 , attached with J L^\dag L^-
    int subtype;// subtype
    bool dag;// if dag
    
    op& operator=(const op & b)
    {
        t=b.t;
        type = b.type;
        subtype = b.subtype;
        dag = b.dag;
        return *this;
    }
    
    void set_value( double t_, bool dag_, int type_, int subtype_)
    {
        t = t_; dag=dag_; type = type_; subtype = subtype_;
        
    }
    
};
inline bool operator==( const op& a, const op& b){ return ( a.t==b.t && a.type==b.type && a.subtype==b.subtype && a.dag==b.dag) ;}
inline bool operator!=( const op& a, const op& b){ return !(a==b) ;}
struct comp_op
{
    bool operator()( const op& a, const op& b){return a.t<b.t;}
};

#endif /* operator_hpp */
