//
//  mat_vec.hpp
//  cthyb
//
//  Created by 胡昊昱 on 2019/7/7.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#ifndef mat_vec_hpp
#define mat_vec_hpp


#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>




typedef double fermi_type;
typedef boost::numeric::ublas::matrix<fermi_type> fermi_mat;
typedef boost::numeric::ublas::vector<fermi_type> fermi_vec;


#endif /* mat_vec_hpp */
