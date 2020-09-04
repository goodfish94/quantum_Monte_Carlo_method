//
//  random.hpp
//  cthyb
//
//  Created by 胡昊昱 on 2019/7/8.
//  Copyright © 2019 胡昊昱. All rights reserved.
//

#ifndef random_hpp
#define random_hpp

#include <ctime>
#include <cstdlib>

struct rand_gen
{
    rand_gen(int seed ){ srand( seed );}
//    rand_gen(int seed ){ srand( seed );}
    
    double random_double( double x)
    {//[0,x)
        return (x*rand())/(RAND_MAX + 1.0);
        
    }
    int random_int(int N)
    {// [0,N)
        return static_cast<int>(random_double(1.0)*N);
    }
    
};
#endif /* random_hpp */
