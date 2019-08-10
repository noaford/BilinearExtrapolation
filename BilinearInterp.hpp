//
//  BilinearInterp.hpp
//  BilinearExtrapolation
//
//  Created by Noah Ford on 9/29/18.
//  Copyright © 2018 Noah Ford. All rights reserved.
//

#ifndef BilinearInterp_hpp
#define BilinearInterp_hpp

#include <stdio.h>
#include <math.h>
#include <complex>
#include <algorithm>
#include <queue>

void bilinearInterp(int N, int M, double dx, double dy, double* phi, double* F, int* havevalue);

#endif /* BilinearInterp_hpp */
