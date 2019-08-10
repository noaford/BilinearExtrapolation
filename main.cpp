//
//  main.cpp
//  BilinearExtrapolation
//
//  Created by Noah Ford on 9/29/18.
//  Copyright © 2018 Noah Ford. All rights reserved.
//
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "BilinearInterp.hpp"

int main(int argc, const char * argv[]) {
    //Length of Grid
    double Lx = 1.;
    double Ly = 1.;
    
    //Number of Points on Grid
    int N = 500;
    int M = 500;
    
    //Size of grid spacing
    double dx = Lx/N;
    double dy = Ly/M;
    
    //Max height of hemisphere
    double maxheight = .1;
    
    //x and y grid
    double* x = new double[N+1];
    double* y = new double[M+1];
    
    for(int i=0;i<N+1;i++)
        x[i] = i*dx;
    for(int j=0;j<M+1;j++)
        y[j] = j*dy;
    
    std::cout << "Initialize Surface" << std::endl;
    //Initialize value
    double* Phi = new double[(N+1)*(M+1)]; //Height Function
    double* F = new double[(N+1)*(M+1)]; //Speed Function
    int* havevalue = new int[(N+1)*(M+1)]; //Boolian to indicate weather we already have value for point
    
    for(int i=0;i<N+1;i++)
        for(int j=0;j<M+1;j++){
            double heightatlocation = maxheight-(x[i]-x[N/2])*(x[i]-x[N/2])-(y[j]-y[M/2])*(y[j]-y[M/2]);
            //double heightatlocation = 2.*(x[i]+y[j]-1+x[i]*y[i]);
            if(heightatlocation>0.){
                Phi[j*(N+1)+i] = heightatlocation;
                havevalue[j*(N+1)+i] = 1;
            }else{
                Phi[j*(N+1)+i] = 0.;
                havevalue[j*(N+1)+i] = 0;
            }
        }
    
    std::cout << "Perform Bilinear Interpolation" << std::endl;
    //Interpolate to the Rest of Grid
    bilinearInterp(N,M,dx,dy,Phi,F, havevalue);
    
    std::cout << "Outputting to File" << std::endl;
    //Output to File
    FILE *file = fopen("data.bin", "wb");
    if(file == NULL)
    {
        std::cout << "Error" << std::endl;
    }
    fwrite(&N,sizeof(int),1,file);
    fwrite(x,sizeof(double),(N+1),file);
    fwrite(&M,sizeof(int),1,file);
    fwrite(y,sizeof(double),(M+1),file);
    fwrite(Phi,sizeof(double),(M+1)*(N+1),file);
    fwrite(F,sizeof(double),(M+1)*(N+1),file);
    fclose(file);
    
    
    
    // Delete Variables
    delete[] x;
    delete[] y;
    delete[] Phi;
    delete[] F;
    delete[] havevalue;
    
    
    
    return 0;
}
