//
//  BilinearInterp.cpp
//  BilinearExtrapolation
//
//  Created by Noah Ford on 9/29/18.
//  Copyright © 2018 Noah Ford. All rights reserved.
//

#include "BilinearInterp.hpp"

//Function to perform complex square root
static std::complex<double> complex_sqrt(const std::complex<double> & z)
{
    return pow(z, 1. / 2.);
}
//function to perform complex cubed root
static std::complex<double> complex_cbrt(const std::complex<double> & z)
{
    return pow(z, 1. / 3.);
}
//Function to solve a quartic equation
void solve_quartic(const std::complex<double> coefficients[5], std::complex<double> roots[4])
{
    // The algorithm below was derived by solving the quartic in Mathematica, and simplifying the resulting expression by hand.
    
    const std::complex<double> a = coefficients[4];
    const std::complex<double> b = coefficients[3] / a;
    const std::complex<double> c = coefficients[2] / a;
    const std::complex<double> d = coefficients[1] / a;
    const std::complex<double> e = coefficients[0] / a;
    
    const std::complex<double> Q1 = c * c - 3. * b * d + 12. * e;
    const std::complex<double> Q2 = 2. * c * c * c - 9. * b * c * d + 27. * d * d + 27. * b * b * e - 72. * c * e;
    const std::complex<double> Q3 = 8. * b * c - 16. * d - 2. * b * b * b;
    const std::complex<double> Q4 = 3. * b * b - 8. * c;
    
    const std::complex<double> Q5 = complex_cbrt(Q2 / 2. + complex_sqrt(Q2 * Q2 / 4. - Q1 * Q1 * Q1));
    const std::complex<double> Q6 = (Q1 / Q5 + Q5) / 3.;
    const std::complex<double> Q7 = 2. * complex_sqrt(Q4 / 12. + Q6);
    
    roots[0] = (-b - Q7 - complex_sqrt(4. * Q4 / 6. - 4. * Q6 - Q3 / Q7)) / 4.;
    roots[1] = (-b - Q7 + complex_sqrt(4. * Q4 / 6. - 4. * Q6 - Q3 / Q7)) / 4.;
    roots[2] = (-b + Q7 - complex_sqrt(4. * Q4 / 6. - 4. * Q6 + Q3 / Q7)) / 4.;
    roots[3] = (-b + Q7 + complex_sqrt(4. * Q4 / 6. - 4. * Q6 + Q3 / Q7)) / 4.;
}



void bilinearInterp(int N, int M, double dx, double dy, double* phi, double* F, int* havevalue){
    
    //Calculate Gradients which will be the speed functions
    for(int i = 0; i<N+1;i++)
        for(int j=0; j<M+1; j++){
            //If surface value is defined at point, take the gradient
            if(havevalue[j*(N+1)+i]){
                double rightderiv, leftderiv, upderiv, downderiv;
                int rightleft = 0; //used to count how many sides we take derivative info from
                int updown = 0; //used to count how many sides we take derivative info from
                //Take left and right derivative
                if(i<N && havevalue[j*(N+1)+i+1]){
                    rightderiv = (phi[j*(N+1)+i+1] - phi[j*(N+1)+i])/dx;
                    rightleft += 1;}
                else
                    rightderiv = 0.;
                if(i>0 && havevalue[j*(N+1)+i-1]){
                    leftderiv = (phi[j*(N+1)+i] - phi[j*(N+1)+i-1])/dx;
                    rightleft += 1;}
                else
                    leftderiv = 0.;
                //Take up and down derivative
                if(j<M && havevalue[(j+1)*(N+1)+i]){
                    upderiv = (phi[(j+1)*(N+1)+i] - phi[j*(N+1)+i])/dy;
                    updown += 1;}
                else
                    upderiv = 0.;
                if(j>0 && havevalue[(j-1)*(N+1)+i]){
                    downderiv = (phi[j*(N+1)+i] - phi[(j-1)*(N+1)+i])/dy;
                    updown += 1;}
                else
                    downderiv = 0.;
                
                //Calculate full left/right and up/down derivative
                double rightleftderiv = 0.;
                double updownderiv = 0.;
                
                //Calculate full gradient
                if(rightleft>0 && updown>0){
                    rightleftderiv = (rightderiv+leftderiv)/rightleft;
                    updownderiv = (upderiv+downderiv)/updown;
                    double gradient = sqrt(pow(rightleftderiv,2)+pow(updownderiv,2));
                    //F is inverse of the derivative
                    if(gradient>0.)
                        F[j*(N+1)+i] = 1./gradient;
                    else
                        F[j*(N+1)+i] = 100.;
                }else //If we don't have derivative info from one side, we will fill it in with algorithm
                    F[j*(N+1)+i] = -1.;
            }else //If we don't have value at the point, we will also fill it in with algorithm
                F[j*(N+1)+i] = -1.;
        }
    
    //While we still have values to calculate, indicated by F at point < 0, we continue to calculate more grid points
    //Calculation done for grid points that are adjacent to points that are already accepted
    //Calculation continued in an outward fashion starting at interface
    
    while( *std::min_element(F,F + (N+1)*(M+1)) < 0.){

    std::priority_queue<std::tuple<double,double,int>> q3;
    
    //Find points on boundary and calculate their potential values
    for(int i=0; i<N+1; i++){
        for(int j=0;j<M+1;j++){
            double sidetoside = false; //This tracks whether we have done a diagonal calculation to find the grid point, which takes precidence of on-sided slope calculations
            
            //Check grid to the left for nonzero Biofilm
            if(F[j*(N+1)+i] < 0. && i>0){
                if (F[j*(N+1)+i-1]>0.){
                    double Ftentative = dx/(phi[(j)*(N+1)+i-2]-phi[(j)*(N+1)+i-1]);//F[j*(N+1)+i-1];
                    //double Ftentative = F[j*(N+1)+i-1];
                    double tentativeheight = phi[j*(N+1)+i-1]-dx/Ftentative;
                    Ftentative = F[j*(N+1)+i-1];
                    //q3.emplace(tentativeheight,Ftentative,j*(N+1)+i);
                    
                    //Check if grid below has biofilm too
                    if(j>1)
                        if (F[(j-1)*(N+1)+i]>0.){
                            double phii = phi[j*(N+1)+i-1];
                            double phij = phi[(j-1)*(N+1)+i];
                            double Fi = F[(j)*(N+1)+i-1];
                            double Fj = F[(j-1)*(N+1)+i];
                            
                            std::complex<double> r[4] = {0.,0.,0.,0.};
                            if(abs(phii-phij)<1.e-10){
                                double a = dx*dx+dy*dy;//pow(dx*dx+dy*dy,2);
                                double b = (dx*dx*Fj + dy*dy*Fi);//- 2*(dx*dx+dy*dy)*(dx*dx*Fj + dy*dy*Fi);
                                //double c = pow(dx*dx*Fj + dy*dy*Fi,2);
                                r[0] = (b)/(a);
                                r[1] = (b)/(a);
                                r[2] = -1.;
                                r[3] = -1.;
                                
                            }else{
                                
                                const std::complex<double> a = pow(phii-phij,2)*(dx*dx+dy*dy);
                                const std::complex<double> b = -2*pow(phii-phij,2)*(dx*dx*Fj + dy*dy*Fi);
                                const std::complex<double> c = pow(phii-phij,2)*(dx*dx*Fj*Fj + dy*dy*Fi*Fi) - pow(dx*dx+dy*dy,2);
                                const std::complex<double> d = 2*(dx*dx+dy*dy)*(dx*dx*Fj + dy*dy*Fi);
                                const std::complex<double> e = -pow(dx*dx*Fj + dy*dy*Fi,2);
                                
                                const std::complex<double> coefficients[5] = {e, d, c, b, a};
                                
                                solve_quartic(coefficients, r);
                            }
                            
                            double phinew = -1.e10;
                            double Fnew = 1;
                            for(int rootindx = 0; rootindx<4; rootindx++){
                                if(std::imag(r[rootindx]) <1.e-10 && std::real(r[rootindx])>0.){
                                    double Ftry = std::real(r[rootindx]);
                                    //double phitry = ((dx*dx)*(Fj-Ftry)*phij+dy*dy*(Fi-Ftry)*phii)/(dx*dx*(Fj-Ftry)+dy*dy*(Fi-Ftry));
                                    double determinant = pow(Ftry,2.)*(dx*dx+dy*dy-Ftry*Ftry*pow(phii-phij,2));
                                    double phitry = -1e10;
                                    if(determinant>=0.)
                                        phitry = (dy*dy*pow(Ftry,2.)*phii+dx*dx*pow(Ftry,2)*phij - dx*dy*sqrt(determinant))/((dx*dx+dy*dy)*pow(Ftry,2.));
                                    if(phitry>phinew && phitry < std::min(phii,phij)){
                                        phinew = phitry;
                                        Fnew = Ftry;
                                    }
                                }
                            }
                            tentativeheight = phinew;
                            Ftentative = Fnew;
                            //q3.emplace(tentativeheight,Ftentative,j*(N+1)+i);
                            sidetoside = true;
                        }
                    //Check if grid above has biofilm too
                    if(j<M-1)
                        if (F[(j+1)*(N+1)+i]>0.){
                            double phii = phi[j*(N+1)+i-1];
                            double phij = phi[(j+1)*(N+1)+i];
                            double Fi = F[(j)*(N+1)+i-1];
                            double Fj = F[(j+1)*(N+1)+i];
                            
                            std::complex<double> r[4] = {0.,0.,0.,0.};
                            if(abs(phii-phij)<1.e-10){
                                double a = dx*dx+dy*dy;//pow(dx*dx+dy*dy,2);
                                double b = (dx*dx*Fj + dy*dy*Fi);//- 2*(dx*dx+dy*dy)*(dx*dx*Fj + dy*dy*Fi);
                                //double c = pow(dx*dx*Fj + dy*dy*Fi,2);
                                r[0] = (b)/(a);
                                r[1] = (b)/(a);
                                r[2] = -1.;
                                r[3] = -1.;
                                
                            }else{
                                
                                const std::complex<double> a = pow(phii-phij,2)*(dx*dx+dy*dy);
                                const std::complex<double> b = -2*pow(phii-phij,2)*(dx*dx*Fj + dy*dy*Fi);
                                const std::complex<double> c = pow(phii-phij,2)*(dx*dx*Fj*Fj + dy*dy*Fi*Fi) - pow(dx*dx+dy*dy,2);
                                const std::complex<double> d = 2*(dx*dx+dy*dy)*(dx*dx*Fj + dy*dy*Fi);
                                const std::complex<double> e = -pow(dx*dx*Fj + dy*dy*Fi,2);
                                
                                const std::complex<double> coefficients[5] = {e, d, c, b, a};
                                
                                solve_quartic(coefficients, r);
                            }
                            
                            double phinew = -1.e10;
                            double Fnew = 1;
                            for(int rootindx = 0; rootindx<4; rootindx++){
                                if(std::imag(r[rootindx]) <1.e-10 && std::real(r[rootindx])>0.){
                                    double Ftry = std::real(r[rootindx]);
                                    //double phitry = ((dx*dx)*(Fj-Ftry)*phij+dy*dy*(Fi-Ftry)*phii)/(dx*dx*(Fj-Ftry)+dy*dy*(Fi-Ftry));
                                    double determinant = pow(Ftry,2.)*(dx*dx+dy*dy-Ftry*Ftry*pow(phii-phij,2));
                                    double phitry = -1e10;
                                    if(determinant>=0.)
                                        phitry = (dy*dy*pow(Ftry,2.)*phii+dx*dx*pow(Ftry,2)*phij - dx*dy*sqrt(determinant))/((dx*dx+dy*dy)*pow(Ftry,2.));
                                    if(phitry>phinew && phitry < std::min(phii,phij)){
                                        phinew = phitry;
                                        Fnew = Ftry;
                                    }
                                }
                            }
                            tentativeheight = phinew;
                            Ftentative = Fnew;
                            //q3.emplace(tentativeheight,Ftentative,j*(N+1)+i);
                            sidetoside = true;
                        }
                    q3.emplace(tentativeheight,Ftentative,j*(N+1)+i);
                }
            }
            
            //Check grid to the right for nonzero Biofilm
            if(F[j*(N+1)+i] < 0. && i<N){
                if (F[j*(N+1)+i+1]>0.){
                    double Ftentative = dx/(phi[(j)*(N+1)+i+2]-phi[(j)*(N+1)+i+1]);//F[j*(N+1)+i+1];
                    //double Ftentative = F[j*(N+1)+i+1];
                    double tentativeheight = phi[j*(N+1)+i+1]-dx/Ftentative;
                    Ftentative =F[j*(N+1)+i+1];
                    //q3.emplace(tentativeheight,Ftentative,j*(N+1)+i);
                    
                    //Check if grid below has biofilm too
                    if(j>1)
                        if (F[(j-1)*(N+1)+i]>0.){
                            double phii = phi[j*(N+1)+i+1];
                            double phij = phi[(j-1)*(N+1)+i];
                            double Fi = F[j*(N+1)+i+1];
                            double Fj = F[(j-1)*(N+1)+i];
                            
                            std::complex<double> r[4] = {0.,0.,0.,0.};
                            if(abs(phii-phij)<1.e-1){
                                double a = dx*dx+dy*dy;//pow(dx*dx+dy*dy,2);
                                double b = (dx*dx*Fj + dy*dy*Fi);//- 2*(dx*dx+dy*dy)*(dx*dx*Fj + dy*dy*Fi);
                                //double c = pow(dx*dx*Fj + dy*dy*Fi,2);
                                r[0] = (b)/(a);
                                r[1] = (b)/(a);
                                r[2] = -1.;
                                r[3] = -1.;
                                
                            }else{
                                
                                const std::complex<double> a = pow(phii-phij,2)*(dx*dx+dy*dy);
                                const std::complex<double> b = -2*pow(phii-phij,2)*(dx*dx*Fj + dy*dy*Fi);
                                const std::complex<double> c = pow(phii-phij,2)*(dx*dx*Fj*Fj + dy*dy*Fi*Fi) - pow(dx*dx+dy*dy,2);
                                const std::complex<double> d = 2*(dx*dx+dy*dy)*(dx*dx*Fj + dy*dy*Fi);
                                const std::complex<double> e = -pow(dx*dx*Fj + dy*dy*Fi,2);
                                
                                const std::complex<double> coefficients[5] = {e, d, c, b, a};
                                
                                solve_quartic(coefficients, r);
                            }
                            
                            double phinew = -1.e10;
                            double Fnew = 1;
                            for(int rootindx = 0; rootindx<4; rootindx++){
                                if(std::imag(r[rootindx]) <1.e-10 && std::real(r[rootindx])>0.){
                                    double Ftry = std::real(r[rootindx]);
                                    //double phitry = ((dx*dx)*(Fj-Ftry)*phij+dy*dy*(Fi-Ftry)*phii)/(dx*dx*(Fj-Ftry)+dy*dy*(Fi-Ftry));
                                    double determinant = pow(Ftry,2.)*(dx*dx+dy*dy-Ftry*Ftry*pow(phii-phij,2));
                                    double phitry = -1e10;
                                    if(determinant>=0.)
                                        phitry = (dy*dy*pow(Ftry,2.)*phii+dx*dx*pow(Ftry,2)*phij - dx*dy*sqrt(determinant))/((dx*dx+dy*dy)*pow(Ftry,2.));
                                    if(phitry>phinew && phitry < std::min(phii,phij)){
                                        phinew = phitry;
                                        Fnew = Ftry;
                                    }
                                }
                            }
                            
                            tentativeheight = phinew;
                            Ftentative = Fnew;
                            //q3.emplace(tentativeheight,Ftentative,j*(N+1)+i);
                            sidetoside = true;
                        }
                    
                    //Check if grid above has biofilm too
                    if(j<M-1)
                        if (F[(j+1)*(N+1)+i]>0.){
                            double phii = phi[j*(N+1)+i+1];
                            double phij = phi[(j+1)*(N+1)+i];
                            double Fi = F[j*(N+1)+i+1];
                            double Fj = F[(j+1)*(N+1)+i];
                            
                            std::complex<double> r[4] = {0.,0.,0.,0.};
                            if(abs(phii-phij)<1.e-10){
                                double a = dx*dx+dy*dy;//pow(dx*dx+dy*dy,2);
                                double b = (dx*dx*Fj + dy*dy*Fi);//- 2*(dx*dx+dy*dy)*(dx*dx*Fj + dy*dy*Fi);
                                //double c = pow(dx*dx*Fj + dy*dy*Fi,2);
                                r[0] = (b)/(a);
                                r[1] = (b)/(a);
                                r[2] = -1.;
                                r[3] = -1.;
                                
                            }else{
                                
                            const std::complex<double> a = pow(phii-phij,2)*(dx*dx+dy*dy);
                            const std::complex<double> b = -2*pow(phii-phij,2)*(dx*dx*Fj + dy*dy*Fi);
                            const std::complex<double> c = pow(phii-phij,2)*(dx*dx*Fj*Fj + dy*dy*Fi*Fi) - pow(dx*dx+dy*dy,2);
                            const std::complex<double> d = 2*(dx*dx+dy*dy)*(dx*dx*Fj + dy*dy*Fi);
                            const std::complex<double> e = -pow(dx*dx*Fj + dy*dy*Fi,2);
                            
                            const std::complex<double> coefficients[5] = {e, d, c, b, a};
                            
                            solve_quartic(coefficients, r);
                            }
                            
                            double phinew = -1.e10;
                            double Fnew = 1;
                            for(int rootindx = 0; rootindx<4; rootindx++){
                                if(std::imag(r[rootindx]) <1.e-10 && std::real(r[rootindx])>0.){
                                    double Ftry = std::real(r[rootindx]);
                                    //double phitry = ((dx*dx)*(Fj-Ftry)*phij+dy*dy*(Fi-Ftry)*phii)/(dx*dx*(Fj-Ftry)+dy*dy*(Fi-Ftry));
                                    double determinant = pow(Ftry,2.)*((dx*dx+dy*dy)-Ftry*Ftry*pow(phii-phij,2));
                                    double phitry = -1e10;
                                    if(determinant>=0.)
                                        phitry = (dy*dy*pow(Ftry,2.)*phii+dx*dx*pow(Ftry,2)*phij - dx*dy*sqrt(determinant))/((dx*dx+dy*dy)*pow(Ftry,2.));
                                    if(phitry>phinew && phitry < std::min(phii,phij)){
                                        phinew = phitry;
                                        Fnew = Ftry;
                                    }
                                }
                            }
                            tentativeheight = phinew;
                            Ftentative = Fnew;
                            //q3.emplace(tentativeheight,Ftentative,j*(N+1)+i);
                            sidetoside = true;
                        }
                    q3.emplace(tentativeheight,Ftentative,j*(N+1)+i);
                }
            }
            
            //Check grid below for nonzero Biofilm
            if(F[j*(N+1)+i] < 0. && j>0 &&!sidetoside){
                if (F[(j-1)*(N+1)+i]>0.){
                    double Ftentative = dy/(phi[(j-2)*(N+1)+i]-phi[(j-1)*(N+1)+i]);//F[(j-1)*(N+1)+i];
                    //double Ftentative =F[(j-1)*(N+1)+i];
                    double tentativeheight = phi[(j-1)*(N+1)+i]-dy/Ftentative;
                    Ftentative = F[(j-1)*(N+1)+i];
                    q3.emplace(tentativeheight,Ftentative,j*(N+1)+i);
                }
            }
            
            //Check grid above for nonzero Biofilm
            if(F[j*(N+1)+i] < 0. && j<M &&!sidetoside){
                if (F[(j+1)*(N+1)+i]>0.){
                    double Ftentative = dy/(phi[(j+2)*(N+1)+i]-phi[(j+1)*(N+1)+i]);//F[(j+1)*(N+1)+i];
                    //double Ftentative = F[(j+1)*(N+1)+i];
                    double tentativeheight = phi[(j+1)*(N+1)+i]-dy/Ftentative;
                    Ftentative = F[(j+1)*(N+1)+i];
                    q3.emplace(tentativeheight,Ftentative,j*(N+1)+i);
                }
            }
            
            
        }
    }
    
    //Put the calculated values into phi and F grids
    while (!q3.empty()){
        std::tuple<double,double,int> point = q3.top(); //Get next largest biofilm heigh
        q3.pop();
        
        int index = std::get<2>(point);
        if(F[index]<0.){            //If we haven't already put data into grid here, but height and growth in now
            phi[index] = std::get<0>(point);
            F[index] = std::get<1>(point);
        }
    }
    }
    
    
    //delete[] F;
    
    
    
    
}
