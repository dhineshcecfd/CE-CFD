#include"LBM.h"
 
const uint_t num_dir = 9;
 
real_t t_alpha[9] =
{
        4.0/9.0,
        1.0/9.0,
        1.0/9.0,
        1.0/9.0,
        1.0/9.0,
        1.0/36.0,
        1.0/36.0,
        1.0/36.0,
        1.0/36.0,
};
real_t c_alpha[9][2] =
{
        { 0.0,0.0 }, // row 0
        { 1.0,0.0 }, // row 1
        { 0.0,1.0 }, // row 2
        { -1.0,0.0 }, // row 3
        { 0.0,-1.0 }, // row 4
        { 1.0,1.0 }, // row 5
        { -1.0,1.0 }, // row 6
        { -1.0,-1.0 }, // row 7
        { 1.0,-1.0 } // row 8
};
 
void Lbm::initialize(PDF_Field &f,DensityField &d,VelocityField &v){
    uint_t x,y;
    x = f.get_xsize();
    y = f.get_ysize();
    for(uint_t i=1; i < x-1; ++i){
        for(uint_t j=1; j < y-1; ++j){
            d(i,j,0) =1.0;
            if(j==y-2 &&((i>=1)&&(i<=x-2))){
                v(i,j,0) = 0.08;
                v(i,j,1) = 0.0;
            }
            else{
                v(i,j,0) = 0.0;
                v(i,j,1) = 0.0;
            }
            for(uint_t k=0; k < num_dir ; ++k){
                f(i,j,k) = t_alpha[k] ;
            }
        }
    }
}
 
void Lbm::calc_density(PDF_Field &f,DensityField &d,uint_t i, uint_t j){
    real_t tmp1 = 0.0;
    for(uint_t k=0; k < num_dir ; ++k){
        tmp1+=f(i,j,k);
    }
    d(i,j,0)=tmp1;
}
 
void Lbm::calc_velocity(PDF_Field &f,VelocityField &v,DensityField &d,uint_t i, uint_t j){
    real_t tmp1 = 0.0;
    real_t tmp2 = 0.0;
    for(uint_t k=0; k < num_dir ; ++k){
        tmp1+=f(i,j,k)*c_alpha[k][0];
        tmp2+=f(i,j,k)*c_alpha[k][1];
    }
    v(i,j,0) = tmp1;
    v(i,j,1) = tmp2;
}
 
void Lbm::collide(PDF_Field &f,real_t omega,VelocityField &v,DensityField &d,uint_t i,uint_t j){
    real_t tmp1,f_eq;
    real_t tmp2 = v(i,j,0)*v(i,j,0)+ v(i,j,1)*v(i,j,1);
    for(uint_t k=0; k < num_dir ; ++k){
        tmp1 = c_alpha[k][0]*v(i,j,0)+c_alpha[k][1]*v(i,j,1);
        f_eq = t_alpha[k]*(d(i,j,0) + 3.0*tmp1 + 4.5*tmp1*tmp1- 1.5*tmp2);
        f(i,j,k) = (1.0 - omega)*f(i,j,k) + omega * f_eq;
    }
}
