#ifndef LBM_H_
#define LBM_H_
#include"Grid.h"
using namespace lbm;
 
class Lbm{
public:
    void initialize(PDF_Field &f,DensityField &d,VelocityField &v);
    void calc_density(PDF_Field &f,DensityField &d,uint_t i, uint_t j);
    void calc_velocity(PDF_Field &f,VelocityField &v,DensityField &d,uint_t i, uint_t j);
    void collide(PDF_Field &f,real_t omega,VelocityField &v,DensityField &d,uint_t i,uint_t j);
};
 
#endif /* LBM_H_ */
