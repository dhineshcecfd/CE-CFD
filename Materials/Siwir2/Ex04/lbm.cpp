/*
 *  steps:
 * create class for lattice ----- a 3D array to store the data--- 2D for cell space and 1D for distribution func
 * Implement stream step
 * Boundary conditions
 * implement collide step
 * use GrayScaleImage class to write out density and velocity fields
 * Normalize the image 0v/rho = 0(white), max gray value 255 (black)
 */
#include<iostream>
#include<vector>

using namespace std;

vector <double> omega(9, 0.0);
struct values{
    double x;
    double y;
    double func;
};

void fill_omega(){

    for (int i=0; i<9; i++){
        if (i==0 || i==2 || i==6 || i==8){
            omega[i] = 0.027777778;
        }

        if (i==1 || i==3 || i==5 || i==7){
            omega[i] = 0.111111111;
        }

        if (i==4){
            omega[i] = 0.444444444;
        }
    }
}

namespace lbm{
enum lattice_direction{
    C = 0,
    E = 1,
    N = 2,
    W = 3,
    S = 4,
    NE = 5,
    NW = 6,
    SW = 7,
    SE = 8
};

template < typename Type, size_t Cellsize >

class Grid{
public:

    Grid():xsize_(0), ysize_(0), data_(0){}
    Grid( size_t xsize, size_t ysize ): xsize_(xsize), ysize_(ysize),
            data_(new Type[Cellsize * xsize * ysize]){}
    ~Grid(){
        delete []data_;
        data_ = 0;
    }
    Type& operator()( size_t x, size_t y, size_t f ){
        assert(x < xsize_ && y < ysize_ && f < Cellsize);
        return data_[y*xsize_*Cellsize + x*Cellsize + f];
    }
    Type operator()( size_t x, size_t y, size_t f ) const{
        assert(x < xsize_ && y < ysize_ && f < Cellsize);
        return data_[y*xsize_*Cellsize + x*Cellsize + f];
    }
    void swap(Grid& grid){
        ysize_ = grid.ysize_ ;
        xsize_ = grid.xsize_ ;
        std::swap(data_,grid.data_ );
    }
    size_t get_xsize(){
        return xsize_;
    }
    size_t get_ysize(){
        return ysize_;
    }

private:
    uint_t xsize_;		// Number of nodes in x-dimension
    uint_t ysize_;		// Number of nodes in y-dimension
    Type* data_;		// Linearized, 1-dimensional representation
                        // of the 2D data grid

};

typedef Grid<real_t, 9> DistFunc;
typedef Grid<real_t, 2> VelField;
typedef Grid<real_t, 1> DensField;
}

class LBM{
public:
    void initial_setup(DistFunc &func,DensField &d,VelField &v);
    void calc_density(PDF_Field &f,DensityField &d,uint_t i, uint_t j);
    void calc_velocity(DistFunc &f, VelField &v, DensField &d, int i, int j);
    void collide(DistFunc &func, double relax_fac, VelField &vel, DensField &dens, int i, int j);
};


class lattice{

public:
    void set_value(double x_val, double y_val, double func_val);
    void set_omega();
    void set_relaxation(double nu_l);
    void display_value(double val);
    void initial_setup();
//    void calc_density();
//    void calc_velocity();
    void calc_dist_function();
    void calc_dist_function_equli();


public:
    //3D array for stroing data
    vector < vector < vector <double> > > cell;
    double density;
    double velocity;
    vector <double> omega;
    static double const c[9][2];
    const double c_s = 0.577350269;
};

void lattice::display_value(double val){
    cout << "value of : " << val << endl;
}

void lattice::calc_dist_function(){

}

void LBM::initial_setup(DistFunc &func,DensField &d,VelField &v){
    size_t x,y;
        x = func.get_xsize();
        y = func.get_ysize();

        for(size_t i=1; i < x-1; ++i){
            for(size_t j=1; j < y-1; ++j){

                //make density as 1.0 at initial setup
                d(i,j,0) =1.0;

                //make initial velocity as "zero"
//                if(j==y-2 &&((i>=1)&&(i<=x-2))){
//                    v(i,j,0) = 0.08;
//                    v(i,j,1) = 0.0;
//                }
//                else{
                    v(i,j,0) = 0.0;
                    v(i,j,1) = 0.0;
//                }
                for(uint_t k=0; k < 9 ; ++k){
                    func(i,j,k) = omega[k] ;
                }
            }
        }
}

void LBM::calc_density(DistFunc &func,DensField &dens,int i, int j){

    double temp = 0.0;
    for(int k=0; k < 9 ; ++k){
        temp += func(i,j,k);
    }
    dens(i,j,0)=temp;
}

void LBM::calc_velocity(DistFunc &func,VelField &vel,DensField &dens,int i, int j){

    double temp1 = 0.0;
    double temp2 = 0.0;

    for(int k=0; k < 9 ; ++k){
        temp1 += func(i,j,k) * c[k][0];
        temp2 += func(i,j,k) * c[k][1];
    }

    vel(i,j,0) = temp1;
    vel(i,j,1) = temp2;
}

void LBM::collide(DistFunc &func, double relax_fac, VelField &vel, DensField &dens, int i,int j){

    double temp1, temp2, func_eq;

    temp2 = vel(i,j,0) * vel(i,j,0) +  vel(i,j,1) * vel(i,j,1);

    for(int k=0; k < 9 ; ++k){

        temp1 = c[k][0] * vel(i,j,0) + c[k][1] * vel(i,j,1);
        func_eq = omega[k] * (dens(i,j,0) + (3.0 * temp1) + (4.5 * temp1 * temp1) - (1.5 * temp2));
        func(i,j,k) = (1.0 - relax_fac) * func(i,j,k) + relax_fac * func_eq;

    }
}

void lattice::set_relaxation(double nu_l){

    double relax_fac;
    relax_fac = 1/((3*nu_l)+0.5);

    display_value(relax_fac); //lies between 0 < w < 2

}

void lattice::set_value(double x_val, double y_val, double func_val){
    values v;
    v.x = x_val;
    v.y = y_val;
    v.func = func_val;
    std::cout << "values of : " << v.x << v.y << v.func << endl;
}



double const lattice::c[9][2] = {
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 1.0},
        {-1.0, 0.0},
        {0.0, -1.0},
        {1.0, 1.0},
        {-1.0, 1.0},
        {-1.0, -1.0},
        {1.0, -1.0}
};


int main(int argc, char *argv[]){

    if(argc != 2){
        cout << "Not enough arguments" << endl;
    }

    fill_omega();
    if (argv[1] == "scenario1"){

        double vis = 10e-6;
        double end_t = 3;
        double accl = 0.01;
        double res_cyl_dia = 30;
    }

    if (argv[1] == "scenario2"){
        double vis = 10e-6;
        double end_t = 5;
        double accl = 0.016;
        double res_cyl_dia = 60;
    }
    lattice l;
    values vl;
    l.set_value(0.5, 0.6, 2.33);

}
