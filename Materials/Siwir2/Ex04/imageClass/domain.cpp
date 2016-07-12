#include<iostream>
#include<vector>
#include <cstdlib>
#include "GrayScaleImage.h"
#include <math.h>
#include<fstream>
#include <algorithm>

using namespace std;
struct values{
    double x;
    double y;
};

class domain{
public:
    void set_initial(double xval, double yval);
    void set_global(double x_length, double y_length, double diameter, int resolution, double viscosity, double acceleration, int end_time, double obstacle_x, double obstacle_y);
    void set_domain ();
    void fill_domain();
    void set_omega();
    void set_cords();
    void find_lattice_units();
    void set_relax_factor();
    void calc_density(int j);
    void calc_velocity(int j);
    void calc_funcq_eq(int j);
    void initialize_2D_vector();
    void initialize_vector();
    void calc_func_q(int i);
    void find_metric_units();
    void display(double valu);
    void find_lattice_velocities();
    void stream_step();
    void initial_setup();
    void swapping_grid();
    void set_vector_zero();
    void find_max_min_vel();
    double Normalizer(int i);

public:
    int nx, ny, cellsize, t_end;
    double x_len, y_len, dia, resol, relax_fac, nu, accel, ux = 0.0, uy = 0.0, obst_x, obst_y;
    vector < vector < vector <double> > > vec;
    vector < vector <double> > func_q;
    vector < vector <double> > func_q_eq;
    vector < vector <double> > dupli_func_q;
    vector <double> omega;
    double cords[9][2];
    //variables in lattice units
    double del_x, del_t, nu_l, ax_l, ay_l=0.0, max_vel, min_vel, norm;
    vector <double> ux_l, uy_l, density;
    vector <double> norm_ux_l;

};

void domain::swapping_grid(){

    func_q = dupli_func_q;
}
void domain::find_max_min_vel(){



    min_vel = *min_element(ux_l.begin(), ux_l.end());
//    display(min_vel);

//    std::cout << "The smallest element is " << min_vel << '\n';
//    std::cout << "The largest element is "  << max_vel << '\n';

    for (int i=0; i<cellsize; i++){
        norm_ux_l[i] = ux_l[i] - min_vel;
    }

    max_vel = *max_element(norm_ux_l.begin(), norm_ux_l.end());
//    display(max_vel);

}

double domain::Normalizer(int i){


    norm = norm_ux_l[i] / max_vel;

    return norm;

}

void domain::initial_setup(){

    int counter = 1;
    for (int i=0; i<cellsize; i++){
        density[i] = 1.0;
        for (int j=0; j<9; j++){
            func_q[i][j]= omega[j];/*func_q[i][j] - relax_fac * (func_q[i][j]-omega[j]) +
                                       (3.0 * omega[j] * density[i] * ((cords[j][0] * ax_l) + (cords[j][1]*ay_l)));*/

//            display(func_q[i][j]);
        ++counter;
        }
    }

}

void domain::stream_step(){

//    nx = 5;
//    ny = 5;
    int count = 1;
//    double x_cell = (obst_x - (dia/2)) * resol / dia;
//    double y_cell = (obst_y - (dia/2)) * resol / dia;
//    double start; // = y_cell * nx + x_cell;

////    display(x_cell);
////    display(y_cell);
//    for(int i=y_cell; i<=(y_cell+resol); i++){
//        for(int j=x_cell; j<=(x_cell+resol); j++){
//            start = i * nx + j;
////            display(start);
//            count++;
//        }
//    }

////    display(count);

//    int theta = 0;
    double x_cord, y_cord;
//    for (theta=0; theta < 360; theta++){
//        x_cord = (dia/2) * cos(theta * M_PI/180);
//        y_cord = (dia/2) * sin(theta * M_PI/180);
////        display(x_cord);
////        display(y_cord);
//    }



    for (int i=0; i<ny; i++){
        for (int k=0; k<nx; k++){
//            display(del_x);
            y_cord = (i+1)*(del_x);
            x_cord = (k+1)*(del_x);
//            display(sqrt(pow((obst_x-x_cord),2) + pow((obst_y-y_cord),2)));
//            display(x_cord);
//            display(0.5*dia);

            if( sqrt(pow((obst_x-x_cord),2) + pow((obst_y-y_cord),2)) <= (0.5*dia/100)){
                dupli_func_q[i*nx+k][0] = 1e-13;
                dupli_func_q[i*nx+k][1] = 1e-13;
                dupli_func_q[i*nx+k][2] = 1e-13;
                dupli_func_q[i*nx+k][3] = 1e-13;
                dupli_func_q[i*nx+k][4] = 1e-13;
                dupli_func_q[i*nx+k][5] = 1e-13;
                dupli_func_q[i*nx+k][6] = 1e-13;
                dupli_func_q[i*nx+k][7] = 1e-13;
                dupli_func_q[i*nx+k][8] = 1e-13;
//                display (i*nx+k);
//                count++;

            }else if((i*nx+k) < nx-1 && (i*nx+k) > 0){ //bottom No-slip boundary

//                display((i*nx+k)+nx-1);
//                display((i*nx+k)+nx);
//                display((i*nx+k)+(nx+1));
//                display(i*nx+k-1);
//                display(i*nx+k);
//                display(i*nx+k+1);
//                display(i*nx+k);
//                display(i*nx+k);
//                display(i*nx+k);

                dupli_func_q [(i*nx+k)+nx-1][6]       = func_q [i*nx+k][6];
                dupli_func_q [(i*nx+k)+nx][2]       = func_q [i*nx+k][2];
                dupli_func_q [(i*nx+k)+(nx+1)][5]   = func_q [i*nx+k][5];
                dupli_func_q [i*nx+k-1][3]          = func_q [i*nx+k][3];
                dupli_func_q [i*nx+k][0]            = func_q [i*nx+k][0];
                dupli_func_q [i*nx+k+1][1]          = func_q [i*nx+k][1];
                dupli_func_q [i*nx+k][5]            = func_q [i*nx+k][7];
                dupli_func_q [i*nx+k][2]            = func_q [i*nx+k][4];
                dupli_func_q [i*nx+k][6]            = func_q [i*nx+k][8];
            }
            else if ((i*nx+k) > (nx*(ny-1)) && (i*nx+k) < (nx*ny)-1 ){ //top No-slip boundary

//                display((i*nx+k)-nx-1);
//                display((i*nx+k)-nx);
//                display((i*nx+k)-nx+1);
//                display(i*nx+k-1);
//                display(i*nx+k);
//                display(i*nx+k+1);
//                display(i*nx+k);
//                display(i*nx+k);
//                display(i*nx+k);

                dupli_func_q [(i*nx+k)-(nx+1)][7]   = func_q [i*nx+k][7];
                dupli_func_q [(i*nx+k)-nx][4]       = func_q [i*nx+k][4];
                dupli_func_q [(i*nx+k)-nx+1][8]     = func_q [i*nx+k][8];
                dupli_func_q [i*nx+k-1][3]          = func_q [i*nx+k][3];
                dupli_func_q [i*nx+k][0]            = func_q [i*nx+k][0];
                dupli_func_q [i*nx+k+1][1]          = func_q [i*nx+k][1];
                dupli_func_q [i*nx+k][7]            = func_q [i*nx+k][5];
                dupli_func_q [i*nx+k][4]            = func_q [i*nx+k][2];
                dupli_func_q [i*nx+k][8]            = func_q [i*nx+k][6];
            }
            else if ((i*nx+k) == (i*nx) && (i*nx+k) > 0 && (i*nx+k) < nx*(ny-1)){//left peroidic boundary

//                display(i*nx+k);
//                display((i*nx+k)+(2*nx-1));
//                display(i*nx+k+1);
//                display((i+1)*nx+k);
//                display((i+1)*nx+(k+1));
//                display((i*nx+k)+(nx-1));
//                display((i*nx+k)-1);
//                display((i*nx+k)-nx);
//                display((i*nx+k)-nx+1);

                dupli_func_q [i*nx+k][0]            = func_q [i*nx+k][0];
                dupli_func_q [(i*nx+k)+(2*nx-1)][6] = func_q [i*nx+k][6];
                dupli_func_q [i*nx+k+1][1]          = func_q [i*nx+k][1];
                dupli_func_q [(i+1)*nx+k][2]        = func_q [i*nx+k][2];
                dupli_func_q [(i+1)*nx+(k+1)][5]    = func_q [i*nx+k][5];
                dupli_func_q [(i*nx+k)+(nx-1)][3]   = func_q [i*nx+k][3];
                dupli_func_q [(i*nx+k)-1][7]        = func_q [i*nx+k][7];
                dupli_func_q [(i*nx+k)-nx][4]       = func_q [i*nx+k][4];
                dupli_func_q [(i*nx+k)-nx+1][8]     = func_q [i*nx+k][8];

            }
            else if ((i*nx+k) == (i*nx)+(nx-1) && (i*nx+k) > nx-1 && (i*nx+k) < (nx*ny)-1){//right peroidic boundary

//                display(i*nx+k);
//                display((i*nx+k)+(nx-1));
//                display(i*nx+k-(nx-1));
//                display((i+1)*nx+k);
//                display(i*nx+(k+1));
//                display((i*nx+k)-1);
//                display((i-1)*nx+(k-1));
//                display((i*nx+k)-nx);
//                display((i*nx+k)-(2*nx-1));

                dupli_func_q [i*nx+k][0]            = func_q [i*nx+k][0];
                dupli_func_q [(i*nx+k)+(nx-1)][6]   = func_q [i*nx+k][6];
                dupli_func_q [i*nx+k-(nx-1)][1]     = func_q [i*nx+k][1];
                dupli_func_q [i*nx+k+nx][2]         = func_q [i*nx+k][2];
                dupli_func_q [i*nx+(k+1)][5]        = func_q [i*nx+k][5];
                dupli_func_q [(i*nx+k)-1][3]        = func_q [i*nx+k][3];
                dupli_func_q [(i-1)*nx+(k-1)][7]    = func_q [i*nx+k][7];
                dupli_func_q [(i*nx+k)-nx][4]       = func_q [i*nx+k][4];
                dupli_func_q [(i*nx+k)-(2*nx-1)][8] = func_q [i*nx+k][8];
            }
            else if ((i*nx+k) == 0){ //SW - boundary

//                display(i*nx+k);
//                display((i*nx+k)+(2*nx-1));
//                display(i*nx+k+1);
//                display((i+1)*nx+k);
//                display((i+1)*nx+(k+1));
//                display((i*nx+k)+(nx-1));
//                display(i*nx+k);
//                display(i*nx+k);
//                display(i*nx+k);

                dupli_func_q [i*nx+k][0]            = func_q [i*nx+k][0];
                dupli_func_q [(i*nx+k)+(2*nx-1)][6] = func_q [i*nx+k][6];
                dupli_func_q [i*nx+k+1][1]          = func_q [i*nx+k][1];
                dupli_func_q [(i+1)*nx+k][2]        = func_q [i*nx+k][2];
                dupli_func_q [(i+1)*nx+(k+1)][5]    = func_q [i*nx+k][5];
                dupli_func_q [(i*nx+k)+(nx-1)][3]   = func_q [i*nx+k][3];
                dupli_func_q [i*nx+k][5]            = func_q [i*nx+k][7];
                dupli_func_q [i*nx+k][2]            = func_q [i*nx+k][4];
                dupli_func_q [i*nx+k][6]            = func_q [i*nx+k][8];
            }
            else if ((i*nx+k) == (nx*(ny-1))){ //NW - boundary

//                display(i*nx+k);
//                display((i*nx+k)-nx+1);
//                display(i*nx+k+1);
//                display((i-1)*nx+k);
//                display(i*nx+k-1);
//                display((i*nx+k)+(nx-1));
//                display(i*nx+k);
//                display(i*nx+k);
//                display(i*nx+k);

                dupli_func_q [i*nx+k][0]            = func_q [i*nx+k][0];
                dupli_func_q [(i*nx+k)-nx+1][8]     = func_q [i*nx+k][8];
                dupli_func_q [i*nx+k+1][1]          = func_q [i*nx+k][1];
                dupli_func_q [(i-1)*nx+k][4]        = func_q [i*nx+k][4];
                dupli_func_q [i*nx+k-1][7]          = func_q [i*nx+k][7];
                dupli_func_q [(i*nx+k)+(nx-1)][3]   = func_q [i*nx+k][3];
                dupli_func_q [(i*nx+k)-1][7]        = func_q [i*nx+k][5];
                dupli_func_q [(i*nx+k)-nx][4]       = func_q [i*nx+k][2];
                dupli_func_q [(i*nx+k)-nx+1][8]     = func_q [i*nx+k][6];
            }
            else if ((i*nx+k) == (nx-1)){ //SE - boundary

//                display(i*nx+k);
//                display((i*nx+k)+(nx-1));
//                display(i*nx+k-(nx-1));
//                display(i*nx+k+nx);
//                display(i*nx+(k+1));
//                display((i*nx+k)-1);
//                display(i*nx+k);
//                display(i*nx+k);
//                display(i*nx+k);

                dupli_func_q [i*nx+k][0]            = func_q [i*nx+k][0];
                dupli_func_q [(i*nx+k)+(nx-1)][6]   = func_q [i*nx+k][6];
                dupli_func_q [i*nx+k-(nx-1)][1]     = func_q [i*nx+k][1];
                dupli_func_q [i*nx+k+nx][2]         = func_q [i*nx+k][2];
                dupli_func_q [i*nx+(k+1)][5]        = func_q [i*nx+k][5];
                dupli_func_q [(i*nx+k)-1][3]        = func_q [i*nx+k][3];
                dupli_func_q [i*nx+k][5]            = func_q [i*nx+k][7];
                dupli_func_q [i*nx+k][2]            = func_q [i*nx+k][4];
                dupli_func_q [i*nx+k][6]            = func_q [i*nx+k][8];
            }
            else if ((i*nx+k) == (nx*ny)-1){ //NE - boundary

//                display(i*nx+k);
//                display((i*nx+k)-(2*nx-1));
//                display(i*nx+k-(nx-1));
//                display(i*nx+k-nx);
//                display(i*nx+k-(nx+1));
//                display((i*nx+k)-1);
//                display(i*nx+k);
//                display((i*nx+k));
//                display((i*nx+k));

                dupli_func_q [i*nx+k][0]            = func_q [i*nx+k][0];
                dupli_func_q [(i*nx+k)-(2*nx-1)][8] = func_q [i*nx+k][8];
                dupli_func_q [i*nx+k-(nx-1)][1]     = func_q [i*nx+k][1];
                dupli_func_q [i*nx+k-nx][4]         = func_q [i*nx+k][4];
                dupli_func_q [i*nx+k-(nx+1)][7]     = func_q [i*nx+k][7];
                dupli_func_q [(i*nx+k)-1][3]        = func_q [i*nx+k][3];
                dupli_func_q [i*nx+k][8]            = func_q [i*nx+k][6];
                dupli_func_q [i*nx+k][4]            = func_q [i*nx+k][2];
                dupli_func_q [i*nx+k][7]            = func_q [i*nx+k][5];
//            } else if(){


                }else
                { //default streaming

//                display(i*nx+k);
//                display((i*nx)+(k+1));
//                display(((i+1)*nx)+k);
//                display((i*nx)+k-1);
//                display((i-1)*nx+k);
//                display((i+1)*nx+k+1);
//                display((i+1)*nx+k-1);
//                display((i-1)*nx+k-1);
//                display((i-1)*nx+k+1);

                dupli_func_q [(i*nx)+k][0]      =   func_q [i*nx+k][0];
                dupli_func_q [(i*nx)+(k+1)][1]  =   func_q [i*nx+k][1];
                dupli_func_q [((i+1)*nx)+k][2]  =   func_q [i*nx+k][2];
                dupli_func_q [(i*nx)+k-1][3]    =   func_q [i*nx+k][3];
                dupli_func_q [(i-1)*nx+k][4]    =   func_q [i*nx+k][4];
                dupli_func_q [(i+1)*nx+k+1][5]  =   func_q [i*nx+k][5];
                dupli_func_q [(i+1)*nx+k-1][6]  =   func_q [i*nx+k][6];
                dupli_func_q [(i-1)*nx+k-1][7]  =   func_q [i*nx+k][7];
                dupli_func_q [(i-1)*nx+k+1][8]  =   func_q [i*nx+k][8];
            }

//            display(dupli_func_q[i*nx+k][0]);
//            display(dupli_func_q[i*nx+k][1]);
//            display(dupli_func_q[i*nx+k][2]);
//            display(dupli_func_q[i*nx+k][3]);
//            display(dupli_func_q[i*nx+k][4]);
//            display(dupli_func_q[i*nx+k][5]);
//            display(dupli_func_q[i*nx+k][6]);
//            display(dupli_func_q[i*nx+k][7]);
//            display(dupli_func_q[i*nx+k][8]);
        }
    }
//    display(count);
}

void domain::initialize_2D_vector(){

    func_q.resize( cellsize , vector<double>( 9 , 0.0 ) );
    func_q_eq.resize( cellsize , vector<double>( 9 , 0.0 ) );
    dupli_func_q.resize( cellsize, vector <double> (9, 0.0));
    omega.resize(9, 0.0);

}

void domain::set_vector_zero(){
    ux_l.resize(cellsize, 0.0);
    uy_l.resize(cellsize, 0.0);
    density.resize(cellsize, 0.0);
    norm_ux_l.resize(cellsize, 0.0);
}

void domain::find_metric_units(){
    x_len = del_x * (nx*100);
//    del_t = 2;

    //find metric viscosity
    nu = nu_l * del_x * del_x / del_t;

    //find metric velocity
//    ux = ux_l * del_x / del_t;
//    uy = uy_l * del_x / del_t;

    //find matric acceleration
    accel = ax_l * del_x / (del_t * del_t);
}

void domain::calc_funcq_eq(int j){

//    display(func_q_eq[j][0]);
//    display(func_q_eq[j][1]);
//    display(func_q_eq[j][2]);
//    display(func_q_eq[j][3]);
//    display(func_q_eq[j][4]);
//    display(func_q_eq[j][5]);
//    display(func_q_eq[j][6]);
//    display(func_q_eq[j][7]);
//    display(func_q_eq[j][8]);

//    display(density[j]);
//    display(ux_l[j]);
//    display(uy_l[j]);

    for (int i=0; i<9; i++){
        func_q_eq[j][i] = omega[i] * density[j] * (1.0 + (3.0*((cords[i][0]*ux_l[j]) + (cords[i][1]*uy_l[j]))) +
                       (4.5*((cords[i][0]*ux_l[j]) + (cords[i][1]*uy_l[j])) * ((cords[i][0]*ux_l[j]) + (cords[i][1]*uy_l[j]))) -
                        (1.5*((ux_l[j]*ux_l[j])+(uy_l[j]*uy_l[j]))));

//    cout << func_q_eq[j][i] << endl;
    }
}

void domain::calc_func_q(int i){

//    display(dupli_func_q[i][0]);
//    display(dupli_func_q[i][1]);
//    display(dupli_func_q[i][2]);
//    display(dupli_func_q[i][3]);
//    display(dupli_func_q[i][4]);
//    display(dupli_func_q[i][5]);
//    display(dupli_func_q[i][6]);
//    display(dupli_func_q[i][7]);
//    display(dupli_func_q[i][8]);

//    display(func_q[i][0]);
//    display(func_q[i][1]);
//    display(func_q[i][2]);
//    display(func_q[i][3]);
//    display(func_q[i][4]);
//    display(func_q[i][5]);
//    display(func_q[i][6]);
//    display(func_q[i][7]);
//    display(func_q[i][8]);

//    ax_l = 0.0;
    for (int j=0; j<9; j++){
        dupli_func_q[i][j] = dupli_func_q[i][j] - relax_fac * (dupli_func_q[i][j]-func_q_eq[i][j]) +
                                       (3.0 * omega[j] * density[i] * ((cords[j][0] * ax_l) + (cords[j][1]*ay_l)));
//        display(dupli_func_q[i][j]);

    }
}

void domain::calc_density(int j){

    density[j] = 0.0;
    for (int i=0; i<9; i++){
        density[j] += dupli_func_q[j][i];
    }
//    display(density[j]);
}

void domain::calc_velocity(int j){

    ux_l[j] = 0.0;
    uy_l[j] = 0.0;

//    display(dupli_func_q[j][0]);
//    display(dupli_func_q[j][1]);
//    display(dupli_func_q[j][2]);
//    display(dupli_func_q[j][3]);
//    display(dupli_func_q[j][4]);
//    display(dupli_func_q[j][5]);
//    display(dupli_func_q[j][6]);
//    display(dupli_func_q[j][7]);
//    display(dupli_func_q[j][8]);

    ux_l[j] += (dupli_func_q[j][1]*cords[1][0] + dupli_func_q[j][5]*cords[5][0] + dupli_func_q[j][8]*cords[8][0] +
                dupli_func_q[j][6]*cords[6][0] + dupli_func_q[j][3]*cords[3][0] + dupli_func_q[j][7]*cords[7][0])/density[j];

    uy_l[j] += (dupli_func_q[j][6]*cords[6][0] + dupli_func_q[j][5]*cords[5][0] + dupli_func_q[j][2]*cords[2][0] +
                dupli_func_q[j][7]*cords[7][0] + dupli_func_q[j][4]*cords[4][0] + dupli_func_q[j][8]*cords[8][0])/density[j];

//        display(ux_l[j]);
//        display(uy_l[j]);

}

void domain::set_cords(){
    cords[0][0] = 0.0;
    cords[0][1] = 0.0;  //center

    cords[1][0] = 1.0;
    cords[1][1] = 0.0;  //East

    cords[2][0] = 0.0;
    cords[2][1] = 1.0;  //North

    cords[3][0] = -1.0;
    cords[3][1] = 0.0;  //west

    cords[4][0] = 0.0;
    cords[4][1] = -1.0; //south

    cords[5][0] = 1.0;
    cords[5][1] = 1.0;  //North-East

    cords[6][0] = -1.0;
    cords[6][1] = 1.0;  //North-West

    cords[7][0] = -1.0;
    cords[7][1] = -1.0; //South-West

    cords[8][0] = 1.0;
    cords[8][1] = -1.0; //South-East

}

void domain::display(double valu){
    cout << "Value is : " << valu << endl;
}

void domain::find_lattice_units(){

    del_x = x_len / (nx*100);
//    display(del_x);
    del_t = 1e-4;

    //find lattice viscosity
    nu_l = nu * del_t / (del_x*del_x);
//    display(nu_l);

    //find lattice velocity
//    ux_l = ux * del_t / del_x;
//    uy_l = uy * del_t / del_x;

    //find lattice acceleration
    ax_l = accel * del_t * del_t / del_x;
//    display(ax_l);

}

void domain::find_lattice_velocities(){

    //find lattice velocity
//    ux_l = ux * del_t / del_x;
//    uy_l = uy * del_t / del_x;
}

void domain::set_relax_factor(){
    relax_fac = 1/((3*nu_l)+0.5);
//    display(relax_fac);
}

void domain::set_omega(){

    omega[0] = 4.0/9.0;
    omega[1] = 1.0/9.0;
    omega[2] = 1.0/9.0; //NW
    omega[3] = 1.0/9.0; //NE
    omega[4] = 1.0/9.0; //W
    omega[5] = 1.0/36.0; //S
    omega[6] = 1.0/36.0; //E
    omega[7] = 1.0/36.0; //N
    omega[8] = 1.0/36.0; //C
}



void domain::set_global(double x_length, double y_length, double diameter, int resolution, double viscosity, double acceleration, int end_time, double obstacle_x, double obstacle_y){
    x_len = x_length;
    y_len = y_length;
    dia = diameter;
    resol = resolution;
    nu = viscosity;
    accel = acceleration;
    t_end = end_time;
    obst_x = obstacle_x;
    obst_y = obstacle_y;
}

void domain::set_initial(double xval, double yval){
    values v;
    v.x = xval;
    v.y = yval;
}

void domain::initialize_vector(){

    int counter = 1;

    for (int i=0; i < cellsize; i++){
        vec.push_back(vector<vector<double> >());

        for (int j=0; j<9; j++){
            vec[i].push_back(vector<double>());

            for (int k=0; k<2; k++){
                vec[i][j].push_back(0.0);
//                cout << vec[i][j][k] << endl;
                counter++;
            }
        }
    }
//    cout << "counter : " << vec.size() << endl;std::cout << "The smallest element is " << min_vel << '\n';
    std::cout << "The largest element is "  << max_vel << '\n';
}

//void domain::initial_vector(){
//        typedef vector<double> v1d;
//        typedef vector<v1d> v2d;
//        typedef vector<v2d> v3d;
//        v3d cell_pos(cellsize, v2d(9, v1d(9, 0.0)));
//}

void domain::set_domain (){

    double temp;

    temp = resol/dia;
    nx = x_len * temp;
    ny = y_len * temp;

    cellsize = nx*ny;

//    cout << "Cell size : " << cellsize << endl;
}

void domain::fill_domain(){

    double temp1, temp2;
    temp1 = x_len / nx;
    temp2 = y_len / ny;

    for (int i=0; i<cellsize; i++){


    }
}

int main(){

    double obst_x = 0.02;
    double obst_y = 0.008;

    double x_len = 6;  //in meters
    double y_len = 2;    //in m
    double dia = 0.5;
    double visco = 1e-6;
    double resol = 30;
    double accel = 0.01;
    int t_end = 3;

    domain d;
    d.set_global(x_len, y_len, dia, resol, visco, accel, t_end, obst_x, obst_y);
    d.set_domain();
    d.initialize_2D_vector();
    d.set_vector_zero();
    d.set_omega();
    d.fill_domain();
    d.set_cords();
    d.find_lattice_units();
    d.set_relax_factor();
    d.initial_setup();
//    cout << "summa" << endl;
    for (int j=0; j<100; j++){
        d.stream_step();


        for (int i=0; i<d.cellsize; i++){
            d.calc_density(i);
            d.calc_velocity(i);
            d.calc_funcq_eq(i);
            d.calc_func_q(i);

        }
        d.swapping_grid();

}
    d.find_max_min_vel();
    std::ofstream solution("solution.txt", std::ofstream::out);
    // for loop to write the output
    for (int i=0; i < d.ny; ++i){
        for (int j=0; j<d.nx; ++j ){
//            if ((i*d.nx+j) == (i*d.nx))
                solution << d.ux_l[i*d.nx+j] << std::endl;

        }
    }
    //    solution << std::endl;
    solution.close();
    double i = 0;

    GrayScaleImage g(360,120);

    for (int i=0; i<d.cellsize; i++){
        d.display(d.Normalizer(i));
        g.setElement(d.nx-1, d.ny-1, d.Normalizer(i));
  }
    g.save("scenario.png");
//    g.GrayScaleImage(3, 4);


//    d.calc_density(0);
//    d.calc_velocity(0);
//

//    for (int i=(d.nx+1); i<d.cellsize-(d.nx+1); i++){
//        d.stream_step(i);
//        d.calc_funcq_eq(i);
//        d.calc_func_q(i);
////    int i=1;
//        d.calc_density(i);
//        d.calc_velocity(i);
////        d.calc_funcq_eq(i);
////        d.calc_func_q(i);
////        d.find_metric_units();
//    }
//    d.initialize_vector();
}
