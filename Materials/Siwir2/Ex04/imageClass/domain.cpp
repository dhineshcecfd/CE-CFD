#include<iostream>
#include<vector>
#include <cstdlib>
#include "GrayScaleImage.h"
#include <math.h>
#include<fstream>
#include <algorithm>
#include <string>

using namespace std;

class domain{
public:
    void set_global(double x_length, double y_length, double diameter, int resolution, double viscosity, double acceleration, int end_time, double obstacle_x, double obstacle_y);
    void set_domain ();
    void set_omega();
    void set_cords();
    void find_lattice_units();
    void set_relax_factor();
    void calc_density(int j);
    void calc_velocity(int j);
    void calc_funcq_eq(int j);
    void initialize_2D_vector();
    void calc_func_q(int i);
    void find_metric_units();
    void display(double valu);
    void find_lattice_velocities();
    void stream_step();
    void initial_setup();
    void swapping_grid();
    void set_vector_zero();
    void find_max_min_vel();
    void Normalizer();

public:
    int nx, ny, cellsize, t_end;
    double x_len, y_len, dia, resol, relax_fac, nu, accel, ux = 0.0, uy = 0.0, obst_x, obst_y;
    vector < vector <double> > func_q;
    vector < vector <double> > func_q_eq;
    vector < vector <double> > dupli_func_q;
    vector <double> omega;
    double cords[9][2];
    double del_x, del_t, nu_l, ax_l, ay_l=0.0, max_vel, min_vel, norm;
    vector <double> ux_l, uy_l, density;
};

void domain::swapping_grid(){

    func_q = dupli_func_q;
}
void domain::find_max_min_vel(){

    for (int j=0; j<cellsize; j++){
        if (ux_l[j] < 0.0){
            ux_l[j] = 0.0;
        }
    }
    min_vel = *min_element(ux_l.begin(), ux_l.end());
    max_vel = *max_element(ux_l.begin(), ux_l.end());
    display(max_vel);
}

void domain::Normalizer(){

    for (int i=0;i<cellsize; i++){
        ux_l[i] /= max_vel;
    }
}

void domain::initial_setup(){

    for (int i=0; i<cellsize; i++){
        density[i] = 1.0;
        ux_l[i] = 0.0;
        uy_l[i] = 0.0;

        for (int j=0; j<9; j++){
            func_q[i][j]= omega[j];
        }
    }
}

void domain::stream_step(){

    double x_cord, y_cord;

    for (int i=0; i<ny; i++){
        for (int k=0; k<nx; k++){
            y_cord = (i+1)*(del_x);
            x_cord = (k+1)*(del_x);

            if( sqrt(pow((obst_x-x_cord),2) + pow((obst_y-y_cord),2)) <= (0.5*dia/100)){

                }else if((i*nx+k) < nx-1 && (i*nx+k) > 0){ //bottom No-slip boundary

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

                dupli_func_q [i*nx+k][0]            = func_q [i*nx+k][0];
                dupli_func_q [(i*nx+k)-(2*nx-1)][8] = func_q [i*nx+k][8];
                dupli_func_q [i*nx+k-(nx-1)][1]     = func_q [i*nx+k][1];
                dupli_func_q [i*nx+k-nx][4]         = func_q [i*nx+k][4];
                dupli_func_q [i*nx+k-(nx+1)][7]     = func_q [i*nx+k][7];
                dupli_func_q [(i*nx+k)-1][3]        = func_q [i*nx+k][3];
                dupli_func_q [i*nx+k][8]            = func_q [i*nx+k][6];
                dupli_func_q [i*nx+k][4]            = func_q [i*nx+k][2];
                dupli_func_q [i*nx+k][7]            = func_q [i*nx+k][5];

                }else
                { //default streaming

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
        }
    }
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
}

void domain::calc_funcq_eq(int j){

    for (int i=0; i<9; i++){
        func_q_eq[j][i] = omega[i] * density[j] * (1.0 + (3.0*((cords[i][0]*ux_l[j]) + (cords[i][1]*uy_l[j]))) +
                       (4.5*((cords[i][0]*ux_l[j]) + (cords[i][1]*uy_l[j])) * ((cords[i][0]*ux_l[j]) + (cords[i][1]*uy_l[j]))) -
                        (1.5*((ux_l[j]*ux_l[j])+(uy_l[j]*uy_l[j]))));
    }
}

void domain::calc_func_q(int i){

    for (int j=0; j<9; j++){
        dupli_func_q[i][j] = dupli_func_q[i][j] - relax_fac * (dupli_func_q[i][j]-func_q_eq[i][j]) +
                                       (3.0 * omega[j] * density[i] * ((cords[j][0] * ax_l) + (cords[j][1]*ay_l)));
    }
}

void domain::calc_density(int j){

    density[j] = 0.0;
    for (int i=0; i<9; i++){
        density[j] += dupli_func_q[j][i];
    }
}

void domain::calc_velocity(int j){

    ux_l[j] = 0.0;
    uy_l[j] = 0.0;

    if (density[j] == 0.0){
        ux_l[j] += 0;
        uy_l[j] += 0;
    }else{

    ux_l[j] += (dupli_func_q[j][1]*cords[1][0] + dupli_func_q[j][5]*cords[5][0] + dupli_func_q[j][8]*cords[8][0] +
                dupli_func_q[j][6]*cords[6][0] + dupli_func_q[j][3]*cords[3][0] + dupli_func_q[j][7]*cords[7][0])/density[j];

    uy_l[j] += (dupli_func_q[j][6]*cords[6][0] + dupli_func_q[j][5]*cords[5][0] + dupli_func_q[j][2]*cords[2][0] +
                dupli_func_q[j][7]*cords[7][0] + dupli_func_q[j][4]*cords[4][0] + dupli_func_q[j][8]*cords[8][0])/density[j];
    }

    display(ux_l[j]);

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
    del_t = 1e-4;

    //find lattice viscosity
    nu_l = nu * del_t / (del_x*del_x);

    //find lattice velocity
//    ux_l = ux * del_t / del_x;
//    uy_l = uy * del_t / del_x;

    //find lattice acceleration
    ax_l = accel * del_t * del_t / del_x;
}

void domain::set_relax_factor(){

    relax_fac = 1/((3*nu_l)+0.5);
}

void domain::set_omega(){

    omega[0] = 4.0/9.0;
    omega[1] = 1.0/9.0;
    omega[2] = 1.0/9.0;
    omega[3] = 1.0/9.0;
    omega[4] = 1.0/9.0;
    omega[5] = 1.0/36.0;
    omega[6] = 1.0/36.0;
    omega[7] = 1.0/36.0;
    omega[8] = 1.0/36.0;
}

void domain::set_global(double x_length, double y_length, double diameter, int resolution, double viscosity, double acceleration, int end_time, double obstacle_x, double obstacle_y){

    x_len   = x_length;
    y_len   = y_length;
    dia     = diameter;
    resol   = resolution;
    nu      = viscosity;
    accel   = acceleration;
    t_end   = end_time;
    obst_x  = obstacle_x;
    obst_y  = obstacle_y;
}

void domain::set_domain (){

    double temp;
    temp = resol/dia;
    nx = x_len * temp;
    ny = y_len * temp;

    cellsize = nx*ny;
}

int main(int argc, char* argv[]){

    if (argc <= 1){
        cout << "Too few arguments" << endl;
    }

    double obst_x = 0.02;
    double obst_y = 0.008;

    double x_len = 6;
    double y_len = 2;
    double dia = 0.5;
    double visco = 1e-6;
    double resol, accel, t_end;

    domain d;
    if (string(argv[1]) == "scenario1"){
        resol = 30;
        accel = 0.01;
        t_end = 0.0003;

        d.set_global(x_len, y_len, dia, resol, visco, accel, t_end, obst_x, obst_y);
        d.set_domain();
        d.initialize_2D_vector();
        d.set_vector_zero();
        d.set_omega();
        d.set_cords();
        d.find_lattice_units();
        d.set_relax_factor();

        d.initial_setup();
        cout << "Time step   -------------------> " << d.del_t << endl;
        cout << "Lattice spacing ---------------> " << d.del_x << endl;
        cout << "Acceleration in lattice units -> " << d.ax_l << endl;
        cout << "Relaxation rate ---------------> " << d.relax_fac << endl;

        for (double j=0.0; j<2; j+=1){
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
        d.Normalizer();

        std::ofstream solution("solution.txt", std::ofstream::out);
        for (int i=0; i < d.ny; ++i){
            for (int j=0; j<d.nx; ++j ){
                solution << d.ux_l[i*d.nx+j] << std::endl;
            }
        }
        solution.close();

        double max = *max_element(d.ux_l.begin(), d.ux_l.end());
        double min = *min_element(d.ux_l.begin(), d.ux_l.end());
        d.display(max);
        d.display(min);

        GrayScaleImage g(d.nx,d.ny);

        for (int i=0; i<d.ny; i++){
            for (int j=0; j<d.nx; j++){
                g.setElement(j, i, d.ux_l[i*d.nx+j]);
            }
        }
        g.save("scenario1.png");

        }else{

            resol = 60;
            accel = 0.016;
            t_end = 5;
            d.set_global(x_len, y_len, dia, resol, visco, accel, t_end, obst_x, obst_y);
            d.set_domain();
            d.initialize_2D_vector();
            d.set_vector_zero();
            d.set_omega();
            d.set_cords();
            d.find_lattice_units();
            d.set_relax_factor();

            d.initial_setup();
            cout << "Time step   -------------------> " << d.del_t << endl;
            cout << "Lattice spacing ---------------> " << d.del_x << endl;
            cout << "Acceleration in lattice units -> " << d.ax_l << endl;
            cout << "Relaxation rate ---------------> " << d.relax_fac << endl;

            for (int j=0; j<t_end; j=j+d.del_t){
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
            d.Normalizer();

            std::ofstream solution("solution.txt", std::ofstream::out);
            for (int i=0; i < d.ny; ++i){
                for (int j=0; j<d.nx; ++j ){
                    solution << d.ux_l[i*d.nx+j] << std::endl;
                }
            }
            solution.close();

            double max = *max_element(d.ux_l.begin(), d.ux_l.end());
            double min = *min_element(d.ux_l.begin(), d.ux_l.end());
            d.display(max);
            d.display(min);

            GrayScaleImage g(d.nx,d.ny);

            for (int i=0; i<d.ny; i++){
                for (int j=0; j<d.nx; j++){
                    g.setElement(j, i, d.ux_l[i*d.nx+j]);
                }
          }
            g.save("scenario2.png");
        }
}
