#include<iostream>
#include<vector>
#include <cstdlib>


using namespace std;
struct values{
    double x;
    double y;
};

class domain{
public:
    void set_initial(double xval, double yval);
    void set_global(double x_length, double y_length, double diameter, int resolution, double viscosity, double acceleration, int end_time);
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
    void stream_step(int i);

public:
    int nx, ny, cellsize, t_end;
    double x_len, y_len, dia, resol, relax_fac, nu, accel, ux = 0.0, uy = 0.0, density = 1.0;
    vector < vector < vector <double> > > vec;
    vector < vector <double> > func_q;
    vector < vector <double> > func_q_eq;
    double omega[9];
    double cords[9][2];
    //variables in lattice units
    double del_x, del_t, nu_l, ux_l = 0.0, uy_l = 0.0, ax_l, ay_l=0.0;

};

void domain::stream_step(int i){

//    for (int i=1; i<ny-1; i++){
//        for (int k=1; k<nx-1; k++){

//            func_q [(i*nx)+k][0]        = func_q[i*nx+k][0];
//            func_q [(i*nx)+(k+1)][1]    = func_q[i*nx+k][1];
//            func_q [((i+1)*nx)+k][2]    = func_q[i*nx+k][2];
//            func_q [(i*nx)+k-1][3]      = func_q[i*nx+k][3];
//            func_q [(i-1)*nx+k][4]      = func_q[i*nx+k][4];
//            func_q [(i+1)*nx+k+1][5]    = func_q[i*nx+k][5];
//            func_q [(i+1)*nx+k-1][6]    = func_q[i*nx+k][6];
//            func_q [(i-1)*nx+k-1][7]    = func_q[i*nx+k][7];
//            func_q [(i-1)*nx+k+1][8]    = func_q[i*nx+k][8];
//        }
//    }
    func_q [i][0]           = func_q[i][0];
    func_q [i+1][1]         = func_q[i][1];
    func_q [i+nx][2]        = func_q[i][2];
    func_q [i-1][3]         = func_q[i][3];
    func_q [i-nx][4]        = func_q[i][4];
    func_q [i+nx+1][5]      = func_q[i][5];
    func_q [i+nx-1][6]      = func_q[i][6];
    func_q [i-nx-1][7]      = func_q[i][7];
    func_q [i-nx+1][8]      = func_q[i][8];

}

void domain::initialize_2D_vector(){

    func_q.resize( cellsize , vector<double>( 9 , 0.0 ) );
    func_q_eq.resize( cellsize , vector<double>( 9 , 0.0 ) );
}

void domain::find_metric_units(){
    x_len = del_x * (nx*100);
    del_t = 2;

    //find metric viscosity
    nu = nu_l * del_x * del_x / del_t;

    //find metric velocity
    ux = ux_l * del_x / del_t;
    uy = uy_l * del_x / del_t;

    //find matric acceleration
    accel = ax_l * del_x / (del_t * del_t);
}

void domain::calc_funcq_eq(int j){

    for (int i=0; i<9; i++){
        func_q_eq[j][i] = omega[i] * density * (1.0 + (3.0*((cords[i][0]*ux_l) + (cords[i][1]*uy_l))) +
                       (4.5*((cords[i][0]*ux_l) + (cords[i][1]*uy_l)) * ((cords[i][0]*ux_l) + (cords[i][1]*uy_l))) -
                        (1.5*((ux_l*ux_l)+(uy_l*uy_l))));

//    cout << func_q_eq[j][i] << endl;
    }
}

void domain::calc_func_q(int i){

    for (int j=0; j<9; j++){
        func_q[i][j] = func_q[i][j] - (relax_fac * (func_q[i][j]-func_q_eq[i][j]) +
                                       (3.0 * omega[i] * density * (cords[j][0] * ax_l) * (cords[j][1]*ay_l)));
//        display(func_q[i][j]);
    }
}

void domain::calc_density(int j){
    for (int i=0; i<9; i++){
        density += func_q[j][i];
    }
    display(density);
}

void domain::calc_velocity(int j){

//    display(ux_l);
//    display(uy_l);

        ux_l += (func_q[j][1]*cords[1][0] + func_q[j][5]*cords[5][0] + func_q[j][8]*cords[8][0] -
                func_q[j][6]*cords[6][0] - func_q[j][3]*cords[3][0] - func_q[j][7]*cords[7][0])/density;

        uy_l += (func_q[j][6]*cords[6][0] + func_q[j][5]*cords[5][0] + func_q[j][2]*cords[2][0] -
                func_q[j][7]*cords[7][0] - func_q[j][4]*cords[4][0] - func_q[j][7]*cords[7][0])/density;

        display(ux_l);
        display(uy_l);
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
    del_t = 0.00002;

    //find lattice viscosity
    nu_l = nu * del_t / (del_x*del_x);

    //find lattice velocity
//    ux_l = ux * del_t / del_x;
//    uy_l = uy * del_t / del_x;

    //find lattice acceleration
    ax_l = accel * del_t * del_t / del_x;

}

void domain::find_lattice_velocities(){

    //find lattice velocity
    ux_l = ux * del_t / del_x;
    uy_l = uy * del_t / del_x;
}

void domain::set_relax_factor(){
    relax_fac = 1/((3*nu_l)+0.5);
//    display(relax_fac);
}

void domain::set_omega(){
    omega[0] = 0.027777778; //SE
    omega[2] = 0.027777778; //SW
    omega[6] = 0.027777778; //NW
    omega[8] = 0.027777778; //NE
    omega[1] = 0.111111111; //W
    omega[3] = 0.111111111; //S
    omega[5] = 0.111111111; //E
    omega[7] = 0.111111111; //N
    omega[4] = 0.444444444; //C
}



void domain::set_global(double x_length, double y_length, double diameter, int resolution, double viscosity, double acceleration, int end_time){
    x_len = x_length;
    y_len = y_length;
    dia = diameter;
    resol = resolution;
    nu = viscosity;
    accel = acceleration;
    t_end = end_time;
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
//    cout << "counter : " << vec.size() << endl;
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


    double x_len = 6;  //in meters
    double y_len = 2;    //in m
    double dia = 0.5;
    double visco = 10e-6;
    double resol = 30;
    double accel = 0.01;
    int t_end = 3;

    domain d;
    d.set_global(x_len, y_len, dia, resol, visco, accel, t_end);
    d.set_domain();
    d.set_omega();
    d.fill_domain();
    d.set_cords();

    d.set_relax_factor();
    d.initialize_2D_vector();
    d.calc_density(0);
    d.calc_velocity(0);
    d.find_lattice_units();

    for (int i=(d.nx+1); i<d.cellsize-(d.nx+1); i++){
        d.stream_step(i);
        d.calc_funcq_eq(i);
        d.calc_func_q(i);
//    int i=1;
        d.calc_density(i);
        d.calc_velocity(i);
//        d.calc_funcq_eq(i);
//        d.calc_func_q(i);
//        d.find_metric_units();
    }
//    d.initialize_vector();
}
