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


struct values{
    double x;
    double y;
    double func;
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

void lattice::initial_setup(){
    density = 1.0;
    x_pos = 0;
    y_pos = 0;

    func = func_eq(density, x_pos, y_pos);
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

void lattice::omega(){

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


void lattice::lbm(){

    for (int t=0; t<t_end; t=t+det_t){
        for (size_t i=0; i<cell.size(); i++){

            //calc density and velocity
            density += func[i];
            velocity += (1/density) * func[i]*c[i];

            //collide
            func_eq[i] = omega[i] * density ( 1.0 + ( ( c[i] * velocity ) / c_s * c_s ) +
                                             ( pow ( ( c[i] * velocity ), 2.0) / ( 2.0 * pow( c_s, 4.0) ) ) -
                                              ( velocity * velocity / ( 2.0 * c_s * c_s ) ) );

            func[i] = func[i] - ((1/tow)*(func[i]) - func_eq[i]);

            //stream
            func(k,t+del_t) = func[i];

        }
    }

}

int main(int argc, char *argv[]){

    if(argc != 2){
        cout << "Not enough arguments" << endl;
    }

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
