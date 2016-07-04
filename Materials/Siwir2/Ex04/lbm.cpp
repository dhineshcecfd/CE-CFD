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
    void set_c_c_s();


public:
    //3D array for stroing data
    vector < vector < vector <double> > > cell;
    double density;
    double velocity;
    vector <double> omega;
    vector <double> c;
    double c_s;
};

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

void lattice::set_c_c_s(){
    c[0] =
}

void lattice::lbm(){

    for (int t=0; t<t_end; t=t+det_t){
        for (size_t i=0; i<cell.size(); i++){

            //calc density and velocity
            density += func[i];
            velocity += func[i]*c[i];

            //collide
            func_eq[i] = omega[i] * density (1+ ((c[i]*velocity)/c_s*c_s) + (pow((c[i]*velocity),2) / (2*pow(c_s,4))) - (velocity * velocity / (2*c_s*c_s)));

            func[i] = func[i] - ((1/tow)*(func[i]) - func_eq[i]);

            //stream
            func(k,t+del_t) = func[i];

        }
    }

}

int main(){
    lattice l;
    values vl;
    l.set_value(0.5, 0.6, 2.33);

}
