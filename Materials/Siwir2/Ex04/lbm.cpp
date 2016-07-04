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

public:
    //3D array for stroing data
    vector < vector < vector <double> > > cell;
};

void lattice::set_value(double x_val, double y_val, double func_val){
    values v;
    v.x = x_val;
    v.y = y_val;
    v.func = func_val;
    std::cout << "values of : " << v.x << v.y << v.func << endl;
}

int main(){
    lattice l;
    values vl;
    l.set_value(0.5, 0.6, 2.33);

}
