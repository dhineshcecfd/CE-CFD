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
    void set_global(double x_length, double y_length, double diameter, int resolution);
    void set_doamin ();
    void fill_domain();
//    void initial_vector();
    void initialize_vector();
    void calc_func_q();

public:
    int nx, ny, cellsize;
    double x_len, y_len, dia, resol;
    vector < vector < vector <double> > > vec;
    vector < vector <double> > func_q;
};

//    typedef vector<int> v1d;
//    typedef vector<v1d> v2d;
//    typedef vector<v2d> v3d;
//    v3d v(3, v2d(3, v1d(2, 4)));

void domain::calc_func_q(){
    for (size_t i=0; i<cellsize; i++){
        for (int j=0; j<9; j++){
//            for (int k=0; k<2; k++){
                func_q[i][j] = func_q[i][j] - (relax_fac * (func_q[i][j]-func_q_eq[i][j]) + (3.0 * omega_c * density * c_c * accel));
//            }
        }

    }
}

void domain::set_global(double x_length, double y_length, double diameter, int resolution){
    x_len = x_length;
    y_len = y_length;
    dia = diameter;
    resol = resolution;
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
    cout << "counter : " << vec.size() << endl;
}

//void domain::initial_vector(){
//        typedef vector<double> v1d;
//        typedef vector<v1d> v2d;
//        typedef vector<v2d> v3d;
//        v3d cell_pos(cellsize, v2d(9, v1d(9, 0.0)));
//}

void domain::set_doamin (){

    double temp;

    temp = resol/dia;
    nx = x_len * temp;
    ny = y_len * temp;

    cellsize = nx*ny;

    cout << "Cell size : " << cellsize << endl;
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
    double resol = 30;

    domain d;
    d.set_global(x_len, y_len, dia, resol);
    d.set_doamin();
    d.fill_domain();
    d.initialize_vector();


}
