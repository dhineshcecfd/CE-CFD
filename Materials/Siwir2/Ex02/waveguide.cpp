// reading a text file
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include <map>
#include "./myColsamm/Source/Colsamm.h"
#include <iomanip>

using namespace std;
using namespace ::_COLSAMM_;

ELEMENTS::Triangle my_element;

double delta, eps;

void call_function(double del, double epsilon){
    delta = del;
    eps = epsilon;

}

//myColsamm my_func()
double my_func(double x, double y){

    double k;
    double temp;

        temp = -50 * ((x*x) + (y*y));
        k = ((100+delta) * exp(temp)) - 100;
        return (k);


}

void local_matrix (vector<double>&corners, vector<double>&vertex0, vector<double>&vertex1, vector<double>&vertex2, vector<double>&x_cord, vector<double>&y_cord, vector< vector <double> >&my_local_matrix){

    vector < vector<double> >local_mass_matrix;

    //for assembling global matrix
    vector < map<int, double> > global_stiff(x_cord.size());
    vector < map<int, double> > global_mass(x_cord.size());

    //myCollsamm corners
    for (size_t i=0; i<vertex0.size(); i++){

        corners[0] = x_cord[vertex0[i]];
        corners[1] = y_cord[vertex0[i]];
        corners[2] = x_cord[vertex1[i]];
        corners[3] = y_cord[vertex1[i]];
        corners[4] = x_cord[vertex2[i]];
        corners[5] = y_cord[vertex2[i]];

        my_element(corners);

        //local stiffness matrix
        my_local_matrix = my_element.integrate(grad(v_()) * grad(w_()) - func<double>(my_func) * v_() * w_());

        //local mass matrix
        local_mass_matrix = my_element.integrate(v_() * w_());



            global_mass[vertex0[i]][vertex0[i]] += local_mass_matrix[0][0];
            global_mass[vertex0[i]][vertex1[i]] += local_mass_matrix[0][1];
            global_mass[vertex0[i]][vertex2[i]] += local_mass_matrix[0][2];
            global_mass[vertex1[i]][vertex0[i]] += local_mass_matrix[1][0];
            global_mass[vertex1[i]][vertex1[i]] += local_mass_matrix[1][1];
            global_mass[vertex1[i]][vertex2[i]] += local_mass_matrix[1][2];
            global_mass[vertex2[i]][vertex0[i]] += local_mass_matrix[2][0];
            global_mass[vertex2[i]][vertex1[i]] += local_mass_matrix[2][1];
            global_mass[vertex2[i]][vertex2[i]] += local_mass_matrix[2][2];



            global_stiff[vertex0[i]][vertex0[i]] += my_local_matrix[0][0];
            global_stiff[vertex0[i]][vertex1[i]] += my_local_matrix[0][1];
            global_stiff[vertex0[i]][vertex2[i]] += my_local_matrix[0][2];
            global_stiff[vertex1[i]][vertex0[i]] += my_local_matrix[1][0];
            global_stiff[vertex1[i]][vertex1[i]] += my_local_matrix[1][1];
            global_stiff[vertex1[i]][vertex2[i]] += my_local_matrix[1][2];
            global_stiff[vertex2[i]][vertex0[i]] += my_local_matrix[2][0];
            global_stiff[vertex2[i]][vertex1[i]] += my_local_matrix[2][1];
            global_stiff[vertex2[i]][vertex2[i]] += my_local_matrix[2][2];
}

    //Write M.txt
    std::cout << "Writing M.txt ......" << std::endl;
    std::ofstream solution("M.txt", std::ofstream::out);
    for(size_t h = 0; h < global_mass.size(); h++){
        for (map<int,double>::iterator ii=global_mass[h].begin(); ii!=global_mass[h].end(); ++ii){
            solution << h << " ";
            solution << (*ii).first << " ";
            solution << (*ii).second << " " << endl;
        }

    }
    solution.close();


    //Write A.txt
         std::cout << "Writing A.txt ......" << std::endl;
    std::ofstream solution1("A.txt", std::ofstream::out);
    for(size_t h = 0; h < global_mass.size(); h++){
        for (map<int,double>::iterator ii=global_stiff[h].begin(); ii!=global_stiff[h].end(); ++ii){
            solution1 << h << " ";
            solution1 << (*ii).first << " ";
            solution1 << (*ii).second << " " << endl;
        }

    }
    solution1.close();

    //Run Inverse iteration
    void Inverse_power_iteration(vector < map<int, double> > &global_mass, vector < map<int, double> > &global_stiff, vector<double>&x_cord, vector<double>&y_cord);
    Inverse_power_iteration(global_mass, global_stiff, x_cord, y_cord);

}


//display desired vector
void display(std::vector<double>&obj){

    for(int i=0;i<5;i++){
        cout<<"element of vec..."<<obj[i]<<endl;
    }
}


//display desired value
void display_value(double obj){

    cout<<"element of object..."<<obj << endl;
}

//Inverse iteration function
void Inverse_power_iteration(vector < map<int, double> > &global_mass, vector < map<int, double> > &global_stiff, vector<double>&x_cord, vector<double>&y_cord){

    double k = 1.0;
    double lamda_old = 1.0;
    double lamda = 2.0;
    int counter = 0;

    //initialize the u, r vectors
    vector <double> u(x_cord.size(), (1.0/sqrt(x_cord.size())));
    vector <double> r(x_cord.size(), 0.0);

    while(k > 1e-10){
        lamda_old = lamda;

        vector <double> fun (x_cord.size(), 0.0);
        for (size_t i = 0; i<x_cord.size(); i++){
            for (map<int,double>::iterator ii=global_mass[i].begin(); ii!=global_mass[i].end(); ++ii){
                fun[i] += (*ii).second * u[(*ii).first];
            }
        }

        vector <double> A_u(x_cord.size(), 0.0);
        for (size_t i = 0; i<x_cord.size(); i++){
            for (map<int,double>::iterator ii=global_stiff[i].begin(); ii!=global_stiff[i].end(); ++ii){
                A_u[i] += (*ii).second * u[(*ii).first];
            }
        }

        // r=Au-b
        for (size_t i=0; i<fun.size(); i++){
        r[i] = A_u[i] - fun[i];

        }

        //call cg method
        void conjugate_gradient(vector<double>&function, vector<map<int, double>>&g_stiff, vector<double>&u,vector <double> &r);
        conjugate_gradient(fun, global_stiff, u,r);

        double u_temp = 0.0;
        double u_norm = 0.0;

        //||u||
        for (size_t i=0; i<u.size(); i++){
            u_temp += u[i]*u[i];
        }
        u_norm = sqrt(u_temp);

        for (size_t i=0; i<u.size(); i++){
            u[i] = u[i]/u_norm;
        }

        vector <double> temp_A_u(x_cord.size());
        for (size_t i = 0; i<x_cord.size(); i++){
            for (map<int,double>::iterator ii=global_stiff[i].begin(); ii!=global_stiff[i].end(); ++ii){
                temp_A_u[i] += (*ii).second * u[(*ii).first];
            }
        }

        double lamda_up = 0.0;
        for (size_t i=0; i<u.size(); i++){
            lamda_up += u[i] * temp_A_u[i];
        }

        vector <double> temp_M_u(x_cord.size());
        for (size_t i = 0; i<x_cord.size(); i++){
            for (map<int,double>::iterator ii=global_mass[i].begin(); ii!=global_mass[i].end(); ++ii){
                temp_M_u[i] += (*ii).second * u[(*ii).first];
            }
        }

        double lamda_down = 0.0;
        for (size_t i=0; i<u.size(); i++){
            lamda_down += u[i] * temp_M_u[i];
        }

        lamda = lamda_up / lamda_down;

        k = abs((lamda - lamda_old)/lamda_old);
        ++counter;

    }

    //Write eigenmode.txt
         std::cout << "Writing eigenmode.txt ......" << std::endl;
    std::ofstream solution1("eigenmode.txt", std::ofstream::out);
    for(size_t i=0; i<u.size(); i++){
            solution1 << x_cord[i] << " ";
            solution1 << y_cord[i] << " ";
            solution1 << u[i] << " " << endl;
        }
    solution1.close();

    cout << "lamda : " << setprecision(10) << lamda << endl;
    cout << "No of iterations to converge : " << counter << endl;
}





//cg function
void conjugate_gradient(vector<double>&function, vector<map<int, double>>&g_stiff, vector<double>&u,vector <double> &r){

    vector <double> d(function.size(), 0.0);
    for (size_t i=0; i<function.size(); i++){
    }

    int counter = 1;
    double beta = 0.0;
    double alpha = 0.0;
    double l2_norm = 1.0;

    // d = -r
    for(size_t i = 0; i < function.size(); i++ ){
        d[i] = -1.0*r[i];
    }

    while(l2_norm > eps){

            vector<double>ad(function.size(), 0.0);
            for(size_t i=0; i<function.size(); i++){
                for (map<int,double>::iterator aditr=g_stiff[i].begin(); aditr!=g_stiff[i].end(); ++aditr){
                ad[i] += (*aditr).second * d[(*aditr).first];
            }
        }

            double rho = 0.0;
            for (size_t i=0; i<function.size(); i++){
                rho += d[i] * ad[i];
            }

            double al_num = 0.0;
            for (size_t i=0; i<function.size(); i++){
                al_num += r[i]*d[i];
            }

            alpha = (al_num / (-1.0*rho));

            for (size_t i=0; i<function.size(); i++){
                u[i] += (alpha * d[i]);
            }

            for (size_t i=0; i<function.size(); i++){
                r[i] += alpha * ad[i];
            }

            double beta_num = 0.0;
            for(size_t i=0; i<function.size(); i++){
                beta_num += r[i] * ad[i];
            }

            beta = beta_num / rho;

            for(size_t i=0; i<function.size(); i++){
                d[i] = (-1.0*r[i]) + (beta * d[i]);
            }

            double temp1 = 0.0;
            for (size_t i=0; i<function.size(); i++){
                temp1 += r[i] * r[i];
            }

            l2_norm = sqrt(temp1);
            counter++;
    }

}

//k_square
void k_square(std::vector<double>&k_sq, std::vector<double>&x_cord, std::vector<double>&y_cord){
    double temp;
    for(size_t i=0; i<x_cord.size(); i++){
        temp = -50 * ((x_cord[i]*x_cord[i]) + (y_cord[i]*y_cord[i]));
        k_sq[i] = ((100+delta) * exp(temp)) - 100;
    }
}


int main (int argc , char *argv[]) {

    if (argc <= 2)
    {
        cout << "Too few arguments" << std::endl;
    }

    double delt = atof(argv[1]);
    double epsi = atof(argv[2]);

    call_function(delt, epsi);

    vector<double>x_cord, y_cord;
    vector<double>vertex0, vertex1, vertex2;
    vector< vector <double> > my_local_matrix;
    vector<double>corners(6, 0.0);
    string line;
    int option;

    cout << "please enter 1-unit_circle.txt or 2-unit_circle_fine.txt to proceed : ";
    cin >> option;

    if (option == 1){

    ifstream myfile ("unit_circle.txt");
    if (myfile.is_open())
        {
        int line_no = 0;
        while (getline (myfile,line))
            {
            //leaving string lines
            if (1 < line_no && line_no < 1041)
                {
                double col1 = 0.0;
                double col2 = 0.0;
                double col3 = 0.0;
                //read column wise
                istringstream ss(line);
                ss >> col1 >> col2 >> col3;
                x_cord.push_back(col2);
                y_cord.push_back(col3);

//                waveguide(col2, col3);
                }
            else if(line_no > 1042){
                double ver0 = 0.0;
                double ver1 = 0.0;
                double ver2 = 0.0;
                //read column wise
                istringstream ss1(line);
                ss1 >> ver0 >> ver1 >> ver2;
                vertex0.push_back(ver0);
                vertex1.push_back(ver1);
                vertex2.push_back(ver2);
            }
            line_no = line_no+1;
            }
        myfile.close();
        }
    else cout << "Unable to open file";
    }else{
        ifstream myfile ("unit_circle_fine.txt");
        if (myfile.is_open())
            {
            int line_no = 0;
            while ( getline (myfile,line) )
                {
                //leaving string lines
                if (1 < line_no && line_no < 12184)
                    {
                    double col1 = 0.0;
                    double col2 = 0.0;
                    double col3 = 0.0;
                    //read column wise
                    istringstream ss(line);
                    ss >> col1 >> col2 >> col3;
                    x_cord.push_back(col2);
                    y_cord.push_back(col3);

    //                waveguide(col2, col3);
                    }
                else if(line_no > 12185){
                    double ver0 = 0.0;
                    double ver1 = 0.0;
                    double ver2 = 0.0;
                    //read column wise
                    istringstream ss1(line);
                    ss1 >> ver0 >> ver1 >> ver2;
                    vertex0.push_back(ver0);
                    vertex1.push_back(ver1);
                    vertex2.push_back(ver2);
                }
                line_no = line_no+1;
                }
            myfile.close();
            }
        else cout << "Unable to open file";
    }


    vector<double>k_sq(x_cord.size());

    //k_square function
    k_square(k_sq, x_cord, y_cord);

    //calling local_matrix function
    local_matrix(corners, vertex0, vertex1, vertex2, x_cord, y_cord, my_local_matrix);

    //Write k_sq.txt
         std::cout << "Writing ksq.txt ......" << std::endl;
    std::ofstream solution("ksq.txt", std::ofstream::out);

    // for loop to write the output
    for (size_t i=0; i < x_cord.size(); i++){
            solution << x_cord[i] << "\t";
            solution << y_cord[i] << "\t";
            solution << k_sq[i] << std::endl;

    //    solution << std::endl;
    }
    solution.close();



    return 0;
}
