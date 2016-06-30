#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <ctype.h>
using namespace std;
int main()
{
    vector<double>x_cord, y_cord;
    string line;
    ifstream myfile ("unit_circle.txt");
//    myfile >> count;
    if (myfile.is_open())
        {
        getline (myfile,line);
                getline (myfile,line);
        //while (myfile.good())
            //{
            while(istringstream ss(line)){
            //leaving string lines
//            if (1 < line_no && line_no < 1041)

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

        }
            for (size_t i=0; i<x_cord.size(); i++)
                cout << "X_cord : " << x_cord[i] << endl;
            cout<< x_cord.size() << endl;

            //Write A.txt
                 std::cout << "Writing .txt ......" << std::endl;
            std::ofstream solution1("x_cord.txt", std::ofstream::out);
            for(size_t h = 0; h < x_cord.size(); h++){
                    solution1 << x_cord[h] << " " << endl;
                }


            solution1.close();
    return 0;
}
