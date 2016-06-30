#include<iostream>
#include<fstream>
#include<map>
#include<string>
#include<sstream>
#include<vector>
#include<math.h>
#include <stdexcept>
#include <sys/stat.h>

using namespace std;

class ParameterReader{

public:
    ParameterReader(){}
    bool readParameters(const std::string& filename);
    inline bool IsDefined(const std::string& key)const ;
    template< typename D >
    inline void GetParameter(const std::string& key, D &value)const ;
    map<std::string, string >param_list;

};

bool ParameterReader :: readParameters(const std::string &filename) {

    std::ifstream file_name(filename);
    if (file_name.is_open()){

        string line;

        while (getline(file_name, line)){
            string col1;
            string col2;
            stringstream ss(line);
            ss >> col1 >> col2;
            param_list.insert (pair <string, string> (col1, col2));
        }
    }
    file_name.close();
    return true;
}

template< typename D >
inline void ParameterReader ::  GetParameter(const std::string& key, D &value) const {
    stringstream ss;
    string dummy = param_list.at(key);
    ss << dummy;
          ss >> value ;

}

inline bool ParameterReader:: IsDefined(const std::string& key)const {
    size_t ct = 0;
    for(map<string, string> :: const_iterator itr = param_list.cbegin() ;itr!= param_list.cend(); itr++ ){
        if( key == (*itr).first){
            return true;
            break;
        }
        else{
            ++ct;
            if(ct == param_list.size()){
                return false;
            }
        }
    }
}

struct triples{
    double x;
    double y;
    double z;
};

struct parameters{
    double rcut;
    double dt;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    double eps;
    double sigma;
    int vs;
    double ts;
    double te;
};

class LinkedCell
{
public:
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    LinkedCell() {}
    size_t nc[3];
    vector<triples>pos;
    vector<triples>vel;
    vector<triples>force;
    vector<triples>pf;
    vector<triples>l;
    vector<parameters>p;
    vector<double>mass;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Getdata(const std::string& filename);
    void export_parameters(double a,double b,double c,double d,double e,double f,double g,double h,double i,double j,double k,double l1, int v );
    void setpos(double a,double b, double c);
    void setvel(double a,double b, double c);
    void setforce(double a,double b, double c);
    void setpf(double a,double b, double c);
    void force_cal(int a,int b);
    void pos_update(size_t a);
    void vel_update(size_t a);
    void integerator();
    void display(vector<triples>p);
    void setf_zero(size_t i);
    void setpf_zero(size_t i);
    void initiate_f();
    void simulation_space();
    void cell_filling();
    void force_calc(int i, vector<int> &c);
    void neighbour_cells(int i);
    void print_position(std::ofstream &);
    void print_force(std::ofstream &);
    void print_velocity(std::ofstream &);
    void print_mass(std::ofstream &);
    void printcell();
    void resize_cell();
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
private:
    
    double cell_space[3];
    vector<vector<int>>cell;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//freeing the particles inforamtion in cells after every timestep
void LinkedCell :: resize_cell(){
    cell.resize(0);
    cell.resize(nc[0]*nc[1]*nc[2]);
}

//finding the neighbour cells
void LinkedCell::neighbour_cells(int i){

    int x_cell,y_cell,z_cell,x_right, x_left, y_up, y_down, z_front, z_back;

    x_cell = (pos[i].x / cell_space[0]);
    y_cell = (pos[i].y / cell_space[1]);
    z_cell = (pos[i].z / cell_space[2]);

    x_right = (x_cell+1)%nc[0];
    x_left  = (x_cell+nc[0]-1)%nc[0];
    y_up    = (y_cell+1)%nc[1];
    y_down  = (y_cell+nc[1]-1)%nc[1];
    z_front = (z_cell+1)%nc[2];
    z_back  = (z_cell+nc[2]-1)%nc[2];

    for (int ne=0; ne<3; ne++){
        if(ne == 1){
            z_cell = z_front;
        }if(ne == 2){
            z_cell = z_back;
        }

        int ic = x_cell + (y_cell)*nc[0] + (z_cell)*nc[0]*nc[1];
        force_calc(i, cell[ic]);

        int ie = x_right + (y_cell)*nc[0] + (z_cell)*nc[0]*nc[1];
        force_calc(i, cell[ie]);

        int iw = x_left + (y_cell)*nc[0] + (z_cell)*nc[0]*nc[1];
        force_calc(i, cell[iw]);

        int in = x_cell + (y_up)*nc[0] + (z_cell)*nc[0]*nc[1];
        force_calc(i, cell[in]);

        int is = x_cell + (y_down)*nc[0] + (z_cell)*nc[0]*nc[1];
        force_calc(i, cell[is]);

        int inw = x_left + (y_up)*nc[0] + (z_cell)*nc[0]*nc[1];
        force_calc(i, cell[inw]);

        int ine = x_right + (y_up)*nc[0] + (z_cell)*nc[0]*nc[1];
        force_calc(i, cell[ine]);

        int isw = x_left + (y_down)*nc[0] + (z_cell)*nc[0]*nc[1];
        force_calc(i, cell[isw]);

        int ise = x_right + (y_down)*nc[0] + (z_cell)*nc[0]*nc[1];
        force_calc(i, cell[ise]);
    }
}

//calculating force for concecutive neighbour particels
void LinkedCell :: force_calc(int i, vector<int> &c){
    
    for(size_t it = 0; it < c.size();it++){
        if(i!=c[it])force_cal(i,c[it]);  
    }
     
}

//finding the cells of corresponding particles
void LinkedCell::cell_filling(){

    int cell_no,x_cell,y_cell,z_cell;
    for (size_t i =0; i<pos.size(); i++){
        x_cell = (pos[i].x / cell_space[0]);
        y_cell = (pos[i].y / cell_space[1]);
        z_cell = (pos[i].z / cell_space[2]);
        cell_no = x_cell + (y_cell)*nc[0] + (z_cell)*nc[0]*nc[1];
        cell[cell_no].push_back(i);
    }
                 

}
 
//calculating simulation space or domain
void LinkedCell :: simulation_space(){

   for(int i = 0; i<3;i++)nc[i] = 0;
   nc[0] = (p[0].xmax - p[0].xmin) / p[0].rcut;
   nc[1] = (p[0].ymax - p[0].ymin) / p[0].rcut;
   nc[2] = (p[0].zmax - p[0].zmin) / p[0].rcut;
   
   cell_space[0] = ((p[0].xmax - p[0].xmin) /nc[0]);
   cell_space[1] = ((p[0].ymax - p[0].ymin) /nc[1]);
   cell_space[2] = ((p[0].zmax - p[0].zmin) /nc[2]);
   cell.resize(nc[0]*nc[1]*nc[2]);
}

void LinkedCell :: export_parameters(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j, double k, double l1, int v){
    parameters o;
    o.eps = a;
    o.sigma = b;
    o.rcut = c;
    o.dt = d;
    o.xmax = e;
    o.xmin = f;
    o.ymax = g;
    o.ymin = h;
    o.zmax = i;
    o.zmin = j;
    o.ts = k;
    o.te = l1;
    o.vs = v;
    p.push_back(o);

}

//calculating initialize force
void LinkedCell :: initiate_f(){
    for(size_t i = 0; i < pos.size();i++){
            for(size_t j = 0; j < pos.size();j++){
                if(i==j){}
                else{
                    force_cal(i,j);
                }
            }
        }

}

//set force to "zero"
void LinkedCell::setf_zero(size_t i){

        force[i].x = 0;
        force[i].y = 0;
        force[i].z = 0;

}

//set old force to "zero"
void LinkedCell::setpf_zero(size_t i){

        pf[i].x = 0;
        pf[i].y = 0;
        pf[i].z = 0;

}

//write particel position into vtk file
void LinkedCell::print_position(std::ofstream &file){

    for(size_t i=0; i < pos.size(); ++i){
        file<<pos[i].x<<" "<<pos[i].y<<" "<<pos[i].z<<"\n";
    }
}

//write particel force into vtk file
void LinkedCell::print_force(std::ofstream &file){

    for(size_t i=0; i < pos.size(); ++i){
        file<<force[i].x<<" "<<force[i].y<<" "<<force[i].z<<"\n";
    }
}

//write particel velocity into vtk file
void LinkedCell::print_velocity(std::ofstream &file){

    for(size_t i=0; i < pos.size(); ++i){
        file<<vel[i].x<<" "<<vel[i].y<<" "<<vel[i].z<<"\n";
    }
}

//write particel mass into vtk file
void LinkedCell::print_mass(std::ofstream &file){

    for(size_t i=0; i < pos.size(); ++i){
        file<< mass[i] <<"\n";
    }
}

//function for MD simulation
void LinkedCell :: integerator(){
    
    //updating particle position
    for(size_t i = 0; i< pos.size();i++){
        pos_update(i);
    }

    //updating neighbour cells
    for(size_t i = 0; i < pos.size();i++){
        setf_zero(i);
        neighbour_cells(i);
    }

    //updating particle velocities
    for(size_t i = 0; i< pos.size();i++){
        vel_update(i);
    }
        
}

//reading .dat file and store into vector
void LinkedCell ::Getdata(const std::string &filename){
    std::ifstream file_name(filename);
    if (file_name.is_open()){

        string line;

        while (getline(file_name, line)){
            double col1, col2, col3, col4, col5, col6, col7;
            stringstream ss(line);
            ss >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7;
            setpos(col2,col3,col4);
            setvel(col5,col6,col7);
            setforce(0.0,0.0,0.0);
            setpf(0.0,0.0,0.0);
            mass.push_back(col1);

        }
    }

}

//set particle position from .dat file
void LinkedCell ::setpos(double a,double b, double c){
    triples t;
    t.x = a;
    t.y = b;
    t.z = c;
    pos.push_back(t);

}

//set particle velocity from .dat file
void LinkedCell ::setvel(double a, double b, double c){
    triples t;
    t.x = a;
    t.y = b;
    t.z = c;
    vel.push_back(t);

}

//set initial force to "zero"
void LinkedCell ::setforce(double a, double b, double c){
    triples t;
    t.x = a;
    t.y = b;
    t.z = c;
    force.push_back(t);

}

//set initial old force to "zero"
void LinkedCell ::setpf(double a, double b, double c){
    triples t;
    t.x = a;
    t.y = b;
    t.z = c;
    pf.push_back(t);

}

//updating particle position for every particle
void LinkedCell::pos_update(size_t a){

    double len_dom = p[0].xmax - p[0].xmin;
        pos[a].x += p[0].dt*vel[a].x + (0.5*p[0].dt*p[0].dt*force[a].x)/mass[a];

            if(pos[a].x < p[0].xmin){
                double diff = fabs(p[0].xmin - pos[a].x);
                pos[a].x = p[0].xmax - fmod(diff,len_dom);
            }

            if(pos[a].x > p[0].xmax){
                double diff = fabs(pos[a].x - p[0].xmax);
                pos[a].x = p[0].xmin + fmod(diff,len_dom);
            }

            pos[a].y += p[0].dt*vel[a].y + (0.5*p[0].dt*p[0].dt*force[a].y)/mass[a];

            if(pos[a].y < p[0].ymin){
                double diff = fabs(p[0].ymin - pos[a].y);
                pos[a].y = p[0].ymax - fmod(diff,len_dom);
            }

            if(pos[a].y > p[0].ymax){
                double diff = fabs(pos[a].y - p[0].ymax);
                pos[a].y = p[0].ymin + fmod(diff,len_dom);
            }

            pos[a].z += p[0].dt*vel[a].z + (0.5*p[0].dt*p[0].dt*force[a].z)/mass[a];

            if(pos[a].z < p[0].zmin){
                double diff = fabs(p[0].zmin - pos[a].z);
                pos[a].z = p[0].zmax - fmod(diff,len_dom);
            }

            if(pos[a].z > p[0].zmax){
                double diff = fabs(pos[a].z - p[0].zmax);
                pos[a].z = p[0].zmin + fmod(diff,len_dom);
            }

            pf[a].x = force[a].x;
            pf[a].y = force[a].y;
            pf[a].z = force[a].z;
}

//updating force on each particle
void LinkedCell::force_cal(int a, int b){

    //finding distance between two particles
    double fx,fy,fz,rx,ry,rz,r_s,r,lx, ly, lz;
    rx = pos[b].x - pos[a].x;
    ry = pos[b].y - pos[a].y;
    rz = pos[b].z - pos[a].z;
    r_s = rx*rx + ry*ry + rz*rz;
    r = sqrt(r_s);
    lx = p[0].xmax - p[0].xmin;
    ly = p[0].ymax - p[0].ymin;
    lz = p[0].zmax - p[0].zmin;

    //finding distance between two boundary particles
    if(r > p[0].rcut){

//        rx = (rx/fabs(rx)) * ((p[0].xmax - p[0].xmin) - fabs(rx));
//        ry = (ry/fabs(ry)) * ((p[0].xmax - p[0].xmin) - fabs(ry));
//        rz = (rz/fabs(rz)) * ((p[0].xmax - p[0].xmin) - fabs(rz));

        if(rx < 0){
            rx = -(lx - fabs(rx));
        }
        if(rx > 0){
            rx = (lx - fabs(rx));
        }
        if(ry < 0){
            ry = -(ly - fabs(ry));
        }
        if(ry > 0){
            ry = (ly - fabs(ry));
        }
        if(rz < 0){
            rz = -(lz - fabs(rz));
        }
        if(rz > 0){
            rz = (lz - fabs(rz));
        }

        r_s = rx*rx + ry*ry + rz*rz;
        r = sqrt(r_s);
    }

    //calculaitng force between particles that distance less than r_cut
    if(r <= p[0].rcut){

            double tempx = pow(p[0].sigma/r,6);
            fx = 24*p[0].eps*((1/(r*r))*tempx*(1 - 2.0*tempx))*rx;
            force[a].x += fx;
            fx = 0.0;

            double tempy = pow(p[0].sigma/r,6);
            fy = 24*p[0].eps*((1/(r*r))*tempy*(1 - 2.0*tempy))*ry;
            force[a].y += fy;
            fy = 0.0;

            double tempz = pow(p[0].sigma/r,6);
            fz = 24*p[0].eps*((1/(r*r))*tempz*(1 - 2.0*tempz))*rz;
            force[a].z += fz ;
            fz = 0.0;
    }
}

//updating particle velocity on each time steps
void LinkedCell::vel_update(size_t a){
    vel[a].x += (0.5*p[0].dt*(force[a].x+pf[a].x))/mass[a];
    vel[a].y += (0.5*p[0].dt*(force[a].y+pf[a].y))/mass[a];
    vel[a].z += (0.5*p[0].dt*(force[a].z+pf[a].z))/mass[a];
}


int main(int argc , char *argv[]) {

    ParameterReader p;
    p.readParameters(argv[1]);
    double eps,sigma,r_cut,dt,ts,te,xmax,xmin,ymax,ymin,zmax,zmin;
    string name;
    int v;

    //reading parameter data with key
    const string e = "epsilon";
    p.GetParameter(e,eps);
    const string f = "sigma";
    p.GetParameter(f,sigma);
    const string g = "r_cut";
    p.GetParameter(g,r_cut);
    const string h = "delta_t";
    p.GetParameter(h,dt);
    const string nam = "name";
    p.GetParameter(nam, name);
    const string a = "x_max";
    const string b = "x_min";
    const string c = "y_min";
    const string d = "y_max";
    const string i = "z_min";
    const string j = "z_max";
    const string k = "t_start";
    const string l = "t_end";
    const string m = "vis_space";
    p.GetParameter(a,xmax);
    p.GetParameter(b,xmin);
    p.GetParameter(c,ymin);
    p.GetParameter(d,ymax);
    p.GetParameter(i,zmin);
    p.GetParameter(j,zmax);
    p.GetParameter(k,ts);
    p.GetParameter(l,te);
    p.GetParameter(m,v);

    //creating object
    LinkedCell dr;

    dr.export_parameters(eps,sigma,r_cut,dt,xmax,xmin,ymax,ymin,zmax,zmin,ts,te,v);
    dr.Getdata(argv[2]);

    //writing vtk files
    std::ostringstream os;
    os <<name<<0<<".vtk";
    std::ofstream file( os.str().c_str() );
    os.seekp(0);
    file<<"# vtk DataFile Version 3.0"<<"\n"
            <<"SiwiRVisFile"<<"\n"
            <<"ASCII\n"
            <<"DATASET UNSTRUCTURED_GRID\n"
            <<"POINTS "<< dr.pos.size() <<" DOUBLE\n";

    //writing particle position in vtk
    dr.print_position(file);

    file<<"POINT_DATA "<< dr.pos.size() <<"\n"
            <<"SCALARS mass double\n"
            <<"LOOKUP_TABLE default\n";

    //writing particle mass in vtk
    dr.print_mass(file);

    file<<"VECTORS force double\n";

    //writing force in vtk
    dr.print_force(file);

    file<<"VECTORS velocity double\n";

    //write velocity in vtk
    dr.print_velocity(file);

    file.close();
/////////////////////////////////////////////////////////////////////////////////////

    //call initialize force function
    dr.initiate_f();

    //call simulation domain calculation function
    dr.simulation_space();

    //call cell filling function in Linkedcell class
    dr.cell_filling();

    //writing vtk files for every time steps
    size_t vtk_t (1);
    for(double t = ts; t < te; t+=dt,++vtk_t){

        //call integerator function
        dr.integerator();

        if(v != 0){
            if((vtk_t%v==0)){
                std::ostringstream os1;
                os1 <<name<<vtk_t<<".vtk";
                std::ofstream file1( os1.str().c_str() );
                os1.seekp(0);
                file1<<"# vtk DataFile Version 3.0"<<"\n"
                        <<"SiwiRVisFile"<<"\n"
                        <<"ASCII\n"
                        <<"DATASET UNSTRUCTURED_GRID\n"
                        <<"POINTS "<<dr.pos.size()<<" DOUBLE\n";
                dr.print_position(file1);
                file1<<"POINT_DATA "<<dr.pos.size()<<"\n"
                        <<"SCALARS mass double\n"
                        <<"LOOKUP_TABLE default\n";
                dr.print_mass(file1);
                file1<<"VECTORS force double\n";
                dr.print_force(file1);
                file1<<"VECTORS velocity double\n";
                dr.print_velocity(file1);
                file1.close();
                }
            }
        dr.resize_cell();
        dr.cell_filling();
    }

    return 0;
}
































