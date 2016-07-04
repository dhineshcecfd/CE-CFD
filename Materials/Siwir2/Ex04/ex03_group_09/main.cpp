
#include"LBM.h"
#include"Fileread.h"
#include"Timing.h"


int main( int argc, char *argv[]) {
	if(argc!=2){
			std::cerr<<"Not enough arguments"<<std::endl;
			throw std::runtime_error("Not enough arguments ");
	}
	
	std::ifstream file(argv[1]);

	if(file.is_open())
		std::cout<<"File open"<<std::endl;
	else{
		std::cout<<"Error opening file"<<std::endl;
		throw std::runtime_error("Error opening file");
	}
	uint_t x,y,timesteps,vtk_step;
	real_t omega;
	FileReader file_rdr;

	//Reading data from params

	file_rdr.readParameterFile(file);
	x = file_rdr.getParameter<uint_t>("sizex")+2;
	y = file_rdr.getParameter<uint_t>("sizey")+2;
	timesteps = file_rdr.getParameter<uint_t>("timesteps");
	vtk_step = file_rdr.getParameter<uint_t>("vtk_step");
	omega = file_rdr.getParameter<real_t>("omega");

	if(omega < 0.5 || omega > 1.95){
		std::cout<<"Omega not in the stablity range"<<std::endl;
		throw std::runtime_error("Omega not in the stablity range");
	}
	PDF_Field src(x,y),dst(x,y);
	VelocityField v(x,y);
	DensityField d(x,y);

	Lbm solver;

	real_t t_alpha[9] =
	{
			4.0/9.0,
			1.0/9.0,
			1.0/9.0,
			1.0/9.0,
			1.0/9.0,
			1.0/36.0,
			1.0/36.0,
			1.0/36.0,
			1.0/36.0,
	};
	real_t c_alpha[9][2] =
	{
			{ 0.0,0.0 }, // row 0
			{ 1.0,0.0 }, // row 1
			{ 0.0,1.0 }, // row 2
			{ -1.0,0.0 }, // row 3
			{ 0.0,-1.0 }, // row 4
			{ 1.0,1.0 }, // row 5
			{ -1.0,1.0 }, // row 6
			{ -1.0,-1.0 }, // row 7
			{ 1.0,-1.0 } // row 8
	};

	solver.initialize(src,d,v);

	double beg,end,ct,MLUPs;
	timing(beg,ct);

	for(uint_t t = 0; t <= timesteps; ++t){
		if(vtk_step !=0){
			if((t%vtk_step==0)){
				std::ostringstream os;
				os << "ldc"<<t<<".vtk";
				std::ofstream file( os.str().c_str() );
				os.seekp(0);
				file<<"# vtk DataFile Version 4.0"<<"\n"
						<<"SiwiRVisFile"<<"\n"
						<<"ASCII\n"
						<<"DATASET STRUCTURED_POINTS\n"
						<<"DIMENSIONS "<<x-2<<" "<<y-2<<" 1\n"
						<<"ORIGIN 0 0 0\n"
						<<"SPACING 1 1 1\n"
						<<"POINT_DATA "<<(x-2)*(y-2)<<"\n\n"
						<<"SCALARS density double 1\n"
						<<"LOOKUP_TABLE default\n";
				for(uint_t j = 1; j < y-1 ; ++j){
					for(uint_t i = 1; i < x-1; ++i){
						file<<d(i,j,0)<<"\n";
					}
				}
				file<<"\n"
						<<"VECTORS velocity double\n";
				for(uint_t j = 1; j < y-1 ; ++j){
					for(uint_t i = 1; i < x-1; ++i){
						file<<v(i,j,0)<<" "<<v(i,j,1)<<" "<<"0"<<"\n";
					}
				}
			}
		}
		//Streaming step
		for(uint_t j = 1; j < y-1 ; ++j){
			for(uint_t i = 1; i < x-1; ++i){
				//SW corner point
				if(i==1 && j==1){
					dst(i,j,0) = src(i,j,0);
					dst(i,j,1) = src(i,j,3);	//bounce back
					dst(i,j,2) = src(i,j,4);	//bounce back
					dst(i,j,3) = src(i+1,j,3);
					dst(i,j,4) = src(i,j+1,4);
					dst(i,j,5) = src(i,j,7);	//bounce back
					dst(i,j,6) = src(i,j,8);	//bounce back
					dst(i,j,7) = src(i+1,j+1,7);
					dst(i,j,8) = src(i,j,6);	//bounce back
				}
				//South boundary
				else if(j==1 && (i>1 && i<x-2)){
					dst(i,j,0) = src(i,j,0);
					dst(i,j,1) = src(i-1,j,1);
					dst(i,j,2) = src(i,j,4);	//bounce back
					dst(i,j,3) = src(i+1,j,3);
					dst(i,j,4) = src(i,j+1,4);
					dst(i,j,5) = src(i,j,7);	//bounce back
					dst(i,j,6) = src(i,j,8);	//bounce back
					dst(i,j,7) = src(i+1,j+1,7);
					dst(i,j,8) = src(i-1,j+1,8);
				}
				//SE corner
				else if(j==1 && i==x-2){
					dst(i,j,0) = src(i,j,0);
					dst(i,j,1) = src(i-1,j,1);
					dst(i,j,2) = src(i,j,4);	//bounce back
					dst(i,j,3) = src(i,j,1);	//bounce back
					dst(i,j,4) = src(i,j+1,4);
					dst(i,j,5) = src(i,j,7);	//bounce back
					dst(i,j,6) = src(i,j,8);	//bounce back
					dst(i,j,7) = src(i,j,5);	//bounce back
					dst(i,j,8) = src(i-1,j+1,8);
				}
				//East
				else if(i==x-2 && (j>1 && j<y-2)){
					dst(i,j,0) = src(i,j,0);
					dst(i,j,1) = src(i-1,j,1);
					dst(i,j,2) = src(i,j-1,2);
					dst(i,j,3) = src(i,j,1);	//bounce back
					dst(i,j,4) = src(i,j+1,4);
					dst(i,j,5) = src(i-1,j-1,5);
					dst(i,j,6) = src(i,j,8);	//bounce back
					dst(i,j,7) = src(i,j,5);	//bounce back
					dst(i,j,8) = src(i-1,j+1,8);
				}
				//NE corner
				else if(i==x-2 && j==y-2){
					dst(i,j,0) = src(i,j,0);
					dst(i,j,1) = src(i-1,j,1);
					dst(i,j,2) = src(i,j-1,2);
					dst(i,j,3) = src(i,j,1);	//bounce back
					dst(i,j,4) = src(i,j,2) - (2.0 * t_alpha[2] * 3.0 * (c_alpha[2][0]*0.08));
					dst(i,j,5) = src(i-1,j-1,5);
					dst(i,j,6) = src(i,j,8);	//bounce back
					dst(i,j,7) = src(i,j,5) - (2.0 * t_alpha[5] * 3.0 * (c_alpha[5][0]*0.08));
					dst(i,j,8) = src(i,j,6) - (2.0 * t_alpha[6] * 3.0 * (c_alpha[6][0]*0.08));
				}
				//North boundary
				else if(j==y-2 &&(i>1 && i <x-2)){
					dst(i,j,0) = src(i,j,0);
					dst(i,j,1) = src(i-1,j,1);
					dst(i,j,2) = src(i,j-1,2);
					dst(i,j,3) = src(i+1,j,3);
					dst(i,j,4) = src(i,j,2) - (2.0 * t_alpha[2] * 3.0 * (c_alpha[2][0]*0.08));
					dst(i,j,5) = src(i-1,j-1,5);
					dst(i,j,6) = src(i+1,j-1,6);
					dst(i,j,7) = src(i,j,5) - (2.0 * t_alpha[5] * 3.0 * (c_alpha[5][0]*0.08));
					dst(i,j,8) = src(i,j,6) - (2.0 * t_alpha[6] * 3.0 * (c_alpha[6][0]*0.08));
				}
				//NW corner
				else if(i==1 && j==y-2){
					dst(i,j,0) = src(i,j,0);
					dst(i,j,1) = src(i,j,3);	//bounce back
					dst(i,j,2) = src(i,j-1,2);
					dst(i,j,3) = src(i+1,j,3);
					dst(i,j,4) = src(i,j,2) -( 2.0 * t_alpha[2] * 3.0 * (c_alpha[2][0]*0.08));
					dst(i,j,5) = src(i,j,7);	//bounce back
					dst(i,j,6) = src(i+1,j-1,6);
					dst(i,j,7) = src(i,j,5) - (2.0 * t_alpha[5] * 3.0 * (c_alpha[5][0]*0.08));
					dst(i,j,8) = src(i,j,6) - (2.0 * t_alpha[6] * 3.0 * (c_alpha[6][0]*0.08));
				}
				//West boundary
				else if(i==1 &&(j>1 && j<y-2)){
					dst(i,j,0) = src(i,j,0);
					dst(i,j,1) = src(i,j,3);	//bounce back
					dst(i,j,2) = src(i,j-1,2);
					dst(i,j,3) = src(i+1,j,3);
					dst(i,j,4) = src(i,j+1,4);
					dst(i,j,5) = src(i,j,7);	//bounce back
					dst(i,j,6) = src(i+1,j-1,6);
					dst(i,j,7) = src(i+1,j+1,7);
					dst(i,j,8) = src(i,j,6);	//bounce back
				}
				//Normal cells
				else {
					dst(i,j,0) = src(i,j,0);
					dst(i,j,1) = src(i-1,j,1);
					dst(i,j,2) = src(i,j-1,2);
					dst(i,j,3) = src(i+1,j,3);
					dst(i,j,4) = src(i,j+1,4);
					dst(i,j,5) = src(i-1,j-1,5);
					dst(i,j,6) = src(i+1,j-1,6);
					dst(i,j,7) = src(i+1,j+1,7);
					dst(i,j,8) = src(i-1,j+1,8);
				}

				solver.calc_density(dst,d,i,j);
				solver.calc_velocity(dst,v,d,i,j);
				//std::cout<<"Density before collision = "<<d(i,j,0)<<std::endl;
				//std::cout<<"Velocity before collision = "<<v(i,j,0)<<std::endl;
				solver.collide(dst,omega,v,d,i,j);
				//std::cout<<"Density after collision = "<<d(i,j,0)<<std::endl;
				//std::cout<<"Velocity after collision = "<<v(i,j,0)<<std::endl;
			}
		}

		
		src.swap(dst);
	}
	timing(end,ct);
	std::cout<<"Time taken = "<<(end-beg)<<" sec"<<std::endl;
	MLUPs = (x-2)*(y-2)*timesteps*1e-6/(end-beg);
	std::cout<<"Performance = "<<MLUPs<<" MLUP/s"<<std::endl;
	std::cout<<"Exit success"<<std::endl;
	return 0;
}
