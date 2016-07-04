
#ifndef GRID_H_
#define GRID_H_
#include <iostream>
#include<assert.h>
#include<fstream>
#include<sstream>
#include<string>


namespace lbm{
enum lat_dir{
	C = 0,
	E = 1,
	N = 2,
	W = 3,
	S = 4,
	NE = 5,
	NW = 6,
	SW = 7,
	SE = 8
};

typedef unsigned int uint_t;
typedef double real_t;

template< typename Type, uint_t Cellsize >
class Grid
{
public:
	Grid():xsize_(0), ysize_(0), data_(0){}
	Grid( uint_t xsize, uint_t ysize ): xsize_(xsize), ysize_(ysize),
			data_(new Type[Cellsize * xsize * ysize]){}
	~Grid(){
		delete []data_;
		data_ = 0;
	}
	Type& operator()( uint_t x, uint_t y, uint_t f ){
		assert(x < xsize_ && y < ysize_ && f < Cellsize);
		return data_[y*xsize_*Cellsize + x*Cellsize + f];
	}
	Type operator()( uint_t x, uint_t y, uint_t f ) const{
		assert(x < xsize_ && y < ysize_ && f < Cellsize);
		return data_[y*xsize_*Cellsize + x*Cellsize + f];
	}
	void swap(Grid& grid){
		ysize_ = grid.ysize_ ;
		xsize_ = grid.xsize_ ;
		std::swap(data_,grid.data_ );
	}
	uint_t get_xsize(){
		return xsize_;
	}
	uint_t get_ysize(){
		return ysize_;
	}

private:
	uint_t xsize_;		// Number of nodes in x-dimension
	uint_t ysize_;		// Number of nodes in y-dimension
	Type* data_;		// Linearized, 1-dimensional representation
						// of the 2D data grid


};

typedef Grid<real_t, 9> PDF_Field;
typedef Grid<real_t, 2> VelocityField;
typedef Grid<real_t, 1> DensityField;
}

#endif /* GRID_H_ */
