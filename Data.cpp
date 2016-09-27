#include "Data.h"
#include <iostream>
#include <fstream>
#include <cmath>


#include <assert.h>
#include <CCfits/CCfits>

using namespace std;
using namespace CCfits;

Data Data::instance;

Data::Data()
{

}


void Data::load(const char* metadata_file, const char* image_file, const char* sigma_file )
{

     std::auto_ptr<FITS> pInfile(new FITS(image_file,Read,true));
     PHDU& image_=pInfile->pHDU(); 

     std::valarray<double>  contents;      
     std::vector<String> name;

     name.push_back("MAGZP");
//     name.push_back("M01_MODMAG");

//     // read all user-specifed, coordinate, and checksum keys in the image
     image_.readAllKeys();
     image_.read(contents);

//     // this doesn't print the data, just header info.
      std::cout << image_ << std::endl;

//     // create an vector of vectors
     ni = image_.axis(0);
     nj = image_.axis(1);

      image_.keyWord(name[0]).value(magzp);
  //   image_.keyWord(name[1]).value(magzero);  
//      
      x_min = 0;
      y_min = 0;
      y_max = nj;
      x_max = ni;
     
    
      int counter =0; 
      image.assign(ni, vector<double>(nj));
      for (long j = 0;j< nj ; j++)
      {
           for (long i = 0;i < ni ; i++)
           {     
	         image[i][j] = double(contents[counter]);                         // pixels (x,y)
                 counter++;
            }
     }

     max_mag = contents.max();
     min_mag = contents.min();
     std::cout << "Max flux"  << contents.max() << std::endl;
     std::cout << "Min flux"  << contents.min() << std::endl;
     std::cout << "magzp"     << magzp << std::endl;
 
//	fstream fin(metadata_file, ios::in);
//	if(!fin)
//		cerr<<"# ERROR: couldn't open file "<<metadata_file<<"."<<endl;
//	fin>>ni>>nj;
//	fin>>x_min>>x_max>>y_min>>y_max;
//	fin.close();

	// Make sure maximum > minimum
	if(x_max <= x_min || y_max <= y_min)
		cerr<<"# ERROR: strange input in "<<metadata_file<<"."<<endl;

	// Compute pixel widths
	dx = (x_max - x_min)/nj;
	dy = (y_max - y_min)/ni;

	// Check that pixels are square
	if(abs(log(dx/dy)) >= 1E-3)
		cerr<<"# ERROR: pixels aren't square."<<endl;

	compute_ray_grid();





     std::auto_ptr<FITS> pInfile2(new FITS(sigma_file,Read,true));
     PHDU& sigma_=pInfile2->pHDU(); 

     std::valarray<double>  contents_sigma;      
     std::vector<String> name_sigma;

//     // read all user-specifed, coordinate, and checksum keys in the image
     sigma_.readAllKeys();
     sigma_.read(contents_sigma);

//     // this doesn't print the data, just header info.
      std::cout << sigma_ << std::endl;

//     // create an vector of vectors
     ni = sigma_.axis(0);
     nj = sigma_.axis(1);

      x_min = 0;
      y_min = 0;
      y_max = nj;
      x_max = ni;
     
    
      counter =0; 
      sig.assign(ni, vector<double>(nj));
      for (long j = 0;j< nj ; j++)
      {
           for (long i = 0;i < ni ; i++)
           {     
	         sig[i][j] = double(contents_sigma[counter]);                         // pixels (x,y)
                 counter++;
            }
     }
	
}

void Data::compute_ray_grid()
{
	// Make vectors of the correct size
	x_rays.assign(ni, vector<double>(nj));
	y_rays.assign(ni, vector<double>(nj));

	// Distance between adjacent rays
	double L = dx;

	for(size_t i=0; i<x_rays.size(); i++)
	{
		for(size_t j=0; j<x_rays[i].size(); j++)
		{
			x_rays[i][j] = x_min + (j + 0.5)*L;
			y_rays[i][j] = y_max - (i + 0.5)*L;
		}
	}
}

