#include "Galaxy.h"
#include "DNest4/code/DNest4.h"
#include "Data.h"
#include "DH_Util.h"
#include <string> 

#include <sstream>
#include <fstream>                 // includes
#include <iostream>


#define ARMA_DONT_USE_CXX11 
#include <armadillo>

/* Mandatory */
 #include <profit/profit.h>
//#include <"/home/dhui890/Documents/Research/Profit/libprofit/profit/profit.h>

/* Depending on which profiles you will use... */
#include <profit/sersic.h>
#include <profit/ferrer.h>
#include <profit/moffat.h>
#include <profit/sky.h>
#include <profit/psf.h>


// use:   more sample_info.txt 

#include <mutex>          // std::mutex



using namespace std;
using namespace DNest4;
using namespace arma;

Galaxy::Galaxy()
:params(15)
{
    set_limit();
    global_index = 0;
    calculate_fft(Data::get_instance().get_ni(), Data::get_instance().get_nj());

    
}


//double Galaxy::draw_tdistribution(double gamma, double location, double shape)
//{         
//      double a = rng.randn();        // rng.rand(); uniform 
//      double b = rng.rand();         // rng.randn(); gaussian
//      x =   a / sqrt(-ln(b));   
//      x = shape*x + location;
//      return x;
//}


//double Galaxy::prob_tdistribution(double x, double gamma, double location, double shape)
//{
//  double p = (p-location)/shape;

//}
//      		


void Galaxy::set_limit()
{ 
    double delta = 1e-1;
    lower_limit.assign(15,0);
    upper_limit.assign(15,0);
    
//	   xc    = components[k][0];     yc = components[k][1];     
//           I_gal = components[k][2];  R_gal = components[k][3]; n_gal = components[k][4]; q _gal = components[k][5]; theta_gal = components[k][6]; c_gal = components[k][7];
//	   I_disc= components[k][8];  R_disc = components[k][9]; q_disc   = components[k][10]; theta_disc = components[k][11]; c_disc = components[k][12];	 
//	   I_bar = components[k][13]; Rout_bar = components[k][14]; rout_a = components[k][15]; rout_b = components[k][16];
//                                       q_bar   = components[k][17]; theta_bar = components[k][18]; c_bar = components[k][19]
    lower_limit[0] = 0; //x
    lower_limit[1] = 0; //y
    lower_limit[2] = Data::get_instance().get_magzp() - 18; // I
    lower_limit[3] = 1;  // R
    lower_limit[4] = 0.5; // n 
    lower_limit[5] = 0;  // q
    lower_limit[6] = 0;    // theta
    lower_limit[7] = -1.0;  // boxi

    lower_limit[8] = Data::get_instance().get_magzp() - 18; // I
    lower_limit[9] = delta;  // rout
    lower_limit[10] =  delta; // a                      // alpha [ 0, 10]     (Tayuen, Kim paper 2014) 
    lower_limit[11] =  delta;   // b                    // beta [0 ,2] 
    lower_limit[12] = 0;    // q
    lower_limit[13] = 0 ;  // theta  
    lower_limit[14] = -1.0;  // boxi


// here A rat is  the  minor  to  major  axis  ratio  (so  always  a number between 0 and 1, where 0 is an infinitely thin line and 1 is a circle or dis


    upper_limit[0] = 100; //x
    upper_limit[1] = 100; //y
    upper_limit[2] = Data::get_instance().get_magzp(); // I
    upper_limit[3] = 50;  // R
    upper_limit[4] = 10; // n 
    upper_limit[5] = 1;  // q
    upper_limit[6] = 180;    // theta
    upper_limit[7] = 1.0;  // boxi

    upper_limit[8] = Data::get_instance().get_magzp(); // I
    upper_limit[9] = 3;  // rout
    upper_limit[10] = 10; // a
    upper_limit[11] = 2.0 - delta;   // b
    upper_limit[12] = 1;    // q
    upper_limit[13] = 180 ;  // theta  
    upper_limit[14] = 1.0;  // boxi

}






void Galaxy::from_prior(RNG& rng)
{
    //    std::cout << std::endl << " Started drawing from prior " << std::endl;
    //    std::cout << std::endl <<  params.size() << std::endl;  
//      params[0] = x_min + (x_max - x_min)*rng.rand();
      // and so on...

      for(size_t i=0; i<params.size(); ++i)
      {
       if ((i == 2) || (i == 8))
       {
               params[i] = 5*rt(2) + 25;   // shape=2, location=25, scale=5
        }
        else
        {
              params[i] = lower_limit[i] + (upper_limit[i] - lower_limit[i])*rng.rand();
        }
       } 
	calculate_image();
   //     std::cout << " Finished image " << std::endl; 
}


double Galaxy::perturb(RNG& rng)
{      
	int which = rng.rand_int(params.size());
        double log_H = 0.0;
        if ((which == 2) || (which == 8))
        {
                log_H -= dt((params[which]-25)/5, 2,0);
         	params[which] += 5*rng.randh();
                log_H += dt((params[which]-25)/5, 2,0);

         }
         else
         {

         	params[which] += (upper_limit[which] - lower_limit[which])*rng.randh();
	        wrap(params[which], lower_limit[which], upper_limit[which]);                   // wraps?
        }
	calculate_image();
 	return log_H;

}



double Galaxy::log_likelihood() const
{
  //      std::cout << std::endl << " Started eval loglikelihood " << std::endl;   
	const vector< vector<double> >& data = Data::get_instance().get_image();
	const vector< vector<double> >& sig =  Data::get_instance().get_sigma();

	double logL = 0.;
	double var;
	for(size_t i=0; i<data.size(); i++)
	{
		for(size_t j=0; j<data[i].size(); j++)
		{
			var = sigma*sigma + sig[i][j]*sig[i][j];
			logL += -0.5*log(2.*M_PI*var)
				-0.5*pow(data[i][j] - image[i][j], 2)/var;
		}
	}

	if(std::isnan(logL) || std::isinf(logL))
		logL = -1E300;

        //std::cout << std::endl << " Finished eval loglikelihood sigma:" <<  var << std::endl;    
        //std::cout << data.size() << " " << data[0].size() << std::endl;    
	return logL;

}


string Galaxy::description() const
{
	return string("objects");
        return string("Each column is one of the 20 parameters.");
}



void Galaxy::print(std::ostream& out) const
{
	for(const double& x: params)
		out<<x<<' ';
        for (long i=0; i < long(image.size()); i++)
	  {
	     for (long j=0;  j < long(image[0].size()); j++)
	     {
	       out << image[i][j] <<" "; 
	     } //end-for-j	   	     
	  } //end for -i       
} // 

void Galaxy::calculate_image()
{

       // static std::mutex mtx2;           // mutex for critical section
      //  mtx2.lock();
     //   std::cout << std::endl << " Started making image " << std::endl;    
        global_index++;
	// Get coordinate stuff from data
	//const vector< vector<double> >& x = Data::get_instance().get_x_rays();
//	const vector< vector<double> >& y = Data::get_instance().get_y_rays();
        int ni= Data::get_instance().get_ni();
        int nj= Data::get_instance().get_nj();

	image.assign(Data::get_instance().get_ni(), vector<long double>(Data::get_instance().get_nj(), 0.));
	
        // the model consist of two general coordinates of a bulge, disk and a bar 
        //  xc(0), yc(1)         (general)
        //  sersic: I_gal(2), R_gal(3), n_gal(4), q_gal(5), theta_gal(6), c_gal(7)  (galaxy)
        //  sersic: I_disc(8), R_disc(9), q_disc=(10), theta_disc(11), c_disc(12) (disc, n is fixed )
        //  ferro: I_bar(13),rout_bar(14), rout_a(15),rout_b(16), q_bar=(17), theta_bar(18), c_bar(19) 

//         double xc, yc;
//         double I_gal, R_gal, n_gal, q_gal, theta_gal, c_gal;
//         double I_disc, R_disc, q_disc, theta_disc, c_disc;
//         double I_bar,rout_bar, rout_a,rout_b, q_bar, theta_bar, c_bar; 

//	   xc    = components[k][0];     yc = components[k][1];     
//           I_gal = components[k][2];  R_gal = components[k][3]; n_gal = components[k][4]; q _gal = components[k][5]; theta_gal = components[k][6]; c_gal = components[k][7];
//	   I_disc= components[k][8];  R_disc = components[k][9]; q_disc   = components[k][10]; theta_disc = components[k][11]; c_disc = components[k][12];	 
//	   I_bar = components[k][13]; Rout_bar = components[k][14]; rout_a = components[k][15]; rout_b = components[k][16];
//                                       q_bar   = components[k][17]; theta_bar = components[k][18]; c_bar = components[k][19];



//So, you'd have to go like this:

//// create the profile
//profit::Model *model = new profit::Model();
//model->width = metadata[0];
//model->height = metadata[1];

//// add a sersic profile it and tune it
//profit::Profile *sersic_profile = model->add_profile("sersic");
//profit::SersicProfile *sp = static_cast<profit::SersicProfile*>(sersic_profile);
//sp->xcen = ...

//// add a ferrer profile and tune it
//profit::Profile *ferrer_profile = model->add_profile("ferrer");
//profit::FerrerProfile *fp = static_cast<profit::FerrerProfile*>(ferrer_profile);
//fp->xcen = ...



//    rout: The outer truncation radius.
//    a: The global power-law slope to the profile center
//    b: The strength of truncation as the radius approaches rout.


//// finally:
//model->evaluate();

     profit::Model *model = new profit::Model();
     profit::Profile *sersic_profile = model->add_profile("sersic");

     model->width = ni; 
     model->height = nj;
     model->magzero = Data::get_instance().get_magzp();

     profit::SersicProfile *sp = static_cast<profit::SersicProfile *>(sersic_profile);
     sp->xcen = params[0];
     sp->ycen = params[1];
     sp->mag = params[2];
     sp->re =  params[3];
     sp->nser = params[4];
     sp->axrat = params[5];
     sp->ang = params[6];
     sp->box = params[7];


     profit::Profile *ferrer_profile = model->add_profile("ferrer");
     profit::FerrerProfile *fp = static_cast<profit::FerrerProfile*>(ferrer_profile);
     fp->xcen = params[0];
     fp->ycen = params[1];
     fp->mag = params[8];
     fp->rout = params[9]*params[3];
     fp->a = params[10];
     fp->b = params[11];
     fp->axrat = params[12]*params[5];  // ax_bar < ax_disc
     fp->ang = params[13];
     fp->box = params[14];

     double x, y;
     double half_xbin = model->scale_x/2.;
     double half_ybin = model->scale_y/2.;
    // std::cout << "rout= " << params[9]*params[3] << " a= " << params[10] << " b= " << params[11] << std::endl;

     model->evaluate();

     for (int j=0; j < nj; j++) 
     {
		x += half_xbin;
		y = 0;
		for (int k=0; k < ni; k++) 
                {
		       y += half_ybin;                       
                       image[j][k] = image[j][k]  +   model->image[j + k*nj];
         	       y += half_ybin;
		} // end-for
		x += half_xbin;
      } //end-for
      delete model;

     //  PSF_Object.blur_image2(image);
     
//          profit::Model *bulge = new profit::Model();
//          profit::Profile *sersic_profile = bulge->add_profile("sersic");

//          bulge->width = ni; 
//          bulge->height = ni;

//          profit::SersicProfile *sp_bulge = static_cast<profit::SersicProfile *>(sersic_profile);
//          sp_bulge->xcen = params[0];
//          sp_bulge->ycen = params[1];
//          sp_bulge->mag = params[2];
//          sp_bulge->re = params[3];
//          sp_bulge->nser = params[4];
//          sp_bulge->axrat = params[5];
//          sp_bulge->ang = params[6];
//          sp_bulge->box = params[7];
////          bulge->magzero = metadata[3];
//          double x; double y;
//          double half_xbin = bulge->scale_x/2.;
//          double half_ybin = bulge->scale_y/2.;
//          bulge->evaluate();

//     for (int j=0; j < nj; j++) 
//     {
//		x += half_xbin;
//		y = 0;
//		for (int i=0; i < ni; i++) 
//                {
//		       y += half_ybin;                       
//                       image[j][i] = image[j][i]  +   bulge->image[j + i*ni];
//         	       y += half_ybin;
//		} // end-for
//		x += half_xbin;
//      } //end-for
//      delete bulge;
        blur_image2(image);
   //     mtx2.unlock();
  
       if (global_index % 500 == 0)
        {
             writefile(image); 
        }
   //   std::cout << std::endl << " Ended making image " << std::endl;    

//        mtx.unlock(); 

}




void Galaxy::writefile(const std::vector< std::vector<long double> >  &image)
{

     int sub2 = int(global_index / 500); 
 //   std::stringstream result;
  //  result << sub2;

//     std::string sub = result.str();
//      std::string name="file_" + sub  +".txt";        

      std::string name="file_" +      std::to_string(sub2)  +".txt";        
      std::ofstream myout; // This is the object that sends data to the file. It is used like cout, but now writes in file. 
      myout.open(name.c_str());  // ofstream means Output File Stream.
    //  std::cout << "This is an " << image.size()  << " by " << image[0].size() << " image" << std::endl;
	  for (long i=0; i < long(image.size()); i++)
	  {
	     for (long j=0;  j < long(image[0].size()); j++)
	     {
	       myout << image[i][j] <<" "; 
	     } //end-for-j	   
	     myout << std::endl;      // preferable for c++  (but gives readin problems in read-in in IDL)     	     
	  }; //end for -i  
// create a new starting array, and filling it with random values
myout.close();  
}





void Galaxy::calculate_fft(int Ni, int Nj)
{
   const vector< vector<double> >& psf_image =  Data::get_instance().get_psf();

// Make the psf the same size as the image
mat psf(Ni, Nj);
psf.zeros();

int ni = psf_image.size();
int nj = psf_image[0].size();

int m, n;
for(int i=0; i<ni; i++)
{
m = mod(i - ni/2, Ni);
for(int j=0; j<nj; j++)
{
n =  mod(j - nj/2, Nj);
psf(m, n) = psf_image[i][j];
}
}

fft_of_psf = fft2(psf);
fft_ready = true;
}

void Galaxy::blur_image2(std::vector< std::vector<long double> >  &img) 
{
 
if(!fft_ready)
cerr<<"# Blurring failed."<<endl;

// Copy the image into an Armadillo matrix
mat A(img.size(), img[0].size());
for(size_t i=0; i<img.size(); i++)	
for(size_t j=0; j<img[0].size(); j++)
A(i, j) = img[i][j];

// Do the fft of it
cx_mat B = fft2(A);

// Multiply the two ffts
for(size_t i=0; i<img.size(); i++)
for(size_t j=0; j<img[0].size(); j++)
B(i, j) *= fft_of_psf(i, j);

// Do the inverse fft
B = ifft2(B);

// Put back in img
for(size_t i=0; i<img.size(); i++)
for(size_t j=0; j<img[0].size(); j++)
img[i][j] = real(B(i, j));
}


