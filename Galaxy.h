#ifndef DNest4_Galaxy
#define DNest4_Galaxy

#include "DNest4/code/DNest4.h"
#include <valarray>
#include <ostream>
#include <string> 

#include <sstream>
#include <fstream>                 // includes
#include <iostream>                // includes:  cerr

#define ARMA_DONT_USE_CXX11 
#include <armadillo>


using namespace std;
class Galaxy
{
	private:
	       std::valarray<double> params;

		// The model image
		std::vector< std::vector<long double> > image;
		std::vector< std::vector<long double> > psf_image;
		void calculate_image();
         	double sigma;
                double magzp;
        

                // FFT of the PSF
                arma::cx_mat fft_of_psf;
                bool fft_ready;
                int Ni, Nj;

	public:
                int global_index;
		// Constructor only gives size of params
		Galaxy();                 
		// Generate the point from the prior
		void from_prior(DNest4::RNG& rng);

		// Metropolis-Hastings proposals
		double perturb(DNest4::RNG& rng);

		// Likelihood function
		double log_likelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;
 
                // Set Limits
                void set_limit(); 

                // Limits
                std::vector<double> upper_limit;
                std::vector<double> lower_limit;
 
                void writefile(const std::vector< std::vector<long double> >  &image);
 
	        std::string IntToStr(int n);
           
                // PSF
                void load_psf();
                void calculate_fft(int Ni, int Nj);
                void blur_image2(std::vector< std::vector<long double> >  &img);

};

#endif
