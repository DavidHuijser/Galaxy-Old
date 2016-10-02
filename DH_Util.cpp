
// gcc -Wall -c DH_Util.cpp

// #include "RJ_Util5.h"
#include <sstream>
#include <string> 
#include <fstream>                 // includes
#include <iostream>
#include <CCfits/CCfits>
#include <math.h>       /* tgamma */
#define MATHLIB_STANDALONE 1



#include <Rmath.h>
#define PI 3.14159265352982





/*
   Student t Distibution 
   90     inline double dt(double x, double n, int lg)                        { return ::Rf_dt(x, n, lg); }
   91     inline double pt(double x, double n, int lt, int lg)                { return ::Rf_pt(x, n, lt, lg); }
   92     inline double qt(double p, double n, int lt, int lg)                { return ::Rf_qt(p, n, lt, lg); }
   93     inline double rt(double n)     
 */

// [[Rcpp::export]] 
double logsubstCpp(double x, double y)
{ if( x >= y){ return x + log(1.0 - exp(y - x)); }else{ return NAN; } }




// Probability according to Uniform probability
double DH_prandom(double val, double lower_limit, double upper_limit, bool logp)
{
  double in_range;
  if (logp == FALSE)
  {
     in_range = 1.0;
     if  ((val < lower_limit) || (val > upper_limit))
     {
       in_range = 0.0;
     }
  }
  else
  {
     in_range = 0.0;
     if  ((val < lower_limit) || (val > upper_limit))
     {
       in_range = -1e300;
     }
  }
  return in_range;
}


// Draw from uniform probability
double DH_drandom(double lower_limit, double upper_limit)
{
  double val  = (upper_limit - lower_limit)*(rand()/double(RAND_MAX)) + lower_limit; 
  return val;
}


//  Normal Probability 
double DH_pnorm(double x, double mean, double sd)    // fraction = numerator / denominator
{
    double numerator = exp(-0.5*pow((x-mean)/sd,2) );
    double denominator = pow(2.0*M_PI,0.5)*sd;
    return (numerator/denominator);    
}


// Draw from Normal
double DH_dnorm(void)
{
  double u1  = DH_drandom(0,1); 
  double u2  = DH_drandom(0,1); 
  double z = sqrt(-2.0*log(u1))*cos(2.0*PI*u2);
  return z;
}

//// [[Rcpp::export]]
//double DH_dt(double gamma, bool logp)
//{
//    double x;
//    if (logp == FALSE)
//    {
//     x = rt(gamma);
//    }
//    else 
//    {
//      x = rt(gamma);  
//    }
//    return x;
//}

//// [[Rcpp::export]]
//double DH_pt(double x, double gamma, bool logp)
//{
//    double p;
//    if (logp == FALSE)
//    {
//     p =  pt(x,gamma, 0, 0); 
//    }
//    else 
//    {
//     p = pt(x, gamma, 0, 1);
//    }
// //  p =  tgamma(0.5*(gamma+1.0))*pow(1.0+x*x/gamma,-0.5*(gamma+1.0) ) / ( pow(gamma*M_PI   , 0.5)*tgamma(0.5*gamma) );
////  p = pt(x, gamma);
//  return p;
//}


double DH_pTruncNormal(double x, double mean, double sigma, double a, double b, bool logp)
{
    double p = 0;
    double log_p = -1e300;
    double results = 0;
    double denominator = 0;
    if ( ((x > a)  && ( x < b)) and (logp == FALSE))
    {
          results = dnorm4((x-mean)/sigma, 0.0, 1.0, 0)/sigma;
          denominator = pnorm5((b-mean)/sigma,0, 1,1,0)  - pnorm5((a-mean)/sigma,0, 1,1,0);
          p = results/denominator;
    } //end-if

    // This not seem to work properly
    if (((x > a)  && ( x < b)) and (logp == TRUE))
    {  
         results = dnorm4((x-mean)/sigma, 0.0, 1.0, TRUE) - log(sigma);
 //        ------- method 1 - works well
 //        denominator = exp(pnorm5((b-mean)/sigma,0, 1,1,TRUE))  - exp(pnorm5((a-mean)/sigma,0, 1,1,TRUE));
 //        log_p = results - log(denominator);
 //       ------- method 2 - works well  too
         denominator = logsubstCpp(pnorm5((b-mean)/sigma,0, 1,1,1),  pnorm5((a-mean)/sigma,0, 1,1,1));     
         log_p = results - denominator;
    } //end-if
    if (logp == TRUE )
    {        
      return(log_p);
    }
    else
    {        
      return(p);
    }
}




double DH_dTruncNormal(double mean, double sd, double a, double b)
{
    int len = 300;
//    int burn = 500; 
    std::vector<double>  chain;
//    std::vector<double>  new_chain;
    chain.assign(len,0);
  //  new_chain.assign(500,0);
    chain[1] =  mean;
    double proposal, probab;  
    double u;
    for (int i=0; i < len; i++)
    {
      proposal =  3.0*DH_dnorm() + chain[i];
     probab = DH_pTruncNormal(proposal,mean, sd,a,b,0) / DH_pTruncNormal(chain[i],mean, sd,a,b,0);
     u = DH_drandom(0,1); 
     if (u < probab){
       chain[i+1] = proposal;
     }else{
       chain[i+1] = chain[i];
     }
   }
   return(chain[len-1]);
 }


// Store in file 
/*
void writefile(string filename, const std::vector< std::vector<double> > &image)
{
      string flname=filename;
      std::ofstream myout; // This is the object that sends data to the file. It is used like cout, but now writes in file. 
      myout.open(flname.c_str());  // ofstream means Output File Stream.
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


double mod(double y, double x)
{
	if(x <= 0)
		std::cerr<<"Warning in mod(double, double) (Utils.cpp)"<<std::endl;
	return (y/x - floor(y/x))*x;
}

// wrap around:  if outside range, it wraps aroud and enters from the other side
void wrap(double& x, double min, double max)
{
	x = mod(x - min, max - min) + min;
}

// modeluo function
int mod(int y, int x)
{
	if(x <= 0)
		std::cerr<<"Warning in mod(int, int) (Utils.cpp)"<<std::endl;
	if(y >= 0)
		return y - (y/x)*x;
	else
		return (x-1) - mod(-y-1, x);
}

*/
std::string IntToStr(int n) 
{
    std::stringstream result;
    result << n;
    return result.str();
}


// FUNCTION pcauchy calc the probability density of the given point.
double DH_pcauchy(double x, double location, double scale )
{
  double prob =  scale/ (M_PI*(scale*scale  + ( x  - location )*( x  - location ))  );
  return prob;
}
 
// FUNCTION dcauchy (draw from a cauchy distribution) 
  double DH_dcauchy(double location, double scale){
    double u  = DH_drandom(0,1); 
    double x;   
    x = location + scale*std::tan(M_PI*( u - 0.5 ) );
    return x;
}
