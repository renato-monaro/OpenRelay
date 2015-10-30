/*
# Copyright (C) 2008-2015 Renato Machado Monaro
# This file is part of OpenRelay.
#
# OpenRelay is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenRelay is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
*/
#ifndef FILTER_H
#define FILTER_H

#include <iostream>
#include <fstream>
#include <math.h>
#include "types.h"

using namespace std;

#define BANDPASS 1
#define DIFFERENTIATOR 2
#define HILBERT 3
#define NEGATIVE_F 0
#define POSITIVE_F 1
#define GRIDDENSITY 16
#define MAXITERATIONS 40



/*
class Filter{
	public:
		virtual double Apply()=0;
	};
class FIR: public Filter{
	public:
		double Apply(){return 0.0;};
	protected:
		vector<double> DemCoeff;
};

class IIR: public Filter{
	public:
		double Apply(){return 0.0;};
	protected:
		vector<double> NumCoeff;
		vector<double> DemCoeff;
};

class Butterworth: public IIR{
	public: Butterworth(unsigned Order, double Cut_off_Frequency);
};

class Remez: public FIR{
	public: 
		Remez(unsigned NumTaps,vector<double> Bands, vector<double> Desired, vector<double> Weights, unsigned Type);
};
*/


template <class Numeric> class iir_1st
{
    protected:   
    	double gain;                    
    	Numeric out;
        Numeric previous_out;
		Numeric previous_in;
        
    public:
        iir_1st(double A=0) : gain(A) {
        	previous_in = previous_out = out = 0 ; }
		void set_coeff(double A) { gain=A;}
		//! Constructor reading coefficient from a file.
		iir_1st(const char* file)
		  {
			FILE *iirf = fopen(file,"r"); 
			fscanf(iirf,"%lf",&gain);
			fclose(iirf);
			previous_in = previous_out = out = 0;
		  }             
		//! Print out coefficients
		void print() {
			cout << "IIR Coefficient gain = " << gain << "\n";
		}
		//! Input new sample and calculate output
		Numeric clock(Numeric input) {
		  // Shift previous outputs and calculate new output */
		  //	out = gain*previous_out + (1-gain)*input;
		  out = gain*previous_out + (input+previous_in);
		  previous_out = out;
		  previous_in = input;
		  return(out);
		}
		//! Reset
		void reset() {
		  previous_in = previous_out = out = 0;
		}
};
template <class Numeric> class iir_2nd
{
    protected:   
    	Numeric b0,b1,b2;                    
    	Numeric a1,a2;
        Numeric in[3];
        Numeric out[3]; 
        
    public:
        iir_2nd(Numeric B0, Numeric B1, Numeric B2, Numeric A1, Numeric A2) : b0(B0), b1(B1), b2(B2), a1(A1), a2(A2) {
        	in[0] = in[1] = in[2] = out[2] = out[1] = out[0] = 0 ; }
        iir_2nd(Numeric A1=0, Numeric A2=0) : b0(1), b1(2), b2(1), a1(A1), a2(A2) {
        	in[0] = in[1] = in[2] = out[2] = out[1] = out[0] = 0 ; } 
        void reset() {
        	in[0] = in[1] = in[2] = out[2] = out[1] = out[0] = 0 ; } 
		void set_a(Numeric A1, Numeric A2) { a1=A1; a2=A2;}
		void set_b(Numeric A1, Numeric A2) { b1=A1; b2=A2;}
		void set_coeff(Numeric A1, Numeric A2) { a1=A1; a2=A2;}

//! Constructor reading coefficients from a file.
iir_2nd(const char* file)
{
	FILE *iirf = fopen(file,"r"); 
	fscanf(iirf,"%lf %lf %lf %lf %lf",&b0,&b1,&b2,&a1,&a2);
	fclose(iirf);
   	in[0] = in[1] = in[2] = out[2] = out[1] = out[0] = 0;
}             
//! Print out coefficients
void print() {
    printf("IIR Coefficients B0 = %lf, B1 =  %lf, B2 = %lf",b0,b1,b2);
    printf(" A0 = 1, A1 = %lf, A2 = %lf\n",a1,a2);
}
//! Input new sample and calculate output
Numeric clock(Numeric input) {
	// Shift inputs by one time sample and place new sample into array
	in[0] = in[1];
	in[1] = in[2];
	in[2] = input;
	// Shift previous outputs and calculate new output */
	out[0] = out[1];
	out[1] = out[2];
	out[2] = b0*in[2] + b1*in[1] + b2*in[0] - a1*out[1] - a2*out[0];
	return(out[2]);
}
};            


template <class Numeric> class butterworth 
{
	private:
	long order;
	long odd;
	long n2;
	complex<double>* roots;
	iir_2nd< Numeric>* iir;
	iir_1st< Numeric>* iir_1;
	double gain;
	double* coeff;

	public:
	//! Constructor, fcd = cut-off (1=sampling rate)
	//! ord = Filter order
	//! amax = attenuation at cut-off
	butterworth(double fcd, long ord=1, double amax=3.0) {
		// amax - attenuation at cutoff
		gain = 1;
		order = ord;
		double epi = pow( (pow(10.0,(amax/10.)) - 1.0) ,(1./(2.0*ord)));
		// fcd - desired cutoff frequency (normalized)
		double wca = 2.0*tan(M_PI*fcd)/epi;
		// wca - pre-warped angular frequency
	    n2 = (order+1)/2;
		odd = (order%2);
		roots = new complex<double>[n2];
		if (odd) iir_1 = new iir_1st< Numeric >;
		get_roots(wca,order,n2);
		bilinear(n2);
		iir = new iir_2nd< Numeric >[order/2];
		get_coeff(n2);
	}
	//! Destructor
	~butterworth() {
	  delete [] roots;
	  if (odd) delete iir_1;
	  delete [] iir;
	}
	//! Reset history
	void reset() {
	  for (int j=odd;j<n2;j++) { 
		iir[j-odd].reset();
		if (odd) iir_1->reset();
	  }
	} 
	//! print coefficients
	void print() {
	  for (int j=odd;j<n2;j++) { 
		iir[j-odd].print();
	  }
	} 
	//! Clock in sample and get output.
	Numeric clock(Numeric in) {
		Numeric tmp = in;
		for (int i=odd;i<n2;i++) {
			tmp = iir[i-odd].clock(tmp);
		}
		if (odd) tmp = iir_1->clock(tmp);
		return(gain*tmp);
	}
	private:
	//! Calculate roots
	void get_roots(double wp, long n, long n2) {
       long l = 0;
       if (n%2 == 0) l = 1;                                                 
	   double arg;
	   for (int j=0;j<n2;j++) {
		 arg = -0.5*M_PI*l/float(n);
		 roots[j] = wp*expj(arg);
		 l += 2;
	  }
	}
	//! Do bilinear transformation
	void bilinear(long n2) {
      double td;
	  int j;
	  const double a = 2.;
      if (odd) {
		// For s-domain zero (or pole) at infinity we
		// get z-domain zero (or pole) at z = -1
		double tmp = (roots[0].real()-a)/(a+roots[0].real());
		iir_1->set_coeff(-tmp);
		gain *= 0.5*(1+tmp);
	  }
      for (j=odd;j<n2;j++) {                       
		td = a*a - 2*a*roots[j].real() + magsq(roots[j]);
		roots[j] = complex<double>((a*a - magsq(roots[j])),
					2.0*a*roots[j].imag())/td;
	  }
	}
	//! Get 2nd order IIR coefficients
	void get_coeff(long n2) {
        for (int j=odd;j<n2;j++) { 
		  iir[j-odd].set_coeff(-2*roots[j].real()/magsq(roots[j]),1.0/magsq(roots[j]));
		  gain *= (magsq(roots[j]) - 2*roots[j].real() + 1.0)/(4*magsq(roots[j])); 
	  }
	}
};                               


template <class Numeric> class chebyshev 
{
	private:
	long order;
	long odd;
	long n2;
	complex<double>* roots;
	iir_2nd< Numeric>* iir;
	iir_1st< Numeric>* iir_1;
	double gain;
	double* coeff;
	double epi;

	public:
	//! Constructor, fcd = cut-off (1=sampling rate)
	//! ord = Filter order
	//! ripple = passband ripple in dB
	chebyshev(double fcd, long ord=1, double ripple=3.0) {
		gain = 1;
		order = ord;
		epi = pow(10.0,(ripple/10.)) - 1.0;
		epi = pow(epi,(1./(1.0*ord)));
		double wca = tan(0.5*M_PI*fcd);
		//! wca - pre-warped angular frequency
	    n2 = (order+1)/2;
		odd = (order%2);
		roots = new complex<double>[n2];
		if (odd) iir_1 = new iir_1st< Numeric >;
		get_roots(wca,order,n2);
		bilinear(n2);
		iir = new iir_2nd< Numeric >[order/2];
		get_coeff(n2);
	}
	//! Destructor
	~chebyshev() {
	  delete [] roots;
	  if (odd) delete iir_1;
	  delete [] iir;
	}
	//! Reset history
	void reset() {
	  for (int j=odd;j<n2;j++) { 
		iir[j-odd].reset();
		if (odd) iir_1->reset();
	  }
	} 
	//! print coefficients
	void print() {
		int j;
		for (j=odd;j<n2;j++) { 
			iir[j-odd].print();
		}
		if (odd) iir_1->print();
	} 
	//! Clock in sample and get output.
	Numeric clock(Numeric in) {
		Numeric tmp = in;
		for (int i=odd;i<n2;i++) {
			tmp = iir[i-odd].clock(tmp);
		}
		if (odd) tmp = iir_1->clock(tmp);
		return(gain*tmp);
	}
	private:
	//! Calculate roots (chebyshev)
	void get_roots(double wp, long n, long n2) {
       long l = 0;
       if (n%2 == 0) l = 1;                                                 
	   double arg;
	   double x = 1/epi;
	   double asinh = log(x + sqrt(1.0+x*x));
	   double v0 = asinh/(double(n));
	   double sm = sinh(v0);
	   double cm = cosh(v0);
	   for (int j=0;j<n2;j++) {
		 arg = -0.5*M_PI*l/float(n);
		 roots[j] = wp*complex<double>(-sm*cos(arg),cm*sin(arg));
		 l += 2;
	  }
	}
	//! Do bilinear transformation
	void bilinear(long n2) {
      double td;
	  int j;
	  const double a = 1.;
      if (odd) {
		// For s-domain zero (or pole) at infinity we
		// get z-domain zero (or pole) at z = -1
		double tmp = (roots[0].real()-a)/(a+roots[0].real());
		iir_1->set_coeff(-tmp);
		gain *= 0.5*(1+tmp);
	  }
      for (j=odd;j<n2;j++) {                       
		td = a*a - 2*a*roots[j].real() + magsq(roots[j]);
		roots[j] = complex<double>((a*a - magsq(roots[j])),
					2.0*a*roots[j].imag())/td;
	  }
	}
	//! Get 2nd order IIR coefficients
	void get_coeff(long n2) {
        for (int j=odd;j<n2;j++) { 
		  iir[j-odd].set_coeff(-2*roots[j].real(),magsq(roots[j]));
		  gain *= (magsq(roots[j]) - 2*roots[j].real() + 1.0)/(4*magsq(roots[j])); 
	  }
	}
};


class remez_fir {

 public:
	remez_fir() {}  

 private:
	static void createDenseGrid(int r, int numtaps, int numband, double bands[],
                       double des[], double weight[], int gridSize,
                       double grid[], double d[], double w[],
                       int symmetry);
	static void initialGuess(int r, int ext[], int gridSize);
	static void calcParms(int r, int ext[], double grid[], double d[], double w[],
				 double ad[], double x[], double y[]);
	static double computeA(double freq, int r, double ad[], double x[], double y[]);
	static void calcError(int r, double ad[], double x[], double y[],
				 int gridSize, double grid[],
				 double d[], double w[], double e[]);
	static void search(int r, int ext[], int gridSize, double e[]);
	static double* freqSample(double A[],int numtaps, int symm);
	static bool isDone(int r, int ext[], double e[]);
 public:
	static double* remez(int n, int numband,
						 double bands[], double des[], double weight[], int type);
};


template <class Numeric> class fir
{
 public: 
  long num_taps;
  Numeric* coeff;
  //      protected:
  Numeric* z; 
  Numeric output;
  
 public: 
  //! Set tap weights
  void settap(long i, Numeric tap) { coeff[i] = tap; }  
  //! Reset
  void reset() { 	
	for (int i=0;i<num_taps;i++) z[i] = 0;  
	output = 0;
  }
  //! Get sum of coefficients
  Numeric coeff_sum() { 
	  int i;
	  Numeric s;
	  for (s=0,i=0;i<num_taps;i++) s += coeff[i];
	  return(s);
  }
  //! Get current output
  Numeric out() { return(output); }
  //! Clock in new sample & compute current output
  Numeric check(long i) { return(z[i]); }
  ~fir(void) {
	if (num_taps>0) {
	  delete [] z;
	  delete [] coeff;
	}
  }
  //! Constructor
  fir(void) : coeff(NULL), z(NULL) { ;};
  //! Constructor
  fir(long n) : coeff(NULL), z(NULL), num_taps(n)
	{
	  int i;
	  if (n>0) {
		coeff = new Numeric[n];	
		z = new Numeric[n];
		for (i=0;i<n;i++) z[i] = coeff[i] = 0;  
	  }
	}
  //! Set size of Filter 
  void set_size(long n) 
	{
	  int i;
	  num_taps = n;
	  if (n>0) {
		coeff = new Numeric[n];
		z = new Numeric[n];
		for (i=0;i<n;i++) z[i] = coeff[i] = 0;
	  }
	}
  long get_size(void) { return(num_taps); }
  //!  Constructor that gets coefficients from file (requires fir.cpp) 
  fir(const char* file) { read_taps(file); };               
  //! Update filter by inputting 1 sample and returning convolved output sample.
  Numeric clock(Numeric in) { return(update(in)); }
  Numeric update(Numeric in) {
	int i;                                                       
	// Update history of inputs
	for (i=num_taps-1;i>0;i--) z[i] = z[i-1];  
	// Add new input
	z[0] = in;   
	// Perform FIR
	for (output=0,i=0;i<num_taps;i++) output += coeff[i]*z[i];
	return(output);
  }
  // Tapped delay line uses previous outputs (to behave like an IIR)
  Numeric iir(Numeric in) {
	  int i;
	  for (output=0,i=0;i<num_taps;i++) output += coeff[i]*z[i];
	  // Update history of outputs
	  for (i=num_taps-1;i>0;i--) z[i] = z[i-1];  
	  output += in;
	  // Add new output to delay line
	  z[0] = output;   
	  return(output);
  }
  int read_taps(const char* file);
  int read_complex_taps(const char* file);
  void print(void);

  template <class N> friend vector<N> Vector_taps(fir<N> x);
  template <class N> friend vector<N> Vector_input(fir<N> y);
  void settap(vector<Numeric> z) {
	  for (int i=0;i<num_taps;i++) coeff[i] = z[i]; 
  }  
};    

template <class T> vector<T> Vector_taps(fir<T> f) {
	long N = f.num_taps;
	vector<T> V(N);
	for (int i=0;i<N;i++) V[i] = f.coeff[i];
	return(V);
}
template <class T> vector<T> Vector_input(fir<T> f) {
	long N = f.num_taps;
	vector<T> V(N);
	for (int i=0;i<N;i++) V[i] = f.z[i];
	return(V);
}

void create_remez_lpfir(fir<double>& remezfir, double* edge, double* fx, double* wtx);

#endif /*!FILTER_H*/
