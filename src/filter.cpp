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
#include "filter.h"



double* remez_fir::remez(int numtaps, int numband,
						 double bands[], double des[], double weight[], int type) {

  double c;
  int i;
  int symmetry = (type == BANDPASS)? POSITIVE_F : NEGATIVE_F;
  int r = numtaps/2;   // number of extrema
  if ((numtaps%2 != 0) && (symmetry == POSITIVE_F)) r++;
  
  // Predict dense grid size in advance for array sizes
  int gridSize = 0;
  for (i=0; i<numband; i++) {
	gridSize += (int)floor(0.5+2*r*GRIDDENSITY*(bands[2*i+1] - bands[2*i]));
  }
  if (symmetry == NEGATIVE_F) gridSize--;
  
  double* grid = new double[gridSize];
  double* d    = new double[gridSize];
  double* w    = new double[gridSize];
  double* e    = new double[gridSize];
  double* x    = new double[r+1];
  double* y    = new double[r+1];
  double* ad   = new double[r+1];
  double* taps = new double[r+1];
  int*    ext  = new int[r+1];
  
  // Create dense frequency grid
  createDenseGrid(r, numtaps, numband, bands, des, weight,
				  gridSize, grid, d, w, symmetry);
  initialGuess(r, ext, gridSize);
  
  // For Differentiator: (fix grid)
  if (type == DIFFERENTIATOR) {
	for (i=0; i<gridSize; i++) {
	  if (d[i] > 0.0001) w[i] /= grid[i];
	}
  }
  
  // For odd or NEGATIVE_F symmetry filters, alter the
  // d[] and w[] according to Parks McClellan
  if (symmetry == POSITIVE_F) {
	if (numtaps % 2 == 0) {
	  for (i=0; i<gridSize; i++) {
		c = cos(M_PI * grid[i]);
		d[i] /= c;
		w[i] *= c;
	  }
	}
  }
  else {
	if (numtaps % 2 != 0) {
	  for (i=0; i<gridSize; i++) {
		c = sin(2*M_PI * grid[i]);
		d[i] /= c;
		w[i] *= c;
	  }
	}
	else {
	  for (i=0; i<gridSize; i++) {
		c = sin(M_PI * grid[i]);
		d[i] /= c;
		w[i] *= c;
	  }
	}
  }

  // Perform the Remez Exchange algorithm
  int iter;
  for (iter=0; iter<MAXITERATIONS; iter++) {
	calcParms(r, ext, grid, d, w, ad, x, y);
	calcError(r, ad, x, y, gridSize, grid, d, w, e);
	search(r, ext, gridSize, e);
	if (isDone(r, ext, e)) break;
  }
  if (iter == MAXITERATIONS) {
	cout << "Reached maximum iteration count.\n";
	cout << "Results may be bad\n";
  }
  
  calcParms(r, ext, grid, d, w, ad, x, y);
  
  // Find the 'taps' of the filter for use with Frequency
  // Sampling.  If odd or NEGATIVE_F symmetry, fix the taps
  // according to Parks McClellan
  for (i=0; i<=numtaps/2; i++) {
	if (symmetry == POSITIVE_F) {
	  if (numtaps%2 != 0) c = 1;
	  else c = cos(M_PI * (double)i/numtaps);
	}
	else {
	  if (numtaps%2 != 0) c = sin(2*M_PI * (double)i/numtaps);
	  else c = sin(M_PI * (double)i/numtaps);
	}
	taps[i] = computeA((double)i/numtaps, r, ad, x, y)*c;
  }
  // Frequency sampling design with calculated taps
  return freqSample(taps, numtaps, symmetry);
}
/*******************
 * createDenseGrid
 *=================
 * Creates the dense grid of frequencies from the specified bands.
 * Also creates the Desired Frequency Response function (d[]) and
 * the Weight function (w[]) on that dense grid
 *
 *
 * INPUT:
 * ------
 * int      r        - 1/2 the number of filter coefficients
 * int      numtaps  - Number of taps in the resulting filter
 * int      numband  - Number of bands in user specification
 * double[] bands    - User-specified band edges [2*numband]
 * double[] des      - Desired response per band [numband]
 * double[] weight   - Weight per band [numband]
 * int      symmetry - Symmetry of filter - used for grid check
 *
 * OUTPUT:
 * -------
 * int    gridSize   - Number of elements in the dense frequency grid
 * double[] grid     - Frequencies (0 to 0.5) on the dense grid [gridSize]
 * double[] d        - Desired response on the dense grid [gridSize]
 * double[] w        - Weight function on the dense grid [gridSize]
 *******************/
void remez_fir::createDenseGrid(int r, int numtaps, int numband, double bands[],
									 double des[], double weight[], int gridSize,
									 double grid[], double d[], double w[],
									 int symmetry) {
  double lowf, highf;
  double delf = 0.5/(GRIDDENSITY*r);
  
  // For differentiator, hilbert,
  //  symmetry is odd and grid[0] = max(delf, band[0])
  if ((symmetry == NEGATIVE_F) && (delf > bands[0]))
	bands[0] = delf;
  
    int j=0;
    int k,i;
    for (int band=0; band < numband; band++) {
      grid[j] = bands[2*band];
      lowf = bands[2*band];
      highf = bands[2*band + 1];
      k = (int)floor((highf - lowf)/delf+0.5);
      for (i=0; i<k; i++) {
        d[j] = des[band];
        w[j] = weight[band];
        grid[j] = lowf;
        lowf += delf;
        j++;
      }
      grid[j-1] = highf;
    }

    // Similar to above, if odd symmetry, last grid point can't be .5
    // - but, if there are even taps, leave the last grid point at .5
    if ((symmetry == NEGATIVE_F) &&
        (grid[gridSize-1] > (0.5 - delf)) && ((numtaps % 2) != 0)) {
          grid[gridSize-1] = 0.5-delf;
    }
  }

  /********************
   * initialGuess
   *==============
   * Places Extremal Frequencies evenly throughout the dense grid.
   *
   *
   * INPUT:
   * ------
   * int r        - 1/2 the number of filter coefficients
   * int gridSize - Number of elements in the dense frequency grid
   *
   * OUTPUT:
   * -------
   * int ext[]    - Extremal indexes to dense frequency grid [r+1]
   ********************/

  void remez_fir::initialGuess(int r, int ext[], int gridSize) {
	int i;
    for (i=0; i<=r; i++) ext[i] = i * (gridSize-1) / r;
  }


  /***********************
   * calcParms
   *===========
   *
   *
   * INPUT:
   * ------
   * int      r    - 1/2 the number of filter coefficients
   * int[]    ext  - Extremal indexes to dense frequency grid [r+1]
   * double[] grid - Frequencies (0 to 0.5) on the dense grid [gridSize]
   * double[] d    - Desired response on the dense grid [gridSize]
   * double[] w    - Weight function on the dense grid [gridSize]
   *
   * OUTPUT:
   * -------
   * double[] ad   - 'b' in Oppenheim & Schafer [r+1]
   * double[] x    - [r+1]
   * double[] y    - 'C' in Oppenheim & Schafer [r+1]
   ***********************/

  void remez_fir::calcParms(int r, int ext[], double grid[], double d[], double w[],
                  double ad[], double x[], double y[]) {
    double sign, xi, delta, denom, numer;
	int i;

    // Find x[]
    for (i=0; i<=r; i++) x[i] = cos(2*M_PI * grid[ext[i]]);

    // Calculate ad[]  - Oppenheim & Schafer eq 7.132
    int ld = (r-1)/15 + 1;         // Skips around to avoid round errors
    for (i=0; i<=r; i++) {
      denom = 1.0;
      xi = x[i];
      for (int j=0; j<ld; j++) {
        for (int k=j; k<=r; k+=ld)
          if (k != i) denom *= 2.0*(xi - x[k]);
       }
       if (fabs(denom)<0.00001) denom = 0.00001;
       ad[i] = 1.0/denom;
    }

    // Calculate delta  - Oppenheim & Schafer eq 7.131
    numer = denom = 0;
    sign = 1;
    for (i=0; i<=r; i++) {
      numer += ad[i] * d[ext[i]];
      denom += sign * ad[i]/w[ext[i]];
      sign = -sign;
    }
    delta = numer/denom;
    sign = 1;

    // Calculate y[]  - Oppenheim & Schafer eq 7.133b
    for (i=0; i<=r; i++) {
      y[i] = d[ext[i]] - sign * delta/w[ext[i]];
      sign = -sign;
    }
  }


  /*********************
   * computeA
   *==========
   * Using values calculated in CalcParms, ComputeA calculates the
   * actual filter response at a given frequency (freq).  Uses
   * eq 7.133a from Oppenheim & Schafer.
   *
   *
   * INPUT:
   * ------
   * double   freq - Frequency (0 to 0.5) at which to calculate A
   * int      r    - 1/2 the number of filter coefficients
   * double[] ad   - 'b' in Oppenheim & Schafer [r+1]
   * double[] x    - [r+1]
   * double[] y    - 'C' in Oppenheim & Schafer [r+1]
   *
   * OUTPUT:
   * -------
   * Returns double value of A[freq]
   *********************/

  double remez_fir::computeA(double freq, int r, double ad[], double x[], double y[]) {

	int i;
    double c;
    double numer = 0;
    double denom = 0;
    double xc = cos(2*M_PI * freq);
    for (i=0; i<=r; i++) {
      c = xc - x[i];
      if (fabs(c) < 1.0e-7) {
        numer = y[i];
        denom = 1;
        break;
      }
      c = ad[i]/c;
      denom += c;
      numer += c*y[i];
    }
    return numer/denom;
  }


  /************************
   * calcError
   *===========
   * Calculates the Error function from the desired frequency response
   * on the dense grid (d[]), the weight function on the dense grid (w[]),
   * and the present response calculation (A[])
   *
   *
   * INPUT:
   * ------
   * int      r        - 1/2 the number of filter coefficients
   * double[] ad       - [r+1]
   * double[] x        - [r+1]
   * double[] y        - [r+1]
   * int      gridSize - Number of elements in the dense frequency grid
   * double[] grid     - Frequencies on the dense grid [gridSize]
   * double[] d        - Desired response on the dense grid [gridSize]
   * double[] w        - Weight function on the desnse grid [gridSize]
   *
   * OUTPUT:
   * -------
   * double[] e        - Error function on dense grid [gridSize]
   ************************/

  void remez_fir::calcError(int r, double ad[], double x[], double y[],
                   int gridSize, double grid[],
                   double d[], double w[], double e[]) {
    double A;
	int i;
    for (i=0; i<gridSize; i++) {
      A = computeA(grid[i], r, ad, x, y);
      e[i] = w[i] * (d[i] - A);
    }
  }

  /************************
   * search
   *========
   * Searches for the maxima/minima of the error curve.  If more than
   * r+1 extrema are found, it uses the following heuristic (thanks
   * Chris Hanson):
   * 1) Adjacent non-alternating extrema deleted first.
   * 2) If there are more than one excess extrema, delete the
   *    one with the smallest error.  This will create a non-alternation
   *    condition that is fixed by 1).
   * 3) If there is exactly one excess extremum, delete the smaller
   *    of the first/last extremum
   *
   *
   * INPUT:
   * ------
   * int      r        - 1/2 the number of filter coefficients
   * int[]    ext      - Indexes to grid[] of extremal frequencies [r+1]
   * int      gridSize - Number of elements in the dense frequency grid
   * double[] e        - Array of error values.  [gridSize]
   *
   * OUTPUT:
   * -------
   * int[]    ext      - New indexes to extremal frequencies [r+1]
   ************************/

  void remez_fir::search(int r, int ext[], int gridSize, double e[]) {
    bool up, alt;
    int* foundExt = new int[gridSize];  /* Array of found extremals */
    int k = 0;
	int i;
    // Check for extremum at 0.
    if (((e[0]>0.0) && (e[0]>e[1])) || ((e[0]<0.0) && (e[0]<e[1])))
      foundExt[k++] = 0;

    // Check for extrema inside dense grid
    for (i=1; i<gridSize-1; i++) {
      if (((e[i]>=e[i-1]) && (e[i]>e[i+1]) && (e[i]>0.0)) ||
          ((e[i]<=e[i-1]) && (e[i]<e[i+1]) && (e[i]<0.0)))
            foundExt[k++] = i;
    }

    // Check for extremum at 0.5
    int j = gridSize-1;
    if (((e[j]>0.0) && (e[j]>e[j-1])) || ((e[j]<0.0) && (e[j]<e[j-1])))
      foundExt[k++] = j;

    // Remove extra extremals
    int extra = k - (r+1);
    int l;
    while (extra > 0) {
      up = e[foundExt[0]] > 0.0;
      // up = true -->  first one is a maximum
      // up = false --> first one is a minimum
      l=0;
      alt = true;
      for (j=1; j<k; j++) {
        if (fabs(e[foundExt[j]]) < fabs(e[foundExt[l]]))
          l = j;               // new smallest error.
        if (up && (e[foundExt[j]] < 0.0))
            up = false;             // switch to a minima
        else if (!up && (e[foundExt[j]] > 0.0))
            up = true;             // switch to a maxima
        else {
          alt = false;
          break;              // Ooops, found two non-alternating
        }                     // extrema.  Delete smallest of them
      }  // if the loop finishes, all extrema are alternating

      // If there's only one extremal and all are alternating,
      // delete the smallest of the first/last extremals.
      if (alt && (extra == 1)) {
        if (fabs(e[foundExt[k-1]]) < fabs(e[foundExt[0]]))
          l = foundExt[k-1];   // Delete last extremal
        else
          l = foundExt[0];     // Delete first extremal
      }

      for (j=l; j<k; j++)      // Loop that does the deletion
        foundExt[j] = foundExt[j+1];
      k--;
      extra--;
    }
    //  Copy found extremals to ext[]
    for (i=0; i<=r; i++) ext[i] = foundExt[i];
  }


  /*********************
   * freqSample
   *============
   * Simple frequency sampling algorithm to determine the impulse
   * response h[] from A's found in ComputeA
   *
   *
   * INPUT:
   * ------
   * int      N        - Number of filter coefficients
   * double[] A        - Sample points of desired response [N/2]
   * int      symmetry - Symmetry of desired filter
   *
   * OUTPUT:
   * -------
   * double[] h        - Impulse Response of final filter [N]
   *********************/
  double* remez_fir::freqSample(double A[], int numtaps, int symm) {
    double x, val;
    int N = numtaps;
    double M = (N-1.0)/2.0;
    double* h = new double[N];
    if (symm == POSITIVE_F) {
      if (N%2 != 0) {
         for (int n=0; n<N; n++) {
           val = A[0];
           x = 2*M_PI * (n - M)/N;
           for (int k=1; k<=M; k++) val += 2.0 * A[k] * cos(x*k);
           h[n] = val/N;
         }
      }
      else {
        for (int n=0; n<N; n++) {
          val = A[0];
          x = 2*M_PI * (n - M)/N;
          for (int k=1; k<=(N/2-1); k++) val += 2.0 * A[k] * cos(x*k);
          h[n] = val/N;
        }
      }
    }
    else {
      if (N%2 != 0) {
        for (int n=0; n<N; n++) {
          val = 0;
          x = 2*M_PI * (n - M)/N;
          for (int k=1; k<=M; k++) val += 2.0 * A[k] * sin(x*k);
          h[n] = val/N;
        }
      }
      else {
        for (int n=0; n<N; n++) {
          val = A[N/2] * sin(M_PI * (n - M));
          x = 2*M_PI * (n - M)/N;
          for (int k=1; k<=(N/2-1); k++) val += 2.0 * A[k] * sin(x*k);
          h[n] = val/N;
        }
      }
    }
    return h;
  }

  /*******************
   * isDone
   *========
   * Checks to see if the error function is small enough to consider
   * the result to have converged.
   *
   * INPUT:
   * ------
   * int      r   - 1/2 the number of filter coeffiecients
   * int[]    ext - Indexes to extremal frequencies [r+1]
   * double[] e   - Error function on the dense grid [gridSize]
   *
   * OUTPUT:
   * -------
   * Returns true if the result converged
   * Returns false if the result has not converged
   ********************/

  bool remez_fir::isDone(int r, int ext[], double e[]) {

	int i;
    double min, max, current;
    min = max = fabs(e[ext[0]]);
    for (i=1; i<=r; i++){
      current = fabs(e[ext[i]]);
      if (current < min) min = current;
      if (current > max) max = current;
    }
    if (((max-min)/max) < 0.0001) return true;
    return false;
  }

void create_remez_lpfir(fir<double>& remezfir, double* edge, double* fx, double* wtx) {
  
  long nfilt = remezfir.num_taps;
  remez_fir Remz;
  double* fir_coef;
  fir_coef = Remz.remez(nfilt,2,edge,fx,wtx,1);
  for (int i=0;i<nfilt;i++) remezfir.settap(i,fir_coef[i]);
}


template <> int fir< complex<double> >::read_complex_taps(const char* file)
{
// Assumes coeficients are complex.
	int i=0;
	double tmp;
	double tmpi,tmpr;
	num_taps = 0;

	ifstream firf(file);
	if (!firf) {
	  cout << "Could not open file " << file << "\n";
	  return(-1);
	}
	while (!firf.eof()) {
			firf >>	tmp;
			firf >> tmp;
			num_taps++;
	}
	firf.close();

	coeff = new complex<double>[num_taps];
	z = new complex<double>[num_taps];

	firf.open(file);
	while (!firf.eof()) {
			firf >> tmpr;
			firf >> tmpi;
			coeff[i++] = complex<double>(tmpr,tmpi);
	}							
	firf.close();
   
	for (i=0;i<num_taps;i++) z[i] = 0;  
	return(0);
}
template <> int fir< complex<long> >::read_taps(const char* file)
{
// Assumes coeficients are real ONLY.
	int i=0;
	long tmp;
	num_taps = 0;

	ifstream firf(file);
	if (!firf) {
	  cout << "Could not open file " << file << "\n";
	  return(-1);
	}
	while (!firf.eof()) {
	  firf >>	tmp;
	  num_taps++;
	}
	firf.close();

	coeff = new complex<long>[num_taps];
	z = new complex<long>[num_taps];

	firf.open(file);
	while (!firf.eof()) {
	  firf >> tmp;
	  coeff[i++] = tmp;
	}							
	firf.close();
   
	for (i=0;i<num_taps;i++) z[i] = 0;  
	return(0);
}
template <> int fir< complex<double> >::read_taps(const char* file)
{
// Assumes coeficients are real ONLY.
	int i=0;
	double tmp;
	num_taps = 0;

	ifstream firf(file);
	if (!firf) {
	  cout << "Could not open file " << file << "\n";
	  return(-1);
	}
	while (!firf.eof()) {
	  firf >>	tmp;
	  num_taps++;
	}
	firf.close();

	coeff = new complex<double>[num_taps];
	z = new complex<double>[num_taps];

	ifstream firfx(file);
	while (!firfx.eof()) {
	  firfx >> tmp;
	  coeff[i++] = tmp;
	}							
	firfx.close();
	
	for (i=0;i<num_taps;i++) z[i] = 0;  
	return(0);
}
template <> int fir<long>::read_taps(const char* file)
{
	int i=0;
	long tmp;
	num_taps = 0;

	ifstream firf(file);
	if (!firf) {
	  cout << "Could not open file " << file << "\n";
	  return(-1);
	}
	while (!firf.eof()) {
			firf >>	tmp;
			num_taps++;
	}
	firf.close();

	coeff = new long[num_taps];
	z = new long[num_taps];

	firf.open(file);
	while (!firf.eof()) firf >> coeff[i++];
	firf.close();
   
	for (i=0;i<num_taps;i++) z[i] = 0;  
	return(0);
}
template <> int fir<double>::read_taps(const char* file)
{
	int i=0;
	double tmp;
	num_taps = 0;

	ifstream firf(file);
	if (! firf) {
	  cout << "Error opening file " << file << "\n";
	  return(-1);
	}

	while (!firf.eof()) {
	  firf >>	tmp;
	  num_taps++;
	}
	firf.close();
	
	coeff = new double[num_taps];
	z = new double[num_taps];

	ifstream firfx(file);
	if (! firfx) {
	  cout << "Error opening file " << file << "\n";
	}
	while (!firfx.eof()) {
	  firfx >> coeff[i++];
	}
	firfx.close();
   
	for (i=0;i<num_taps;i++) z[i] = 0;  
	return(0);
}   
template <> void fir<double>::print() {
	cout << "FIR filter coefficients" << '\n';
	for (long i=0;i<num_taps;i++) {
		cout << coeff[i] << cout.width(10) << ' ';
		if ((i+1)%6 == 0) cout << '\n';
	}
	cout << '\n';
	cout.flush();
}
template <> void fir<complex<double> >::print() {
	long i;
	cout << "Real FIR filter coefficients" << '\n';
	for (i=0;i<num_taps;i++) {
		cout << real(coeff[i]) << cout.width(10) << ' ';
		if ((i+1)%6 == 0) cout << '\n';
	}
	cout << '\n';
	cout << "Imaginary FIR filter coefficients" << '\n';
	for (i=0;i<num_taps;i++) {
		cout << imag(coeff[i]) << cout.width(10) << ' ';
		if ((i+1)%6 == 0) cout << '\n';
	}	
	cout << '\n';
	cout.flush();
}    
template <> void fir<complex<long> >::print() {
	long i;
	cout << "Real FIR filter coefficients" << '\n';
	for (i=0;i<num_taps;i++) {
		cout << (long)real(coeff[i]) << cout.width(10) << ' ';
		if ((i+1)%6 == 0) cout << '\n';
	}
	cout << '\n';
	cout << "Imaginary FIR filter coefficients" << '\n';
	for (i=0;i<num_taps;i++) {
		cout << (long)imag(coeff[i]) << cout.width(10) << ' ';
		if ((i+1)%6 == 0) cout << '\n';
	}	
	cout << '\n';
	cout.flush();
}
