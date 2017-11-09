// Programming for making running mean.
// Input should be a single column of numbers.
// Command line input: file name, halfwindow size n (window is 2n + 1).

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include<cstring>
using namespace std;
#define maxlength 500000

int main(int argc, char* argv[])
{
  if (3 != argc) {
    cerr << "Aborting; use '" << argv[0] << " logfile " 
	 << "window_half_size'.\n";
    exit(EXIT_FAILURE);
  }

  int j, k, nlines, hfw, n;
  double col1[maxlength], sum1;

  int verbose = 0;

  hfw = atoi(argv[2]);

  // Initialize arrays
  for (k = 0; k < maxlength; ++k) {
    col1[k] = 0.;
  }

  // Open log file.
  ifstream logfile(argv[1]); 
  if (!logfile) {
    cerr << "Aborting: can't find log file '" << argv[1] << "'.\n";
    exit(EXIT_FAILURE);
  }

  // Number of lines in log file
  nlines = -1;
  while (logfile.ignore(1000, '\n')) {
    ++nlines;
  }
  logfile.clear();
  logfile.seekg(0);
  if (verbose) {
    cout << "There are " << nlines << " data lines.\n\n";
  }
  if (nlines > maxlength) {
    cerr << "Input file is too large. Max data lines = " << maxlength
	 << ".\n";
        exit(EXIT_FAILURE);
  }

  // Read in logfile
  for (k = 0; k < nlines; ++k) {
    logfile >> col1[k];
    logfile.ignore(1000, '\n');
  }

  // Running mean
  for (k = 0; k < nlines; ++k) {
    sum1 = 0.;
    n = 0;
    for (j = k - hfw; j <= k + hfw; ++j) { 
      if (j >= 0 && j < nlines) {
	sum1 += col1[j];
	n += 1;
      }
    }
    cout << col1[k] << " " << sum1/float(n) << "\n";
  }

}
