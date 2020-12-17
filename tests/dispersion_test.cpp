#include <synth_turb/SynthTurb3d_periodic_box.hpp>
#include <synth_turb/SynthTurb3d_all_waves.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>

#define _NMODES 200
#define _NWAVES 5//50

#define NPairs 1000 // number of test particle pairs
#define InitSep 1 // initial separation [in units of the Kolmogorov length]
#define LKol 1e-3 // Kolmogorov length [m]
#define Lmax 1 // integral length [m]
#define DT 0.01 // [s]
#define T 10 // [s]
#define EPS 1e-2 // TKE diss rate [m2/s3]


// TODO: 3 tests: all_waves, periodic_box, GA17

int main()
{
//  SynthTurb::SynthTurb3d_periodic_box<double, _NMODES, _NWAVES> rm_d(EPS, Lmax, LKol); // eps [m2/s3?], Lmax [m], Lmin[m] (Lmin has no role in the periodic version)
  SynthTurb::SynthTurb3d_all_waves<double, _NMODES, _NWAVES> rm_d(EPS, Lmax, LKol); // eps [m2/s3?], Lmax [m], Lmin[m] (Lmin has no role in the periodic version)

  rm_d.generate_random_modes();
  std::random_device rd;

  std::default_random_engine rand_eng(rd());
  std::uniform_real_distribution<double> x_d(0, Lmax),
                                         h_d(-1, std::nextafter(1, std::numeric_limits<double>::max())), // uniform in [-1,1]
                                         th_d(0, std::nextafter(2 * M_PI, std::numeric_limits<double>::max()));  // uniform in [0,2*Pi]
                                         
  std::ofstream ofs ("dispersion_tests.dat", std::ofstream::out);
  ofs << "#time\t<r>\tsig(r)" << std::endl;
  ofs << std::setprecision(5) << 0 << "\t" << InitSep*LKol << "\t" << 0 << std::endl;



  double x[NPairs][2][2][3]; // positions [Npairs][NdropsInPair][n/n+1(pred-corr)][NDims]
  double v[NPairs][2][2][3]; // velocities
  double r[NPairs]; // separations
  double e[3]; // unit vector

  // random position of the first droplet in a pair
  for(int i=0; i<3; ++i)
    for(int p=0; p<NPairs; ++p)
      x[p][0][0][i] = x_d(rand_eng);

  // positions of the second droplets seprated randomly by InitSep * LKol
  for(int p=0; p<NPairs; ++p)
  {
    // random unit vector
    double h  = h_d(rand_eng);
    double th = th_d(rand_eng);
    e[0] = sqrt(1. - h*h) * cos(th);
    e[1] = sqrt(1. - h*h) * sin(th);
    e[2] = h;
    for(int i=0; i<3; ++i)
      x[p][1][0][i] = x[p][0][0][i] + e[i] * InitSep * LKol;
  }

  double mean_r, sig_r;

  // time loop
  //  const int n_step = T / DT;
  for(double t=0; t<=T; t+=DT)
  {
    // update positions
    #pragma omp parallel for
    for(int p=0; p<NPairs; ++p)
    {
      for(int d=0; d<2; ++d)
      {
        // predictor
        rm_d.calculate_velocity_field(v[p][d][0][0], v[p][d][0][1], v[p][d][0][2], x[p][d][0], t);
        for(int i=0; i<3; ++i)
        {
          x[p][d][1][i] = x[p][d][0][i] + v[p][d][0][i] * DT; 
        }
        // corrector
        rm_d.calculate_velocity_field(v[p][d][1][0], v[p][d][1][1], v[p][d][1][2], x[p][d][1], t+DT);
        for(int i=0; i<3; ++i)
        {
          x[p][d][0][i] += 0.5 * (v[p][d][0][i] + v[p][d][1][i]) * DT; 
        }
      }
      // calculate separation
      r[p] = sqrt(
        (x[p][1][0][0] - x[p][0][0][0]) * (x[p][1][0][0] - x[p][0][0][0]) + 
        (x[p][1][0][1] - x[p][0][0][1]) * (x[p][1][0][1] - x[p][0][0][1]) + 
        (x[p][1][0][2] - x[p][0][0][2]) * (x[p][1][0][2] - x[p][0][0][2]) 
      );
    }

    // calculate dispersion statistics
    mean_r=0;
    sig_r=0;

    // mean
    for(int p=0; p<NPairs; ++p)
    {
      mean_r += r[p];
    }
    mean_r /= NPairs;

    // std dev
    for(int p=0; p<NPairs; ++p)
    {
      sig_r += (r[p] - mean_r) * (r[p] - mean_r);
    }
    sig_r = sqrt(sig_r / NPairs);

    std::cout << "t: " << t+DT << " <r>: " << mean_r << " sig(r): " << sig_r << std::endl;
    ofs << std::setprecision(5) << t+DT << "\t" << mean_r << "\t" << sig_r << std::endl;
  }
}
