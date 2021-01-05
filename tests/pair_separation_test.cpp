#include <synth_turb/SynthTurb3d_periodic_box.hpp>
#include <synth_turb/SynthTurb3d_all_waves.hpp>
#include <synth_turb/RandTurbGA17.hpp>
#include <synth_turb/SynthTurb3d_periodic_box_multiwave.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include <chrono>
#include <omp.h>

//#define _NMODES 200
//#define _NWAVES 5//50

#define NPairs 100 // 1000 // number of test particle pairs
#define InitSep 1 // initial separation [in units of the Kolmogorov length]
#define LKol 1e-3 // Kolmogorov length [m]
#define Lmax 1 // integral length [m]
#define DT 0.1 // [s]
#define T 100 // [s]
#define EPS 1e-1 // TKE diss rate [m2/s3]


class tester_common
{
  protected:

  std::random_device rd;
  std::default_random_engine rand_eng;//(rd());
  std::uniform_real_distribution<double> x_d,//(0, Lmax),
                                         h_d,//(-1, std::nextafter(1, std::numeric_limits<double>::max())), // uniform in [-1,1]
                                         th_d;//(0, std::nextafter(2 * M_PI, std::numeric_limits<double>::max()));  // uniform in [0,2*Pi]
                                         
  std::ofstream ofs;// (outfile, std::ofstream::out);

  double x[NPairs][2][2][3]; // positions [Npairs][NdropsInPair][n/n+1(pred-corr)][NDims]
  double v[NPairs][2][2][3]; // velocities
  double r[NPairs]; // separations
  double e[3]; // unit vector

  double mean_r, sig_r;

  virtual void predictor(const int&, const double&)=0;
  virtual void corrector(const int&, const double&)=0;
  virtual void update_time(const double&)=0;

  public:
  //ctor
  tester_common(const std::string &outfile):
    rand_eng(rd()),
    x_d(0, Lmax),
    h_d(-1, std::nextafter(1, std::numeric_limits<double>::max())),
    th_d(0, std::nextafter(2 * M_PI, std::numeric_limits<double>::max())),
    ofs(outfile, std::ofstream::out)
  {
    ofs << "#time\t<r>\tsig(r)" << std::endl;
    ofs << std::setprecision(5) << std::setw(5) << 0 << "\t" << InitSep*LKol << "\t" << 0 << std::endl;

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
  }

  void test()
  {
    auto t1 = std::chrono::high_resolution_clock::now();
    // time loop
    //  const int n_step = T / DT;
    double t=0;
    while(t<=T)
    {
      #pragma omp parallel for
      for(int p=0; p<NPairs; ++p)
      {
        this->predictor(p, t);
      }

      t+=DT;

      this->update_time(t);

      #pragma omp parallel for
      for(int p=0; p<NPairs; ++p)
      {
        this->corrector(p, t);
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
  
  //    std::cout << "t: " << t+DT << " <r>: " << mean_r << " sig(r): " << sig_r << std::endl;
      ofs << std::setprecision(5) << std::setw(5) << t+DT << "\t" << mean_r << "\t" << sig_r << std::endl;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    std::cout << "wall time: " << duration << " [ms] OpenMP threads: " << omp_get_max_threads() <<  std::endl;
  }
};

template<template<class,int,int> class SynthTurb_t, int NModes, int NWaves>
class tester_synth_turb_common : public tester_common
{
  using parent_t = tester_common;

  virtual void calculate_velocity(const int &p, const int &d, const int &n, const double &t)=0;

  void predictor(const int &p, const double &t) override
  {
    for(int d=0; d<2; ++d)
    {
      // predictor
      this->calculate_velocity(p,d,0,t); 
      for(int i=0; i<3; ++i)
      {
        x[p][d][1][i] = x[p][d][0][i] + v[p][d][0][i] * DT; 
      }
    }
  }

  void corrector(const int &p, const double &t) override
  {
    for(int d=0; d<2; ++d)
    {
      // corrector
      this->calculate_velocity(p,d,1,t);
      for(int i=0; i<3; ++i)
      {
        x[p][d][0][i] += 0.5 * (v[p][d][0][i] + v[p][d][1][i]) * DT; 
      }
    }
  }

  protected:

  SynthTurb_t<double, NModes, NWaves> st;//(EPS, Lmax, LKol); // eps [m2/s3?], Lmax [m], Lmin[m] (Lmin has no role in the periodic version)

  public:
  //ctor
  tester_synth_turb_common(const std::string &outfile):
    parent_t(outfile),
    st(EPS, Lmax, LKol)
  {
    std::cout << outfile << " test, NModes: " << NModes << " NWaves: " << NWaves << std::endl;
  }
};

template<template<class,int,int> class SynthTurb_t, int NModes, int NWaves>
class tester_synth_turb : public tester_synth_turb_common<SynthTurb_t, NModes, NWaves>
{
  using parent_t = tester_synth_turb_common<SynthTurb_t, NModes, NWaves>;
  using parent_t::parent_t;

  void calculate_velocity(const int &p, const int &d, const int &n, const double &t) override
  {
    this->st.calculate_velocity(this->v[p][d][n][0], this->v[p][d][n][1], this->v[p][d][n][2], this->x[p][d][n], t);
  }

  void update_time(const double &t) override
  {}
};


template<template<class,int,int> class SynthTurb_t, int NModes, int NWaves>
class tester_synth_turb_multiwave : public tester_synth_turb_common<SynthTurb_t, NModes, NWaves>
{
  using parent_t = tester_synth_turb_common<SynthTurb_t, NModes, NWaves>;
  using parent_t::parent_t;

  void calculate_velocity(const int &p, const int &d, const int &n, const double &t) override
  {
    this->st.calculate_velocity(this->v[p][d][n][0], this->v[p][d][n][1], this->v[p][d][n][2], this->x[p][d][n]);
  }

  void update_time(const double &t) override
  {
    this->st.update_time(t);
  }
};


template<template<class> class RandTurb_t>
class tester_rand_turb : public tester_common
{
  using parent_t = tester_common;

  void predictor(const int &p, const double &t) override
  {
    for(int d=0; d<2; ++d)
    {
      for(int i=0; i<3; ++i)
        rt.update_sgs_velocity(v[p][d][0][i], DT);
    }
  }

  void corrector(const int &p, const double &t) override
  {
    for(int d=0; d<2; ++d)
    {
      for(int i=0; i<3; ++i)
      {
        v[p][d][1][i] = v[p][d][0][i];
        rt.update_sgs_velocity(v[p][d][1][i], DT);
        v[p][d][0][i] = 0.5 * (v[p][d][0][i] + v[p][d][1][i]);
        x[p][d][0][i] += v[p][d][0][i] * DT; 
      }
    }
  }

  void update_time(const double &dt) override
  {}

  private:
  RandTurb_t<double> rt;

  public:
  //ctor
  tester_rand_turb(const std::string &outfile):
    parent_t(outfile),
    rt(EPS, Lmax)
  {
    std::cout << "rand_turb test" << std::endl;
  }
};

int main()
{
  std::cout <<
    "NPairs: " << NPairs <<
    " InitSep: " << InitSep <<
    " Lkol: " << LKol <<
    " Lmax: " << Lmax <<
    " DT: " << DT <<
    " T: " << T <<
    " EPS: " << EPS <<
    std::endl;

  // synth turb with periodic box flow
//  {
//    constexpr int NModes=1000,
//                  NWaves=6;
//    tester_synth_turb<SynthTurb::SynthTurb3d_periodic_box, NModes, NWaves> periodic_box("pair_separation_periodic_box.dat");
//    periodic_box.test();
//  }
//  // synth turb with periodic box flow with more waves
  {
    constexpr int NModes=1000,//200,
                  NWaves=50;
    tester_synth_turb_multiwave<SynthTurb::SynthTurb3d_periodic_box_multiwave, NModes, NWaves> periodic_box_multiwave("pair_separation_periodic_box_multiwave.dat");
    periodic_box_multiwave.test();
  }
//  // synth turb with all waves
//  {
//    constexpr int NModes=200,//200,
//                  NWaves=50;//50;
//    tester_synth_turb<SynthTurb::SynthTurb3d_all_waves, NModes, NWaves> all_waves("pair_separation_all_waves.dat");
//    all_waves.test();
//  }
//  // GA17 SGS model
//  {
//    tester_rand_turb<RandTurb::RandTurb_GA17> GA17("pair_separation_GA17.dat");
//    GA17.test();
//  }
}
