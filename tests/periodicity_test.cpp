#include <synth_turb/SynthTurb3d_periodic_box.hpp>
#include <synth_turb/SynthTurb3d_all_waves.hpp>
#include <synth_turb/SynthTurb3d_periodic_box_multiwave.hpp>
#include <iostream>
#include <chrono>
#include <cassert>
#include <omp.h>

#define NX 60
#define DX 1
#define _NMODES NX-1
#define _NWAVES 6//3
#define TIME 10
#define EPSILON 1e-10

bool AreClose(double a, double b)
{
    return fabs(a - b) < EPSILON;
}

template<template<class, int, int> class SynthTurb_t>
void test(const bool check_periodicity)
{

  // Lmax = (nx-1)*dx zeby sprawdzic warunek na periodycznosc
  SynthTurb_t<double, _NMODES, _NWAVES> st(1e-4, (NX-1)*DX, DX); // eps [m2/s3], Lmax [m], Lmin[m]

  double u[NX][NX][NX],
        v[NX][NX][NX],
        w[NX][NX][NX];

  auto t1 = std::chrono::high_resolution_clock::now();
  for(double time = 1; time <= TIME; time+=1)
  {
    st.generate_velocity_field(u,v,w,DX,time);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  std::cout << "generate_velocity_field wall time: " << duration << " [ms] OpenMP threads: " << omp_get_max_threads() <<  std::endl;

  if(check_periodicity)
  {
    for(int i=0; i<NX; ++i)
      for(int j=0; j<NX; ++j)
      {
        assert(AreClose(u[i][j][0], u[i][j][NX-1]));
        assert(AreClose(v[i][j][0], v[i][j][NX-1]));
        assert(AreClose(w[i][j][0], w[i][j][NX-1]));
      }
    for(int i=0; i<NX; ++i)
      for(int k=0; k<NX; ++k)
      {
        assert(AreClose(u[i][0][k], u[i][NX-1][k]));
        assert(AreClose(v[i][0][k], v[i][NX-1][k]));
        assert(AreClose(w[i][0][k], w[i][NX-1][k]));
      }
    for(int j=0; j<NX; ++j)
      for(int k=0; k<NX; ++k)
      {
        assert(AreClose(u[0][j][k], u[NX-1][j][k]));
        assert(AreClose(v[0][j][k], v[NX-1][j][k]));
        assert(AreClose(w[0][j][k], w[NX-1][j][k]));
      }
  }
}

int main()
{
  std::cout << std::endl << "PERIODIC MULTIWAVE TEST" << std::endl;
  test<SynthTurb::SynthTurb3d_periodic_box_multiwave>(true);

  std::cout << std::endl << "PERIODIC TEST" << std::endl;
  test<SynthTurb::SynthTurb3d_periodic_box>(true);

  std::cout << std::endl << "ALL WAVES TEST" << std::endl;
  test<SynthTurb::SynthTurb3d_all_waves>(false);
}
