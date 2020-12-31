#include <synth_turb/SynthTurb3d_periodic_box.hpp>
#include <synth_turb/SynthTurb3d_all_waves.hpp>
#include <iostream>

#define NX 3
//#define NX 50
#define DX 1
#define _NMODES NX
//#define _NMODES 200
#define _NWAVES 6//3
//#define _NWAVES 50

template<template<class, int, int> class SynthTurb_t>
void test()
{
  SynthTurb_t<double, _NMODES, _NWAVES> rm_d(1e-4, (NX-1)*DX, 1e-3); // eps [m2/s3], Lmax [m], Lmin[m] (Lmin has no role in the periodic version)
  // Lmax = (nx-1)*dx zeby sprawdzic warunek na periodycznosc
  rm_d.generate_random_modes();

  double u[NX][NX][NX],
        v[NX][NX][NX],
        w[NX][NX][NX];

  double t=0;

  rm_d.generate_velocity_field(u,v,w,DX,t);

  std::cout << "u:" << std::endl;
  for(int i=0; i<NX; ++i)
  {
    for(int j=0; j<NX; ++j)
    {
      for(int k=0; k<NX; ++k)
        std::cout << u[i][j][k] << " ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  std::cout << "v:" << std::endl;
  for(int i=0; i<NX; ++i)
  {
    for(int j=0; j<NX; ++j)
    {
      for(int k=0; k<NX; ++k)
        std::cout << v[i][j][k] << " ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  std::cout << "w:" << std::endl;
  for(int i=0; i<NX; ++i)
  {
    for(int j=0; j<NX; ++j)
    {
      for(int k=0; k<NX; ++k)
        std::cout << w[i][j][k] << " ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  std::cout << "u[0][0][0]: " << u[0][0][0] << std::endl;
  std::cout << "u[NX-1][0][0]: " << u[NX-1][0][0] << std::endl;
}

int main()
{
  std::cout << std::endl << "PERIODIC TEST" << std::endl << std::endl;
  test<SynthTurb::SynthTurb3d_periodic_box>();

  std::cout << std::endl << "ALL WAVES TEST" << std::endl << std::endl;
  test<SynthTurb::SynthTurb3d_all_waves>();
}
