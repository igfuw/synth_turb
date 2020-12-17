// Generating a velocity field from a prescribed energy spectrum,
// Sources: Sidin et al. in Physics of Fluids (2009), Zhou et al. Phys. of Fluids (2018), Fortran implementation from Gustavo C. Abade

// TODO: uzyc poprawki do spektrum energii dla malych skal?

#pragma once

#include <random>
#include <cmath>

#include <iostream>

namespace SynthTurb
{
  template<class real_t, int Nmodes, int Nwaves> // Nmodes - number of Fourier modes, Nwaves - number of wave vectors in each mode
  class SynthTurb3d_common
  {
    private:

    real_t AA[3],
           BB[3]; // magnitues of A and B coefficients

    real_t dk[Nmodes],  // differences between norms of subsequent wave vectors
           E[Nmodes],   // kinetic energy in a mode
           var[Nmodes], // variances 
           wn[Nmodes];   // frequencies 

    const real_t lambda = 1; // unsteadiness parameter, within [0,1]; see Sidin et al. 2009

    real_t Anm[3][Nmodes][Nwaves],
           Bnm[3][Nmodes][Nwaves],
           knm[3][Nmodes][Nwaves];

    protected:
    std::default_random_engine rand_eng;
    std::uniform_real_distribution<real_t> h_d, th_d; // random numbers used in generating e[3]

    real_t k[Nmodes];   // norms of wave vectors
    real_t e[3]; // random unit vector

    virtual void generate_wavenumbers(const real_t &Lmax, const real_t &Lmin)=0;
    virtual void generate_unit_wavevectors(const int &m)=0;

    // init, can't be in the constructor because it calls virtual functions
    void init(
      const real_t &eps,        // TKE dissipation rate [m2/s3]
      const real_t &Lmax = 100, // maximum length scale [m]
      const real_t &Lmin = 1e-3 // Kolmogorov length scale [m]
    )
    {
      this->generate_wavenumbers(Lmax, Lmin);

      // Energy spectrum; first equation in Appendix of Sidin et al. 2009, but without the corrections f_L and f_eta
      for(int n=0; n<Nmodes; ++n)
        E[n] = 1.44 * pow(eps, 2./3.) * pow(k[n], -5./3.);
      
      // Wave vector differences
      dk[0] =         0.5 *  (k[1]        - k[0]);
      dk[Nmodes-1] =  0.5 *  (k[Nmodes-1] - k[Nmodes-2]);
      
      for(int n=1; n<Nmodes-1; ++n)
        dk[n] = 0.5 * (k[n+1] - k[n-1]);

      // Variances
      for(int n=0; n<Nmodes; ++n)
        var[n] = E[n] * dk[n] / Nwaves;

      // Frequencies; Eq. A7 in Sidin et al. 2009
      for(int n=0; n<Nmodes; ++n)
        wn[n] = lambda * sqrt(E[n] * k[n] * k[n] * k[n]);


      // std::cerr << "SynthTurb3d ctor debug output:" << std::endl;
     // for(int n=0; n<Nmodes; ++n)
     // {
     //    std::cerr << "n: " << n << " k[n]: " << k[n] << std::endl;
     //    std::cerr << "n: " << n << " E[n]: " << E[n] << std::endl;
     //    std::cerr << "n: " << n << " dk[n]: " << dk[n] << std::endl;
     //    std::cerr << "n: " << n << " var[n]: " << var[n] << std::endl;
     //    std::cerr << "n: " << n << " wn[n]: " << wn[n] << std::endl;
     //    std::cerr << std::endl;
     // }
    }

    public:
 
    // ctor
    SynthTurb3d_common():
      rand_eng(std::random_device()()),
      h_d(-1, std::nextafter(1, std::numeric_limits<real_t>::max())), // uniform in [-1,1]
      th_d(0, std::nextafter(2 * M_PI, std::numeric_limits<real_t>::max()))  // uniform in [0,2*Pi]
      {}

    // fills the kn, An and Bn arrays
    void generate_random_modes()
    {
      for(int n=0; n<Nmodes; ++n)
      {
        std::normal_distribution<real_t> G_d(0, sqrt(var[n]));

        for(int m=0; m<Nwaves; ++m)
        {
          this->generate_unit_wavevectors(m);

          // knm = random vector * magnitude
          for(int i=0; i<3; ++i)
            knm[i][n][m] = e[i] * k[n];

          // calculate coefficients Anm and Bnm - see Zhou et al.
          real_t AA[3];
          for(int i=0; i<3; ++i)
            AA[i] = G_d(rand_eng);

          Anm[0][n][m] = AA[1] * e[2] - AA[2] * e[1];
          Anm[1][n][m] = AA[2] * e[0] - AA[0] * e[2];
          Anm[2][n][m] = AA[0] * e[1] - AA[1] * e[0];

          real_t BB[3];
          for(int i=0; i<3; ++i)
            BB[i] = G_d(rand_eng);

          Bnm[0][n][m] = BB[1] * e[2] - BB[2] * e[1];
          Bnm[1][n][m] = BB[2] * e[0] - BB[0] * e[2];
          Bnm[2][n][m] = BB[0] * e[1] - BB[1] * e[0];
        }
      }
    }

//    // generate velocity field assuming uniform grid size and spacing, Arakawa-C staggering 
//    // TODO: periodic bcond (velocities at opposite edges are equal)
//    template<int nx>
//    void generate_velocity_field_ArakawaC(real_t u[nx+1][nx][nx], real_t v[nx][nx+1][nx], real_t w[nx][nx][nx+1], const real_t &dx, const real_t &t)
//    {
//      #pragma omp parallel for
//      for(int i=0; i<nx; ++i)
//        for(int j=0; j<nx; ++j)
//          for(int k=0; k<nx; ++k)
//          {
//            u[i][j][k] = 0;
//            v[i][j][k] = 0;
//            w[i][j][k] = 0;
//
//            for(int n=0; n<Nmodes; ++n)
//            {
//              const real_t wnt = wn[n] * t;
//              //// std::cerr << i << " " << j << " " << k << " wnt: " << wnt << std::endl;
//              for(int m=0; m<Nwaves; ++m)
//              {
//                const real_t xu = (knm[0][n][m] * i * dx + knm[1][n][m] * (j+0.5) * dx + knm[2][n][m] * (k+0.5) * dx) + wnt;
//                u[i][j][k] += Anm[0][n][m]*cos(xu) + Bnm[0][n][m]*sin(xu); 
//                // std::cerr << i << " " << j << " " << k << " x: " << x << std::endl;
//
//                const real_t xv = (knm[0][n][m] * (i+0.5) * dx + knm[1][n][m] * j * dx + knm[2][n][m] * (k+0.5) * dx) + wnt;
//                v[i][j][k] += Anm[1][n][m]*cos(xv) + Bnm[1][n][m]*sin(xv); 
//
//                const real_t xw = (knm[0][n][m] * (i+0.5) * dx + knm[1][n][m] * (j+0.5) * dx + knm[2][n][m] * k * dx) + wnt;
//                w[i][j][k] += Anm[2][n][m]*cos(xw) + Bnm[2][n][m]*sin(xw); 
//              }
//            }
//          }
//    };


    void calculate_velocity_field(real_t &u, real_t &v, real_t &w, const real_t x[3], const real_t &t)
    {
      {
        u = 0;
        v = 0;
        w = 0;

        for(int n=0; n<Nmodes; ++n)
        {
          const real_t wnt = wn[n] * t;
          //// std::cerr << i << " " << j << " " << k << " wnt: " << wnt << std::endl;
          for(int m=0; m<Nwaves; ++m)
          {
            const real_t r = (knm[0][n][m] * x[0] + knm[1][n][m] * x[1] + knm[2][n][m] * x[2]) + wnt;
            u += Anm[0][n][m]*cos(r) + Bnm[0][n][m]*sin(r); 
            v += Anm[1][n][m]*cos(r) + Bnm[1][n][m]*sin(r); 
            w += Anm[2][n][m]*cos(r) + Bnm[2][n][m]*sin(r); 
          }
        }
      }
    };

    template<int nx>
    void generate_velocity_field(real_t u[nx][nx][nx], real_t v[nx][nx][nx], real_t w[nx][nx][nx], const real_t &dx, const real_t &t)
    {
      #pragma omp parallel for
      for(int i=0; i<nx; ++i)
        for(int j=0; j<nx; ++j)
          for(int k=0; k<nx; ++k)
          {
            const real_t x[3] = {i*dx, j*dx, k*dx};
            calculate_velocity_field(u[i][j][k], v[i][j][k], w[i][j][k], x, t);
          }
    };
  };
};
