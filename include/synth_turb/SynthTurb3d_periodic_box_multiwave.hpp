// Generating a velocity field from a prescribed energy spectrum, periodic in a box
// Sources: document by Gustavo C. Abade

// TODO: uzyc poprawki do spektrum energii dla malych skal?

#pragma once

#include <random>
#include <cmath>
#include <array>
#include <algorithm>
#include <iostream>

namespace SynthTurb
{
  // TODO: a co z vectors.size==0?
  int degeneracy_generator(const int &n, std::vector<std::array<int,3>> &vectors) // returns degeneracy, fills vectors with half of the degeneracy (opposite vectors are not added)
  {
    int degeneracy = 0;
    vectors.clear();
    const double max = sqrt(n);
    for(int nx=-max; nx<=max; ++nx)
      for(int ny=-max; ny<=max; ++ny)
        for(int nz=-max; nz<=max; ++nz)
        {
          if(nx*nx + ny*ny + nz*nz == n)
          {
 //           std::cout << "(" << nx << "," << ny << "," << nz << ")" << std::endl;
            ++degeneracy;
            if(std::find(vectors.begin(), vectors.end(), std::array<int,3>{-nx,-ny,-nz}) == vectors.end()) // add only if the opposite vector has not been added yet
              vectors.push_back({nx,ny,nz});
          }
        }
 //   std::cout << "degeneracy: " << degeneracy << std::endl;
    return degeneracy;
  }


  template<class real_t, int Nmodes, int Nwaves_max> // Nmodes - number of Fourier modes, Nwaves_max - number of wave vectors in each mode (maximum, could be less in periodic domain)
  class SynthTurb3d_periodic_box_multiwave
  {
    private:
    real_t dk[Nmodes],  // differences between norms of subsequent wave vectors
           E[Nmodes],   // kinetic energy in a mode
           std_dev[Nmodes], // sqrt(variances)
           wn[Nmodes];   // frequencies 

    real_t time;

    int Nwaves[Nmodes]; // actual number of wave vectors in thath mode


    // wavevectors in the form k = (nx,ny,nz) * 2 PI / L, where n is integer to get periodic flow 
    int nn[Nmodes]; // = nx^2 + ny^2 + nz^2
    real_t enm[3][Nmodes][Nwaves_max]; // unit vectors along wavevectors

    const real_t lambda = 1; // unsteadiness parameter, within [0,1]; see Sidin et al. 2009

    real_t Anm[3][Nmodes][Nwaves_max],
           Bnm[3][Nmodes][Nwaves_max],
           knm[3][Nmodes][Nwaves_max];

    real_t k[Nmodes];   // norms of wave vectors


    public:

    void update_time(const real_t &totime)
    {
      if(totime <= time)
        throw std::runtime_error("time cannot go back (or stand still!)");

      const real_t dt = totime - time;

      #pragma omp parallel for
      for(int n=0; n<Nmodes; ++n)
      {
        std::normal_distribution<real_t> normal_d(0,1);
        std::default_random_engine local_rand_eng(std::random_device{}());
        real_t relax = exp(-wn[n] * dt);

        for(int m=0; m<Nwaves[n]; m+=2)
        {
          for(int i=0; i<3; ++i)
          {
            Anm[i][n][m] = relax * Anm[i][n][m] + std_dev[n] * sqrt(1. - relax * relax) * normal_d(local_rand_eng);
            Anm[i][n][m+1] = -Anm[i][n][m];

            Bnm[i][n][m] = relax * Bnm[i][n][m] + std_dev[n] * sqrt(1. - relax * relax) * normal_d(local_rand_eng);
            Bnm[i][n][m+1] = Bnm[i][n][m];
          }
        }
      }
      time = totime;
    };
 
    // ctor
    SynthTurb3d_periodic_box_multiwave(
      const real_t &eps,        // TKE dissipation rate [m2/s3]
      const real_t &Lmax = 100, // maximum length scale [m]
      const real_t &Lmin = 1e-3 // Kolmogorov length scale [m]
    ) : time(0)
    {
      if(Nwaves_max % 2 != 0) throw std::runtime_error("Nwaves_max needs to be even, because we need to include opposites of all wavevectors.");

      // nn = nx^2 + ny^2 + nz^2 

      // --- linear distribution of nn (nn = 1, 2, 3, 4, ..., Nmodes) ---
   //   for(int n=0; n<Nmodes; ++n)
     //   k[n] = sqrt(n+1) * (2. * M_PI / Lmax);

      // --- geometric distribution of nn ---
      if(Nmodes > Lmax / Lmin)
        throw std::runtime_error("too many modes: Nmodes is greater than Lmax / Lmin");

      k[0] = 2. * M_PI / Lmax;
      nn[0]=1;

      real_t alpha = pow(Lmax / Lmin, 1. / (Nmodes - 1));
      while(1)
      {
        for(int n=1; n<Nmodes; ++n)
        {
          nn[n] = -1;
          int exponent = n;
          while(nn[n] <= nn[n-1])
          {
            nn[n] = std::round(std::pow(alpha, exponent++));
          }
//          std::cerr << "alpha: " << alpha << " nn[" << n << "]: " << nn[n] << std::endl;
          if(nn[n] > Lmax / Lmin) break;
        }
        if(nn[Nmodes-1] <= Lmax / Lmin && nn[Nmodes-1]!=0)
          break;
        else
          alpha /= 1.001;
      }

      for(int n=1; n<Nmodes; ++n)
      {
  //      std::cerr << "nn[" << n << "]: " << nn[n] << std::endl;
        k[n] = k[0] * sqrt(real_t(nn[n]));
      }


      std::vector<std::array<int,3>> vectors;
      for(int n=0; n<Nmodes; ++n)
      {
        Nwaves[n] = degeneracy_generator(nn[n], vectors);

        if(Nwaves[n] > Nwaves_max) // random shuffle, because not all possible degeneracies will be used
        {
          std::default_random_engine local_rand_eng(std::random_device{}());
          std::shuffle(std::begin(vectors), std::end(vectors), local_rand_eng);
          Nwaves[n] = Nwaves_max;
        }
      //  if(Nwaves_max != 6) throw std::runtime_error("nwaves max needs to be 6 for this test");
      //  vectors = {{1,0,0},{0,1,0},{0,0,1}};

        for(int m=0; m<Nwaves[n]; m+=2)
        {
          enm[0][n][m] = vectors.at(m/2)[0] / sqrt(real_t(nn[n]));
          enm[1][n][m] = vectors.at(m/2)[1] / sqrt(real_t(nn[n]));
          enm[2][n][m] = vectors.at(m/2)[2] / sqrt(real_t(nn[n]));
          // opposite vector
          enm[0][n][m+1] = -vectors.at(m/2)[0] / sqrt(real_t(nn[n]));
          enm[1][n][m+1] = -vectors.at(m/2)[1] / sqrt(real_t(nn[n]));
          enm[2][n][m+1] = -vectors.at(m/2)[2] / sqrt(real_t(nn[n]));
        }
      }

      // Energy spectrum; first equation in Appendix of Sidin et al. 2009, but without the corrections f_L and f_eta
      for(int n=0; n<Nmodes; ++n)
        E[n] = 1.44 * pow(eps, 2./3.) * pow(k[n], -5./3.);
      
      // Wave vector differences
      dk[0] =         0.5 *  (k[1]        - k[0]);
      dk[Nmodes-1] =  0.5 *  (k[Nmodes-1] - k[Nmodes-2]);
      
      for(int n=1; n<Nmodes-1; ++n)
        dk[n] = 0.5 * (k[n+1] - k[n-1]);

      // std deviations 
      for(int n=0; n<Nmodes; ++n)
      {
        std_dev[n] = sqrt(E[n] * dk[n] / Nwaves[n]);
      }

      // Frequencies; Eq. A7 in Sidin et al. 2009
      for(int n=0; n<Nmodes; ++n)
        wn[n] = lambda * sqrt(E[n] * k[n] * k[n] * k[n]);


      // std::cerr << "SynthTurb3d ctor debug output:" << std::endl;
     // for(int n=0; n<Nmodes; ++n)
     // {
     //    std::cerr << "n: " << n << " k[n]: " << k[n] << std::endl;
     //    std::cerr << "n: " << n << " E[n]: " << E[n] << std::endl;
     //    std::cerr << "n: " << n << " dk[n]: " << dk[n] << std::endl;
     //    std::cerr << "n: " << n << " wn[n]: " << wn[n] << std::endl;
     //    std::cerr << std::endl;
     // }

      for(int n=0; n<Nmodes; ++n)
      {
        for(int m=0; m<Nwaves[n]; ++m)
        {
          // knm = unit vector * magnitude
          for(int i=0; i<3; ++i)
            knm[i][n][m] = enm[i][n][m] * k[n];

          // init random coefficients
          for(int i=0; i<3; ++i)
          {
            Anm[i][n][m] = 0;
            Bnm[i][n][m] = 0;
          }
        }
      }
    }

    template <int dim>
    real_t calculate_velocity_dir(const real_t x[3])
    {
      real_t u = 0;
      static thread_local constexpr int dim1 = (dim + 1) % 3;
      static thread_local constexpr int dim2 = (dim + 2) % 3;

      for(int n=0; n<Nmodes; ++n)
      {
        for(int m=0; m<Nwaves[n]; ++m)
        {
          const real_t r = (knm[0][n][m] * x[0] + knm[1][n][m] * x[1] + knm[2][n][m] * x[2]);
          u += (Anm[dim1][n][m] * enm[dim2][n][m] - Anm[dim2][n][m] * enm[dim1][n][m])*cos(r) - (Bnm[dim1][n][m] * enm[dim2][n][m] - Bnm[dim2][n][m] * enm[dim1][n][m])*sin(r); 
        }
      }
      return u;
    };

    void calculate_velocity(real_t &u, real_t &v, real_t &w, const real_t x[3])
    {
      u = 0;
      v = 0;
      w = 0;

      for(int n=0; n<Nmodes; ++n)
      {
        //// std::cerr << i << " " << j << " " << k << " wnt: " << wnt << std::endl;
        for(int m=0; m<Nwaves[n]; ++m)
        {
//          std::cerr << "mode_idx: " << n << " wave_idx: " << m << " knm[0]: " << knm[0][n][m] << " knm[1]: " << knm[1][n][m] << " knm[2]: " << knm[2][n][m] << std::endl;

          const real_t r = (knm[0][n][m] * x[0] + knm[1][n][m] * x[1] + knm[2][n][m] * x[2]);

  //        std::cerr << "mode_idx: " << n << " wave_idx: " << m << " cos coeff: " << (Anm[1][n][m] * e[2] - Anm[2][n][m] * e[1]) << " sin coeff: " << (Bnm[1][n][m] * e[2] - Bnm[2][n][m] * e[1]) << " Anm[1][n][m]: " << Anm[1][n][m] << " e[2]: std::endl; 

          u += (Anm[1][n][m] * enm[2][n][m] - Anm[2][n][m] * enm[1][n][m])*cos(r) - (Bnm[1][n][m] * enm[2][n][m] - Bnm[2][n][m] * enm[1][n][m])*sin(r); 
          v += (Anm[2][n][m] * enm[0][n][m] - Anm[0][n][m] * enm[2][n][m])*cos(r) - (Bnm[2][n][m] * enm[0][n][m] - Bnm[0][n][m] * enm[2][n][m])*sin(r); 
          w += (Anm[0][n][m] * enm[1][n][m] - Anm[1][n][m] * enm[0][n][m])*cos(r) - (Bnm[0][n][m] * enm[1][n][m] - Bnm[1][n][m] * enm[0][n][m])*sin(r); 
        }
      }
    };

    template<int nx>
    void generate_velocity_field(real_t u[nx][nx][nx], real_t v[nx][nx][nx], real_t w[nx][nx][nx], const real_t &dx, const real_t &t)
    {
      update_time(t);

      #pragma omp parallel for
      for(int i=0; i<nx; ++i)
        for(int j=0; j<nx; ++j)
          for(int k=0; k<nx; ++k)
          {
            const real_t x[3] = {i*dx, j*dx, k*dx};
            calculate_velocity(u[i][j][k], v[i][j][k], w[i][j][k], x);
          }
    };
  };
};
