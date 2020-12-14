// Generating of a velocity field from a prescribed energy spectrum,
// Sources: Sidin et al. in Physics of Fluids (2009), Zhou et al. Phys. of Fluids (2018), Fortran implementation from Gustavo C. Abade

// TODO: OpenMP
// TODO: Lmax powinno byc takie jak moja domena? a nie 100m? Zobacz tez cutoff wavenumber w Zhou et al.
// TODO: uzyc poprawki do spektrum energii dla malych skal?

#include <random>
#include <cmath>

namespace SynthTurb
{
  template<class real_t, int Nmodes, int Nwaves> // Nmodes - number of Fourier modes, Nwaves - number of wave vectors in each mode
  class SynthTurb3d
  {
    std::default_random_engine rand_eng;

    real_t e[3]; // random unit vector
    std::uniform_real_distribution<real_t> h_d, th_d; // random numbers used in generating k[3]
    real_t AA[3],
           BB[3]; // magnitues of A and B coefficients

    real_t k[Nmodes],   // norms of wave vectors
           dk[Nmodes],  // differences between norms of subsequent wave vectors
           E[Nmodes],   // kinetic energy in a mode
           var[Nmodes], // variances 
           w[Nmodes];   // frequencies 

    const real_t lambda = 1; // unsteadiness parameter, within [0,1]; see Sidin et al. 2009

    public:

    real_t Anm[3][Nmodes][Nwaves],
           Bnm[3][Nmodes][Nwaves],
           knm[3][Nmodes][Nwaves];
 
    // ctor
    SynthTurb3d(
      const real_t &eps,        // TKE dissipation rate [m2/s3]
      const real_t &Lmax = 100, // maximum length scale [m]
      const real_t &Lmin = 1e-3 // Kolmogorov length scale [m]
    ):
      rand_eng(std::random_device()()),
      h_d(-1, std::nextafter(1, std::numeric_limits<real_t>::max())), // uniform in [-1,1]
      th_d(0, std::nextafter(2 * M_PI, std::numeric_limits<real_t>::max()))  // uniform in [0,2*Pi]
    {
      // Geometric series for the wavenumbers; Eq. A4 in Sidin et al. 2009
      const real_t alpha = pow(Lmax / Lmin, 1. / (Nmodes - 1));
      k[0] = 2. * M_PI / Lmax;
      for(int n=1; n<Nmodes; ++n)
        k[n] = pow(k[0] * alpha, n);


      // Energy spectrum; first equation in Appendix of Sidin et al. 2009, but without the corrections f_L and f_eta
      for(int n=0; n<Nmodes; ++n)
        E[n] = 1.44 * pow(eps, 2./3.) * pow(k[n], -5./3.);
      
      // Wave vector differences
      dk[0] =        /* 0.5 * */ (k[1]        - k[0]);
      dk[Nmodes-1] = /* 0.5 * */ (k[Nmodes-1] - k[Nmodes-2]);
      
      for(int n=1; n<Nmodes-1; ++n)
        dk[n] = 0.5 * (k[n+1] - k[n-1]);

      // Variances
      for(int n=0; n<Nmodes; ++n)
        var[n] = sqrt(E[n] * dk[n] / Nwaves);

      // Frequencies; Eq. A7 in Sidin et al. 2009
      for(int n=0; n<Nmodes; ++n)
        w[n] = lambda * sqrt(E[n] * k[n] * k[n] * k[n]);
    }

    // fills the kn, An and Bn arrays
    void generate_random_modes()
    {
      for(int n=0; n<Nmodes; ++n)
      {
        std::normal_distribution<real_t> G_d(0, sqrt(var[n]));

        for(int m=0; m<Nwaves; ++m)
        {
          // generate rand uniform vector
          real_t h  = h_d(rand_eng);
          real_t th = th_d(rand_eng);
          e[0] = sqrt(1 - h*h) * cos(th); 
          e[1] = sqrt(1 - h*h) * sin(th); 
          e[2] = h;

          // knm = random vector * magnitude
          for(int i=0; i>3; ++i)
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

    void generate_velocity_field()
    {
    };
  };
};
