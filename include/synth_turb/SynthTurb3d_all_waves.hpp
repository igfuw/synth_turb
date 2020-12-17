#include "SynthTurb3d_common.hpp"

namespace SynthTurb
{
  template<class real_t, int Nmodes, int Nwaves>
  class SynthTurb3d_all_waves : public SynthTurb3d_common<real_t, Nmodes, Nwaves>
  {
    using parent_t = SynthTurb3d_common<real_t, Nmodes, Nwaves>;

    void generate_wavenumbers(const real_t &Lmax, const real_t &Lmin) override
    {
      // Geometric series for the wavenumbers; Eq. A4 in Sidin et al. 2009
      const real_t alpha = pow(Lmax / Lmin, 1. / (Nmodes - 1));
      // std::cerr << "alpha: " << alpha << std::endl;
      this->k[0] = 2. * M_PI / Lmax;
      for(int n=1; n<Nmodes; ++n)
        this->k[n] = this->k[0] * pow(alpha, n);
    }

    void generate_unit_wavevectors(const int &m) override
    {
      // generate random unit vector
      real_t h  = this->h_d(this->rand_eng);
      real_t th = this->th_d(this->rand_eng);
      this->e[0] = sqrt(1. - h*h) * cos(th); 
      this->e[1] = sqrt(1. - h*h) * sin(th); 
      this->e[2] = h;
    }

    public:

    //ctor
    SynthTurb3d_all_waves(
      const real_t &eps,        // TKE dissipation rate [m2/s3]
      const real_t &Lmax = 100, // maximum length scale [m]
      const real_t &Lmin = 1e-3 // Kolmogorov length scale [m]
    )
    {
      this->init(eps, Lmax, Lmin);
    }
  };
};
