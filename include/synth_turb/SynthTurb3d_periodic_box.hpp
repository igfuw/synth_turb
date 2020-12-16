#include "SynthTurb3d_common.hpp"

namespace SynthTurb
{
  template<class real_t, int Nmodes, int Nwaves>
  class SynthTurb3d_periodic_box : public SynthTurb3d_common<real_t, Nmodes, Nwaves>
  {
    using parent_t = SynthTurb3d_common<real_t, Nmodes, Nwaves>;

    void generate_wavenumbers(const real_t &Lmax, const real_t &Lmin) override
    {
      // wavenumbers in the form k = n * 2 PI / L, where n=1,2,3,...,Nmodes to get periodic flow 
      for(int n=0; n<Nmodes; ++n)
        this->k[n] = (n+1) * (2. * M_PI / Lmax);
    }

    void generate_unit_wavevectors(const int &m) override
    {
      // generate unit vector along xyz
      if(Nwaves!=3) throw std::runtime_error("Nwaves needs to be 3 for periodic flow");
      this->e[0]=0;
      this->e[1]=0;
      this->e[2]=0;
      this->e[m]=1;
    }

    public:

    //ctor
    SynthTurb3d_periodic_box(
      const real_t &eps,        // TKE dissipation rate [m2/s3]
      const real_t &Lmax = 100, // maximum length scale [m]
      const real_t &Lmin = 1e-3 // Kolmogorov length scale [m]
    )
    {
      this->init(eps, Lmax, Lmin);
    }
  };
};
