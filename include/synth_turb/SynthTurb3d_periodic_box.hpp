#include "SynthTurb3d_common.hpp"

namespace SynthTurb
{
  template<class real_t, int Nmodes, int Nwaves>
  class SynthTurb3d_periodic_box : public SynthTurb3d_common<real_t, Nmodes, Nwaves>
  {
    using parent_t = SynthTurb3d_common<real_t, Nmodes, Nwaves>;

    void generate_wavenumbers(const real_t &Lmax, const real_t &Lmin) override
    {
    //  if(Nmodes!=2) throw std::runtime_error("Nmodes needs to be 2 for periodic flow");

      // wavenumbers in the form k = n * 2 PI / L, where n=1,2,3,...,Nmodes to get periodic flow 
      for(int n=0; n<Nmodes; ++n)
        this->k[n] = (n+1) * (2. * M_PI / Lmax);
//      this->k[0] = (1) * (2. * M_PI / Lmax);
  //    this->k[1] = sqrt(2) * (2. * M_PI / Lmax);
    }

    void generate_unit_wavevectors(const int &wave_idx) override
    {
      // generate unit vector 
      if(Nwaves!=6) throw std::runtime_error("Nwaves needs to be 6 in the non-multiwave periodic box");
   //   if(mode_idx == 0)
   //   {
        this->e[0]=0;
        this->e[1]=0;
        this->e[2]=0;
        this->e[wave_idx%3]= 1 * (int(wave_idx/3) == 0 ? 1 : -1);
   /*   }
      else
      {
        this->e[0]= int(wave_idx/3) == 0 ? 1 : -1;
        this->e[1]= int(wave_idx/3) == 0 ? 1 : -1;
        this->e[2]= int(wave_idx/3) == 0 ? 1 : -1;
        this->e[wave_idx%3]=0;
      }

      std::cerr << "mode_idx: " << mode_idx << " wave_idx: " << wave_idx << " e[0]: " << this->e[0] << " e[1]: " << this->e[1] << " e[2]: " << this->e[2] << std::endl;
      */
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
