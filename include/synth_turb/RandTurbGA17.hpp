// calculating SGS velocity from a Langevin equation as in Grabowski & Abade JAS 2017

#pragma once

#include <random>
#include <cmath>

namespace RandTurb
{
  template<class real_t>
  class RandTurb_GA17
  {
    private:
    const real_t CE=0.845,
                 Ctau=1.5,
                 TKE,
                 tau,
                 std_dev;

    std::default_random_engine rand_eng;
    std::normal_distribution<real_t> normal_dist;

    public:
    //ctor
    RandTurb_GA17(const real_t eps, const real_t Lmax):
      TKE(pow(Lmax*eps/CE, real_t(2./3.))),
      tau(Lmax / pow(real_t(2 * M_PI), real_t(1./3.)) * sqrt(Ctau / TKE)),
      std_dev(sqrt(2./3. * TKE)),
      rand_eng(std::random_device()()),
      normal_dist(0,1)
      {}

    void update_sgs_velocity(real_t &u, const real_t &dt)
    {
      u = u * exp(-dt / tau) + sqrt(1 - exp(-2 * dt / tau)) * std_dev * normal_dist(rand_eng); 
    }
  };
}

