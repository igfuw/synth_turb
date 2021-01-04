import numpy as np
import matplotlib.pyplot as plt

#datadir = '/home/piotr/praca/coal_fluctu_dim/synth_turb_tests/separation_data/DT0.1_L1_EPS0.1/'
datadir = 'build/'

for sim_type in ["all_waves", "periodic_box", "GA17"]:
  time = np.loadtxt(str(datadir)+'pair_separation_'+str(sim_type)+'.dat', usecols=0)
  mean_r = np.loadtxt(str(datadir)+'pair_separation_'+str(sim_type)+'.dat', usecols=1)
  sig_r = np.loadtxt(str(datadir)+'pair_separation_'+str(sim_type)+'.dat', usecols=2)
  plt.plot(time, mean_r, label=sim_type)
  plt.fill_between(time, mean_r-sig_r, mean_r+sig_r, alpha=0.1)

plt.xlabel('time [s]')
plt.ylabel('pair separation [m] (mean+-std dev)')
plt.legend()
plt.show()
