import numpy as np
import matplotlib.pyplot as plt

for sim_type in ["all_waves", "periodic_box"]:
  time = np.loadtxt('build/pair_separation_'+str(sim_type)+'.dat', usecols=0)
  mean_r = np.loadtxt('build/pair_separation_'+str(sim_type)+'.dat', usecols=1)
  sig_r = np.loadtxt('build/pair_separation_'+str(sim_type)+'.dat', usecols=2)
  plt.errorbar(time, mean_r, sig_r, label=sim_type)

plt.xlabel('time [s]')
plt.ylabel('mean separation [m]')
plt.legend()
plt.show()
