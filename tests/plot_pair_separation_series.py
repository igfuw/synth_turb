import numpy as np
import matplotlib.pyplot as plt

#datadir = '/home/piotr/praca/coal_fluctu_dim/synth_turb_tests/separation_data/DT0.1_L1_EPS0.1/'
datadir = '/home/piotr/praca/coal_fluctu_dim/synth_turb_tests/separation_data/synth_turb_separate_code/Feb_2022/multiwave_nwave50_nmode10_and_others/'
#datadir = 'build/'

labels = {
  "all_waves" : "non-periodic",
  "periodic_box" : "periodic",
  "periodic_box_multiwave" : "periodic"
}

fig, axs = plt.subplots(3,1, figsize=(6.3,6.3))

#for eps in [0.01,0.1,1,10,100,1000]:
for i, eps in enumerate([1,10,100]):
#  plt.clf()
  #for sim_type in ["all_waves", "periodic_box", "periodic_box_multiwave", "GA17"]:
  for sim_type in ["all_waves", "periodic_box"]:
    simname = str(sim_type)+'_EPS'+str(eps)
    time = np.loadtxt(str(datadir)+'pair_separation_'  +simname+'.dat', usecols=0)
    mean_r = np.loadtxt(str(datadir)+'pair_separation_'+simname+'.dat', usecols=1)
    sig_r = np.loadtxt(str(datadir)+'pair_separation_' +simname+'.dat', usecols=2)
  
  #  time /= 0.06976
  #  mean_r /= 1e-3
  #  sig_r /= 1e-3
  
    axs[i].plot(time, mean_r, label=labels[sim_type])
    axs[i].fill_between(time, mean_r-sig_r, mean_r+sig_r, alpha=0.1)
    axs[i].set_title('$\epsilon = '+str(eps)+'\mathrm{cm}^2 \mathrm{s}^{-3}$')
  
  #plt.xlim(0,150)
  #plt.ylim(0,1500)
  axs[i].set_ylabel('pair separation [m]')

axs[0].tick_params(labelbottom=False)
axs[1].tick_params(labelbottom=False)

axs[2].set_xlabel('time [s]')
axs[0].legend()
plt.tight_layout()
plt.savefig('pair_separation_series.png')
plt.savefig('pair_separation_series.pdf')
