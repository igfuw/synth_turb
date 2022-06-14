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

for eps in [0.01,0.1,1,10,100,1000]:
  plt.clf()
  #for sim_type in ["all_waves", "periodic_box", "periodic_box_multiwave", "GA17"]:
  for sim_type in ["all_waves", "periodic_box"]:
    simname = str(sim_type)+'_EPS'+str(eps)
    time = np.loadtxt(str(datadir)+'pair_separation_'  +simname+'.dat', usecols=0)
    mean_r = np.loadtxt(str(datadir)+'pair_separation_'+simname+'.dat', usecols=1)
    sig_r = np.loadtxt(str(datadir)+'pair_separation_' +simname+'.dat', usecols=2)
  
  #  time /= 0.06976
  #  mean_r /= 1e-3
  #  sig_r /= 1e-3
  
    plt.plot(time, mean_r, label=labels[sim_type])
    plt.fill_between(time, mean_r-sig_r, mean_r+sig_r, alpha=0.1)
  
  plt.xlabel('time [s]')
  #plt.xlim(0,150)
  #plt.ylim(0,1500)
  plt.ylabel('pair separation [m] (mean+-std dev)')
  plt.legend()
  plt.savefig('pair_separation_series_EPS'+str(eps)+'.png')
