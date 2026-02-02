from scipy.stats import qmc
import pandas as pd
import matplotlib.pyplot as plt

plotthis = True

sampler = qmc.Sobol(d=4)
sample = sampler.random_base2(m=11)
#print(f'{sample}')
#print(f'{qmc.discrepancy(sample)}')

#lower parameter bounds for sig_C800, sig_br730, eps_C800, eps_br730
l_bounds = [0.15, 0.20, 0.01, 1.00]
#upper parameter bounds for sig_C800, sig_br730, eps_C800, eps_br730
u_bounds = [0.75, 0.75, 0.99, 3.00]

sample = qmc.scale(sample, l_bounds, u_bounds)
#print(f'{sample}')

df_sample = pd.DataFrame({'SigC800': sample[:, 0], 'SigBr730': sample[:, 1], 'EpsC800': sample[:, 2], 'EpsBr730': sample[:, 3]})
print(f'{df_sample}')
df_sample.to_csv('parameter_sample.csv')

if plotthis:
  plt.plot(sample[:, 0], sample[:, 1], '.')
  plt.xlabel("SigC800")
  plt.ylabel("SigBr730")
  
  plt.figure()
  plt.plot(sample[:, 0], sample[:, 2], '.')
  plt.xlabel("SigC800")
  plt.ylabel("EpsC800")
  
  plt.figure()
  plt.plot(sample[:, 0], sample[:, 3], '.')
  plt.xlabel("SigC800")
  plt.ylabel("EpsBr730")
  
  plt.figure()
  plt.plot(sample[:, 1], sample[:, 2], '.')
  plt.xlabel("SigBr730")
  plt.ylabel("EpsC800")
  
  plt.figure()
  plt.plot(sample[:, 1], sample[:, 3], '.')
  plt.xlabel("SigBr730")
  plt.ylabel("EpsBr730")
  
  plt.figure()
  plt.plot(sample[:, 2], sample[:, 3], '.')
  plt.xlabel("EpsC800")
  plt.ylabel("EpsBr730")
  
  plt.show()
