import pandas as pd
import numpy as np
import os

targetFile1 = "00_1-bromobutane_RCE.target"
targetFile2 = "00_2-bromobutane_RCE.target"

energyFile1OPLS = os.path.join("1-bromobutane_OPLS", "output","Energy.rel.kcal.txt")
energyFile1octOpt = os.path.join("1-bromobutane_octOpt", "output","Energy.rel.kcal.txt")
energyFile2OPLS = os.path.join("2-bromobutane_OPLS", "output","Energy.rel.kcal.txt")
energyFile2octOpt = os.path.join("2-bromobutane_octOpt", "output","Energy.rel.kcal.txt")

def  calcMAPE(targets, sims):
  #print(f'{targets}')
  #print(f'{sims}')
  if not len(targets) == len(sims):
    print(f'len(targets) ({len(targets)}) != len(sims) ({len(sims)}). Exit.')
    exit()
  mape = 0
  for i in range(len(targets)):
    if targets[i] == 0 and sims[i] == 0:
      pass
    elif targets[i] == 0 and not sims[i] == 0:
      print(f'targets[i] == 0, but sims[i] != 0. Add 1e-6 to targets.')
      mape = mape + np.absolute((1e-6 - sims[i])/1e-6)
    else:
      mape = mape + np.absolute((targets[i] - sims[i])/targets[i])
  return (mape/len(targets))*100


def csv2arr(filename):
  try:
    df = pd.read_csv(filename, header=None, delimiter=" ")
  except:
    print(f'could not save file \"{filename}\" to a pandas df. Exit.')
    exit()
  df = df.drop(df.columns[0], axis=1)  
  arr = df.to_numpy()
  return arr

if __name__ == "__main__":
  targets1 = csv2arr(targetFile1)
  targets2 = csv2arr(targetFile2)
  
  energies1OPLS = csv2arr(energyFile1OPLS)
  energies1octOpt = csv2arr(energyFile1octOpt)
  energies2OPLS = csv2arr(energyFile2OPLS)
  energies2octOpt = csv2arr(energyFile2octOpt)
  
  mape1OPLS = calcMAPE(targets1, energies1OPLS)
  print(f'1-bromobutane OPLS:   {mape1OPLS}')
  mape1octOpt = calcMAPE(targets1, energies1octOpt)
  print(f'1-bromobutane octOpt: {mape1octOpt}\n')
  mape2OPLS = calcMAPE(targets2, energies2OPLS)
  print(f'2-bromobutane OPLS:   {mape2OPLS}')
  mape2octOpt = calcMAPE(targets2, energies2octOpt)
  print(f'2-bromobutane octOpt: {mape2octOpt}')
