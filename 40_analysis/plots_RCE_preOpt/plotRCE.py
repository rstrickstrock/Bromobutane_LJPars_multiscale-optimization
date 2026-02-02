import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from pandas.plotting import parallel_coordinates
import plotly.graph_objects as go
from matplotlib.lines import Line2D
import sys

figDPI = 100

RCE1TargetFile = '01_1-bromobutane_RCE.target'
RCE2TargetFile = '01_2-bromobutane_RCE.target'

RCE1OctOptFile = '01_1-bromobutane_RCE.octOpt'
RCE2OctOptFile = '01_2-bromobutane_RCE.octOpt'
RCE1OPLSFile = '01_1-bromobutane_RCE.OPLS'
RCE2OPLSFile = '01_2-bromobutane_RCE.OPLS'

try:
  saveOrShow = sys.argv[1]
except:
  saveOrShow = "show"
if saveOrShow == "save":
  pass
else:
  saveOrShow = "show"
  
  
def readRCEFile(pathToFile):
  try:
    dfRCEs = pd.read_csv(pathToFile)
  except:
    print(f'Could not read csv file. Exit.')
    exit()
  else:
    #print(f'{dfRCEs}')
    try:
      dfRCEs = dfRCEs.drop(columns=["Unnamed: 0"])
    except:
      pass
    else:
      #print(f'{dfRCEs}')
      pass
  
  return dfRCEs.RCE.to_numpy()


if __name__ == "__main__":
  RCE1Targets = readRCEFile(RCE1TargetFile)
  idx1BB = np.argsort(RCE1Targets)
  #print(f'{RCE1Targets}')
  RCE2Targets = readRCEFile(RCE2TargetFile)
  idx2BB = np.argsort(RCE2Targets)
  
  RCE1OctOpt = readRCEFile(RCE1OctOptFile)
  RCE2OctOpt = readRCEFile(RCE2OctOptFile)
  RCE1OPLS = readRCEFile(RCE1OPLSFile)
  RCE2OPLS = readRCEFile(RCE2OPLSFile)
  
  srtTargets1BB = []
  srtTargets2BB = []
  srtOctOpt1BB = []
  srtOctOpt2BB = []
  srtOPLS1BB = []
  srtOPLS2BB = []
  
  for idx in idx1BB:
    srtTargets1BB.append(RCE1Targets[idx])
    srtOctOpt1BB.append(RCE1OctOpt[idx])
    srtOPLS1BB.append(RCE1OPLS[idx])
  for idx in idx2BB:
    srtTargets2BB.append(RCE2Targets[idx])
    srtOctOpt2BB.append(RCE2OctOpt[idx])
    srtOPLS2BB.append(RCE2OPLS[idx])
  
  #print(f'{np.min([np.min(srtTargets1BB), np.min(srtOctOpt1BB), np.min(srtOPLS1BB)])}')
  minRCE = np.min([np.min(srtTargets1BB), np.min(srtOctOpt1BB), np.min(srtOPLS1BB)])
  minRCE = minRCE - abs(minRCE*0.05)
  maxRCE = np.max([np.max(srtTargets1BB), np.max(srtOctOpt1BB), np.max(srtOPLS1BB)])
  maxRCE = maxRCE + abs(maxRCE*0.05)
  plt.figure(figsize=(5, 5))
  plt.plot([minRCE, maxRCE], [minRCE, maxRCE], '-', color='#000000')
  plt.plot(srtTargets1BB, srtTargets1BB, '-X', color='#000000', label="Targets")
  plt.plot(srtTargets1BB, srtOctOpt1BB, 'o', ls='dotted', color='#FF0000', label='Octane Opt')
  plt.plot(srtTargets1BB, srtOPLS1BB, '^', ls='dashed', color='#8B0000', label='OPLS')
  plt.xlabel('target RCE [kJ/mol]', fontsize=15, fontweight='bold')
  plt.ylabel('reproduced RCE [kJ/mol]', fontsize=15, fontweight='bold')
  ticks = [-2.0, -1.0, 0.0, 1.0, 2.0]
  plt.xticks(ticks, ["-2.0", "-1.0", "0.0", "1.0", "2.0"], fontsize=13)
  plt.yticks(ticks, ["-2.0", "-1.0", "0.0", "1.0", "2.0"], fontsize=13)
  plt.xlim([minRCE, maxRCE])
  plt.ylim([minRCE, maxRCE])
  plt.legend()
  plt.tight_layout()
  if saveOrShow == 'save':
    plt.savefig(f'RCE_targets-vs-preOpts_1BB.png', dpi=100)
  
  minRCE = np.min([np.min(srtTargets2BB), np.min(srtOctOpt2BB), np.min(srtOPLS2BB)])
  minRCE = minRCE - abs(minRCE*0.05)
  maxRCE = np.max([np.max(srtTargets2BB), np.max(srtOctOpt2BB), np.max(srtOPLS2BB)])  
  maxRCE = maxRCE + abs(maxRCE*0.05)
  plt.figure(figsize=(5, 5))
  plt.plot([minRCE, maxRCE], [minRCE, maxRCE], '-', color='#000000')
  plt.plot(srtTargets2BB, srtTargets2BB, '-X', color='#000000', label="Targets")
  plt.plot(srtTargets2BB, srtOctOpt2BB, 'o', ls='dotted', color='#FF0000', label='Octane Opt')
  plt.plot(srtTargets2BB, srtOPLS2BB, '^', ls='dashed', color='#8B0000', label='OPLS')
  plt.xlabel('target RCE [kJ/mol]', fontsize=15, fontweight='bold')
  plt.ylabel('reproduced RCE [kJ/mol]', fontsize=15, fontweight='bold')
  ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
  plt.xticks(ticks, fontsize=13)
  plt.yticks(ticks, fontsize=13)
  plt.xlim([minRCE, maxRCE])
  plt.ylim([minRCE, maxRCE])
  plt.legend()
  plt.tight_layout()
  if saveOrShow == 'save':
    plt.savefig(f'RCE_targets-vs-preOpts_2BB.png', dpi=100)  
  else:
    plt.show()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
