import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

pwd = "/home/rstric2s/current_sim/Application/20_1-bromobutane_opt/maxR2/"
topValues = 3
diffType = "relTargSimDiff"

try:
  saveOrShow = sys.argv[1]
except:
  saveOrShow = "show"
if saveOrShow == "save":
  pass
else:
  saveOrShow = "show"

if diffType == "relTargSimDiff":
  cbarLabel = f'rel. diff Target vs. Simulation'
elif diffType == "absTargSimDiff":
  cbarLabel = f'abs. diff Target vs. Simulation'
elif diffType == "relPredSimDiff":
  cbarLabel = f'rel. diff Prediction vs. Simulation'
elif diffType == "absPredSimDiff":
  cbarLabel = f'abs. diff Prediction vs. Simulation'
else:
  print(f'diffType: {diffType} not valid. Must be one of "relTargSimDiff", "absTargSimDiff", "relPredSimDiff", "absPredSimDiff". Exit.')
  exit()


thisDir = os.path.basename(pwd)
#print(f'{thisDir}')

statsFile = os.path.join(pwd, "StatsDensity_withSimResults.csv")
stats = pd.read_csv(statsFile)
stats = stats.drop(stats.columns[0], axis=1)
#print(f'{stats}')

SigC = stats["SigC136"].to_numpy()
#print(f'{SigC136}')
SigBr = stats["SigBr730"].to_numpy()
EpsC = stats["EpsC136"].to_numpy()
EpsBr = stats["EpsBr730"].to_numpy()

diff = stats[diffType].to_numpy()

sortedDiff = stats[diffType].abs().sort_values()
#print(f'{type(sortedDiff)}')
topSigC = []
topSigBr = []
topEpsC = []
topEpsBr = []
topDiff = []
for i in range(0, topValues):
  #print(f'{sortedDiff.iloc[i]}')
  thisDiff = sortedDiff.iloc[i]
  #print(f'{thisDiff}')
  thisIndex = sortedDiff[sortedDiff == thisDiff].index[0]
  #print(f'{thisIndex}')
  #print(f'{stats["SigC136"][thisIndex]}')
  topSigC.append(stats["SigC136"][thisIndex])
  topSigBr.append(stats["SigBr730"][thisIndex])
  topEpsC.append(stats["EpsC136"][thisIndex])
  topEpsBr.append(stats["EpsBr730"][thisIndex])
  topDiff.append(stats[diffType][thisIndex])

gs_kw = dict(width_ratios=[1, 1, 1], height_ratios=[1, 1])
fig, axd = plt.subplot_mosaic([['SigCvsSigBr', 'SigCvsEpsC', 'SigCvsEpsBr'], 
                               ['SigBrvsEpsC', 'SigBrvsEpsBr', 'EpsCvsEpsBr']], 
                               gridspec_kw=gs_kw, figsize=(18.0, 12.0))

p1 = axd["SigCvsSigBr"].scatter(SigC, SigBr, c=diff, marker='o')#, edgecolor="#44AA99")
axd["SigCvsSigBr"].set(xlabel=f'SigC', ylabel=f'SigBr')
axd["SigCvsSigBr"].set_title("SigC vs SigBr", fontweight='bold')
fig.colorbar(p1, ax=axd["SigCvsSigBr"], label=f'{cbarLabel}')
axd["SigCvsSigBr"].scatter(topSigC, topSigBr, label='top 3 results', marker='x', color="#FF0000")
axd["SigCvsSigBr"].legend()

p2 = axd["SigCvsEpsC"].scatter(SigC, EpsC, c=diff, marker='o')#, edgecolor="#44AA99")
axd["SigCvsEpsC"].set(xlabel=f'SigC', ylabel=f'EpsC')
axd["SigCvsEpsC"].set_title("SigC vs EpsC", fontweight='bold')
fig.colorbar(p2, ax=axd["SigCvsEpsC"], label=f'{cbarLabel}')
axd["SigCvsEpsC"].scatter(topSigC, topEpsC, label='top 3 results', marker='x', color="#FF0000")
axd["SigCvsEpsC"].legend()

p3 = axd["SigCvsEpsBr"].scatter(SigC, EpsBr, c=diff, marker='o')#, edgecolor="#44AA99")
axd["SigCvsEpsBr"].set(xlabel=f'SigC', ylabel=f'EpsBr')
axd["SigCvsEpsBr"].set_title("SigC vs EpsBr", fontweight='bold')
fig.colorbar(p3, ax=axd["SigCvsEpsBr"], label=f'{cbarLabel}')
axd["SigCvsEpsBr"].scatter(topSigC, topEpsBr, label='top 3 results', marker='x', color="#FF0000")
axd["SigCvsEpsBr"].legend()

p4 = axd["SigBrvsEpsC"].scatter(SigBr, EpsC, c=diff, marker='o')#, edgecolor="#44AA99")
axd["SigBrvsEpsC"].set(xlabel=f'SigBr', ylabel=f'EpsC')
axd["SigBrvsEpsC"].set_title("SigBr vs EpsC", fontweight='bold')
fig.colorbar(p4, ax=axd["SigBrvsEpsC"], label=f'{cbarLabel}')
axd["SigBrvsEpsC"].scatter(topSigBr, topEpsC, label='top 3 results', marker='x', color="#FF0000")
axd["SigBrvsEpsC"].legend()

p5 = axd["SigBrvsEpsBr"].scatter(SigBr, EpsBr, c=diff, marker='o')#, edgecolor="#44AA99")
axd["SigBrvsEpsBr"].set(xlabel=f'SigBr', ylabel=f'EpsBr')
axd["SigBrvsEpsBr"].set_title("SigBr vs EpsBr", fontweight='bold')
fig.colorbar(p5, ax=axd["SigBrvsEpsBr"], label=f'{cbarLabel}')
axd["SigBrvsEpsBr"].scatter(topSigBr, topEpsBr, label='top 3 results', marker='x', color="#FF0000")
axd["SigBrvsEpsBr"].legend()

p6 = axd["EpsCvsEpsBr"].scatter(EpsC, EpsBr, c=diff, marker='o')#, edgecolor="#44AA99")
axd["EpsCvsEpsBr"].set(xlabel=f'EpsC', ylabel=f'EpsBr')
axd["EpsCvsEpsBr"].set_title("EpsC vs EpsBr", fontweight='bold')
fig.colorbar(p6, ax=axd["EpsCvsEpsBr"], label=f'{cbarLabel}')
axd["EpsCvsEpsBr"].scatter(topEpsC, topEpsBr, label='top 3 results', marker='x', color="#FF0000")
axd["EpsCvsEpsBr"].legend()

plt.tight_layout()

if saveOrShow == "save":
  plt.savefig(f'optParams-vs-{diffType}_{thisDir}.png', dpi=300, format='png')
  #break
if saveOrShow == "show":
  plt.show()

#TODO:
# write top results in File





