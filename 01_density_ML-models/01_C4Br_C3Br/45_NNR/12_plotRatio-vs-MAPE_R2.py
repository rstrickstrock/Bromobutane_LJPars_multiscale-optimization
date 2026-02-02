import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import os
import glob
import sys

datasetname = "1-Bromobutane"
datasetname = "2-Bromobutane"

if datasetname == "1-Bromobutane":
  statisticsFile = 'Stats_1-bromobutane.csv'
elif datasetname == "2-Bromobutane":
  statisticsFile = 'Stats_2-bromobutane.csv'
else:
  print(f'Check datasetname. Exit.')
  exit()

mapeColor = "#005AB5"
R2Color = "#FFA54A"

figDPI = 100

minmaxvals = True

if minmaxvals:
  minMAPE = 0.00
  maxMAPE = 0.20
  xAxisTol = 0.005
  MAPETicks = [0.0, 0.04, 0.08, 0.12, 0.16, 0.2]
  MAPETickLabels = []
  for tick in MAPETicks:
    MAPETickLabels.append(f'{tick*100:.1f}')
  minR2 = 0.0
  maxR2 = 1.0
  yAxisTol = 0.05
  R2Ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
  R2TickLabels = []
  for tick in R2Ticks:
    R2TickLabels.append(f'{tick:.1f}')

try:
  saveOrShow = sys.argv[1]
except:
  saveOrShow = "show"
if saveOrShow == "save":
  pass
else:
  saveOrShow = "show"

  
if not os.path.isfile(statisticsFile):
  print(f'Can not find and open \'{statisticsFile}\'. Exit.')
  exit()
else:
  thisDF = pd.read_csv(statisticsFile)
  #print(f'{thisDF}')
  try:
    thisDF = thisDF.drop(columns=["Unnamed: 0"])
  except:
    print(f'Something went wrong with\'thisDF = thisDF.drop(columns=["Unnamed: 0"])\'.')
  else:
    #print(f'{thisDF}')
    pass

## plots
datasetSubsets = thisDF["ratio"].unique()
#print(f'{datasetSubsets}')
#print(f'{len(datasetSubsets)}')

MAPE = []
MAPEerr = []

R2 = []
R2err = []

Ratios = []

for nSubset in range(0, len(datasetSubsets)):
  #print(f'nSubset: {nSubset}')
  thisRatio = datasetSubsets[nSubset]
  Ratios.append(thisRatio)
  
  dfThisSubset = thisDF[(thisDF["ratio"] == thisRatio)]
  #print(f'{dfThisSubset}')
  thisMAPE = dfThisSubset["mape"].to_numpy()
  thisR2 = dfThisSubset["r2"].to_numpy()
  #print(f'{thisMAPE}')
  #print(f'{thisR2}')
  
  thisAvgMAPE = np.mean(thisMAPE)
  MAPE.append(thisAvgMAPE)
  thisUpperMAPEerr = np.std(thisMAPE)
  if thisAvgMAPE - thisUpperMAPEerr <= 0.0:
    thisLowerMAPEerr = thisAvgMAPE
  else:
    thisLowerMAPEerr = thisUpperMAPEerr
  MAPEerr.append([thisLowerMAPEerr, thisUpperMAPEerr])
    
  thisAvgR2 = np.mean(thisR2)
  R2.append(thisAvgR2)
  thisLowerR2err = np.std(thisR2)
  if thisAvgR2 + thisLowerR2err <= 1.0:
    thisUpperR2err = thisLowerR2err
  else:
    thisUpperR2err = 1.0 - thisAvgR2
  R2err.append([thisLowerR2err, thisUpperR2err])
  #break
  
#print(f'{np.array(MAPE)}')
#print(f'{np.array(MAPEerr)}')
#print(f'{R2}')
#print(f'{R2err}')
#print(f'{Ratios}')
#exit()
fig, ax1 = plt.subplots()

ax1.errorbar(np.array(Ratios), np.array(MAPE), yerr=np.array(MAPEerr).T, color=mapeColor, ls='dotted', capsize=4.0, marker='x')
ax1.set_xlabel("Splitting Ratio [%]", fontsize=20, fontweight='bold')
ax1.set_ylabel("MAPE [%]", fontsize=20, fontweight='bold')
#ax1.set_ylabel("MAPE [%]", color=mapeColor, fontsize=14, fontweight='bold')
#ax1.tick_params(axis='y', labelcolor=mapeColor)

ax1.set_xticks([0.2, 0.4, 0.6, 0.8])
ax1.set_xticklabels(["20", "40", "60", "80"], fontsize=15)
plt.grid(color='lightgray', linestyle='dotted')

if minmaxvals:
  ax1.set_ylim([np.array(MAPETicks).min()-xAxisTol, np.array(MAPETicks).max()+xAxisTol])
  ax1.set_yticks(MAPETicks)
  ax1.set_yticklabels(MAPETickLabels, fontsize=15)

plt.tight_layout()
# Show or save
if saveOrShow == "show":
    plt.show()
if saveOrShow == "save":
    plt.savefig(f'Ratio-vs-MAPE_NNR_{datasetname}.png', dpi=figDPI)


fig, ax1 = plt.subplots()
ax1.errorbar(np.array(Ratios), np.array(R2), yerr=np.array(R2err).T, color=R2Color, ls='dotted', capsize=4.0, marker='x')
ax1.set_xlabel("Splitting Ratio [%]", fontsize=20, fontweight='bold')
ax1.set_ylabel("R$^2$", fontsize=20, fontweight='bold')
#ax1.set_ylabel("R$^2$", color=R2Color, fontsize=14, fontweight='bold')
#ax1.tick_params(axis='y', labelcolor=R2Color)

ax1.set_xticks([0.2, 0.4, 0.6, 0.8])
ax1.set_xticklabels(["20", "40", "60", "80"], fontsize=15)
plt.grid(color='lightgray', linestyle='dotted')

if minmaxvals:
  ax1.set_ylim([np.array(R2Ticks).min()-yAxisTol, np.array(R2Ticks).max()+yAxisTol])
  ax1.set_yticks(R2Ticks)
  ax1.set_yticklabels(R2TickLabels, fontsize=15)

plt.tight_layout()
# Show or save
if saveOrShow == "show":
    plt.show()
if saveOrShow == "save":
    plt.savefig(f'Ratio-vs-R2_NNR_{datasetname}.png', dpi=figDPI)









