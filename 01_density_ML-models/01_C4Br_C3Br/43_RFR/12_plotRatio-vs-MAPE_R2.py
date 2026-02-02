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
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
markers = ['1', '3', 'x', '|', 'v', '<', '*', 'o', 'd', 's']
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

#print(f'{thisDF}')
#print(f'{thisDF["#trees"].unique()}')
#print(f'{len(thisDF["#trees"].unique())}')
#exit()

## plots
numTrees = thisDF["#trees"].unique()
#print(f'{numTrees}')
#print(f'{len(numTrees)}')
Ratios = thisDF["ratio"].unique()
#print(f'{Ratios}')

MAPE = []
MAPEerr = []
R2 = []
R2err = []
for nTree in numTrees:
  #print(f'nTree: {nTree}')  
  dfThisSubset = thisDF[(thisDF["#trees"] == nTree)]
  #print(f'{dfThisSubset}')
  
  thisAvgMAPE = []
  thisMAPEerr = []
  thisAvgR2 = []
  thisR2err = []
  for thisRatio in Ratios:
    dfThisRatioSubset = dfThisSubset[(dfThisSubset["ratio"] == thisRatio)]
    thisRatioMAPE = dfThisRatioSubset["mape"].to_numpy()
    thisRatioR2 = dfThisRatioSubset["r2"].to_numpy()
    #print(f'{thisRatioMAPE}')
    #print(f'{thisRatioR2}')

    thisRatioAvgMAPE = np.mean(thisRatioMAPE)
    thisAvgMAPE.append(thisRatioAvgMAPE)
    thisRatioUpperMAPEerr = np.std(thisRatioMAPE)
    if thisRatioAvgMAPE - thisRatioUpperMAPEerr <= 0.0:
      thisRatioLowerMAPEerr = thisRatioAvgMAPE
    else:
      thisRatioLowerMAPEerr = thisRatioUpperMAPEerr
    thisMAPEerr.append([thisRatioLowerMAPEerr, thisRatioUpperMAPEerr])
    
    thisRatioAvgR2 = np.mean(thisRatioR2)
    thisAvgR2.append(thisRatioAvgR2)
    thisRatioLowerR2err = np.std(thisRatioR2)
    if thisRatioAvgR2 + thisRatioLowerR2err <= 1.0:
      thisRatioUpperR2err = thisRatioLowerR2err
    else:
      thisRatioUpperR2err = 1.0 - thisRatioAvgR2
    thisR2err.append([thisRatioLowerR2err, thisRatioUpperR2err])
    #break
  #print(f'{thisAvgMAPE}')
  #print(f'{len(thisAvgMAPE)}')
  #print(f'{thisMAPEerr}')
  #print(f'{thisAvgR2}')
  #print(f'{thisR2err}')
  MAPE.append(thisAvgMAPE)
  MAPEerr.append(thisMAPEerr)
  R2.append(thisAvgR2)
  R2err.append(thisR2err)
  
#print(f'{np.array(MAPE)}')
#print(f'{np.array(MAPEerr)}')
#print(f'{R2}')
#print(f'{R2err}')
#print(f'{Ratios}')
#exit()

fig, ax1 = plt.subplots()
for n in range(0,len(numTrees)):
  #print(f'{n}')
  nTree = np.array(numTrees)[n]
  ax1.errorbar(np.array(Ratios), np.array(MAPE[n]), yerr=np.array(MAPEerr[n]).T, label=f'num. of trees: {int(nTree)}', color=colors[n], ls='dotted', capsize=4.0, marker=markers[n])

ax1.set_xlabel("Splitting Ratio [%]", fontsize=20, fontweight='bold')
ax1.set_ylabel("MAPE [%]", fontsize=20, fontweight='bold')
#ax1.set_ylabel("MAPE [%]", color=mapeColor, fontsize=14, fontweight='bold')
#ax1.tick_params(axis='y', labelcolor=mapeColor)
plt.legend(loc="upper left")

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
    plt.savefig(f'Ratio-vs-MAPE_numTrees_{datasetname}.png', dpi=figDPI)


fig, ax1 = plt.subplots()
for n in range(0,len(numTrees)):
  #print(f'{n}')
  nTree = numTrees[n]
  ax1.errorbar(np.array(Ratios), np.array(R2[n]), yerr=np.array(R2err[n]).T, label=f'num. of trees: {int(nTree)}', color=colors[n], ls='dotted', capsize=4.0, marker=markers[n])

ax1.set_xlabel("Splitting Ratio [%]", fontsize=20, fontweight='bold')
ax1.set_ylabel("R$^2$", fontsize=20, fontweight='bold')
#ax1.set_ylabel("R$^2$", color=R2Color, fontsize=14, fontweight='bold')
#ax1.tick_params(axis='y', labelcolor=R2Color)
plt.legend(loc="lower left")

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
    plt.savefig(f'Ratio-vs-R2_numTrees_{datasetname}.png', dpi=figDPI)









