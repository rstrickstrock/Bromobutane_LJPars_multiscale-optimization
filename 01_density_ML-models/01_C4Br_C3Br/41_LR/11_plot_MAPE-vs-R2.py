import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
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
  
figDPI = 100

xLabel = "MAPE"
yLabel = "R²"

facecolor = "#005AB5"
minmaxvals = True

if minmaxvals:
  Filter = True
else:
  Filter = False

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
  dfStatistics = pd.read_csv(statisticsFile)
  #print(f'{dfStatistics}')
  try:
    dfStatistics = dfStatistics.drop(columns=["Unnamed: 0"])
  except:
    print(f'Something went wrong with\'dfStatistics = dfStatistics.drop(columns=["Unnamed: 0"])\'.')
  else:
    #print(f'{dfStatistics}')
    pass

if Filter:
  # Data filtering
  thisDF = dfStatistics[
      (dfStatistics["mape"] >= minMAPE) &
      (dfStatistics["mape"] <= maxMAPE) &
      (dfStatistics["r2"] >= minR2) &
      (dfStatistics["r2"] <= maxR2)
  ]
else:
  thisDF = dfStatistics


g = sns.jointplot(data=thisDF, x="mape", y="r2", color=facecolor)
# Clear default histograms
g.ax_marg_x.cla()
g.ax_marg_y.cla()
# X marginal KDE
sns.kdeplot(
    data=thisDF,
    x="mape",
    ax=g.ax_marg_x,
    fill=True,           # optional: fill under curve
    color=facecolor
)
# Y marginal KDE
sns.kdeplot(
    data=thisDF,
    y="r2",
    ax=g.ax_marg_y,
    fill=True,
    color=facecolor
)

g.ax_joint.set_xlabel('MAPE [%]', fontweight='bold', fontsize=20)
g.ax_joint.set_ylabel('R²', fontweight='bold', fontsize=20)
g.ax_joint.grid(color='lightgray', linestyle='dotted')
g.ax_marg_x.axis("off")
g.ax_marg_y.axis("off")

if minmaxvals:
  g.ax_joint.set_xlim([minMAPE-xAxisTol, maxMAPE+xAxisTol])
  g.ax_joint.set_ylim([minR2-yAxisTol, maxR2+yAxisTol])
  g.ax_joint.set_xticks(MAPETicks)
  g.ax_joint.set_xticklabels(MAPETickLabels, fontsize=15)
  g.ax_joint.set_yticks(R2Ticks)
  g.ax_joint.set_yticklabels(R2TickLabels, fontsize=15)


 
plt.tight_layout()
# Show or save
if saveOrShow == "show":
    plt.show()
if saveOrShow == "save":
    plt.savefig(f'MAPE-vs-R2_LR_{datasetname}_n{len(thisDF)}.png', dpi=figDPI)



