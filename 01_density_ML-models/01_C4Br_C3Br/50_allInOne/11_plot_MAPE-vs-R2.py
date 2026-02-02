import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys


datasetname = "1-Bromobutane"
datasetname = "2-Bromobutane"
top = 20

marker = "o"
nbins = 4
xBinwidth = 0.0005
yBinwidth = 0.005
if datasetname == "1-Bromobutane":
  statisticsFiles = [f'Stats_LR_1-bromobutane_top{top}.csv', f'Stats_PR_1-bromobutane_top{top}.csv', f'Stats_RFR_1-bromobutane_top{top}.csv', f'Stats_GPR_1-bromobutane_top{top}.csv', f'Stats_NNR_1-bromobutane_top{top}.csv']
elif datasetname == "2-Bromobutane":
  statisticsFiles = [f'Stats_LR_2-bromobutane_top{top}.csv', f'Stats_PR_2-bromobutane_top{top}.csv', f'Stats_RFR_2-bromobutane_top{top}.csv', f'Stats_GPR_2-bromobutane_top{top}.csv', f'Stats_NNR_2-bromobutane_top{top}.csv']
else:
  print(f'Check datasetname. Exit.')
  exit()
  
figDPI = 100

xLabel = "MAPE"
yLabel = "R²"

facecolors = ['#f781bf', '#ff7f00', '#984ea3', '#4daf4a', '#377eb8']
minmaxvals = True

if minmaxvals:
  minMAPE = 0.00
  maxMAPE = 0.18
  xAxisTol = 0.001
  MAPETicks = [0.00, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18]
  MAPETickLabels = []
  for tick in MAPETicks:
    MAPETickLabels.append(f'{tick*100:.1f}')
  minR2 = 0.5
  maxR2 = 1.0
  yAxisTol = 0.01
  R2Ticks = [0.5, 0.6, 0.7, 0.80, 0.90, 1.0]
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

first = True
for n in range(0, len(statisticsFiles)):
  if not os.path.isfile(statisticsFiles[n]):
    print(f'Can not find and open \'{statisticsFiles[n]}\'. Exit.')
    exit()
  else:
    thisDF = pd.read_csv(statisticsFiles[n])
    #print(f'{thisDF}')
    try:
      thisDF = thisDF.drop(columns=["Unnamed: 0"])
    except:
      print(f'Something went wrong with\'thisDF = thisDF.drop(columns=["Unnamed: 0"])\'.')
    else:
      #print(f'{thisDF}')
      pass
  thisLabel = statisticsFiles[n].split("_")[1]
  if thisLabel == "LR":
    thisLabel = "linear regression"
  elif thisLabel == "PR":
    thisLabel = "polynomial regression"
  elif thisLabel == "RFR":
    thisLabel = "random forest regression"
  elif thisLabel == "GPR":
    thisLabel = "gaussian process regression"
  elif thisLabel == "NNR":
    thisLabel = "neural network regression"
    
  if first:
    g = sns.jointplot(data=thisDF,
                      x="mape",
                      y="r2",
                      marginal_ticks=True,
                      kind="scatter",
                      color=facecolors[n],
                      joint_kws={"marker": marker,
                                 "label": thisLabel}
                      )

    g.ax_marg_x.cla()
    g.ax_marg_y.cla()
    # X marginal KDE
    sns.kdeplot(
        data=thisDF,
        x="mape",
        ax=g.ax_marg_x,
        fill=True,           # optional: fill under curve
        color=facecolors[n]
    )
    # Y marginal KDE
    sns.kdeplot(
        data=thisDF,
        y="r2",
        ax=g.ax_marg_y,
        fill=True,
        color=facecolors[n]
    )
    first = False
  else:
    sns.scatterplot(data=thisDF,
                    x="mape",
                    y="r2",
                    ax=g.ax_joint,
                    color=facecolors[n],
                    marker=marker,
                    label=thisLabel
                    )
    sns.kdeplot(
        data=thisDF,
        x="mape",
        ax=g.ax_marg_x,
        fill=True,           # optional: fill under curve
        color=facecolors[n]
    )
    # Y marginal KDE
    sns.kdeplot(
        data=thisDF,
        y="r2",
        ax=g.ax_marg_y,
        fill=True,
        color=facecolors[n]
    )

#g.ax_marg_x.set_visible(False)
#g.ax_marg_y.set_visible(False)

g.ax_joint.set_xlabel('MAPE [%]', fontweight='bold', fontsize=18)
g.ax_joint.set_ylabel('R²', fontweight='bold', fontsize=18)
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


g.ax_joint.legend(loc="lower left")
plt.tight_layout()
# Show or save
if saveOrShow == "show":
    plt.show()
if saveOrShow == "save":
    plt.savefig(f'MAPE-vs-R2_allInOne_{datasetname}_top{top}.png', dpi=figDPI)



