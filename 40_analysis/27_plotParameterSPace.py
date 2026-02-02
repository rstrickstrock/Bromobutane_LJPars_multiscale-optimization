import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from pandas.plotting import parallel_coordinates
import plotly.graph_objects as go
from matplotlib.lines import Line2D
import sys

DensFile1BB = "20_1-bromobutane_density_trainingsdata.csv"
DensFile2BB = "20_2-bromobutane_density_trainingsdata.csv"

RCEFile1BB = "20_1-bromobutane_RCE_trainingsdata.csv"
RCEFile2BB = "20_2-bromobutane_RCE_trainingsdata.csv"

targetDensFile = "20_bromobutane_density.target"
targetRCEFile = "20_bromobutane_RCE.target"

figDPI = 100

try:
  saveOrShow = sys.argv[1]
except:
  saveOrShow = "show"
if saveOrShow == "save":
  pass
else:
  saveOrShow = "show"


def mkPlotPC(df1, df2, prop1, prop2, maxMAPE1, maxMAPE2):
  colors = []
  if prop1 == 'Dens1BB':
    thisMAPE1Type = 'mape1BBDens'
    thisLBL1 = r'MAPE($\rho_{(\mathrm{1BB})}$) $\leq$ ' + f'{maxMAPE1} [%]'
    colors.append('#8E44AD')
  elif prop1 == 'Dens2BB':
    thisMAPE1Type = 'mape2BBDens'
    thisLBL1 = r'MAPE($\rho_{(\mathrm{2BB})}$) $\leq$ ' + f'{maxMAPE1} [%]'
    colors.append('#2ECC71')
  elif prop1 == 'RCE1BB':
    thisMAPE1Type = 'mape1BBRCE'
    thisLBL1 = r'MAPE(RCE$_{(\mathrm{1BB})}$) $\leq$ ' + f'{maxMAPE1} [%]'
    colors.append('#D55E00')
  elif prop1 == 'RCE2BB':
    thisMAPE1Type = 'mape2BBRCE'
    thisLBL1 = r'MAPE(RCE$_{(\mathrm{2BB})}$) $\leq$ ' + f'{maxMAPE1} [%]'
    colors.append('#1B4F72')
  else:
    print(f'check prop1. Exit()')
    exit()
    
  if prop2 == 'Dens1BB':
    thisMAPE2Type = 'mape1BBDens'
    thisLBL2 = r'MAPE($\rho_{(\mathrm{1BB})}$) $\leq$ ' + f'{maxMAPE2} [%]'
    colors.append('#8E44AD')
  elif prop2 == 'Dens2BB':
    thisMAPE2Type = 'mape2BBDens'
    thisLBL2 = r'MAPE($\rho_{(\mathrm{2BB})}$) $\leq$ ' + f'{maxMAPE2} [%]'
    colors.append('#2ECC71')
  elif prop2 == 'RCE1BB':
    thisMAPE2Type = 'mape1BBRCE'
    thisLBL2 = r'MAPE(RCE$_{(\mathrm{1BB})}$) $\leq$ ' + f'{maxMAPE2} [%]'
    colors.append('#D55E00')
  elif prop2 == 'RCE2BB':
    thisMAPE2Type = 'mape2BBRCE'
    thisLBL2 = r'MAPE(RCE$_{(\mathrm{2BB})}$) $\leq$ ' + f'{maxMAPE2} [%]'
    colors.append('#1B4F72')
  else:
    print(f'check prop2. Exit()')
    exit()
  
  #print(f'{df1}')  
  thisDF1 = df1[(df1[thisMAPE1Type] <= maxMAPE1)]
  #print(f'{thisDF1}')
  #print(f'{df2}')  
  thisDF2 = df2[(df2[thisMAPE2Type] <= maxMAPE2)]
  #print(f'{thisDF1}')
  thisDFs = [thisDF1, thisDF2]

  value_cols = ['SigC800', 'SigBr730', 'EpsC800', 'EpsBr730']
  y_limits = {
    'SigC800': (0.00, 0.75),
    'SigBr730': (0.00, 0.75),
    'EpsC800': (0.00, 1.00),
    'EpsBr730': (0.00, 3.00),
    }
    
  x = np.arange(len(value_cols))

  plt.figure(figsize=(12, 6))
  for k in range(0, len(thisDFs)):
    df = thisDFs[k]
    if k == 0:
      color1 = colors[k]
    elif k == 1:
      color2 = colors[k]
    else:
      print(f'{k}')
      color = '#000000'
    for _, row in df.iterrows():
      y = [
          (row[col] - y_limits[col][0]) /
          (y_limits[col][1] - y_limits[col][0])
          for col in value_cols
      ]
      plt.plot(x, y, color=colors[k], alpha=0.7)
    
  # Add vertical lines at each x-position
  for xi in x:
    plt.axvline(x=xi, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
  #plt.xticks(x, value_cols)
  xlbls = [r'$\sigma_{\mathrm{C}_{\mathrm{Br}}}$', r'$\sigma_{\mathrm{Br}}$', r'$\varepsilon_{\mathrm{C}_{\mathrm{Br}}}$', r'$\varepsilon_{\mathrm{Br}}$']
  plt.xticks(x, xlbls, fontsize=15, fontweight='bold')
  plt.tick_params(axis='x', pad=15)
  plt.ylim(0, 1)
  plt.yticks([])
  
  # Annotate per-axis scales
  tick_length = 0.03  # how long the ticks should be
  for i, col in enumerate(value_cols):
    # Draw the main vertical axis line
    plt.plot([i, i], [0, 1], color='black', linestyle='--', linewidth=1, alpha=0.5)
    
    # Draw ticks at normalized positions (e.g., bottom, middle, top)
    num_ticks = 6
    for t in np.linspace(0, 1, num_ticks):
        plt.plot([i - tick_length/2, i + tick_length/2], [t, t], color='black', linewidth=1)
        
    lo, hi = y_limits[col]
    plt.text(i, -0.015, f"{lo}", ha='center', va='top', fontsize=15, fontweight='bold')
    plt.text(i, 1.02, f"{hi}", ha='center', va='bottom', fontsize=15, fontweight='bold')
  
  legend_elements = [
    Line2D([0], [0], color=color1, lw=2, label=thisLBL1),
    Line2D([0], [0], color=color2, lw=2, label=thisLBL2)
  ]
  plt.legend(handles=legend_elements, loc='lower right', fontsize=12)
  plt.ylabel(r'$\sigma_{\mathrm{C}_{\mathrm{Br}}}$ and $\sigma_{\mathrm{Br}}$ in [nm]; $\varepsilon_{\mathrm{C}_{\mathrm{Br}}}$ and $\varepsilon_{\mathrm{Br}}$ in [kJ/mol]', fontsize=15)
  #plt.title(f'{substance} - {sorting}')
  
  if saveOrShow == 'save':
    plt.savefig(f'PC_MAPEs_{prop1}-vs-{prop2}.png', dpi=figDPI)
  
  return plt
  

def printTbl(df):
  #print(f'{df}')
  dfNNR = df[(df["method"] == 'NNR')]
  dfPR = df[(df["method"] == 'PR')]
  thisDFs = [dfNNR, dfPR]
  
  for thisDF in thisDFs:
    print(f'{thisDF["method"].unique()}')
    print("$\\sigma_{\\mathrm{C}}$ & $\\sigma_{\\mathrm{Br}}$ & $\\varepsilon_{\\mathrm{C}}$ & $\\varepsilon_{\\mathrm{Br}}$ \\\\")
    for i, row in thisDF.iterrows():
      if i%2 == 1:
        print("    \\rowcolor[HTML]{D3EEFA}")
      print(f'    {row.SigC800:.6f} & {row.SigBr730:.6f} & {row.EpsC800:.6f} & {row.SigBr730:.6f} \\\\')
    print(f'')  
    

def readThisCSV(pathToFile):
  try:
    thisDF = pd.read_csv(pathToFile)
  except:
    print(f'Could not read csv file \"{pathToFile}\". Exit().')
    exit()
  else:
    #print(f'{thisDF}')
    try:
      thisDF = thisDF.drop(columns=["Unnamed: 0"])
    except:
      pass
    else:
      #print(f'{thisDF}')
      pass
      
    return thisDF


def getDensMAPE(dfDens, target):
  densities = dfDens.density.to_numpy()
  #print(f'{densities}') 
  mapes = []
  for density in densities:
    #print(f'{density}') 
    mape = float((np.absolute(target - density)/np.absolute(target))*100)
    #print(f'{mape}')
    mapes.append(mape)
  return np.array(mapes)


def getDensRCE(dfRCE, targets):
  mapes = []
  RCENames = [rce for rce in dfRCE.columns if rce.startswith("RCE")] 
  for i, row in dfRCE.iterrows():
    mape = 0
    for k in range(0,len(RCENames)):
      #print(f'{row[RCENames[k]]}')
      if targets[k] == 0:
        pass
      else:
        thisMape = float((np.absolute(targets[k] - row[RCENames[k]])/np.absolute(targets[k]))*100)
        #print(f'{thisMape}')
        mape = mape + thisMape
    mape = float(mape/len(RCENames))
    mapes.append(mape)
  #print(f'{mapes}')
  return mapes 


if __name__ == "__main__":
  pwd = os.getcwd()
  #print(f'{pwd}')
  DensFile1BB = os.path.join(pwd, DensFile1BB)
  DensFile2BB = os.path.join(pwd, DensFile2BB)
  RCEFile1BB = os.path.join(pwd, RCEFile1BB)
  RCEFile2BB = os.path.join(pwd, RCEFile2BB)
  targetDensFile = os.path.join(pwd, targetDensFile)
  targetRCEFile = os.path.join(pwd, targetRCEFile)
  
  dfDens1BB = readThisCSV(DensFile1BB)
  dfDens2BB = readThisCSV(DensFile2BB)
  dftargetDens = readThisCSV(targetDensFile)
  target1BBDens = dftargetDens[(dftargetDens["substance"] == "1-bromobutane")].density.to_numpy()
  target2BBDens = dftargetDens[(dftargetDens["substance"] == "2-bromobutane")].density.to_numpy()
  mape1BBdens = getDensMAPE(dfDens1BB, target1BBDens)
  mape2BBdens = getDensMAPE(dfDens2BB, target2BBDens)
  
  dfSummaryDens1BB = pd.DataFrame({"SigC800": dfDens1BB["SigC800"],
                                   "SigBr730": dfDens1BB["SigBr730"],
                                   "EpsC800": dfDens1BB["EpsC800"],
                                   "EpsBr730": dfDens1BB["EpsBr730"],
                                   "mape1BBDens": mape1BBdens})
  #print(f'{dfSummaryDens1BB}')
  dfSummaryDens2BB = pd.DataFrame({"SigC800": dfDens2BB["SigC800"],
                                   "SigBr730": dfDens2BB["SigBr730"],
                                   "EpsC800": dfDens2BB["EpsC800"],
                                   "EpsBr730": dfDens2BB["EpsBr730"],
                                   "mape2BBDens": mape2BBdens})
  #print(f'{dfSummaryDens2BB}')
  
  dfRCE1BB = readThisCSV(RCEFile1BB)
  dfRCE2BB = readThisCSV(RCEFile2BB)
  dftargetRCE = readThisCSV(targetRCEFile)
  target1BBRCE = dftargetRCE[(dftargetRCE["Substance"] == "1-bromobutane")].RCE.to_numpy()
  target2BBRCE = dftargetRCE[(dftargetRCE["Substance"] == "2-bromobutane")].RCE.to_numpy()
  mape1BBRCE = getDensRCE(dfRCE1BB, target1BBRCE)
  mape2BBRCE = getDensRCE(dfRCE2BB, target2BBRCE)
  
  dfSummaryRCE1BB = pd.DataFrame({"SigC800": dfRCE1BB["SigC800"],
                                  "SigBr730": dfRCE1BB["SigBr730"],
                                  "EpsC800": dfRCE1BB["EpsC800"],
                                  "EpsBr730": dfRCE1BB["EpsBr730"],
                                  "mape1BBRCE": mape1BBRCE})
  #print(f'{dfSummaryRCE1BB}')
  dfSummaryRCE2BB = pd.DataFrame({"SigC800": dfRCE2BB["SigC800"],
                                  "SigBr730": dfRCE2BB["SigBr730"],
                                  "EpsC800": dfRCE2BB["EpsC800"],
                                  "EpsBr730": dfRCE2BB["EpsBr730"],
                                  "mape2BBRCE": mape2BBRCE})
  #print(f'{dfSummaryRCE2BB}')

  plt = mkPlotPC(dfSummaryDens1BB, dfSummaryDens2BB, 'Dens1BB', 'Dens2BB', 1.5, 1.5)
  plt = mkPlotPC(dfSummaryRCE1BB, dfSummaryRCE2BB, 'RCE1BB', 'RCE2BB', 250, 30)
  plt = mkPlotPC(dfSummaryDens1BB, dfSummaryRCE1BB, 'Dens1BB', 'RCE1BB', 1.5, 250)
  plt = mkPlotPC(dfSummaryDens2BB, dfSummaryRCE2BB, 'Dens2BB', 'RCE2BB', 1.5, 30)
  if saveOrShow == 'show':
    plt.show()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
