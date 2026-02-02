import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from pandas.plotting import parallel_coordinates
import plotly.graph_objects as go
from matplotlib.lines import Line2D
import sys

dirs = ["30_1+2-Bromobutane_opt", "31_1-Bromobutane_opt", "31_2-Bromobutane_opt", "32_1+2-Bromobutane_opt_PR", "33_1-Bromobutane_opt_PR", "33_2-Bromobutane_opt_PR"]

optResultsFileName = "optResults_withDensitySims.csv"
figDPI = 100

try:
  saveOrShow = sys.argv[1]
except:
  saveOrShow = "show"
if saveOrShow == "save":
  pass
else:
  saveOrShow = "show"


def concatDFs(df1, df2, method):
  n = len(df2)
  dfThisSummary = pd.DataFrame({"method": [f'{method}'] * n,
                                "Dens1SimMAPE": df2["Dens1SimMAPE"],
                                "Dens2SimMAPE": df2["Dens2SimMAPE"],
                                "RCE1MAPE": df2["RCE1MAPE"],
                                "RCE2MAPE": df2["RCE2MAPE"],
                                "SigC800": df2["SigC800"],
                                "SigBr730": df2["SigBr730"],
                                "EpsC800": df2["EpsC800"],
                                "EpsBr730": df2["EpsBr730"],
                                "#": df2["#"]})
  #print(f'{dfThisSummary}')
  return pd.concat([df1, dfThisSummary], ignore_index=True)


def mkPlotPC2(df, substance):
  dfNNR = df[(df["method"] == 'NNR')]
  #dfNNR = df[(df["sort"] == 'Dens2SimMAPE') & (df["method"] == 'NNR')]
  #dfR1NNR = df[(df["sort"] == 'RCE1MAPE') & (df["method"] == 'NNR')]
  #dfR2NNR = df[(df["sort"] == 'RCE2MAPE') & (df["method"] == 'NNR')]
  
  dfPR = df[(df["method"] == 'PR')]
  #dfD2PR = df[(df["sort"] == 'Dens2SimMAPE') & (df["method"] == 'PR')]
  #dfR1PR = df[(df["sort"] == 'RCE1MAPE') & (df["method"] == 'PR')]
  #dfR2PR = df[(df["sort"] == 'RCE2MAPE') & (df["method"] == 'PR')]
  thisDFs = [dfNNR, dfPR]
  
  maxDens1SimMAPE = int(np.ceil(df["Dens1SimMAPE"].max()))
  maxDens2SimMAPE = int(np.ceil(df["Dens2SimMAPE"].max()))
  maxRCE1MAPE = df["RCE1MAPE"].max()
  maxRCE1MAPE = int((np.ceil(maxRCE1MAPE/100))*100)
  maxRCE2MAPE = df["RCE2MAPE"].max()
  maxRCE2MAPE = int((np.ceil(maxRCE2MAPE/10))*10)
  value_cols = ['Dens1SimMAPE', 'Dens2SimMAPE', 'RCE1MAPE', 'RCE2MAPE']
  y_limits = {
    'Dens1SimMAPE': (0, maxDens1SimMAPE),
    'Dens2SimMAPE': (0, maxDens2SimMAPE),
    'RCE1MAPE': (0, maxRCE1MAPE),
    'RCE2MAPE': (0, maxRCE2MAPE),
    }
    
  x = np.arange(len(value_cols))

  plt.figure(figsize=(12, 6))
  for df in thisDFs:
    for _, row in df.iterrows():
      y = [
          (row[col] - y_limits[col][0]) /
          (y_limits[col][1] - y_limits[col][0])
          for col in value_cols
      ]
      if row.method == "NNR":
        color = '#377eb8'
      elif row.method == "PR":
        color = '#ff7f00'
      else:
        print(f'{row}')
        color = '#000000'
      plt.plot(x, y, color=color, alpha=0.7)
    
  # Add vertical lines at each x-position
  for xi in x:
    plt.axvline(x=xi, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
  #plt.xticks(x, value_cols)
  xlbls = [r'MAPE$(\rho_{(1BB)})$',r'MAPE$(\rho_{(2BB)})$',r'MAPE$(\mathrm{RCE}_{(1BB)})$',r'MAPE$(\mathrm{RCE}_{(2BB)})$']
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
    Line2D([0], [0], color='#377eb8', lw=2, label='NNR surrog. Models'),
    Line2D([0], [0], color='#ff7f00', lw=2, label='PR surrog. Models')
  ]
  plt.legend(handles=legend_elements, loc='lower right', fontsize=12)
  plt.ylabel("MAPE [%]", fontsize=15, fontweight='bold')
  #plt.title(f'{substance} - {sorting}')
  
  if saveOrShow == 'save':
    plt.savefig(f'PC_NNR-vs-PR_{substance}.png', dpi=100)
  
  return plt
  

if __name__ == "__main__":
  pwd = os.path.dirname(os.getcwd())
  #print(f'{pwd}')
  try:
    densTargetFile = os.path.join(os.getcwd(), densTargetFile)
  except:
    densTargetFile = None
  #print(f'{densTargetFile}')
  try:
    rceTargetFile = os.path.join(os.getcwd(), rceTargetFile)
  except:
    rceTargetFile = None
  #print(f'{rceTargetFile}')
  
  dfSummaryBoth = pd.DataFrame({"method": [],
                                "Dens1SimMAPE": [],
                                "Dens2SimMAPE": [],
                                "RCE1MAPE": [],
                                "RCE2MAPE": [],
                                "SigC800": [],
                                "SigBr730": [],
                                "EpsC800": [],
                                "EpsBr730": [],
                                "#": []})
  dfSummary1Bro = pd.DataFrame({"method": [],
                                "Dens1SimMAPE": [],
                                "Dens2SimMAPE": [],
                                "RCE1MAPE": [],
                                "RCE2MAPE": [],
                                "SigC800": [],
                                "SigBr730": [],
                                "EpsC800": [],
                                "EpsBr730": [],
                                "#": []})
  dfSummary2Bro = pd.DataFrame({"method": [],
                                "Dens1SimMAPE": [],
                                "Dens2SimMAPE": [],
                                "RCE1MAPE": [],
                                "RCE2MAPE": [],
                                "SigC800": [],
                                "SigBr730": [],
                                "EpsC800": [],
                                "EpsBr730": [],
                                "#": []})
  for thisDir in dirs:
    if '1+2' in thisDir:
      substance = 'both'
    elif '1-Bromobutane' in thisDir:
      substance = '1-bromobutane'
    elif '2-Bromobutane' in thisDir:
      substance = '2-bromobutane'
    else:
      print(f'Substance could not be set. Exit.')
      exit()
    #print(f'main(): {substance}')
    
    thisDir = os.path.join(pwd, thisDir, '01_C4Br_C3Br')
    if substance == 'both':
      thisDirs = [os.path.join(thisDir, '1-bromo_initPars')]#, os.path.join(thisDir, '2-bromo_initPars')]
    else:
      thisDirs = [thisDir]
      
    for thisDir in thisDirs:
      #print(f'substance: {substance}')
      #print(f'thisDir: {thisDir}')
      thisResultsFilePath = os.path.join(thisDir, optResultsFileName)
      try:
        dfthisResults = pd.read_csv(thisResultsFilePath)
      except:
        print(f'Could not read {thisResultsFilePath} as csv file.')
      else:
        #print(f'Check.')
        try:
          dfthisResults = dfthisResults.drop(columns=["Unnamed: 0"])
        except:
          pass
        else:
          #print(f'{dfthisResults}')
          pass
      
      if "_PR" in thisDir:
        method = "PR"
      else:
        method = "NNR"
      #print(f'{method}')
      
      if substance == 'both':
        dfSummaryBoth = concatDFs(dfSummaryBoth, dfthisResults, method)
      elif substance == '1-bromobutane':
        dfSummary1Bro = concatDFs(dfSummary1Bro, dfthisResults, method)
      elif substance == '2-bromobutane':
        dfSummary2Bro = concatDFs(dfSummary2Bro, dfthisResults, method)

  #print(f'{dfSummaryBoth}')
  #print(f'{dfSummary1Bro}')
  #print(f'{dfSummary2Bro}')
 
  plt = mkPlotPC2(dfSummaryBoth, 'both')
  if saveOrShow == 'show':
    plt.show()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
