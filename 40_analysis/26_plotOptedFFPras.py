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


def mkPlotPC(df):
  dfNNR = df[(df["method"] == 'NNR')]
  #dfNNR = df[(df["sort"] == 'Dens2SimMAPE') & (df["method"] == 'NNR')]
  #dfR1NNR = df[(df["sort"] == 'RCE1MAPE') & (df["method"] == 'NNR')]
  #dfR2NNR = df[(df["sort"] == 'RCE2MAPE') & (df["method"] == 'NNR')]
  
  dfPR = df[(df["method"] == 'PR')]
  #dfD2PR = df[(df["sort"] == 'Dens2SimMAPE') & (df["method"] == 'PR')]
  #dfR1PR = df[(df["sort"] == 'RCE1MAPE') & (df["method"] == 'PR')]
  #dfR2PR = df[(df["sort"] == 'RCE2MAPE') & (df["method"] == 'PR')]
  thisDFs = [dfNNR, dfPR]
  
  maxSigC = df["SigC800"].max()
  maxSigBr = df["SigBr730"].max()
  maxEpsC = df["EpsC800"].max()
  maxEpsBr = df["EpsBr730"].max()
  
  value_cols = ['SigC800', 'SigBr730', 'EpsC800', 'EpsBr730']
  y_limits = {
    'SigC800': (0.00, 0.75),
    'SigBr730': (0.00, 0.75),
    'EpsC800': (0.00, 1.00),
    'EpsBr730': (0.00, 3.00),
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
    Line2D([0], [0], color='#377eb8', lw=2, label='NNR surrog. Models'),
    Line2D([0], [0], color='#ff7f00', lw=2, label='PR surrog. Models')
  ]
  plt.legend(handles=legend_elements, loc='lower right', fontsize=12)
  plt.ylabel(r'$\sigma_{\mathrm{C}_{\mathrm{Br}}}$ and $\sigma_{\mathrm{Br}}$ in [nm]; $\varepsilon_{\mathrm{C}_{\mathrm{Br}}}$ and $\varepsilon_{\mathrm{Br}}$ in [kJ/mol]', fontsize=15)
  #plt.title(f'{substance} - {sorting}')
  
  if saveOrShow == 'save':
    plt.savefig(f'PC_optimizedFFPars_simultaneous.png', dpi=100)
  
  return plt
  

def mkPlotPC2(df1, df2):
  #print(f'{df1["Dens1SimMAPE"].max()}')
  #print(f'{df1["RCE1MAPE"].max()}')
  #print(f'{df2["Dens2SimMAPE"].max()}')
  #print(f'{df2["RCE2MAPE"].max()}')
  #exit()
  thisDFs = [df1, df2]
  
  #maxSigC = df["SigC800"].max()
  #maxSigBr = df["SigBr730"].max()
  #maxEpsC = df["EpsC800"].max()
  #maxEpsBr = df["EpsBr730"].max()
  
  value_cols = ['SigC800', 'SigBr730', 'EpsC800', 'EpsBr730']
  y_limits = {
    'SigC800': (0.00, 0.75),
    'SigBr730': (0.00, 0.75),
    'EpsC800': (0.00, 1.00),
    'EpsBr730': (0.00, 3.00),
    }
    
  x = np.arange(len(value_cols))

  plt.figure(figsize=(12, 6))
  for k in range(0,len(thisDFs)):
    df = thisDFs[k]
    for _, row in df.iterrows():
      y = [
          (row[col] - y_limits[col][0]) /
          (y_limits[col][1] - y_limits[col][0])
          for col in value_cols
      ]
      if k == 0:
        color = '#377eb8'
      elif k == 1:
        color = '#ff7f00'
      else:
        print(f'{row}')
        color = '#000000'
      plt.plot(x, y, color=color, alpha=0.7)
    
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
    Line2D([0], [0], color='#377eb8', lw=2, label='opt. 1BB separately'),
    Line2D([0], [0], color='#ff7f00', lw=2, label='opt. 2BB separately')
  ]
  plt.legend(handles=legend_elements, loc='lower right', fontsize=12)
  plt.ylabel(r'$\sigma_{\mathrm{C}_{\mathrm{Br}}}$ and $\sigma_{\mathrm{Br}}$ in [nm]; $\varepsilon_{\mathrm{C}_{\mathrm{Br}}}$ and $\varepsilon_{\mathrm{Br}}$ in [kJ/mol]', fontsize=15)
  #plt.title(f'{substance} - {sorting}')
  
  if saveOrShow == 'save':
    plt.savefig(f'PC_optimizedFFPars_separately.png', dpi=100)
  
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
 
  plt = mkPlotPC(dfSummaryBoth)
  #printTbl(dfSummaryBoth)
  
  plt = mkPlotPC2(dfSummary1Bro, dfSummary2Bro)
  #printTbl(dfSummary1Bro)
  #printTbl(dfSummary2Bro)
  if saveOrShow == 'show':
    plt.show()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
