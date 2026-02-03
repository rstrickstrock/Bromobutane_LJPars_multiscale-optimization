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
top = 5
figDPI = 100

try:
  saveOrShow = sys.argv[1]
except:
  saveOrShow = "show"
if saveOrShow == "save":
  pass
else:
  saveOrShow = "show"

def concatDFs(df1, df2, method, sorting, n):
  dfThisSummary = pd.DataFrame({"method": [f'{method}'] * n,
                                "sort": [f'{sorting}'] * n,
                                "Dens1SimMAPE": df2.head(top)["Dens1SimMAPE"],
                                "Dens2SimMAPE": df2.head(top)["Dens2SimMAPE"],
                                "RCE1MAPE": df2.head(top)["RCE1MAPE"],
                                "RCE2MAPE": df2.head(top)["RCE2MAPE"],
                                "SigC800": df2.head(top)["SigC800"],
                                "SigBr730": df2.head(top)["SigBr730"],
                                "EpsC800": df2.head(top)["EpsC800"],
                                "EpsBr730": df2.head(top)["EpsBr730"],
                                "#": df2.head(top)["#"]})
  #print(f'{dfThisSummary}')
  return pd.concat([df1, dfThisSummary], ignore_index=True)


def mkPlotPC1(df, prop, saveOrShow, substance):
  legend_elements = []
  if prop == 'density' and substance == "1-bromobutane":
    dfNNR = df[(df["sort"] == 'Dens1SimMAPE') & (df["method"] == 'NNR')]
    dfPR = df[(df["sort"] == 'Dens1SimMAPE') & (df["method"] == 'PR')]
    legend_elements.append(Line2D([0], [0], color='#8E44AD', lw=2, label=r'best MAPE$(\rho_{(1BB)})$'))
    xlbls = [r'MAPE$(\rho_{(1BB)})$',r'MAPE$(\mathrm{RCE}_{(1BB)})$']
    value_cols = ['Dens1SimMAPE', 'RCE1MAPE']
    y_limits = {
    'Dens1SimMAPE': (0, 5),
    'RCE1MAPE': (0, 350),
    }
  elif prop == 'density' and substance == "2-bromobutane":
    dfNNR = df[(df["sort"] == 'Dens2SimMAPE') & (df["method"] == 'NNR')]
    dfPR = df[(df["sort"] == 'Dens2SimMAPE') & (df["method"] == 'PR')]
    legend_elements.append(Line2D([0], [0], color='#2ECC71', lw=2, label=r'best MAPE$(\rho_{(2BB)})$'))
    xlbls = [r'MAPE$(\rho_{(2BB)})$',r'MAPE$(\mathrm{RCE}_{(2BB)})$']
    value_cols = ['Dens2SimMAPE', 'RCE2MAPE']
    y_limits = {
    'Dens2SimMAPE': (0, 5),
    'RCE2MAPE': (0, 100),
    }
  elif prop == 'RCE' and substance == "1-bromobutane":
    dfNNR = df[(df["sort"] == 'RCE1MAPE') & (df["method"] == 'NNR')]
    dfPR = df[(df["sort"] == 'RCE1MAPE') & (df["method"] == 'PR')]
    legend_elements.append(Line2D([0], [0], color='#D55E00', lw=2, label=r'MAPE$(\mathrm{RCE}_{(1BB)})$'))
    xlbls = [r'MAPE$(\rho_{(1BB)})$',r'MAPE$(\mathrm{RCE}_{(1BB)})$']
    value_cols = ['Dens1SimMAPE', 'RCE1MAPE']
    y_limits = {
    'Dens1SimMAPE': (0, 5),
    'RCE1MAPE': (0, 350),
    }
  elif prop == 'RCE' and substance == "2-bromobutane":
    dfNNR = df[(df["sort"] == 'RCE2MAPE') & (df["method"] == 'NNR')]
    dfPR = df[(df["sort"] == 'RCE2MAPE') & (df["method"] == 'PR')]
    legend_elements.append(Line2D([0], [0], color='#1B4F72', lw=2, label=r'MAPE$(\mathrm{RCE}_{(2BB)})$'))
    xlbls = [r'MAPE$(\rho_{(2BB)})$',r'MAPE$(\mathrm{RCE}_{(2BB)})$']
    value_cols = ['Dens2SimMAPE', 'RCE2MAPE']
    y_limits = {
    'Dens2SimMAPE': (0, 5),
    'RCE2MAPE': (0, 100),
    }
  else:
    print(f'set prop or substance correctly. Exit.')
    exit()
    
  thisDFs = [dfNNR, dfPR]

  x = np.arange(len(value_cols))

  plt.figure(figsize=(12, 6))
  for df in thisDFs:
    for _, row in df.iterrows():
      y = [
          (row[col] - y_limits[col][0]) /
          (y_limits[col][1] - y_limits[col][0])
          for col in value_cols
      ]
      if row.method == "NNR" and row.sort == "Dens1SimMAPE":
        color = '#8E44AD'
      elif row.method == "NNR" and row.sort == "Dens2SimMAPE":
        color = '#2ECC71'
      elif row.method == "NNR" and row.sort == "RCE1MAPE":
        color = '#D55E00'
      elif row.method == "NNR" and row.sort == "RCE2MAPE":
        color = '#1B4F72'
      elif row.method == "PR" and row.sort == "Dens1SimMAPE":
        color = '#8E44AD'
      elif row.method == "PR" and row.sort == "Dens2SimMAPE":
        color = '#2ECC71'
      elif row.method == "PR" and row.sort == "RCE1MAPE":
        color = '#D55E00'
      elif row.method == "PR" and row.sort == "RCE2MAPE":
        color = '#1B4F72'
      else:
        print(f'{row}')
        color = '#000000'
      plt.plot(x, y, color=color, alpha=0.7)
    
  # Add vertical lines at each x-position
  for xi in x:
    plt.axvline(x=xi, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
  plt.xticks(x, xlbls, fontsize=15, fontweight='bold')
  plt.tick_params(axis='x', pad=15)
  plt.ylim(0, 1)
  plt.yticks([])
  # Annotate per-axis scales
  for i, col in enumerate(value_cols):
    lo, hi = y_limits[col]
    plt.text(i, -0.015, f"{lo}", ha='center', va='top', fontsize=14)
    plt.text(i, 1.02, f"{hi}", ha='center', va='bottom', fontsize=14)

  plt.legend(handles=legend_elements, loc='upper right', fontsize=12)
  plt.ylabel("MAPE [%]", fontsize=15, fontweight='bold')
  #plt.title(f'{substance} - {sorting}')
  
  if saveOrShow == 'save':
    plt.savefig(f'PC_NNR+PR_sorted-by-{prop}_{substance}.png', dpi=100)
  
  return plt  


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
  
  value_cols = ['Dens1SimMAPE', 'Dens2SimMAPE', 'RCE1MAPE', 'RCE2MAPE']
  y_limits = {
    'Dens1SimMAPE': (0, 5),
    'Dens2SimMAPE': (0, 5),
    'RCE1MAPE': (0, 350),
    'RCE2MAPE': (0, 100),
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
  plt.legend(handles=legend_elements, loc='upper right', fontsize=12)
  plt.ylabel("MAPE [%]", fontsize=15, fontweight='bold')
  #plt.title(f'{substance} - {sorting}')
  
  if saveOrShow == 'save':
    plt.savefig(f'PC_NNR-vs-PR_{substance}.png', dpi=100)
  
  return plt


def mkPlotPC3(df, prop, saveOrShow, substance):
  plotSort = 1
  if prop == 'density':
    df1NNR = df[(df["sort"] == 'Dens1SimMAPE') & (df["method"] == 'NNR')]
    df2NNR = df[(df["sort"] == 'Dens2SimMAPE') & (df["method"] == 'NNR')]
    df1PR = df[(df["sort"] == 'Dens1SimMAPE') & (df["method"] == 'PR')]
    df2PR = df[(df["sort"] == 'Dens2SimMAPE') & (df["method"] == 'PR')]
    legend_elements = [
      Line2D([0], [0], color='#8E44AD', lw=2, label=r'best MAPE$(\rho_{(1BB)})$'),
      Line2D([0], [0], color='#2ECC71', lw=2, label=r'best MAPE$(\rho_{(2BB)})$')
    ]
  elif prop == 'RCE':
    df1NNR = df[(df["sort"] == 'RCE1MAPE') & (df["method"] == 'NNR')]
    df2NNR = df[(df["sort"] == 'RCE2MAPE') & (df["method"] == 'NNR')]
    df1PR = df[(df["sort"] == 'RCE1MAPE') & (df["method"] == 'PR')]
    df2PR = df[(df["sort"] == 'RCE2MAPE') & (df["method"] == 'PR')]
    legend_elements = [
      Line2D([0], [0], color='#D55E00', lw=2, label=r'MAPE$(\mathrm{RCE}_{(1BB)})$'),
      Line2D([0], [0], color='#1B4F72', lw=2, label=r'MAPE$(\mathrm{RCE}_{(2BB)})$')
    ]
  else:
    print(f'set prop correctly. Exit.')
    exit()
    
  thisDFs = [df1NNR, df2NNR, df1PR, df2PR]
  
  if plotSort == 1:
    xlbls = [r'MAPE$(\rho_{(1BB)})$',r'MAPE$(\rho_{(2BB)})$',r'MAPE$(\mathrm{RCE}_{(1BB)})$',r'MAPE$(\mathrm{RCE}_{(2BB)})$']
    value_cols = ['Dens1SimMAPE', 'Dens2SimMAPE', 'RCE1MAPE', 'RCE2MAPE']
    y_limits = {
      'Dens1SimMAPE': (0, 5),
      'Dens2SimMAPE': (0, 5),
      'RCE1MAPE': (0, 350),
      'RCE2MAPE': (0, 100),
      }
  elif plotSort == 2:
    xlbls = [r'MAPE$(\rho_{(1BB)})$',r'MAPE$(\mathrm{RCE}_{(1BB)})$',r'MAPE$(\rho_{(2BB)})$',r'MAPE$(\mathrm{RCE}_{(2BB)})$']
    value_cols = ['Dens1SimMAPE', 'RCE1MAPE', 'Dens2SimMAPE', 'RCE2MAPE']
    y_limits = {
      'Dens1SimMAPE': (0, 5),
      'RCE1MAPE': (0, 350),
      'Dens2SimMAPE': (0, 5),
      'RCE2MAPE': (0, 100),
    }
  else:
    print(f'check plotSort. Exit.')
    exit()
    
  x = np.arange(len(value_cols))

  plt.figure(figsize=(12, 6))
  for df in thisDFs:
    for _, row in df.iterrows():
      y = [
          (row[col] - y_limits[col][0]) /
          (y_limits[col][1] - y_limits[col][0])
          for col in value_cols
      ]
      if row.method == "NNR" and row.sort == "Dens1SimMAPE":
        color = '#8E44AD'
      elif row.method == "NNR" and row.sort == "Dens2SimMAPE":
        color = '#2ECC71'
      elif row.method == "NNR" and row.sort == "RCE1MAPE":
        color = '#D55E00'
      elif row.method == "NNR" and row.sort == "RCE2MAPE":
        color = '#1B4F72'
      elif row.method == "PR" and row.sort == "Dens1SimMAPE":
        color = '#8E44AD'
      elif row.method == "PR" and row.sort == "Dens2SimMAPE":
        color = '#2ECC71'
      elif row.method == "PR" and row.sort == "RCE1MAPE":
        color = '#D55E00'
      elif row.method == "PR" and row.sort == "RCE2MAPE":
        color = '#1B4F72'
      else:
        print(f'{row}')
        color = '#000000'
      plt.plot(x, y, color=color, alpha=0.7)
    
  # Add vertical lines at each x-position
  for xi in x:
    plt.axvline(x=xi, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
    
  #plt.xticks(x, value_cols)
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
    

  plt.legend(handles=legend_elements, loc='upper right', fontsize=12)
  plt.ylabel("MAPE [%]", fontsize=15, fontweight='bold')
  #plt.title(f'{substance} - {sorting}')
  
  if saveOrShow == 'save':
    plt.savefig(f'PC_NNR+PR_sorted-by-{prop}_{substance}.png', dpi=100)
  
  return plt


def mkPlotScatter1(df, saveOrShow, substance):
  #legend_elements = []
  if substance == "1-bromobutane":
    dfDens = df[(df["sort"] == 'Dens1SimMAPE')]
    dfRCE = df[(df["sort"] == 'RCE1MAPE')]
    
    xlbl = r'MAPE$(\rho_{(1BB)})$ [%]'
    ylbl = r'MAPE$(\mathrm{RCE}_{(1BB)})$ [%]'
    value_cols = ['Dens1SimMAPE', 'RCE1MAPE']
    xlims = [0, 5]
    ylims = [0, 350]
    colors = ['#8E44AD', '#D55E00']
    thisLBLs = [r'best MAPE$(\rho_{(1BB)})$', r'best MAPE$(\mathrm{RCE}_{(1BB)})$']

  elif substance == "2-bromobutane":
    dfDens = df[(df["sort"] == 'Dens2SimMAPE')]
    dfRCE = df[(df["sort"] == 'RCE2MAPE')]
    
    xlbl = r'MAPE$(\rho_{(2BB)})$ [%]'
    ylbl = r'MAPE$(\mathrm{RCE}_{(2BB)})$ [%]'
    value_cols = ['Dens2SimMAPE', 'RCE2MAPE']
    xlims = [0, 5]
    ylims = [0, 125]
    colors = ['#2ECC71', '#1B4F72']
    thisLBLs = [r'best MAPE$(\rho_{(2BB)})$', r'best MAPE$(\mathrm{RCE}_{(2BB)})$']
    
  else:
    print(f'set prop or substance correctly. Exit.')
    exit()


  thisDFs = [dfDens, dfRCE]
  plt.figure(figsize=(6, 6))

  for k in range(0, len(thisDFs)):
    df = thisDFs[k]
    thisX = df[value_cols[0]].to_numpy()
    thisY = df[value_cols[1]].to_numpy()

    plt.scatter(thisX, thisY, color=colors[k], label=thisLBLs[k])

  
  plt.xlabel(xlbl, fontsize=20, fontweight='bold')
  plt.ylabel(ylbl, fontsize=20, fontweight='bold')
  plt.xlim(xlims)
  plt.ylim(ylims)
  plt.xticks(fontsize=15)
  plt.yticks(fontsize=15)
  plt.grid(True, linestyle=':', alpha=0.6)
  plt.legend(loc='best', fontsize=16)
  plt.tight_layout()
  
  if saveOrShow == 'save':
    plt.savefig(f'ScatterPlot_MAPEDens-vs-MAPERCE_{substance}', dpi=100)
  
  return plt


def mkPlotScatter2(df, saveOrShow, substance):
  #color = '#377eb8' #nnr // color = '#ff7f00' #pr
  if substance == "1-bromobutane":
    dfNNR = df[(df["method"] == 'NNR')]
    dfPR = df[(df["method"] == 'PR')] 
    xlbl = r'MAPE$(\rho_{(1BB)})$ [%]'
    ylbl = r'MAPE$(\mathrm{RCE}_{(1BB)})$ [%]'
    value_cols = ['Dens1SimMAPE', 'RCE1MAPE']
    xlims = [0, 5]
    ylims = [0, 350]
    colors = ['#377eb8', '#ff7f00']
    thisLBLs = ['NNR surrog. Models', 'PR surrog. Models']

  elif substance == "2-bromobutane":
    dfNNR = df[(df["method"] == 'NNR')]
    dfPR = df[(df["method"] == 'PR')] 
    xlbl = r'MAPE$(\rho_{(2BB)})$ [%]'
    ylbl = r'MAPE$(\mathrm{RCE}_{(2BB)})$ [%]'
    value_cols = ['Dens2SimMAPE', 'RCE2MAPE']
    xlims = [0, 5]
    ylims = [0, 125]
    colors = ['#377eb8', '#ff7f00']
    thisLBLs = ['NNR surrog. Models', 'PR surrog. Models']
    
  else:
    print(f'set prop or substance correctly. Exit.')
    exit()

  thisDFs = [dfNNR, dfPR]
  plt.figure(figsize=(6, 6))

  for k in range(0, len(thisDFs)):
    df = thisDFs[k]
    thisX = df[value_cols[0]].to_numpy()
    thisY = df[value_cols[1]].to_numpy()

    plt.scatter(thisX, thisY, color=colors[k], label=thisLBLs[k])

  
  plt.xlabel(xlbl, fontsize=20, fontweight='bold')
  plt.ylabel(ylbl, fontsize=20, fontweight='bold')
  plt.xlim(xlims)
  plt.ylim(ylims)
  plt.xticks(fontsize=15)
  plt.yticks(fontsize=15)
  plt.grid(True, linestyle=':', alpha=0.6)
  plt.legend(loc='best', fontsize=16)
  plt.tight_layout()
  
  if saveOrShow == 'save':
    plt.savefig(f'ScatterPlot_NNR-vs-PR_{substance}', dpi=100)
  
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
                                "sort": [],
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
                                "sort": [],
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
                                "sort": [],
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
        df1BromoDens = dfthisResults.sort_values("Dens1SimMAPE")
        dfSummaryBoth = concatDFs(dfSummaryBoth, df1BromoDens, method, 'Dens1SimMAPE', top)
        df2BromoDens = dfthisResults.sort_values("Dens2SimMAPE")
        dfSummaryBoth = concatDFs(dfSummaryBoth, df2BromoDens, method, 'Dens2SimMAPE', top)
        df1BromoRCE = dfthisResults.sort_values("RCE1MAPE")
        dfSummaryBoth = concatDFs(dfSummaryBoth, df1BromoRCE, method, 'RCE1MAPE', top)
        df2BromoRCE = dfthisResults.sort_values("RCE2MAPE")
        dfSummaryBoth = concatDFs(dfSummaryBoth, df2BromoRCE, method, 'RCE2MAPE', top)
      elif substance == '1-bromobutane':
        df1BromoDens = dfthisResults.sort_values("Dens1SimMAPE")
        dfSummary1Bro = concatDFs(dfSummary1Bro, df1BromoDens, method, 'Dens1SimMAPE', top)
        df1BromoRCE = dfthisResults.sort_values("RCE1MAPE")
        dfSummary1Bro = concatDFs(dfSummary1Bro, df1BromoRCE, method, 'RCE1MAPE', top)
      elif substance == '2-bromobutane':
        df2BromoDens = dfthisResults.sort_values("Dens2SimMAPE")
        dfSummary2Bro = concatDFs(dfSummary2Bro, df2BromoDens, method, 'Dens2SimMAPE', top)
        df2BromoRCE = dfthisResults.sort_values("RCE2MAPE")
        dfSummary2Bro = concatDFs(dfSummary2Bro, df2BromoRCE, method, 'RCE2MAPE', top)

  #print(f'{dfSummaryBoth}')
  #print(f'{dfSummary1Bro}')
  #print(f'{dfSummary2Bro}')
  
 
  #plt = mkPlotPC3(dfSummaryBoth, 'density', saveOrShow, 'both')
  #plt = mkPlotPC3(dfSummaryBoth, 'RCE', saveOrShow, 'both')
  #plt = mkPlotPC2(dfSummaryBoth, 'both')
  
  plt = mkPlotScatter1(dfSummary1Bro, saveOrShow, '1-bromobutane')
  plt = mkPlotScatter2(dfSummary1Bro, saveOrShow, '1-bromobutane')
  plt = mkPlotScatter1(dfSummary2Bro, saveOrShow, '2-bromobutane')
  plt = mkPlotScatter2(dfSummary2Bro, saveOrShow, '2-bromobutane')
  if saveOrShow == 'show':
    plt.show()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
