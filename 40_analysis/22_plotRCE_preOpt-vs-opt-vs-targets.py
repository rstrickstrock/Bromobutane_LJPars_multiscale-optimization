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
targetFileName = "20_bromobutane_RCE.target"
preOptFileName = "20_bromobutane_RCE.preOpt.csv"

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


def readRCEFile(filePath):
  try:
    dfThisRCEs = pd.read_csv(filePath)
  except:
    print(f'can not pd.read_csv({filePath}). Exit.')
    exit()
  else:
    try:
      dfThisRCEs = dfThisRCEs.drop(columns=["Unnamed: 0"])
    except:
      pass
      
  return dfThisRCEs

def concatDFs(df1, df2, method, sorting, n):
  dfThisSummary = pd.DataFrame({"method": [f'{method}'] * n,
                                "sort": [f'{sorting}'] * n,
                                "RCE1MAPE": df2.head(n)["RCE1MAPE"],
                                "RCE2MAPE": df2.head(n)["RCE2MAPE"],
                                "RCE1": df2.head(n)["RCE1"],
                                "RCE2": df2.head(n)["RCE2"],
                                "RCE3": df2.head(n)["RCE3"],
                                "RCE4": df2.head(n)["RCE4"],
                                "RCE5": df2.head(n)["RCE5"],
                                "RCE6": df2.head(n)["RCE6"],
                                "RCE7": df2.head(n)["RCE7"],
                                "RCE8": df2.head(n)["RCE8"],
                                "RCE9": df2.head(n)["RCE9"],
                                "RCE21": df2.head(n)["RCE21"],
                                "RCE22": df2.head(n)["RCE22"],
                                "RCE23": df2.head(n)["RCE23"],
                                "#": df2.head(n)["#"]})
  #print(f'{dfThisSummary}')
  return pd.concat([df1, dfThisSummary], ignore_index=True)
  

def mkScatterPlot(df, RCEtargets, preOptRCE, sort, substance):
  targets1BB = RCEtargets["RCE"].to_numpy()[:8]
  #print(f'{targets1BB}')
  idx1BB = np.argsort(targets1BB)
  targets2BB = RCEtargets["RCE"].to_numpy()[9:]
  #print(f'{targets2BB}')
  idx2BB = np.argsort(targets2BB)
  
  preOpt1BB = preOptRCE["RCE"].to_numpy()[:8]
  #print(f'{preOpt1BB}')
  preOpt2BB = preOptRCE["RCE"].to_numpy()[9:]
  #print(f'{preOpt2BB}')
  
  srtTargets1BB = []
  srtTargets2BB = []
  srtpreOpt1BB = []
  srtpreOpt2BB = []
  
  for idx in idx1BB:
    srtTargets1BB.append(targets1BB[idx])
    srtpreOpt1BB .append(preOpt1BB[idx])
  for idx in idx2BB:
    srtTargets2BB.append(targets2BB[idx])
    srtpreOpt2BB .append(preOpt2BB[idx])
    
  if substance == 'both':
    dfNNR = df[(df["method"] == 'NNR')]
    #print(f'{dfNNR}')
    #exit()
    RCE1MAPENNR = dfNNR[(dfNNR["sort"] == sort)]
    #print(f'{RCE1MAPENNR}')
    RCE2MAPENNR = dfNNR[(dfNNR["sort"] == sort)]
    
    dfPR = df[(df["method"] == 'PR')]
    RCE1MAPEPR = dfPR[(dfPR["sort"] == sort)]
    #print(f'{RCE1MAPENNR}')
    RCE2MAPEPR = dfPR[(dfPR["sort"] == sort)]

  plt.figure(figsize=(5, 5))
  plt.plot(srtTargets1BB, srtTargets1BB, '-', color='#000000', label="Targets")
  thisLBL = False
  for i, row in RCE1MAPENNR.iterrows():
    #print(f'{row}')
    thisRCEs = np.array([row.RCE1, row.RCE2, row.RCE3, row.RCE4, row.RCE5, row.RCE6, row.RCE7, row.RCE8, row.RCE9])
    srtThisRCEs = []
    for idx in idx1BB:
      srtThisRCEs.append(thisRCEs[idx])
    if not thisLBL:
      plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#377eb8', label="NNR")
      thisLBL = True
    else:
      plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#377eb8')
  thisLBL = False
  for i, row in RCE1MAPEPR.iterrows():
    #print(f'{row}')
    thisRCEs = np.array([row.RCE1, row.RCE2, row.RCE3, row.RCE4, row.RCE5, row.RCE6, row.RCE7, row.RCE8, row.RCE9])
    srtThisRCEs = []
    for idx in idx1BB:
      srtThisRCEs.append(thisRCEs[idx])
    if not thisLBL:
      plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#ff7f00', label="PR")
      thisLBL = True
    else:
      plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#ff7f00')
  plt.plot(srtTargets1BB, srtpreOpt1BB, '-x', color='#ff0000', label="pre optimization")
  plt.xlabel('rel. conf. Energy [kJ/mol]', fontsize=15, fontweight='bold')
  plt.ylabel('rel. conf. Energy [kJ/mol]', fontsize=15, fontweight='bold')
  #plt.xlim([min(targets1BB.min(), preOpt1BB.min()), max(targets1BB.max(), preOpt1BB.max())])
  #plt.ylim([min(targets1BB.min(), preOpt1BB.min()), max(targets1BB.max(), preOpt1BB.max())])
  plt.legend()
  plt.tight_layout()
  if saveOrShow == 'save':
    plt.savefig(f'RCE_targets-vs-preOpt-vs-Opt-{substance}_{sort}_1BB.png', dpi=100)
  
  plt.figure(figsize=(5, 5))
  plt.plot(srtTargets2BB, srtTargets2BB, '-', color='#000000', label="Targets")
  thisLBL = False
  for i, row in RCE1MAPENNR.iterrows():
    #print(f'{row}')
    thisRCEs = np.array([row.RCE21, row.RCE22, row.RCE23])
    srtThisRCEs = []
    for idx in idx2BB:
      srtThisRCEs.append(thisRCEs[idx])
    if not thisLBL:
      plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#377eb8', label="NNR")
      thisLBL = True
    else:
      plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#377eb8')
  thisLBL = False
  for i, row in RCE2MAPEPR.iterrows():
    #print(f'{row}')
    thisRCEs = np.array([row.RCE21, row.RCE22, row.RCE23])
    srtThisRCEs = []
    for idx in idx2BB:
      srtThisRCEs.append(thisRCEs[idx])
    if not thisLBL:
      plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#ff7f00', label="PR")
      thisLBL = True
    else:
      plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#ff7f00')
  plt.plot(srtTargets2BB, srtpreOpt2BB, '-x', color='#ff0000', label="pre optimization")
  plt.xlabel('rel. conf. Energy [kJ/mol]', fontsize=15, fontweight='bold')
  plt.ylabel('rel. conf. Energy [kJ/mol]', fontsize=15, fontweight='bold')
  plt.xticks(fontsize=13)
  plt.yticks(fontsize=13)
  plt.legend()
  plt.tight_layout()
 
  if saveOrShow == 'save':
    plt.savefig(f'RCE_targets-vs-preOpt-vs-Opt-{substance}_{sort}_2BB.png', dpi=100)
  
  return plt


def mkScatterPlot1(df, RCEtargets, preOptRCE, substance):
  targets1BB = RCEtargets["RCE"].to_numpy()[:8]
  #print(f'{targets1BB}')
  idx1BB = np.argsort(targets1BB)
  targets2BB = RCEtargets["RCE"].to_numpy()[9:]
  #print(f'{targets2BB}')
  idx2BB = np.argsort(targets2BB)
  
  preOpt1BB = preOptRCE["RCE"].to_numpy()[:8]
  #print(f'{preOpt1BB}')
  preOpt2BB = preOptRCE["RCE"].to_numpy()[9:]
  #print(f'{preOpt2BB}')
  
  srtTargets1BB = []
  srtTargets2BB = []
  srtpreOpt1BB = []
  srtpreOpt2BB = []
  
  for idx in idx1BB:
    srtTargets1BB.append(targets1BB[idx])
    srtpreOpt1BB .append(preOpt1BB[idx])
  for idx in idx2BB:
    srtTargets2BB.append(targets2BB[idx])
    srtpreOpt2BB .append(preOpt2BB[idx])
    
  dfNNR = df[(df["method"] == 'NNR')]
  dfPR = df[(df["method"] == 'PR')]
  
  RCE1MAPENNR = dfNNR[(dfNNR["sort"] == 'Dens1SimMAPE')]
  RCE1MAPEPR = dfPR[(dfPR["sort"] == 'Dens1SimMAPE')]
  RCE2MAPEPR = dfPR[(dfPR["sort"] == 'Dens2SimMAPE')]
  RCE2MAPENNR = dfNNR[(dfNNR["sort"] == 'Dens2SimMAPE')]
  
  minRCE1BB = np.min([np.min(targets1BB), np.min(preOpt1BB), RCE1MAPENNR[["RCE1", "RCE2", "RCE3", "RCE4", "RCE5", "RCE6", "RCE7", "RCE8", "RCE9"]].min(numeric_only=True).min(), RCE1MAPEPR[["RCE1", "RCE2", "RCE3", "RCE4", "RCE5", "RCE6", "RCE7", "RCE8", "RCE9"]].min(numeric_only=True).min()])
  minRCE1BB = minRCE1BB - abs(minRCE1BB*0.05)
  maxRCE1BB = np.max([np.max(targets1BB), np.max(preOpt1BB), RCE1MAPENNR[["RCE1", "RCE2", "RCE3", "RCE4", "RCE5", "RCE6", "RCE7", "RCE8", "RCE9"]].max(numeric_only=True).max(), RCE1MAPEPR[["RCE1", "RCE2", "RCE3", "RCE4", "RCE5", "RCE6", "RCE7", "RCE8", "RCE9"]].max(numeric_only=True).max()])
  maxRCE1BB = maxRCE1BB + abs(maxRCE1BB*0.05) 
  
  minRCE2BB = np.min([np.min(targets2BB), np.min(preOpt2BB), RCE1MAPENNR[["RCE21", "RCE22", "RCE23"]].min(numeric_only=True).min(), RCE1MAPEPR[["RCE21", "RCE22", "RCE23"]].min(numeric_only=True).min()])
  minRCE2BB = minRCE2BB - abs(minRCE2BB*0.05)
  maxRCE2BB = np.max([np.max(targets2BB), np.max(preOpt2BB), RCE1MAPENNR[["RCE21", "RCE22", "RCE23"]].max(numeric_only=True).max(), RCE1MAPEPR[["RCE21", "RCE22", "RCE23"]].max(numeric_only=True).max()])
  maxRCE2BB = maxRCE2BB + abs(maxRCE2BB*0.05)   
  
  if substance == 'both':
    plt.figure(figsize=(5, 5))
    plt.plot([minRCE1BB, maxRCE1BB], [minRCE1BB, maxRCE1BB], '-', color='#000000')
    plt.plot(srtTargets1BB, srtTargets1BB, '-X', color='#000000', label="Targets")
    thisLBL = False
    for i, row in RCE1MAPENNR.iterrows():
      #print(f'{row}')
      thisRCEs = np.array([row.RCE1, row.RCE2, row.RCE3, row.RCE4, row.RCE5, row.RCE6, row.RCE7, row.RCE8, row.RCE9])
      srtThisRCEs = []
      for idx in idx1BB:
        srtThisRCEs.append(thisRCEs[idx])
      if not thisLBL:
        plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#8E44AD', label=r'best MAPE$(\rho_{(1BB)})$')
        thisLBL = True
      else:
        plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#8E44AD')
    thisLBL = True
    for i, row in RCE1MAPEPR.iterrows():
      #print(f'{row}')
      thisRCEs = np.array([row.RCE1, row.RCE2, row.RCE3, row.RCE4, row.RCE5, row.RCE6, row.RCE7, row.RCE8, row.RCE9])
      srtThisRCEs = []
      for idx in idx1BB:
        srtThisRCEs.append(thisRCEs[idx])
      if not thisLBL:
        plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#8E44AD', label=r'best MAPE$(\rho_{(1BB)})$')
        thisLBL = True
      else:
        plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#8E44AD')
    thisLBL = False
    for i, row in RCE2MAPENNR.iterrows():
      #print(f'{row}')
      thisRCEs = np.array([row.RCE1, row.RCE2, row.RCE3, row.RCE4, row.RCE5, row.RCE6, row.RCE7, row.RCE8, row.RCE9])
      srtThisRCEs = []
      for idx in idx1BB:
        srtThisRCEs.append(thisRCEs[idx])
      if not thisLBL:
        plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#2ECC71', label=r'best MAPE$(\rho_{(2BB)})$')
        thisLBL = True
      else:
        plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#2ECC71')
    thisLBL = True
    for i, row in RCE2MAPEPR.iterrows():
      #print(f'{row}')
      thisRCEs = np.array([row.RCE1, row.RCE2, row.RCE3, row.RCE4, row.RCE5, row.RCE6, row.RCE7, row.RCE8, row.RCE9])
      srtThisRCEs = []
      for idx in idx1BB:
        srtThisRCEs.append(thisRCEs[idx])
      if not thisLBL:
        plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#2ECC71', label=r'best MAPE$(\rho_{(2BB)})$')
        thisLBL = True
      else:
        plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#2ECC71')
  
    plt.plot(srtTargets1BB, srtpreOpt1BB, '-o', color='#ff0000', label="pre optimization")
  
    plt.xlabel('target RCE [kJ/mol]', fontsize=15, fontweight='bold')
    plt.ylabel('reproduced RCE [kJ/mol]', fontsize=15, fontweight='bold')
    #plt.xlim([min(targets1BB.min(), preOpt1BB.min()), max(targets1BB.max(), preOpt1BB.max())])
    #plt.ylim([min(targets1BB.min(), preOpt1BB.min()), max(targets1BB.max(), preOpt1BB.max())])
    plt.xticks(ticks=[-3.0, -2.0, -1.0, 0.0, 1.0, 2.0], labels=["-3.0", "-2.0", "-1.0", "0.0", "1.0", "2.0"],fontsize=13)
    plt.yticks(ticks=[-3.0, -2.0, -1.0, 0.0, 1.0, 2.0], labels=["-3.0", "-2.0", "-1.0", "0.0", "1.0", "2.0"],fontsize=13)
    plt.xlim([minRCE1BB, maxRCE1BB])
    plt.ylim([minRCE1BB, maxRCE1BB])
    plt.legend()
    plt.tight_layout()
    if saveOrShow == 'save':
      plt.savefig(f'RCE_targets-vs-preOpt-vs-Opt-{substance}_sort-by-density_1BB.png', dpi=100)
    
    plt.figure(figsize=(5, 5))
    plt.plot([minRCE2BB, maxRCE2BB], [minRCE2BB, maxRCE2BB], '-', color='#000000')
    plt.plot(srtTargets2BB, srtTargets2BB, '-X', color='#000000', label="Targets")
    thisLBL = False
    for i, row in RCE1MAPENNR.iterrows():
      #print(f'{row}')
      thisRCEs = np.array([row.RCE21, row.RCE22, row.RCE23])
      srtThisRCEs = []
      for idx in idx2BB:
        srtThisRCEs.append(thisRCEs[idx])
      if not thisLBL:
        plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#8E44AD', label=r'best MAPE$(\rho_{(1BB)})$')
        thisLBL = True
      else:
        plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#8E44AD')
    thisLBL = True
    for i, row in RCE1MAPEPR.iterrows():
      #print(f'{row}')
      thisRCEs = np.array([row.RCE21, row.RCE22, row.RCE23])
      srtThisRCEs = []
      for idx in idx2BB:
        srtThisRCEs.append(thisRCEs[idx])
      if not thisLBL:
        plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#8E44AD', label=r'best MAPE$(\rho_{(1BB)})$')
        thisLBL = True
      else:
        plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#8E44AD')
    thisLBL = False
    for i, row in RCE2MAPENNR.iterrows():
      #print(f'{row}')
      thisRCEs = np.array([row.RCE21, row.RCE22, row.RCE23])
      srtThisRCEs = []
      for idx in idx2BB:
        srtThisRCEs.append(thisRCEs[idx])
      if not thisLBL:
        plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#2ECC71', label=r'best MAPE$(\rho_{(2BB)})$')
        thisLBL = True
      else:
        plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#2ECC71')
    thisLBL = True
    for i, row in RCE2MAPEPR.iterrows():
      #print(f'{row}')
      thisRCEs = np.array([row.RCE21, row.RCE22, row.RCE23])
      srtThisRCEs = []
      for idx in idx2BB:
        srtThisRCEs.append(thisRCEs[idx])
      if not thisLBL:
        plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#2ECC71', label=r'best MAPE$(\rho_{(2BB)})$')
        thisLBL = True
      else:
        plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#2ECC71')
    plt.plot(srtTargets2BB, srtpreOpt2BB, '-o', color='#ff0000', label="pre optimization")
    plt.xlabel('target RCE [kJ/mol]', fontsize=15, fontweight='bold')
    plt.ylabel('reproduced RCE [kJ/mol]', fontsize=15, fontweight='bold')
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.xlim([minRCE2BB, maxRCE2BB])
    plt.ylim([minRCE2BB, maxRCE2BB])
    plt.legend()
    plt.tight_layout()
 
    if saveOrShow == 'save':
      plt.savefig(f'RCE_targets-vs-preOpt-vs-Opt-{substance}_sort-by-density_2BB.png', dpi=100)
  
  return plt


def mkScatterPlot2(df, RCEtargets, preOptRCE, substance):
  targets1BB = RCEtargets["RCE"].to_numpy()[:8]
  #print(f'{targets1BB}')
  idx1BB = np.argsort(targets1BB)
  targets2BB = RCEtargets["RCE"].to_numpy()[9:]
  #print(f'{targets2BB}')
  idx2BB = np.argsort(targets2BB)
  
  preOpt1BB = preOptRCE["RCE"].to_numpy()[:8]
  #print(f'{preOpt1BB}')
  preOpt2BB = preOptRCE["RCE"].to_numpy()[9:]
  #print(f'{preOpt2BB}')
  
  srtTargets1BB = []
  srtTargets2BB = []
  srtpreOpt1BB = []
  srtpreOpt2BB = []
  
  for idx in idx1BB:
    srtTargets1BB.append(targets1BB[idx])
    srtpreOpt1BB .append(preOpt1BB[idx])
  for idx in idx2BB:
    srtTargets2BB.append(targets2BB[idx])
    srtpreOpt2BB .append(preOpt2BB[idx])
    
  #dfNNR = df[(df["method"] == 'NNR')]
  #dfPR = df[(df["method"] == 'PR')]  
  if substance == '1-bromobutane':
    DensMAPENNR = df[(df["sort"] == 'Dens1SimMAPE')]
    #DensMAPEPR = dfPR[(dfPR["sort"] == 'Dens1SimMAPE')]
    RCEMAPENNR = df[(df["sort"] == 'RCE1MAPE')]
    #RCEMAPEPR = dfPR[(dfPR["sort"] == 'RCE1MAPE')]
    minRCE1BB = np.min([np.min(targets1BB), np.min(preOpt1BB), DensMAPENNR[["RCE1", "RCE2", "RCE3", "RCE4", "RCE5", "RCE6", "RCE7", "RCE8", "RCE9"]].min(numeric_only=True).min(), RCEMAPENNR[["RCE1", "RCE2", "RCE3", "RCE4", "RCE5", "RCE6", "RCE7", "RCE8", "RCE9"]].min(numeric_only=True).min()])
    minRCE1BB = minRCE1BB - abs(minRCE1BB*0.05)
    maxRCE1BB = np.max([np.max(targets1BB), np.max(preOpt1BB), DensMAPENNR[["RCE1", "RCE2", "RCE3", "RCE4", "RCE5", "RCE6", "RCE7", "RCE8", "RCE9"]].max(numeric_only=True).max(), RCEMAPENNR[["RCE1", "RCE2", "RCE3", "RCE4", "RCE5", "RCE6", "RCE7", "RCE8", "RCE9"]].max(numeric_only=True).max()])
    maxRCE1BB = maxRCE1BB + abs(maxRCE1BB*0.05)
    
    plt.figure(figsize=(5, 5))
    plt.plot([minRCE1BB, maxRCE1BB], [minRCE1BB, maxRCE1BB], '-', color='#000000')
    plt.plot(srtTargets1BB, srtTargets1BB, '-X', color='#000000', label="Targets")
    thisLBL = False
    for i, row in DensMAPENNR.iterrows():
      #print(f'{row}')
      thisRCEs = np.array([row.RCE1, row.RCE2, row.RCE3, row.RCE4, row.RCE5, row.RCE6, row.RCE7, row.RCE8, row.RCE9])
      srtThisRCEs = []
      for idx in idx1BB:
        srtThisRCEs.append(thisRCEs[idx])
      if not thisLBL:
        plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#8E44AD', label=r'best MAPE$(\rho_{(1BB)})$')
        thisLBL = True
      else:
        plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#8E44AD')
    thisLBL = False
    for i, row in RCEMAPENNR.iterrows():
      #print(f'{row}')
      thisRCEs = np.array([row.RCE1, row.RCE2, row.RCE3, row.RCE4, row.RCE5, row.RCE6, row.RCE7, row.RCE8, row.RCE9])
      srtThisRCEs = []
      for idx in idx1BB:
        srtThisRCEs.append(thisRCEs[idx])
      if not thisLBL:
        plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#D55E00', label=r'best MAPE$(\mathrm{RCE}_{(1BB)})$')
        thisLBL = True
      else:
        plt.plot(srtTargets1BB, srtThisRCEs, 'x', ls='dotted', color='#D55E00')

    plt.plot(srtTargets1BB, srtpreOpt1BB, '-o', color='#ff0000', label="pre optimization")
    plt.xlabel('target RCE [kJ/mol]', fontsize=15, fontweight='bold')
    plt.ylabel('reproduced RCE [kJ/mol]', fontsize=15, fontweight='bold')
    #plt.xlim([min(targets1BB.min(), preOpt1BB.min()), max(targets1BB.max(), preOpt1BB.max())])
    #plt.ylim([min(targets1BB.min(), preOpt1BB.min()), max(targets1BB.max(), preOpt1BB.max())])
    plt.xticks(ticks=[-2.0, -1.0, 0.0, 1.0, 2.0], labels=["-2.0", "-1.0", "0.0", "1.0", "2.0"],fontsize=13)
    plt.yticks(ticks=[-2.0, -1.0, 0.0, 1.0, 2.0], labels=["-2.0", "-1.0", "0.0", "1.0", "2.0"],fontsize=13)
    plt.xlim([minRCE1BB, maxRCE1BB])
    plt.ylim([minRCE1BB, maxRCE1BB])
    plt.legend()
    plt.tight_layout()
    if saveOrShow == 'save':
      plt.savefig(f'RCE_targets-vs-preOpt-vs-Opt-{substance}.png', dpi=100)
    
  if substance == '2-bromobutane':
    DensMAPENNR = df[(df["sort"] == 'Dens2SimMAPE')]
    #DensMAPEPR = dfPR[(dfPR["sort"] == 'Dens2SimMAPE')]
    RCEMAPENNR = df[(df["sort"] == 'RCE2MAPE')]
    #RCEMAPEPR = dfPR[(dfPR["sort"] == 'RCE2MAPE')]
    minRCE2BB = np.min([np.min(targets2BB), np.min(preOpt2BB), DensMAPENNR[["RCE21", "RCE22", "RCE23"]].min(numeric_only=True).min(), RCEMAPENNR[["RCE21", "RCE22", "RCE23"]].min(numeric_only=True).min()])
    minRCE2BB = minRCE2BB - abs(minRCE2BB*0.05)
    maxRCE2BB = np.max([np.max(targets2BB), np.max(preOpt2BB), DensMAPENNR[["RCE21", "RCE22", "RCE23"]].max(numeric_only=True).max(), RCEMAPENNR[["RCE21", "RCE22", "RCE23"]].max(numeric_only=True).max()])
    maxRCE2BB = maxRCE2BB + abs(maxRCE2BB*0.05)
    
    plt.figure(figsize=(5, 5))
    plt.plot([minRCE2BB, maxRCE2BB], [minRCE2BB, maxRCE2BB], '-', color='#000000')
    plt.plot(srtTargets2BB, srtTargets2BB, '-X', color='#000000', label="Targets")    
    thisLBL = False
    for i, row in DensMAPENNR.iterrows():
      #print(f'{row}')
      thisRCEs = np.array([row.RCE21, row.RCE22, row.RCE23])
      srtThisRCEs = []
      for idx in idx2BB:
        srtThisRCEs.append(thisRCEs[idx])
      if not thisLBL:
        plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#2ECC71', label=r'best MAPE$(\rho_{(2BB)})$')
        thisLBL = True
      else:
        plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#2ECC71')
    thisLBL = False
    for i, row in RCEMAPENNR.iterrows():
      #print(f'{row}')
      thisRCEs = np.array([row.RCE21, row.RCE22, row.RCE23])
      srtThisRCEs = []
      for idx in idx2BB:
        srtThisRCEs.append(thisRCEs[idx])
      if not thisLBL:
        plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#1B4F72', label=r'best MAPE$(\mathrm{RCE}_{(2BB)})$')
        thisLBL = True
      else:
        plt.plot(srtTargets2BB, srtThisRCEs, 'x', ls='dotted', color='#1B4F72')
    plt.plot(srtTargets2BB, srtpreOpt2BB, '-o', color='#ff0000', label="pre optimization")
    plt.xlabel('target RCE [kJ/mol]', fontsize=15, fontweight='bold')
    plt.ylabel('reproduced RCE [kJ/mol]', fontsize=15, fontweight='bold')
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.xlim([minRCE2BB, maxRCE2BB])
    plt.ylim([minRCE2BB, maxRCE2BB])
    plt.legend()
    plt.tight_layout()
 
    if saveOrShow == 'save':
      plt.savefig(f'RCE_targets-vs-preOpt-vs-Opt-{substance}.png', dpi=100)
  
  return plt


if __name__ == "__main__":
  pwd = os.path.dirname(os.getcwd())
  #print(f'{pwd}')
  try:
    targetFilePath = os.path.join(os.getcwd(), targetFileName)
  except:
    targetFilePath = None
  #print(f'{targetFilePath}')
  if not targetFilePath is None:
    RCEtargets = readRCEFile(targetFilePath)
    #print(f'{RCEtargets}')
  
  try:
    preOptFilePath = os.path.join(os.getcwd(), preOptFileName)
  except:
    preOptFilePath = None
  #print(f'{preOptFilePath}')
  if not preOptFilePath is None:
    preOptRCE = readRCEFile(preOptFilePath)
    #print(f'{preOptRCE}')
  
  dfSummaryBoth = pd.DataFrame({"method": [],
                                "sort": [],
                                "RCE1MAPE": [],
                                "RCE2MAPE": [],
                                "RCE1": [],
                                "RCE2": [],
                                "RCE3": [],
                                "RCE4": [],
                                "RCE5": [],
                                "RCE6": [],
                                "RCE7": [],
                                "RCE8": [],
                                "RCE9": [],
                                "RCE21": [],
                                "RCE22": [],
                                "RCE23": [],
                                "#": []})
  dfSummary1Bro = pd.DataFrame({"method": [],
                                "sort": [],
                                "RCE1MAPE": [],
                                "RCE2MAPE": [],
                                "RCE1": [],
                                "RCE2": [],
                                "RCE3": [],
                                "RCE4": [],
                                "RCE5": [],
                                "RCE6": [],
                                "RCE7": [],
                                "RCE8": [],
                                "RCE9": [],
                                "RCE21": [],
                                "RCE22": [],
                                "RCE23": [],
                                "#": []})
  dfSummary2Bro = pd.DataFrame({"method": [],
                                "sort": [],
                                "RCE1MAPE": [],
                                "RCE2MAPE": [],
                                "RCE1": [],
                                "RCE2": [],
                                "RCE3": [],
                                "RCE4": [],
                                "RCE5": [],
                                "RCE6": [],
                                "RCE7": [],
                                "RCE8": [],
                                "RCE9": [],
                                "RCE21": [],
                                "RCE22": [],
                                "RCE23": [],
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
          #exit()
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
  
  
  #plt = mkScatterPlot1(dfSummaryBoth, RCEtargets, preOptRCE, 'both')
  plt = mkScatterPlot2(dfSummary1Bro, RCEtargets, preOptRCE, '1-bromobutane')
  plt = mkScatterPlot2(dfSummary2Bro, RCEtargets, preOptRCE, '2-bromobutane')
  if saveOrShow == 'show':
    plt.show()
  
  
  
  
  
  
  
  
  
  
  
