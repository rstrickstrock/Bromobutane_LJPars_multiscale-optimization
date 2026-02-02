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


def calcMAPE(npArray1, npArray2):
  """
  returns MAPE in [%]
  """
  MAPE = 0
  for i in range(len(npArray1)):
    #print(f'{i}')
    if npArray1[i] == 0:
      pass
    else:
      #print(f'old MAPE: {MAPE}')
      #print(f'add: abs(({npArray1[i]} - {npArray2[i]})/{npArray1[i]}) = abs({npArray1[i]-npArray2[i]}/{npArray1[i]}) = abs({(npArray1[i]-npArray2[i])/npArray1[i]})')
      MAPE = MAPE + abs((npArray1[i] - npArray2[i])/npArray1[i])
      #print(f'new MAPE: {MAPE}\n')
  return (MAPE*100)/(len(npArray1))
  
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
    #print(f'{RCEtargets["RCE"][0]}')

  try:
    preOptFilePath = os.path.join(os.getcwd(), preOptFileName)
  except:
    preOptFilePath = None
  #print(f'{preOptFilePath}')
  if not preOptFilePath is None:
    preOptRCE = readRCEFile(preOptFilePath)
    #print(f'{preOptRCE}')
  
  dfTargets = pd.DataFrame({"method": f'T',
                            "sort": f'None',
                            "RCE1MAPE": f'None',
                            "RCE2MAPE": f'None',
                            "RCE1": [RCEtargets["RCE"][0]],
                            "RCE2": [RCEtargets["RCE"][1]],
                            "RCE3": [RCEtargets["RCE"][2]],
                            "RCE4": [RCEtargets["RCE"][3]],
                            "RCE5": [RCEtargets["RCE"][4]],
                            "RCE6": [RCEtargets["RCE"][5]],
                            "RCE7": [RCEtargets["RCE"][6]],
                            "RCE8": [RCEtargets["RCE"][7]],
                            "RCE9": [RCEtargets["RCE"][8]],
                            "RCE21": [RCEtargets["RCE"][9]],
                            "RCE22": [RCEtargets["RCE"][10]],
                            "RCE23": [RCEtargets["RCE"][11]],
                            "#": f'None'})
  #print(f'{dfTargets}')                          
  targets1BB = RCEtargets["RCE"][:9].to_numpy()
  #print(f'{targets1BB}')
  targets2BB = RCEtargets["RCE"][9:].to_numpy()
  #print(f'{targets2BB}')
  preOpts1BB = preOptRCE["RCE"][:9].to_numpy()
  #print(f'{preOpts1BB}')
  preOpts2BB = preOptRCE["RCE"][9:].to_numpy()
  #print(f'{preOpts2BB}')
  preOpt1MAPE = calcMAPE(targets1BB, preOpts1BB)
  #print(f'{preOpt1MAPE}')
  preOpt2MAPE = calcMAPE(targets2BB, preOpts2BB)
  #print(f'{preOpt2MAPE}')

  dfPreOpts = pd.DataFrame({"method": f'preOpt',
                            "sort": f'None',
                            "RCE1MAPE": [preOpt1MAPE],
                            "RCE2MAPE": [preOpt2MAPE],
                            "RCE1": [preOptRCE["RCE"][0]],
                            "RCE2": [preOptRCE["RCE"][1]],
                            "RCE3": [preOptRCE["RCE"][2]],
                            "RCE4": [preOptRCE["RCE"][3]],
                            "RCE5": [preOptRCE["RCE"][4]],
                            "RCE6": [preOptRCE["RCE"][5]],
                            "RCE7": [preOptRCE["RCE"][6]],
                            "RCE8": [preOptRCE["RCE"][7]],
                            "RCE9": [preOptRCE["RCE"][8]],
                            "RCE21": [preOptRCE["RCE"][9]],
                            "RCE22": [preOptRCE["RCE"][10]],
                            "RCE23": [preOptRCE["RCE"][11]],
                            "#": f'None'})
  #print(f'{dfPreOpts}')
  #exit()
  
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
  dfSummaryBoth = pd.concat([dfSummaryBoth, dfTargets], ignore_index=True)
  dfSummaryBoth = pd.concat([dfSummaryBoth, dfPreOpts], ignore_index=True)
  
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
  dfSummary1Bro = pd.concat([dfSummary1Bro, dfTargets], ignore_index=True)
  dfSummary1Bro = pd.concat([dfSummary1Bro, dfPreOpts], ignore_index=True)
  
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
  dfSummary2Bro = pd.concat([dfSummary2Bro, dfTargets], ignore_index=True)
  dfSummary2Bro = pd.concat([dfSummary2Bro, dfPreOpts], ignore_index=True)
  
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
  if False:
    for i, row in dfSummaryBoth.iterrows():
      if row.sort == "Dens1SimMAPE":
        #print(f'1')
        thisSorting = '$\\rho_{(1\\mathrm{BB})}$'
      elif row.sort == "Dens2SimMAPE":
        #print(f'2')
        thisSorting = '$\\rho_{(2\\mathrm{BB})}$'
      else:
        #print(f'3')
        thisSorting = 'miau'
      
      if thisSorting == 'miau':
        #print(f'4')
        if row.method == "T":
          print(f'    {row.RCE1:.2f} & {row.RCE2:.2f} & {row.RCE3:.2f} & {row.RCE4:.2f} & {row.RCE5:.2f} & {row.RCE6:.2f} & {row.RCE7:.2f} & {row.RCE8:.2f} & {row.RCE9:.2f} & - & - & - \\\\')
        elif row.method == "preOpt":
          print(f'    {row.RCE1:.2f} & {row.RCE2:.2f} & {row.RCE3:.2f} & {row.RCE4:.2f} & {row.RCE5:.2f} & {row.RCE6:.2f} & {row.RCE7:.2f} & {row.RCE8:.2f} & {row.RCE9:.2f} & {row.RCE1MAPE:.2f} & - & - \\\\')
        else:
          pass
      else:
        #print(f'5')
          print(f'    {row.RCE1:.2f} & {row.RCE2:.2f} & {row.RCE3:.2f} & {row.RCE4:.2f} & {row.RCE5:.2f} & {row.RCE6:.2f} & {row.RCE7:.2f} & {row.RCE8:.2f} & {row.RCE9:.2f} & {row.RCE1MAPE:.2f} & {row.method} & {thisSorting} \\\\')
      
    for i, row in dfSummaryBoth.iterrows():
      if row.sort == "Dens1SimMAPE":
        #print(f'1')
        thisSorting = '$\\rho_{(1\\mathrm{BB})}$'
      elif row.sort == "Dens2SimMAPE":
        #print(f'2')
        thisSorting = '$\\rho_{(2\\mathrm{BB})}$'
      else:
        #print(f'3')
        thisSorting = 'miau'
      
      if thisSorting == 'miau':
        #print(f'4')
        if row.method == "T":
          print(f'    {row.RCE21:.2f} & {row.RCE22:.2f} & {row.RCE23:.2f} & - & - & - \\\\')
        elif row.method == "preOpt":
          print(f'    {row.RCE21:.2f} & {row.RCE22:.2f} & {row.RCE23:.2f} & {row.RCE2MAPE:.2f} & - & - \\\\')
        else:
          pass
      else:
        #print(f'5')
        print(f'    {row.RCE21:.2f} & {row.RCE22:.2f} & {row.RCE23:.2f} & {row.RCE2MAPE:.2f} & {row.method} & {thisSorting} \\\\') 
             
  #print(f'{dfSummary1Bro}')
  if False: 
    for i, row in dfSummary1Bro.iterrows():
      if row.sort == "Dens1SimMAPE":
        #print(f'1')
        thisSorting = '$\\rho_{(1\\mathrm{BB})}$'
      elif row.sort == "RCE1MAPE":
        #print(f'2')
        thisSorting = '$\\mathrm{RCE}_{(1\\mathrm{BB})}$'
      else:
        #print(f'3')
        thisSorting = 'miau'
      
      if thisSorting == 'miau':
        #print(f'4')
        if row.method == "T":
          print(f'    {row.RCE1:.2f} & {row.RCE2:.2f} & {row.RCE3:.2f} & {row.RCE4:.2f} & {row.RCE5:.2f} & {row.RCE6:.2f} & {row.RCE7:.2f} & {row.RCE8:.2f} & {row.RCE9:.2f} & - & - & - \\\\')
        elif row.method == "preOpt":
          print(f'    {row.RCE1:.2f} & {row.RCE2:.2f} & {row.RCE3:.2f} & {row.RCE4:.2f} & {row.RCE5:.2f} & {row.RCE6:.2f} & {row.RCE7:.2f} & {row.RCE8:.2f} & {row.RCE9:.2f} & {row.RCE1MAPE:.2f} & - & - \\\\')
        else:
          pass
      else:
        #print(f'5')
          print(f'    {row.RCE1:.2f} & {row.RCE2:.2f} & {row.RCE3:.2f} & {row.RCE4:.2f} & {row.RCE5:.2f} & {row.RCE6:.2f} & {row.RCE7:.2f} & {row.RCE8:.2f} & {row.RCE9:.2f} & {row.RCE1MAPE:.2f} & {row.method} & {thisSorting} \\\\')
  
  #print(f'{dfSummary2Bro}')
  if True:
    for i, row in dfSummary2Bro.iterrows():
      if row.sort == "Dens2SimMAPE":
        #print(f'1')
        thisSorting = '$\\rho_{(2\\mathrm{BB})}$'
      elif row.sort == "RCE2MAPE":
        #print(f'2')
        thisSorting = '$\\mathrm{RCE}_{(2\\mathrm{BB})}$'
      else:
        #print(f'3')
        thisSorting = 'miau'
      
      if thisSorting == 'miau':
        #print(f'4')
        if row.method == "T":
          print(f'    {row.RCE21:.2f} & {row.RCE22:.2f} & {row.RCE23:.2f} & - & - & - \\\\')
        elif row.method == "preOpt":
          print(f'    {row.RCE21:.2f} & {row.RCE22:.2f} & {row.RCE23:.2f} & {row.RCE2MAPE:.2f} & - & - \\\\')
        else:
          pass
      else:
        #print(f'5')
        print(f'    {row.RCE21:.2f} & {row.RCE22:.2f} & {row.RCE23:.2f} & {row.RCE2MAPE:.2f} & {row.method} & {thisSorting} \\\\') 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
