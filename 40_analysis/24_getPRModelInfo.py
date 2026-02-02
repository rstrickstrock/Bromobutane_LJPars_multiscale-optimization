import pandas as pd
import numpy as np
import os

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import pickle

dirs = ["30_1+2-Bromobutane_opt", "31_1-Bromobutane_opt", "31_2-Bromobutane_opt", "32_1+2-Bromobutane_opt_PR", "33_1-Bromobutane_opt_PR", "33_2-Bromobutane_opt_PR"]

optResultsFileName = "optResults_withDensitySims.csv"
top = 5

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


def getPolynInfo(df, substance):
  if substance == 'both':
    pathToModels1BB = f'/work/rstric2s/current_sim/bromobutane.app/01_density_ML-models/01_C4Br_C3Br/32_1-bromobutane_PR'
    pathToModels2BB = f'/work/rstric2s/current_sim/bromobutane.app/01_density_ML-models/01_C4Br_C3Br/32_2-bromobutane_PR'
  elif substance == '1-bromobutane':
    pathToModels1BB = f'/work/rstric2s/current_sim/bromobutane.app/01_density_ML-models/01_C4Br_C3Br/32_1-bromobutane_PR'
    pathToModels2BB = f'None'
  elif substance == '2-bromobutane':
    pathToModels1BB = f'None'
    pathToModels2BB = f'/work/rstric2s/current_sim/bromobutane.app/01_density_ML-models/01_C4Br_C3Br/32_2-bromobutane_PR'
  else:
    print(f'Check \"substance\". Exit.')
    exit()
  pathToModels = [pathToModels1BB, pathToModels2BB]
  
  PRruns = df[df["method"] == "PR"]["#"].unique()
  for run in PRruns:
    if "maxR2" in run:
      modelDir = f'bestTrainedModels_maxR2'
    elif "minMAPE" in run:
      modelDir = f'bestTrainedModels_minMAPE'
    
    for thisPath in pathToModels:
      if thisPath == f'None':
        pass
      else:
        thisBestModelsLogFile = os.path.join(thisPath, modelDir, f'bestModels.log')
        #print(f'{thisBestModelsLogFile}')
        #print(f'{run}')
        with open(thisBestModelsLogFile, 'r') as f:
          thisBestModelsLog = f.readlines()
          for line in thisBestModelsLog:
            if run in line:
              thisLine = line
              break
        #print(f'{line}')
        
        thisModelSpecs = line.split(" ")[0]
        thisModelName = f'{line.split(" ")[2].split(".")[0]}.sav'
        thisModelPath = os.path.join(thisPath, modelDir, thisModelName)
        #print(f'{thisModelSpecs}')
        #print(f'{thisModelName}')
        
        with open(thisModelPath, "rb") as thisModelInput:
          thisModel = pickle.load(thisModelInput)
        #print(f'{thisModel.named_steps}') 
        poly = thisModel.named_steps["polynomialfeatures"]
        linreg = thisModel.named_steps["linearregression"]

        degree = poly.degree
        coefs = linreg.coef_
        intercept = linreg.intercept_

        terms = poly.get_feature_names_out()

        print("Polynomial degree:", degree)
        print("Intercept:", intercept)

        for t, c in zip(terms, coefs):
            print(f"{t}: {c}")
        exit()

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
        if method == "PR":
          dfInfo = getPolynInfo(dfSummaryBoth, substance)
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
  
 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
