import pandas as pd
import numpy as np
import os
import glob
from natsort import natsorted

#dirs = ["30_1+2-Bromobutane_opt", "31_1-Bromobutane_opt", "31_2-Bromobutane_opt", "32_1+2-Bromobutane_opt_PR", "33_1-Bromobutane_opt_PR", "33_2-Bromobutane_opt_PR"]
#dirs = ["33_2-Bromobutane_opt_PR"]
#dirs = ["31_2-Bromobutane_opt"]
dirs = ["30_1+2-Bromobutane_opt"]

optResultsFileName = "optResults.csv"
optResultsSortedDensityMAPEFileName = "optResults_sortedDensityMAPE.csv"
optResultsSortedRCEMAPEFileName = "optResults_sortedRCEMAPE.csv"


if __name__ == "__main__":
  pwd = os.path.dirname(os.getcwd())
  #print(f'{pwd}')
  
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
      thisDirs = [os.path.join(thisDir, '1-bromo_initPars'), os.path.join(thisDir, '2-bromo_initPars')]
    else:
      thisDirs = [thisDir]
      
    for thisDir in thisDirs:
      #print(f'main(): {thisDir}')
      thisOptResultsFile = os.path.join(thisDir, optResultsFileName)
      try:
        dfThisOptResults = pd.read_csv(thisOptResultsFile)
      except:
        print(f'Could not read {thisOptResultsFile}.')
      else:
        #print(f'{dfThisOptResults}')
        try:
          dfThisOptResults = dfThisOptResults.drop(dfThisOptResults.columns[0], axis=1)
        except:
          print(f'Could not drop column.')
          break
      
      #print(f'{dfThisOptResults["MAPE_dens1"]}')
      #print(f'{dfThisOptResults.sort_values("MAPE_dens1")}')
      #print(f'{substance}')
      if substance == 'both':
        if '1-bromo' in os.path.basename(thisDir):
          #print(f'{dfThisOptResults.sort_values("MAPE_dens1")}')
          dfThisOptResultsDensityMAPE = dfThisOptResults.sort_values("MAPE_dens1")
          dfThisOptResultsRCEMAPE = dfThisOptResults.sort_values("MAPE_RCE1")
          optResultsSortedDensityMAPEFileName = "optResults_sortedDensityMAPE_1-bromo.csv"
          optResultsSortedRCEMAPEFileName = "optResults_sortedRCEMAPE_1-bromo.csv"
        if '2-bromo' in os.path.basename(thisDir):
          #print(f'{dfThisOptResults.sort_values("MAPE_dens2")}')
          dfThisOptResultsDensityMAPE = dfThisOptResults.sort_values("MAPE_dens2")
          dfThisOptResultsRCEMAPE = dfThisOptResults.sort_values("MAPE_RCE2")
          optResultsSortedDensityMAPEFileName = "optResults_sortedDensityMAPE_2-bromo.csv"
          optResultsSortedRCEMAPEFileName = "optResults_sortedRCEMAPE_2-bromo.csv"
      elif substance == '1-bromobutane':
        #print(f'{dfThisOptResults.sort_values("MAPE_dens1")}')
        dfThisOptResultsDensityMAPE = dfThisOptResults.sort_values("MAPE_dens1")
        dfThisOptResultsRCEMAPE = dfThisOptResults.sort_values("MAPE_RCE1")
        optResultsSortedDensityMAPEFileName = "optResults_sortedDensityMAPE_1-bromo.csv"
        optResultsSortedRCEMAPEFileName = "optResults_sortedRCEMAPE_1-bromo.csv"
      elif substance == '2-bromobutane':
        #print(f'{dfThisOptResults.sort_values("MAPE_dens2")}')
        dfThisOptResultsDensityMAPE = dfThisOptResults.sort_values("MAPE_dens2")
        dfThisOptResultsRCEMAPE = dfThisOptResults.sort_values("MAPE_RCE2")
        optResultsSortedDensityMAPEFileName = "optResults_sortedDensityMAPE_2-bromo.csv"
        optResultsSortedRCEMAPEFileName = "optResults_sortedRCEMAPE_2-bromo.csv"
      
      dfThisSummaryDensity = dfThisOptResultsDensityMAPE[["#", "MAPE_dens1", "MAPE_dens2", "MAPE_RCE1", "MAPE_RCE2", "SigC800", "SigBr730", "EpsC800", "EpsBr730"]].copy()
      dfThisSummaryRCE = dfThisOptResultsRCEMAPE[["#", "MAPE_dens1", "MAPE_dens2", "MAPE_RCE1", "MAPE_RCE2", "SigC800", "SigBr730", "EpsC800", "EpsBr730"]].copy()
      #print(f'{dfThisSummaryDensity}')
      #print(f'{dfThisSummaryRCE}')
      
      thisOptResultsSortedDensityMAPEFileName = os.path.join(thisDir, optResultsSortedDensityMAPEFileName)
      thisOptResultsSortedRCEMAPEFileName = os.path.join(thisDir, optResultsSortedRCEMAPEFileName)
      #print(f'{thisOptResultsSortedDensityMAPEFileName}')
      #print(f'{thisOptResultsSortedRCEMAPEFileName}')
      
      if os.path.exists(thisOptResultsSortedDensityMAPEFileName):
        os.remove(thisOptResultsSortedDensityMAPEFileName)
        print(f'Removed existing optimization results file: \'{thisOptResultsSortedDensityMAPEFileName}\'.')
      dfThisSummaryDensity.to_csv(thisOptResultsSortedDensityMAPEFileName)
      print(f'Wrote optimization results summary to file: \'{thisOptResultsSortedDensityMAPEFileName}\'.')
      
      if os.path.exists(thisOptResultsSortedRCEMAPEFileName):
        os.remove(thisOptResultsSortedRCEMAPEFileName)
        print(f'Removed existing optimization results file: \'{thisOptResultsSortedRCEMAPEFileName}\'.')
      dfThisSummaryRCE.to_csv(thisOptResultsSortedRCEMAPEFileName)
      print(f'Wrote optimization results summary to file: \'{thisOptResultsSortedRCEMAPEFileName}\'.')
      print(f'')











