import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

resultsFileName, sortedBy = "optResultsSortedDensityMAPEWithDensitySims.csv", "density" # sorted by the density MAPE
resultsFileName , sortedBy = "optResultsSortedRCEMAPEWithDensitySims.csv", "RCE"     # sorted by the RCE MAPE

plotRuns = []
plotRuns = ["minMAPE-6"]

preOptEnergiesFileName = "10_RCE_1-bromobutane.preOpt.csv"
targetEnergiesFileName = "10_1-bromobutane_RCE.target"
substance = "1-bromobutane"

colors = ["#577590", "#9bc53d", "#ffba40", "#ff8b94", "#ffee55"]

def CSVFile2DataFrame(csvFile):
  try:
    dfFromCSVFile = pd.read_csv(csvFile)
  except:
    print(f'Can not read CSV File \"{csvFile}\". Exit')
    exit()
  else:
    dfFromCSVFile = dfFromCSVFile.drop(dfFromCSVFile.columns[0], axis=1)  
  return dfFromCSVFile


if __name__ == "__main__":
  try:
    saveOrShow = sys.argv[1]
  except:
    saveOrShow = "show"
  if saveOrShow == "save":
    pass
  else:
    saveOrShow = "show"
    
  pwd = os.path.dirname(os.getcwd())
  
  resultsFile = os.path.join(pwd, resultsFileName)
  dfResults = CSVFile2DataFrame(resultsFile)
  #print(f'{dfResults}')
  
  targetEnergeisFile = os.path.join(os.getcwd(), targetEnergiesFileName)
  dfTargetEnergies = CSVFile2DataFrame(targetEnergeisFile)
  #print(f'{dfTargetEnergies}')
  sortIDX = dfTargetEnergies["RCE"].to_numpy().argsort()
  #print(f'{sortIDX}')
  
  preOptEnergeisFile = os.path.join(os.getcwd(), preOptEnergiesFileName)
  dfPreOptEnergies = CSVFile2DataFrame(preOptEnergeisFile)
  #print(f'{dfPreOptEnergies}')
  
  sortedTargets = []
  sortedPreOpts = []
  prevOptMAPE = 0
  for idx in sortIDX:
    thisTargetE = dfTargetEnergies["RCE"].iloc[idx]
    thisPreOptE = dfPreOptEnergies["RCE"].iloc[idx]
    sortedTargets.append(thisTargetE)
    sortedPreOpts.append(thisPreOptE)
    if thisTargetE == 0 and thisPreOptE == 0:
      pass
    elif thisTargetE == 0 and not thisPreOptE == 0:
      print(f'thisTargetE = 0, but thisPreOptE != 0. Set thisTargetE = 1e-6.')
      prevOptMAPE = prevOptMAPE + np.absolute((0.000001 - thisPreOptE)/0.000001)
    else:
      prevOptMAPE = prevOptMAPE + np.absolute((thisTargetE - thisPreOptE)/thisTargetE)
  prevOptMAPE = prevOptMAPE/len(sortIDX)
  #print(f'{sortedTargets}')
  #print(f'{sortedPreOpts}')
  #print(f'{prevOptMAPE*100:.2f}%')
  
  plt.plot(sortedTargets, sortedTargets, '-', color='#000000')
  plt.plot(sortedTargets, sortedPreOpts, 'o', linestyle='dotted', color='#D4282F', label=f'pre opt, MAPE: {prevOptMAPE*100:.2f}%')
  #plt.show()
  
  if len(plotRuns) > 0:
    for i in range(len(plotRuns)):
      dfThisOptRun = dfResults[dfResults["#"] == plotRuns[i]]
      #print(f'{dfThisOptRun}')
      
      thisSortedOptedE = []
      for idx in sortIDX:
        thisSortedOptedE.append(dfThisOptRun[f'RCE{idx+1}'])
      thisSortedOptedE = np.array(thisSortedOptedE)
      #print(f'{thisSortedOptedE}')
      thisRunName = dfThisOptRun["#"].to_string().split(" ")[-1]
      #print(f'{thisRunName}')
      thisRunMAPE = dfThisOptRun["RCEMAPE"].to_numpy()[0]
      #print(f'{thisRunMAPE}')
      plt.plot(sortedTargets, thisSortedOptedE, 'o', linestyle='dotted', color=colors[i], label=f'{thisRunName}, MAPE: {thisRunMAPE:.2f}%')
      #break
 
  plt.legend()  
  plt.xlabel(f'target rel. conf. energy [kcal/mol]')
  plt.ylabel(f'calculated rel. conf. energy [kcal/mol]')
  plt.suptitle(f'{substance}')
  plt.tight_layout()

  if saveOrShow == "save":
    figFilePath = os.path.join(pwd, f'renameThisFigurePLS.png')
    plt.savefig(figFilePath, dpi=300, format='png')
  if saveOrShow == "show":
    plt.show()


