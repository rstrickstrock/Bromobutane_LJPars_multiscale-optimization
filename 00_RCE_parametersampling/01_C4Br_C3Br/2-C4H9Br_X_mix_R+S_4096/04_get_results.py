import os
import pandas as pd
import glob
import numpy as np

resultsFileName = '2-Bromobutane_RCE.csv'
targetEnergiesFile = '00_2-bromobutane_target_RCE.csv'
cwd = os.getcwd()
#print(f'{cwd}')
simdirs = glob.glob(os.path.join(cwd, "experiments", "*"))

dfTargetEnergies = pd.read_csv(targetEnergiesFile)
#print(f'{dfTargetEnergies.RCE[0]}')

dfTargets = pd.DataFrame({'SigC800': [0],
                          'SigBr730': [0],
                          'EpsC800': [0],
                          'EpsBr730': [0],
                          'RCE1': [dfTargetEnergies.RCE[0]],
                          'RCE2': [dfTargetEnergies.RCE[1]],
                          'RCE3': [dfTargetEnergies.RCE[2]],
                          'AvgMAPE': [0]})
#print(f'{dfTargets}')
targets = [float(dfTargetEnergies.RCE[0]), float(dfTargetEnergies.RCE[1]), float(dfTargetEnergies.RCE[2])]
#print(f'{targets}')

dfResults = pd.DataFrame({'SigC800': [],
                          'SigBr730': [],
                          'EpsC800': [],
                          'EpsBr730': [],
                          'RCE1': [],
                          'RCE2': [],
                          'RCE3': [],
                          'AvgMAPE': []})

def calcAVG(targets, thisRCEs):
  if not len(targets) == len(thisRCEs):
    print(f'len(targets) [{len(targets)}] != len(thisRCEs) [{len(thisRCEs)}]. Exit()')
    exit()
  thisAvgMAPE = 0  
  for i in range(len(targets)):
    thisTarget = targets[i]
    thisRCE = thisRCEs[i]
    if thisTarget == 0.0:
      ## dont divide by 0
      thisTarget = 0.01
      thisRCE = 0.01
    thisAvgMAPE = thisAvgMAPE + np.absolute((thisTarget - thisRCE)/thisTarget)
  thisAvgMAPE = thisAvgMAPE/len(targets)
  return thisAvgMAPE


for simdir in simdirs:
  ffPars = os.path.basename(simdir).split("_")
  #print(f'{ffPars}')

  resultFile = os.path.join(simdir, "output", "Energy.rel.kcal.txt")
  #print(f'{resultFile}')
  if os.path.isfile(resultFile):
    with open(resultFile, "r") as f:
      lines = f.readlines() 
    #print(f'{lines}')
    thisRCEs = []
    for line in lines:
      #print(f'{line}')
      thisRCEs.append(float(line.split(" ")[1]))
    #print(f'{thisRCEs}')
    thisAvgMAPE = calcAVG(targets, thisRCEs)
    
    dfThisResult = pd.DataFrame({'SigC800': [ffPars[0]],
                                 'SigBr730': [ffPars[1]],
                                 'EpsC800': [ffPars[2]],
                                 'EpsBr730': [ffPars[3]],
                                 'RCE1': [thisRCEs[0]],
                                 'RCE2': [thisRCEs[1]],
                                 'RCE3': [thisRCEs[2]],
                                 'AvgMAPE': [thisAvgMAPE]})
   # print(f'{dfThisResult}')
    dfResults = pd.concat([dfResults, dfThisResult], ignore_index=True)
  
  
if os.path.exists(resultsFileName):
  os.remove(resultsFileName)
  print(f'Removed existing statistics file: \'{resultsFileName}\'.')
dfResults.to_csv(resultsFileName)
print(f'Wrote results to file: \'{resultsFileName}\'.')
print(f'{dfResults}')
