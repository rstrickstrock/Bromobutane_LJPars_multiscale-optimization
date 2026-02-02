import os
import glob
import pandas as pd
import numpy as np
from natsort import natsorted

optResultsFileName = "optResults.csv"
densTargetFile = "00_bromobutane_density.target"
evalDir = "evalDensities_withOptParams"

#dirs = ["30_1+2-Bromobutane_opt", "31_1-Bromobutane_opt", "31_2-Bromobutane_opt", "32_1+2-Bromobutane_opt_PR", "33_1-Bromobutane_opt_PR", "33_2-Bromobutane_opt_PR"]
dirs = ["31_2-Bromobutane_opt", "33_2-Bromobutane_opt_PR"]

newOptResultsFileName = "optResults_withDensitySims.csv"


def  get_sim_results(slurmFile, propertyName="Density"):
  try:
    with open(slurmFile, 'r') as f:
      lines = f.readlines()
  except:
    print(f'Could not read slurmFile \"{slurmFile}\". Exit.')
    exit()
  
  simRes = []
  simErr = []  
  for line in lines:
    if line.startswith(propertyName):
      #print(f'{line}')
      try:
        simRes.append(float(line.split()[1]))
      except:
        simRes.append(f'None')
        #print(f'Failed to execute \"simRes = float(line.split()[1])\". Exit.')
        #print(f'Check slurmFile \"{slurmFile}\"')
        #exit()
        
      try:
        simErr.append(float(line.split()[2]))
      except:
        simErr.append(f'None')
        #print(f'Failed to execute \"simErr = float(line.split()[2])\". Exit.')
        #print(f'Check slurmFile \"{slurmFile}\"')
        #exit()
      #print(f'{simRes}')
      #print(f'{simErr}')
  if len(simRes) > 0:
    return simRes, simErr, True
  else:
    return False, None, False


def get_targets(targetFile):
  """
  return optimization targets from provided target file
  """
  if targetFile is None:
    return None
    
  try:
    with open(targetFile, 'r') as f:
      lines = f.readlines()
  except:
    print(f'Could not open targetFile: \"{targetFile}\". Exit.')
    exit()
  else:
    #print(f'{lines}')
    targets = []
    for line in lines:
      targets.append(float(line.split(" ")[1]))
    #print(f'{targets}')    
    return np.array(targets)
  

def get_opt_results(optResultsFile):
  try:
    dfOptResults = pd.read_csv(optResultsFile)
  except:
    print(f'Can not read optResultsFile \"{optResultsFile}\". Exit')
    exit()
  else:
    dfOptResults = dfOptResults.drop(dfOptResults.columns[0], axis=1)
  
  n = len(dfOptResults["MAPE_dens1"])
  dfThisOptResults = pd.DataFrame({"Dens1SimMAPE": [None] * n,
                                   "Dens2SimMAPE": [None] * n,
                                   "RCE1MAPE": dfOptResults["MAPE_RCE1"],
                                   "RCE2MAPE": dfOptResults["MAPE_RCE2"],
                                   "SigC800": dfOptResults["SigC800"],
                                   "SigBr730": dfOptResults["SigBr730"],
                                   "EpsC800": dfOptResults["EpsC800"],
                                   "EpsBr730": dfOptResults["EpsBr730"],
                                   "Dens1PredMAPE": dfOptResults["MAPE_dens1"],
                                   "Dens1Pred": dfOptResults["prediction_1-bromobutane"],
                                   "Dens2PredMAPE": dfOptResults["MAPE_dens2"],
                                   "Dens2Pred": dfOptResults["prediction_2-bromobutane"],
                                   "Dens1Sim": [None] * n,
                                   "Dens1SimErr": [None] * n,
                                   "Dens2Sim": [None] * n,
                                   "Dens2SimErr": [None] * n,
                                   "RCE1": dfOptResults["RCE1"],
                                   "RCE2": dfOptResults["RCE2"],
                                   "RCE3": dfOptResults["RCE3"],
                                   "RCE4": dfOptResults["RCE4"],
                                   "RCE5": dfOptResults["RCE5"],
                                   "RCE6": dfOptResults["RCE6"],
                                   "RCE7": dfOptResults["RCE7"],
                                   "RCE8": dfOptResults["RCE8"],
                                   "RCE9": dfOptResults["RCE9"],
                                   "RCE21": dfOptResults["RCE21"],
                                   "RCE22": dfOptResults["RCE22"],
                                   "RCE23": dfOptResults["RCE23"],
                                   "#": dfOptResults["#"]})

  return dfThisOptResults


def calc_mape(targets, simulationResults):
  """
  return resulting MAPE of optimized property !! in [%] !!
  """
  if simulationResults is None or simulationResults == f'None':
    return f'None'
    
  try:
    n = len(targets)
  except:
    try:
      target = float(targets)
    except:
      print(f'Handle this error (1). Exit.')
      exit()
    else:
      targets = []
      targets.append(target)
  
  try:
    n = len(simulationResults)
  except:
    try:
      simulationResult = float(simulationResults)
    except:
      print(f'Handle this error (2). Exit.')
      exit()
    else:
      simulationResults = []
      simulationResults.append(simulationResult)
      
  if not len(targets) is len(simulationResults):
    #print(f'Dimension missmatch len(targets) = {len(targets)}, len(simulationResults) = {len(simulationResults)}. Skip')
    return f'None'
    
  MAPE = 0
  for i in range(len(targets)):
    #print(f'{targets[i]}')
    if targets[i] == 0:
      if simulationResults[i] == 0:
        pass
      else:
        print(f'target = {targets[i]} and simulationResult = {simulationResult[i]}. Added 0.000001 to each.')
        MAPE = MAPE + np.abs((targets[i] - simulationResults[i])/0.000001)
    else:
      #print(f'{targets[i]}')
      #print(f'{simulationResults[i]}')
      MAPE = MAPE + np.abs((targets[i] - simulationResults[i])/targets[i])
  MAPE = MAPE/len(targets)
  return MAPE*100


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
      thisDirs = [os.path.join(thisDir, '1-bromo_initPars', evalDir), os.path.join(thisDir, '2-bromo_initPars', evalDir)]
    else:
      thisDirs = [os.path.join(thisDir, evalDir)]
    
    for thisDir in thisDirs:
      #print(f'main(): {thisDir}')
      optResultsFile = os.path.join(os.path.dirname(thisDir), optResultsFileName)
      #print(f'main(): {optResultsFile}')
      
      dfOptResults = get_opt_results(optResultsFile)
      #print(f'{dfOptResults}')

      densTargets = get_targets(densTargetFile)
      #print(f'{densTargets}')
      if substance == "both":
        densTargets1 = densTargets[0]
        densTargets2 = densTargets[1]
      elif substance == "1-bromobutane":
        densTargets1 = densTargets[0]
        densTargets2 = f'None'
      elif substance == "2-bromobutane":
        densTargets1 = f'None'
        densTargets2 = densTargets[1]
      #print(f'{densTargets1}')
      #print(f'{densTargets2}')
      
      simDirs = glob.glob(os.path.join(thisDir, "*"))
      #print(f'{simDirs}')

      for simDir in simDirs:
        slurmFiles = natsorted(glob.glob(os.path.join(simDir, "slurm-*")))
        #print(f'{slurmFiles}')
        for slurmFile in slurmFiles:
          simRes, simErr, finish = get_sim_results(slurmFile)
          #print(f'{simRes}, {simErr}, {finish}')
          if finish:
            simRes1 = simRes[0]
            simRes2 = simRes[1]
            simErr1 = simErr[0]
            simErr2 = simErr[1]
            simRun = os.path.basename(os.path.dirname(slurmFile))
            break
            
        if not simRes:
          print(f'Could not get simRes for \"{simDir}\". Set to \"None\".')
          simRes1 = f'None'
          simRes2 = f'None'
          simErr1 = f'None'
          simErr2 = f'None'
          
        simMAPE1 = calc_mape(densTargets1, simRes1)
        simMAPE2 = calc_mape(densTargets2, simRes2)
          
        idx = dfOptResults.index[dfOptResults["#"] == simRun].tolist()[0]
        #print(f'{idx}')
        dfOptResults["Dens1SimMAPE"].iloc[idx] = simMAPE1
        dfOptResults["Dens1Sim"].iloc[idx] = simRes1
        dfOptResults["Dens1SimErr"].iloc[idx] = simErr1
          
        dfOptResults["Dens2SimMAPE"].iloc[idx] = simMAPE2
        dfOptResults["Dens2Sim"].iloc[idx] = simRes2
        dfOptResults["Dens2SimErr"].iloc[idx] = simErr2
      
      if substance == 'both':
        if '1-bromo' in os.path.basename(os.path.dirname(thisDir)):
          #print(f'{dfThisOptResults.sort_values("MAPE_dens1")}')
          dfThisOptResultsDensityMAPE = dfOptResults.sort_values("Dens1SimMAPE")
          dfThisOptResultsRCEMAPE = dfOptResults.sort_values("RCE1MAPE")
          newOptResultsSortedDensityMAPEFileName = "optResults_sortedDensityMAPE_1-bromo_withDensitySims.csv"
          newOptResultsSortedRCEMAPEFileName = "optResults_sortedRCEMAPE_1-bromo_withDensitySims.csv"
        if '2-bromo' in os.path.basename(os.path.dirname(thisDir)):
          #print(f'{dfThisOptResults.sort_values("MAPE_dens2")}')
          dfThisOptResultsDensityMAPE = dfOptResults.sort_values("Dens2SimMAPE")
          dfThisOptResultsRCEMAPE = dfOptResults.sort_values("RCE2MAPE")
          newOptResultsSortedDensityMAPEFileName = "optResults_sortedDensityMAPE_2-bromo_withDensitySims.csv"
          newOptResultsSortedRCEMAPEFileName = "optResults_sortedRCEMAPE_2-bromo_withDensitySims.csv"
      elif substance == '1-bromobutane':
        #print(f'{dfThisOptResults.sort_values("MAPE_dens1")}')
        dfThisOptResultsDensityMAPE = dfOptResults.sort_values("Dens1SimMAPE")
        dfThisOptResultsRCEMAPE = dfOptResults.sort_values("RCE1MAPE")
        newOptResultsSortedDensityMAPEFileName = "optResults_sortedDensityMAPE_1-bromo_withDensitySims.csv"
        newOptResultsSortedRCEMAPEFileName = "optResults_sortedRCEMAPE_1-bromo_withDensitySims.csv"
      elif substance == '2-bromobutane':
        #print(f'{dfThisOptResults.sort_values("MAPE_dens2")}')
        dfThisOptResultsDensityMAPE = dfOptResults.sort_values("Dens2SimMAPE")
        dfThisOptResultsRCEMAPE = dfOptResults.sort_values("RCE2MAPE")
        newOptResultsSortedDensityMAPEFileName = "optResults_sortedDensityMAPE_2-bromo_withDensitySims.csv"
        newOptResultsSortedRCEMAPEFileName = "optResults_sortedRCEMAPE_2-bromo_withDensitySims.csv"
      
      newOptResultsFile = os.path.join(os.path.dirname(thisDir), newOptResultsFileName)
      #print(f'{newOptResultsFile}')
      if os.path.exists(newOptResultsFile):
        os.remove(newOptResultsFile)
        print(f'Removed existing optimization results file: \'{newOptResultsFile}\'.')
      dfOptResults.to_csv(newOptResultsFile)
      print(f'Wrote optimization results summary to file: \'{newOptResultsFile}\'.')
      print(f'{dfOptResults}')
      
      #print(f'{newOptResultsSortedDensityMAPEFileName}')
      thisNewOptResultsSortedDensityMAPEFileName = os.path.join(os.path.dirname(thisDir), newOptResultsSortedDensityMAPEFileName)
      thisNewOptResultsSortedRCEMAPEFileName = os.path.join(os.path.dirname(thisDir), newOptResultsSortedRCEMAPEFileName)
      #print(f'{thisNewOptResultsSortedDensityMAPEFileName}')
      #print(f'{thisNewOptResultsSortedRCEMAPEFileName}')
      
      if os.path.exists(thisNewOptResultsSortedDensityMAPEFileName):
        os.remove(thisNewOptResultsSortedDensityMAPEFileName)
        print(f'Removed existing optimization results file: \'{thisNewOptResultsSortedDensityMAPEFileName}\'.')
      dfThisOptResultsDensityMAPE.to_csv(thisNewOptResultsSortedDensityMAPEFileName)
      print(f'Wrote optimization results summary to file: \'{thisNewOptResultsSortedDensityMAPEFileName}\'.')
      
      if os.path.exists(thisNewOptResultsSortedRCEMAPEFileName):
        os.remove(thisNewOptResultsSortedRCEMAPEFileName)
        print(f'Removed existing optimization results file: \'{thisNewOptResultsSortedRCEMAPEFileName}\'.')
      dfThisOptResultsRCEMAPE.to_csv(thisNewOptResultsSortedRCEMAPEFileName)
      print(f'Wrote optimization results summary to file: \'{thisNewOptResultsSortedRCEMAPEFileName}\'.')
      print(f'')
        
  
      
      
      
      

    
    
    
    
