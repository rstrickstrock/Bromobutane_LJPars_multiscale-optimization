import pandas as pd
import numpy as np
import os
import glob
from natsort import natsorted

dirs = ["30_1+2-Bromobutane_opt", "31_1-Bromobutane_opt", "31_2-Bromobutane_opt", "32_1+2-Bromobutane_opt_PR", "33_1-Bromobutane_opt_PR", "33_2-Bromobutane_opt_PR"]
#dirs = ["33_2-Bromobutane_opt_PR"]
#dirs = ["31_2-Bromobutane_opt"]

optResultsFileName = "optResults.csv"
optResultsSortedDensityFileName1 = "optResultsSortedDensityMAPE1.csv"
optResultsSortedRCEFileName1 = "optResultsSortedRCEMAPE1.csv"
optResultsSortedDensityFileName2 = "optResultsSortedDensityMAPE2.csv"
optResultsSortedRCEFileName2 = "optResultsSortedRCEMAPE2.csv"

densTargetFile = "00_bromobutane_density.target"
rceTargetFile = "00_bromobutane_RCE.target"


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


def get_last_opt_cycle_dir(optDir):
  """
  return directory names of the last optimization cycle
  """
  try:
    ppDirs = natsorted(glob.glob(os.path.join(optDir, 'PhysProp', 'g.*.0')))
  except:
    print(f'Optimization not started?! Could not get ppDirs for \"{optDir}\".')
    ppDir = None
  else:
    #print(f'\t\tget_last_opt_cycle_dir(): {ppDirs}')
    if len(ppDirs) >= 1:
      ppDir = ppDirs[-1]
    elif len(ppDirs) == 0:
      ppDir = None
  
  try:
    qmDirs = natsorted(glob.glob(os.path.join(optDir, 'QMMM', 'g.*.0')))
  except:
    print(f'Optimization not started?! Could not get qmDirs for \"{optDir}\".')
    qmDir = None
  else:
    #print(f'\t\tget_last_opt_cycle_dir(): {qmDirs}')
    if len(qmDirs) >= 1:
      qmDir = qmDirs[-1]
    elif len(qmDirs) == 0:
      qmDir = None
      
  #print(f'{ppDir}')
  #print(f'{qmDir}')
  #if ppDir is None and qmDir is None:
  #  print(f'No Optimization results to analyse. Exit.')
  #  exit()
  return ppDir, qmDir


def get_density_predictions(ppDir):
  """
  return density predictions
  """
  try:
    with open(os.path.join(ppDir, "this_prediction.csv"), 'r') as f:
      lines = f.readlines()
  except:
    print(f'Could not read density prediction file \"{os.path.join(ppDir, "this_prediction.csv")}\". Exit')
    exit()

  densities = []
  for line in lines:
    densities.append(float(line.split(" ")[1]))
  
  if len(densities) == 0:
    print(f'No density predictions found. Exit.')
    exit()
  return np.array(densities)


def get_RCEs(qmDir, substance):
  """
  return RCEs
  """
  try:
    with open(os.path.join(qmDir, "output", "Energy.rel.kcal.txt"), 'r') as f:
      lines = f.readlines()
  except:
    #print(f'Could not read RCE result file \"{os.path.join(qmDir, "output", "Energy.rel.kcal.txt")}\". Exit')
    if substance == 'both':
      return [f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None']
    elif substance == '1-bromobutane':
      return [f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None']
    elif substance == '2-bromobutane':
      return [f'None', f'None', f'None']
    else:
      print(f'Could not read RCE result file \"{os.path.join(qmDir, "output", "Energy.rel.kcal.txt")}\". Exit')
      exit()
    
  RCEs = []
  for line in lines:
    RCEs.append(float(line.split(" ")[1]))
    
  if len(RCEs) == 0:
    print(f'No RCE results found. Exit.')
    exit()
  return np.array(RCEs)


def get_ff_params(ppDir):
  """
  return optimized ff parameters
  """
  try:
    with open(os.path.join(ppDir, "this_parameters.csv"), 'r') as f:
      lines = f.readlines()
  except:
    print(f'Could not read this parameters file \"{os.path.join(ppDir, "this_parameters.csv")}\". Exit')
    exit()
    
  params = lines[1].split(",")
  ffParams = []
  for ffPar in params:
    ffParams.append(float(ffPar))
  return np.array(ffParams)


def calc_mape(targets, simulationResults):
  """
  return resulting MAPE of optimized property !! in [%] !!
  """
    
  if not type(targets) is type(simulationResults):
    if simulationResults is None or simulationResults == f'None':
      #print(f'simulationResults is \"None\". Set mape = \"None\".')
      return f'None'
    if isinstance(simulationResults , list) and simulationResults[0] is None or simulationResults[0] == f'None':
      #print(f'simulationResults is \"None\". Set mape = \"None\".')
      return f'None'
    print(f'type missmatch type(targets) = {type(targets)}, type(simulationResults) = {type(simulationResults)}. Exit')
    print(f'\t\tcalc_mape(): {targets}')
    print(f'\t\tcalc_mape(): {simulationResults}')
    exit()

  #print(f'\t\tcalc_mape(): {targets}')
  #print(f'\t\tcalc_mape(): {simulationResults}')
  if not isinstance(simulationResults, np.float64):
    if not len(targets) == len(simulationResults):
      print(f'len missmatch: len(targets) = {len(targets)}, len(simulationResults) = {len(simulationResults)}. Set mape = \"None\".')
      return f'None'
    
  MAPE = 0
  try:
    n = len(targets)
  except:
    n = 1
    try:
      targets = np.array(float(targets))
    except:
      print(f'do something about this error.')
  
  if n == 1:
    if targets == 0:
      if simulationResults == 0:
        pass
      else:
        print(f'target = {targets} and simulationResult = {simulationResult}. Added 0.000001 to each.')
        MAPE = MAPE + np.abs((targets - simulationResults)/0.000001)
    else:
      MAPE = MAPE + np.abs((targets - simulationResults)/targets)
  else:
    for i in range(n):
      #print(f'{targets}')
      if targets[i] == 0:
        if simulationResults[i] == 0:
          pass
        else:
          print(f'target = {targets[i]} and simulationResult = {simulationResult[i]}. Added 0.000001 to each.')
          MAPE = MAPE + np.abs((targets[i] - simulationResults[i])/0.000001)
      else:
        MAPE = MAPE + np.abs((targets[i] - simulationResults[i])/targets[i])
  MAPE = MAPE/n
  return MAPE*100

  
def get_results(optDirs, densTargets1, densTargets2, rceTargets1, rceTargets2, substance):
  if densTargets1 is None and densTargets2 is None:
    print(f'\"densTargets1\" is None AND \"densTargets2\" is None. Exit.')
    exit()
  if rceTargets1 is None and rceTargets2 is None:
    print(f'\"rceTargets1\" is None AND \"rceTargets2\" is None. Exit.')
    exit()
    
  dfResults = pd.DataFrame({"#": [],
                            "substance": [],
                            "SigC800": [],
                            "SigBr730": [],
                            "EpsC800": [],
                            "EpsBr730": [],
                            "MAPE_dens1": [],
                            "MAPE_dens2": [],
                            "MAPE_RCE1": [],
                            "MAPE_RCE2": [],
                            "prediction_1-bromobutane": [],
                            "prediction_2-bromobutane": [],
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
                            "RCE23": []})
  
  for optDir in optDirs:                              
    #print(f'\tget_results(): {optDir}')
    ppDir, qmDir = get_last_opt_cycle_dir(optDir)
    #print(f'\tget_results(): {ppDir}')
    #print(f'\tget_results(): {qmDir}')
    
    if not ppDir == None:
      thisDensityPredictions = get_density_predictions(ppDir)
    else:
      thisDensityPredictions = [f'None', f'None']
    #print(f'\tget_results(): {thisDensityPredictions}')
    
    if not qmDir == None:
      thisRCEs = get_RCEs(qmDir, substance)
    else:
      thisRCEs = [f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None']
    #print(f'\tget_results(): {thisRCEs}')
      
    if densTargets1 is None and rceTargets1 is None:
      #print(f'\tget_results(): 1 == None')
      thisDensityPredictions1 = f'None'
      thisDensityPredictions2 = thisDensityPredictions[0]
      thisRCEs1 = [f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None', f'None']
      thisRCEs2 = thisRCEs[0:2]
    elif densTargets2 is None and rceTargets2 is None:
      #print(f'\tget_results(): 2 == None')
      thisDensityPredictions1 = thisDensityPredictions[0]
      thisDensityPredictions2 = f'None'
      thisRCEs1 = thisRCEs[:9]
      thisRCEs2 = [f'None', f'None', f'None']
    else:
      #print(f'\tget_results(): both')
      thisDensityPredictions1 = thisDensityPredictions[0]
      thisDensityPredictions2 = thisDensityPredictions[1] 
      thisRCEs1 = thisRCEs[:9]
      thisRCEs2 = thisRCEs[9:]
    #print(f'\tget_results(): {thisDensityPredictions1}')
    #print(f'\tget_results(): {thisDensityPredictions1}')
    #print(f'\tget_results(): {thisRCEs1}')
    #print(f'\tget_results(): {thisRCEs2}')
    #print(f'\tget_results(): {ppDir}')
    #print(f'\tget_results(): {qmDir}')
    while len(thisRCEs2) < 3:
      #print(f'{len(thisRCEs2)}')
      thisRCEs2 = np.append(thisRCEs2, 0)
      
    if ppDir == None and qmDir == None:
      thisFFParams = [f'None', f'None', f'None', f'None']
    elif not ppDir is None:
      #print(f'\tget_results(): {ppDir}')
      thisFFParams = get_ff_params(ppDir)
    elif not qmDir is None:
      #print(f'\tget_results(): {qmDir}')
      thisFFParams = get_ff_params(qmDir)  
    #print(f'\tget_results(): {thisFFParams}')

    if not densTargets1 is None:
      thisMAPEDens1 = calc_mape(densTargets1, thisDensityPredictions1)
    else:
      thisMAPEDens1 = f'None'
    if not densTargets2 is None:
      thisMAPEDens2 = calc_mape(densTargets2, thisDensityPredictions2)
    else:
      thisMAPEDens2 = f'None'
    if not rceTargets1 is None:
      thisMAPERCEs1 = calc_mape(rceTargets1, thisRCEs1)
    else:
      thisMAPERCEs1 = f'None'
    if not rceTargets2 is None:
      try:
        thisMAPERCEs2 = calc_mape(rceTargets2, thisRCEs2)
      except:
        thisMAPERCEs2 = f'None'
    else:
      thisMAPERCEs2 = f'None'
    #print(f'\tget_results(): {thisMAPEDens1}')
    #print(f'\tget_results(): {thisMAPERCEs1}')
    #print(f'\tget_results(): {thisMAPEDens2}')
    #print(f'\tget_results(): {thisMAPERCEs2}')
    
    #print(f'{os.path.basename(optDir)}')
    #print(f'\tget_results(): {thisRCEs2}')
    dfThisResults = pd.DataFrame({"#": [os.path.basename(optDir)],
                                  "substance": [substance],
                                  "SigC800": [thisFFParams[0]],
                                  "SigBr730": [thisFFParams[1]],
                                  "EpsC800": [thisFFParams[2]],
                                  "EpsBr730": [thisFFParams[3]],
                                  "MAPE_dens1": [thisMAPEDens1],
                                  "MAPE_dens2": [thisMAPEDens2],
                                  "MAPE_RCE1": [thisMAPERCEs1],
                                  "MAPE_RCE2": [thisMAPERCEs2],
                                  "prediction_1-bromobutane": [thisDensityPredictions1],
                                  "prediction_2-bromobutane": [thisDensityPredictions2],
                                  "RCE1": [thisRCEs1[0]],
                                  "RCE2": [thisRCEs1[1]],
                                  "RCE3": [thisRCEs1[2]],
                                  "RCE4": [thisRCEs1[3]],
                                  "RCE5": [thisRCEs1[4]],
                                  "RCE6": [thisRCEs1[5]],
                                  "RCE7": [thisRCEs1[6]],
                                  "RCE8": [thisRCEs1[7]],
                                  "RCE9": [thisRCEs1[8]],
                                  "RCE21": [thisRCEs2[0]],
                                  "RCE22": [thisRCEs2[1]],
                                  "RCE23": [thisRCEs2[2]]})
    #print(f'{dfThisResults}')
    dfResults = pd.concat([dfResults, dfThisResults], ignore_index=True)
  #print(f'{dfResults}')
  return dfResults


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
      thisDirs = [os.path.join(thisDir, '1-bromo_initPars'), os.path.join(thisDir, '2-bromo_initPars')]
    else:
      thisDirs = [thisDir]
      
    for thisDir in thisDirs:
      #print(f'main(): {thisDir}')
      optDirsMaxR2 = natsorted(glob.glob(os.path.join(thisDir, 'maxR2', 'maxR2-*')))
      optDirsMinMAPE = natsorted(glob.glob(os.path.join(thisDir, 'minMAPE', 'minMAPE-*')))
      #print(f'main(): {optDirsMaxR2}')
      #print(f'main(): {optDirsMinMAPE}')
  
      densTargets = get_targets(densTargetFile)
      #print(f'main(): {densTargets}')
      rceTargets = get_targets(rceTargetFile)
      #print(f'main(): {rceTargets}')
      if substance == 'both':
        densTargets1 = densTargets[0]
        densTargets2 = densTargets[1]
        #print(f'main(): {densTargets1}')
        #print(f'main(): {densTargets2}')
        rceTargets1 = rceTargets[:9]
        rceTargets2 = rceTargets[9:]
        #print(f'main(): {rceTargets1}')
        #print(f'main(): {rceTargets2}')
        dfResultsMaxR2 = get_results(optDirsMaxR2, densTargets1, densTargets2, rceTargets1, rceTargets2, substance)
        dfResultsMinMAPE = get_results(optDirsMinMAPE, densTargets1, densTargets2, rceTargets1, rceTargets2, substance)

      elif substance == '1-bromobutane':
        densTargets1 = densTargets[0]
        #print(f'main(): {densTargets1}')
        rceTargets1 = rceTargets[:9]
        #print(f'main(): {rceTargets1}')
        dfResultsMaxR2 = get_results(optDirsMaxR2, densTargets1, None, rceTargets1, None, substance)
        dfResultsMinMAPE = get_results(optDirsMinMAPE, densTargets1, None, rceTargets1, None, substance)
      elif substance == '2-bromobutane':
        densTargets2 = densTargets[1]
        #print(f'main(): {densTargets2}')
        rceTargets2 = rceTargets[9:]
        #print(f'main(): {rceTargets2}')
        dfResultsMaxR2 = get_results(optDirsMaxR2, None, densTargets2, None, rceTargets2, substance)
        dfResultsMinMAPE = get_results(optDirsMinMAPE, None, densTargets2, None, rceTargets2, substance)
      #print(f'main(): {dfResultsMaxR2}')
      #print(f'main(): {dfResultsMinMAPE}')
   
      dfResults = pd.concat([dfResultsMaxR2, dfResultsMinMAPE], ignore_index=True)
      print(f'{dfResults}')
  
      optResultsFile = os.path.join(thisDir, optResultsFileName)
      #print(f'main(): {optResultsFile}')
      if os.path.exists(optResultsFile):
        os.remove(optResultsFile)
        print(f'Removed existing optimization results file: \'{optResultsFile}\'.')
      dfResults.to_csv(optResultsFile)
      print(f'Wrote optimization results summary to file: \'{optResultsFile}\'.')











