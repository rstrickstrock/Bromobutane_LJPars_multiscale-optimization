import pandas as pd
import os
import glob
import shutil
import subprocess


optResultsFileName = "optResults.csv"
#dirs = ["30_1+2-Bromobutane_opt", "31_1-Bromobutane_opt", "31_2-Bromobutane_opt", "32_1+2-Bromobutane_opt_PR", "33_1-Bromobutane_opt_PR", "33_2-Bromobutane_opt_PR"]
#dirs = ["33_2-Bromobutane_opt_PR"]
dirs = ["31_2-Bromobutane_opt"]

def get_all_optimized_ff_parameters(optResultsFile):
  optResults = pd.read_csv(optResultsFile)
  optResults = optResults.drop(optResults.columns[0], axis=1)
  #print(f'{optResults}')
  optedFFPars = pd.DataFrame({"#": [],
                              "SigC800": [],
                              "SigBr730": [],
                              "EpsC800": [],
                              "EpsBr730": []})
  for i in range(len(optResults)):
    thisDFFFPars = pd.DataFrame({"#": [optResults.iloc[i]["#"]],
                                 "SigC800": [optResults.iloc[i]["SigC800"]],
                                 "SigBr730": [optResults.iloc[i]["SigBr730"]],
                                 "EpsC800": [optResults.iloc[i]["EpsC800"]],
                                 "EpsBr730": [optResults.iloc[i]["EpsBr730"]]})
    #print(f'{thisDFFFPars}')
    optedFFPars = pd.concat([optedFFPars, thisDFFFPars], ignore_index=True)
  #print(f'{optedFFPars}')
  return optedFFPars


def run_density_eval(pwd, evalDir, dfFFPars):
  #print(f'{pwd}')
  #print(f'{evalDir}')
  #print(f'{dfFFPars}')
  
  densityTemplateDir = os.path.join(pwd, "density.template")
  #print(f'{densityTemplateDir}')

  if not os.path.isdir(densityTemplateDir):
    print(f'template directory \"{densityTemplateDir}\" for density calculations does not exist. Exit.')
    exit()
  thisSimDir = os.path.join(evalDir, dfFFPars["#"])
  #print(f'{thisSimDir}')
  thisCommand = f'sh {os.path.join(densityTemplateDir, "physprop_setup.sh")} {thisSimDir} {dfFFPars["SigC800"]} {dfFFPars["SigBr730"]} {dfFFPars["EpsC800"]} {dfFFPars["SigBr730"]} {dfFFPars["#"]}'
  #print(f'{thisCommand}')
  p = subprocess.Popen(thisCommand, stdout=subprocess.PIPE, shell=True)
  (output, err) = p.communicate()
  p_status = p.wait()
  #print(f'Command output: {output}')



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
      optResultsFile = os.path.join(thisDir, optResultsFileName)
      if os.path.isfile(optResultsFile):
        ffPars = get_all_optimized_ff_parameters(optResultsFile)
      else:
        print(f'optResultsFile \"{optResultsFile}\" does not exist. Exit')
        exit()
      #print(f'{ffPars}')
      
      evalDir = os.path.join(thisDir, "evalDensities_withOptParams")
      #print(f'{evalDir}')
      if os.path.exists(evalDir):
        shutil.rmtree(evalDir)
      os.mkdir(evalDir)

      for i in range(len(ffPars)):
        #print(f'{ffPars.iloc[i]}')
        run_density_eval(os.getcwd(), evalDir, ffPars.iloc[i])

    
  





