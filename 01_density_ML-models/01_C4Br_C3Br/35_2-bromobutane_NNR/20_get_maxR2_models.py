import pandas as pd
import shutil
import os

statisticsFile = 'Stats.csv'
numberOfModels = 20
modelDir = 'bestTrainedModels_maxR2'
logFile = os.path.join(modelDir, 'bestModels.log')

if not os.path.isfile(statisticsFile):
  print(f'Can not find and open \'{statisticsFile}\'. Exit.')
  exit()
else:
  dfStatistics = pd.read_csv(statisticsFile)
  #print(f'{dfStatistics}')
  try:
    dfStatistics = dfStatistics.drop(columns=["Unnamed: 0"])
  except:
    print(f'Something went wrong with\'dfStatistics = dfStatistics.drop(columns=["Unnamed: 0"])\'.')
  else:
    #print(f'{dfStatistics}')
    pass

if os.path.isdir(modelDir):
  shutil.rmtree(modelDir)
os.mkdir(modelDir) 
with open(logFile, "w") as f:
  f.write(f'Logfile for best models:\n')

for i in range(0, numberOfModels):
  idxMaxR2 = dfStatistics['r2'].idxmax() 
  maxR2Row = dfStatistics.loc[idxMaxR2]
  #print(f'{maxR2Row}')
  dfStatistics = dfStatistics.drop([idxMaxR2])
  
  thisRatio = maxR2Row.ratio
  if str(thisRatio).endswith("5"):
    thisRatio = f'{thisRatio:.2f}'
  else:
    thisRatio = f'{thisRatio:.1f}'
  #print(f'{thisRatio}')
  thisRndInt = str(int(maxR2Row.rndint))
  #print(f'{thisRndInt}')

  trainedModelFile = f'trainedModel_{thisRatio}_{thisRndInt}.pth'
  pathToModelFile = os.path.join("trainedModels", thisRatio, trainedModelFile)
  pathToModelDest = os.path.join(modelDir, f'2-bromobutane_maxR2-{i}.pth')
  if not os.path.isfile(pathToModelFile):
    print(f'trained model file \"{pathToModelFile}\" does not exist. Exit.')
    exit()
  else:
    print(f'copy \"{pathToModelFile}\" to \"{pathToModelDest}\"')
    with open(logFile, "a") as f:
      f.write(f'{trainedModelFile} -> 2-bromobutane_maxR2-{i}.pth\n')
    shutil.copyfile(pathToModelFile, pathToModelDest)
