import os
import pandas as pd
import glob

resultsFileName = '1-Bromobutane_Trainingsdata.tmp.csv'
cwd = os.getcwd()
#print(f'{cwd}')
densdirs = glob.glob(os.path.join(cwd, "experiments", "*"))

dfResults = pd.DataFrame({'SigC800': [],
                          'SigBr730': [],
                          'EpsC800': [],
                          'EpsBr730': [],
                          'density': [],
                          'density_err': []})

for densdir in densdirs:
  #print(f'{os.path.basename(densdir)}')
  ffPars = os.path.basename(densdir).split("_")
  #print(f'{ffPars}')

  densdir = os.path.join(densdir, "density")
  #print(f'{densdir}')
  slurmFiles = glob.glob(os.path.join(densdir, "slurm-*"))
  #print(f'{slurmFiles}')
  #print(f'{slurmFiles[-1]}')
  densityLine = []
  if len(slurmFiles) >= 1:
    lastSlurmFile = slurmFiles[-1]
    with open(lastSlurmFile, "r") as f:
      lines = f.readlines() 
    densityLine = [line.strip() for line in lines if line.lstrip().startswith("Density")]
    #print(f'{densityLine}')
  
  #print(f'{densityLine}')
  if len(densityLine) < 1:
    ## no density found - probably not all sims finished yet
    pass
  else:
    Results = []
    for item in densityLine[0].split(" "):
      if len(item) > 1:
        Results.append(item)
    #print(f'{Results}')
  
    dfThisEntry = pd.DataFrame({'SigC800': [ffPars[0]],
                                'SigBr730': [ffPars[1]],
                                'EpsC800': [ffPars[2]],
                                'EpsBr730': [ffPars[3]],
                                'density': [Results[1]],
                                'density_err': [Results[2]]})
    dfResults = pd.concat([dfResults, dfThisEntry], ignore_index=True)
  
if os.path.exists(resultsFileName):
  os.remove(resultsFileName)
  print(f'Removed existing statistics file: \'{resultsFileName}\'.')
dfResults.to_csv(resultsFileName)
print(f'Wrote results to file: \'{resultsFileName}\'.')
print(f'{dfResults}')
