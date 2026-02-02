import os
import glob

startedSimulationsFile = "started_simulations.txt"
batchRerunSimsFile = "04_batch_rerunSims.sh"
cwd = os.getcwd()
#print(f'{cwd}')
densdirs = glob.glob(os.path.join(cwd, "experiments", "*"))

f = open(startedSimulationsFile, "r")
lines = f.readlines()
f.close()
nSims = len(lines)
#print(f'{nSims}')

nDens = 0
nFinished = 0
rerunSims = 0
batchCommands = []
for densdir in densdirs:
  densdir = os.path.join(densdir, "density")
  #print(f'{densdir}')
  slurmFiles = glob.glob(os.path.join(densdir, "slurm-*"))
  #print(f'{slurmFiles}')
  #print(f'{len(slurmFiles)}')
  if os.path.isfile(os.path.join(densdir, "density.xvg")):
    ## density already calculated
    nDens = nDens + 1
    nFinished = nFinished + 1
  elif len(slurmFiles) >= 4:
    ## density not calculated, but MD sim finished
    nFinished = nFinished + 1
  else:
    ## MD sim not finished, rerun
    rerunSims = rerunSims + 1
    rerunSimsFile = f'04_rerun_sims-{rerunSims}.sh'
    batchCommands.append(f'/bin/bash {rerunSimsFile}\n')
    command = f'rm {os.path.join(densdir, "1*")} && rm {os.path.join(densdir, "slurm-")} && cd {densdir} && python run_sim_chain.py && cd {cwd}'
    #print(f'{command}')
    if os.path.isfile(rerunSimsFile):
      os.remove(rerunSimsFile)
    with open(rerunSimsFile, 'w') as f:
      f.write(f'{command}')


if rerunSims > 0:
  if os.path.isfile(batchRerunSimsFile):
    os.remove(batchRerunSimsFile)
  with open(batchRerunSimsFile, "w") as f:
    f.write(f'#!/bin/bash\n')
    f.write(f'#SBATCH --partition=hpc,hpc1,hpc3\n')
    f.write(f'#SBATCH --nodes=1\n')
    f.write(f'#SBATCH --mem 10G\n')
    f.write(f'#SBATCH --time=72:00:00\n')
    f.write(f'#SBATCH --job-name=rerunScript\n\n')
    for command in batchCommands:
      f.write(command)
  
print(f'number of started sims: {nSims}')
print(f'number of finished sims: {nFinished}')
print(f'number of calculated densities: {nDens}')
if nDens == nSims:
  print(f'  -> run 04_get_results.py')
if rerunSims > 0:
  print(f'number of not finished sims: {rerunSims}')
  print(f'  -> run {batchRerunSimsFile}.')
