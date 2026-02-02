import os
import glob

startedSimulationsFile = "started_simulations.txt"
batchRerunSimsFile = "04_batch_rerunSims.sh"
cwd = os.getcwd()
#print(f'{cwd}')
simdirs = glob.glob(os.path.join(cwd, "experiments", "*"))

f = open(startedSimulationsFile, "r")
lines = f.readlines()
f.close()
nSims = len(lines)
#print(f'{nSims}')

nRCE = 0
nFinished = 0
rerunSims = 0
batchCommands = []
for simdir in simdirs:
  enerFile = os.path.join(simdir, "output", "Energy.rel.kcal.txt")
  #print(f'{enerFile}')
  if os.path.isfile(enerFile):
    nRCE = nRCE + 1
    nFinished = nFinished + 1
  else:
    ## MD sim not finished, rerun
    rerunSims = rerunSims + 1
    rerunSimsFile = f'04_rerun_sims-{rerunSims}.sh'
    batchCommands.append(f'/bin/bash {rerunSimsFile}\n')
    command = f'rm {os.path.join(simdir, "output")} && rm {os.path.join(simdir, "slurm-")} && cd {simdir} && python run_sim_chain.py && cd {cwd}'
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
print(f'number of calculated densities: {nRCE}')
if nRCE == nSims:
  print(f'  -> run 04_get_results.py')
if rerunSims > 0:
  print(f'number of not finished sims: {rerunSims}')
  print(f'  -> run {batchRerunSimsFile}.')
