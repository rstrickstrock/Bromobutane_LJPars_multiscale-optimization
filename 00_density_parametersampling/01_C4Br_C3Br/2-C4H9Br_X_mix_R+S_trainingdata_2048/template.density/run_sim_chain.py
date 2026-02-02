import subprocess
import time
#import os
#from glob import glob
#import threading

time_check_intervall = 60 #seconds
commands =[['sbatch', '01_batch_slurm_emin.sh'],
           ['sbatch', '02_batch_slurm_nvt.sh'],
           ['sbatch', '03_batch_slurm_npt.sh'],
           ['sbatch', '04_batch_slurm_prod.sh'],
           ['sbatch', '05_batch_slurm_eval_density.sh']]

for command in commands:
  ## start simulation     
  simulation = subprocess.Popen(command, stdout=subprocess.PIPE)
  #print(f'started {command[1]}')
  #print(f'in cwd: {os.getcwd()}')
  ## get std out to determine job ID
  simulation_stdout = simulation.stdout.read()
  batch_job_id = [int(s) for s in simulation_stdout.split() if s.isdigit()]
  batch_job_id = batch_job_id[0]

  # check if batch scrip terminated and wait if not
  while 1:
    if simulation.poll() == 0:
      break
    else:
      time.sleep(1)
    
  # wait for job to finish
  while 1:
    check_job_status = subprocess.Popen(['squeue', '-j', str(batch_job_id)], stdout=subprocess.PIPE)
    check_job_status = check_job_status.stdout.read()
    if str(check_job_status).find(str(batch_job_id)) > 0:
      ## job is still running, check again in time_check_intervall seconds
      time.sleep(time_check_intervall)
    else:
      ## job has finished
      break
