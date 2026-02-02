import os
import glob
import shutil
import subprocess
from natsort import natsorted

#dirs = ["30_1+2-Bromobutane_opt", "31_1-Bromobutane_opt", "31_2-Bromobutane_opt", "32_1+2-Bromobutane_opt_PR", "33_1-Bromobutane_opt_PR", "33_2-Bromobutane_opt_PR"]
dirs = ["31_2-Bromobutane_opt", "33_2-Bromobutane_opt_PR"]

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
      evalDir = os.path.join(thisDir, "evalDensities_withOptParams")
      simDirs = glob.glob(os.path.join(evalDir, "*"))
      #print(f'{simDirs}')
      thisCWD = os.getcwd()
      for simDir in simDirs:
        #print(f'{simDir}')
        os.chdir(simDir)
        thisCommand = f'sbatch batch_slurm_eval_PhysProp.sh'
        p = subprocess.Popen(thisCommand, stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        p_status = p.wait()
        os.chdir(thisCWD)

  
  
  
  
  
    




