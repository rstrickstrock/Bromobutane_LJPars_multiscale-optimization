import os
from natsort import natsorted # natsorted(glob.glob(os.path.join(simDir, "slurm-*")))
import glob
import shutil


def removeDirs(dirs):
  for thisDir in dirs:
    try:
      shutil.rmtree(thisDir)
    except:
      print(f'Could not remove \"{thisDir}\". Continue.')


def remove_temp_opt_files(listOfDirs):
  """
  removes all intermediate files created during optimization
  only keeps the last gradient iteration
  """
  for thisDir in listOfDirs:
    #print(f'{thisDir}')
    lastDir = natsorted(glob.glob(os.path.join(thisDir, "QMMM", "g.*.0")))[-1]
    #print(f'{lastDir}')
    maxIteration = os.path.basename(lastDir).split(".")[1]
    try:
      maxIteration = int(maxIteration)
    except:
      print(f'Could not cast maxIteration = {maxIteration} to int. Exit.')
    #print(f'{maxIteration}')
    
    ### remove dirs:
    ## a.* dirs for QMMM
    removeTheseDirs = glob.glob(os.path.join(thisDir, "QMMM", "a.*"))
    #print(f'{removeTheseDirs}')
    removeDirs(removeTheseDirs)
    ## a.* dirs for PhysProp
    removeTheseDirs = glob.glob(os.path.join(thisDir, "PhysProp", "a.*"))
    #print(f'{removeTheseDirs}')
    removeDirs(removeTheseDirs)
    ## g.* dirs
    if maxIteration == 1:
      print(f'maxIteration == 1, will keep every g.* directory.')
    else:
      removeTheseDirs = []
      for i in range(1, maxIteration):
        tmpRemoveDirs = glob.glob(os.path.join(thisDir, "QMMM", f'g.{i}.*'))
        #print(f'{tmpRemoveDirs}')
        for tmpRemoveDir in tmpRemoveDirs:
          removeTheseDirs.append(tmpRemoveDir)
        tmpRemoveDirs = glob.glob(os.path.join(thisDir, "PhysProp", f'g.{i}.*'))
        for tmpRemoveDir in tmpRemoveDirs:
          removeTheseDirs.append(tmpRemoveDir)
      #print(f'{removeTheseDirs}')
      removeDirs(removeTheseDirs)


if __name__ == "__main__":
  pwd = os.path.dirname(os.getcwd())
  
  optDirsMaxR2 = glob.glob(os.path.join(pwd, "maxR2", "maxR2-*"))
  optDirsMinMAPE = glob.glob(os.path.join(pwd, "minMAPE", "minMAPE-*"))
  #print(f'{optDirsMaxR2}')
  #print(f'{optDirsMinMAPE}')
  
  remove_temp_opt_files(optDirsMaxR2)
  remove_temp_opt_files(optDirsMinMAPE)












