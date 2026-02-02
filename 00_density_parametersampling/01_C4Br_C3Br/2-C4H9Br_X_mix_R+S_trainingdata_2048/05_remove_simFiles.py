import os
import glob
import shutil

densdirs = glob.glob(os.path.join(os.getcwd(), "experiments", "*"))
for densdir in densdirs:
  #print(f'{densdir}')
  eminFiles = glob.glob(os.path.join(densdir, "density", "11_emin.*"))
  #print(f'{eminFiles}')        
  for eminFile in eminFiles:
    if os.path.isfile(eminFile):
      os.remove(eminFile)
  nvtFiles = glob.glob(os.path.join(densdir, "density", "12_nvt.*"))
  #print(f'{nvtFiles}')        
  for nvtFile in nvtFiles:
    if os.path.isfile(nvtFile):
      os.remove(nvtFile)
  nptFiles = glob.glob(os.path.join(densdir, "density", "13_npt.*"))
  #print(f'{eminFiles}')        
  for nptFile in nptFiles:
    if os.path.isfile(nptFile):
      os.remove(nptFile)
  prodFiles = glob.glob(os.path.join(densdir, "density", "14_prod.*"))
  #print(f'{prodFiles}')        
  for prodFile in prodFiles:
    if prodFile.endswith(".edr"):
      pass
    elif os.path.isfile(prodFile):
      os.remove(prodFile)
  slurmFiles = glob.glob(os.path.join(densdir, "density", "slurm-*"))
  #print(f'{slurmFiles}')        
  for slurmFile in slurmFiles:
    if os.path.isfile(slurmFile):
      os.remove(slurmFile)
  
  tmpFiles = glob.glob(os.path.join(densdir, "density", "#*"))
  #print(f'{tmpFiles}')
  for tmpFile in tmpFiles:
    #print(f'tmpFile: {tmpFile}')
    if os.path.isfile(tmpFile):
      os.remove(tmpFile)
      
  ffDir = os.path.join(densdir, "density", "oplsaa_robin.ff")
  if os.path.isdir(ffDir):
    shutil.rmtree(ffDir)
    
  

