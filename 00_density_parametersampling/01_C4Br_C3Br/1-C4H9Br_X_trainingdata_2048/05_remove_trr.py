import os
import glob

densdirs = glob.glob(os.path.join(os.getcwd(), "experiments", "*"))
for densdir in densdirs:
  #print(f'{densdir}')
  trrFiles = glob.glob(os.path.join(densdir, "density", "*.trr"))
  #print(f'{trrFiles}')        
  for trrFile in trrFiles:
    if os.path.isfile(trrFile):
      os.remove(trrFile)
  
  tmpFiles = glob.glob(os.path.join(densdir, "density", "#*"))
  #print(f'{tmpFiles}')
  for tmpFile in tmpFiles:
    #print(f'tmpFile: {tmpFile}')
    if os.path.isfile(tmpFile):
      os.remove(tmpFile)

