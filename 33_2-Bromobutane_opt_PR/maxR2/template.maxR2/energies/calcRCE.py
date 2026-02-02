import os
import pandas as pd

substance = "2-bromobutane"

absEnergyFiles = os.path.join("output", "Energy.abs.kcal.txt")
#print(f'{absEnergyFiles}')

relEnergyFile = os.path.join("output", "Energy.rel.kcal.txt")
if os.path.exists(relEnergyFile):
  os.remove(relEnergyFile)

if __name__ == "__main__":
  with open(absEnergyFiles, 'r') as f:
    lines = f.readlines()
    
  absRCEs = []
  for line in lines:
    absRCEs.append(float(line.split(" ")[1]))
  
  minRCE1 = None
  minRCE2 = None
  minRCEboth = None
  if substance == "1-bromobutane":
    minRCE1 = absRCEs[3]
    dfData = pd.DataFrame({"molecule": ["1-bromobutane_conformer_01", "1-bromobutane_conformer_02", "1-bromobutane_conformer_03", "1-bromobutane_conformer_04", "1-bromobutane_conformer_05", "1-bromobutane_conformer_06", "1-bromobutane_conformer_07", "1-bromobutane_conformer_08", "1-bromobutane_conformer_09"],
                           "absEnergy [kJ/mol]": [None] * 9,
                           "relEnergy [kJ/mol]": [None] * 9,
                           "relEnergy [kcal/mol]": [None] * 9})
  elif substance == "2-bromobutane":
    minRCE2 = absRCEs[0]
    dfData = pd.DataFrame({"molecule": ["2-bromobutane_conformer_01", "2-bromobutane_conformer_02", "2-bromobutane_conformer_03"],
                           "absEnergy [kJ/mol]": [None] * 3,
                           "relEnergy [kJ/mol]": [None] * 3,
                           "relEnergy [kcal/mol]": [None] * 3})
  elif substance == "both":
    minRCEboth = [absRCEs[3],  absRCEs[0]]
    dfData = pd.DataFrame({"molecule": ["1-bromobutane_conformer_01", "1-bromobutane_conformer_02", "1-bromobutane_conformer_03", "1-bromobutane_conformer_04", "1-bromobutane_conformer_05", "1-bromobutane_conformer_06", "1-bromobutane_conformer_07", "1-bromobutane_conformer_08", "1-bromobutane_conformer_09", "2-bromobutane_conformer_01", "2-bromobutane_conformer_02", "2-bromobutane_conformer_03"],
                           "absEnergy [kJ/mol]": [None] * 12,
                           "relEnergy [kJ/mol]": [None] * 12,
                           "relEnergy [kcal/mol]": [None] * 12})
  else:
    print(f'\"substance\" needs to be \"1-bromobutane\" or \"2-bromobutane\" or \"both\" (is: \"{substance}\"). Exit.')
    exit()
  #print(f'{dfData}')

  relRCEs = []
  if not minRCE1 is None:
    for i in range(len(absRCEs)):
      relRCEs.append((minRCE1 - absRCEs[i])/4.184)
  if not minRCE2 is None:
    for i in range(len(absRCEs)):
      relRCEs.append((minRCE2 - absRCEs[i])/4.184)
  if not minRCEboth is None:
    for i in range(0,9):
      relRCEs.append((minRCE1 - absRCEs[i])/4.184)
    for i in range(9,12):
      relRCEs.append((minRCE2 - absRCEs[i])/4.184)
  #print(f'{relRCEs}')

  f = open(relEnergyFile, 'w')
  for i in range(0, len(relRCEs)):
    if substance == "1-bromobutane":
      confName = f'1-bromobutane_conformer_0{i+1}'
    elif substance == "2-bromobutane":
      confName = f'2-bromobutane_conformer_0{i+1}'
    elif substance == "both":
      if i < 9: # I know :(
        confName = f'1-bromobutane_conformer_0{i+1}'
      else:
        confName = f'2-bromobutane_conformer_0{i+1}'
    # print(f'{confName}')
    f.write(f'{confName} {relRCEs[i]}\n')
  f.close()
