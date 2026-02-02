import glob
import os
import pandas as pd

thisEnergyFiles = glob.glob(os.path.join("output", "*totalEnergy*"))
#print(f'{thisEnergyFiles}')

outFile = os.path.join("output", "Energy.rel.kcal.txt")
if os.path.exists(outFile):
  os.remove(outFile)
  
dfData = pd.DataFrame({"molecule": ["2-bromobutane_conformer_01", "2-bromobutane_conformer_02", "2-bromobutane_conformer_03"],
                       "absEnergy [kJ/mol]": [None] * 3,
                       "relEnergy [kJ/mol]": [None] * 3,
                       "relEnergy [kcal/mol]": [None] * 3})
#print(f'{dfData}')

for energyFile in thisEnergyFiles:
  #print(f'{energyFile}')
  n = os.path.basename(energyFile).split("_")[0]
  if n.startswith('1'):
    n = n[1]
    #print(f'{n}')
    molecName = f'1-bromobutane_conformer_0{n}'
  elif n.startswith('2'):
    n = n[1]
    #print(f'{n}')
    molecName = f'2-bromobutane_conformer_0{n}'
  
  with open(energyFile, 'r') as f:
    lines = f.readlines()
  
  last_line = lines[-1]
  #print(f'File: {energyFile}\nlast_line: {last_line}\n')
  for item in last_line.split(" "):
    if len(item) > 2 and item != "Total" and item != "Energy":
      dfData.loc[dfData["molecule"] == molecName, "absEnergy [kJ/mol]"] = float(item)
      break
   
#print(f'{dfData}')

#conf4AbsEner_1bromo = dfData.loc[dfData["molecule"] == "1-bromobutane_conformer_04", "absEnergy [kJ/mol]"].to_numpy()[0]
conf1AbsEner_2bromo = dfData.loc[dfData["molecule"] == "2-bromobutane_conformer_01", "absEnergy [kJ/mol]"].to_numpy()[0]
f = open(outFile, 'w')
for n in range(1,len(dfData["molecule"])+1):
  molecName = dfData.iloc[n-1].molecule
  #print(f'{molecName}')
  confNAbsEner = dfData.loc[dfData["molecule"] == molecName, "absEnergy [kJ/mol]"].to_numpy()[0]
  if molecName.startswith('1-bromobutane'):
    dfData.loc[dfData["molecule"] == molecName, "relEnergy [kJ/mol]"] = confNAbsEner - conf4AbsEner_1bromo
    dfData.loc[dfData["molecule"] == molecName, "relEnergy [kcal/mol]"] = (confNAbsEner - conf4AbsEner_1bromo)/4.184
    f.write(f'{molecName} {(confNAbsEner - conf4AbsEner_1bromo)/4.184}\n')
  elif molecName.startswith('2-bromobutane'):
    dfData.loc[dfData["molecule"] == molecName, "relEnergy [kJ/mol]"] = confNAbsEner - conf1AbsEner_2bromo
    dfData.loc[dfData["molecule"] == molecName, "relEnergy [kcal/mol]"] = (confNAbsEner - conf1AbsEner_2bromo)/4.184
    f.write(f'{molecName} {(confNAbsEner - conf1AbsEner_2bromo)/4.184}\n')
f.close()
print(f'{dfData}')

 
