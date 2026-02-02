import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

resultsFileName, sortedBy = "optResultsSortedDensityMAPEWithDensitySims.csv", "density" # sorted by the density MAPE
#resultsFileName , sortedBy = "optResultsSortedRCEMAPEWithDensitySims.csv", "RCE"     # sorted by the RCE MAPE

preOptDensityFileName = "10_Density_1-bromobutane.preOpt.csv"
substance = "1-bromobutane"
#substance = "2-bromobutane"

colors = ["#577590", "#9bc53d", "#ffba40", "#ff8b94", "#ffee55"]

def CSVFile2DataFrame(csvFile):
  try:
    dfFromCSVFile = pd.read_csv(csvFile)
  except:
    print(f'Can not read CSV File \"{csvFile}\". Exit')
    exit()
  else:
    dfFromCSVFile = dfFromCSVFile.drop(dfFromCSVFile.columns[0], axis=1)  
  return dfFromCSVFile


if __name__ == "__main__":
  try:
    saveOrShow = sys.argv[1]
  except:
    saveOrShow = "show"
  if saveOrShow == "save":
    pass
  else:
    saveOrShow = "show"
  
  if substance == "1-bromobutane":
    targetDensity = 1275.8
  elif substance == "2-bromobutane":
    targetDensity =1258.5
  else:
    print(f'substance needs to be \"1-bromobutane\" or \"2-bromobutane\" (is \"{substance}\"). Exit.')
    exit()
  
  pwd = os.path.dirname(os.getcwd())
  
  resultsFile = os.path.join(pwd, resultsFileName)
  dfResults = CSVFile2DataFrame(resultsFile)
  #print(f'{dfResults}')
  
