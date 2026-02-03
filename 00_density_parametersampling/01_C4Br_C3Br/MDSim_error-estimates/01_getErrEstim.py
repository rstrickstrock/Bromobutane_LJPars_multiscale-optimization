import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import sys

datasetnames = ["1-bromobutane", "2-bromobutane"]
statisticsFiles = ["1-Bromobutane_Trainingsdata.csv", "2-Bromobutane_Trainingsdata.csv"]
colors = ['#0ba1e2', '#e24c0b']

figDPI = 100

try:
  saveOrShow = sys.argv[1]
except:
  saveOrShow = "show"
if saveOrShow == "save":
  pass
else:
  saveOrShow = "show"


for n in range(len(statisticsFiles)):
  if not os.path.isfile(statisticsFiles[n]):
    print(f"Can not find and open '{statisticsFiles[n]}'. Exit.")
    exit()
  else:
    thisDF = pd.read_csv(statisticsFiles[n])
    try:
      thisDF = thisDF.drop(columns=["Unnamed: 0"])
    except:
      print("Column 'Unnamed: 0' not found.")

  densities = thisDF["density"].to_numpy()
  densityErrs = thisDF["density_err"].to_numpy()
  relErrs = (densityErrs / densities) * 100  # [%]

  avgRelErr = np.mean(relErrs)
  ErrErrEstim = np.std(relErrs)
  print(f"{datasetnames[n]} - avg Error: ({avgRelErr:.3f} +/- {ErrErrEstim:.3f})%")

  plt.figure(figsize=(5, 5))
  plt.hist(relErrs,
       bins=50,
       density=True,
       alpha=0.5,
       label="uncertainties",
       color=colors[0])

  plt.axvline(avgRelErr,
        color=colors[1],
        linestyle="--",
        linewidth=2,
        label="avg. uncertainty")
        
  lowerBoundary = avgRelErr - 1 * ErrErrEstim
  if lowerBoundary < 0:
    lowerBoundary = 0
  upperBoundary = avgRelErr + 1 * ErrErrEstim
  plt.axvspan(lowerBoundary,
        upperBoundary,
        color=colors[1],
        alpha=0.2,
        label="variance")

  plt.xlabel("Rel. Uncertainty [%]", fontweight='bold', fontsize=18)
  plt.ylabel("Occurence", fontweight='bold', fontsize=18)
  plt.xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
  plt.ylim([0, 9])
  plt.yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
  
  plt.tick_params(axis="x", labelsize=15)
  plt.tick_params(axis="y", labelsize=0)
  #plt.title("Relative Error Distributions")
  plt.legend()
  plt.grid(True, linestyle=':', alpha=0.6)
  plt.tight_layout()
  if saveOrShow == "show":
    plt.show()
  if saveOrShow == "save":
    plt.savefig(f'MDSim-uncertainties_{datasetnames[n]}.png', dpi=figDPI)

#1-bromobutane - avg Error: (0.241 +/- 0.453)%
#2-bromobutane - avg Error: (0.189 +/- 0.272)%

