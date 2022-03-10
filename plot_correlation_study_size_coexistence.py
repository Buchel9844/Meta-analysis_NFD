import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("../Results/NFD_meta.csv")

specs = set(data[["species.i", "species.j"]].values.flatten())


for spec in specs:
    ind = (data["species.i"] == spec) | (data["species.j"] == spec)
    studies = set(data.loc[ind, "First.Author"])
    if len(studies)>1:
        print(spec, studies)
        
studies = sorted(set(data["First.Author"]))

comp_outcome = np.empty((len(studies), 3))

for s, study in enumerate(studies):
    for i in range(3):
        comp_outcome[s, i] = np.sum((data["First.Author"] == study) 
                                    & np.isfinite(data["r_i"])
                                    & (data["comp_outcome"] == i))
    
size = np.sum(comp_outcome, axis = 1)
comp_outcome /= size[:,np.newaxis]

fig = plt.figure()

label = ["Priority effect", "Competitive exclusion", "Coexistence"]
for i in range(3):
    plt.plot(size, comp_outcome[:,i], '.', label = label[i])
    
plt.xlabel("Number of communities in study")
plt.ylabel("Proportion of outcomes")
plt.legend()

fig.savefig("Figure_ap_correlation_study_size_coexistence.pdf")