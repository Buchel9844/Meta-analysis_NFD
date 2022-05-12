import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.cm import Set1 as cmap

data = pd.read_csv("../Results/NFD_metaanalysis.csv")
data = data[data.Definition == "Spaak"]
#data.loc[data.Kingdom == "phytoplankton", "Kingdom"] = "Phytoplankton"
data.loc[data["lab.OR.field"] == "field", "lab.OR.field"] = "Field"

keys = ["Kingdom", "Model.use", "Initial.Definition", "lab.OR.field"]
sub_keys = {key:list(set(data[key])) for key in keys}
sub_keys["Kingdom"] = ["Phytoplankton", "Bacteria", "Annual plant", "Perennial plant"]
sub_keys["Model.use"] = ["Lotka Volterra", "APM", "AP2", "None"]
sub_keys["lab.OR.field"] = ["Field", "Mesocosm", "Greenhouse", "lab"]
sub_keys["Initial.Definition"] =  sorted(sub_keys["Initial.Definition"])

ratios = pd.DataFrame(data = np.nan, columns = keys,
            index = np.arange(max([len(sub_keys[sk]) for sk in sub_keys])))

for key in keys:
    for i, sk in enumerate(sub_keys[key]):
        ratios.loc[i, key] = np.sum(data[key] == sk)/len(data)
        
fig, ax = plt.subplots(2,2, sharex = True, sharey = True, figsize =(7,7))
ax[0,0].set_ylabel("Fitness differences")
ax[1,0].set_ylabel("Fitness differences")
ax[1,0].set_xlabel("Niche differences")
ax[1,1].set_xlabel("Niche differences")
ax = ax.flatten()

w = 0.5
colors = cmap(np.linspace(0,1,len(ratios)))
names = ["System", "Community\nmodel", "Initial\nDefinition", "Setting"]

txt_ind = np.append(np.zeros((1,len(keys))), ratios, axis = 0)
txt_ind = np.cumsum(txt_ind, axis = 0)
txt_ind = (txt_ind[:-1] + txt_ind[1:])/2
ha = ["right", "left"]
err = w/10
x_loc = [0-err,w+err]
x_loc_marker = [0,w]
x = 2.5
ax[0].set_xlim(np.nanpercentile(data["ND"], [x,100-x]))

for i, key in enumerate(keys):
    print(key)
    for j,subkey in enumerate(sub_keys[key]):
        ind = data[key] == subkey
        color = np.repeat(colors[[j]], sum(ind), axis = 0)
        color[:,-1] = 0.5
        color[:, -1][data.loc[ind, "Model.use"] == "None"] /= 4
        ax[i].scatter(data.loc[ind, "ND"], data.loc[ind, "FD"],
                      s = 5, color = color, label = subkey)
    ax[i].legend(ncol = 2 if len(sub_keys[key])>4 else 1)
    ax[i].set_yticks([0,1])
    ax[i].axvline(0, color = "grey")
    ax[i].axvline(1, color = "grey")
    ax[i].axhline(0, color = "grey")
    ax[i].plot(ax[i].get_xlim(), ax[i].get_xlim(), 'k--')
    ax[i].set_title(names[i])



ax[0].set_ylim([-0.75,1])
fig.savefig("Figure_ap_NFD_empirical_setting.pdf")
"""

dataframes = {key: pd.DataFrame(None, index = sub_keys["Kingdom"],
                                columns = sub_keys[key])
              for key in keys}

for key in keys:
    for kingdom in sub_keys["Kingdom"]:
        dataframes[key].loc[kingdom, "Total"] =  np.sum(
                data["Kingdom"] == kingdom)
        for subkey in sub_keys[key]:
            dataframes[key].loc[kingdom, subkey] = np.sum(
                (data["Kingdom"] == kingdom) & (data[key] == subkey))
    
    # reduce phytoplankton counts due to uncertainty of NFD
    dataframes[key].loc["Phytoplankton",:] /= 10
    
for key in keys:
    if key == "Kingdom":
        continue
    dataframes[key].to_csv("Table_{}.csv".format(key))
    
from scipy.stats import ttest_ind
king = "Perennial plant"
a = data.ND[(data.Kingdom == king)
            & (data["Model.use"] == "Lotka Volterra")]
b = data.ND[(data.Kingdom == king)
            & (data["Model.use"] == "APM")]
plt.figure()
print(ttest_ind(a,b, equal_var = False, nan_policy="omit"))
bins = np.linspace(*np.nanpercentile(np.append(a,b), [1,99]), 30)
plt.hist(a, density = True, alpha = .5, bins = bins, label = "LV")
plt.hist(b, density = True, alpha = .5, bins = bins, label = "APM")
plt.legend()"""