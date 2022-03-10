import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from matplotlib.cm import Set1 as cmap
from sklearn import mixture
from sklearn.metrics import adjusted_rand_score

from clustering_number_of_clusters import compute_NFD_for_clusters, data

data_all = pd.read_csv("../Results/NFD_metaanalysis.csv")

data_car = data_all[data_all["Definition"] == "Carroll"]
data_car = data_car[np.isfinite(data_car[["ND", "FD"]].values).all(axis = 1)]
data_zhao = data_all[data_all["Definition"] == "Zhao"]
data_zhao = data_zhao[np.isfinite(data_zhao[["ND", "FD"]].values).all(axis = 1)]

n_rep = 9
data_dict = {"Carroll": n_rep*[data_car[["ND", "FD"]].values],
             "Zhao": n_rep*[data_zhao[["ND", "FD"]].values],
             "Spaak": [compute_NFD_for_clusters(data) for i in range(n_rep)]}

defs = ["Spaak", "Carroll", "Zhao"]
clfs = {key: np.empty(n_rep, dtype = "object") for key in defs}

n_clusts = np.arange(1, 11)
aic, bic, loglike = np.empty((3, len(defs), n_rep, len(n_clusts)))
rand_metrics = np.full((len(defs), n_rep, n_rep, len(n_clusts)), np.nan)

n_pref = np.nan

for d, definition in enumerate(defs):
    for j, n_clus in enumerate(n_clusts):
        if j == n_pref:
            fig, ax = plt.subplots(3,3, sharex = True, sharey = True)
            ax = ax.flatten()
            fig.suptitle(definition)
        predicts = np.empty((n_rep, len(data_dict[definition][0])))
        for i in range(n_rep):
            clf = mixture.GaussianMixture(n_components = n_clus, covariance_type="full")
            clf.fit(data_dict[definition][i])
            # compute metrics of fit
            predicts[i] = clf.predict(data_dict[definition][i])
            bic[d, i, j] = clf.bic(data_dict[definition][i])
            aic[d, i, j] = clf.aic(data_dict[definition][i])
            loglike[d,i,j] = clf.score(data_dict[definition][i])
            if j == n_pref:
                ax[i].scatter(*data_dict[definition][i].T, s = 2,
                              c = clf.predict(data_dict[definition][i]))
                ax[i].set_xlim(np.nanpercentile(data_dict[definition][i][:,0], [1,99]))
                ax[i].set_ylim(np.nanpercentile(data_dict[definition][i][:,1], [1,99]))
            for k in range(i):
                rand_metrics[d,i,k,j] = adjusted_rand_score(predicts[k],
                                                           predicts[i])
            
            

fig, ax = plt.subplots(2,2, sharex = True)
color = ["r","b", "g"]
for d, definition in enumerate(["Spaak"]):
    ax[0,0].plot(n_clusts, np.mean(aic[d], axis = 0), '.',
                 label = definition, color = color[d])
    ax[0,1].plot(n_clusts, np.mean(bic[d], axis = 0), '.',
                 label = definition, color = color[d])
    ax[1,0].plot(n_clusts, np.mean(loglike[d], axis = 0), '.',
                 label = definition, color = color[d])
    ax[1,1].plot(n_clusts[1:]+d/10, np.nanmean(rand_metrics[d].reshape(-1, len(n_clusts)), axis = 0).T[1:],
                 '.', color = color[d], label = definition)
    
ax[0,0].set_ylabel("AIC")
ax[0,1].set_ylabel("BIC")
ax[1,0].set_ylabel("Log-likelyhood")
ax[1,1].set_ylabel("rand metric")
#ax[1,1].legend()
ax[1,1].set_xlabel("Number of clusters")
ax[1,0].set_xlabel("Number of clusters")

fig.tight_layout()


fig.savefig("Figure_ap_number_of_clusters.pdf")
    
