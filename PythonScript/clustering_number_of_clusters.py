from sklearn import mixture
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#np.random.seed(1234)

# load empirical data
data = pd.read_csv("../Results/NFD_meta.csv")
data.loc[data.Kingdom == "phytoplankton", "Kingdom"] = "Phytoplankton"
data.loc[data["lab.OR.field"] == "field", "lab.OR.field"] = "Field"
apm2 = "Annual plant\nmodel\nwith exponential\nspecies interactions"
apm2 = "Exponential Annual plant"
data.loc[data["Model.use"] == "AP2", "Model.use"] = apm2
data.loc[data["Model.use"] == "APM", "Model.use"] = "Annual plant\nmodel"

no_model = (data["Model.use"] == "None").values


# add intrinsic growth rates and invasion growth rates for cases where no model is used
no_S_i = np.isnan(data.S_i)
data.loc[no_S_i, "S_i"] = ((data.g_ii - data.g_ij)/data.g_ii)[no_S_i]
data.loc[no_S_i, "S_j"] = ((data.g_jj - data.g_ji)/data.g_jj)[no_S_i]

data.loc[no_model, "fi_0"] = 1
data.loc[no_model, "r_i"] = 1 - data.loc[no_model, "S_i"]
data.loc[no_model, "r_j"] = 1 - data.loc[no_model, "S_j"]
data.loc[no_model, "fj_0"] = 1
data["eta_i"] = data["fi_0"]*data["FS_i"]/(data["FS_i"] - 1)

def compute_NFD_for_clusters(data_org = data.copy(), include_nan = False):
    # create estimates for niche and fitness differences where these are not known
    
    data = data_org.copy()
    # used to estimate eta_i for cases where eta_i is not known
    p = np.random.uniform(0,1, len(data))
    
    # given niche overlap compute no-niche growth rate
    data.loc[no_model, "eta_i"] = p*(data.fi_0 - data.fj_0*np.abs((data.fi_0-data.r_i)/(data.fj_0-data.r_j)))
    NO = np.abs(data.fi_0-data.r_i)/(data.fi_0 - data.eta_i)
    data.loc[no_model, "eta_j"] = (data.fj_0 - np.abs(data.fj_0-data.r_j)/NO)[no_model]
    
    # given eta, compute niche and fitness differences
    data.loc[no_model, "NS_i"] = ((data.r_i-data.eta_i)/(data.fi_0-data.eta_i))[no_model]
    data.loc[no_model, "FS_i"] = -(data.eta_i/(data.fi_0 - data.eta_i))[no_model]
    data.loc[no_model, "NS_j"] = ((data.r_j-data.eta_j)/(data.fj_0-data.eta_j))[no_model]
    data.loc[no_model, "FS_j"] = -(data.eta_j/(data.fj_0 - data.eta_j))[no_model]
    
    # check for correctness, NS_1 should be limited by 1-(fi_0 -r_i)/ri_0
    test1 = np.abs((data.fi_0 - data.r_i)/data.fi_0)
    test2 = np.abs((data.fj_0 - data.r_j)/data.fj_0)
    test1, test2 = np.sort([test1, test2], axis = 0)
    if np.any(np.isfinite(data["NS_i"]) &
              ((np.abs(1-data["NS_i"]) < test1 - 1e-8)
               | (np.abs(1-data["NS_i"]) > test2 + 1e-8))):
        raise

    if include_nan:
        ind = np.argmax(data[["FS_i", "FS_j"]].values, axis = 1)
        data["FS"] = data[["FS_i", "FS_j"]].values[np.arange(len(data)),ind]
        data["NS"] = data[["NS_i", "NS_j"]].values[np.arange(len(data)),ind]
        
        return data[["NS", "FS"]].values
    
    # return only data necessary for clusterin  
    ND = data[["NS_i", "NS_j"]].values
    FD = data[["FS_i", "FS_j"]].values
    ind = np.argmax(FD, axis = 1)
    NFD = np.array([ND[np.arange(len(ND)), ind],
                    FD[np.arange(len(ND)), ind]]).T
    
    return NFD[np.isfinite(NFD).all(axis = 1)]

# compute clustering and plot
np.random.seed(0)
NFD = compute_NFD_for_clusters()
n_clus = 3
clf = mixture.GaussianMixture(n_components = n_clus, covariance_type="full")
clf.fit(NFD)

# 20 estimates for each community where NFD is not known exactly
NFD_estimates = np.array([compute_NFD_for_clusters(include_nan = True)
                          for i in range(10)])

ind_no_nan = np.isfinite(NFD_estimates[0,:,1])
# apply clustering
probs = np.array([clf.predict_proba(NFD[ind_no_nan]) for NFD in NFD_estimates])

# determine competitive outcome
data["comp_outcome"] = np.sum(data[["r_i", "r_j"]]>0, axis = 1)
data.loc[np.isnan(data["r_i"]), "comp_outcome"] = np.nan

# compute clustering based on estimates of niche and fitness differences
for i, clust in enumerate(["low_ND_clust", "high_ND_clust", "var_clust"]):
    data[clust] = np.nan
    data.loc[ind_no_nan, clust] = np.nanmean(probs[...,i], axis = 0)
    
# frequency dependency estimates
bounds = [[-np.inf, 0], [0, 1], [1, np.inf]]
for i, freq in enumerate(["prob_prio", "prob_neg_freq", "prob_fac"]):
    data[freq] = np.mean((NFD_estimates[...,0] >= bounds[i][0])
                         & (NFD_estimates[...,0] < bounds[i][1]), axis = 0)
    data.loc[~ind_no_nan, freq] = np.nan
    
# save data
data.to_csv("NFD_clustering.csv", index = False)

if __name__ == "__main__":
    fig, ax = plt.subplots(3,3, figsize = (9,6), sharex = True, sharey = True)
    ax = ax.flatten()
    n_clus = 10
    n_rep = 2
    aic = np.empty((n_rep, n_clus))
    
    for i in range(n_clus):
        for j in range(n_rep):
            NFD = compute_NFD_for_clusters()
            NFD = NFD[:,[0]]
            clf = mixture.GaussianMixture(n_components=i+1, covariance_type="full")
            clf.fit(NFD)
            aic[j,i] = clf.bic(NFD)
            NFD = compute_NFD_for_clusters(data)
        if i<len(ax):
            ax[i].scatter(*NFD.T, c = clf.fit_predict(NFD), s = 2)
            ax[i].set_xlim([-1.2,1.5])
            ax[i].set_ylim([0, 1])
            
            cov = clf.covariances_
            ax[i].plot(*clf.means_.T, 'o', color = "r")
            ax[i].set_title("{} clusters".format(i+1))            
    
    fig.tight_layout()
    fig.savefig("Figure_ap_multiple_clusters.png")
    
    fig = plt.figure()
    plt.plot(1 + np.arange(n_clus), np.mean(aic, axis = 0), 'o')
    plt.xlabel("Number of clusters")
    plt.ylabel("Goodnes of fit")
    fig.savefig("Figure_ap_number_of_clusters.png")
