import numpy as np
import matplotlib.pyplot as plt

from clustering_number_of_clusters import NFD_estimates, clf, probs

fig = plt.figure(figsize = (9,9))
coex_col = ["grey", "blue", "red"]
labels = [["Competitive\nexclusion", "Coexistence", "Priority Effect"],
            ["Low ND", "High ND", "variable"]]

ax_ND = fig.add_axes([0.35, 0.08, 0.55, 0.2])
ax_FD = fig.add_axes([0.08, 0.35, 0.2, 0.55])
ax_coex = fig.add_axes([0.35, 0.35,0.55,0.55])

ax_coex.scatter(NFD_estimates[...,0], NFD_estimates[...,1], 
                c = "k", alpha = 1/len(NFD_estimates), s = 5)

ND_lim = [-1,2]
ND_bins, d_ND = np.linspace(*ND_lim, 31, retstep = True)
FD_lim = [0,1]
FD_bins, d_FD = np.linspace(*FD_lim, 21, retstep = True)


##############################################################################
# compute clustering and plot
clust_col = ["orange", "green", "purple"]

scal_ND = np.sum(probs, axis = (0,1))/np.sqrt(2*np.pi*clf.covariances_[:,0,0])
scal_FD = np.sum(probs, axis = (0,1))/np.sqrt(2*np.pi*clf.covariances_[:,1,1])
#clf_n_tot = 1/np.sqrt(2*np.pi*clf.covariances_[:,0,0])
alpha = 0.7
x_ND = np.linspace(*ND_lim, 1001)
x_FD = np.linspace(*FD_lim, 1001)
for i in range(3):
    ax_ND.hist(NFD_estimates[:,np.isfinite(NFD_estimates[0,:,0]),0].flatten(),
                bins = ND_bins, alpha = alpha, color = clust_col[i],
                weights = probs[...,i].flatten(), label = labels[1][i])
    ax_FD.hist(NFD_estimates[:,np.isfinite(NFD_estimates[0,:,0]),1].flatten(),
                bins = FD_bins, alpha = alpha, color = clust_col[i],
                weights = probs[...,i].flatten(), label = labels[1][i],
                orientation = "horizontal")
    
    # multiply with bin width to have same integral on display
    ax_ND.plot(x_ND, d_ND*scal_ND[i]*np.exp(-(x_ND-clf.means_[i,0])**2/clf.covariances_[i,0,0]/2),
                color = clust_col[i], linewidth = 3)
    ax_FD.plot(d_FD*scal_FD[i]*np.exp(-(x_FD-clf.means_[i,1])**2/clf.covariances_[i,1,1]/2),x_FD, 
                color = clust_col[i], linewidth = 3)
 
    
ax_FD.legend()
# add cluster ellipses in NFD graph
theta = np.linspace(0, 2*np.pi, 1000);
ls = ["-", "--", ":"]
for i in range(3):
    eigv, eig_vec = np.linalg.eigh(clf.covariances_[i])
    ellipsis = (np.sqrt(eigv[None,:]) * eig_vec).dot([np.sin(theta), np.cos(theta)])
    for j in range(3):
        plt.plot(clf.means_[i,0] + (j+1)*ellipsis[0,:],
                 clf.means_[i,1] + (j+1)*ellipsis[1,:],
             color = clust_col[i], linestyle = ls[j], lw = 3, zorder = 0)


###############################################################################
# add layout

###############################################################################
# add layout
fs = 16
ax_coex.set_xlim(ND_lim)
ax_coex.set_ylim(FD_lim)

ax_FD.set_ylabel("Fitness differences", fontsize = fs)
ax_FD.set_xlabel("Frequency", fontsize = fs)
#ax_FD2.set_ylabel("Frequency")
#ax_FD.set_ylabel("Frequency")
ax_ND.set_ylabel("Frequency", fontsize = fs)
ax_ND.set_ylabel("Frequency", fontsize = fs)
ax_ND.set_xlabel("Niche differences", fontsize = fs)

for i, ax in enumerate([ax_FD, ax_coex, ax_ND]):
    ax.set_title("ABCDEF"[i], loc = "left")
    

ax_ND.set_xticks(np.arange(-0.5, 2, 0.5))
ax_ND.set_xticklabels(["Positive\nfrequency dependence", 0, "Negative\nfrequency dependence", 1, "Facilitation"])
ax_ND.set_xlim([-1,2])

ax_FD.set_ylim([0,1])
ax_FD.set_yticks([0,0.5, 1])

ax_coex.set_xticks([0,1])
ax_coex.set_yticks([0,0.5, 1])


ax_ND.legend(fontsize = 16, bbox_to_anchor = [-0.15, 1])
#ax_coex.set_xlabel("Niche differences")

fig.savefig("Figure_clusters.pdf")
