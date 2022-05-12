import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kruskal

from clustering_number_of_clusters import data, NFD_estimates

fig = plt.figure(figsize = (9,9))
coex_col = ["grey", "blue", "red"]
labels = [["Competitive\nexclusion", "Coexistence", "Priority Effect"],
            ["Low ND", "High ND", "variable"]]

ax_ND = fig.add_axes([0.35, 0.08, 0.55, 0.2])
ax_FD = fig.add_axes([0.08, 0.35, 0.2, 0.55])
ax_coex = fig.add_axes([0.35, 0.35,0.55,0.55])
#ax_coex = fig.add_subplot(1,2,1, sharex = ax_ND)

for i, outcome in enumerate([1, 2, 0]):
    ind = (data["comp_outcome"] == outcome).values
    ax_coex.scatter(NFD_estimates[:, ind, 0], NFD_estimates[:, ind, 1], 
                c = coex_col[i], alpha = 1/len(NFD_estimates), s = 5)

ND_lim = [-1,2]
ND_bins, d_ND = np.linspace(*ND_lim, 31, retstep = True)
FD_lim = [0,1]
FD_bins, d_FD = np.linspace(*FD_lim, 21, retstep = True)

# histograms of niche differences
inds = [np.isfinite(data.comp_outcome.values), # all
        data.comp_outcome.values != 1, # coexistence + priority
        data.comp_outcome.values == 0] # priority
weights = np.full(NFD_estimates.shape[:2], 1/len(NFD_estimates))
for i in range(3):
    ax_ND.hist(NFD_estimates[:,inds[i], 0].flatten(),
                bins = ND_bins, color = coex_col[i], label = labels[0][i])
    ax_FD.hist(NFD_estimates[:,inds[i], 1].flatten(),
                bins = FD_bins, color = coex_col[i], label = labels[0][i], 
                orientation ="horizontal")


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

ax_coex.axvline(0, color = "k")
ax_coex.axvline(1, color = "k")
ax_coex.plot([0,1], [0,1], color = "b", linewidth = 2, linestyle = "--")
ND_prio = np.linspace(-1, 0, 1000)
ax_coex.plot(ND_prio, 1-1/(1-ND_prio), color = "r", linestyle = "--")


ax_ND.legend(fontsize = 16, bbox_to_anchor = [-0.15, 1])
#ax_coex.set_xlabel("Niche differences")
#ax_coex.set_ylabel("Fitness differences")
"""
ax_lab1 = fig.add_subplot(4,4,10)
ax_lab2 = fig.add_subplot(4,4,13)

cols = [coex_col, clust_col]
for j, ax in enumerate([ax_lab1, ax_lab2]):
    for i in range(3):
        ax.hist(np.ones(5), color = cols[j][i], label = labels[j][i])
    
    ax.set_frame_on(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim([-1,0])
    ax.legend()"""

#fig.tight_layout()

fig.savefig("Figure_comp_outcome.pdf")

##############################################################################
# compute kruskal wallis difference between distributions

no_nan = np.isfinite(NFD_estimates).all(axis = (0,2))
comp = data.loc[no_nan, "comp_outcome"]
NFD_estimates = NFD_estimates[:,no_nan]
for i in range(3):
    j = (i+1)%3
    
    print("ND", "comp_outcome = {} and {}".format(i,j),
          kruskal(NFD_estimates[0, comp == i, 0].flatten(),
                  NFD_estimates[0, comp == j, 0].flatten())[1])
    print("FD", "comp_outcome = {} and {}".format(i,j),
          kruskal(NFD_estimates[1, comp == i, 1].flatten(),
                  NFD_estimates[1, comp == j, 1].flatten())[1])
    
print("FD_medians", [np.round(np.nanmedian(NFD_estimates[:,comp == i, 1]), 3) for i in range(3)])
print("ND_medians", [np.round(np.nanmedian(NFD_estimates[:,comp == i, 0]), 3) for i in range(3)])
