import json
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt


# Break the y-axis into two parts
SPLITS = list(reversed([
    (0, 15+2),
    (20-4, 60),
]))


COLORS = {
    "gray"        : [118/255, 118/255, 118/255],
    "red"         : [219/255,  58/255,  53/255],
    "blue"        : [ 23/255, 114/255, 240/255],
    "green"       : [  8/255, 139/255,   5/255],
    "lightgray"   : [171/255, 171/255, 171/255],
    "lightred"    : [255/255, 138/255, 118/255],
    "lightblue"   : [145/255, 167/255, 255/255],
    "lightgreen"  : [118/255, 188/255,  99/255],
}


def plot_and_fit(fits, ax_splits, k, runtime, label, marker, color, fit_color):
     # BBHash is only available for k <= 12
    fit = linregress(k[:len(runtime)], runtime)
    fit_dict = {
        "slope": fit.slope,
        "slope_err": fit.stderr,
        "intercept": fit.intercept,
        "intercept_err": fit.intercept_stderr,
        "coeff": fit.rvalue**2
    }
    fits[label] = fit_dict

    for ax in ax_splits:
        ax.plot(
            k[:len(runtime)], runtime,
            linestyle="None", marker=marker, color=color, label=label
        )
        limits = np.array([0, k[-1]+1])
        ax.plot(
            limits, fit.slope * limits + fit.intercept,
            linestyle="--", marker="None", color=fit_color
        )


with open("benchmark.json") as file:
    runtimes = json.load(file)
k = np.array(runtimes["k"])
runtime_naive    = np.array(runtimes["naive"]   ) * 1e-3
runtime_fast     = np.array(runtimes["fast"]    ) * 1e-3
runtime_bbhash   = np.array(runtimes["bbhash"]  ) * 1e-3


fig, ax_splits = plt.subplots(
    len(SPLITS), 1, sharex=True, figsize=(6.0, 4.0)
)
fig.subplots_adjust(hspace=0.05)

fits = {}
plot_and_fit(fits, ax_splits, k, runtime_naive,  "naive",  "o", COLORS["gray"], COLORS["lightgray"])
plot_and_fit(fits, ax_splits, k, runtime_fast,   "fast",   "d", COLORS["red"],  COLORS["lightred"])
plot_and_fit(fits, ax_splits, k, runtime_bbhash, "bbhash", "s", COLORS["blue"], COLORS["lightblue"])

for ax, (begin, end) in zip(ax_splits, SPLITS):
    ax.set_ylim(begin, end)
upper_ax, lower_ax = ax_splits
lower_ax.set_xlim(0, k[-1]+1)
lower_ax.set_xlabel("k")
# Place a common y-label between the two axes
lower_ax.set_ylabel("run time (µs)")
lower_ax.yaxis.set_label_coords(-0.1, 1.0)
# Remove axis between the two splits
upper_ax.spines.bottom.set_visible(False)
upper_ax.tick_params(bottom=False)
lower_ax.spines.top.set_visible(False)
lower_ax.xaxis.tick_bottom()
# Top split gets the legend
upper_ax.legend(loc="upper right")

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
upper_ax.plot([0, 1], [0, 0], transform=upper_ax.transAxes, **kwargs)
lower_ax.plot([0, 1], [1, 1], transform=lower_ax.transAxes, **kwargs)

fig.savefig("benchmark.svg")


with open("fits.json", "w") as file:
    json.dump(fits, file, indent=4)