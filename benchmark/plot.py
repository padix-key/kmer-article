import json
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt


COLORS = {
    # L = 50
    "gray"        : [118/255, 118/255, 118/255],
    "red"         : [219/255,  58/255,  53/255],
    "blue"        : [ 23/255, 114/255, 240/255],
    "green"       : [  8/255, 139/255,   5/255],
    "violet"      : [203/255,  56/255, 170/255],
    # L = 70
    "lightgray"   : [171/255, 171/255, 171/255],
    "lightred"    : [255/255, 138/255, 118/255],
    "lightblue"   : [145/255, 167/255, 255/255],
    "lightgreen"  : [118/255, 188/255,  99/255],
    "lightviolet" : [241/255, 135/255, 211/255],
    "yellow"      : [200/255, 168/255,  67/255],
}


def plot_and_fit(fits, ax, k, runtime, label, marker, color, fit_color):
    fit = linregress(k, runtime)
    fit_dict = {
        "slope": fit.slope,
        "slope_err": fit.stderr,
        "intercept": fit.intercept,
        "intercept_err": fit.intercept_stderr,
        "coeff": fit.rvalue**2
    }
    fits[label] = fit_dict

    ax.plot(k, runtime, linestyle="None", marker=marker, color=color, label=label)
    ax.plot(k, fit.slope * k + fit.intercept, linestyle="--", marker="None", color=fit_color)


with open("benchmark.json") as file:
    runtimes = json.load(file)
k = np.array(runtimes["k"])
runtime_naive    = np.array(runtimes["naive"]   ) * 1e-3
runtime_fast     = np.array(runtimes["fast"]    ) * 1e-3
runtime_bitshift = np.array(runtimes["bitshift"]) * 1e-3


fig, ax = plt.subplots(figsize=(6.0, 4.0), constrained_layout=True)

fits = {}
plot_and_fit(fits, ax, k, runtime_naive,    "naive",    "o", COLORS["gray"], COLORS["lightgray"])
plot_and_fit(fits, ax, k, runtime_fast,     "fast",     "v", COLORS["blue"], COLORS["lightblue"])
plot_and_fit(fits, ax, k, runtime_bitshift, "bitshift", "^", COLORS["red"],  COLORS["lightred"])

ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.set_xlabel("k")
ax.set_ylabel("runtime (µs)")
ax.legend()

fig.savefig("benchmark.svg")


with open("fits.json", "w") as file:
    json.dump(fits, file, indent=4)