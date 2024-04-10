import os
import json
import matplotlib.pyplot as plt
import stylia
import shutil
import numpy as np
import blitzgsea as blitz
from stylia.figure.figure import stylize
import pandas as pd
from stylia.colors.colors import NamedColors
import warnings
warnings.filterwarnings("ignore")

ROOT = os.path.abspath(os.path.dirname(__file__))
PATH = os.path.abspath(os.path.join(ROOT, "../results/proteins/"))
DATA = os.path.abspath(os.path.join(ROOT, "../data"))
CACHE = os.path.abspath(os.path.join(ROOT, "../cache"))

COLORS = stylia.colors.colors.NamedColors()

top_plots_number = 99999

def iterate_all_cache_folders():
    R = []
    for profile_subfolder in os.listdir(os.path.join(CACHE, "detailed")):
        for annotations_subfolder in os.listdir(os.path.join(CACHE, "detailed", profile_subfolder)):
            R += [(profile_subfolder, annotations_subfolder)]

    R = sorted(R, key=lambda x: x[0])
    return R

for profile_subfolder, annotations_subfolder in iterate_all_cache_folders():

    profile = profile_subfolder

    profile_number = int(profile[1:])
    if profile_number < 78:
        continue

    cache_folder = os.path.join(CACHE, "detailed", profile_subfolder, annotations_subfolder)

    def truncate_string(value, max_length=43, suffix="..."):
        string_value = str(value)
        if string_value.endswith("]"):
            string_value = string_value.split("[")[0]
        string_value = string_value.title()
        string_truncated = string_value[
            : min(len(string_value), (max_length - len(suffix)))
        ]
        suffix = suffix if len(string_value) > max_length else ""
        return string_truncated + suffix
    
    def running_sum(signature, geneset, library, result):
        
        signature = signature.sort_values(1, ascending=False).set_index(0)
        signature = signature[~signature.index.duplicated(keep="first")]

        signature_map = {}
        for i, h in enumerate(signature.index):
            signature_map[h] = i

        gs = set(library[geneset])
        hits = [i for i, x in enumerate(signature.index) if x in gs]

        running_sum, es = blitz.enrichment_score(
            np.array(np.abs(signature.iloc[:, 0])), signature_map, gs
        )
        running_sum = list(running_sum)

        cmap_types = stylia.colors.ContinuousColorMap(
            stylia.colors.NamedColorMaps().spectral, transformation=None
        )
        cmap_types.fit([-5,5])

        cmap_nes = stylia.colors.ContinuousColorMap(
            stylia.colors.NamedColorMaps().spectral, transformation=None
        )
        cmap_nes.fit([0, 4])

        fig = plt.figure(figsize=(2, 2))

        ax = fig.add_gridspec(12, 11, wspace=0, hspace=0.5)
        ax1 = fig.add_subplot(ax[0:7, 0:11])
        ax1 = stylize(ax1)

        color_ranksum = cmap_nes.transform([4-result.loc[geneset, "nes"]])[0]

        ax1.plot(list(running_sum), color=color_ranksum, lw=1)
        ax1.fill_between([i for i in range(len(running_sum))], list(running_sum), color=color_ranksum, alpha=0.2)
        
        ax1.tick_params(labelsize=6)
        plt.xlim([0, len(running_sum)])

        ylim = ax1.get_ylim()
        nn = np.where(np.abs(running_sum) == np.max(np.abs(running_sum)))[0][0]
        ax1.vlines(
            x=nn,
            ymin=np.min(ylim[0]),
            ymax=np.max(running_sum),
            linestyle=":",
            lw = 0.5,
            color=color_ranksum,
        )
        ax1.set_ylim(ylim)

        if es > 0:
            ax1.text(
                len(running_sum)*0.7,
                np.max(running_sum)*0.95,
                "NES = {0:.3f}\nP = {1:.2e}".format(result.loc[geneset, "nes"], result.loc[geneset, "pval"]),
                size=6,
                ha="left",
                va="top",
                zorder=100,
            )
        else:
            ax1.text(
                len(running_sum) / 30,
                0,
                "NES=" + "{:.3f}".format(result.loc[geneset, "nes"]),
                size=6,
                bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 1},
                ha="left",
                va="top",
                zorder=100,
            )

        ax1.grid(True, which="both")
        ax1.set(xticks=[])
        plt.title("{0} / {1}".format(profile_subfolder, truncate_string(geneset)), fontsize=6)

        plt.ylabel("Enrichment Score (ES)", fontsize=6)
        plt.xlabel("")
        
        rank_vec = signature[1]

        ax1 = fig.add_subplot(ax[7:8, 0:11])
        ax1 = stylize(ax1)
        stylia.label(ax1, title="", xlabel="", ylabel="")
        colors = cmap_types.transform([-rank_vec[i] for i in hits])
        ax1.vlines(x=hits, ymin=-1, ymax=1, color=colors, lw=0.7)
        plt.xlim([0, len(running_sum)])
        plt.ylim([-1, 1])
        ax1.set(yticks=[])
        ax1.set(xticks=[])
        ax1 = fig.add_subplot(ax[8:11, 0:11])
        ax1 = stylize(ax1)
        x = np.arange(0.0, len(rank_vec), 20).astype("int")
        x = np.append(x, signature.shape[0] - 1)
        ax1.fill_between(x, np.array(rank_vec)[x], color=NamedColors().gray, alpha=0.5)
        ax1.plot(x, np.array(rank_vec)[x], lw=0.5, color=NamedColors().gray)
        ax1.hlines(y=0, xmin=0, xmax=len(rank_vec), color="black", zorder=100, lw=0.5)
        plt.xlim([0, len(running_sum)])
        plt.ylim([-5.8, 5.8])
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        plt.xlabel("Ranked proteins", fontsize=6)
        plt.ylabel("Z-score", fontsize=6)
        ax1.tick_params(labelsize=6)
        stylia.label(ax1, title="")
        plt.ion()
        fig.patch.set_facecolor("white")
        return fig

    def gsea_plot(signature, term, annotations, result):
        ax = running_sum(signature, term, annotations, result=result)
        return ax

    i = 0
    j = 0
    
    MIN_EDGE_SIZE = 5
    MAX_EDGE_SIZE = 5000

    result = pd.read_csv(os.path.join(cache_folder, "result.tsv"), delimiter="\t")
    result = result[result["leading_edge_size"] >= MIN_EDGE_SIZE]
    result = result[result["leading_edge_size"] <= MAX_EDGE_SIZE]
    df = result.copy()
    df = df.sort_values("nes", ascending=False).reset_index(drop=True)
    df = df.head(50)

    print(profile_subfolder, annotations_subfolder, df.shape[0])

    with open(os.path.join(cache_folder, "paths.json"), "r") as f:
        paths = json.load(f)

    results_path = os.path.join(ROOT, "..", paths["results_path"])

    signature_ = pd.read_csv(
        os.path.join(results_path, "signature.tsv"), delimiter="\t", header=None
    )

    with open(os.path.join(results_path, "annotations.json"), "r") as f:
        annotations_ = json.load(f)

    counts = 0
    idx2term = {}

    if os.path.exists(os.path.join(cache_folder, "basic")):
        shutil.rmtree(os.path.join(cache_folder, "basic"))
    os.makedirs(os.path.join(cache_folder, "basic"), exist_ok=True)

    for term in list(df["Term"]):

        result = pd.read_csv(
            os.path.join(results_path, "result.tsv"), delimiter="\t"
        )

        result_ = result.copy()
        result_ = result_.set_index("Term")

        if j == top_plots_number:
            break
        fig, ax = plt.subplots(1,1,figsize=(5,5))
        fig = gsea_plot(signature_, term, annotations_, result=result_)
        i += 1
        j += 1

        stylia.save_figure(os.path.join(cache_folder, "basic", "plot_{0}.png".format(counts)))
        
        idx2term[counts] = term
        counts += 1

        plt.close(fig)

    with open(os.path.join(cache_folder, "basic", "idx2term.json"), "w") as f:
        json.dump(idx2term, f)