import os
import sys
import csv
import json
import collections
import pandas as pd
from tqdm import tqdm
import matplotlib
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Arial"
import numpy as np
import matplotlib.pyplot as plt
import stylia
import matplotlib.pyplot as plt
import stylia
import joblib

proteins_set = None

ROOT = os.path.abspath(os.path.dirname(__file__))
TMP = os.path.join(ROOT, "tmp")
if not os.path.exists(TMP):
    os.mkdir(TMP)

MIN_SET_SIZE = 1
PROFILE_TYPE = "Fragment"
OVERVIEW_PVALUE_CUTOFF = 0.05

# relative imports
sys.path.append(os.path.join(ROOT, "../src/"))
from util import listdir_util

# import metadata
from proteome_meta import global_profiles_dict
from proteome_meta import task_suf
from proteome_meta import annotation_type_dict
from proteome_meta import annotation_dict
from proteome_meta import universe_dict

# path to results and original data
PATH = os.path.abspath(os.path.join(ROOT, "../results/proteins/"))
DATA = os.path.abspath(os.path.join(ROOT, "../data"))
CACHE = os.path.abspath(os.path.join(ROOT, "../cache"))

# generic inputs

# protein id to gene name
def get_protein_id_to_gene_name():
    pid2name = {}
    name2pid = {}
    with open(os.path.join(DATA, "general", "pid2name_primary.tsv"), "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for r in reader:
            pid2name[r[0]] = r[1]
            name2pid[r[1]] = r[0]
    return pid2name, name2pid

pid2name, name2pid = get_protein_id_to_gene_name()

def pid2gene(x):
    if x in pid2name:
        return pid2name[x]
    else:
        return x


def gene2pid(x):
    if x in name2pid:
        return name2pid[x]
    else:
        return x


def pretty_term(x):
    x = x.title()
    if x.endswith("]"):
        x = x.split(" [")[0]
    return x

profile_type = PROFILE_TYPE
profile_type_subfolder = profile_type.lower()

def get_sorted_fids():
    fids = []
    for fid in listdir_util(os.path.join(DATA, "signatures", "proteins", "fragment")):
        fids += [fid]
    fids = sorted(fids)
    return fids

fids = get_sorted_fids()
for profile in tqdm(fids):
    profile_subfolder = profile
    all_cases = fids
    draw_fragment = True

    def universe_selector():
        preselected="HEK293T Core"
        universe = preselected
        universe_subfolder = universe_dict[universe]
        return universe, universe_subfolder

    def task_finder(signature_data):
        tasks = []
        task_filenames = {}
        for k, v in task_suf.items():
            for task in listdir_util(signature_data):
                if task.endswith(k + ".tsv"):
                    tasks += [v]
                    task_filenames[v] = task
        return tasks, task_filenames
    
    universe, universe_subfolder = universe_selector()
    universe_file = os.path.abspath(
        os.path.join(DATA, "universes", "proteins", universe_subfolder + ".tsv")
    )

    annotations_ = []

    annotation_type_labels = []
    annotation_keys = []
    for at, v in annotation_type_dict.items():
        for x in v:
            an = annotation_dict[x]
            annotation_type_labels += [(at, x, an)]
            annotation_keys += ["{0}---{1}".format(at, an)]

    signature_files = []
    results_folders = []
    annotation_files = []
    for an in annotation_type_labels:
        an = an[2]
        signature_data = os.path.join(
            DATA, "signatures", "proteins", profile_type_subfolder, profile_subfolder
        )
        tasks, task_filenames = task_finder(signature_data)
        for t in task_filenames.items():
            signature_files += [os.path.join(signature_data, t[1])]
            results_folders += [
                os.path.join(
                    PATH,
                    universe_subfolder,
                    an,
                    profile_type_subfolder,
                    profile_subfolder,
                    t[1].split(".tsv")[0],
                )
            ]
            annotation_files += [
                os.path.abspath(
                    os.path.join(DATA, "annotations", "proteins", an + ".tsv")
                )
            ]

    for sf, rf, af in zip(signature_files, results_folders, annotation_files):
        if os.path.exists(os.path.join(rf, "result.tsv")):
            pass
        else:
            print("Missing", rf)

    p_value = OVERVIEW_PVALUE_CUTOFF

    protein_label = "Gene Name"
    if protein_label == "Gene Name":
        convert_to_gene = True
    else:
        convert_to_gene = False

    direction = "Up"

    annotation_items = collections.defaultdict(list)
    for an in annotation_type_labels:
        with open(
            os.path.join(DATA, "annotations", "proteins", an[-1] + ".tsv"), "r"
        ) as f:
            reader = csv.reader(f, delimiter="\t")
            for r in reader:
                annotation_items[(an[2], r[1])] += [r[0]]

    item_annotations = collections.defaultdict(list)
    for k, v in annotation_items.items():
        for x in list(v):
            item_annotations[x] += [k]

    all_proteins = sorted(item_annotations.keys())

    proteins_in_universe = []
    with open(
        os.path.join(DATA, "universes", "proteins", universe_subfolder + ".tsv"), "r"
    ) as f:
        reader = csv.reader(f, delimiter="\t")
        for r in reader:
            proteins_in_universe += [r[0]]

    annotation_type_inv = {}
    for k, v in annotation_type_dict.items():
        for x in v:
            if x in annotation_dict:
                annotation_type_inv[x] = k

    class Overviewer(object):
        
        def __init__(self, profile_subfolder, annotations):
            self.profile_subfolder = profile_subfolder
            self.annotations = annotations
            self.anno_subtype_dict = dict((v,k) for k,v in annotation_dict.items())

        def get(self):
            R = []
            anno_type = []
            anno_subtype = []
            anno_subtype_subfolder = []
            ann_idxs = []
            ann_vals = []
            sig_sizes = []
            edge_idxs = []
            for k in self.annotations:
                results_folder = os.path.join(PATH, "hek293t_core", k.split("---")[1], "fragment", self.profile_subfolder, self.profile_subfolder+"_val_log2fc")
                path = os.path.join(results_folder, "harmon.tsv")
                df = pd.read_csv(path, delimiter="\t")
                df = df[df["pval"] < 0.05]
                df = df[df["score"] > 0]
                df = df[df["overlap"] > 5]
                df = df[df["setsize"] < 500]
                anno_type += [k.split("---")[0]]*df.shape[0]
                anno_subtype += [self.anno_subtype_dict[k.split("---")[1]]]*df.shape[0]
                anno_subtype_subfolder += [k.split("---")[1]]*df.shape[0]
                terms = list(df["term"])
                dr = pd.read_csv(os.path.join(results_folder, "result.tsv"), delimiter="\t")
                dr = dr[dr["Term"].isin(terms)][["Term", "leading_edge"]]
                dr = dr.rename(columns={"Term": "term", "leading_edge": "edge"})
                df = pd.merge(df, dr, on="term", how="left")
                R += [df]
                with open(os.path.join(results_folder, "signature.tsv"), "r") as f:
                    reader = csv.reader(f, delimiter="\t")
                    sig_val = {}
                    sig_idx = {}
                    for i, r in enumerate(reader):
                        sig_val[r[0]] = float(r[1])
                        sig_idx[r[0]] = i
                with open(os.path.join(results_folder, "annotations.json"), "r") as f:
                    ann = json.load(f)
                for r in df[["term", "edge"]].values:
                    t = r[0]
                    edge = r[1].split(",")
                    pids = ann[t]
                    ann_idxs += [",".join([str(sig_idx[pid]) for pid in pids])]
                    ann_vals += [",".join([str(sig_val[pid]) for pid in pids])]
                    edge_idxs += [",".join([str(sig_idx[e]) for e in edge])]
                sig_sizes += [len(sig_val)]*df.shape[0]
            df = pd.concat(R, axis=0)
            df["type"] = anno_type
            df["subtype"] = anno_subtype
            df["subtype_key"] = anno_subtype_subfolder
            df["edge_idxs"] = edge_idxs
            df["ann_idxs"] = ann_idxs
            df["ann_vals"] = ann_vals
            df["sig_size"] = sig_sizes
            df = df.sort_values(by="pval").reset_index(drop=True)
            return df

    ov = Overviewer(profile_subfolder, annotation_keys)
    df = ov.get()
    df = df.head(200)

    df.to_csv(os.path.join(CACHE, "overview", "{0}.tsv".format(profile)), sep="\t", index=False)

    # Now the plots

    cmap_types = stylia.colors.ContinuousColorMap(
        stylia.colors.NamedColorMaps().spectral, transformation=None
    )
    cmap_types.fit([0,1,2,3,4])
    annotation_type_colors = {}
    i = 0
    for k, v in annotation_type_dict.items():
        annotation_type_colors[k] = cmap_types.transform([i])[0]
        i += 1

    def many_ticks_plots(ax, data, fitted_cmap):
        sig_size = np.max(data["sig_size"])
        for j, r in enumerate(data[["term", "ann_vals", "ann_idxs", "edge", "edge_idxs", "subtype", "type", "score", "pval"]].values):
            idxs = [int(x) for x in r[2].split(",")]
            vals = [float(x) for x in r[1].split(",")]
            colors = fitted_cmap.transform(vals)
            for i in range(len(idxs)):
                ax.plot([idxs[i], idxs[i]], [j-0.4, j+0.4], color=colors[i], lw=0.5)
            ax.text(sig_size*0.6, j, pretty_term(r[0]), ha="center", va="center", fontsize=6)
            eidx = int(np.max([int(x) for x in r[4].split(",")]))
            ax.scatter([eidx], [j+0.4], color="black", marker="^", s=10, zorder=100000)
            edge = [pid2gene(x) for x in r[3].split(",")]
            le_prop = "{0} / {1}".format(len(edge), eidx+1)
            edge = " ".join(edge[:5]) + " ..."
            sc = "{:.2f}".format(r[7])
            pv = "{:.1e}".format(r[8])
            s = "[{0}] {1} ({2})".format(le_prop, sc, pv)
            ax.text(sig_size*0.01, j, s=edge + " " + "{0}".format(s), ha="left", va="center", fontsize=6)
            ax.text(sig_size*0.98, j, s=r[5], fontsize=6, va="center", ha="right")
            ax.scatter([sig_size*0.99], [j], color=annotation_type_colors[r[6]], s=30, edgecolor="black", lw=0.2, zorder=1000)
            ax.text(-sig_size*0.01, j, s=(j+1), fontsize=6, va="center", ha="right")

        ax.set_xlim((0, sig_size))
        ax.set_ylim(data.shape[0], -1)
        ax.set_axis_off()
        ax.set_title("Enrichment leaderboard for fragment {0}".format(profile_subfolder), fontsize=8)

    prot2idx = collections.defaultdict(list)
    for i,r in enumerate(list(df["edge"])):
        for x in r.split(","):
            gn = pid2gene(x)
            prot2idx[gn] += [i]
    all_proteins_ = sorted(prot2idx.keys())
    ann2idx = collections.defaultdict(list)
    for i,r in enumerate(df["term"]):
        ann2idx[r] += [i]
    all_annotations_ = sorted(ann2idx.keys())

    type2idx = collections.defaultdict(list)
    for i,r in enumerate(list(df["type"])):
        type2idx[r] += [i]
    all_types_ = sorted(type2idx.keys())

    subtype2idx = collections.defaultdict(list)
    for i,r in enumerate(list(df["subtype"])):
        subtype2idx[r] += [i]
    all_subtypes_ = sorted(subtype2idx.keys())

    selected_proteins = all_proteins_
    selected_annotations = all_annotations_
    selected_subtypes = all_subtypes_
    selected_types = all_types_

    keep_idxs = []
    if selected_proteins is not None:
        for x in selected_proteins:
            for idx in prot2idx[x]:
                keep_idxs += [idx]
    
    if selected_annotations is not None:
        for x in selected_annotations:
            for idx in ann2idx[x]:
                keep_idxs += [idx]

    if selected_subtypes is not None:
        for x in selected_subtypes:
            for idx in subtype2idx[x]:
                keep_idxs += [idx]
    
    if selected_types is not None:
        for x in selected_types:
            for idx in type2idx[x]:
                keep_idxs += [idx]
    
    if keep_idxs:
        keep_idxs = sorted(set(keep_idxs))
        df = df.iloc[keep_idxs]

    df["edge_genes"] = [" ".join([pid2gene(x) for x in r.split(",")]) for r in list(df["edge"])]

    df_view = df[["term", "overlap", "setsize", "score", "pval", "edge_genes", "subtype", "type"]]
    df_view = df_view.rename(columns = {
        "term": "Term",
        "overlap": "Edge size",
        "setsize": "Set size",
        "score": "Score",
        "pval": "P-value",
        "edge_genes": "Leading edge",
        "subtype": "Category subtype",
        "type": "Category type"
    })
    df_view["rank"] = [i+1 for i in range(df_view.shape[0])]
    df_view = df_view.set_index("rank")

    cmap = stylia.colors.ContinuousColorMap(
        stylia.colors.NamedColorMaps().coolwarm, transformation=None
    )
    cmap.fit([-5, 5])

    number_to_show = df_view.shape[0]

    fig, ax = plt.subplots(1,1,figsize=(10,28*number_to_show/100))
    many_ticks_plots(ax, df.head(number_to_show), cmap)
    plt.tight_layout()
    plt.savefig(os.path.join(CACHE, "overview", "{0}.png".format(profile)), dpi=300)