# regular imports
import os
import sys
import csv
import json
import pandas as pd
import collections
import numpy as np
import shutil
from sklearn.preprocessing import MinMaxScaler
import streamlit as st
import h5py
import matplotlib.pyplot as plt
import stylia
from stylia.colors.colors import NamedColors
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from griddify import Cloud2Grid


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
from enrichers import Enricher, Harmonizer, PreSummarizer, Summarizer
from util import listdir_util

# import metadata
from proteome_meta import global_profiles_dict
from proteome_meta import task_suf
from proteome_meta import annotation_type_dict
from proteome_meta import annotation_dict
from proteome_meta import universe_dict

# set page layout
st.set_page_config(layout="wide", page_title="Ligand Discovery Protein Set Enrichment Analysis")

# path to results and original data
PATH = os.path.abspath(os.path.join(ROOT, "../results/proteins/"))
DATA = os.path.abspath(os.path.join(ROOT, "../data"))

# generic inputs

# protein id to gene name
pid2name = {}
name2pid = {}
with open(os.path.join(DATA, "general", "pid2name_primary.tsv"), "r") as f:
    reader = csv.reader(f, delimiter="\t")
    for r in reader:
        pid2name[r[0]] = r[1]
        name2pid[r[1]] = r[0]


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

# side bar

st.sidebar.title("Ligand Discovery Proteome Set Enrichment Analysis")

# signatures (aka profiles)
st.sidebar.header("Select a fragment")

profile_type = PROFILE_TYPE
profile_type_subfolder = profile_type.lower()

if profile_type == "Global":
    options = ["Detectability", "Promiscuity"]
    profile = st.sidebar.radio("Global protein profile", options)
    profile_subfolder = global_profiles_dict[profile]
    all_cases = [global_profiles_dict[o] for o in options]
    draw_fragment = False
else:
    fids = []
    for fid in listdir_util(os.path.join(DATA, "signatures", "proteins", "fragment")):
        fids += [fid]
    fids = sorted(fids)
    profile = st.sidebar.selectbox("Fragment identifier", options=fids)
    profile_subfolder = profile
    all_cases = fids
    draw_fragment = True

st.sidebar.header("Choose a type of analysis")

type_of_analysis = st.sidebar.radio(
    "Type of analysis", options=["Overview", "Detailed"]
)

# Universe

def universe_selector(preselected="HEK293T Core"):
    if preselected is None:
        st.sidebar.header("Background proteomes")
        universes = [
            "Human Proteome",
            "HEK293T Core",
            "Bind Degs Detected",
            "Bind Degs Enriched",
        ]
        if profile_type == "Fragment":
            universes += ["Pulldown"]
        universe = st.sidebar.radio("Protein universe (background)", universes, index=1)
    else:
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


if type_of_analysis == "Overview":

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
            st.error("Enrichment not found! Please run the analysis first.")
            enricher = Enricher(
                signature_file=sf,
                annotations_file=af,
                universe_file=universe_file,
                min_set_size=MIN_SET_SIZE,
            )
            enricher.calculate()
            if not os.path.exists(rf):
                os.makedirs(rf, exist_ok=True)
            enricher.save(rf)

    st.header("Enrichment overview for {0} {1}".format(profile_type.lower(), profile))

    p_value = OVERVIEW_PVALUE_CUTOFF

    columns = st.columns(3)

    protein_label = "Gene Name"
    if protein_label == "Gene Name":
        convert_to_gene = True
    else:
        convert_to_gene = False

    direction = "Up"

    view = columns[0].radio("View", options=["Table", "Plot"])

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

    named_colors = NamedColors()

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

    ov = Overviewer(profile_subfolder, annotation_keys)
    df = ov.get()
    df = df.head(200)

    columns = st.columns(4)

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

    selected_proteins = columns[0].multiselect("Filter by proteins in leading edge ({0} unique proteins)".format(len(all_proteins_)), options=all_proteins_)
    selected_annotations = columns[1].multiselect("Select annotations", options=all_annotations_)
    selected_subtypes = columns[2].multiselect("Filter by annotation subtype", options=all_subtypes_)
    selected_types = columns[3].multiselect("Filter by annotation type", options=all_types_)
    
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
    if view == "Table":
        st.dataframe(df_view.reset_index(drop=True), height=2000)
    else:
        cmap = stylia.colors.ContinuousColorMap(
            stylia.colors.NamedColorMaps().coolwarm, transformation=None
        )
        cmap.fit([-5, 5])

        number_to_show = df_view.shape[0]

        fig, ax = plt.subplots(1,1,figsize=(10,28*number_to_show/100))
        many_ticks_plots(ax, df.head(number_to_show), cmap)
        st.pyplot(fig)
    
else:

    def signature_data_finder(profile_subfolder):
        signature_data = os.path.join(
            DATA, "signatures", "proteins", profile_type_subfolder, profile_subfolder
        )
        return signature_data

    def task_finder(signature_data):
        tasks = []
        task_filenames = {}
        for k, v in task_suf.items():
            for task in listdir_util(signature_data):
                if task.endswith(k + ".tsv"):
                    tasks += [v]
                    task_filenames[v] = task
        return tasks, task_filenames

    signature_data = signature_data_finder(profile_subfolder)
    tasks, task_filenames = task_finder(signature_data)

    task_type = "Consensus"
    if task_type == "Only continuous":
        task_type = "Continuous"
    if task_type == "Only binary":
        task_type = "Binary"

    def annotations_selector():
        st.sidebar.header("Select protein annotation category")

        annotation_types = [
            "Sequence",
            "Functions",
            "Processes and pathways",
            "Localization",
            "Drugs and Diseases",
        ]
        annotation_type = st.sidebar.radio("Type of annotation", annotation_types)

        annotations = annotation_type_dict[annotation_type]

        annotation = st.sidebar.selectbox("Annotation source", options=annotations)
        annotation_subfolder = annotation_dict[annotation]

        return annotation, annotation_subfolder, annotation_type, annotations

    if task_type == "Consensus":

        class EnrichmentMeasuresConsensus(object):
            def __init__(self, root_data, results_path, pvalue_cutoff=0.05):
                self.root_data = os.path.abspath(root_data)
                self.results_path = os.path.abspath(results_path)
                self.pvalue_cutoff = pvalue_cutoff

            def subfolder_tasks(self):
                sf = {}
                for l in listdir_util(self.results_path):
                    meta_path = os.path.join(self.results_path, l, "meta.json")
                    if not os.path.exists(meta_path):
                        continue
                    with open(meta_path, "r") as f:
                        data = json.load(f)
                        if data["status"] == "empty":
                            continue
                        sf[l] = data["task"]
                return sf

            def find_ranskum_subfolders(self):
                subfolders = []
                for k, v in self.subfolder_tasks().items():
                    if v == "ranksum":
                        if "log2fc" in k or "gauss" in k:
                            subfolders += [k]
                assert len(subfolders) == 1
                return subfolders

            def find_hypergeometric_subfolders(self):
                subfolders = []
                for k, v in self.subfolder_tasks().items():
                    if v == "hypergeometric":
                        subfolders += [k]
                return self.sort_hypergeometric_subfolders(subfolders)

            def sort_hypergeometric_subfolders(self, subfolders):
                hit = []
                for sf in subfolders:
                    if sf.endswith("hit"):
                        hit += [sf]
                with_digits = []
                for i, sf in enumerate(subfolders):
                    suf = sf.split("_")[-1]
                    if suf.isdigit():
                        with_digits += [(i, int(suf))]
                with_digits = sorted(with_digits, key=lambda x: x[1])
                digit = [subfolders[x[0]] for x in with_digits]
                return hit + digit

            def find_up_hypergeometric_subfolders(self):
                subfolders = [
                    x
                    for x in self.find_hypergeometric_subfolders()
                    if "_top_" in x or "_bin_" in x
                ]
                return subfolders

            def find_dw_hypergeometric_subfolders(self):
                subfolders = [
                    x for x in self.find_hypergeometric_subfolders() if "_bottom_" in x
                ]
                return subfolders

            def enricher_if_not_done(self, path):
                if not os.path.exists(os.path.join(path, "meta.json")):
                    path_ = path.split("/")
                    annotations = path_[-4]
                    universe = path_[-5]
                    entity = path_[-6]
                    profile_type = path_[-3]
                    profile = path_[-2]
                    task = path_[-1]
                    universe_file = os.path.join(
                        self.root_data, "universes", entity, universe + ".tsv"
                    )
                    signature_file = os.path.join(
                        self.root_data,
                        "signatures",
                        entity,
                        profile_type,
                        profile,
                        task + ".tsv",
                    )
                    annotations_file = os.path.join(
                        self.root_data, "annotations", entity, annotations + ".tsv"
                    )
                    enricher = Enricher(
                        signature_file=signature_file,
                        universe_file=universe_file,
                        annotations_file=annotations_file,
                        min_set_size=MIN_SET_SIZE,
                    )
                    enricher.calculate()
                    if not os.path.exists(path):
                        os.makedirs(path, exist_ok=True)
                    enricher.save(path)

            def harmonizer_if_not_done(self, path):
                if not os.path.exists(os.path.join(path, "harmon.tsv")):
                    h = Harmonizer(path)
                    h.harmonize()
                    h.save()

            def summarizer_if_not_done(self, path):
                sm = Summarizer(path)
                if not sm.is_summary_available():
                    sm.summarize()
                    sm.save()

            def get_enricher_file_path(self, path):
                return os.path.join(path, "result.tsv")

            def get_harmonizer_file_path(self, path):
                return os.path.join(path, "harmon.tsv")

            def get_summarizer_file_path(self, path):
                sm = Summarizer(path)
                return sm.get_summary_file()

            def prepare(self, subfolder):
                path = os.path.join(self.results_path, subfolder)
                self.enricher_if_not_done(path)
                self.summarizer_if_not_done(path)

            def read_data(self, subfolder):
                path = os.path.join(self.results_path, subfolder)
                self.enricher_if_not_done(path)
                self.harmonizer_if_not_done(path)
                self.summarizer_if_not_done(path)
                ds = pd.read_csv(self.get_summarizer_file_path(path), delimiter="\t")
                dr = pd.read_csv(self.get_harmonizer_file_path(path), delimiter="\t")
                df = pd.merge(dr, ds, on="term", how="left")
                return df

            def all_expected_subfolders(
                self, strings=["_gauss", "_val_log2fc", "_bin_", "_bottom_", "_top_"]
            ):
                path_ = self.results_path.split("/")
                profile = path_[-1]
                profile_type = path_[-2]
                entity = path_[-5]
                path = os.path.join(
                    self.root_data, "signatures", entity, profile_type, profile
                )
                subfolders = []
                path = os.path.abspath(path)
                for signature_file in listdir_util(path):
                    for s in strings:
                        if s in signature_file:
                            subfolders += [signature_file.replace(".tsv", "")]
                return subfolders

            def _merge(self, direction):
                subfolders = self.all_expected_subfolders()
                available_subfolders = []
                for sf in subfolders:
                    self.prepare(sf)
                    if direction == "up":
                        if "_bottom_" in sf:
                            continue
                    else:
                        if "_top_" in sf or "_bin_" in sf:
                            continue
                    available_subfolders += [sf]
                ranksum_subfolders = self.find_ranskum_subfolders()
                if direction == "up":
                    hypergeometric_subfolders = self.find_up_hypergeometric_subfolders()
                else:
                    hypergeometric_subfolders = self.find_dw_hypergeometric_subfolders()
                subfolders = ranksum_subfolders + hypergeometric_subfolders
                reference = subfolders[0]
                ref_data = self.read_data(reference)
                data = collections.OrderedDict()
                terms = pd.DataFrame(ref_data["term"])
                data[reference] = ref_data
                for sf in subfolders[1:]:
                    d = self.read_data(sf)
                    data[sf] = pd.merge(terms, d, on="term", how="left").reset_index(
                        drop=True
                    )
                return data

            def _assemble(self, direction, data):
                df = None
                for k, v in data.items():
                    if df is None:
                        df = v.copy()
                        if direction == "up":
                            df = df[df["score"] >= 0].reset_index(drop=True)
                        else:
                            df = df[df["score"] <= 0].reset_index(drop=True)
                        df = df[df["pval"] <= self.pvalue_cutoff].reset_index(drop=True)
                        df = df[df["setsize"] >= 5].reset_index(drop=True)
                    else:
                        suffix = k.split("_")[-1]
                        v = v[["term", "pval"]]
                        v.columns = ["term", "pval_{0}".format(suffix)]
                        df = df.merge(v, on="term", how="left")
                return df

            def _consensus(self):
                consensus_path = os.path.join(self.results_path, "consensus.json")

                # return if exists
                if os.path.exists(consensus_path):
                    with open(consensus_path, "r") as f:
                        data = json.load(f)
                    return data

                # hypergeometric subfolders
                up_hypergeometric_subfolders = self.find_up_hypergeometric_subfolders()
                hg_up_dh = collections.OrderedDict()
                for sf in up_hypergeometric_subfolders:
                    tag = "bin_{0}".format(sf.split("_")[-1])
                    path = os.path.join(self.results_path, sf)
                    df = pd.read_csv(self.get_harmonizer_file_path(path), delimiter="\t")
                    d = {}
                    for r in df.values:
                        d[r[0]] = r[1:]
                    hg_up_dh[tag] = d
                dw_hypergeometric_subfolders = self.find_dw_hypergeometric_subfolders()
                hg_dw_dh = collections.OrderedDict()
                for sf in dw_hypergeometric_subfolders:
                    tag = "bin_{0}".format(sf.split("_")[-1])
                    path = os.path.join(self.results_path, sf)
                    df = pd.read_csv(self.get_harmonizer_file_path(path), delimiter="\t")
                    d = {}
                    for r in df.values:
                        d[r[0]] = r[1:]
                    hg_dw_dh[tag] = d

                # ranksum subfolder

                ranksum_subfolders = self.find_ranskum_subfolders()

                path = os.path.join(self.results_path, ranksum_subfolders[0])
                dr = pd.read_csv(self.get_enricher_file_path(path), delimiter="\t")
                ds = pd.read_csv(self.get_harmonizer_file_path(path), delimiter="\t")
                dh = pd.read_csv(self.get_summarizer_file_path(path), delimiter="\t")
                df = pd.merge(ds, dh, on="term", how="left")
                sig = []
                with open(os.path.join(path, "signature.tsv"), "r") as f:
                    reader = csv.reader(f, delimiter="\t")
                    for r in reader:
                        sig += [(r[0], float(r[1]))]
                sig = sorted(sig, key=lambda x: -x[1])
                sig_idxs = {}
                for i, r in enumerate(sig):
                    sig_idxs[r[0]] = i
                sig_vals = {}
                for r in sig:
                    sig_vals[r[0]] = r[1]
                edges = collections.OrderedDict()
                for r in dr[["Term", "leading_edge"]].values:
                    v = str(r[1])
                    if v == "nan":
                        v = []
                    else:
                        v = v.split(",")
                    edges[r[0]] = v

                with open(os.path.join(path, "annotations.json"), "r") as f:
                    annotations = json.load(f)

                data = collections.OrderedDict()
                data["signature"] = {"size": len(sig), "values": sig}

                data["terms"] = collections.OrderedDict()

                # pvalue cutoff
                df = df[df["pval"] <= self.pvalue_cutoff].reset_index(
                    drop=True
                )  # TODO parametrize

                col_idxs = dict((k, i) for i, k in enumerate(list(df.columns)))

                for v in df.values:
                    v = dict((k, v[i]) for k, i in col_idxs.items())

                    t = v["term"]

                    d = collections.OrderedDict()
                    # direction
                    if v["score"] >= 0:
                        d["direction"] = "up"
                        hg_dh = hg_up_dh
                    else:
                        d["direction"] = "dw"
                        hg_dh = hg_dw_dh

                    # overlap
                    d["overlap"] = collections.OrderedDict()
                    d["overlap"]["value"] = v["overlap"]
                    d["overlap"]["mean"] = v["overlap_mean"]
                    d["overlap"]["std"] = v["overlap_std"]
                    d["overlap"]["bin"] = collections.OrderedDict()
                    for tag, w in hg_dh.items():
                        d["overlap"]["bin"][tag] = w[t][0]

                    # setsize
                    d["setsize"] = collections.OrderedDict()
                    d["setsize"]["value"] = v["setsize"]
                    d["setsize"]["mean"] = v["setsize_mean"]
                    d["setsize"]["std"] = v["setsize_std"]
                    d["setsize"]["bin"] = collections.OrderedDict()
                    for tag, w in hg_dh.items():
                        d["setsize"]["bin"][tag] = w[t][1]

                    # score
                    d["score"] = collections.OrderedDict()
                    d["score"]["value"] = v["score"]
                    d["score"]["mean"] = v["score_mean"]
                    d["score"]["std"] = v["score_std"]
                    d["score"]["bin"] = collections.OrderedDict()
                    for tag, w in hg_dh.items():
                        d["score"]["bin"][tag] = w[t][2]

                    # pval
                    d["pval"] = collections.OrderedDict()
                    d["pval"]["value"] = v["pval"]
                    d["pval"]["mean"] = v["pval_mean"]
                    d["pval"]["std"] = v["pval_std"]
                    d["pval"]["bin"] = collections.OrderedDict()
                    for tag, w in hg_dh.items():
                        d["pval"]["bin"][tag] = w[t][3]

                    # signature idxs
                    d["signature_idxs"] = collections.OrderedDict()
                    d["signature_idxs"]["edge"] = [sig_idxs[k] for k in edges[t]]
                    d["signature_idxs"]["full"] = [sig_idxs[k] for k in annotations[t]]

                    # include
                    data["terms"][t] = d

                with open(self.get_consensus_path(), "w") as f:
                    json.dump(data, f)

                return data

            def run_merge(self, direction):
                data = self._merge(direction)
                data = self._assemble(direction, data)
                return data

            def run_consensus(self):
                data = self._consensus()
                return data

            def get_consensus_path(self):
                return os.path.join(self.results_path, "consensus.json")

            def clean(self):
                path = self.get_consensus_path()
                if os.path.exists(path):
                    os.remove(self.get_consensus_path())

        annotation, annotation_subfolder, annotation_type, annotations = (
            annotations_selector()
        )
        universe, universe_subfolder = universe_selector()
        
        path_ = os.path.join(
            "results",
            "proteins",
            universe_subfolder,
            annotation_subfolder,
            profile_type_subfolder,
            profile_subfolder,
        )
        ec = EnrichmentMeasuresConsensus("data", path_)
        ec.clean()
        df_merge = ec.run_merge("up")

        def read_annotations(path):
            found = False
            for p in os.listdir(path):
                if p.endswith("log2fc") or p.endswith("gauss"):
                    found = True
                    break
            assert found
            with open(os.path.join(path, p, "annotations.json"), "r") as f:
                ann = json.load(f)
            return ann

        def read_ranksum_harmonized(path):
            found = False
            for p in os.listdir(path):
                if p.endswith("log2fc") or p.endswith("gauss"):
                    found = True
                    break
            assert found
            return pd.read_csv(os.path.join(path, p, "harmon.tsv"), delimiter="\t")

        def read_hit_signature(path):
            found = False
            for p in os.listdir(path):
                if p.endswith("_hit"):
                    found = True
                    break
            if not found:
                return None
            else:
                path = os.path.join(path, p, "signature.tsv")
                if not os.path.exists(path):
                    return None
                with open(path, "r") as f:
                    reader = csv.reader(f, delimiter="\t")
                    prots = []
                    for r in reader:
                        prots += [r[0]]
                return set(prots)

        def read_group_data(path):
            for p in os.listdir(path):
                if p.endswith("log2fc"):
                    break
            path = os.path.join(path, p, "group_data.json")
            if os.path.exists(path):
                with open(path, "r") as f:
                    group_data = json.load(f)
            else:
                group_data = None
            return group_data

        ann = read_annotations(path_)
        harm = read_ranksum_harmonized(path_)
        hit_signature = read_hit_signature(path_)
        group_data = read_group_data(path_)
        
        from sklearn.feature_extraction.text import TfidfVectorizer
        from umap import UMAP
        from griddify import Cloud2Grid

        @st.cache_data
        def annotations_grid_as_dict(ann, harm):
            df = harm.sort_values(["setsize", "pval", "overlap"], ascending=[False, True, False]).reset_index(drop=True)
            df = df.head(2048)
            terms = list(df["term"])
            corpus = []
            for t in terms:
                corpus += [" ".join(ann[t])]
            vec = TfidfVectorizer(max_features=10000)
            X = vec.fit_transform(corpus).todense()
            red = UMAP()
            X = red.fit_transform(X, y=list(np.clip(df["score"], -5, 5)))
            X = MinMaxScaler().fit_transform(X)
            cg = Cloud2Grid(max_side=16)
            cg.fit(X)
            X = cg.transform(X)
            term_locs = {}
            for i in range(len(terms)):
                term = terms[i]
                x,y = X[i]
                term_locs[term] = (x,y)
            side = cg._side
            return term_locs, side
                    
        term_locs, term_grid_side = annotations_grid_as_dict(ann, harm)

        def read_universe_proteins(universe_subfolder):
            universe_file = os.path.join(DATA, "universes", "proteins", universe_subfolder+".tsv")
            with open(universe_file, "r") as f:
                reader = csv.reader(f, delimiter="\t")
                prots = []
                for r in reader:
                    prots += [r[0]]
                return set(prots)

        universe_proteins = read_universe_proteins(universe_subfolder)

        consensus = ec.run_consensus()

        @st.cache_data
        def proteome_grid_as_dict(grid=True):
            def read_projection(h5_file):
                with h5py.File(h5_file, "r") as f:
                    proteins = [x.decode("utf-8") for x in f["Keys"][:]]
                    X = f["Values"][:]
                return X, proteins

            X, proteins = read_projection(
                os.path.join(
                    ROOT, "../data/general/projections/esm1b_sequence.h5"
                )
            )

            idxs = [i for i,p in enumerate(proteins) if p in universe_proteins]
            X = X[idxs]
            proteins = list(np.array(proteins)[idxs])

            X = MinMaxScaler().fit_transform(X)

            if grid:
                cg = Cloud2Grid(max_side=16)
                cg.fit(X)
                X = cg.transform(X)
            else:
                G = []
                for i in np.linspace(0,1,64):
                    for j in np.linspace(0,1,64):
                        G += [[i,j]]
                G = np.array(G)

                from sklearn.neighbors import NearestNeighbors

                nn = NearestNeighbors(n_neighbors=1)
                nn.fit(G)
                I = nn.kneighbors(X, return_distance=False)

                X = G[I[:,0]]

            prot_locs = {}
            for i in range(len(proteins)):
                pid = proteins[i]
                x,y = X[i]
                prot_locs[pid] = (x,y)
            
            return prot_locs

        prot_locs = proteome_grid_as_dict(grid=True)

        def annotations_grid_plot(ax, term, consensus, data, term_locs, term_grid_side):
            scores = {}
            for r in data[["term", "score"]].values:
                scores[r[0]] = float(r[1])

            n = term_grid_side
            n = (n-2)*2+1
            radius = (1/(n+1))

            loc_terms = collections.defaultdict(list)
            for k,v in term_locs.items():
                x,y = v
                loc_terms[(x,y)] += [k]

            cmap = stylia.colors.ContinuousColorMap(
                stylia.colors.NamedColorMaps().coolwarm, transformation=None
            )
            cmap.fit([-5, 5])
            for k,v in loc_terms.items():
                sc = np.array([scores[x] for x in v])
                sc = np.mean(sc)
                color = cmap.transform([sc])[0]
                patch = Circle((k[0], k[1]), radius, facecolor=color, edgecolor="none")
                ax.add_patch(patch)

            sc = consensus["terms"][term]["score"]["value"]
            x,y = term_locs[term]
            ax.plot([x, x], [-10, 10], lw=0.5, color="black")
            ax.plot([-10, 10], [y, y], lw=0.5, color="black")
            ax.scatter([x], [y], color=cmap.transform([sc])[0], s=10, zorder=1000, lw=0.5, edgecolor="white")
            labels = [i for i in range(term_grid_side)]
            labels_norm = dict((k, i) for i,k in enumerate(np.linspace(0,1,term_grid_side)))
            ax.set_xticks([x])
            ax.set_xticklabels([labels[labels_norm[x]]])
            ax.set_yticks([y])
            ax.set_yticklabels([labels[labels_norm[y]]])

            ax.set_xlim(0-radius-0.05,1+radius+0.05)
            ax.set_ylim(0-radius-0.05,1+radius+0.05)

            ax.set_title("Annotation")

            ax.grid(False)


        def leading_edge_promiscuity_analysis(ax, term, consensus):
            data = consensus["terms"][term]
            sig = consensus["signature"]["values"]
            le_idxs = data["signature_idxs"]["edge"]
            le_pids = set([sig[i][0] for i in le_idxs])
            prom_prots_sorted = [x[0] for x in sorted(protein_promiscuity.items(), key=lambda x: -x[1])]

            cmap = stylia.colors.ContinuousColorMap(
                stylia.colors.NamedColorMaps().spectral, transformation=None
            )
            cmap.fit([-3, 3])

            def add_bar(top, i):
                prom = set(prom_prots_sorted[:top])
                x = prom.intersection(le_pids)
                v = protein_promiscuity[prom_prots_sorted[top-1]]
                color = cmap.transform([v])[0]
                x = len(x)/len(le_pids)
                ax.scatter([x], [i], color=color)
                ax.plot([0,x], [i, i], color=color)
                
            
            tops = [25, 100, 250, 500]
            for i, top in enumerate(tops):
                add_bar(top, i)

            ax.set_xlim(-0.1, 1.1)
            ax.set_ylim(3.5, -0.5)

            ax.set_yticks([0,1,2,3])
            ax.set_yticklabels([25, 100, 250, 500])

            ax.set_xticks([0, 0.5, 1])
            ax.set_xlabel("Proportion")
            ax.set_title("Edge prom.")
        

        def consensus_ticks_plot(ax, term, consensus):
            sig = consensus["signature"]["values"]
            data = consensus["terms"][term]
            idxs = data["signature_idxs"]["full"]
            cmap = stylia.colors.ContinuousColorMap(
                stylia.colors.NamedColorMaps().coolwarm, transformation=None
            )
            cmap.fit([-5, 5])
            values = [float(sig[i][1]) for i in idxs]
            colors = cmap.transform(values)
            for i, idx in enumerate(idxs):
                ax.plot([idx, idx], [0, 1], color=colors[i], lw=1)

            # add binary dots
            
            for k,v in data["pval"]["bin"].items():
                loc = k.split("_")[-1]
                if loc.isnumeric():
                    loc = int(loc)-1
                    v = float(v)
                    if v < 0.05: # TODO parametrize
                        color = "black"
                    else:
                        color = "white"
                    ax.scatter([loc], [0.5], color=color, edgecolor="black", zorder=100000, lw = 0.5)

            # ranksum leading edge
            
            le = np.max(data["signature_idxs"]["edge"])-1             
            ax.scatter([le], [0.01], marker="^", color="black", edgecolor="black", zorder=100000)

            pad = 0.01
            
            xlim = (0, len(sig))
            
            xrng = xlim[1] - xlim[0]
            ylim = (0,1)
            yrng = ylim[1] - ylim[0]
            ylim = (ylim[0]-yrng*0.1, ylim[1]+yrng*0.1)
            xlim = (xlim[0]-xrng*pad, xlim[1]+xrng*pad)

            sc = "{:.2f}".format(data["score"]["value"])
            pv = "{:.1e}".format(data["pval"]["value"])

            s = "NES = {0} (P = {1})".format(sc, pv)
            ax.text(x=np.mean(xlim), y=np.mean(ylim), s=s, ha="center", va="center")

            ax.axes.yaxis.set_visible(False)
            ax.set_title(pretty_term(term)[:60])
            
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax.grid(False)

            ax.set_xlabel("Ranked list")

            p = -xlim[1]/25

            e = len(data["signature_idxs"]["edge"])
            n = len(data["signature_idxs"]["full"])

            ax.text(p, 0.7, e, va="center", ha="center", fontsize=6)
            ax.plot([p*0.7, p*1.3], [0.5, 0.5], lw = 0.5, color="black", clip_on=False)
            ax.text(p, 0.25, n, va="center", ha="center", fontsize=6)

            return ax

        def expected_value_plot(ax, data, term, groups=False):

            cmap = stylia.colors.ContinuousColorMap(
                stylia.colors.NamedColorMaps().coolwarm, transformation=None
            )
            cmap.fit([-5, 5])

            v = data["score"]["value"]
            m = data["score"]["mean"]
            s = data["score"]["std"]

            rect = Rectangle(xy=(m-s,-1), width=s*2, height=10, facecolor=cmap.transform([m])[0], alpha=0.5)

            if groups and group_data is not None:
                if term in group_data:
                    d = group_data[term]["other"]
                    osc = []
                    for _,x in d.items():
                        osc += [x[0]]
                    osc = np.clip(osc, -5, 5)
                    ax.scatter(osc, [0.3]*len(osc), color="black", s=2, zorder=10000)

            ax.scatter(np.clip([v], -5, 5), [0.5], color=cmap.transform([v]), zorder=10000, s=5)
            ax.add_patch(rect)

            ax.plot([m,m], [-10,10], color=cmap.transform([m])[0], zorder=1000)

            ax.axes.yaxis.set_visible(False)
            
            ax.set_xlim(-5.5,5.5)
            ax.set_ylim(-0.1,1.1)

            ax.set_title("Exp. NES")
            ax.grid(False)

            colors = cmap.transform(np.linspace(-5,5,10))
            for i,y in enumerate(np.linspace(0,1,len(colors))):
                ax.scatter([-7.5], [y], color=colors[i], clip_on=False, s=5)

            for y, l in zip([0, 0.5, 1], ["-5", "0", "+5"]):
                ax.text(-9.2, y, l, va="center", ha="right", fontsize=6)

            return ax

        def get_promiscuity():
            protein_promiscuity = {}
            with open(os.path.join(DATA, "signatures", "proteins", "global", "promiscuity", "promiscuity_gauss.tsv"), "r") as f:
                reader = csv.reader(f, delimiter="\t")
                for r in reader:
                    protein_promiscuity[r[0]] = np.clip(float(r[1]), -3, 3)
            return protein_promiscuity

        def get_detectability():
            protein_detectability = {}
            with open(os.path.join(DATA, "signatures", "proteins", "global", "detectability", "detectability_gauss.tsv"), "r") as f:
                reader = csv.reader(f, delimiter="\t")
                for r in reader:
                    protein_detectability[r[0]] = np.clip(float(r[1]), -3, 3)
            return protein_detectability

        protein_promiscuity = get_promiscuity()
        protein_detectability = get_detectability()

        def leading_edge(ax, term, consensus):
            data = consensus["terms"][term]
            sig = consensus["signature"]["values"]
            le_idxs = set(data["signature_idxs"]["edge"])
            fl_idxs = set(data["signature_idxs"]["full"])
            sel = []
            for i, s in enumerate(sig):
                if i in fl_idxs:
                    if i in le_idxs:
                        l = 1
                    else:
                        l = 0
                    sel += [(s[0], pid2gene(s[0]), float(s[1]), l)]

            top = 10    
            for i in range(top):
                ax.scatter([i], [0.5], s=1, color="black")

            sel = sel[:top]
            
            cmap = stylia.colors.ContinuousColorMap(
                stylia.colors.NamedColorMaps().coolwarm, transformation=None
            )
            cmap.fit([-5, 5])

            cmap_prom = stylia.colors.ContinuousColorMap(
                stylia.colors.NamedColorMaps().spectral, transformation=None
            )

            cmap_prom.fit([-3, 3])

            for i, s in enumerate(sel):

                pid = s[0]
                if pid not in protein_promiscuity:
                    prom = -3
                else:
                    prom = protein_promiscuity[pid]
                
                color_prom = cmap_prom.transform([prom])[0]
                
                el = i-0.4
                er = i+0.4
                rng = er-el

                ax.plot([el, er], [0.5, 0.5], lw=0.5, color=color_prom, zorder=-1)
                
                prom_scaled = (prom - (-3))/(3 - (-3))

                x = prom_scaled*rng + el
                
                if s[3] == 1:
                    ax.scatter([x], [0.5], color=color_prom, edgecolor=color_prom, lw=0.5)
                else:
                    ax.scatter([x], [0.5], color="white", edgecolor=color_prom, lw=0.5)
                if int(i/2) == i/2:
                    y = 1.5
                else:
                    y = -1
                ax.text(x=i, y=y, s=s[1], ha="center", va="center", rotation=0, fontsize=6)
            
            ax.set_ylim(-1,3)
            ax.set_xlim(-0.5,9.5)
            ax.set_xlabel("Top 10 proteins")
            ax.axes.yaxis.set_visible(False)
            ax.axes.xaxis.set_visible(False)
            ax.set_axis_off()

            footer = "{0} {1} / {2}".format(profile_type, profile, annotation)
            ax.text(4.5,-3.7, footer, ha="center", va="center", fontsize=6)

        def proteome_grid_plot(ax, term, consensus, protein_locations):

            idxs = consensus["terms"][term]["signature_idxs"]["edge"]
            sig = [x[0] for x in consensus["signature"]["values"]]
            prots = [sig[i] for i in idxs]
            
            prot_locs = protein_locations
            
            loc_prots = collections.defaultdict(list)
            for k,v in prot_locs.items():
                loc_prots[v] += [k]

            cmap_prom = stylia.colors.ContinuousColorMap(
                stylia.colors.NamedColorMaps().spectral, transformation=None
            )

            cmap_prom.fit([-3, 3])

            cmap = stylia.colors.ContinuousColorMap(
                plt.get_cmap("Greys"), transformation=None
            )
            cmap.fit([len(v) for k,v in loc_prots.items()])

            locs = sorted(loc_prots.keys())

            for loc in locs: 
                color = cmap.transform([len(loc_prots[loc])])
                ax.scatter([loc[0]], [loc[1]], s=0.5, color=color, alpha=0.5)

            my_coords = collections.defaultdict(list)
            for p in prots:
                x,y = prot_locs[p]
                my_coords[(x,y)] += [p]
            
            def sizer(n):
                factor = 8
                return n*factor

            xs = []
            ys = []
            ss = []
            cs = []
            for k,v in my_coords.items():
                s = sizer(len(v))
                proms = []
                for x in v:
                    if x in protein_promiscuity:
                        proms += [protein_promiscuity[x]]
                    else:
                        proms += [-3]
                prom = np.mean(proms)
                color = cmap_prom.transform([prom])[0]
                xs += [k[0]]
                ys += [k[1]]
                ss += [s]
                cs += [color]
            
            idxs = np.argsort(ss)[::-1]
            xs = [xs[i] for i in idxs]
            ys = [ys[i] for i in idxs]
            ss = [ss[i] for i in idxs]
            cs = [cs[i] for i in idxs]

            ax.scatter(xs, ys, s=ss, color=cs, edgecolor="white", lw=0.5, clip_on=False)

            ax.set_xlim(-0.05, 1.05)
            ax.set_ylim(-0.05, 1.05)

            pos = [0.1, 0.3, 0.55, 0.9]
            for i, n in enumerate([1, 2, 5, 10]):
                s = sizer(n)
                ax.scatter([pos[i]], [-0.6], s=s, color="none", edgecolor="black", clip_on=False, lw=0.5)

            ax.axes.yaxis.set_visible(False)
            ax.tick_params(axis=u'both', which=u'both',length=0)
            ax.set_xticks(pos)
            ax.set_xticklabels(["1", "2", "5", "10"], color="black")

            colors = cmap_prom.transform(np.linspace(-3,3,10))
            for i,y in enumerate(np.linspace(0,1,len(colors))):
                ax.scatter([-0.25], [y], color=colors[i], clip_on=False, s=5)

            for y, l in zip([0, 0.5, 1], ["-3", "0", "+3"]):
                ax.text(-0.4, y, l, va="center", ha="right", fontsize=6)

            ax.set_title("Sequence")
            
            ax.grid(False)

        def rankscore_promiscuity_scatter(ax, term, consensus):
            prot_vals_list = consensus["signature"]["values"]
            prot_vals = dict((r[0], float(r[1])) for r in prot_vals_list)
            avail_prots = set([k for k in prot_vals.keys()])

            proms = collections.OrderedDict()
            for k,v in protein_promiscuity.items():
                if k not in avail_prots:
                    continue
                proms[k] = v
            for k,v in protein_detectability.items():
                if k not in avail_prots:
                    continue
                if k in proms:
                    continue
                proms[k] = -3
            for k in list(avail_prots):
                if k in proms:
                    continue
                proms[k] = -3

            proteins = list(avail_prots)
            prom_rnks = collections.OrderedDict()
            i = 0
            for k,v in proms.items():
                prom_rnks[k] = i
                i += 1
            
            xs = []
            ys = []
            zs = []
            for p in proteins:
                xs += [prot_vals[p]]
                ys += [prom_rnks[p]]
                zs += [proms[p]]

            from sklearn.preprocessing import QuantileTransformer
            tr = QuantileTransformer(output_distribution="normal")
            ys = tr.fit_transform(-np.array(ys).reshape(-1,1))[:,0]

            ax.scatter(xs, ys, s=1, color="lightgray", alpha=0.2)

            # map proteins of interest

            edge = [prot_vals_list[i][0] for i in consensus["terms"][term]["signature_idxs"]["edge"]]
            full = [prot_vals_list[i][0] for i in consensus["terms"][term]["signature_idxs"]["full"]]

            protein_idxs = dict((k,i) for i,k in enumerate(proteins))

            cmap_prom = stylia.colors.ContinuousColorMap(
                stylia.colors.NamedColorMaps().spectral, transformation=None
            )
            cmap_prom.fit([-3, 3])
            
            for p in edge:
                idx = protein_idxs[p]
                x = xs[idx]
                y = ys[idx]
                z = zs[idx]
                color = cmap_prom.transform([z])[0]
                ax.scatter([x], [y], color=color, edgecolor=color, lw=0.5, zorder=1000)

            for p in full:
                idx = protein_idxs[p]
                x = xs[idx]
                y = ys[idx]
                z = zs[idx]
                color = cmap_prom.transform([z])[0]
                ax.scatter([x], [y], color="white", edgecolor=color, lw=0.5, zorder=900)                            

            ax.set_xlim(-5.5, 5.5)
            ax.set_ylim(-3.3, 3.3)

            ax.set_xticks([-5, -2.5, 0, 2.5, 5])
            ax.set_yticks([-3, -1.5, 0, 1.5, 3])

            if profile_type_subfolder == "fragment":
                ax.set_xlabel("Log2FC")
            else:
                ax.set_xlabel("Z-score")
            ax.set_ylabel("Promiscuity (z-score)")

            ax.set_title("Promiscuity")

        def consensus_plot(term, consensus):
            
            fig = plt.figure(constrained_layout=True, figsize=(stylia.TWO_COLUMNS_WIDTH, stylia.TWO_COLUMNS_WIDTH/3.8), dpi=600)
            spec = gridspec.GridSpec(ncols=10, nrows=2)
            
            ax = fig.add_subplot(spec[0,2:8])
            ax = stylia.figure.figure.stylize(ax)
            stylia.label(ax, title="", xlabel="", ylabel="")
            consensus_ticks_plot(ax, term, consensus)

            ax = fig.add_subplot(spec[0,8])
            ax = stylia.figure.figure.stylize(ax)
            stylia.label(ax, title="", xlabel="", ylabel="")
            expected_value_plot(ax, consensus["terms"][term], term, groups=True)
            
            ax = fig.add_subplot(spec[0,9])
            ax = stylia.figure.figure.stylize(ax)
            stylia.label(ax, title="", xlabel="", ylabel="")
            annotations_grid_plot(ax, term, consensus, harm, term_locs, term_grid_side)

            ax = fig.add_subplot(spec[1,2:8])
            ax = stylia.figure.figure.stylize(ax)
            stylia.label(ax, title="", xlabel="", ylabel="")
            leading_edge(ax, term, consensus)

            ax = fig.add_subplot(spec[1,8])
            ax = stylia.figure.figure.stylize(ax)
            stylia.label(ax, title="", xlabel="", ylabel="")
            leading_edge_promiscuity_analysis(ax, term, consensus)
            
            ax = fig.add_subplot(spec[1,9])
            ax = stylia.figure.figure.stylize(ax)
            stylia.label(ax, title="", xlabel="", ylabel="")
            proteome_grid_plot(ax, term, consensus, prot_locs)
            
            ax = fig.add_subplot(spec[:2,:2])
            ax = stylia.figure.figure.stylize(ax)
            rankscore_promiscuity_scatter(ax, term, consensus)

            return fig

        df_view = df_merge[["term", "overlap", "setsize", "score", "pval", "score_mean"]]
        df_view = df_view.rename(columns={
            "term": "Term",
            "overlap": "Edge size",
            "setsize": "Set size",
            "score": "Score",
            "pval": "P-value",
            "score_mean": "Expected score"
        })

        st.header("Fragment: {0} & Category: {2} ({1})".format(profile_subfolder, annotation_type, annotation))
        
        if True:

            binary_tasks = set()
            for t in tasks:
                tfn = task_filenames[t].split("/")[-1]
                if "_top_" in tfn or "_bottom_" in tfn or "_bin_" in tfn:
                    binary_tasks.update([t])

            task_type = "Continuous"
            if task_type == "Continuous":
                task_options = [t for t in tasks if t not in binary_tasks]
            else:
                task_options = [t for t in tasks if t in binary_tasks]
            task = "Log2FC"
            task_filename = task_filenames[task]
            universe, universe_subfolder = universe_selector()

            # Define relevant paths

            results_path = os.path.abspath(
                os.path.join(
                    PATH,
                    universe_subfolder,
                    annotation_subfolder,
                    profile_type_subfolder,
                    profile_subfolder,
                    task_filename.split(".tsv")[0],
                )
            )
            signature_file = os.path.abspath(
                os.path.join(
                    DATA,
                    "signatures",
                    "proteins",
                    profile_type_subfolder,
                    profile_subfolder,
                    task_filename,
                )
            )
            annotations_file = os.path.abspath(
                os.path.join(DATA, "annotations", "proteins", annotation_subfolder + ".tsv")
            )
            universe_file = os.path.abspath(
                os.path.join(DATA, "universes", "proteins", universe_subfolder + ".tsv")
            )

            # Check if results exist

            def do_calculations(results_path, signature_file, min_set_size):
                with st.spinner(
                    "Results for {} not available. Calculating. Please be patient...".format(
                        os.path.basename(results_path)
                    )
                ):
                    enricher = Enricher(
                        signature_file, annotations_file, universe_file, min_set_size
                    )
                    enricher.calculate()
                    if not os.path.exists(results_path):
                        os.makedirs(results_path, exist_ok=True)
                    enricher.save(results_path)

            if not os.path.exists(os.path.join(results_path, "meta.json")):
                do_calculations(results_path, signature_file, min_set_size=MIN_SET_SIZE)

            if not os.path.exists(os.path.join(results_path, "meta.json")):
                st.error(
                    "Results could not be calculated for your chosen configuration. Please report issue to mduran-frigola@cemm.at"
                )

            with open(os.path.join(results_path, "meta.json"), "r") as f:
                meta = json.load(f)

            if meta["status"] == "empty":
                st.warning("Profile was empty, no enrichments were calculated")

            if meta["status"] == "done":

                type_of_task = meta["task"]

                # read main results
                result = pd.read_csv(
                    os.path.join(results_path, "result.tsv"), delimiter="\t"
                )

                result_ = result.copy()
                result_ = result_.set_index("Term")

                # defining good ranges

                if type_of_task == "ranksum":
                    leading_edge_sizes = []
                    for l in list(result["leading_edge"]):
                        if str(l) == "nan":
                            leading_edge_sizes += [0]
                        else:
                            leading_edge_sizes += [len(l.split(","))]
                    result["leading_edge_size"] = leading_edge_sizes
                    set_sizes = result["geneset_size"]
                    common_sizes = result["leading_edge_size"]
                else:
                    set_sizes = result["geneset_size"]
                    common_sizes = result["overlap"]

                lim_set_size = (0, int(np.max(set_sizes)))
                default_min_set_size = lim_set_size[0]
                default_max_set_size = lim_set_size[1]
                default_set_size = (default_min_set_size, default_max_set_size)

                lim_common_size = (0, int(np.max(common_sizes)))
                default_min_common_size = lim_common_size[0]
                default_max_common_size = lim_common_size[1]
                default_common_size = (default_min_common_size, default_max_common_size)

                # Make sure that we have enough to do a summary (i.e. pre-summarize)
                signatures_folder = "/".join(signature_file.split("/")[:-1])
                ps = PreSummarizer(
                    results_folder=results_path,
                    signatures_folder=signatures_folder,
                    min_set_size=MIN_SET_SIZE,
                )
                ps.run(overwrite=True, minimum_number_of_cases=30)

                # check if summary is done
                if ps.is_summary_available():
                    sm = Summarizer(results_path)
                    sm.summarize()
                    sm.save()
                else:
                    print("Summary exists")

                # read signature
                signature_ = pd.read_csv(
                    os.path.join(results_path, "signature.tsv"), delimiter="\t", header=None
                )
                signature_size = signature_.shape[0]

                # read annotations
                with open(os.path.join(results_path, "annotations.json"), "r") as f:
                    annotations_ = json.load(f)
                annotations_size = len(annotations_)

                # signature and annotations metric
                metric_cols = st.columns(3)
                metric_cols[0].metric(
                    "{0} profile: {1}".format(profile_type, profile),
                    value="{0} proteins".format(signature_size),
                )
                metric_cols[1].metric(
                    "{0}: {1}".format(annotation_type, annotation),
                    value="{0} categories".format(annotations_size),
                )
                if type_of_task == "hypergeometric":
                    universe_size = pd.read_csv(
                        os.path.join(results_path, "universe.tsv"),
                        delimiter="\t",
                        header=None,
                    ).shape[0]
                    metric_cols[2].metric(
                        "Background size", value="{0} proteins".format(universe_size)
                    )

                else:
                    average_annotations_per_term = np.mean([len(v) for k,v in annotations_.items()])
                    metric_cols[2].metric(
                        "Average number of proteins per term", value="{0:.1f} proteins".format(average_annotations_per_term)
                    )

                columns = st.columns(6)
                view = columns[0].radio("View", options=["Tables", "Basic plots", "Advanced plots"])

                p_value_cutoff = columns[2].number_input("P-value cutoff", value=0.05, min_value=0., max_value=1., format="%.3f")
                min_edge_size = columns[3].number_input("Minimum leading edge size", value=5, min_value=0, max_value=10000)
                max_edge_size = columns[4].number_input("Maximum leading edge size", value=5000, min_value=1, max_value=10000)
                protein_label = "Gene Name"
                if protein_label == "Gene Name":
                    convert_to_gene = True
                else:
                    convert_to_gene = False

                # Annotations of interest and Proteins of interest

                all_annotations = sorted(annotations_.keys())

                select_columns = st.columns(3)
                selected_annotations = select_columns[2].multiselect(
                    "Select annotation categories", options=all_annotations
                )

                proteins_in_signature = set(signature_[0])
                available_proteins = sorted(
                    set(
                        [
                            x
                            for k, v in annotations_.items()
                            for x in v
                            if x in proteins_in_signature
                        ]
                    )
                )

                convert_to_gene = True
                if convert_to_gene:
                    available_proteins = sorted([pid2gene(x) for x in available_proteins])
                else:
                    available_proteins = available_proteins
                selected_proteins = select_columns[0].multiselect(
                    "Filter by proteins found in at least one annotation term ({0})".format(
                        len(available_proteins)
                    ),
                    options=available_proteins,
                )
                if selected_proteins:
                    if convert_to_gene:
                        selected_proteins = [gene2pid(x) for x in selected_proteins]
                    selected_proteins = set(selected_proteins)
                    if not selected_annotations:
                        for k, v in annotations_.items():
                            if len(selected_proteins.intersection(v)) > 0:
                                selected_annotations += [k]
                    if not selected_annotations:
                        st.warning(
                            "No available annotations for any of your proteins of interest..."
                        )

                result = result[result["nes"] > 0]
                result = result[result["leading_edge_size"] >= min_edge_size]
                result = result[result["leading_edge_size"] <= max_edge_size]
                result = result.reset_index(drop=True)
                
                prot2idx = collections.defaultdict(list)
                for i, r in enumerate(list(result["leading_edge"])):
                    if str(r) == "nan":
                        continue
                    for x in r.split(","):
                        prot2idx[pid2gene(x)] += [i]

                selected_leading_proteins = select_columns[1].multiselect(
                    "Filter by proteins found in at least one leading edge ({0})".format(len(prot2idx.keys())),
                    options = sorted(prot2idx.keys())
                )

                if selected_leading_proteins:
                    idxs = []
                    for v in selected_leading_proteins:
                        for x in prot2idx[v]:
                            idxs += [x]
                    idxs = sorted(set(idxs))
                    result = result.iloc[idxs]
                     
                # Type specific buttons

                if type_of_task == "ranksum":

                    sort_by = "NES"
                    if sort_by == "NES":
                        sort_by_nes = True
                    else:
                        sort_by_nes = False

                    direction = "Up"
                    if direction == "Up":
                        is_up = True
                    else:
                        is_up = False

                    df = result.copy()
                    df = df.rename(columns = {"Term": "term"})

                    df_merge = df_merge[["term", "score_mean"]]

                    df = df.merge(df_merge, how="left", on="term")

                    df = df[df["leading_edge"].notnull()]

                    df["edge_genes"] = [" ".join([pid2gene(x) for x in r.split(",")]) for r in list(df["leading_edge"])]

                    df = df[["term","leading_edge_size",  "geneset_size", "nes", "pval", "fdr", "score_mean", "edge_genes", "leading_edge"]]

                    if selected_annotations:
                        df = df[df["term"].isin(selected_annotations)]

                    if is_up:
                        df = df[df["nes"] >= 0]
                    else:
                        df = df[df["nes"] < 0]
                    if sort_by_nes:
                        if is_up:
                            df = df.sort_values(by="nes", ascending=False)
                        else:
                            df = df.sort_values(by="nes", ascending=True)
                    else:
                        df = df.sort_values(by="pval")

                    df = df.reset_index(drop=True)

                    df = df.rename(columns = {
                        "term": "Term",
                        "leading_edge_size": "Edge size",
                        "geneset_size": "Set size",
                        "nes": "Score",
                        "pval": "P-value",
                        "fdr": "FDR",
                        "score_mean": "Mean score",
                        "edge_genes": "Leading edge",
                    })
                    
                    st.write("Rank-Sum enrichment results")
                    if view == "Tables":
                        st.dataframe(df[[c for c in list(df.columns)[:-1] if c != "Mean score"]].reset_index(drop=True))
                        
                        term = st.selectbox("Explore term...", df["Term"])

                        if term is not None:
                            # Explore term

                            t_values = {}
                            for r in signature_.values:
                                t_values[r[0]] = r[1]
                            o_values = {}
                            signature_original = pd.read_csv(
                                signature_file, delimiter="\t", header=None
                            )
                            for r in signature_original.values:
                                o_values[r[0]] = r[1]

                            cols = st.columns([0.15, 1])

                            col = cols[0]

                            annotations_size = len(annotations_[term])
                            signature_size = len(signature_)

                            df_filt = df[df["Term"] == term]
                            leading_edge = list(df_filt["leading_edge"])[0]
                            if str(leading_edge) == "nan":
                                leading_edge = []
                            else:
                                leading_edge = leading_edge.split(",")
                            display_proteins = col.radio(
                                "Display proteins",
                                [
                                    "Leading edge ({0})".format(len(leading_edge)),
                                    "In category ({0})".format(annotations_size),
                                    "Full profile ({0})".format(signature_size),
                                ],
                            )
                            if "Leading" in display_proteins:
                                proteins = leading_edge
                            elif "category" in display_proteins:
                                proteins = annotations_[term]
                            else:
                                proteins = signature_[0]
                            o_values = [o_values[pid] for pid in proteins]
                            t_values = [t_values[pid] for pid in proteins]

                            proteins_set = set(proteins)
                            if convert_to_gene:
                                genes = [pid2gene(x) for x in proteins]
                                label = "Gene Name"
                            else:
                                label = "UniProtAC"
                            dl = pd.DataFrame(
                                {"Gene Name": genes, "UniProt AC": proteins, "Log2FC": o_values, "Z-score": t_values}
                            )

                            sort_by = col.radio(
                                "Sort proteins", ["By Z-score", "Alphabetically"]
                            )
                            if sort_by != "Alphabetically":
                                if is_up:
                                    dl = dl.sort_values("Z-score", ascending=False)
                                else:
                                    dl = dl.sort_values("Z-score", ascending=True)
                            else:
                                dl = dl.sort_values(label)
                            dl = dl.reset_index(drop=True)

                            col = cols[1]
                            col.dataframe(dl.reset_index(drop=True))

                def remove_redundancy(df, overlap_min=0.5):
                    leading_edges = []
                    for x in list(df["leading_edge"]):
                        if str(x) == "nan":
                            le = set()
                        else:
                            le = set(x.split(","))
                        leading_edges += [le]
                    idxs = []
                    done = []
                    for i, le in enumerate(leading_edges):
                        is_found = False
                        for d in done:
                            n_int = len(le.intersection(d))
                            n_min = min(len(le), len(d))
                            if n_min == 0:
                                is_found = True
                            else:
                                overlap = float(n_int) / n_min
                                if overlap >= overlap_min:
                                    is_found = True
                        if not is_found:
                            idxs += [i]
                            done += [le]
                    st.write(df.shape)
                    return df.iloc[idxs, :]

        if view == "Advanced plots":

            files_folder = os.path.join(TMP, "files")
            if os.path.exists(files_folder):
                shutil.rmtree(files_folder)
            os.mkdir(files_folder)

            files_zip = os.path.join(TMP, "files")
            if os.path.exists(files_zip+".zip"):
                os.remove(files_zip+".zip")

            top_plots_number = columns[1].number_input("Maximum number of plots", value=5, min_value=1, max_value=50)

            i = 0
            counts = 0
            for term in list(df["Term"]):
                if term not in consensus["terms"].keys():
                    continue
                if i == top_plots_number:
                    break
                fig = consensus_plot(term, consensus)
                plt.tight_layout(w_pad=0.01, h_pad=0.01)
                st.pyplot(fig, dpi=300)

                plt.savefig(os.path.join(files_folder, "file_{0}.png".format(counts)))
                plt.savefig(os.path.join(files_folder, "file_{0}.pdf".format(counts)))
                counts += 1

                i += 1
            
            if i == 0:
                st.warning("No summary plots of sufficient quality could be created...")
            else:
                shutil.make_archive(files_zip, "zip", files_folder)
                with open(files_zip+".zip", "rb") as fp:
                    btn = st.download_button(
                        label = "Download plots",
                        data = fp,
                        file_name = files_zip+".zip",
                        mime = "application/zip"
                    )


        if view == "Basic plots":

            top_plots_number = columns[1].number_input("Maximum number of plots", value=12, min_value=1, max_value=50)

            import blitzgsea as blitz
            from stylia.figure.figure import stylize

            COLORS = stylia.colors.colors.NamedColors()

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

            columns = st.columns(4)
            i = 0
            j = 0

            files_folder = os.path.join(TMP, "files")
            if os.path.exists(files_folder):
                shutil.rmtree(files_folder)
            os.mkdir(files_folder)

            files_zip = os.path.join(TMP, "files")
            if os.path.exists(files_zip+".zip"):
                os.remove(files_zip+".zip")
            
            counts = 0
            for term in list(df["Term"]):
                if j == top_plots_number:
                    break
                fig, ax = plt.subplots(1,1,figsize=(5,5))
                fig = gsea_plot(signature_, term, annotations_, result=result_)
                if i == len(columns):
                    i = 0
                col = columns[i]
                col.pyplot(fig)
                i += 1
                j += 1
                stylia.save_figure(os.path.join(files_folder, "file_{0}.png".format(counts)))
                stylia.save_figure(os.path.join(files_folder, "file_{0}.pdf".format(counts)))
                counts += 1
            
            shutil.make_archive(files_zip, "zip", files_folder)

            with open(files_zip+".zip", "rb") as fp:
                btn = st.download_button(
                    label = "Download plots",
                    data = fp,
                    file_name = files_zip+".zip",
                    mime = "application/zip"
                )


