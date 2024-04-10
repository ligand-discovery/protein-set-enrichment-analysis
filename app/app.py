# regular imports
import os
import sys
import csv
import collections
import pandas as pd
import streamlit as st
import json

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
from proteome_meta import task_suf
from proteome_meta import annotation_type_dict
from proteome_meta import annotation_dict
from proteome_meta import universe_dict

# set page layout
st.set_page_config(layout="wide", page_title="Ligand Discovery Protein Set Enrichment Analysis")

# path to results and original data
PATH = os.path.abspath(os.path.join(ROOT, "../results/proteins/"))
DATA = os.path.abspath(os.path.join(ROOT, "../data"))
CACHE = os.path.abspath(os.path.join(ROOT, "../cache"))

# generic inputs

# protein id to gene name
@st.cache_data
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

# side bar

st.sidebar.title("Ligand Discovery Proteome Set Enrichment Analysis")

# signatures (aka profiles)
st.sidebar.header("Select a fragment")

profile_type = PROFILE_TYPE
profile_type_subfolder = profile_type.lower()

@st.cache_data
def get_sorted_fids():
    fids = []
    for fid in listdir_util(os.path.join(DATA, "signatures", "proteins", "fragment")):
        fids += [fid]
    fids = sorted(fids)
    return fids

fids = get_sorted_fids()
profile = st.sidebar.selectbox("Fragment identifier", options=fids)
profile_subfolder = profile
all_cases = fids
draw_fragment = True

st.sidebar.header("Choose a type of analysis")

type_of_analysis = st.sidebar.radio(
    "Type of analysis", options=["Overview", "Detailed"]
)

@st.cache_data
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

# OVERVIEW TYPE OF ANALSYS

if type_of_analysis == "Overview":

    universe, universe_subfolder = universe_selector()
    universe_file = os.path.abspath(
        os.path.join(DATA, "universes", "proteins", universe_subfolder + ".tsv")
    )

    @st.cache_data
    def get_annotation_keys_and_labels():
        annotation_type_labels = []
        annotation_keys = []
        for at, v in annotation_type_dict.items():
            for x in v:
                an = annotation_dict[x]
                annotation_type_labels += [(at, x, an)]
                annotation_keys += ["{0}---{1}".format(at, an)]
        return annotation_keys, annotation_type_labels
    
    annotation_keys, annotation_type_labels = get_annotation_keys_and_labels()

    @st.cache_data
    def get_annotation_files():
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
        return tasks, annotation_files, signature_files, results_folders
    
    tasks, annotation_files, signature_files, results_folders = get_annotation_files()

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

    @st.cache_data
    def get_annotation_items_and_all_proteins():
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
        return annotation_items, item_annotations, all_proteins
    
    annotation_items, item_annotations, all_proteins = get_annotation_items_and_all_proteins()

    @st.cache_data
    def get_proteins_in_universe():
        proteins_in_universe = []
        with open(
            os.path.join(DATA, "universes", "proteins", universe_subfolder + ".tsv"), "r"
        ) as f:
            reader = csv.reader(f, delimiter="\t")
            for r in reader:
                proteins_in_universe += [r[0]]
        return proteins_in_universe
    
    proteins_in_universe = get_proteins_in_universe()

    annotation_type_inv = {}
    for k, v in annotation_type_dict.items():
        for x in v:
            if x in annotation_dict:
                annotation_type_inv[x] = k

    df = pd.read_csv(os.path.join(CACHE, "overview", "{0}.tsv".format(profile)), sep="\t")

    if view == "Table":

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

        st.dataframe(df_view.reset_index(drop=True), height=2000)
    
    else:

        st.image(os.path.join(CACHE, "overview", "{0}.png".format(profile)))


## DETAILED TYPE OF ANALYSIS

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
    

    annotation, annotation_subfolder, annotation_type, annotations = (
        annotations_selector()
    )
    universe, universe_subfolder = universe_selector()

    cache_folder = os.path.join(CACHE, "detailed", profile_subfolder, annotation_subfolder)

    path_ = os.path.join(
        "results",
        "proteins",
        universe_subfolder,
        annotation_subfolder,
        profile_type_subfolder,
        profile_subfolder,
    )

    st.header("Fragment: {0} & Category: {2} ({1})".format(profile_subfolder, annotation_type, annotation))

    # read metrics

    with open(os.path.join(cache_folder, "metrics.json"), "r") as f:
        metrics = json.load(f)

    metric_cols = st.columns(3)
    metric_cols[0].metric(
        "{0} profile: {1}".format(profile_type, profile),
        value="{0} proteins".format(metrics["signature_size"]),
    )
    metric_cols[1].metric(
        "{0}: {1}".format(annotation_type, annotation),
        value="{0} categories".format(metrics["annotations_size"]),
    )
    metric_cols[2].metric(metrics["title"], value=round(metrics["value"], 2))

    columns = st.columns(6)
    view = columns[0].radio("View", options=["Tables", "Basic plots", "Advanced plots"])

    if view == "Tables":

        p_value_cutoff = columns[2].number_input("P-value cutoff", value=0.05, min_value=0., max_value=1., format="%.3f")
        min_edge_size = columns[3].number_input("Minimum leading edge size", value=5, min_value=0, max_value=10000)
        max_edge_size = columns[4].number_input("Maximum leading edge size", value=5000, min_value=1, max_value=10000)
        protein_label = "Gene Name"
        if protein_label == "Gene Name":
            convert_to_gene = True
        else:
            convert_to_gene = False

        available_selections = json.load(open(os.path.join(cache_folder, "selections.json"), "r"))

        all_annotations = available_selections["all_annotations"]
        available_proteins = available_selections["available_proteins"]

        select_columns = st.columns(3)
        selected_annotations = select_columns[2].multiselect(
            "Select annotation categories", options=available_proteins
        )

        selected_proteins = select_columns[0].multiselect(
            "Filter by proteins found in at least one annotation term ({0})".format(
                len(available_proteins)
            ),
            options=available_proteins,
        )

        task = "Log2FC"
        task_filename = task_filenames[task]

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

        with open(os.path.join(results_path, "annotations.json"), "r") as f:
            annotations_ = json.load(f)

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

        result = pd.read_csv(os.path.join(cache_folder, "result.tsv"), sep="\t")
        result = result[result["leading_edge_size"] >= min_edge_size]
        result = result[result["leading_edge_size"] <= max_edge_size]
        result = result.reset_index(drop=True)

        leading_proteins = available_selections["leading_proteins"]

        selected_leading_proteins = select_columns[1].multiselect(
            "Filter by proteins found in at least one leading edge",
            options = leading_proteins)

        if selected_leading_proteins:

            prot2idx = collections.defaultdict(list)
            for i, r in enumerate(list(result["leading_edge"])):
                if str(r) == "nan":
                    continue
                for x in r.split(","):
                    prot2idx[pid2gene(x)] += [i]

            idxs = []
            for v in selected_leading_proteins:
                for x in prot2idx[v]:
                    idxs += [x]
            idxs = sorted(set(idxs))
            result = result.iloc[idxs]

        df_merge = pd.read_csv(os.path.join(cache_folder, "df_merge.tsv"), sep="\t")

        type_of_task = metrics["type_of_task"]
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

        st.dataframe(df[[c for c in list(df.columns)[:-1] if c != "Mean score"]].reset_index(drop=True))

        term = st.selectbox("Explore term...", df["Term"])

        if term is not None:

            signature_ = pd.read_csv(
                os.path.join(results_path, "signature.tsv"), delimiter="\t", header=None
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

    if view == "Basic plots":
        top_plots_number = columns[1].number_input("Maximum number of plots", value=12, min_value=1, max_value=50)
        plot_columns = st.columns(4)

        with open(os.path.join(cache_folder, "basic", "idx2term.json"), "r") as f:
            idx2term = json.load(f)

        idxs = [i for i in range(len(idx2term))]

        i = 0
        j = 0

        for idx in idxs:

            if i == len(plot_columns):
                i = 0
            col = plot_columns[i]

            if j == top_plots_number:
                break

            col.image(os.path.join(cache_folder, "basic", "plot_{0}.png".format(idx)))
            i += 1
            j += 1


    if view == "Advanced plots":
        top_plots_number = columns[1].number_input("Maximum number of plots", value=5, min_value=1, max_value=10)

        with open(os.path.join(cache_folder, "advanced", "idx2term.json"), "r") as f:
            idx2term = json.load(f)

        idxs = [i for i in range(len(idx2term))]

        j = 0
        for idx in idxs:
            if j == top_plots_number:
                break

            st.image(os.path.join(cache_folder, "advanced", "plot_{0}.png".format(idx)))
            j += 1
