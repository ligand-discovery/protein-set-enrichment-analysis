import os
import sys
import csv
import json
import collections
import pandas as pd
import numpy as np
import joblib

ROOT = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(ROOT, "../src/"))
from enrichers import Summarizer, Harmonizer, PreSummarizer

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

def universe_selector():
    preselected="HEK293T Core"
    universe = preselected
    universe_subfolder = universe_dict[universe]
    return universe, universe_subfolder

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


def annotations_iterator():
    annotation_types = [
        "Sequence",
        "Functions",
        "Processes and pathways",
        "Localization",
        "Drugs and Diseases",
    ]
    R = []
    for annotation_type in annotation_types:
        annotations = annotation_type_dict[annotation_type]
        for annotation in annotations:
            annotation_subfolder = annotation_dict[annotation]
            R += [(annotation, annotation_subfolder, annotation_type, annotations)]
    return R
        

iter_count = 0
for profile in fids:
    profile_subfolder = profile
    all_cases = fids

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

    universe, universe_subfolder = universe_selector()

    for annotation, annotation_subfolder, annotation_type, annotations in annotations_iterator():
        print(iter_count, profile_subfolder, annotation_subfolder)
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

        consensus = ec.run_consensus()

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

        df_view = df_merge[["term", "overlap", "setsize", "score", "pval", "score_mean"]]
        df_view = df_view.rename(columns={
            "term": "Term",
            "overlap": "Edge size",
            "setsize": "Set size",
            "score": "Score",
            "pval": "P-value",
            "score_mean": "Expected score"
        })

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

        with open(os.path.join(results_path, "meta.json"), "r") as f:
            meta = json.load(f)

        if meta["status"] == "empty":
            continue

        if meta["status"] == "done":

            type_of_task = meta["task"]

            # read main results
            result = pd.read_csv(
                os.path.join(results_path, "result.tsv"), delimiter="\t"
            )

            result_ = result.copy()
            result_ = result_.set_index("Term")

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

            metrics = {
                "type_of_task": type_of_task,
                "signature_size": signature_size,
                "annotations_size": annotations_size,
            }

            if type_of_task == "hypergeometric":
                universe_size = pd.read_csv(
                    os.path.join(results_path, "universe.tsv"),
                    delimiter="\t",
                    header=None,
                ).shape[0]
                metrics["title"] = "Background size"
                metrics["value"] = universe_size
  
            else:
                average_annotations_per_term = np.mean([len(v) for k,v in annotations_.items()])
                metrics["title"] = "Average number of proteins per term"
                metrics["value"] = average_annotations_per_term

            folder_path = os.path.join(CACHE, "detailed", profile_subfolder, annotation_subfolder)
            if not os.path.exists(folder_path):
                os.makedirs(folder_path, exist_ok=True)
            
            with open(os.path.join(folder_path, "metrics.json"), "w") as f:
                json.dump(metrics, f, indent=4)

            all_annotations = sorted(annotations_.keys())

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

            MIN_EDGE_SIZE = 1
            MAX_EDGE_SIZE = 1000

            result = result[result["nes"] > 0]
            result = result[result["leading_edge_size"] >= MIN_EDGE_SIZE]
            result = result[result["leading_edge_size"] <= MAX_EDGE_SIZE]
            result = result.reset_index(drop=True)
            
            prot2idx = collections.defaultdict(list)
            for i, r in enumerate(list(result["leading_edge"])):
                if str(r) == "nan":
                    continue
                for x in r.split(","):
                    prot2idx[pid2gene(x)] += [i]

            selections = {
                "all_annotations": all_annotations,
                "available_proteins": available_proteins,
                "leading_proteins": sorted(prot2idx.keys())
            }

            with open(os.path.join(folder_path, "selections.json"), "w") as f:
                json.dump(selections, f)

            paths = {
                "results_path": results_path.split("/protein-set-enrichment-analysis/")[1]
            }
            with open(os.path.join(folder_path, "paths.json"), "w") as f:
                json.dump(paths, f)

            result.to_csv(os.path.join(folder_path, "result.tsv"), sep="\t", index=False)
            df_merge.to_csv(os.path.join(folder_path, "df_merge.tsv"), sep="\t", index=False)
            
            joblib.dump(consensus, os.path.join(folder_path, "consensus.joblib"))

        iter_count += 1

