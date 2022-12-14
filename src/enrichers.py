import os
import shutil
import csv
import collections
import pandas as pd
import numpy as np
import random
import json
from sklearn.preprocessing import PowerTransformer
from statsmodels.stats.multitest import multipletests
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(ROOT, "../tools/blitzgsea"))
import blitzgsea as blitz
from scipy.stats import fisher_exact

from util import listdir_util


class BaseEnricher(object):
    def __init__(self, signature_file, annotations_file, universe_file, min_set_size):
        self.annotations_file = annotations_file
        self.signature_file = signature_file
        self.univers_file = universe_file
        self.min_set_size = min_set_size
        self.universe = None
        self.signature = None
        self.annotations = None

    def read_annotations_sets(self):
        annotations = collections.defaultdict(set)
        with open(self.annotations_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for r in reader:
                annotations[r[1]].update([r[0]])
        return annotations

    def read_universe_standard(self):
        if self.universe is not None:
            return self.universe
        universe = []
        with open(self.univers_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for r in reader:
                universe += [r[0]]
        self.universe = universe
        return self.universe

    def read_universe_from_pulldown(self):
        if self.universe is not None:
            return self.universe
        signature_name = self.signature_file.split("/")[-1].split("_")[0]
        df = pd.read_csv(
            os.path.join(ROOT, "../data/general/cemm_primary_data.tsv"), delimiter="\t"
        )  # TODO remove hardcoding
        fids = set(df["FragID"])
        if signature_name in fids:
            exp_id = list(set(df[df["FragID"] == signature_name]["expID"]))[0]
            universe = sorted(set(df[df["expID"] == exp_id]["UniProtID"]))
        else:
            universe = sorted(set(df["UniProtID"]))
        self.universe = universe
        return self.universe

    def read_universe(self):
        if "pulldown" in self.univers_file.split("/t")[-1]:
            return self.read_universe_from_pulldown()
        else:
            return self.read_universe_standard()

    def read_signature(self):
        signature = []
        with open(self.signature_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for r in reader:
                signature += [r]
        return signature

    def get_signature(self):
        if self.signature is not None:
            return self.signature
        R = self.read_signature()
        universe = self.read_universe()
        signature = []
        for r in R:
            if r[0] not in universe:
                continue
            signature += [r]
        self.signature = signature
        return signature

    def get_annotations(self, consider_signature):
        if self.annotations is not None:
            return self.annotations
        D = self.read_annotations_sets()
        if consider_signature:
            universe = set([x[0] for x in self.get_signature()])
        else:
            universe = self.read_universe()
        annotations = {}
        for k, v in D.items():
            v = v.intersection(universe)
            if len(v) < self.min_set_size:
                continue
            annotations[k] = v
        self.annotations = annotations
        return annotations


class RankSumEnricher(BaseEnricher):
    def __init__(self, signature_file, annotations_file, universe_file, min_set_size):
        BaseEnricher.__init__(
            self,
            signature_file=signature_file,
            annotations_file=annotations_file,
            universe_file=universe_file,
            min_set_size=min_set_size,
        )
        self.get_signature()
        self.get_annotations(consider_signature=True)
        self.reformat_signature()
        self.transform_signature()
        self.reformat_annotations()

    def reformat_signature(self):
        signature = []
        for s in self.signature:
            signature += [[s[0], float(s[1])]]
        self.signature = pd.DataFrame(signature)

    def reformat_annotations(self):
        d = {}
        for k, v in self.annotations.items():
            d[k] = list(v)
        self.annotations = d

    def transform_signature(self):
        transf = PowerTransformer()
        R = []
        for r in list(self.signature[1]):
            R += [[r]]
        X = np.array(R)
        X = transf.fit_transform(X)
        self.signature[1] = X[:, 0]

    def calculate(
        self, anchors=20, permutations=2000, processes=4, signature_cache=True
    ):
        result = blitz.gsea(
            self.signature,
            self.annotations,
            anchors=anchors,
            permutations=permutations,
            min_size=self.min_set_size,
            processes=processes,
            signature_cache=signature_cache,
            verbose=True,
        )
        self.result = result
        return result

    def save(self, dirname):
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname, exist_ok=True)
        self.result.to_csv(os.path.join(dirname, "result.tsv"), sep="\t")
        meta = {"status": "done", "task": "ranksum", "min_set_size": self.min_set_size}
        with open(os.path.join(dirname, "meta.json"), "w") as f:
            json.dump(meta, f, indent=4)
        with open(os.path.join(dirname, "annotations.json"), "w") as f:
            json.dump(dict((k, sorted(v)) for k, v in self.annotations.items()), f)
        with open(os.path.join(dirname, "signature.tsv"), "w") as f:
            writer = csv.writer(f, delimiter="\t")
            for p in self.signature.values:
                writer.writerow(p)


class HyperGeometricEnricher(BaseEnricher):
    def __init__(self, signature_file, annotations_file, universe_file, min_set_size):
        BaseEnricher.__init__(
            self,
            signature_file=signature_file,
            annotations_file=annotations_file,
            universe_file=universe_file,
            min_set_size=min_set_size,
        )
        self.get_signature()
        self.get_annotations(consider_signature=False)
        self.reformat_signature()

    def reformat_signature(self):
        self.signature = [s[0] for s in self.signature]

    def calculate(self):
        print("HyperGeometric test")
        universe = set(self.universe)
        signature = set(self.signature)
        U = len(universe)
        S = len(signature)
        R = []
        commons = []
        for k, v in self.annotations.items():
            A = len(v)
            common = signature.intersection(v)
            a = len(common)
            commons += [",".join(sorted(common))]
            b = S - a
            c = A - a
            d = U - (a + b + c)
            odds, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
            R += [[k, a, A, odds, pvalue]]
        result = pd.DataFrame(
            R, columns=["Term", "overlap", "geneset_size", "odds", "pval"]
        )
        result["fdr"] = multipletests(list(result["pval"]), method="fdr_bh")[1]
        result["leading_edge"] = commons
        result = result.sort_values(by="fdr")
        result = result.reset_index(drop=True)
        result = result.set_index("Term")
        self.result = result
        return result

    def save(self, dirname):
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname, exist_ok=True)
        self.result.to_csv(os.path.join(dirname, "result.tsv"), sep="\t")
        meta = {
            "status": "done",
            "task": "hypergeometric",
            "min_set_size": self.min_set_size,
        }
        with open(os.path.join(dirname, "meta.json"), "w") as f:
            json.dump(meta, f, indent=4)
        with open(os.path.join(dirname, "annotations.json"), "w") as f:
            json.dump(dict((k, sorted(v)) for k, v in self.annotations.items()), f)
        with open(os.path.join(dirname, "universe.tsv"), "w") as f:
            writer = csv.writer(f, delimiter="\t")
            for p in sorted(self.universe):
                writer.writerow([p])
        with open(os.path.join(dirname, "signature.tsv"), "w") as f:
            writer = csv.writer(f, delimiter="\t")
            for p in sorted(self.signature):
                writer.writerow([p])


class Enricher(object):
    def __init__(self, signature_file, annotations_file, universe_file, min_set_size):
        if self._signature_is_empty(signature_file):
            self._is_empty = True
        else:
            self._is_empty = False
            if self._signature_is_set(signature_file):
                self.enricher = HyperGeometricEnricher(
                    signature_file, annotations_file, universe_file, min_set_size
                )
            else:
                self.enricher = RankSumEnricher(
                    signature_file, annotations_file, universe_file, min_set_size
                )

    def _signature_is_empty(self, signature_file):
        with open(signature_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for r in reader:
                return False
        return True

    def _signature_is_set(self, signature_file):
        with open(signature_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for r in reader:
                if len(r) == 1:
                    return True
                else:
                    return False

    def calculate(self):
        if self._is_empty:
            return None
        else:
            result = self.enricher.calculate()
            return result

    def _save_empty(self, dirname):
        meta = {"status": "empty"}
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname, exist_ok=True)
        with open(os.path.join(dirname, "meta.json"), "w") as f:
            json.dump(meta, f, indent=4)

    def save(self, dirname):
        if self._is_empty:
            self._save_empty(dirname)
        else:
            self.enricher.save(dirname)


class Harmonizer(object):
    def __init__(self, results_folder):
        self.results_folder = os.path.abspath(results_folder)
        self.df = None

    def is_ranksum(self):
        task_name = self.results_folder.rstrip("/").split("/")[-1]
        for x in ["_bin_", "_top_", "_bottom_"]:
            if x in task_name:
                return False
        return True

    def harmonize(self):
        path = os.path.join(self.results_folder, "result.tsv")
        if not os.path.exists(path):
            return
        df = pd.read_csv(path, sep="\t")
        if self.is_ranksum():
            overlap = []
            for x in list(df["leading_edge"]):
                x = str(x)
                if x == "nan":
                    overlap += [0]
                else:
                    overlap += [len(x.split(","))]
            df["overlap"] = overlap
            self.df = df[["Term", "overlap", "geneset_size", "nes", "pval"]]
            self.df.columns = ["term", "overlap", "setsize", "score", "pval"]
        else:
            self.df = df[["Term", "overlap", "geneset_size", "odds", "pval"]]
            self.df.columns = ["term", "overlap", "setsize", "score", "pval"]

    def save(self):
        if self.df is None:
            return
        self.df.to_csv(
            os.path.join(self.results_folder, "harmon.tsv"), sep="\t", index=False
        )


class PreSummarizer(object):
    def __init__(self, results_folder, signatures_folder, min_set_size):
        self.min_set_size = min_set_size
        results_folder = os.path.abspath(results_folder)
        signatures_folder = os.path.abspath(signatures_folder)
        min_set_size = min_set_size
        if not os.path.exists(results_folder):
            self._is_ready = False
        else:
            self._is_ready = True
        self.results_folder = os.path.abspath(results_folder).rstrip("/")
        self._harmonize_if_not_exists()
        self.signatures_folder = os.path.abspath(signatures_folder).rstrip("/")
        self.base_name = self.results_folder.split("/")[-2]
        self.profile_type = self.results_folder.split("/")[-1].replace(
            self.base_name, ""
        )
        self.signature_group = self.results_folder.split("/")[-3]
        self.annotations = self.results_folder.split("/")[-4]
        self.universe = self.results_folder.split("/")[-5]
        self.entity_type = self.results_folder.split("/")[-6]
        self.bases_folder = "/".join(self.signatures_folder.split("/")[:-1])
        self.results_root = "/".join(self.results_folder.split("/")[:-6])
        self.signatures_root = "/".join(self.signatures_folder.split("/")[:-3])
        self.data_root = "/".join(self.signatures_root.split("/")[:-1])
        self.annotations_root = os.path.join(self.data_root, "annotations")
        self.annotations_file = os.path.join(
            self.annotations_root, self.entity_type, self.annotations + ".tsv"
        )
        self.universes_root = os.path.join(self.data_root, "universes")
        self.universe_file = os.path.join(
            self.universes_root, self.entity_type, self.universe + ".tsv"
        )

    def _harmonize_if_not_exists(self):
        if not os.path.exists(os.path.join(self.results_folder, "harmon.tsv")):
            h = Harmonizer(self.results_folder)
            h.harmonize()
            h.save()

    def iterate_bases(self):
        for base in listdir_util(self.bases_folder):
            if base == self.base_name:
                continue
            yield base

    def get_result_file(self, base):
        path = os.path.join(
            self.results_root,
            self.entity_type,
            self.universe,
            self.annotations,
            self.signature_group,
            base,
            "{0}{1}".format(base, self.profile_type),
            "result.tsv",
        )
        return path

    def get_harmon_file(self, base):
        path = os.path.join(
            self.results_root,
            self.entity_type,
            self.universe,
            self.annotations,
            self.signature_group,
            base,
            "{0}{1}".format(base, self.profile_type),
            "harmon.tsv",
        )
        return path

    def get_signature_file(self, base):
        path = os.path.join(
            self.signatures_root,
            self.entity_type,
            self.signature_group,
            base,
            "{0}{1}.tsv".format(base, self.profile_type),
        )
        return path

    def get_available_results_files(self):
        results_files = []
        for base in self.iterate_bases():
            results_file = self.get_result_file(base)
            if os.path.isfile(results_file):
                results_files += [results_file]
        return results_files

    def get_unavailable_results_folders(self):
        results_folders = []
        signatures_files = []
        for base in self.iterate_bases():
            results_file = self.get_result_file(base)
            if not os.path.isfile(results_file):
                results_folders += ["/".join(results_file.split("/")[:-1])]
                signatures_files += [self.get_signature_file(base)]
        data = [(x, y) for x, y in zip(results_folders, signatures_files)]
        return data

    def get_summary_file(self):
        path = os.path.join(
            self.results_root,
            self.entity_type,
            self.universe,
            self.annotations,
            self.signature_group,
            "_summary_",
            self.profile_type[1:],
            "summary.tsv",
        )
        return path

    def is_summary_available(self):
        path = self.get_summary_file()
        if os.path.exists(path):
            return True
        else:
            return False

    def run(self, overwrite=False, minimum_number_of_cases=30):
        if not overwrite:
            if self.is_summary_available():
                return
        done_results = self.get_available_results_files()
        todo_results = []
        if len(done_results) < minimum_number_of_cases:
            todo_results = self.get_unavailable_results_folders()
            if todo_results:
                n = min(len(todo_results), minimum_number_of_cases - len(done_results))
                todo_results = random.sample(todo_results, n)
        if todo_results:
            for r in todo_results:
                results_folder = r[0]
                signature_file = r[1]
                enricher = Enricher(
                    signature_file=signature_file,
                    annotations_file=self.annotations_file,
                    universe_file=self.universe_file,
                    min_set_size=self.min_set_size,
                )
                enricher.calculate()
                enricher.save(results_folder)
                harmonizer = Harmonizer(results_folder=results_folder)
                harmonizer.harmonize()
                harmonizer.save()


class Summarizer(object):
    def __init__(self, results_folder, overwrite=False):
        if not os.path.exists(results_folder):
            self._is_ready = False
        else:
            self._is_ready = True
        self.results_folder = os.path.abspath(results_folder).rstrip("/")
        self.base_name = self.results_folder.split("/")[-2]
        self.profile_type = self.results_folder.split("/")[-1].replace(
            self.base_name, ""
        )
        self.signature_group = self.results_folder.split("/")[-3]
        self.annotations = self.results_folder.split("/")[-4]
        self.universe = self.results_folder.split("/")[-5]
        self.entity_type = self.results_folder.split("/")[-6]
        self.results_root = "/".join(self.results_folder.split("/")[:-6])
        self._is_done = self.file_exists()
        self._overwrite = overwrite
        if not self._is_done:
            self._todo = True
        else:
            if self._overwrite:
                self._todo = True
            else:
                self._todo = False

    def all_results_folders(self):
        upper_results_folder = "/".join(self.results_folder.split("/")[:-2])
        arfs = []
        for d in listdir_util(upper_results_folder):
            if d.startswith("_"):
                continue
            arf = os.path.join(upper_results_folder, d, d + self.profile_type)
            meta = os.path.join(arf, "meta.json")
            if os.path.exists(meta):
                with open(meta, "r") as f:
                    meta_data = json.load(f)
                if meta_data["status"] == "done":
                    arfs += [arf]
        return arfs

    def harmon_as_dict(self, df):
        rd = collections.OrderedDict()
        for r in df.values:
            rd[r[0]] = r[1:]
        return rd

    def summarize(self):
        if not self._is_ready:
            return
        if not self._todo:
            return
        all_overlap = collections.defaultdict(list)
        all_setsize = collections.defaultdict(list)
        all_score = collections.defaultdict(list)
        all_pvalue = collections.defaultdict(list)
        all_results_sample_size = 0
        keys = set()
        for arf in self.all_results_folders():
            harmon = pd.read_csv(os.path.join(arf, "harmon.tsv"), delimiter="\t")
            rd = self.harmon_as_dict(harmon)
            for k, v in rd.items():
                all_overlap[k] += [int(np.clip(v[0], 0, 1e8))]
                all_setsize[k] += [int(np.clip(v[1], 0, 1e8))]
                all_score[k] += [float(np.clip(v[2], -10, 10))]
                all_pvalue[k] += [float(np.clip(v[3], 0, 1))]
                keys.update([k])
            all_results_sample_size += 1
        R = []
        keys = sorted(keys)
        for k in keys:
            means = []
            stds = []
            means += [np.mean(all_overlap[k])]
            stds += [np.std(all_overlap[k])]
            means += [np.mean(all_setsize[k])]
            stds += [np.std(all_setsize[k])]
            means += [np.mean(all_score[k])]
            stds += [np.std(all_score[k])]
            means += [np.mean(all_pvalue[k])]
            stds += [np.std(all_pvalue[k])]
            R += [[k] + means + stds]
        df = pd.DataFrame(
            R,
            columns=[
                "term",
                "overlap_mean",
                "setsize_mean",
                "score_mean",
                "pval_mean",
                "overlap_std",
                "setsize_std",
                "score_std",
                "pval_std",
            ],
        )
        self.summary = df
        self.all_results_sample_size = all_results_sample_size

    def get_summary_folder(self):
        path = os.path.join(
            self.results_root,
            self.entity_type,
            self.universe,
            self.annotations,
            self.signature_group,
            "_summary_",
            self.profile_type[1:],
        )
        return path

    def get_summary_file(self):
        return os.path.join(self.get_summary_folder(), "summary.tsv")

    def get_summary_meta_file(self):
        return os.path.join(self.get_summary_folder(), "summary_meta.json")

    def is_summary_available(self):
        if os.path.exists(self.get_summary_file()):
            return True
        else:
            return False

    def file_exists(self):
        if not self._is_ready:
            return False
        if os.path.exists(self.get_summary_file()):
            return True
        else:
            return False

    def save(self):
        if not self._is_ready:
            return
        if not self._todo:
            return
        output_folder = self.get_summary_folder()
        os.makedirs(output_folder, exist_ok=True)
        self.summary.to_csv(
            self.get_summary_file(), sep="\t", index=False
        )
        with open(self.get_summary_meta_file(), "w") as f:
            meta = {"sample_size": self.all_results_sample_size}
            json.dump(meta, f, indent=4)


class GroupSummarizer(object):

    def __init__(self, data_path):
        fid2group = {}
        df = pd.read_csv(os.path.join(data_path, "general", "cemm_primary_data.tsv"), delimiter="\t")
        for r in df[["FragID", "expID"]].values:
            fid2group[r[0]] = r[1]
        group2fid = collections.defaultdict(list)
        for k,v in fid2group.items():
            group2fid[v] += [k]
        self.fid2group = fid2group
        self.group2fid = group2fid
        self.fid = None

    def set(self, results_folder):
        self.results_folder = os.path.abspath(results_folder).rstrip("/")
        self.base_name = self.results_folder.split("/")[-2]
        self.profile_type = self.results_folder.split("/")[-1].replace(
            self.base_name, ""
        )
        self.signature_group = self.results_folder.split("/")[-3]
        self.annotations = self.results_folder.split("/")[-4]
        self.universe = self.results_folder.split("/")[-5]
        self.entity_type = self.results_folder.split("/")[-6]
        self.results_root = "/".join(self.results_folder.split("/")[:-6])
        # it only works with fragment ids at the moment
        # TODO
        fid = self.base_name
        self.fid = fid
        self.other_fids = []
        if fid in self.fid2group:
            gid = self.fid2group[fid]
            for k in self.group2fid[gid]:
                if k == fid:
                    continue
                self.other_fids += [k]
        
    def summarize(self):
        assert self.fid is not None
        df = pd.read_csv(os.path.join(self.results_folder, "harmon.tsv"), delimiter="\t")
        df = df[df["pval"] <= 0.05] # TODO
        terms = list(df["term"]) # TODO
        scores = list(df["score"])
        pvals = list(df["pval"])
        upper_results_folder = "/".join(self.results_folder.split("/")[:-2])
        if len(terms) > 0:
            pairs = {}
            for ofid in self.other_fids:
                df = pd.read_csv(os.path.join(upper_results_folder, ofid, ofid + self.profile_type, "harmon.tsv"), delimiter="\t")
                for term in terms:
                    df_ = df[df["term"] == term]
                    if df_.shape[0] == 0:
                        sc_ = None
                        pv_ = None
                    else:
                        sc_ = list(df_["score"])[0]
                        pv_ = list(df_["pval"])[0]
                    pairs[(ofid, term)] = (sc_, pv_)
        else:
            pairs = {}
        data = {}
        for t, s, p in zip(terms, scores, pvals):
            data_ = {}
            data_["this"] = [s, p]
            data_["other"] = {}
            for k,v in pairs.items():
                if k[1] == t:
                    data_["other"][k[0]] = [v[0], v[1]]
            data[t] = data_
        self.data = data
        self.fid = None

    def save(self):
        with open(os.path.join(self.results_folder, "group_data.json"), "w") as f:
            json.dump(self.data, f)

    def clean(self):
        path = os.path.join(self.results_folder, "group_data.json")
        if os.path.exists(path):
            os.remove(path)
