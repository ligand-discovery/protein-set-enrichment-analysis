import os
import sys
import pandas as pd
import json
import csv
import joblib
import shutil
import h5py
import warnings
import collections
import stylia
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from sklearn.preprocessing import QuantileTransformer
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.preprocessing import MinMaxScaler
from sklearn.neighbors import NearestNeighbors
from umap import UMAP
from griddify import Cloud2Grid
warnings.filterwarnings("ignore")

ROOT = os.path.abspath(os.path.dirname(__file__))
PATH = os.path.abspath(os.path.join(ROOT, "../results/proteins/"))
DATA = os.path.abspath(os.path.join(ROOT, "../data"))
CACHE = os.path.abspath(os.path.join(ROOT, "../cache"))

sys.path.append(os.path.join(ROOT, "../src/"))
from proteome_meta import universe_dict
from proteome_meta import annotation_dict

annotation_dict_inv = dict((v,k) for k,v in annotation_dict.items())

top_plots_number = 10

def iterate_all_cache_folders():
    R = []
    for profile_subfolder in os.listdir(os.path.join(CACHE, "detailed")):
        for annotations_subfolder in os.listdir(os.path.join(CACHE, "detailed", profile_subfolder)):
            R += [(profile_subfolder, annotations_subfolder)]

    R = sorted(R, key=lambda x: x[0])
    return R


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

profile_type = "Fragment"
profile_type_subfolder = profile_type.lower()


for profile_subfolder, annotations_subfolder in iterate_all_cache_folders():

    profile = profile_subfolder

    profile_number = int(profile[1:])
    if profile_number < 127:
        continue

    annotation = annotation_dict_inv[annotations_subfolder]
    
    universe = "HEK293T Core"
    universe_subfolder = universe_dict[universe]

    cache_folder = os.path.join(CACHE, "detailed", profile_subfolder, annotations_subfolder)

    path_ = os.path.join(
        "results",
        "proteins",
        universe_subfolder,
        annotations_subfolder,
        profile_type_subfolder,
        profile_subfolder,
    )

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
        X = red.fit_transform(np.asarray(X), y=list(np.clip(df["score"], -5, 5)))
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

        profile_type_subfolder = "fragment"
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

    consensus = joblib.load(os.path.join(cache_folder, "consensus.joblib"))

    i = 0
    j = 0
    
    MIN_EDGE_SIZE = 5
    MAX_EDGE_SIZE = 5000

    result = pd.read_csv(os.path.join(cache_folder, "result.tsv"), delimiter="\t")
    result = result[result["leading_edge_size"] >= MIN_EDGE_SIZE]
    result = result[result["leading_edge_size"] <= MAX_EDGE_SIZE]
    df = result.copy()
    df = df.sort_values("nes", ascending=False).reset_index(drop=True)
    df = df.head(top_plots_number)

    print(profile_subfolder, annotations_subfolder, df.shape[0])

    if os.path.exists(os.path.join(cache_folder, "advanced")):
        shutil.rmtree(os.path.join(cache_folder, "advanced"))
    os.makedirs(os.path.join(cache_folder, "advanced"))

    i = 0
    counts = 0
    idx2term = {}
    for term in list(df["Term"]):
        if term not in term_locs.keys():
            continue
        if term not in consensus["terms"].keys():
            print(term, "not in consensus")
            continue
        if i == top_plots_number:
            break
        fig = consensus_plot(term, consensus)
        plt.tight_layout(w_pad=0.01, h_pad=0.01)
        
        figure_path = os.path.join(cache_folder, "advanced", "plot_{0}.png".format(counts))
        plt.savefig(figure_path, dpi=300)

        plt.close(fig)
        
        idx2term[counts] = term

        counts += 1

        i += 1

    with open(os.path.join(cache_folder, "advanced", "idx2term.json"), "w") as f:
        json.dump(idx2term, f)