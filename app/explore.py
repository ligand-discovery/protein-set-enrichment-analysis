import os
import csv
import streamlit as st
import pandas as pd
import json

st.set_page_config(layout="wide")

st.title("Quick fragment search")

ROOT = os.path.dirname(os.path.abspath(__file__))

RESULTS_FOLDER = os.path.join(os.path.abspath(os.path.join(ROOT, "..", "results")))
RESULTS_FILE = os.path.abspath(os.path.join(ROOT, "..", "results", "all_results.tsv"))

SPEC_FRAG_CUTOFF = 50

columns = ["term", "overlap", "setsize", "score", "pval", "category", "fid"]

def get_frag_counts():
    df = pd.read_csv(os.path.join(ROOT, "../data/general/cemm_primary_hit_data.tsv"), sep="\t")
    df = pd.DataFrame(df.value_counts("FragID"), columns=["counts"])
    df = df.reset_index()
    df = df.rename(columns={"FragID": "fid"})
    fids = set(pd.read_csv(os.path.join(ROOT, "../data/general/cemm_primary_data.tsv"), sep="\t")["FragID"])
    fids = fids.difference(df["fid"])
    dm = pd.DataFrame({"fid": sorted(fids), "counts": [0]*len(fids)})
    df = pd.concat([df, dm], axis=0).reset_index(drop=True)
    return df

spec_frag = get_frag_counts()
to_remove_fids = set(list(spec_frag[spec_frag["counts"] >= SPEC_FRAG_CUTOFF]["fid"]))

def file_filter(p_max, file_name):
    R = []
    with open(RESULTS_FILE, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for r in reader:
            o = float(r[1])
            m = float(r[2])
            s = float(r[3])
            p = float(r[4])
            if m >= 10 and s > 0 and o > 3 and p < p_max:
                R += [r]
    with open(file_name, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(columns)
        for r in R:
            writer.writerow(r)

mc_file = os.path.join(RESULTS_FOLDER, "mc_results.tsv")
if not os.path.exists(mc_file):
    file_filter(0.01, mc_file)

lc_file = os.path.join(RESULTS_FOLDER, "lc_results.tsv")
if not os.path.exists(lc_file):
    file_filter(0.05, lc_file)

hc_file = os.path.join(RESULTS_FOLDER, "hc_file.tsv")
if not os.path.exists(hc_file):
    file_filter(0.001, hc_file)


def terms_json(file_name):
    terms = []
    cats = []
    fids = []
    with open(RESULTS_FILE, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for r in reader:
            terms += [r[0]]
            cats += [r[5]]
            fids += [r[6]]
    with open(file_name, "w") as f:
        data = {
            "terms": sorted(set(terms)),
            "categories": sorted(set(cats)),
            "fids": sorted(set(fids))
        }
        json.dump(data, f)

mc_json = os.path.join(RESULTS_FOLDER, "mc_terms.json")
if not os.path.exists(mc_json):
    terms_json(mc_json)

lc_json = os.path.join(RESULTS_FOLDER, "lc_terms.json")
if not os.path.exists(lc_json):
    terms_json(lc_json)

hc_json = os.path.join(RESULTS_FOLDER, "hc_terms.json")
if not os.path.exists(hc_json):
    terms_json(hc_json)


cols = st.columns(3)
confidence = cols[0].radio(label="Confidence (enrichment)", options=["High", "Medium", "Low"], index=1, horizontal=True)
if confidence == "Medium":
    data_file = mc_file
    json_file = mc_json
if confidence == "High":
    data_file = hc_file
    json_file = hc_json
if confidence == "Low":
    data_file = lc_file
    json_file = lc_json

with open(json_file, "r") as f:
    terms_data = json.load(f)

prom = cols[1].radio(label="Promiscuous fragments", options=["Yes", "No"], index=0, horizontal=True)
if prom == "Yes":
    remove_prom = False
else:
    remove_prom = True

category_text = cols[2].multiselect(label="Search for a category...", options=terms_data["categories"])
annotation_text = st.text_input(label="Search for a term...")

if category_text:
    do_category = True
    category_text = set(category_text)
else:
    do_category = False

if annotation_text:
    with open(data_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        R = []
        for r in reader:
            if remove_prom:
                if r[6] in to_remove_fids:
                    continue
            if do_category:
                if r[5] not in category_text:
                    continue
            if annotation_text in r[0].lower():
                R += [[r[0], int(r[1]), int(r[2]), round(float(r[3]), 3), round(float(r[4]), 5), r[5], r[6]]]

    df = pd.DataFrame(R, columns=["Term", "Edge size", "Set size", "Score", "P-value", "Category subtype", "Fragment"])

    df = df.sort_values("Score", ascending=False).reset_index(drop=True)
    df["rank"] = [i+1 for i in range(df.shape[0])]
    df = df.set_index("rank")

    st.dataframe(df, width=1000, height=1000)

