# Protein Set Enrichment Analysis for Ligand Discovery
Perform and explore enrichment analyses based on Ligand Discovery primary screening data

## Usage

The repository contains a large amount of orthogonal data collected from the public domain. The data corresponds to protein annotations of multiple types, including (but not limited to):

### Sequence
 - Domain, family, philogeny groups
 - Active site, binding site
### Functions
 - Protein class
 - Molecular function
### Processes and pathways
 - Biological processes
 - Pathways
 - Complexes
### Localization
 - Subcellular localization
 - Cellular component
### Drugs and diseases
 - Drug target classes
 - Target druggability
 - Disease category

All annotations are available [here](https://ligand-discovery.s3.eu-central-1.amazonaws.com/protein-set-enrichment-analysis/protein_annotations.zip).

## Primitive plots
For each fragment, we ranked proteins by their Log2FC (z-normalized) and performed a ranksum (GSEA-like) enrichment test across all annotations. The figure below shows some annotations found to be enriched for fragment C001.

<img width="1067" alt="primitive-plots" src="https://user-images.githubusercontent.com/19725330/233401945-a6937ed6-1e8b-43e2-9fdb-3c486e4bc47c.png">
Conventional Ranksum enrichment analysis.

In addition, we performed **hypergeometric tests**, based on binarized data, as well as top-25, 50, 100, 250 and 500. The protein universe used was the basal proteome of HEK293T.
We designed a primitive version of the Streamlit App to navigate the enormous amount of enrichment results. Two limitations became apparent:
* A panel displaying a large number of top enrichment results was necessary in order to extract biological insights.
* Promiscuity of proteins propagates to promiscuity of enrichment results, resulting in frequently occurring annotation terms.

## Advanced plots
To address the above limitations, we provide the following two plot types.

### Leaderboard
The leaderboard below corresponds to fragment C170. Vacuolar proteins (a GO Cellular componet) are enriched for this fragment. In the leading edge of this enrichment result, we find TMEM59, TPP1, etc. The normalized enrichment score (NES) is 5.95, and the P-value is 2.9e-09. In red, we see high Log2FC values, for the vacuolar proteins, and in blue lower Log2FCs. The dot at the right is colored by category (in this case, localization).

<img width="1067" alt="leaderboard" src="https://user-images.githubusercontent.com/19725330/233402017-cd36d89d-0b5a-4cd5-9bc9-b14e0109a67a.png">
Leaderboard plot. The leaderboard can have an arbitrary length (10, 50, 100...).

### In-depth plots
Here we focus on one particular annotation (Vacuolar Lumen) and fragment (C175).

<img width="1091" alt="indepth" src="https://user-images.githubusercontent.com/19725330/233402523-6e22166f-e303-4e01-924b-26f0722e7c78.png">
In depth-plot. (Left) The promiscuity plot highlights proteins in the annotaiton (coloured). Filled circles correspond to the leading edge. Color denotes promiscuity (blue) or specificity (red). (Center) On top, ranksum plot, including circles denoting the result of a hypergeometric test at top-25, top-50, top-250 and top-500. Empty dots denote non-significant result (P-value > 0.05). In the bottom, top-10 proteins in the leading edge, colored and located by promiscuity. (Right) In the upper-left panel, the expected normalized enrichment score (NES) of this annotation across other fragments is shown (mean and standard deviation), along with fragments of the same pull down (in black). In the upper-right panel, a griddified projection of annotations is shown (coloured by promiscuity), in order to geolocate the annotation with respect to the rest of annotations. In the lower-left panel, the number of proteins in the leading edge at different degrees of promiscuity is shown. In the lower-right panel, proteins are projected (and griddified) by sequence similarity, and the leading edge proteins are highlighted (coloured by promiscuity).

The current Streamlit App capitalizes on these two display items to provide informative navigation of the enrichment results. The following is a mockup of the Streamlit App:

<img width="946" alt="mockup" src="https://user-images.githubusercontent.com/19725330/233402593-be189f9c-7365-4840-ad92-b87ad22ca6d8.png">
Mockup of the Streamlit Protein Set Enrichment Analysis App. The two main pages are highlighted. On the left, we sketch the leaderboard page, focused on a given fragment. On the right, we sketch the focus page, specific to a fragment-category pair.

## Installation

We recommend using Conda. Make sure a C++ compiler is installed:

```bash
conda install -c conda-forge cxx-compiler
```

Install the necessary dependencies:
```bash
pip install -r requirements.txt
```

### Download the necessary data

If you are interested in running a local instance of this app, please reach out to [miquel@ersilia.io](mailto:miquel@ersilia.io). Necessary data is >50GB in size.

## Run

You can run the app with:
```bash
streamlit run app/proteome.py
```

You can also run an app for quick exploration:
```bash
streamlit run app/explore.py
```

## About

This project was performed at [Georg Winter Lab](https://www.winter-lab.com/), based at [CeMM](https://cemm.at/), Vienna.
