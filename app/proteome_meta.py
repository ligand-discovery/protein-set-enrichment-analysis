import collections

# Global profiles
global_profiles_dict = {"Detectability": "detectability", "Promiscuity": "promiscuity"}

# Task suffixes
task_suf = collections.OrderedDict()

task_suf["_val_log2fc"] = "Log2FC"
task_suf["_gauss"] = "Gaussianized score"

task_suf["_bin_hit"] = "Protein hits (custom cutoffs)"
task_suf["_bin_50"] = "Top 50 proteins"
task_suf["_bin_100"] = "Top 100 proteins"
task_suf["_bin_250"] = "Top 250 proteins"
task_suf["_bin_500"] = "Top 500 proteins"
task_suf["_top_50"] = "Top 50 proteins"
task_suf["_top_100"] = "Top 100 proteins"
task_suf["_top_250"] = "Top 250 proteins"
task_suf["_top_500"] = "Top 500 proteins"
task_suf["_bottom_50"] = "Bottom 50 proteins"
task_suf["_bottom_100"] = "Bottom 100 proteins"
task_suf["_bottom_250"] = "Bottom 250 proteins"
task_suf["_bottom_500"] = "Bottom 500 proteins"

# Annotation type dict
annotation_type_dict = {
    "Sequence": [
        "InterPro Superfamily",
        "InterPro Family",
        "InterPro Domain",
        "InterPro Active site",
        "InterPro Binding site",
        "InterPro Conserved site",
        "InterPro Repeat",
        "InterPro PTM"
    ],
    "Functions": [
        "Panther Protein class",
        "Protein Atlas Protein class",
        "Protein Atlas Molecular function",
        "GO Molecular function",
    ],
    "Processes and pathways": [
        "Reactome Pathways",
        "KEGG Pathways",
        "WikiPathways",
        "Protein Atlas Biological process",
        "GO Biological process",
    ],
    "Localization": [
        "Protein Atlas Subcellular",
        "GO Cellular component",
        "OpenCell Localization",
    ],
    "Drugs and Diseases": [
        "Pharos IDG Category",
        "Pharos Drug Target Ontology",
        "Protein Atlas Disease involvement",
        "Human Phenotype Ontology",
    ],
}

# Annotation dict
annotation_dict = {
    # Sequence
    "InterPro Superfamily": "interpro_homologous_superfamily",
    "InterPro Family": "interpro_family",
    "InterPro Domain": "interpro_domain",
    "InterPro Active site": "interpro_active_site",
    "InterPro Binding site": "interpro_binding_site",
    "InterPro Conserved site": "interpro_conserved_site",
    "InterPro Repeat": "interpro_repeat",
    "InterPro PTM": "interpro_ptm",
    # Functions
    "Panther Protein class": "panther_protein_class",
    "Protein Atlas Protein class": "protein_atlas_protein_class",
    "Protein Atlas Molecular function": "protein_atlas_molecular_function",
    "GO Molecular function": "msigdb_gomf",
    # Pathways
    "Reactome Pathways": "msigdb_reactome",
    "KEGG Pathways": "msigdb_kegg",
    "WikiPathways": "msigdb_wp",
    "Protein Atlas Biological process": "protein_atlas_biological_process",
    "GO Biological process": "msigdb_gobp",
    # Localization
    "Protein Atlas Subcellular": "protein_atlas_subcellular_location_all", # _main?
    "GO Cellular component": "msigdb_gocc",
    "OpenCell Localization": "opencell_localization",
    # Drugs and Diseases
    "Pharos IDG Category": "pharos_protein_category",
    "Pharos Drug Target Ontology": "pharos_dto",
    "Protein Atlas Disease involvement": "protein_atlas_disease_involvement",
    "Human Phenotype Ontology": "msigdb_hp",
}

# Universe dict
universe_dict = {
    "Human Proteome": "human_proteome",
    "HEK293T Core": "hek293t_core",
    "Bind Degs Detected": "cemm_detected",
    "Bind Degs Enriched": "cemm_enriched",
    "Pulldown": "pulldown",
}
