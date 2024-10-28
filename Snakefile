""" Snake pipeline for running RNAseq analysis """

# Configure shell for all rules 
shell.executable("/bin/bash")
shell.prefix("set -o nounset -o pipefail -o errexit -x; ")
import subprocess # check if necessary
import glob
import os 
import re
import pandas as pd
from itertools import combinations


# Set apptainer bindings NOT A PERMANENT SOLUTION
os.environ["APPTAINER_BIND"] = "/scratch/alpine:/scratch/alpine,/pl/active/Anschutz_BDC:/pl/active/Anschutz_BDC,/tmp:/tmp"

# Parameters from config.yaml
SCRATCH_DIR    = config["SCRATCH_DIR"]
BASE_PATH      = config["BASE_PATH"]
SAMPLE_TABLE      = config["SAMPLE_TABLE"]
GENOME            = config["GENOME"]
COMBINED_GENOME   = config["COMBINED_GENOME"]
GTF               = config["GTF"]
BED               = config["BED"]
COMBINED_GTF      = config["COMBINED_GTF"]
RESULTS           = config["RESULTS"]
ADAPTORS          = config["ADAPTORS"]
TRIM_METHOD       = config["TRIM_METHOD"]
PROJECT           = config["PROJECT"]
FASTQ_SCREEN_DATA = config["FASTQ_SCREEN_DATA"]
RMATS_INFO        = config["RMATS_INFO"]
GENERAL_CONTAINER = config["GENERAL_CONTAINER"]
FASTQ_SCREEN      = config["FASTQ_SCREEN"]
R_CONTAINER       = config["R_CONTAINER"]


# Pull out sample names and fastq files
SAMPLE_LIST = pd.read_table(config["SAMPLE_TABLE"]).set_index("sample", drop = False)
SAMPLES = SAMPLE_LIST.index.values
IS_PAIRED = "fastq2" in SAMPLE_LIST.columns

# Pull out if there are any samples with spike-ins
if "spike_in" in SAMPLE_LIST.columns:
    SPIKE_SAMPLES = SAMPLE_LIST[SAMPLE_LIST.spike_in].index.values
    NON_SPIKE_SAMPLES = SAMPLE_LIST[~SAMPLE_LIST.spike_in].index.values
else:
    NON_SPIKE_SAMPLES = SAMPLE_LIST
    SPIKE_SAMPLES = []

if IS_PAIRED:
    CMD_PARAMS = config["PE"]
else:
    CMD_PARAMS = config["SE"]

# Function to check paths for input files/directories
def _check_path(path):
    if os.path.exists(path):
        return os.path.abspath(path)
    else:
        sys.exit("ERROR: " + path + " does not exist.")

# Make the desired path in scratch
#user_variable = os.environ.get('USER')
user_variable = "kwellswrasman@xsede.org"
scratch_dir = os.path.join(SCRATCH_DIR, user_variable)
current_dir = os.getcwd()
rel_path_to_home = os.path.relpath(current_dir, start=BASE_PATH)

RESULTS2 = os.path.join(scratch_dir, rel_path_to_home, RESULTS)
os.makedirs(RESULTS2, exist_ok = True)

# Set sample/group names
SAMPLES = [x.strip() for x in SAMPLES]


# Check/set directory/file paths
if not os.path.exists(RESULTS):
    os.makedirs(RESULTS)
RESULTS  = _check_path(RESULTS)
GENOME   = _check_path(GENOME)

# Only check the combined genome if combined samples are present in the sample sheet
if len(SPIKE_SAMPLES) > 0:
    COMBINED_GENOME = _check_path(COMBINED_GENOME)

# Build rMats
# Write a function that takes all pairwise groups
# 1. Change the first rule to only make the input of one group (add in a group variable)
# 2. Change the second rule to take the group files as input based on the group variable
# This needs to have the scructure groupx_groupy in the wildcards
RMATS_SAMPLES = pd.read_table(RMATS_INFO).set_index("sample", drop = False)
GROUPS = RMATS_SAMPLES["group"].unique()
pairwise_combinations = list(combinations(GROUPS, 2))
group_comparisons = []
for combination in pairwise_combinations:
    group_comparisons.append(f"{combination[0]}_{combination[1]}")

# Final output files
rule all:
    input:
         # Fastqc summary
        expand(
            "{results}/fastqc_pre_trim_summary_untrimmed.tsv",
            results = RESULTS
            ),
        # Run fastq screen
        expand(
            "{results}/fastq_screen/fastq_screen_barplot.pdf",
            results = RESULTS
        ),
        # Output of adaptor trimming
        expand(
            "{results}/{trim_method}_trim/{sample}.txt",
            results = RESULTS, sample = SAMPLES,
            trim_method = TRIM_METHOD),
        # Fastqc summary
        expand(
            "{results}/fastqc_{trim_method}_summary_trimmed.tsv",
            results = RESULTS, trim_method = TRIM_METHOD
            ),
        # STAR output
        expand(
            "{results}/star_summary_{trim_method}_trim.tsv",
            results = RESULTS, trim_method = TRIM_METHOD
            ),
        # Count table
        expand(
            "{results}/{project}_countTable_{trim_method}_trim.txt",
            results = RESULTS, project = PROJECT,
            trim_method = TRIM_METHOD
            ),

        # Run rmats
        expand(
           "{results}/rmats_{trim_method}_trim/{group}_input.txt",
            results = RESULTS, trim_method = TRIM_METHOD, group = GROUPS
        ),
        expand(
            "{results}/rmats_{trim_method}_trim/rmats_{group_compare}_done.txt",
            results = RESULTS, trim_method = TRIM_METHOD, group_compare = group_comparisons
        ),
        expand(
            "{results}/rmats_{trim_method}_trim/{group_compare}/combined_JCEC.txt",
            results = RESULTS, trim_method = TRIM_METHOD, group_compare = group_comparisons
        ),
        expand(
          "{results}/{project}_quality_{trim_method}_trim.html",
          results = RESULTS, project = PROJECT,
          trim_method = TRIM_METHOD
        )

# Snakes to run
include: "src/rules/fastqc.snake"
include: "src/rules/trimming.snake"
include: "src/rules/star.snake"
include: "src/rules/featurecounts.snake"
include: "src/rules/run_rmats.snake"
include: "src/rules/fastq_screen.snake"
include: "src/rules/qc_plots.snake"
