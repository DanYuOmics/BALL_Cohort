#!/usr/bin/env python3
# ==============================================================================
# Script: 02_pyclone_deconvolution.py
# Description: Execute PyClone-VI variational inference on pooled 
#              cohort inputs, extract latent variable parameters (MAP cluster IDs), 
#              perform Dirichlet normalization to compute Cellular Fraction (CCF),
#              and serialize the clonal architecture table.
# ==============================================================================

import os
import sys
import subprocess
import h5py
import numpy as np
import pandas as pd

# ------------------------------------------------------------------------------
# 1. Environment Setup & Directory Configurations
# ------------------------------------------------------------------------------
INPUT_TSV = os.getenv("PYCLONE_INPUT", "./article/PyClone_VI_Main_Cohort_Input_pooled.tsv")
H5_OUTPUT = os.getenv("PYCLONE_H5_OUT", "./article/Fig2_Upgraded_Clonal_Architecture.tsv")
TXT_OUTPUT = os.getenv("PYCLONE_TXT_OUT", "./article/Fig2_Upgraded_Clonal_Architecture_TXT.tsv")

# Ensure parent directory exists
os.makedirs(os.path.dirname(H5_OUTPUT), exist_ok=True)

# ------------------------------------------------------------------------------
# 2. Step A: Fit PyClone-VI Variational Model
# ------------------------------------------------------------------------------
def run_pyclone_fit(input_file, output_h5):
    """
    Run PyClone-VI variational model with a maximum of 10 mixture components
    and 20 independent random restarts.
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"⚠️ Error: PyClone-VI input file not found at: {input_file}")

    print(f"⚡ Step A: Fitting PyClone-VI variational model on input: {input_file}")
    
    cmd = [
        "pyclone-vi", "fit",
        "-i", input_file,
        "-o", output_h5,
        "-c", "10",
        "-r", "20"
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"✅ Successfully fitted model. Saved HDF5 parameters to: {output_h5}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error executing PyClone-VI command: {e}")
        sys.exit(1)

# ------------------------------------------------------------------------------
# 3. Step B: Mathematical Normalization & Extraction of Clonal Metrics
# ------------------------------------------------------------------------------
def extract_clonal_architecture(h5_path, txt_out_path):
    """
    Extract mutation identifiers, MAP cluster assignments from the z-matrix, 
    and normalize Dirichlet mixture weights (pi) into clean CCF values.
    """
    print(f"🧬 Step B: Parsing HDF5 parameters and extracting subclonal architecture...")

    if not os.path.exists(h5_path):
        raise FileNotFoundError(f"⚠️ Error: PyClone HDF5 output file missing at: {h5_path}")

    with h5py.File(h5_path, 'r') as f:
        # 1. Retrieve unique mutation identifiers
        mutations = [m.decode('utf-8') for m in f['data']['mutations'][:]]
        
        # 2. Extract latent variable matrix 'z' and apply MAP (Maximum A Posteriori)
        z_matrix = f['var_params']['z'][:]
        cluster_ids = np.argmax(z_matrix, axis=1)
        
        # 3. Retrieve raw mixture weights 'pi' (Dirichlet alpha parameters)
        pi_matrix = f['var_params']['pi'][:]
        
        # 4. Perform Dirichlet normalization to map alpha counts into a clean [0, 1] probability space
        ccf_proportions = pi_matrix / np.sum(pi_matrix)
        
        # 5. Map normalized cellular proportions back to corresponding somatic mutations
        cellular_fractions = []
        for idx, c_id in enumerate(cluster_ids):
            val = ccf_proportions[c_id]
            cellular_fractions.append(float(val))

    # 6. Construct human-readable DataFrame and serialize to disk
    df_out = pd.DataFrame({
        'mutation_id': mutations,
        'sample_id': 'BALL_Cohort',
        'cluster_id': cluster_ids,
        'cellular_fraction': cellular_fractions
    })

    df_out.to_csv(txt_out_path, sep='\t', index=False)
    print(f"🎉 Extraction and mathematical Dirichlet normalization successfully completed!")
    print(f"📄 Serialized subclonal architecture table to: {txt_out_path}")

# ------------------------------------------------------------------------------
# 4. Main Driver Execution
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    # Step A: Run PyClone-VI fitting
    run_pyclone_fit(INPUT_TSV, H5_OUTPUT)
    
    # Step B: Post-process HDF5 output into TSV
    extract_clonal_architecture(H5_OUTPUT, TXT_OUTPUT)
