# Overview
Here is code for analyses of ovary ATAC data. The overall plan is to 
a) calculate gen score
b) intergate and label transfer from single cell atlas
c) perform differential accessability analyses for granulosa cells
The analyses is performed using snapATAC2 (s) and pyscenic+ (p) in combination or in parallel. 

# The pipeline
Following pipeline was used:
1. Combine samples and calculate gene scores.
    snapATAC2: see 01.snapatac2_combine.sh (based on https://github.com/cellgeni/nf-atac/tree/snapatac2)
    pycistopic: see 01.pycistopic_combine.sh (to generate cistopic object) + bin/calc_gene_scores_pycistopic.py (to calculate scores)
2. Combine with single cell reference (02.prepare_for_scvi.ipynb)
3. Use scanvi for integration and label transfer with dataset_donor as batch (03.run_scanvi.sh)
4. Use 04.clean_scvi_results.ipynb to clean up transfered annotation. Details: for some reason scanvi label transfer misannotated many reference cells (for examples it annotated all Stroma_InnCor as some Immune cells) so I used bbknn on scanvi embedding to transfer labels bby majority voting of 50 nearest reference neighbors. Integration for pycistopic-generated gene scores worked poorly, so I continued using annotation based on snapATAC gene scores.
5. Call and quantify peaks using celltype annotation (05.snapatac2_peak_calling.sh or 05.pycistopic_peak_calling.sh). snapatac2 found 49186 peaks while pycistopic returned 645326 peaks, 43743 of them are common (more than 50% overlap). 
6. Search for differentially accesable regions within Granulosa cells (only Granulosa_sq, Granulosa_sq_atr, Granulosa_sq_transitioning, Granulosa_AMH_early, Granulosa_AMH_ml) using scanpy (see 06.scanpy_DARs.ipynb). Identified marker regions are in snapatac2_granulosa_dars.csv and pycistopic_granulosa_dars.csv.
