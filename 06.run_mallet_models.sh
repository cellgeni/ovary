#! /bin/bash -e
#BSUB -G cellgeni
#BSUB -J mallet[1-12]
#BSUB -o logs/%J.%I.mallet.log
#BSUB -e logs/%J.%I.mallet.err
#BSUB -n 20
#BSUB -M400000
#BSUB -R "span[hosts=1] select[mem>400000] rusage[mem=400000]"
#BSUB -q long

cd /lustre/scratch127/cellgen/cellgeni/tickets/tic-3942

topics=(2 5 10 15 20 25 30 35 40 45 50 55)
ntopics=${topics[$LSB_JOBINDEX-1]}
singularity run -B /nfs,/lustre /nfs/cellgeni/singularity/images/scenicplus-fa55dae.sif \
 actions/ovary/bin/pycistopic_run_mallet.py \
 --cistopic_pkl_in work/pycistopic/results_pycistopic_call_peaks/combined_cistopic_object.pkl \
 --pkl_out work/pycistopic/results_pycistopic_call_peaks/mallet_models/model_${ntopics}.pkl \
 --n_topics ${ntopics} \
 --ncpu 20 \
 --mem 400