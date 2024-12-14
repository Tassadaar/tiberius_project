
import sys, os
# ADD PATH TO LEARN MSA HERE (if not installed with pip)
learn_msa_path = ''
sys.path.insert(0, learn_msa_path)

sys.path.append("bin")
from bin.genome_anno import Anno
import numpy as np
import tensorflow as tf
import tensorflow.keras as keras
from bin.eval_model_class import PredictionGTF
from Bio import SeqIO
from Bio.Seq import Seq
from bin.tiberius import assemble_transcript, check_in_frame_stop_codons

# CHANGE MODEL PATH IF NEEDED
model_path = "data/tiberius_model_CM.h5"
batch_size = 2
seq_len = 500004
strand = '+'

emb = False
hmm_parallel = 817
trans_lstm = False

inp_data_dir = 'inp/'
out_dir = 'test_train/'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

genome_path = f'data/Ergobibamus_cyprinoides_CL.scaffolds.edited.fasta'
# output gtf file
gtf_out = 'tiberius.out'

pred_gtf = PredictionGTF(
    model_path=model_path,
    seq_len=seq_len,
    batch_size=batch_size,
    hmm=True,
    emb=False,
    num_hmm=1,
    hmm_factor=1,
    genome_path=genome_path,
    softmask=True, strand=strand,
    parallel_factor=hmm_parallel
)

# load model
pred_gtf.load_model()

# load input data x_seq
x_seq, y_seq, coords = pred_gtf.load_inp_data(
    strand=strand,
    chunk_coords=True, softmask=True
)

# generate LSTM and HMM predictions
hmm_pred = pred_gtf.get_predictions(x_seq, hmm_filter=True)

# infer gene structures and write GTF file
anno, tx_id = pred_gtf.create_gtf(y_label=hmm_pred, coords=coords,
                                  out_file=gtf_out, f_chunks=x_seq, strand=strand)

# Filter results and write to gtf
gtf_out = 'tiberius.gtf'
genome = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
anno_outp = Anno('', f'anno')
out_tx = {}
for tx_id, tx in anno.transcripts.items():
    exons = tx.get_type_coords('CDS', frame=False)
    filt=False

    # filter out tx with inframe stop codons
    coding_seq, prot_seq = assemble_transcript(exons, genome[tx.chr], tx.strand )
    if not coding_seq or check_in_frame_stop_codons(prot_seq):
        filt = True
    # filter out transcripts with cds len shorter than args.filter_short
    if not filt and tx.get_cds_len() < 201:
        filt = True

    if not filt:
        out_tx[tx_id] = tx

anno_outp.add_transcripts(out_tx, f'anno')
anno_outp.norm_tx_format()
anno_outp.find_genes()
anno_outp.rename_tx_ids()
anno_outp.write_anno(gtf_out)