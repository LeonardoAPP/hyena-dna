# Finetunning on A. Thaliana chromatine profile sequences

### Requirements
* A. Thaliana genome `.fasta` file (`TAIR10_mod.fna`)
* Chromatin profile annotations 
    - `Arabidopsis_thaliana.TAIR10.57.gff3` --> introns, exons, etc...
    - `TAIR10_mod.Rfam.gff3` --> non-coding RNA
    - `TAIR10_Transposable_Elements.bed` --> transposons

* Create a `.bed` file with annotations `chr_n`, `strart`, `end`, `label`, `split` [0,1,2], `strand` [+,-,.].
* Create a script that loop thru `.bed` and `.fasta` file and return a sequence of length `max_length` and its label from a random entry of the `.bed` file. (**Done**)
* The training dataset have to be a list of nucleotides sequences and a list of labes corresponding to each sequence (**Done**)

* Base line experiment: genes vs inter-genic regions
    - gene seq from tss 
    - inter-genic regions far from centromers