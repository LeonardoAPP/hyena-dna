# Finetunning on A. Thaliana chromatine profile sequences

### Requirements
* A. Thaliana genome `.fasta` file (`TAIR10_mod.fna`)
* Chromatin profile annotations 
    - `Arabidopsis_thaliana.TAIR10.57.gff3` --> introns, exons, RNAs...
    - `TAIR10_mod.Rfam.gff3` --> non-coding RNA
    - `TAIR10_Transposable_Elements.bed` --> transposons

* Split the genome previously (to prevent data leakage)
    - 3/4 train
    - 1/8 validation
    - 1/8 test
* Create a `.bed` file with annotations `chr_n`, `strart`, `end`, `label`, `split` [0,1,2], `strand` [+,-,.]. (***Done***)
* Create a script that loop thru `.bed` and `.fasta` file and return a sequence of length `max_length` and its label from a random entry of the `.bed` file. (***Done***)
    - When the `seq_len` from the `.bed` file is less than `max_length`, the sequence is extended to fill the `max_length` 
        - If `label == 'intergenic'`, the sequence is extended with random data
        - If `label == 'gene'`, the sequence is extended with cromosome data
* The training dataset have to be a list of nucleotides sequences and a list of labes corresponding to each sequence (***Done***)

* Base-line experiment: genes vs intergenic regions
    - gene seq from gff (***Done***)
    - intergenic regions far from centromers (***Done***)
    - Study trainig loss vs validation loss (***Done***)
    - Study performance vs ds_train_length (***Done***)
    - Evaluate performance vs `max_length`
    - Evaluate the effects of random data in `intergenic` on model performance
    - Create a `pipeline(raw_data): -> train_data_set`  