
from pathlib import Path
from pyfaidx import Fasta
import polars as pl
import pandas as pd
import torch
from random import randrange, random, choices
import random as rand
import numpy as np


"""

Dataset for sampling arbitrary intervals from the Arabidobsis Thaliana genome.

"""


# helper functions

def exists(val):
    return val is not None

def coin_flip():
    return random() > 0.5

# augmentations

string_complement_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}

def string_reverse_complement(seq):
    rev_comp = ''
    for base in seq[::-1]:
        if base in string_complement_map:
            rev_comp += string_complement_map[base]
        # if bp not complement map, use the same bp
        else:
            rev_comp += base
    return rev_comp


class FastaInterval():
    def __init__(
        self,
        *,
        fasta_file,
        # max_length = None,
        return_seq_indices = False,
        shift_augs = None,
        rc_aug = False,
        pad_interval = False,
        pad_only = False,  # NUEVO: solo pad con puntos, no extender con datos reales
    ):
        fasta_file = Path(fasta_file)
        assert fasta_file.exists(), 'path to fasta file must exist'

        self.seqs = Fasta(str(fasta_file))
        self.return_seq_indices = return_seq_indices
        # self.max_length = max_length # -1 for adding sos or eos token
        self.shift_augs = shift_augs
        self.rc_aug = rc_aug
        self.pad_interval = pad_interval 
        self.pad_only = pad_only

        # calc len of each chromosome in fasta file, store in dict
        self.chr_lens = {}

        for chr_name in self.seqs.keys():
            # remove tail end, might be gibberish code
            # truncate_len = int(len(self.seqs[chr_name]) * 0.9)
            # self.chr_lens[chr_name] = truncate_len
            self.chr_lens[chr_name] = len(self.seqs[chr_name])


    def __call__(self, chr_name, start, end, max_length, random_aug = False):
        """
        max_length passed from dataset, not from init
        """
        interval_length = end - start
        chromosome = self.seqs[chr_name]
        chromosome_length = self.chr_lens[chr_name]
        rand_seq_left = ''
        rand_seq_right = ''

        if exists(self.shift_augs):
            min_shift, max_shift = self.shift_augs
            max_shift += 1

            min_shift = max(start + min_shift, 0) - start
            max_shift = min(end + max_shift, chromosome_length) - end

            rand_shift = randrange(min_shift, max_shift)
            start += rand_shift
            end += rand_shift

        left_padding = right_padding = 0

        # checks if not enough sequence to fill up the start to end
        if interval_length < max_length:
            extra_seq = max_length - interval_length

            extra_left_seq = extra_seq // 2
            extra_right_seq = extra_seq - extra_left_seq

            # if random_aug is enabled, add random sequence to the left and right
            if random_aug:
                rand_seq_left = str.join('', choices(['A', 'C', 'G', 'T'], k=extra_left_seq))
                rand_seq_right = str.join('', choices(['A', 'C', 'G', 'T'], k=extra_right_seq))
            # if pad_only is enabled, pad with '.' only
            elif self.pad_only:
                left_padding = extra_left_seq
                right_padding = extra_right_seq
            # else, extend sequence left and right with the data from the chromosome
            else:
                start -= extra_left_seq
                end += extra_right_seq

        if start < 0:
            left_padding = -start
            start = 0

        if end > chromosome_length:
            right_padding = end - chromosome_length
            end = chromosome_length

        # Added support!  need to allow shorter seqs
        if interval_length > max_length:
            end = start + max_length

        seq = rand_seq_left + str(chromosome[start:end]) + rand_seq_right

        # if reverse complement augmentation is enabled, 
        # reverse complement the sequence with 50% probability 
        if self.rc_aug and coin_flip():
            seq = string_reverse_complement(seq)

        if self.pad_interval or self.pad_only:
            seq = ('.' * left_padding) + seq + ('.' * right_padding)
        
        # Uppercase normalization
        seq = seq.upper()

        return seq

class a_thalinana_Dataset(torch.utils.data.Dataset):

    '''
    Loop thru bed file, retrieve (chr, start, end, 'label' ,'split', 'strand'), query fasta file for sequence.
    
    '''

    def __init__(
        self,
        split,
        bed_file,
        fasta_file,
        max_length,
        pad_max_length = None,
        tokenizer = None,
        tokenizer_name = None,
        add_eos = False,
        return_seq_indices = False,
        shift_augs = None,
        rc_aug = False,
        rand_aug = True, # random augmentation for intergenic regions
        rc_strand = False, # reverse complement the sequence in strand -
        deterministic_mode = False, # deterministic mode for random augmentations
        replace_N_token = False,  # replace N token with pad token
        pad_interval = False,  # options for different padding
        pad_only = False,
    ):

        self.max_length = max_length
        self.pad_max_length = pad_max_length if pad_max_length is not None else max_length
        self.tokenizer_name = tokenizer_name
        self.tokenizer = tokenizer
        self.add_eos = add_eos
        self.replace_N_token = replace_N_token  
        self.pad_interval = pad_interval 
        self.split_dict = {'train': 0, 'val': 1, 'test': 2}
        self.rc_strand = rc_strand     
        self.rand_aug = rand_aug
        self.deterministic_mode = deterministic_mode


        bed_path = Path(bed_file)
        assert bed_path.exists(), 'path to .bed file must exist'

        # read bed file
        df_raw = pd.read_csv(str(bed_path), sep = '\t', names=['chr_name', 'start', 'end', 'label' ,'split', 'strand'], header=0)        
        # select only split df
        self.df = df_raw[df_raw['split'] == self.split_dict[split]]

        self.fasta = FastaInterval(
            fasta_file = fasta_file,
            # max_length = max_length,
            return_seq_indices = return_seq_indices,
            shift_augs = shift_augs,
            rc_aug = rc_aug,
            pad_interval = pad_interval,
            pad_only = pad_only,  
        )

        classes = self.df['label'].unique().tolist()
        self.label_dict = {}
        for i, label in enumerate(classes):
            self.label_dict[label] = i

    def __len__(self):
        return len(self.df)

    def replace_value(self, x, old_value, new_value):
        return torch.where(x == old_value, new_value, x)

    def __getitem__(self, idx):
        """Returns a sequence of specified len"""
        # Fijar semilla basada en idx para hacer aumentos aleatorios deterministas
        if self.deterministic_mode:
            seed = 42 + idx  # Puedes usar cualquier número base (42) + idx
            torch.manual_seed(seed)
            np.random.seed(seed)
            rand.seed(seed)  # Si usas `random` en FastaInterval

        # sample a random row from df
        row = self.df.iloc[idx]
        # row = (chr, start, end, label ,split, strand)
        chr_name, start, end, seq_label, seq_strand = (row[0], row[1], row[2], row[3], row[5])

        # hack, random augmentation for intergenic regions
        if seq_label == 'intergenic':
            seq = self.fasta(chr_name, start, end, max_length=self.max_length, random_aug = self.rand_aug)
        else:
            seq = self.fasta(chr_name, start, end, max_length=self.max_length)

        if self.rc_strand & (seq_strand == '-'):
            seq = string_reverse_complement(seq)

        if self.tokenizer_name == 'char':

            seq = self.tokenizer(seq,
                add_special_tokens=True if self.add_eos else False,  # this is what controls adding eos
                padding="max_length",
                max_length=self.max_length,
                truncation=True,
            )
            seq = seq["input_ids"]  # get input_ids

        elif self.tokenizer_name == 'bpe':
            seq = self.tokenizer(seq, 
                # add_special_tokens=False, 
                padding="max_length",
                max_length=self.pad_max_length,
                truncation=True,
            ) 
            # get input_ids
            if self.add_eos:
                seq = seq["input_ids"][1:]  # remove the bos, keep the eos token
            else:
                seq = seq["input_ids"][1:-1]  # remove both special tokens
        
        # convert to tensor
        seq = torch.LongTensor(seq)  # hack, remove the initial cls tokens for now

        if self.replace_N_token:
            # replace N token with a pad token, so we can ignore it in the loss
            seq = self.replace_value(seq, self.tokenizer._vocab_str_to_int['N'], self.tokenizer.pad_token_id)

        data = seq.clone()  # the whole sequence

        label = torch.LongTensor([self.label_dict[seq_label]])
        target = label.clone()  # the label

        return {"input_ids": data, "labels": target}
