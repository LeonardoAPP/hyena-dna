from Bio import SeqIO
import random

def create_bed_file(output_file, num_records):
    chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']
    labels = ['gene', 'intergenic']
    splits = [0, 1, 2]
    strands = ['+', '-', '.']

    with open(output_file, 'w') as bed_file:
        bed_file.write(f"{chromosomes[0]}\t{50_000}\t{50_100}\t{labels[0]}\t{2}\t{'-'}\n")
        for _ in range(num_records):
            chr = random.choice(chromosomes)
            start = random.randint(0, 1000000)
            end = start + random.randint(100, 10000)
            label = random.choice(labels)
            split = random.choice(splits)
            strand = random.choice(strands)
            bed_file.write(f"{chr}\t{start}\t{end}\t{label}\t{split}\t{strand}\n")

create_bed_file('A_Thaliana.bed', 1000)