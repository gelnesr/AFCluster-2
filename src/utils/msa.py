from Bio import SeqIO

def load_fasta(fil):
    ''' Read a fasta file and return ids, seqs'''
    IDs, seqs = [], []
    for record in SeqIO.parse(fil, "fasta"):
        IDs.append(record.id)
        seqs.append(str(record.seq))
    return IDs, seqs

def write_fasta(names, seqs, outfile):
    ''' Write a fasta file '''
    with open(outfile, 'w') as f:
        for name, seq in zip(names, seqs):
            f.write(f">{name}\n{seq}\n")

def clean_seqs(seqs):
    return [''.join([x for x in s if x.isupper() or x == '-']) for s in seqs]
