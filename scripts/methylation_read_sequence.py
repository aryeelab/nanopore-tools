import sys
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

def iter_chunk_by_read_name(file):
    csv_reader = pd.read_csv(file, iterator=True, chunksize=1, sep='\t')
    first_chunk = csv_reader.get_chunk()
    id = first_chunk.iloc[0]["read_name"]
    chunk = pd.DataFrame(first_chunk)
    for l in csv_reader:
        if id == l.iloc[0]["read_name"]:
            chunk = chunk.append(l)
            continue
        id = l.iloc[0]["read_name"]
        yield chunk
        chunk = pd.DataFrame(l)
    yield chunk

of = open(outfile, 'w')

read_iter = iter_chunk_by_read_name(infile)
for read in read_iter:
    read_start_pos = min(read["start"])
    read["s"] = read["start"] - read_start_pos
    read["length"] = read["sequence"].str.len()
    read["e"] = read["s"] + read["length"]
    read["pad"] = read["s"] - read["e"].shift()
    read.loc[read.index[0], 'pad'] = 0
    seqs = []
    for index, kmer in read.iterrows():
        seq = kmer['sequence']
        if isinstance(seq, str):
            if kmer['log_lik_ratio'] >= 2.5:
                seq = seq.replace("CG", "MG")
            elif kmer['log_lik_ratio'] <= -2.5:
                seq = seq.replace("CG", "UG")
            else:
                seq = seq.replace("CG", "?G")                
            seqs.append('-' * int(kmer['pad']))
            seqs.append(seq)
    sequence = ''.join(seqs)
    read_end_pos = read_start_pos + len(sequence)
    chromosome = read["chromosome"].iloc[0]
    read_name = read["read_name"].iloc[0]
    of.write('\t'.join([str(chromosome), str(read_start_pos), str(read_end_pos), read_name, sequence]) + "\n")

of.close()