def fasta_list(filepath: str):
        with open(filepath, 'r') as data:
            file = data.read()
            sequences = file.split('>')[1:]
            list_seqs = []
            for sequence in sequences:
                seq_parts = sequence.partition('\n')
                header = seq_parts[0]
                seq = seq_parts[2].replace('\n', '')
                list_seqs.append((header, seq))
            return list_seqs