from pathlib import Path
from typing import List

transitions = {
    'A': ['G'],
    'C': ['T'],
    'T': ['C'],
    'G': ['A'],
}


transversions = {
    'A': ['C', 'T'],
    'C': ['A', 'G'],
    'T': ['A', 'G'],
    'G': ['C', 'T'],
}


def is_transition(n1: str, n2: str) -> bool:
    return n1 in transitions[n2]


def is_transversion(n1: str, n2: str) -> bool:
    return n1 in transversions[n2]


def download_sequence(seq_id: str) -> str:
    seq_filename = f'{seq_id}.fasta'
    if Path(seq_filename).exists():
        return seq_filename

    from Bio import Entrez
    print(f'[*] Downloading sequence {seq_id}')
    Entrez.email = 'A.N.Other@example.com'
    handle = Entrez.efetch(db='nucleotide', id=seq_id, rettype='fasta', retmode='text')
    record = handle.read()
    if not record:
        raise Exception(f'[-] Downloading sequence {seq_id} FAILED')
    with open(seq_filename, 'w') as f:
        f.write(record)
    print(f'[+] Downloaded sequence stored at {seq_filename}')
    return seq_filename


def read_sequences(filename: str) -> List[str]:
    from Bio import SeqIO
    print(f'[*] Reading sequences from {filename}')
    seqs = [str(record.seq) for record in SeqIO.parse(filename, 'fasta')]
    print(f'[+] Read {len(seqs)} sequences')
    return seqs


def plot(data_lists: List[List[int]], data_labels: List[str], /, *,
         title: str = '', xlabel: str = '', ylabel: str = '',
         output_file: str = ''):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for data, label in zip(data_lists, data_labels):
        ax.plot(data, label=label)
    ax.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if output_file:
        fig.savefig(output_file, dpi=fig.dpi)
    plt.show()
