import argparse
import sys
import os
import csv
import re
import textwrap
from re import Pattern
from pathlib import Path
from typing import List, Union, Optional


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, 
                                     usage=f"{sys.argv[0]} -h")
    parser.add_argument('-i', 
                        dest='genome_file', 
                        type=isfile, 
                        required=True, 
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', 
                        dest='min_gene_len', 
                        type=int, 
                        default=50,
                        help="Minimum gene length to consider (default 50)")
    parser.add_argument('-s',
                        dest='max_shine_dalgarno_distance', 
                        type=int, 
                        default=16, 
                        help=("Maximum distance from start codon "
                              "where to look for a Shine-Dalgarno "
                              "motif (default 16)"))
    parser.add_argument('-d', 
                        dest='min_gap', 
                        type=int, 
                        default=40,
                        help=("Minimum gap between two genes - "
                              "shine box not included (default 40)"))
    parser.add_argument('-p', 
                        dest='predicted_genes_file', 
                        type=Path, 
                        default=Path("predict_genes.csv"),
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', 
                        dest='fasta_file', 
                        type=Path,
                        default=Path("genes.fna"),
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file: Path) -> str:
    """Extract genome sequence from fasta files.

    :param fasta_file: (Path) Path to the fasta file.
    :return: (str) Sequence from the genome. 
    """
    sequence = []
    with open(fasta_file, 'r') as file:
        for line in file:
            # Ignore header lines starting with ">"
            if not line.startswith(">"):
                # Remove any newline characters and convert to uppercase
                sequence.append(line.strip().upper())
    return ''.join(sequence)


def find_start(start_regex: Pattern,
               sequence: str,
               start: int,
               stop: int) -> Union[int, None]:
    """Find next start codon before the end position.

    :param start_regexp: A regex object that identifies a start codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :param stop: (int) Stop position of the research
    :return: (int) If exist, position of the start codon. Otherwise None. 
    """
    match = start_regex.search(sequence, start, stop)
    if match:
        return match.start(0)
    return None


def find_stop(stop_regex: Pattern,
              sequence: str,
              start: int) -> Union[int, None]:
    """Find next stop codon in the same reading phase as the start.

    :param stop_regexp: A regex object that identifies a stop codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :return: (int) If exist, position of the stop codon. Otherwise None. 
    """
    # Start looking after the start codon, in steps of 3 to stay in frame
    for match in stop_regex.finditer(sequence, start):
        # Ensure it's in the same reading frame
        if (match.start(0) - start) % 3 == 0:
            return match.start(0)
    return None


def has_shine_dalgarno(shine_regex: Pattern,
                       sequence: str,
                       start: int,
                       max_shine_dalgarno_distance: int) -> bool:
    """Find a shine dalgarno motif before the start codon

    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Position of the start in the genome
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine 
                                              dalgarno to the start position
    :return: (boolean) true -> has a shine dalgarno upstream to the gene, 
                       false -> no
    """
    # Define the region to search for Shine-Dalgarno 
    # (from start-max_shine_dalgarno_distance to start-6)
    search_start = start - max_shine_dalgarno_distance
    search_end = start - 6

    # If the search range is negative, return False
    if search_start < 0:
        return False

    # Search the sequence using the Shine-Dalgarno regex
    match = shine_regex.search(sequence, search_start, search_end)
    
    # Return True if a Shine-Dalgarno sequence is found, otherwise False
    return match is not None


def predict_genes(sequence: str,
                  start_regex: Pattern,
                  stop_regex: Pattern,
                  shine_regex: Pattern, 
                  min_gene_len: int,
                  max_shine_dalgarno_distance: int,
                  min_gap: int) -> List:
    """Predict most probable genes

    :param sequence: (str) Sequence from the genome.
    :param start_regexp: A regex object that identifies a start codon.
    :param stop_regexp: A regex object that identifies a stop codon.
    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param min_gene_len: (int) Minimum gene length.
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine 
                                              dalgarno to the start position.
    :param min_gap: (int) Minimum distance between two genes.
    :return: (list) List of [start, stop] position of each predicted genes.
    """
    predicted_genes = []
    current_position = 0

    # Loop until we have checked the entire sequence with enough gap
    while len(sequence) - current_position >= min_gap:
        # Find the next start codon
        start_codon_pos = find_start(start_regex, 
                                     sequence, 
                                     current_position, 
                                     len(sequence))
        if start_codon_pos is None:
            break  # Stop search if no more start codons are found

        # Find the corresponding stop codon
        stop_codon_pos = find_stop(stop_regex, sequence, start_codon_pos)
        
        # Ensure the gene is long enough
        if stop_codon_pos and (
            stop_codon_pos - start_codon_pos >= min_gene_len
            ):
            # Check for the Shine-Dalgarno sequence upstream of the start codon
            if has_shine_dalgarno(shine_regex, 
                                  sequence, 
                                  start_codon_pos,
                                  max_shine_dalgarno_distance):
                # A probable gene is identified, store its positions
                predicted_genes.append(
                    [
                        start_codon_pos + 1, 
                        stop_codon_pos + 3
                    ]
                    )
                # Move the current position to after the stop codon + min_gap
                current_position = stop_codon_pos + 3 + min_gap
            else:
                # No Shine-Dalgarno found
                current_position += 1
        else:
            # No valid stop codon or gene is too short
            current_position += 1

    return predicted_genes


def write_genes_pos(predicted_genes_file: Path,
                    probable_genes: List[List[int]]) -> None:
    """Write list of gene positions.

    :param predicted_genes_file: (Path) Output file of gene positions.
    :param probable_genes: List of [start, stop] position of predicted genes.
    """
    try:
        with predicted_genes_file.open("wt") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def fill(sequence: str) -> str:
    """Fill the sequence to 60 characters per line

    :param sequence: (str) DNA Sequence.
    :return: (str) Filled sequence.
    """
    return os.linesep.join(textwrap.wrap(sequence, width=60))


def write_genes(fasta_file: Path,
                sequence: str,
                probable_genes: List[List[int]],
                sequence_rc: str, 
                probable_genes_comp: List[List[int]]):
    """Write gene sequence in fasta format

    :param fasta_file: (Path) Output fasta file.
    :param sequence: (str) Sequence of genome file in 5'->3'.
    :param probable_genes: (list) List of [start, stop] position of each 
                                  predicted genes in 5'->3'.
    :param sequence_rc: (str) Sequence of genome file in 3' -> 5'.
    :param probable_genes_comp: (list) List of [start, stop] position of 
                                       each predicted genes in 3' -> 5'.
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for i,gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    i + 1, os.linesep, 
                    fill(sequence[gene_pos[0]-1:gene_pos[1]])))
            i = i + 1
            for j,gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i+1+j, os.linesep,
                            fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement

    :param sequence: (str) DNA Sequence.
    :return: (str) Reverse complemented sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in sequence[::-1]])


#==============================================================
# Main program
#==============================================================
def main() -> None: # pragma: no cover
    """
    Main program function
    """
    
    # Gene detection over genome involves considering thymine instead of uracil
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    
    # Shine-Dalgarno sequence regex (AGGAGGUAA: AGGA or GGAGG)
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    
    # Get arguments from the command line
    args = get_arguments()

    # Read the genome sequence from the provided FASTA file
    sequence = read_fasta(args.genome_file)

    # Predict genes in the 5' to 3' direction
    probable_genes = predict_genes(
        sequence, start_regex, stop_regex, shine_regex,
        args.min_gene_len, args.max_shine_dalgarno_distance, args.min_gap
    )

    # Predict genes in the 3' to 5' direction (reverse complement)
    sequence_rc = reverse_complement(sequence)
    probable_genes_comp = predict_genes(
        sequence_rc, start_regex, stop_regex, shine_regex,
        args.min_gene_len, args.max_shine_dalgarno_distance, args.min_gap
    )

    # Correct reverse complement positions and combine with 5' to 3' prediction
    corrected_genes_comp = [
        [len(sequence) - end + 1, len(sequence) - start + 1]
        for start, end in probable_genes_comp
    ]
    
    # Combine the predicted genes from both directions
    all_predicted_genes = probable_genes + corrected_genes_comp
    all_predicted_genes.sort()

    # Write predicted gene positions to file
    write_genes_pos(args.predicted_genes_file, all_predicted_genes)

    # Write the gene sequences (original and reverse complement)
    write_genes(args.fasta_file,
                sequence,
                probable_genes,
                sequence_rc,
                corrected_genes_comp)


if __name__ == "__main__":
    main()
