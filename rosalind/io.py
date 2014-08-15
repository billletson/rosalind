import urllib2
from .sequences import *


def load_fasta_file(f, s_type="DNA"):
    """
    Reads a FASTA formatted file (or file-like object) and creates a list of
    sequence objects
    Arguments: file f,str s_type
    Returns: Sequence[] (DNA[],RNA[],or Protein[])
    """
    sequences = []
    tmp_string = ""
    name = ""
    for line in f:
        line = line.strip()
        if line[0] == ">":
            if tmp_string:
                if s_type == "RNA":
                    sequences.append(RNA(name, tmp_string))
                elif s_type == "Protein":
                    sequences.append(Protein(name, tmp_string))
                else:
                    sequences.append(DNA(name, tmp_string))
            name = line[1:]
            tmp_string = ""
        else:
            tmp_string += line
    if tmp_string:
        if s_type == "RNA":
            sequences.append(RNA(name, tmp_string))
        elif s_type == "Protein":
            sequences.append(Protein(name, tmp_string))
        else:
            sequences.append(DNA(name, tmp_string))
    return sequences


def load_fasta_uniprot(uniprot_id, s_type="Protein"):
    """
    Given a uniprot id, creates a sequence object
    Arguments: string uniprot_id, str s_type
    Returns: Sequence (DNA,RNA,or Protein)
    """
    response = urllib2.urlopen("http://www.uniprot.org/uniprot/%s.fasta" % uniprot_id)
    response.next()
    sequence = "".join(line.strip() for line in response)
    if s_type == "RNA":
        return RNA(uniprot_id, sequence)
    elif s_type == "DNA":
        return DNA(uniprot_id, sequence)
    else:
        return Protein(uniprot_id, sequence)


def load_scoring_matrix(matrix_type):
    """
    Returns a dict representing a scoring matrix. Inversion is because this package uses modifications of the
    levenshtein distance rather than postive alignment scores under the covers.
    """
    matrix = {}
    with open(os.path.dirname(__file__) + ("/reference_data/%s.txt" % matrix_type), 'rb') as f:
        header = f.next().strip().split()
        for line in f:
            values = line.strip().split()
            for i in xrange(1, len(values)):
                matrix[(header[i - 1], values[0])] = int(values[i])
    return matrix

