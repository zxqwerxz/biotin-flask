"""Fasta class.

This module contains classes and methods used for interacting with the typical
fasta file used in targeted NextGen sequencing. In particular, these fasta
files often contain multiple gene sequences in one fasta file.

"""

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Development'


class Fasta:
    """Represents a FASTA file.

    Attributes:
        genes (dict): A dictionary of genes. See example:

        example_genes_instance_variable = {
            'geneA': 'ACTGTCTCTTTTT',
            'geneB': 'TCTGGGTCTCCAAAATC'
        }

    """

    def __init__(self, filepath):
        """Construct a Fasta object from a fasta file path."""
        self.genes = {}
        try:
            with open(filepath, 'r') as file_h:
                prev_id = ''
                prev_seq = ''
                for line in file_h.readlines():
                    line = line.rstrip()
                    if line[0] == '>':
                        # Add the previous gene
                        if prev_id and prev_seq:
                            if prev_id in self.genes:
                                raise IOError('')
                            else:
                                self.genes[prev_id] = prev_seq

                        # Start the new gene
                        prev_id = line[1:]
                        continue
                    if prev_id:
                        prev_seq = prev_seq + line

                # Don't forget to add the last gene
                if prev_id in self.genes:
                    raise IOError('')
                else:
                    self.genes[prev_id] = prev_seq
        except:
            raise IOError('Unable to read FASTA file. Is it corrupt?')
        if len(self.genes) == 0:
            raise IOError('Unable to read FASTA file. Is it empty?')

    def validate(self, gene_list):
        """Validate if the fasta file contains all of the genes in the list.

        Parameter:
            gene_list (list): A list of string gene names.

        Returns:
            True if the fasta file contains all of the genes.

        """
        for gene in gene_list:
            if gene not in self.genes:
                return False
        return True
