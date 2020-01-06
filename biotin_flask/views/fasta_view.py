"""
This module takes a list of gene names as input,
(1) searches NCBI Nucleotide for the IDs of genes and their sequences
(2) fetches the FASTA files from the nucleotide IDs and gene IDs
and outputs a FASTA file with every gene and its sequence.
Default sequence includes 5kb upstream and 2kb downstream flanking sequence.
"""

import os, io, csv, StringIO, zipfile
from flask import render_template, request, flash, make_response, send_file
from werkzeug import secure_filename
from Bio import Entrez
from xml.dom import minidom

from biotin_flask import app

__author__ = 'Eric Zhou'
__copyright__ = 'Copyright (C) 2019, EpigenDx Inc.'
__credits__ = ['Eric Zhou']
__version__ = '0.0.1'
__status__ = 'Production'

Entrez.email = 'zhou.eric.y@gmail.com'

@app.route('/misc/fasta', methods=['GET', 'POST'])
def fasta():
    """Handle GET or POST requests to the psq URL.

    Form args:
        gene:   (string) a list of gene names

    Returns:
        A CSV file.

    """

    # Render an empty form for GET request
    if request.method == 'GET':
        return render_template('fasta/form.html')

    # Otherwise validate the form on a POST request and process uploaded files
    species = request.form.get('species')
    gene = request.form.get('gene').splitlines()

    # Obtain FASTA file from NCBI
    id_list, gene_success = uid(species, gene)
    document = fetch(id_list, gene_success)

    # Modify the FASTA file

    return render_template('fasta/form.html')

def uid(species, gene):
    """
    Searches NCBI Nucleotide and RefSeqGene for unique IDs of gene names using ESearch.
    This step is necessary because EFetch takes unique IDs as an argument, not gene names.

    :param
        species:        (string) the species
        gene:           (List, string) a list of gene names

    :return:
        id_list         (List, string) a list of unique IDs corresponding to the gene names
        gene_success    (List, string) a list of gene names for which NCBI successfully retrieved IDs
    """

    id_list = []
    gene_success = []

    for gene in gene:
        try:
            response = Entrez.esearch(db='nucleotide', term=species+'[Organism] '+gene+'[GENE] RefSeqGene[KYWD]', )
            doc = minidom.parseString(response.read())

            try:
                id_list.append(doc.getElementsByTagName('Id')[0].firstChild.nodeValue)
                gene_success.append(gene)
            except:
                flash('NCBI did not return any results for ' + species + ' ' + gene +
                      '. It is possible that NCBI was overloaded. ' + species + ' ' + gene + ' skipped.', 'alert-warning')

        except:
            flash('NCBI was unable to search for ' + species + ' ' + gene +
                  '. It is possible that NCBI was overloaded. ' + species + ' ' + gene + ' skipped.', 'alert-warning')

    return id_list, gene_success

def fetch(id_list, gene_success):
    """
    Downloads the FASTA sequences from NCBI Nucleotide and RefSeqGene and concatenates them into one document.
    Determines size of flanking regions, just in case they are not 5kb upstream and 2kb downstream as by default.

    :param
        id_list:        (List, string) a list of unique IDs corresponding to the gene names
        gene_success    (List, string) a list of gene names that were successfully searched

    :return:
        document:       (FASTA File Object) a FASTA file containing the sequences of all the genes
        ref_list:       (List, string) a list of reference numbers for genes (can be used to search NCBI Gene database)
        flanking:       (List, List, string) a list of [upstream, downstream] flanking sizes for each gene
    """

    local = []

    try:
        # Fetch the FASTA file
        document = Entrez.efetch(db='nucleotide', id=id_list, rettype='fasta', retmode='text')
        # Fetch the full record
        detail = Entrez.efetch(db='nucleotide', id=id_list, rettype='native', retmode='xml')
        doc = minidom.parseString(detail.read())
        # Grab the gene IDs, e.g. NG_007392
        ref_list = doc.getElementsByTagName('GBSeq_locus').firstChild.nodeValue
        # Grab the lengths of the entries
        length = doc.getElementsByTagName('GBSeq_length').firstChild.nodeValue
        # Get the local coordinates of the gene
        for i, gene in doc.getElementsByTagName('GBSeq'):
            for element in gene.getElementsByTagName('GBFeature_key'):
                if element.firstChild.nodeValue == 'gene' and gene_success[i] \
                        in element.parentNode.getElementsByTagName('GBQualifier_value').firstChild.nodeValue:
                    # Coordinates usually come in this format '####..####' or 'complement(join(####..####,####..####))'
                    # But since genes are not split up, it should only come in this format '####..####'
                    local.append(element.nextSibling.firstChild.nodeValue.split('..'))
                    break

        # TODO: figure out how to get a list of values with getElementsByTagName
        print ref_list

    except:
        flash('Something went wrong when NCBI was fetching the FASTA file', 'alert')
        return render_template('fasta/form.html')

    return document, ref_list
