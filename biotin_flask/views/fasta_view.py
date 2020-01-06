"""
This module takes a list of gene names as input,
(1) searches NCBI Nucleotide for the IDs of genes and their sequences
(2) fetches the FASTA files from the nucleotide IDs and gene IDs
and outputs a FASTA file with every gene and its sequence.
Default sequence includes 5kb upstream and 2kb downstream flanking sequence.
"""

import os, io, csv, StringIO, zipfile, re
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
API_KEY = 'c9208f5fd7eba42663b68a994027be7bb408'

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

    # Obtain FASTA file from NCBI, list of NCBI Gene IDs and size of flanking sequences
    id_list, gene_success = uid(species, gene)
    document, gene_id, flanking = fetch(id_list, gene_success)

    # Get genomic coordinates of the FASTA file
    coord = coordinates(gene_id, flanking)

    # Modify the FASTA file
    # TODO: Confirm that coordinates are correct and then write code to modify FASTA file

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
            response = Entrez.esearch(db='nucleotide', term=species+'[Organism] '+gene+'[GENE] RefSeqGene[KYWD]',
                                      api_key=API_KEY)
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
        id_list:        (List, string) a list of NCBI unique IDs corresponding to the gene names
        gene_success    (List, string) a list of gene names that were successfully searched

    :return:
        document:       (FASTA File Object) a FASTA file containing the sequences of all the genes
        ref_list:       (List, string) a list of reference numbers for genes (can be used to search NCBI Gene database)
        flanking:       (List, List, string) a list of [upstream, downstream] flanking sizes for each gene
    """

    flanking = []
    ref_list = []

    # Fetch the FASTA file
    document = Entrez.efetch(db='nucleotide', id=id_list, rettype='fasta', retmode='text', api_key=API_KEY)
    # Fetch the full record
    detail = Entrez.efetch(db='nucleotide', id=id_list, rettype='gb', retmode='xml', api_key=API_KEY)
    doc = minidom.parseString(detail.read())
    # Grab the lengths of the entries
    length = list(map(int,getListOfValues(doc.getElementsByTagName('GBSeq_length'))))
    # Get the local coordinates of the gene
    i = 0 # Set up counter because cannot enumerate over Node List :(
    for gene in doc.getElementsByTagName('GBSeq'):
        for element in gene.getElementsByTagName('GBFeature_key'):

            nodeList = getListOfValues(element.parentNode.getElementsByTagName('GBQualifier_value'))

            # If the feature type is a gene
            # And if one of the sibling nodes "GBQualifier_value" contains the name of the gene
            if element.firstChild.nodeValue == 'gene' and gene_success[i] in nodeList:

                # Get the Gene ID
                for nodeValue in nodeList:
                    if 'GeneID:' in nodeValue:
                        ref_list.append(nodeValue.split(':')[1])
                        break

                # Coordinates usually come in this format '####..####' or 'join(####..####,####..####)'
                # But since genes are not split up, it should only come in this format '####..####'
                # Or the complement 'complement(####..####)'
                coordinates = element.parentNode.getElementsByTagName('GBFeature_location')[0].\
                    firstChild.nodeValue
                if '(' in coordinates:
                    # RegEx to search for '(####..####)'
                    start, end = list(map(int,re.search(r'\(\d*\.\.\d*\)', coordinates).group(1).split('..')))
                else:
                    start,end = list(map(int,coordinates.split('..')))
                flanking.append((start-1,length[i]-end))
                break
        i += 1

    return document, ref_list, flanking

def coordinates(ref_list, flanking):
    """
    Fetches list of genomic coordinates of the FASTA files given NCBI gene IDs and size of flanking sequences

    :param
        ref_list:       (List, String) List of NCBI Gene ID Strings to Search NCBI Gene
        flanking:       (List, Tuple, Int) List of size of (left, right) flanking sequences on the FASTA files

    :return:
        coordinates:    (List, Tuple, String) List of genomic coordinates for each of the NCBI Gene IDs
    """

    coordinates = []
    detail = Entrez.efetch(db='gene', id=ref_list, retmode='xml', api_key=API_KEY)
    doc = minidom.parseString(detail.read())
    i = 0 # Can't enumerate over nodeList so setting counter here
    for gene in doc.getElementsByTagName('Entrezgene_locus'):
        # Assuming that the Reference Coordinates are in the first <Gene-commentary>
        # The chromosome number should be in "Chromosome # Reference GRCh38.p13 Primary Assembly"
        chr = gene.getElementsByTagName('Gene-commentary_label')[0].firstChild.nodeValue.split(' ')[1]
        start = int(gene.getElementsByTagName('Seq-interval_from')[0].firstChild.nodeValue) - flanking[i][0]
        end = int(gene.getElementsByTagName('Seq-interval_to')[1].firstChild.nodeValue) + flanking[i][1]
        if chr and start and end:
            coordinates.append('chr'+chr+':'+str(start)+'-'+str(end))
        else:
            coordinates.append('')
        i += 1

    return coordinates

def getListOfValues(nodeList):
    """
    Takes a DOM Node List Object and returns a list of their Node values.

    :param
        nodeList:       (Node List) a DOM Node List object
    :return:
        valueList:      (List, String) a list of values inside the Nodes, usually but not necessarily Strings
    """

    valueList = []
    for node in nodeList:
        valueList.append(node.firstChild.nodeValue)

    return valueList