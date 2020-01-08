"""
This module takes a list of gene names as input,
(1) searches NCBI Nucleotide for the IDs of genes and their sequences
(2) fetches the FASTA files from the nucleotide IDs and gene IDs
and outputs a FASTA file with every gene and its sequence.
Default sequence includes 5kb upstream and 2kb downstream flanking sequence.
"""

import os, io, csv, StringIO, zipfile, re
from flask import render_template, request, flash, send_file
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
        species:    (String) the species
        gene:       (String) a list of gene names
        name:       (String) [optional] name of file output, 'references_genes.fasta' is default

    :return:
        result:     (File Object) a FASTA file

    """

    # Render an empty form for GET request
    if request.method == 'GET':
        return render_template('fasta/form.html')

    # Otherwise validate the form on a POST request and process uploaded files
    species = request.form.get('species')
    gene = request.form.get('gene').splitlines()
    if not species or not gene:
        flash('Fill out required fields', 'alert-warning')
        return render_template('fasta/form.html')
    name = request.form.get('name')
    if not name:
        name = 'reference_genes.fasta'

    # Obtain NCBI Nucleotide IDs and list of genes that were successfully searched
    id_list, gene_success = uid(species, gene)
    if not len(gene_success) == len(gene):
        return render_template('fasta/form.html')

    # Obtain FASTA document, list of NCBI Gene IDs, size of flanking sequences, list of genes that
    # Were successfully searched
    document, gene_id, flanking, gene_success = fetch(id_list, gene_success)
    if not len(gene_success) == len(gene):
        return render_template('fasta/form.html')
    print gene, gene_success, flanking, gene_id, document
    read_fasta = document.read()

    # Get genomic coordinates of the FASTA file
    coord = coordinates(gene_id, flanking)
    if not coord:
        flash('Error in looking up genomic coordinates. Reduce the number of genes and try again.', 'alert-warning')
        return render_template('fasta/form.html')

    # Form headers of the new FASTA file, e.g. "FOXP3.chrX:49245435-49266931"
    headers = [gene_name + '.' + coordinate for gene_name, coordinate in zip(gene_success, coord)]
    print '', gene_id, flanking, gene_success, coord, headers

    # Modify the FASTA file and return the result
    with open(os.path.join(app.config['UPLOAD_FOLDER'], 'fasta_results.fasta'), 'w') as dest:
        counter = 0 # for counting genes
        for i, line in enumerate(read_fasta.split('\n')):
            if line:
                if line[0] == '>':
                    dest.writelines('>' + headers[counter] + '\n')
                    counter += 1
                else:
                    dest.writelines(line + '\n')

    document.close()

    # Send File
    return send_file(os.path.join(app.config['UPLOAD_FOLDER'], 'fasta_results.fasta'),
                     attachment_filename=name,
                     as_attachment=True)

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
                flash('NCBI ESearch did not return any results for ' + species + ' ' + gene, 'alert-warning')

        except:
            flash('NCBI ESearch was unable to search for ' + species + ' ' + gene +
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
        gene_success:   (List, string) a list of genes that passed EFetch
    """

    flanking = []
    ref_list = []

    # Fetch the full record
    detail = Entrez.efetch(db='nucleotide', id=id_list, rettype='gb', retmode='xml', api_key=API_KEY)
    doc = minidom.parseString(detail.read())
    # Grab the lengths of the entries
    length = list(map(int,getListOfValues(doc.getElementsByTagName('GBSeq_length'))))
    # Get the local coordinates of the gene
    i = 0 # Set up counter because cannot enumerate over Node List :(
    for gene in doc.getElementsByTagName('GBSeq'):

        """
        # This was too error-prone
        # Check to see if gene_success corresponds to the gene name in the document
        gene_name = gene.getElementsByTagName('GBSeq_definition')[0].firstChild.nodeValue
        gene_name = re.search(r'\((.*?)\)', gene_name).group(1)
        if not gene_name == gene_success[i]:
            flash('Warning: NCBI did not find ' + gene_success[i] + ', but found '
                  + gene_name + ' instead. Replacing ' + gene_success[i]
                  + ' with ' + gene_name, 'alert-warning')
        """

        gene_id = 0
        for element in gene.getElementsByTagName('GBFeature_key'):

            nodeList = getListOfValues(element.parentNode.getElementsByTagName('GBQualifier_value'))

            # If the feature type is a gene
            # And if one of the sibling nodes "GBQualifier_value" contains the name of the gene
            if element.firstChild.nodeValue == 'gene' and gene_success[i] in nodeList:
                print gene_success[i]
                # Get the Gene ID
                for nodeValue in nodeList:
                    if 'GeneID:' in nodeValue:
                        gene_id = 1
                        ref_list.append(nodeValue.split(':')[1])
                        break

                # assuming gene coordinates are the first thing
                start = int(element.parentNode.getElementsByTagName('GBInterval_from')[0].firstChild.nodeValue)
                end = int(element.parentNode.getElementsByTagName('GBInterval_to')[0].firstChild.nodeValue)
                flanking.append((start-1,length[i]-end))
                break

        if not gene_id:
            flash('NCBI EFetch did not find ' + gene_success[i]
                  + '.', 'alert-warning')
            # Delete the gene from the ID list, so you don't end up fetching the FASTA file for it.
            id_list.pop(i)
            gene_success.pop(i)
        else:
            i += 1

    # Fetch the FASTA file
    document = Entrez.efetch(db='nucleotide', id=id_list, rettype='fasta', retmode='text', api_key=API_KEY)

    return document, ref_list, flanking, gene_success


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

    # Try reducing the size of the request until NCBI gives you data back
    total = len(ref_list) # Last index
    a = 0
    x = 0
    while total - a > 0:
        print total, a, x, [a, total-x]
        detail = Entrez.efetch(db='gene', id=ref_list[a:total-x], retmode='xml', api_key=API_KEY)
        doc = minidom.parseString(detail.read())
        if not doc.getElementsByTagName('Entrezgene_locus'):
            if total-x-a > 1:
                x += 1
            continue
        else:
            a = total - x
            x = 0

        # Can't enumerate over nodeList so setting counter here
        i = 0
        for gene in doc.getElementsByTagName('Entrezgene_locus'):
            # Assuming that the Reference Coordinates are in the first <Gene-commentary>
            # The chromosome number should be in "Chromosome # Reference GRCh38.p13 Primary Assembly"
            chr = gene.getElementsByTagName('Gene-commentary_label')[0].firstChild.nodeValue.split(' ')[1]
            start = int(gene.getElementsByTagName('Seq-interval_from')[0].firstChild.nodeValue) - flanking[i][0]
            end = int(gene.getElementsByTagName('Seq-interval_to')[0].firstChild.nodeValue) + flanking[i][1]
            if chr and start and end:
                coordinates.append('chr'+chr+':'+str(start)+'-'+str(end))
            else:
                coordinates.append('')
            i += 1

    return coordinates
