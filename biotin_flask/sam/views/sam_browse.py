# -*- coding: utf-8 -*-
"""Browse a SAM file."""

from collections import OrderedDict
import os

from flask import escape, flash, session, redirect, render_template, request
from flask import url_for

from biotin_flask import app
from biotin_flask.models.disk_io import basename, get_firstfile, list_dir
from biotin_flask.models.fasta import Fasta
from biotin_flask.models.pysam_ext import get_all_genes, get_refbase, load_sams
from biotin_flask.models.stream_io import csv_response, zip_csv_response

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.2'
__status__ = 'Production'


@app.route('/sam/browse', methods=['GET', 'POST'])
def sam_browse():
    """Render the sam browser.

    Methods:
        GET: Render the SAM browser page with form.
        POST: Process the form upload display browser results.

    """
    # Redirect user if there are no bam files uploaded
    if not session['id']:
        return redirect(url_for('sam_upload'), code=302)
    bampaths = list_dir('sam', session['id'], ('.bam'))
    fastapath = get_firstfile('sam', session['id'], ('.fa', '.fasta'))
    if len(bampaths) == 0:
        return redirect(url_for('sam_upload'), code=302)

    # Load samfiles and get list of genes in samfile
    bamfiles = load_sams(bampaths, 'rb')
    gene_list = get_all_genes(bamfiles)

    # Load Fasta file and validate it
    badFasta = True
    if fastapath:
        try:
            fasta = Fasta(fastapath)
            if fasta.validate(gene_list):
                badFasta = False
                gene_list = sorted(fasta.genes.keys())
        except:
            badFasta = True

    # Render an empty form for a GET request
    if request.method == 'GET':

        # Construct table where (rows, columns) = (genes, samples)
        gene_rows = []
        for gene in gene_list:
            row = []
            for bam in bamfiles:
                if gene in bam.references:
                    row.append(True)
                else:
                    row.append(False)
            gene_rows.append(row)

        # Render an empty form for a GET request
        return render_template(
            'sam/browse.html',
            session_id=escape(session['id']),
            samfiles=basename(bampaths, session['id']),
            fastafile=basename(fastapath, session['id']),
            badFasta=badFasta,
            genes=gene_list,
            rows=gene_rows
        )

    # Process form POST
    if request.method == 'POST':

        # Load common variables
        form = {
            'mode': request.form.get('mode'),
            'gene': request.form.get('gene'),
            'start': request.form.get('start'),
            'end': request.form.get('end'),
            'custom': request.form.get('custom'),
            'output': request.form.get('output'),
        }

        # Validate and parse form
        try:
            parsed_form = _validate_and_parse(form)
        except (KeyError, ValueError), message:
            flash(str(message), 'alert-warning')
            return redirect(url_for('sam_browse'), code=302)

        # Begin main loop
        i = 0
        error_count = 0
        results = []
        for bam in bamfiles:

            # Get coordinate positions of the selected gene as a list
            try:
                positions = _get_positions(bam, parsed_form)
            except ValueError, message:
                error_count = error_count + 1
                continue

            # Make the universal position container
            positionContainer = OrderedDict({parsed_form['gene']: positions})

            # Get the read bases at the requested coordinated positions
            readDict = _bam_fetch_positions(bam, positionContainer)

            # Handle the fasta file if present
            if badFasta:
                fasta_seq = None
            else:
                fasta_seq = _fasta_fetch_positions(fasta, positionContainer)

            results.append({
                'sample': basename(bampaths, session['id'])[i],
                'positions': _flatten_position_container(positionContainer),
                'reads': readDict,
                'fasta': fasta_seq
            })
            i = i + 1

        # Abort if all the bamfiles errored
        if len(bamfiles) == error_count:
            flash('Coordinates not found in sample(s)', 'alert-warning')
            return redirect(url_for('sam_browse'), code=302)

        # Determine output method - HTML
        if parsed_form['output'] == 'html':
            genes, positions = _split_and_add1(results[0]['positions'])
            return render_template(
                'sam/browse_results.html',
                sample=results[0]['sample'],
                gene=genes[0],
                positions=positions,
                reads=results[0]['reads'],
                fasta=results[0]['fasta']
            )

        # Output method - CSV
        if len(bamfiles) == 1:  # Single CSV file, send back a .csv file
            data = _format_csv_file(results[0])
            filename = os.path.splitext(results[0]['sample'])[0] + '.csv'
            return csv_response('result_' + filename, data)

        else:  # Multiple CSV files, send back a .zip file
            data = {}
            for result in results:
                csvdata = _format_csv_file(result)
                filename = os.path.splitext(result['sample'])[0] + '.csv'
                data[filename] = csvdata
            return zip_csv_response('results.zip', data)


def _validate_and_parse(form):
    """Validate and parse forms on POST to sam_browse.

    Note: We convert all inputs to zero-based indexing, because it is pysam's
    convention to use zero-based indexing. Therefore, for consistency, all of
    our code will also utilize zero-based indexing, and we assume the end-user
    is typing everything in using 1-based indexing.

    Parameters:
        form (dict): A dictionary of form inputs.

    Required Form Keys:
        mode (str): One of ('range', 'custom').
        gene (str): Gene name.
        start (str): Start position (1-based).
        end (str): End position (1-based).
        custom (str): A list of custom positions (1-based).
        output (str): Output mode.

    Returns:
        A parsed dictionary of form inputs.

    """
    newform = form
    if not (form['gene'] and form['mode'] and form['output']):
        raise KeyError('Required field is missing.')

    if form['mode'] == 'range':
        if not (form['start'] and form['end']):
            raise KeyError('Required field is missing.')
        newform['start'] = int(form['start']) - 1
        newform['end'] = int(form['end']) - 1
        newform['custom'] = None
        if newform['start'] < 0 or newform['end'] < 0:
            raise ValueError('Coordinate position cannot be negative.')
        if newform['end'] <= newform['start']:
            raise ValueError('`End` cannot be less than `Start`.')
        return newform

    if form['mode'] == 'custom':
        if not form['custom']:
            raise KeyError('Required field is missing.')
        newform['start'] = None
        newform['end'] = None
        try:
            custom = map(int, form['custom'].splitlines())
            custom = [x - 1 for x in custom]
            newform['custom'] = sorted(custom)
        except:
            raise ValueError('Error in custom coordinates provided.')
        return newform


def _get_positions(bam, form):
    """Construct a list of valid coordinate positions for interrogating reads.

    Note 1:
        If start/end coordinates are not present, then we will utilize the
        custom list to formulate a valid list.

    Note 2:
        Pysam AlignmentFile.pileup() uses an inclusive start and and exclusive
        end. In set notation, this is: [start, end). However, we expect our
        end users to enter coordinates inclusively, like: [start, end]. As a
        result we (+1) before running AlignmentFile.pileup().

    Parameters:
        bam (AlignmentFile object): The bam file to check.
        form (dict): A dictionary of form inputs.

    Required Form Keys:
        gene (str): Gene name.
        start (int): Start position (0-based).
        end (str): End position (0-based).
        custom (str): A list of custom positions (0-based).

    Returns:
        A sorted list of valid integer coordinate positions.

    """
    positions = []
    gene, start, end = (form['gene'], form['start'], form['end'])

    # If Start and End is not set - MODE: CUSTOM
    if start is None or end is None:
        start = min(form['custom'])
        end = max(form['custom'])

    try:
        for col in bam.pileup(gene, start, end + 1, truncate=True):
            if form['custom']:  # MODE: CUSTOM
                if col.reference_pos in form['custom']:
                    positions.append(col.reference_pos)
            else:  # MODE: RANGE
                positions.append(col.reference_pos)
        if len(positions) == 0:
            raise ValueError()
    except ValueError:
        raise ValueError('Requested positions not found in this sample.')

    return sorted(positions)


def _bam_fetch_positions(bam, position_container):
    """Construct a table of bases at the requested posistions.

    Note 1:
        We require an OrderedDict in the input because the positions across
        different genes are concatenated horizontally in the output. Therefore,
        order must be maintained in order to know which gene each section
        refers to.

    Note 2:
        Pysam AlignmentFile.fetch() uses inclusive `start` and exclsuive `end`,
        so we (+1) to `end` in order to get an inclusive range. See also:
        method _get_positions (note 2).

    Parameters:
        bam (AlignmentFile object): The bam file to check.
        position_container (OrderedDict): An ordered dict of position lists.

    Returns:
        An OrderedDict of reads of a list of string bases.

    Example_Position_Container_Input = {
        'geneA' = [ 1, 2],
        'geneB' = [ 1, 2]
    }

    Example_Output = {
        'read1' = {
            'strand': False,
            'bases': [ 'A', 'C', '', '']
        },
        'read2' = {
            'strand': True,
            'bases': [ '', '', 'T', 'G']
        }
    }

    """
    # Get the total number of positions
    total_positions = 0
    for positions in position_container.values():
        total_positions = total_positions + len(positions)

    # Begin loop
    reads = OrderedDict()
    buffer_front = 0
    buffer_back = total_positions
    for gene, positions in position_container.items():
        buffer_back = buffer_back - len(positions)
        start, end = min(positions), max(positions)
        for read in bam.fetch(gene, start, end + 1):
            reads[read.query_name] = {
                'strand': read.is_reverse,
                'bases': [''] * buffer_front
            }
            for pos in positions:
                reads[read.query_name]['bases'].append(get_refbase(read, pos))
            if buffer_back:
                reads[read.query_name]['bases'].append([''] * buffer_back)
        buffer_front = len(positions)
    return reads


def _fasta_fetch_positions(fasta, position_container):
    """Fetch requested bases from fasta file.

    Parameters:
        fasta (Fasta object): Fasta file to search.
        position_container (OrderedDict): An ordered dict of position lists.

    Returns:
        A string list of bases at the flattened position_container coordinates.

    """
    result = []
    for gene, positions in position_container.items():
        for pos in positions:
            result.append(fasta.genes[gene][pos].upper())
    return result


def _flatten_position_container(position_container):
    """Flattens the position container structure to a list of positions.

    Arguments:
        position_container (OrderedDict): An ordered dict of position lists.

        Example_Position_Container_Input = {
            'geneA' = [ 1, 2],
            'geneB' = [ 1, 2]
        }

    Returns:
        A string list of 0-based positions concatenated with the gene name.

        Example_Output = ['GeneA:1', 'GeneA:2', 'GeneB:1', 'GeneB:2']

    """
    result = []
    for gene, positions in position_container.items():
        for pos in positions:
            result.append(gene + ':' + str(pos))
    return result


def _split_and_add1(geneposition_list):
    """Split a list of gene:positions and add 1 to prepare for printing.

    Note:
        We add 1 to print the final coordinates because pysam uses 0-based
        coordinates, whereas the end-user expects 1-based coordinates.

    Arguments:
        geneposition_list (list): A list of string `gene:position` values.

        Example_Input = ['GeneA:0', 'GeneA:1', 'GeneB:0', 'GeneB:1']

    Returns:
        A string list of genes and an integer list of positions (+1).

        Example_Output = [
            ['GeneA', 'GeneA', 'GeneB', 'GeneB'],
            [1, 2, 1, 2]
        ]

    """
    result = [[], []]
    for item in geneposition_list:
        gene, pos = item.split(':')
        result[0].append(gene)
        result[1].append(int(pos) + 1)
    return result


def _format_csv_file(data):
    """Format data into csv file ready to be print.

    Arguments:
        data (dict): Data to format. Note that positions are 0-based.

        data_example = {
            'sample': 'sample1.bam',
            'positions': ['GeneA:0', 'GeneA:1', 'GeneB:0', 'GeneB:1'],
            'reads': {
                'read1' = {
                    'strand': False,
                    'bases': [ 'A', 'C', '', '']
                },
                'read2' = {
                    'strand': True,
                    'bases': [ '', '', 'T', 'G']
                }
            },
            'fasta': ['A', 'C', 'T', 'G']
        }

    Returns:
        A list of string lists ready to be print. See example.

        return_example = [
            ['', 'Gene:', 'GeneA', 'GeneA', 'GeneB', 'GeneB'],
            ['Read Strand', 'Read ID', 0, 1, 0, 1],
            ['', 'Fasta Sequence', 'A', 'C', 'T', 'G'],
            ['+', 'read1', 'A', 'C', '', ''],
            ['+', 'read1', '', '', 'T', 'G']
        ]

    """
    genes, positions = _split_and_add1(data['positions'])
    result = []
    result.append(['', 'Gene:'] + genes)
    result.append(['Read Strand', 'Read ID'] + positions)
    if data['fasta']:
        result.append(['', 'Fasta Sequence'] + data['fasta'])
    for read, d in data['reads'].items():
        if d['strand']:
            strand = '-'
        else:
            strand = '+'
        result.append([strand, read] + d['bases'])
    return result
