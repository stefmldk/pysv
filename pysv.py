"""
This script identifies large, structural variants in a bam file containing mapped long-reads. The script can
take a single bam file as input or a folder with one or several bam files. It outputs a "structural_variants" folder
to the folder containing the bam file(s), so write access to this folder is required. If a single bam file is analyzed,
it is possible to provide a sample name, which will be added as a suffix to the output. This option is not available
if multiple bam files are analyzed via the recursive flag. However, if sample names are written in the headers of the
bam-files (requires a non-default setting in minimap2), they will be used as suffix for each bam file's output.
"""

import json
import os
import itertools
import sys
import re
import multiprocessing
import subprocess
import time
import pandas
import zipfile
from bisect import bisect
import random


# Provide user-friendly feedback in case non-standard packages are not installed

def no_trace_error(error_msg):
    """
    Alternative to raise error to avoid the trace and produce a cleaner feedback to the user. Writes to stderr and
    exits with nonzero exit code, so that a pipeline will detect the error and stop.
    :param error_msg:
    :return:
    """
    print('\n' + error_msg, file=sys.stderr)
    exit(-1)


try:
    from bs4 import BeautifulSoup as beautiful_soup
except ModuleNotFoundError:
    no_trace_error('It appears that bs4 is not installed. Please install it and run this script again.')

try:
    import plotly.express as plotly_express
except ModuleNotFoundError:
    no_trace_error('It appears that plotly is not installed. Please install it and run this script again.')

try:
    import pysam
except ModuleNotFoundError:
    no_trace_error('It appears that pysam is not installed. Please install it (available for Linux-based platforms) and run this script again.')

try:
    from scipy import stats
    from scipy.ndimage import median
except ModuleNotFoundError:
    no_trace_error('It appears that scipy is not installed. Please install it and run this script again.')


__author__ = "Steffen Møller Bøttger"
__date__ = "03-06-2025"
__copyright__ = "Copyright 2025, Steffen Møller Bøttger"
__email__ = "Steffen.Moller.Bottger@rsyd.dk"
__license__ = "MIT"
__version__ = '1.8.6'


def zipdir(dir_path, zip_container):
    """
    Inspired by: https://stackoverflow.com/questions/1855095/how-to-create-a-zip-archive-of-a-directory
    :param dir_path:
    :param zip_container:
    :return:
    """
    for path, dirs, files in os.walk(dir_path):
        for _file in files:
            zip_container.write(os.path.join(path, _file))


def iterate_bamfiles(folder=None, recursive=False):
    """
    Traverses the given directory or pwd recursively and returns the path to the first encountered bam file unless
    recursive == True in which case traversal continues to all sub folders and all encountered bam files are returned.

    :param folder:          Folder to start the traversal from
    :param recursive:       If True, traversal continues after first bam file and all encountered bam files are returned
    :return:
    """
    directory = folder or os.getcwd()
    for path, dirs, files in os.walk(directory):
        for file_name in files:

            if file_name.endswith('.bam'):
                yield os.path.join(path, file_name)
                # bam_files.append(os.path.join(path, file_name))
                if not recursive:
                    return
    return False  # Returning False lastly so that we can check for existence of files via Truthy check


def get_bamfiles(fs_node=None, recursive=False):
    """
    Using this because os.walk() does not throw an error for nonexistent folder and because a conditional cannot be
    used in a generator function to print feedback in case the folder is nonexistent.
    :param fs_node:
    :param recursive:
    :return:
    """
    if os.path.exists(fs_node):
        if os.path.isfile(fs_node):
            return [fs_node]
        else:
            return iterate_bamfiles(fs_node, recursive)
    else:
        print('\nIt seems that the provided path does not exist:', fs_node)
        no_trace_error('Please provide a working path to a bam file or a directory containing bam files')



def get_clippings(read):
    """
    Returns a tupple with start/left and end/right clipping of bases derived from the CIGAR string. Both soft- and
    hard-clippings are included.

    This approach assumes that an alignment only has clippings at start and end, which is in accordance with the SAM/BAM
    Format Specification.
    """
    cigar = read.cigarstring
    start_clipping = 0
    end_clipping = 0

    # Start clipping
    for i in range(20):
        if not cigar[i].isnumeric():
            break
    if cigar[i] == 'S' or cigar[i] == 'H':
        start_clipping = int(cigar[:i])

    # End clipping
    if cigar[-1] == 'S' or cigar[-1] == 'H':
        for i in range(2, 20):
            if not cigar[-i].isnumeric():
                break
        end_clipping = int(cigar[-(i - 1):-1])

    return start_clipping, end_clipping


def get_sorted_alphanum_strings_ascending(names_iterable, index=None):
    """
    Sorts a list of alphanumeric strings so that contained numbers will appear in a numerical logic, ascending order. If
    no numbers are part of the strings, they will be sorted lexically. If an index is provided, a list of lists can be
    sorted relative to the values found at that index.
    """
    alphanum_key = lambda key: [int(string) if string.isdigit() else string for string in re.split('([0-9]+)', key[index] if index is not None else key)]
    return sorted(names_iterable, key=alphanum_key)


def get_min_alphanum(names_iterable):
    """
    Returns the lexically smallest string in the provided list of alphanumeric strings
    """
    return get_sorted_alphanum_strings_ascending(names_iterable)[0]


class BamData:
    """
    Class for storing and retrieving bam file data
    """

    def get_contig_lengths(self):
        header = self.bam_obj.header.to_dict()
        contig_data = header['SQ']
        contig_data =  sorted(contig_data, key=lambda x:x['LN'], reverse=True)
        if len(contig_data) >= self.nuclear_genome_contigs:
            contig_data = contig_data[:self.nuclear_genome_contigs]
        bam_contigs = {}
        for element in contig_data:
            bam_contigs[element['SN']] = element['LN']
        return bam_contigs


    def get_pointer_contig(self):
        return self.contig_names[self.pointer['contig']]


    def get_next_contig_chunk(self, chunk_length=None):
        """
        Gets next contig chunk relative to the current pointer and updates the pointer to the base after the chunk
        unless chunk_length == None in which case the rest of the contig_lengths from current pointer is returned.

        :param (int) chunk_length:  If none is given, a chunk corresponding to the rest of the contig_lengths from the current pointer is returned
        :return:
        """
        contig_chunk = []

        if not chunk_length:
            contig_chunk.append((
                self.get_pointer_contig(),
                self.pointer['location'],
                self.contig_lengths[self.get_pointer_contig()]
            ))
            if self.pointer['contig'] < len(self.contig_lengths) - 1:
                for contig_pointer in range(self.pointer['contig'] + 1, len(self.contig_lengths)):
                    contig_chunk.append((
                        self.contig_names[contig_pointer],
                        0,
                        self.contig_lengths[self.contig_names[contig_pointer]]
                    ))
        else:
            remainder = chunk_length
            while remainder:
                rest_of_current_contig = self.contig_lengths[self.get_pointer_contig()] - self.pointer['location']
                if rest_of_current_contig <= remainder:
                    contig_chunk.append((  # Append rest of contig
                        self.get_pointer_contig(),
                        self.pointer['location'],
                        self.contig_lengths[self.get_pointer_contig()]
                    ))
                    remainder -= rest_of_current_contig
                    self.pointer['contig'] += 1
                    self.pointer['location'] = 0
                else:
                    new_pointer_location = self.pointer['location'] + remainder
                    contig_chunk.append((
                        self.get_pointer_contig(),
                        self.pointer['location'],
                        new_pointer_location
                    ))
                    self.pointer['location'] = new_pointer_location
                    remainder = 0

        return contig_chunk


    def get_reference_subdivisions(self, number_of_subdivisions):
        """
        Given a reference dictionary of contig names and their lengths, this function subdivides it into chunks of even
        base lengths according to the given number of subdivisions. A chunk may consist of multiple chromosome parts that
        add up to the chunk length. The last chunk may contain any remainder and so may be a little longer than the rest.
        :param contig_lengths:
        :param number_of_subdivisions:
        :return:
        """
        chunk_length = self.genome_length // number_of_subdivisions

        chunks = []
        for i in range(number_of_subdivisions - 1):
            chunks.append(self.get_next_contig_chunk(chunk_length))
        chunks.append(self.get_next_contig_chunk())  # Gets the last chunk including remainder
        return chunks


    def get_region_depth(self, contig, start, stop, num_samples):
        """
        Gets depth in the region defined by start and stop based on num_samples random samples.

        :param contig:
        :param start:
        :param stop:
        :param num_samples:
        :return:
        """
        start = max(0, start)
        stop = min(stop, self.contig_lengths[contig]) - 1
        depth_measurements = []
        for i in range(num_samples):
            rand_pos = random.randrange(start, stop)
            depth_measurements.append(self.bam_obj.count(contig, rand_pos, rand_pos + 1))
        return median(depth_measurements)


    def get_region_x_depth(self, contig, start, stop, num_samples=10):
        """
        Returns depth as X-count, where 2X is the median depth for the BAM file
        :param contig:
        :param start:
        :param stop:
        :param num_samples:
        :return:
        """
        region_depth = self.get_region_depth(contig, start, stop, num_samples)
        return round(region_depth / self.one_x) if region_depth else 0


    def get_translocation_x_depths(self, breakpoint_1, breakpoint_2, offset_length=5000000):
        """
        Gets x_depths relative to translocation breakpoints

        :param breakpoint_1:        Breakpoint defined as (chr, breakpoint_1)
        :param breakpoint_2:        Breakpoint defined as (chr, breakpoint_1)
        :param offset_length:       Offset from breakpoints, i.e. it defines upstream/downstream regions
        :return:
        """
        if (
                breakpoint_1[1] < 10 or
                breakpoint_2[1] < 10 or
                10 > abs(self.contig_lengths[breakpoint_1[0]] - breakpoint_1[1]) or
                10 > abs(self.contig_lengths[breakpoint_2[0]] - breakpoint_2[1])
        ):
            return None

        return [
            self.get_region_x_depth(breakpoint_1[0], breakpoint_1[1] - offset_length, breakpoint_1[1]),
            self.get_region_x_depth(breakpoint_1[0], breakpoint_1[1], breakpoint_1[1] + offset_length),
            self.get_region_x_depth(breakpoint_2[0], breakpoint_2[1] - offset_length, breakpoint_2[1]),
            self.get_region_x_depth(breakpoint_2[0], breakpoint_2[1], breakpoint_2[1] + offset_length)
        ]


    def get_translocation_with_deletion_x_depths(self, breakpoint_1, del_start, del_stop, offset_length=5000000):
        """
        Gets x_depths relative to breakpoint_1 and start/stop of putative deletion
        :param breakpoint_1:        Breakpoint defined as (chr, breakpoint_1)
        :param del_start:           Start coordinate defined as (chr, start)
        :param del_stop:            Stop coordinate defined as (chr, stop)
        :param offset_length:       Offset from breakpoints, i.e. it defines upstream/downstream regions
        :return:
        """

        if (
                breakpoint_1[1] < 10 or
                10 > abs(self.contig_lengths[breakpoint_1[0]] - breakpoint_1[1]) or
                del_start[1] < 10 or
                10 > abs(self.contig_lengths[del_stop[0]] - del_stop[1])
        ):
            return None, None

        deletion_x_depth = self.get_region_x_depth(del_start[0], del_start[1], del_stop[1], num_samples=4)
        x_depths = [
            self.get_region_x_depth(breakpoint_1[0], breakpoint_1[1] - offset_length, breakpoint_1[1]),
            self.get_region_x_depth(breakpoint_1[0], breakpoint_1[1], breakpoint_1[1] + offset_length),
            self.get_region_x_depth(del_start[0], del_start[1] - offset_length, del_start[1]),
            self.get_region_x_depth(del_stop[0], del_stop[1], del_stop[1] + offset_length),
        ]
        return  deletion_x_depth, x_depths


    def get_translocation_with_double_deletion_x_depths(self, del_1_start, del_1_stop, del_2_start, del_2_stop, offset_length=5000000):
        """
        Gets x_depths relative to start/stop of putative translocation deletions
        :param del_1_start:         Start coordinate defined as (chr, start)
        :param del_1_stop:          Stop coordinate defined as (chr, stop)
        :param del_2_start:         Start coordinate defined as (chr, start)
        :param del_2_stop:          Stop coordinate defined as (chr, stop)
        :param offset_length:       Offset from breakpoints, i.e. it defines upstream/downstream regions
        :return:
        """
        if (
                del_1_start[1] < 10 or
                10 > abs(self.contig_lengths[del_1_stop[0]] - del_1_stop[1]) or
                del_2_start[1] < 10 or
                10 > abs(self.contig_lengths[del_2_stop[0]] - del_2_stop[1])
        ):
            return None, None, None


        del_1_x_depth = self.get_region_x_depth(del_1_start[0], del_1_start[1], del_1_stop[1], num_samples=4)
        del_2_x_depth = self.get_region_x_depth(del_2_start[0], del_2_start[1], del_2_stop[1], num_samples=4)
        x_depths = [
            self.get_region_x_depth(del_1_start[0], del_1_start[1] - offset_length, del_1_start[1]),
            self.get_region_x_depth(del_1_stop[0], del_1_stop[1], del_1_stop[1] + offset_length),
            self.get_region_x_depth(del_2_start[0], del_2_start[1] - offset_length, del_2_start[1]),
            self.get_region_x_depth(del_2_stop[0], del_2_stop[1], del_2_stop[1] + offset_length),
        ]
        return  del_1_x_depth, del_2_x_depth, x_depths


    def has_balanced_translocation_with_deletion_coverage_profile(self, breakpoint_1, del_start, del_stop):
        """
        Boolean. Returns True if X-coverage profile is 1X, 2X, 2X, 2X, 2X. For the deletion, only four samples is used
        as opposed to the default of 10 because it is typically way smaller than the offset_length. This is the only
        coverage assessment function with deletion where the deletion coverage is included. That is because unbalanced
        and malsegregation profiles present uniqueness by themselves without including the deletion.

        :param breakpoint_1:        Breakpoint defined as (chr, breakpoint_1)
        :param del_start:           Start coordinate defined as (chr, start)
        :param del_stop:            Stop coordinate defined as (chr, stop)
        :return:
        """

        deletion_x_depth, x_depths = self.get_translocation_with_deletion_x_depths(breakpoint_1, del_start, del_stop)

        if deletion_x_depth is None:
            return False

        x_depths.append(deletion_x_depth)

        print('Potential balanced translocation with deletion', x_depths)
        if deletion_x_depth > 1:
            print(f'Deletion depth is {deletion_x_depth}. Location is {del_start[0]}:{del_start[1]}-{del_stop[1]}')  # If del depth is > 1 it may indicate a low complexity region with multimapping
        return sorted(x_depths) == [1, 2, 2, 2, 2]
    
    
    def has_balanced_translocation_with_double_deletion_coverage_profile(self, del_1_start, del_1_stop, del_2_start, del_2_stop):
        """

         :param del_1_start:         Start coordinate defined as (chr, start)
        :param del_1_stop:          Stop coordinate defined as (chr, stop)
        :param del_2_start:         Start coordinate defined as (chr, start)
        :param del_2_stop:          Stop coordinate defined as (chr, stop)
        :return:
        """
        deletion_1_x_depth, deletion_2_x_depth,  x_depths = self.get_translocation_with_double_deletion_x_depths(del_1_start, del_1_stop, del_2_start, del_2_stop)

        if x_depths is None:
            return False

        x_depths += [deletion_1_x_depth, deletion_2_x_depth]

        print('Potential balanced translocation with deletion', x_depths)
        if deletion_1_x_depth > 1:
            print(
                f'Deletion depth is {deletion_1_x_depth}. Location is {del_1_start[0]}:{del_1_start[1]}-{del_1_stop[1]}')  # If del depth is > 1 it may indicate a low complexity region with multimapping
        if deletion_2_x_depth > 1:
            print(
                f'Deletion depth is {deletion_2_x_depth}. Location is {del_2_start[0]}:{del_2_start[1]}-{del_2_stop[1]}')  # If del depth is > 1 it may indicate a low complexity region with multimapping
        return sorted(x_depths) == [1, 1, 2, 2, 2, 2]


    def has_unbalanced_translocation_coverage_profile(self, breakpoint_1, breakpoint_2):
        """
        An unbalanced translocation will have region coverages corresponding to two times 2X, one time 1X and one time
        3X. There is oly one deviating chromosome, so only one breakpoint from this chromosome generates informative
        reads. Thus, info about potential deletions during translocation is lost.

        To protect against breakpoints called at or near chromosome ends, we return False if one of the breakpoints
        are within 10 bases from those.


        :param breakpoint_1:        Breakpoint defined as (chr, breakpoint_1)
        :param breakpoint_2:        Breakpoint defined as (chr, breakpoint_1)
        :return:
        """

        x_depths = self.get_translocation_x_depths(breakpoint_1, breakpoint_2)

        if x_depths is None:
            return False

        print('Potential UB translocation', x_depths)
        return sorted(x_depths) == [1, 2, 2, 3]


    def has_translocation_with_malsegregation_coverage_profile(self, breakpoint_1, breakpoint_2):
        """
        Unbalanced translocations with meiotic malsegregation (for instance t(11;22) in Emanuel syndrome) will have
        region coverages corresponding to two times 2X and two times 3X.

        To protect against breakpoints called at or near chromosome ends, we return False if one of the breakpoints
        are within 10 bases from those.


        :param breakpoint_1:        Breakpoint defined as (chr, breakpoint_1)
        :param breakpoint_2:        Breakpoint defined as (chr, breakpoint_1)
        :return:
        """

        x_depths = self.get_translocation_x_depths(breakpoint_1, breakpoint_2)
        if x_depths is None:
            return False

        return sorted(x_depths) == [2, 2, 3, 3]


    def has_translocation_with_malsegregation_and_deletion_coverage_profile(self, breakpoint_1, del_start, del_stop, offset_length=5000000):
        """
        Unbalanced translocations with meiotic malsegregation (for instance t(11;22) in Emanuel syndrome) will have
        region coverages corresponding to two times 2X and two times 3X.

        :param breakpoint_1:              Breakpoint defined as (chr, breakpoint_1)
        :param del_start:               Deletion start defined as (chr, del_start)
        :param del_stop:                Deletion start defined as (chr, del_stop)
        :param offset_length:           Offset from breakpoints, i.e. it defines upstream/downstream regions to find depths from
        :return:
        """

        deletion_x_depth, x_depths = self.get_translocation_with_deletion_x_depths(breakpoint_1, del_start, del_stop, offset_length=offset_length)

        if deletion_x_depth is None:
            return False

        print('Potential translocation with malsegregation and deletion', x_depths[:].append(deletion_x_depth))

        # We do not use the deletion depth because the malsegregation coverage profile in itself is quite unique
        return sorted(x_depths) == [2, 2, 3, 3]


    def has_translocation_with_malsegregation_and_double_deletion_coverage_profile(self, del_1_start, del_1_stop, del_2_start, del_2_stop):
        """
        Boolean
        Unbalanced translocations with meiotic malsegregation (for instance t(11;22) in Emanuel syndrome) will have
        region coverages corresponding to two times 2X and two times 3X.

        :param del_1_start:         Start coordinate defined as (chr, start)
        :param del_1_stop:          Stop coordinate defined as (chr, stop)
        :param del_2_start:         Start coordinate defined as (chr, start)
        :param del_2_stop:          Stop coordinate defined as (chr, stop)
        :param offset_length:       Offset from breakpoints, i.e. it defines upstream/downstream regions
        :return:
        """

        deletion_1_x_depth, deletion_2_x_depth,  x_depths = self.get_translocation_with_double_deletion_x_depths(del_1_start, del_1_stop, del_2_start, del_2_stop)

        if x_depths is None:
            return False

        print('Potential translocation with malsegregation and double deletion', x_depths)

        # We do not use the deletion depth because the malseg coverage profile in itself is quite unique
        return sorted(x_depths) == [2, 2, 3, 3]


    def chrom_number(self, chrom_name):
        """
        Gets chromosome number as integer. If chromosome has no numer (e.g. chrX), the given is returned
        :param chrom_name:
        :return:
        """
        try:
            return int(''.join([s for s in chrom_name if s.isdigit()]))
        except ValueError:  #  Will happen for non-digit chromosomes like X, Y, M etc.
            return chrom_name


    def get_depths(self, samples=50):
        """
        Calculates median depth for each chromosome by random sampling. Thus, measurements differ slightly between runs.
        Acrocentric chromosomes, if defined, are sampled from the last half only to avoid skewed depths in case of too
        many samples from the p-arm.

        :param samples:
        :return:
        """
        print('Analyzing BAM file')
        self.depths = {}
        depth_measurements = {}
        for contig in self.contig_names:
            chr_length = self.contig_lengths[contig]
            start = round(chr_length / 2) if self.chrom_number(contig) in self.acrocentric else 15000
            for i in range(samples):
                rand_pos = max(0, random.randrange(start, chr_length) - 1)
                depth_measurements.setdefault(contig, []).append(self.bam_obj.count(contig, rand_pos, rand_pos + 1))
            self.depths[contig] = median(depth_measurements[contig])

        self.median_depth = median([self.depths[chr] for chr in self.depths if chr not in ['chrX', 'chrY']])
        self.one_x = self.median_depth / 2

        data_frame = pandas.DataFrame.from_dict(depth_measurements, orient='index').T
        figure = plotly_express.box(data_frame, y=list(depth_measurements.keys()), width=1500, height=800)
        figure.update_layout(yaxis_title="Coverage")
        figure.update_layout(yaxis_range=[0, 150])
        figure.update_layout(xaxis_title="Chromosome")
        figure.add_hline(y=self.median_depth, line_width=1, line_dash="dash", line_color="#636EFA")
        self.depth_box_plot_html = figure.to_html(full_html=False, div_id='coverage_plot')


    def read_passes(self, read, filters=None):
        """
        Boolean callback for reads
        :param read:
        :return:
        """
        filters = filters or {
            'mapping_quality': 60,
        }
        for read_filter in filters:
            if not getattr(read, read_filter) == filters[read_filter]:
                return False
        return True


    def regions_have_same_coverage(self, bam_region_1, bam_region_2, sample_size=50):
        """
        Boolean. Returns True if coverage of region 1 and region 2 are not statistically different and False otherwise based
        on a median test of coverage samples made at random points in the provided regions.

        :param bam_region_1:    Tuple/list: (contig, start, stop)
        :param bam_region_2:
        :param sample_size:     How many samples to make in the provided regions
        :return:
        """
        region_1_values = []
        region_2_values = []
        for i in range(sample_size):
            region_1_pos = random.randrange(bam_region_1[1], bam_region_1[2])
            region_2_pos = random.randrange(bam_region_2[1], bam_region_2[2])
            # region_1_values.append(self.bam_obj.count(bam_region_1[0], region_1_pos, region_1_pos + 1, read_callback=read_passes))
            # region_2_values.append(self.bam_obj.count(bam_region_2[0], region_2_pos, region_2_pos + 1, read_callback=read_passes))
            region_1_values.append(self.bam_obj.count(bam_region_1[0], region_1_pos, region_1_pos + 1))
            region_2_values.append(self.bam_obj.count(bam_region_2[0], region_2_pos, region_2_pos + 1))


        res = stats.median_test(region_1_values, region_2_values)
        return res.pvalue > 0.05


    # def bam_file_ok(self):



    def __init__(self, bam_file, nuclear_genome_contigs=24, acrocentric=(13, 14, 15, 21, 22), telomere_length=15000):  # 24 is for the human genome
        """

        :param bam_file:
        :param nuclear_genome_contigs:     Number of contigs in the nuclear genome assembly (autosomes and sex-chromosomes)  == 24 for humans
        :param acrocentric:
        :param telomere_length:            Proximal telomere length at birth (excluded from coverage measurements)
        """
        self.bam_file = bam_file
        self.bam_obj = pysam.AlignmentFile(bam_file, "rb")
        self.nuclear_genome_contigs = nuclear_genome_contigs
        self.acrocentric = acrocentric
        self.telomere_length = telomere_length
        self.contig_lengths = self.get_contig_lengths()
        self.contig_names = get_sorted_alphanum_strings_ascending(list(self.contig_lengths.keys()))
        self.genome_length = sum(self.contig_lengths.values())
        self.depths = None
        self.median_depth = None
        self.depth_box_plot_html = None
        self.one_x = None
        self.get_depths()
        self.pointer = {'contig': 0, 'location': 0}


def get_alignments(bam_obj, contig_areas=None):
    """
    Reads iterator. Iterates through reads in a bam object. If contig areas is specified, fetching is done relative to
    those. Deliberately not in a class in order to be more compatible with multiprocessing.
    :param bam_obj:         Object obtained with pysam.AlignmentFile(<path_to_bam_file>, "rb")
    :param contig_areas:    List of contig area tuples: [('chr1', 0, 1200341),(...)...]
    :return:
    """
    error = None
    if contig_areas:
        for contig_area in contig_areas:
            print(f'Iterating bam file region: {contig_area[0]}:{contig_area[1]}-{contig_area[2]}')
            try:
                for read in bam_obj.fetch(*contig_area):
                    yield read, error
            except ValueError as error:
                yield None, error
    else:
        print('Iterating bam file')
        for read in bam_obj.fetch():
            yield read, error


def get_overlap(a, b):
    """
    Returns the amount of overlap, in bp
    between a and b as well as start/end coordinates.
    If >0, the number of bp of overlap
    If 0,  they are book-ended.
    If <0, the distance in bp between them
    """
    overlap_start = max(a[0], b[0])
    overlap_end = min(a[1], b[1])
    return overlap_end - overlap_start, (overlap_start, overlap_end)


def is_breakpoint_shared(breakpoint_1, breakpoint_2, max_wobble=20):
    """

    :param breakpoint_1:
    :param breakpoint_2:
    :param max_wobble:          Max allowed distance between the breakpoints to infer shared break_value
    :return:
    """
    relative_wobble = max_wobble / 2  # Because max_wobble is allowed for both breakpoints, we half it to ensure that max allowed distance between BPs equals allowed max_wobble
    return get_overlap((breakpoint_1 - relative_wobble, breakpoint_1 + relative_wobble), (breakpoint_2 - relative_wobble, breakpoint_2 + relative_wobble))[0] > 0


def get_breakpoint_clusters(alignments_list, alignment_key, max_wobble):
    """
    Defining function. Identifies clusters based on overlap of reads and shared breakpoints.
    Given a list of alignments (dicts with keys to start/stop tuples) and a corresponding key indicating read name, it
    returns a list of clusters of overlapping reads relative to the given alignment key. The alignment key provides an
    organizational layer to ensure that we only look for overlaps on the same chromosome/contig. For instance, if we
    want to find clusters of reads that overlap in both the "first" and the "second" chromosome (first in this respect
    is typically just the chromosome with the lowest numerical number relative to the second chromosome), we first call
    this function using either alignment key, and then call it again for each resulting cluster with the other key.

    :param (list) alignments_list:  List of alignments, with each item being an alignment dict for a read. For instance:
                                    {
                                        'first': [1844, 7207, 7207, '+', 'l'],
                                        'second': [26056879, 26061816, 26056879, '-', 'l'],
                                        'read': 'cd9436f1-8176-4292-9633-669a3acfaaba'
                                    }, where the first three numbers in first/second are start, stop, and break.
                                    Each key in the dict corresponds to a chromosome/contig. The alignment key provided
                                    in the next param must correspond to one of the keys. That is, finding overlaps is
                                    only relevant between coordinates belonging to the same contig.

    :param (str) alignment_key:     For instance "first" or "second" if those are the keys used in the alignment list.

    :return (list):                 A list of clusters of overlapping reads. Each cluster is a dict with a reads list
                                    and a breakpoint_1.
    """


    # Sort the alignments relative to the breakpoint_1 for the alignment key
    sorted_alignments = sorted(alignments_list, key=lambda x: x[alignment_key][2])

    breakpoint_clusters = []  # Returned
    break_key = f'break_{alignment_key}'

    cluster = {}

    for alignment in sorted_alignments:
        # alignment_added = False

        # Initiate the first cluster with the read/alignment (one-read cluster)
        if not cluster:
            cluster =  {
                    # We use lists in order to accommodate multiple SAs from the same read - duplex reads etc.
                    'alignments': [alignment],
                    'read_names': {alignment['read']},  # Using set of read names to id cluster length as a read may have more alignments with the same breakpoint_1 - for instance duplex reads
                    'breaks': [alignment[alignment_key][2]],  # List of breakpoints for median calculation
                    break_key: alignment[alignment_key][2]  # Breakpoint
            }
        else:
            # Add alignment to cluster if they share breakpoint_1
            if is_breakpoint_shared(alignment[alignment_key][2], cluster[break_key], max_wobble):
                cluster['alignments'].append(alignment)
                cluster['read_names'].add(alignment['read'])
                cluster['breaks'].append(alignment[alignment_key][2])
                cluster[break_key] = round(median(cluster['breaks']))

            else:  # Conclude the cluster and start a new from the read
                if len(cluster['read_names']) > 1:  # only include the cluster if it is supported by more than one read.
                    breakpoint_clusters.append(
                        {
                            'alignments': cluster['alignments'],
                            break_key: cluster[break_key]
                        }
                    )
                cluster = {
                    'alignments': [alignment],
                    'read_names': {alignment['read']},
                    'breaks': [alignment[alignment_key][2]],
                    break_key: alignment[alignment_key][2]
                }

    if len(cluster['read_names']) > 1:  # only include the cluster if it is supported by more than one read.
        breakpoint_clusters.append(
            {
                'alignments': cluster['alignments'],
                break_key: cluster[break_key]
            }
        )
    return breakpoint_clusters


def get_chromosome_to_chromosome_clusters(chromosome_to_different_chromosome_reads, max_wobble):
    """

    :param chromosome_to_different_chromosome_reads:
    :param per_chromosome_coverage:
    :return:
    """

    ##############################################################
    #     Get reads that overlap in the first chromosome        #
    ##############################################################

    first_chrom_clusters = {}

    print('Finding clusters of reads that overlap in the first chromosome')
    for first_chrom in chromosome_to_different_chromosome_reads:
        second_chromosomes = chromosome_to_different_chromosome_reads[first_chrom]

        for second_chrom in second_chromosomes:
            read_entries = second_chromosomes[second_chrom]

            first_chrom_clusters.setdefault(
                first_chrom,
                {}
            )[second_chrom] = get_breakpoint_clusters(read_entries, 'first', max_wobble)

    ##########################################################################
    #     Get reads that overlap in both the first and second chromosome    #
    ##########################################################################

    first_and_second_chrom_clusters = {}
    print('Finding clusters of reads that overlap in both the first and the second chromosome')
    for first_chrom in first_chrom_clusters:
        second_chromosomes = first_chrom_clusters[first_chrom]
        for second_chrom in second_chromosomes:
            alignment_clusters_with_first_chrom_align_overlap = second_chromosomes[second_chrom]

            for alignment_cluster in alignment_clusters_with_first_chrom_align_overlap:
                alignments_with_first_chrom_align_overlap = alignment_cluster['alignments']

                first_and_second_bp_cluster = get_breakpoint_clusters(alignments_with_first_chrom_align_overlap, 'second', max_wobble)

                for cluster in first_and_second_bp_cluster:
                    cluster['break_first'] = alignment_cluster['break_first']

                first_and_second_chrom_clusters.setdefault(
                    first_chrom,
                    {}
                ).setdefault(
                    second_chrom,
                    []
                ).extend(first_and_second_bp_cluster )

    return first_and_second_chrom_clusters


def get_left_clipping(alignment):
    """
    Returns the number of left-clipped bases for an alignment.
    Based on the valid assumption that a CIGAR string always starts with a number.

    :param alignment:
    :return: int
    """
    cigar_string = alignment.cigarstring
    for i, char in enumerate(cigar_string):
        if not char.isnumeric():
            if char in ['S', 'H']:
                return int(cigar_string[:i])
            else:
                return 0


def get_right_clipping(alignment):
    """
    Returns the number of right-clipped bases for an alignment.
    Based on the valid assumption that a CIGAR string always ends with a letter.

    :param alignment:
    :return: int
    """
    cigar_string = alignment.cigarstring
    reverse_cigar = enumerate(reversed(cigar_string))
    if reverse_cigar.__next__()[1] in ['S', 'H']:
        for i, char in reverse_cigar:
            if not char.isnumeric():
                return int(cigar_string[-i:-1])
    else:
        return 0


def is_left_clip_soft(alignment):
    """
    Boolean. Returns True if the left clipping is soft
    :param alignment:
    :return:
    """
    for char in alignment.cigarstring:
        if not char.isnumeric():
            return char == 'S'


def is_right_clip_soft(alignment):
    """
    Boolean. Returns True if the right clipping is soft
    :param alignment:
    :return:
    """
    return alignment.cigarstring.endswith('S')


def is_relative_left_clip_soft(alignment):
    """
    Boolean
    :param alignment:
    :return:
    """
    if alignment.is_reverse:
        return is_left_clip_soft(alignment)
    if alignment.is_forward:
        is_right_clip_soft(alignment)


def is_relative_right_clip_soft(alignment):
    """
    Boolean
    :param alignment:
    :return:
    """
    if alignment.is_forward:
        return is_left_clip_soft(alignment)
    if alignment.is_reverse:
        is_right_clip_soft(alignment)


def get_relative_left_clipping(alignment):
    """
    Returns left-clipping relative to the 5'-end of the query sequence. It returns left-clipping if strand is + and
    right-clipping if strand is -.

    :param alignment:
    :return:
    """
    if alignment.is_forward:
        return get_left_clipping(alignment)
    if alignment.is_reverse:
        return get_right_clipping(alignment)


def get_relative_right_clipping(alignment):
    """
    Returns left-clipping relative to the 5'-end of the query sequence. It returns left-clipping if strand is + and
    right-clipping if strand is -.

    :param alignment:
    :return:
    """
    if alignment.is_reverse:
        return get_left_clipping(alignment)
    if alignment.is_forward:
        return get_right_clipping(alignment)


def is_reverse_direction(alignment):
    """
    Not needed because alignment already has an is_forward attribute.
    Boolean: returns 1 if the reverse complement flag has been set and 0 otherwise.
    Based on decoding the flag value: if the value (int) is written as binary, each number from the right corresponds to
    a specific flag in the SAM bitwise FLAG table going from top to bottom. If a position is 1, the flag is set. If it
    is 0, the flag is not set.

    :param alignment:
    :return: bool
    """
    return int(list(bin(alignment.flag)[:1:-1] + '0000')[4])  # Adding four zeros to ensure index 4 is always in the list range.


def get_alignment_number(alignment):
    """
    Returns the position/index of the current alignment relative to the order of all supplementary alignments the read
    has in a 5'-3' direction. 5'-3'-directionality is necessary because a read may have alignments to both the main
    strand and the complementary strand, which would otherwise switch positions, i.e. the last would become the first.

    Alternative alignments may not always be consecutive because the same stretch of sequence may align to different
    locations. For those cases, alternative alignments of the same sequence or partially the same will receive the same
    alignment number in this function.

    :param alignment:
    :return:
    """
    sa_tag = alignment.get_tag('SA')
    supplementary_alignments = sa_tag.split(';')[:-1]  # deleting the last item because the string always ends with a semicolon
    supplementary_alignments_relative_clippings = []
    for sa in supplementary_alignments:
        sa_data = sa.split(',')
        supplementary_alignments_relative_clippings.append(get_relative_left_clipping(sa_data[3], sa_data[2]))
    sorted_relative_clippings = sorted(supplementary_alignments_relative_clippings)

    relative_left_clipping =  get_relative_left_clipping(alignment.cigarstring, '+') if alignment.is_forward else get_relative_left_clipping(alignment.cigarstring, '-')

    return bisect(sorted_relative_clippings, relative_left_clipping)


def get_number_of_supplementary_alignments(alignment):
    sa_tag = alignment.get_tag('SA')
    return len(sa_tag.split(';')[:-1])


def get_query_start_stop(alignment):
    """
    Returns start and stop coordinates of the query sequence in a 5'-3' direction
    :param alignment:
    :return:
    """
    read_len = alignment.infer_read_length()
    start_clip = get_left_clipping(alignment)
    end_clip = get_right_clipping(alignment)

    if alignment.is_reverse:
        seq_start = end_clip
        seq_end = read_len - start_clip
    else:
        seq_start = start_clip
        seq_end = read_len - end_clip

    return seq_start, seq_end


def alt_get_query_start_stop(alignment):
    """
    Returns start and stop coordinates of the query sequence in a 5'-3' direction
    :param alignment:
    :return:
    """
    read_len = alignment.infer_read_length()
    skip_start = alignment.cigartuples[0][1] if alignment.cigartuples[0][0] in [4, 5] else 0 # start klipping
    skip_end = alignment.cigartuples[-1][1] if alignment.cigartuples[-1][0] in [4, 5] else 0 # end klipping

    if alignment.is_reverse:
        seq_start = skip_end
        seq_end = read_len - skip_start
    else:
        seq_start = skip_start
        seq_end = read_len - skip_end

    return seq_start, seq_end



def get_sa_reads(bam_file, contig_areas=None, min_map_qual=20):
    """
    Identifies reads that contain an SA tag in either a full bam file or specific areas in the bam file.
    :param bam_file:
    :param contig_areas:
    :param min_map_qual:
    :return:

    For each alignment:
    The SA tag is a semicolon-delimited list with each element containing: (contig_name , pos, strand, CIGAR, mapQ, NM;)
    alignment.query_alignment_start     denotes the index (0-based, inclusive) of the first aligned (non-clipped) base.
    alignment.query_alignment_end       denotes the index (0-based, exclusive) of the first base after the last aligned.
    alignment.reference_start
    alignment.reference_end

    """
    reads_dict = {}
    errors = set()

    bam_obj = pysam.AlignmentFile(bam_file, "rb")

    for alignment, error in get_alignments(bam_obj, contig_areas):

        if alignment.has_tag('SA') and alignment.mapping_quality >= min_map_qual:
            query_start, query_stop = get_query_start_stop(alignment)

            reads_dict.setdefault(alignment.query_name, []).append(
                {
                    'chr': alignment.reference_name,
                    'start': alignment.reference_start,
                    'end': alignment.reference_end,
                    'query_start': query_start,
                    'query_stop': query_stop,
                    'orient': '+' if alignment.is_forward else '-',
                    'qual': alignment.mapping_quality,
                }
            )

        if error:
            errors.add(str(error))
    return reads_dict, errors


def pair_alignments(read_name, read_alignments, max_alignments=10, max_consecutive_distance=100, max_relative_overlap=0.5):
    """
    Pairs up alignments from a read and assigns clipping and breakpoints to each alignment.

    Delete after test:
    reads_dict.setdefault(alignment.query_name, []).append(
                {
                    'chr': alignment.reference_name,
                    'start': alignment.reference_start,
                    'end': alignment.reference_end,
                    'query_start': alignment.query_alignment_start,
                    'query_stop': alignment.query_alignment_end,
                    'orient': '+' if alignment.is_forward else '-',
                    'qual': alignment.mapping_quality
                }
            )

            We must end up with:

            first_chrom = get_sorted_alphanum_strings_ascending([align_1['chr'], align_2['chr']])[0]
                    if first_chrom == align_1['chr']:
                        first = align_1
                        second = align_2
                    else:
                        first = align_2
                        second = align_1

                    chromosome_to_different_chromosome_reads.setdefault(
                        first['chr'],
                        {}
                    ).setdefault(
                        second['chr'],
                        []
                    ).append(
                        {
                            'first': (first['start'], first['end'], first['break'], first['orient'], first['clip']),
                            'second': (second['start'], second['end'], second['break'], second['orient'], second['clip']),
                            'read': read_name
                        }
                    )

    :param max_relative_overlap:
    :param read_name:
    :param (str) read_alignments:
    :param (int) max_alignments:              The maximum number of alignments for a read. Setting this low can remove noise and shorten processing time
    :param (int) max_consecutive_distance:    Maximum distance in both directions between stop and start for two alignments to be deemed consecutive. Making this big can potentially increase sensitivity without compromising precision because we still use the more stringent breakpoints for clustering
    :return (dict):
    """
    alignment_pairs = {}
    if len(read_alignments) <= max_alignments:
        sorted_alignments = sorted(read_alignments, key=lambda x:x['query_start'])

        for index in range(len(sorted_alignments)):
            unpaired_alignment = sorted_alignments[index]
            for i in range(index + 1, len(sorted_alignments)):
                next_alignment = sorted_alignments[i]

                quary_alignments_distance = next_alignment['query_start'] - unpaired_alignment['query_stop']  # Distance between alignments on the query - not on the reference
                relative_overlap = -quary_alignments_distance / (unpaired_alignment['query_stop'] - unpaired_alignment['query_start'])
                consecutive_distance = next_alignment['query_start'] - unpaired_alignment['query_stop']

                if relative_overlap <= max_relative_overlap and consecutive_distance <= max_consecutive_distance:
                    # Process clippings and breakpoints
                    unpaired_alignment['break'] = unpaired_alignment['end'] if unpaired_alignment['orient'] == '+' else unpaired_alignment['start']
                    next_alignment['break'] = next_alignment['start'] if next_alignment['orient'] == '+' else next_alignment['end']
                    unpaired_alignment['clip'] = 'r' if unpaired_alignment['orient'] == '+' else 'l'
                    next_alignment['clip'] = 'l' if next_alignment['orient'] == '+' else 'r'

                    if unpaired_alignment['chr'] != next_alignment['chr']:  # No duplex reads in this category
                        sorted_chromosomes = get_sorted_alphanum_strings_ascending([unpaired_alignment['chr'], next_alignment['chr']])
                        alignment_pair = sorted([unpaired_alignment, next_alignment], key=lambda x: sorted_chromosomes.index(x['chr']))
                    else:
                        alignment_pair = sorted([unpaired_alignment, next_alignment], key=lambda x: x['start'])
                    first = alignment_pair[0]
                    second = alignment_pair[1]
                    alignment_pairs.setdefault(
                        first['chr'],
                        {}
                    ).setdefault(
                        second['chr'],
                        []
                    ).append(
                        {
                            'first': (first['start'], first['end'], first['break'], first['orient'], first['clip']),
                            'second': (second['start'], second['end'], second['break'], second['orient'], second['clip']),
                            'read': read_name
                        }
                    )
    return alignment_pairs


def get_chromosome_to_chromosome_reads(bam_data, min_map_qual=20, cores=1):
    """
    Generates and returns two dictionaries:
    - A chromosome_to_chromosome_reads "matrix" layout for reads with supplementary alignments, where one
    alignment is in one chromosome and the other alignment is in another chromosome.
    - A "same chromosome" layout."

    The chromosome_to_chromosome_reads layout has the form:

    {
        'chr1': {
            'chr2': [
                {
                    'first': (first['start'], first['end'], first['break'], first['orient'], first['clip']),
                    'second': (second['start'], second['end'], second['break'], second['orient'], second['clip']),
                    'read': read_name
                },
                {...}
            ],
            'chr3': [...]
        },
        'chr2': {
            'chr3': [...],
            ...
        }
        ...
    }

    For entries where first and second chromosome are equal, "first" and "second" are organized relative to start so
    first['start'] < second['start']. start, end, and break values are integers; orient is '+' or '-'; clip is 'l' or
    'r'.

    :param bam_data:                BamData object
    :param int min_map_qual:        Default is currently 20, which in minimap based mappers will mostly provide reads with a quality = 60
    :param bool development_mode:
    :return:
    """
    
    ####################################################################################################
    #                                    BAM file search for SA reads                                  #
    ####################################################################################################

    t0 = time.time()

    reads_dict = {}  # Contains alignments associated to reads having supplementary alignments

    print('Sorting out reads that have supplementary alignments')
    if cores==1:
        reads_dict = get_sa_reads(bam_data.bam_file)
        t1 = time.time()
    else:
        print('Running multiprocessing')
        reference_subdivisions = bam_data.get_reference_subdivisions(cores)

        # Get SA reads in as many processes as there are cores requested
        with multiprocessing.Pool(cores) as pool:
            processes = [
                pool.apply_async(get_sa_reads, args=(bam_data.bam_file, reference_subdivisions[i], min_map_qual))
                for i in range(cores)
            ]
            process_tuples = [process.get() for process in processes]
            errors = [error for process_tuple in process_tuples for error in process_tuple[1]]
            for process_tuple in process_tuples:
                for read in process_tuple[0]:
                    reads_dict[read] = reads_dict.get(read, []) + process_tuple[0][read]

            if errors:
                # Errors may happen if for instance cancerous cell lines are analyzed. They can sometimes lose entire chromosomes causing ValueErrors if random access queries are used rather than end-to-end iteration.
                print('Some errors were thrown during bam file processing. they have been written to: pysam_fetch_errors.txt')
                errors = sorted(list(set(errors)))
                with open('pysam_fetch_errors.txt', 'r') as error_file:
                    error_file.write('\n'.join(errors))

            t1 = time.time()

    print(f'BAM file processing finished in {t1 - t0} seconds using {cores} cores.')

    ####################################################################################################
    #                              Read organization relative to chromosomes                           #
    ####################################################################################################

    print('Sanitizing and organizing reads into chromosome matrices')

    """
    Here we make pairwise combinations of those reads having more than two alignments and define breakpoints
    """

    chromosome_to_chromosome_reads = {}

    for read_name in reads_dict:
        read_alignments = reads_dict[read_name]
        if len(read_alignments) == 1:
            continue  #  There may be only one alignment if the quality of the other alignments was below the minimum alignment score
        else:
            paired_alignments = pair_alignments(read_name, reads_dict[read_name])

            for first_chrom in paired_alignments:
                for second_chrom in paired_alignments[first_chrom]:
                    chromosome_to_chromosome_reads.setdefault(
                        first_chrom,
                        {}
                    ).setdefault(
                        second_chrom,
                        []
                    ).extend(paired_alignments[first_chrom][second_chrom])

    return chromosome_to_chromosome_reads


def generate_max_cluster_count_heatmap_and_json(cluster_dict, chromosomes, output_folder, sample_suffix, format='jpg'):
    """
    Given a dictionary of chromosome to chromosome overlapping clusters, this function generates a heatmap of the
    number of reads in the cluster that has the most reads. A cluster is a collection of reads that overlap in both
    chromosomes. it also generates a json file of clusters sorted by number of reads per cluster. That way, for each
    chromosome-to-chromosome entry, the cluster with the most reads will be placed first.

    Args:
        self.cluster_dict (dict): Dictionary (chrom to chrom matrix) of overlapping clusters
        file_name (str):     File name for the figure file
        format:              For example pdf / png / svg / jpg
    """

    import plotly.graph_objects as plotly_graph_objects
    import plotly.io as plotly_io
    plotly_io.templates.default = "none"

    print('Generating heatmap...')

    chrom_None_row = {chrom: None for chrom in chromosomes}

    # contains max counts of reads that overlap (with clipping included) in both chromosomes
    chrom_to_chrom_max_overlapping_reads = {}

    for first_chrom in chromosomes:
        if first_chrom in cluster_dict:
            chrom_row = {}
            for second_chrom in chromosomes:
                if second_chrom in cluster_dict[first_chrom] and cluster_dict[first_chrom][second_chrom]:

                    # Sort clusters relative to number of reads in the clusters in descending order
                    sorted_clusters = sorted(cluster_dict[first_chrom][second_chrom],
                                                 key=lambda x: len(x['alignments']), reverse=True)
                    cluster_dict[first_chrom][second_chrom] = sorted_clusters

                    # Since we have sorted, the cluster in the first entry will contain the most reads
                    chrom_row[second_chrom] = len(cluster_dict[first_chrom][second_chrom][0]['alignments'])
                else:
                    chrom_row[second_chrom] = None
            chrom_to_chrom_max_overlapping_reads[first_chrom] = chrom_row
        else:
            chrom_to_chrom_max_overlapping_reads[first_chrom] = chrom_None_row

    # Save cluster dict in a json file to make it easier to inspect positions and reads for those with high counts
    with open(os.path.join(output_folder, f'{sample_suffix}sorted_pairwise_clusters.json'), 'w') as json_file:
        json_file.write(json.dumps(cluster_dict))

    plot_file_name = os.path.join(output_folder, f'{sample_suffix}max_count_clusters_heatmap')

    chrom_to_chrom_alignments_data_frame = pandas.DataFrame.from_dict(chrom_to_chrom_max_overlapping_reads)

    heat = plotly_graph_objects.Heatmap(
        z=chrom_to_chrom_alignments_data_frame,
        x=chromosomes,
        y=chromosomes,
        xgap=1, ygap=1,
        colorscale='greens',
        colorbar_thickness=15,
        colorbar_ticklen=2,
    )


    layout = plotly_graph_objects.Layout(
        title_text='Maximum cluster read counts', title_x=0.5,
        width=600, height=600,
        xaxis_showgrid=False,
        yaxis_showgrid=False,
        yaxis_autorange='reversed',
        font=dict(
            # family="Courier New, monospace",
            # size=18,
            color="#2e2d2d"
        )
    )

    figure = plotly_graph_objects.Figure(data=[heat], layout=layout)
    figure.write_image(f"{plot_file_name}.{format}")


def generate_cluster_read_counts_jsons(cluster_dict, output_folder):
    """
    Creates a range of json files - one for each read counts containing clusters supported by reads of that count and
    returns a html_table with columns: 'Reads per cluster' and 'Cluster counts'
    :param self.cluster_dict:
    :return: (str) html_table
    """

    os.makedirs(output_folder, exist_ok=True)
    per_read_count_clusters = {}

    counter = {}
    for first_chromosome in cluster_dict:
        for second_chromosome in cluster_dict[first_chromosome]:
            entry = cluster_dict[first_chromosome][second_chromosome]
            for cluster in entry:
                number_of_supporting_reads = len(cluster['alignments'])
                counter[number_of_supporting_reads] = counter.get(number_of_supporting_reads, 0) + 1

                per_read_count_clusters.setdefault(number_of_supporting_reads, {}).setdefault(first_chromosome, {}).setdefault(second_chromosome, []).append(cluster)

    for number_of_supporting_reads in per_read_count_clusters:
        with open(os.path.join(output_folder, f'{number_of_supporting_reads}.json'), 'w') as json_file:
            json_file.write(json.dumps(per_read_count_clusters[number_of_supporting_reads]))

    sorted_counter = {}
    for number in sorted(counter):
        sorted_counter[number] = counter[number]

    hist_data = {
        'Reads per cluster': [],
        'Cluster counts': []
    }

    for cluster_read_counts in sorted_counter:
        hist_data['Reads per cluster'].append(cluster_read_counts)
        hist_data['Cluster counts'].append(sorted_counter[cluster_read_counts])

    data_frame = pandas.DataFrame.from_dict(hist_data)
    df_html = data_frame.to_html(col_space=200, index=False)
    # fig = plotly_express.histogram(data_frame, x="Reads per cluster", y='Cluster counts', nbins=max(list(counter.keys())))
    # fig.update_yaxes(range=[0, 100])
    # fig.update_layout(yaxis_title="cluster counts")
    # fig.update_layout(xaxis_title="Reads per cluster")
    # fig.write_image('cluster_read_counts.jpg')
    return df_html


class StuctVarFinder:
    """
    Identifies structural variants based on unique cluster signatures.

    Todo: If an annotation file is given, it will be able to flag SVs that overlap genes of various subtypes and provide
    insertion statistics.

    The operations in this class are quite fast and therefore each SV subtype has its own method even though each method
    must iterate the cluster dictionary independently. This creates more manageable code.

    The strategy for each cluster list that relates to a chromosome pair is to iterate through it once and for each SV
    type, we look for the corresponding event in both an upstream and a downstream configuration when it makes sense.
    """

    def get_cluster_read_clipping_counts(self, cluster, alignment_key):
        """
        Given a cluster it returns a 2-tuple with number of left-clipped reads at index 0 and number of right clipped reads
        at index 1
        :param cluster:
        :param alignment_key:
        :return: (int: left, int: right)
        """
        left = 0
        right = 0
        for read in cluster['alignments']:
            if read[alignment_key][4] == 'l':
                left += 1
            elif read[alignment_key][4] == 'r':
                right += 1
        return left, right


    def only_left_clipped_alignments(self, cluster, alignment_key):
        """
        Boolean
        :param cluster:
        :param alignment_key:
        :return: bool
        """
        left_clipped, right_clipped = self.get_cluster_read_clipping_counts(cluster, alignment_key)
        return left_clipped > right_clipped == 0


    def only_right_clipped_alignments(self, cluster, alignment_key):
        """
        Boolean
        :param cluster:
        :param alignment_key:
        :return: bool
        """
        left_clipped, right_clipped = self.get_cluster_read_clipping_counts(cluster, alignment_key)
        return right_clipped > left_clipped == 0


    def left_and_right_clipped_alignments(self, cluster, alignment_key):
        """
        Boolean
        :param cluster:
        :param alignment_key:
        :return: bool
        """
        left_clipped, right_clipped = self.get_cluster_read_clipping_counts(cluster, alignment_key)
        return left_clipped > 0 and right_clipped > 0


    def only_left_clipped_cluster(self, cluster):
        return self.only_left_clipped_alignments(cluster, 'first') and self.only_left_clipped_alignments(cluster, 'second')


    def only_right_clipped_cluster(self, cluster):
        return self.only_right_clipped_alignments(cluster, 'first') and self.only_right_clipped_alignments(cluster, 'second')


    def polarized_clipping_cluster(self, cluster):
        """
        Boolean. Returns True if both alignments (first and second) are clipped in the same direction
        :param cluster:
        :return:
        """
        return self.only_left_clipped_cluster(cluster) or self.only_right_clipped_cluster(cluster)


    def outward_clipping_cluster(self, cluster):
        """
        Returns True if cluster is left-clipped at first alignment key and right-clipped at second
        :param cluster:
        :return:
        """
        return self.only_left_clipped_alignments(cluster, 'first') and self.only_right_clipped_alignments(cluster, 'second')


    def inward_clipping_cluster(self, cluster):
        """
        Returns True if cluster is right-clipped at first alignment key and left-clipped at second
        :param cluster:
        :return:
        """
        return self.only_right_clipped_alignments(cluster, 'first') and self.only_left_clipped_alignments(cluster, 'second')


    def bidirectional_clipping_cluster(self, cluster):
        """
        Boolean. Returns True if clipping is both to the left and the right at both alignment keys (both qq and pq)
        :param cluster:
        :return:
        """
        return self.left_and_right_clipped_alignments(cluster, 'first') and self.left_and_right_clipped_alignments(cluster, 'second')


    def get_unbalanced_translocation_data(self, cluster, first_chromosome, second_chromosome):
        if self.bam_data.has_unbalanced_translocation_coverage_profile(
                (first_chromosome, cluster['break_first']),
                (second_chromosome, cluster['break_second'])
        ):
            return {
                'chromosome_1': first_chromosome,
                'chromosome_2': second_chromosome,
                'breakpoint_1': cluster['break_first'],
                'breakpoint_2': cluster['break_second'],
                'clusters': cluster,
                'length': 'Nan',
                'igv_1': f"{first_chromosome}:{cluster['break_first']}",
                'igv_2': f"{second_chromosome}:{cluster['break_second']}",
            }
        return None


    def get_translocation_with_malsegregation_data(self, cluster, first_chromosome, second_chromosome):
        if self.bam_data.has_translocation_with_malsegregation_coverage_profile(
                (first_chromosome, cluster['break_first']),
                (second_chromosome, cluster['break_second'])
        ):
            return {
                'chromosome_1': first_chromosome,
                'chromosome_2': second_chromosome,
                'breakpoint_1': cluster['break_first'],
                'breakpoint_2': cluster['break_second'],
                'clusters': cluster,
                'length': 'Nan',
                'igv_1': f"{first_chromosome}:{cluster['break_first']}",
                'igv_2': f"{second_chromosome}:{cluster['break_second']}",
            }
        return None


    def get_balanced_tranlocation_with_deletion_data(self, break_first_clusters, break_second_clusters, first_chromosome, second_chromosome):
        if len(break_first_clusters) == 2 and len(break_second_clusters) == 1:
            sorted_clusters = sorted(break_first_clusters, key=lambda x:x['break_second'])
            if self.bam_data.has_balanced_translocation_with_deletion_coverage_profile(
                    (first_chromosome, sorted_clusters[0]['break_first']),
                    (second_chromosome, sorted_clusters[0]['break_second']),
                    (second_chromosome, sorted_clusters[1]['break_second'])
            ):
                return {
                    'chromosome_1': first_chromosome,
                    'chromosome_2': second_chromosome,
                    'breakpoint': sorted_clusters[0]['break_first'],
                    'del_start': sorted_clusters[0]['break_second'],
                    'del_stop': sorted_clusters[1]['break_second'],
                    'length': sorted_clusters[1]['break_second'] - sorted_clusters[0]['break_second'],
                    'clusters': sorted_clusters,
                    'igv_1': f"{first_chromosome}:{sorted_clusters[0]['break_first']}",
                    'igv_2': f"{second_chromosome}:{sorted_clusters[0]['break_second']}-{sorted_clusters[1]['break_second']}",
                }
        elif len(break_first_clusters) == 1 and len(break_second_clusters) == 2:
            sorted_clusters = sorted(break_second_clusters, key=lambda x: x['break_first'])
            if self.bam_data.has_balanced_translocation_with_deletion_coverage_profile(
                    (second_chromosome, sorted_clusters[0]['break_second']),
                    (first_chromosome, sorted_clusters[0]['break_first']),
                    (first_chromosome, sorted_clusters[1]['break_first'])
            ):
                return {
                    'chromosome_1': first_chromosome,
                    'chromosome_2': second_chromosome,
                    'breakpoint': sorted_clusters[0]['break_second'],
                    'del_start': sorted_clusters[0]['break_first'],
                    'del_stop': sorted_clusters[1]['break_first'],
                    'length': sorted_clusters[1]['break_first'] - sorted_clusters[0]['break_first'],
                    'clusters': sorted_clusters,
                    'igv_1': f"{second_chromosome}:{sorted_clusters[0]['break_second']}",
                    'igv_2': f"{first_chromosome}:{sorted_clusters[0]['break_first']}-{sorted_clusters[1]['break_first']}",
                }
        return None


    def get_translocation_with_malsegregation_data_and_deletion(self, break_first_clusters, break_second_clusters, first_chromosome, second_chromosome):
        if len(break_first_clusters) == 2 and len(break_second_clusters) == 1:
            sorted_clusters = sorted(break_first_clusters, key=lambda x:x['break_second'])
            if self.bam_data.has_translocation_with_malsegregation_and_deletion_coverage_profile(
                    (first_chromosome, sorted_clusters[0]['break_first']),
                    (second_chromosome, sorted_clusters[0]['break_second']),
                    (second_chromosome, sorted_clusters[1]['break_second'])
            ):
                return {
                    'chromosome_1': first_chromosome,
                    'chromosome_2': second_chromosome,
                    'breakpoint': sorted_clusters[0]['break_first'],
                    'del_start': sorted_clusters[0]['break_second'],
                    'del_stop': sorted_clusters[1]['break_second'],
                    'length': sorted_clusters[1]['break_second'] - sorted_clusters[0]['break_second'],
                    'clusters': sorted_clusters,
                    'igv_1': f"{first_chromosome}:{sorted_clusters[0]['break_first']}",
                    'igv_2': f"{second_chromosome}:{sorted_clusters[0]['break_second']}-{sorted_clusters[1]['break_second']}",
                }
        elif len(break_first_clusters) == 1 and len(break_second_clusters) == 2:
            sorted_clusters = sorted(break_second_clusters, key=lambda x: x['break_first'])
            if self.bam_data.has_translocation_with_malsegregation_and_deletion_coverage_profile(
                    (second_chromosome, sorted_clusters[0]['break_second']),
                    (first_chromosome, sorted_clusters[0]['break_first']),
                    (first_chromosome, sorted_clusters[1]['break_first'])
            ):
                return {
                    'chromosome_1': first_chromosome,
                    'chromosome_2': second_chromosome,
                    'breakpoint': sorted_clusters[0]['break_second'],
                    'del_start': sorted_clusters[0]['break_first'],
                    'del_stop': sorted_clusters[1]['break_first'],
                    'length': sorted_clusters[1]['break_first'] - sorted_clusters[0]['break_first'],
                    'clusters': sorted_clusters,
                    'igv_1': f"{second_chromosome}:{sorted_clusters[0]['break_second']}",
                    'igv_2': f"{first_chromosome}:{sorted_clusters[0]['break_first']}-{sorted_clusters[1]['break_first']}",
                }
        return None


    def get_balanced_translocation_with_dual_deletion_data(self, cluster, clusters, first_chromosome, second_chromosome):
        """
        Returns the data if clipping pattern matches the pattern only seen in balanced qq translocations with
        two deletions. To accommodate the deletions, a scanning must be undertaken because the signature will be
        composed of two clusters with no clear link. For the given cluster, the required scanning direction can be
        inferred by looking at the clipping direction since the deletions, if present, will have inwards clipping. To
        follow along with the code, a sketch is recommended.
        If pattern is not found, None is returned.

        Assumptions that must be established by caller:
        - first_chrom != second_chrom
        - cluster is pulled from a position-sorted list of clusters

        :param {} cluster:     Cluster that defines the starting point of the pattern search
        :param [] clusters:    A list of clusters for first and second chromosome (known by caller)
        :return {}:
        """

        # We do not condition on alignment directions in order to accommodate both qq and pq translocations
        if self.only_right_clipped_alignments(cluster, 'first') and self.only_left_clipped_alignments(cluster, 'second'):

            first_chrom_sorted_clusters = sorted(clusters, key=lambda x: x['break_first'])

            # Get index of cluster on the sorted lists
            for index, sorted_cluster in enumerate(first_chrom_sorted_clusters):
                if sorted_cluster['i'] == cluster['i']:
                    first_chrom_cluster_index = index
                    break

            # We should always encounter the start of deletions first in a list of sorted clusters - unless we do not have a deletion
            first_chrom_next_index = first_chrom_cluster_index
            for i in range(5):  # Allow up to five other clusters before scanning is stopped to allow for other SVs within the deletion or noisy environment. Todo: See if this hurts precision
                first_chrom_next_index += 1
                if first_chrom_next_index == len(first_chrom_sorted_clusters):
                    return None

                next_first_chrom_cluster = first_chrom_sorted_clusters[first_chrom_next_index]
                if (self.only_left_clipped_alignments(next_first_chrom_cluster, 'first') and
                        self.only_right_clipped_alignments(next_first_chrom_cluster, 'second') and
                        self.cluster_alignments_are_unidirectional( next_first_chrom_cluster)):
                    first_chrom_start = cluster['break_first']
                    first_chrom_stop = next_first_chrom_cluster['break_first']
                    first_chrom_del_length = first_chrom_stop - first_chrom_start

                    second_chrom_start = next_first_chrom_cluster['break_second']
                    second_chrom_stop = cluster['break_second']
                    second_chrom_del_length = second_chrom_stop - second_chrom_start

                    min_del_length = self.max_breakpoint_wobble * 2  # Less than this is considered same breakpoint_1 and should be picked op elsewhere

                    if first_chrom_del_length > min_del_length and second_chrom_del_length > min_del_length:
                        return {
                            'chromosome_1': first_chromosome,
                            'chromosome_2': second_chromosome,
                            'first_del_start': first_chrom_start,
                            'first_del_stop': first_chrom_stop,
                            'first_del_length': first_chrom_del_length,
                            'second_del_start': second_chrom_start,
                            'second_del_stop': second_chrom_stop,
                            'second_del_length': second_chrom_del_length,
                            'length': max(first_chrom_del_length, second_chrom_del_length),
                            'igv_1': f"{first_chromosome}:{first_chrom_start}-{first_chrom_stop}",
                            'igv_2': f"{second_chromosome}:{second_chrom_start}-{second_chrom_stop}",
                            'clusters': [cluster, next_first_chrom_cluster],
                        }
        return None


    def get_malsegregated_translocation_with_dual_deletion_data(self, cluster, clusters, first_chromosome, second_chromosome):
        """
        Similar to get_balanced_translocation_with_dual_deletion_data, but it also checks coverage pattern.
        Returns the data if clipping pattern matches the pattern only seen in balanced translocations with
        two deletions if coverage pattern matches malsegregation. To accommodate the deletions, a scanning must be
        undertaken because the signature will be composed of two clusters with no clear link. For the given cluster, the
        required scanning direction can be inferred by looking at the clipping direction since the deletions, if
        present, will have inwards clipping. To follow along with the code, a sketch is recommended. If pattern is not
        found, None is returned.

        Assumptions that must be established by caller:
        - first_chrom != second_chrom
        - cluster is pulled from a position-sorted list of clusters

        :param {} cluster:     Cluster that defines the starting point of the pattern search
        :param [] clusters:    A list of clusters for first and second chromosome (known by caller)
        :return {}:
        """

        # We do not condition on alignment directions in order to accommodate both qq and pq translocations
        if self.only_right_clipped_alignments(cluster, 'first') and self.only_left_clipped_alignments(cluster, 'second'):

            first_chrom_sorted_clusters = sorted(clusters, key=lambda x: x['break_first'])

            # Get index of cluster on the sorted lists
            for index, sorted_cluster in enumerate(first_chrom_sorted_clusters):
                if sorted_cluster['i'] == cluster['i']:
                    first_chrom_cluster_index = index
                    break

            # We should always encounter the start of deletions first in a list of sorted clusters - unless we do not have a deletion
            first_chrom_next_index = first_chrom_cluster_index
            for i in range(5):  # Allow up to five other clusters before scanning is stopped to allow for other SVs within the deletion or noisy environment. Todo: See if this hurts precision
                first_chrom_next_index += 1
                if first_chrom_next_index == len(first_chrom_sorted_clusters):
                    return None

                next_first_chrom_cluster = first_chrom_sorted_clusters[first_chrom_next_index]
                if (self.only_left_clipped_alignments(next_first_chrom_cluster,'first') and
                        self.only_right_clipped_alignments(next_first_chrom_cluster, 'second') and
                        self.cluster_alignments_are_unidirectional(next_first_chrom_cluster)):
                    first_chrom_start = cluster['break_first']
                    first_chrom_stop = next_first_chrom_cluster['break_first']
                    first_chrom_del_length = first_chrom_stop - first_chrom_start

                    second_chrom_start = next_first_chrom_cluster['break_second']
                    second_chrom_stop = cluster['break_second']
                    second_chrom_del_length = second_chrom_stop - second_chrom_start

                    min_del_length = self.max_breakpoint_wobble * 2  # Less than this is considered same breakpoint_1 and should be picked op elsewhere

                    if first_chrom_del_length > min_del_length and second_chrom_del_length > min_del_length:
                        if self.bam_data.has_translocation_with_malsegregation_and_double_deletion_coverage_profile(
                                (first_chromosome, first_chrom_start),
                                (first_chromosome, first_chrom_stop),
                                (second_chromosome, second_chrom_start),
                                (second_chromosome, second_chrom_stop)
                        ):
                            return {
                                'chromosome_1': first_chromosome,
                                'chromosome_2': second_chromosome,
                                'first_del_start': first_chrom_start,
                                'first_del_stop': first_chrom_stop,
                                'first_del_length': first_chrom_del_length,
                                'second_del_start': second_chrom_start,
                                'second_del_stop': second_chrom_stop,
                                'second_del_length': second_chrom_del_length,
                                'length': max(first_chrom_del_length, second_chrom_del_length),
                                'igv_1': f"{first_chromosome}:{first_chrom_start}-{first_chrom_stop}",
                                'igv_2': f"{second_chromosome}:{second_chrom_start}-{second_chrom_stop}",
                                'clusters': [cluster, next_first_chrom_cluster],
                            }
        return None


    def cluster_reads_have_inversion_signature(self, cluster):
        """
        Boolean. Tests if pairwise alignments to the first and second chromosomes are clipped in the same direction
        while having bidirectional alignments. This signature is typically seen in inversions or pq translocations.
        :param cluster:
        :return:
        """

        # Alignment direction: index 3, clipping: index 4
        for read in cluster['alignments']:
            if read['first'][3] == read['second'][3] or read['first'][4] != read['second'][4]:
                return False
        return True


    def cluster_alignments_are_unidirectional(self, cluster):
        """
        Boolean. Tests if pairwise alignments to the first and second chromosome have the same direction. If reads are
        not unidirectional, it may indicate an inversion.
        :param cluster:
        :return:
        """
        for read in cluster['alignments']:
            if read['first'][3] != read['second'][3]:
                return False
        return True


    def cluster_alignments_are_bidirectional(self, cluster):
        for read in cluster['alignments']:
            if read['first'][3] == read['second'][3]:
                return False
        return True


    def get_shared_breakpoint_clusters(self, clusters, cluster, break_key):
        """
        Returns a list with clusters that shares the break_value with either break_first or break_second
        :param clusters:
        :param cluster:
        :param break_key:
        :return:
        """
        break_value = cluster[break_key]
        bp_clusters = []
        for i, cluster in enumerate(clusters):
            if len(cluster['alignments']) >= self.min_cluster_reads:
                if is_breakpoint_shared(break_value, cluster[break_key], self.max_breakpoint_wobble):
                    bp_clusters.append(cluster)
        return bp_clusters


    def is_duplex_cluster(self, cluster):
        """
        Boolean. Returns True if the cluster is made up of duplex reads, i.e. reads that contain alignments in
        opposite directions that have the same breakpoint_1.
        :return:    Bool
        """
        for read in cluster['alignments']:
            if not (read['first'][3] != read['second'][3] and is_breakpoint_shared(read['first'][2], read['second'][2], self.max_breakpoint_wobble)):
                return False
        return True


    def cluster_contains_a_duplex_read(self, cluster):
        """
        Boolean. Returns True if cluster contains at least one duplex-like read, i.e. a read containing alignments in
        opposite directions that have the same breakpoint_1. Not as strict as self.is_duplex_cluster and may therefore be
        useful in "noisy" conditions. However, it may also lead to wrongly identifying a cluster as a duplex cluster.
        :return:    Bool
        """
        for read in cluster['alignments']:
            if read['first'][3] != read['second'][3] and is_breakpoint_shared(read['first'][2], read['second'][2], self.max_breakpoint_wobble):
                return True
        return False


    def cluster_contains_read(self, cluster, read_name):
        """
        For debug purposes. Returns True if the cluster contains a read with the name, read_name
        :param cluster:
        :param read_name:
        :return:
        """
        for alignment in cluster['alignments']:
            if alignment['read'] == read_name:
                return True
        return False



    def call_non_inversion_svs(self):
        """
        :return: None   Updates instance
        """

        for first_chromosome in self.cluster_dict:

            for second_chromosome in self.cluster_dict[first_chromosome]:
                clusters = self.cluster_dict[first_chromosome][second_chromosome]

                first_chrom_sorted_clusters = sorted(clusters, key=lambda x: x['break_first'])

                for cluster in first_chrom_sorted_clusters:

                    if len(cluster['alignments']) >= self.min_cluster_reads and self.cluster_alignments_are_unidirectional(cluster):
                        break_first_clusters = self.get_shared_breakpoint_clusters(clusters, cluster, 'break_first')
                        break_second_clusters = self.get_shared_breakpoint_clusters(clusters, cluster, 'break_second')
                        if self.outward_clipping_cluster(cluster):



                            """
                            We have a cluster consisting of only left_clipped reads at the first break_value and only right-clipped reads at the second break_value (outward clipping).
                            Many SVs display this signature, so we nust look for additional, associated signatures to characterize which of several possible SVs we are dealing
                            with. If no additional signatures can reliably be identified, we likely have a tandem duplication as this is the only signature it displays.
                            """

                            if len(break_first_clusters) == len(break_second_clusters) == 1 and first_chromosome == second_chromosome:  # Tandem duplication signature (not including coverage at this point)
                                duplication = {
                                    'chromosome': first_chromosome,
                                    'begin': cluster['break_first'],
                                    'end': cluster['break_second'],
                                    'clusters': cluster,
                                    'length': cluster['break_second'] - cluster['break_first'],
                                    'igv': f'{first_chromosome}:{cluster["break_first"]}-{cluster["break_second"]}'
                                }
                                self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'TANDEM_DUP_{first_chromosome}:{duplication["begin"]}-{duplication["end"]}')
                                self._sv_tandem_duplications.append(duplication)

                            elif len(break_first_clusters) == 1 and len(break_second_clusters) == 2:  # Potential duplication to DOWNSTREAM position on same chromosome or different chromosome

                                break_second_clusters = sorted(break_second_clusters, key=lambda x: x['break_first'])
                                if self.outward_clipping_cluster(break_second_clusters[0]) and self.inward_clipping_cluster(break_second_clusters[1]):
                                    duplication = {
                                        'donor': first_chromosome,
                                        'begin': break_second_clusters[0]['break_first'],
                                        'end': break_second_clusters[1]['break_first'],
                                        'acceptor': second_chromosome,
                                        'insertion_point': cluster['break_second'],
                                        'clusters': break_second_clusters,
                                        'length': break_second_clusters[1]['break_first'] - break_second_clusters[0]['break_first'],
                                        'igv_d': f'{first_chromosome}:{break_second_clusters[0]["break_first"]}-{break_second_clusters[1]["break_first"]}',
                                        'igv_a': f'{second_chromosome}:{cluster["break_second"]}',
                                        'info': 'dup-down first-left second-right clipped',
                                    }
                                    self.cluster_dict[first_chromosome][second_chromosome][break_second_clusters[0]['i']].setdefault('sv', []).append(f'DUP_{first_chromosome}_to_{second_chromosome}')
                                    self.cluster_dict[first_chromosome][second_chromosome][break_second_clusters[1]['i']].setdefault('sv', []).append(f'DUP_{first_chromosome}_to_{second_chromosome}')
                                    self._sv_distant_duplications.append(duplication)

                            elif len(break_first_clusters) == 2 and len(break_second_clusters) == 1:  # Potential duplication to UPSTREAM position on same chromosome or different chromosome
                                break_first_clusters = sorted(break_first_clusters, key=lambda x: x['break_second'])
                                if self.inward_clipping_cluster(break_first_clusters[0]) and self.outward_clipping_cluster(break_first_clusters[1]):

                                    duplication = {
                                        'donor': second_chromosome,
                                        'begin': break_first_clusters[0]['break_second'],
                                        'end': break_first_clusters[1]['break_second'],
                                        'acceptor': first_chromosome,
                                        'insertion_point': cluster['break_first'],
                                        'clusters': break_first_clusters,
                                        'length': break_first_clusters[1]['break_second'] - break_first_clusters[0]['break_second'],
                                        'igv_d': f'{second_chromosome}:{break_first_clusters[0]["break_second"]}-{break_first_clusters[1]["break_second"]}',
                                        'igv_a': f'{first_chromosome}:{cluster["break_first"]}',
                                        'info': 'dup-up first-left second-right clipped',
                                    }
                                    self.cluster_dict[first_chromosome][second_chromosome][break_first_clusters[0]['i']].setdefault('sv', []).append(f'DUP_{second_chromosome}_to_{first_chromosome}')
                                    self.cluster_dict[first_chromosome][second_chromosome][break_first_clusters[1]['i']].setdefault('sv', []).append(f'DUP_{second_chromosome}_to_{first_chromosome}')
                                    self._sv_distant_duplications.append(duplication)

                        elif self.inward_clipping_cluster(cluster):  #
                            if len(break_first_clusters) == len(break_second_clusters) == 1 and first_chromosome == second_chromosome:
                                deletion = {
                                    'chr': first_chromosome,
                                    'begin': cluster['break_first'],
                                    'end': cluster['break_second'],
                                    'length': cluster['break_second'] - cluster['break_first'],
                                    'clusters': cluster,
                                    'igv': f'{first_chromosome}:{cluster["break_first"]}-{cluster["break_second"]}'
                                }
                                self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'DEL_{first_chromosome}:{deletion["begin"]}-{deletion["end"]}')
                                self._sv_deletions.append(deletion)

                            elif len(break_first_clusters) == 2 and len(break_second_clusters) == 1:  # Potential duplication to upstream position on same chromosome or different chromosome
                                break_first_clusters = sorted(break_first_clusters, key=lambda x: x['break_second'])
                                if self.inward_clipping_cluster(break_first_clusters[0]) and self.outward_clipping_cluster(break_first_clusters[1]):
                                    duplication = {
                                        'donor': second_chromosome,
                                        'begin': break_first_clusters[0]['break_second'],
                                        'end': break_first_clusters[1]['break_second'],
                                        'acceptor': first_chromosome,
                                        'insertion_point': cluster['break_first'],
                                        'clusters': break_first_clusters,
                                        'length': break_first_clusters[1]['break_second'] - break_first_clusters[0]['break_second'],
                                        'igv_d': f'{second_chromosome}:{break_first_clusters[0]["break_second"]}-{break_first_clusters[1]["break_second"]}',
                                        'igv_a': f'{first_chromosome}:{cluster["break_first"]}',
                                        'info': 'dup-up first-right second-left clipped',
                                    }
                                    self.cluster_dict[first_chromosome][second_chromosome][break_first_clusters[0]['i']].setdefault('sv', []).append(f'DUP_{second_chromosome}_to_{first_chromosome}')
                                    self.cluster_dict[first_chromosome][second_chromosome][break_first_clusters[1]['i']].setdefault('sv', []).append(f'DUP_{second_chromosome}_to_{first_chromosome}')
                                    self._sv_distant_duplications.append(duplication)

                            elif len(break_first_clusters) == 1 and len(break_second_clusters) == 2:  # Potential duplication to downstream position on same chromosome or different chromosome
                                break_second_clusters = sorted(break_second_clusters, key=lambda x: x['break_first'])
                                if self.outward_clipping_cluster(break_second_clusters[0]) and self.inward_clipping_cluster(break_second_clusters[1]):
                                    duplication = {
                                        'donor': first_chromosome,
                                        'begin': break_second_clusters[0]['break_first'],
                                        'end': break_second_clusters[1]['break_first'],
                                        'acceptor': second_chromosome,
                                        'insertion_point': cluster['break_second'],
                                        'clusters': break_second_clusters,
                                        'length': break_second_clusters[1]['break_first'] - break_second_clusters[0]['break_first'],
                                        'igv_d': f'{first_chromosome}:{break_second_clusters[0]["break_first"]}-{break_second_clusters[1]["break_first"]}',
                                        'igv_a': f'{second_chromosome}:{cluster["break_second"]}',
                                        'info': 'dup-down first-right second-left clipped',
                                    }
                                    self.cluster_dict[first_chromosome][second_chromosome][break_second_clusters[0]['i']].setdefault('sv', []).append(f'DUP_{first_chromosome}_to_{second_chromosome}')
                                    self.cluster_dict[first_chromosome][second_chromosome][break_second_clusters[1]['i']].setdefault('sv', []).append(f'DUP_{first_chromosome}_to_{second_chromosome}')
                                    self._sv_distant_duplications.append(duplication)


    def call_inversion_svs(self):
        """
        Returns SVs that involve _sv_inversions. These include:
        - Clean _sv_inversions
        - Tandem duplications with inversion
        - Distant duplications with _sv_inversions
        :return:
        """

        for first_chromosome in self.cluster_dict:

            for second_chromosome in self.cluster_dict[first_chromosome]:
                clusters = self.cluster_dict[first_chromosome][second_chromosome]

                first_chrom_sorted_clusters = sorted(clusters, key=lambda x: x['break_first'])

                for cluster in first_chrom_sorted_clusters:

                    if len(cluster['alignments']) >= self.min_cluster_reads:

                        if self.only_left_clipped_cluster(cluster) or self.only_right_clipped_cluster(cluster):  # Only seen with inversion SVs

                            break_first_clusters = self.get_shared_breakpoint_clusters(clusters, cluster, 'break_first')
                            break_second_clusters = self.get_shared_breakpoint_clusters(clusters, cluster, 'break_second')

                            for break_first_cluster in break_first_clusters:
                                if break_first_cluster['i'] != cluster['i']: # Two different clusters with shared first breakpoint_1

                                    if (self.only_left_clipped_cluster(cluster) and self.only_right_clipped_cluster(break_first_cluster)) or (self.only_right_clipped_cluster(cluster) and self.only_left_clipped_cluster(break_first_cluster)):  # Hallmark of inversion SVs
                                        if (
                                            is_breakpoint_shared(break_first_cluster['break_second'],cluster['break_second'], self.max_breakpoint_wobble) and
                                            first_chromosome == second_chromosome and
                                            self.cluster_alignments_are_bidirectional(break_first_cluster) and
                                            self.cluster_alignments_are_bidirectional(cluster)
                                        ):
                                            inversion = {
                                                'chromosome': first_chromosome,
                                                'begin': cluster['break_first'],
                                                'end': cluster['break_second'],
                                                'length': cluster['break_second'] - cluster['break_first'],
                                                'clusters': [cluster, break_first_cluster],
                                                'igv': f"{first_chromosome}:{cluster['break_first']}-{cluster['break_second']}"
                                            }
                                            if len(break_first_clusters) > 2 or len(break_second_clusters) > 2:  # Todo add to ambiguous clusters
                                                inversion['multi_cluster_bp'] = 'Several clusters detected at at least one of the breakpoints'
                                            self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'INV_{first_chromosome}')
                                            self.cluster_dict[first_chromosome][second_chromosome][break_first_cluster['i']].setdefault('sv', []).append(f'INV_{first_chromosome}')
                                            self._sv_inversions.append(inversion)

                                        elif self.is_duplex_cluster(break_first_cluster):
                                            if (
                                                first_chromosome == second_chromosome and
                                                self.only_left_clipped_cluster(break_first_cluster) and
                                                self.only_right_clipped_cluster(cluster) and
                                                self.cluster_alignments_are_bidirectional(break_first_cluster) and
                                                self.cluster_alignments_are_bidirectional(cluster)
                                            ): # Tandem upstream dup with inversion (head-to-head):
                                                inverse_tandem_dup = {
                                                    'chromosome': first_chromosome,
                                                    'begin': cluster['break_first'],
                                                    'end': cluster['break_second'],
                                                    'length': cluster['break_second'] - cluster['break_first'],
                                                    'clusters': [cluster, break_first_cluster],
                                                    'orient': 'head-to-head',
                                                    'igv': f"{first_chromosome}:{cluster['break_first']}-{cluster['break_second']}",
                                                }
                                                self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'INV_DUP_{first_chromosome}')
                                                self.cluster_dict[first_chromosome][second_chromosome][break_first_cluster['i']].setdefault('sv', []).append(f'INV_DUP_{first_chromosome}')
                                                self._sv_inverted_tandem_duplications.append(inverse_tandem_dup)

                                        elif self.is_duplex_cluster(cluster):  # This can be both upstream and downstream tandem duplication with inversion
                                            if (
                                                    first_chromosome == second_chromosome and
                                                    self.only_left_clipped_cluster(cluster) and
                                                    self.only_right_clipped_cluster(break_first_cluster) and
                                                    self.cluster_alignments_are_bidirectional(break_first_cluster) and
                                                    self.cluster_alignments_are_bidirectional(cluster)
                                            ):  # Tandem upstream dup with inversion (head-to-head):
                                                inverse_tandem_dup = {
                                                    'chromosome': first_chromosome,
                                                    'begin': break_first_cluster['break_first'],
                                                    'end': break_first_cluster['break_second'],
                                                    'length': break_first_cluster['break_second'] - break_first_cluster['break_first'],
                                                    'clusters': [cluster, break_first_cluster],
                                                    'orient': 'head-to-head',
                                                    'igv': f"{first_chromosome}:{break_first_cluster    ['break_first']}-{cluster['break_second']}",
                                                }
                                                self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'INV_DUP_{first_chromosome}')
                                                self.cluster_dict[first_chromosome][second_chromosome][break_first_cluster['i']].setdefault('sv', []).append(f'INV_DUP_{first_chromosome}')
                                                self._sv_inverted_tandem_duplications.append(inverse_tandem_dup)

                                            elif (
                                                first_chromosome == second_chromosome and
                                                self.only_right_clipped_cluster(cluster) and
                                                self.only_left_clipped_cluster(break_first_cluster) and
                                                self.cluster_alignments_are_bidirectional(break_first_cluster) and
                                                self.cluster_alignments_are_bidirectional(cluster)
                                            ): # Tandem downstream dup with inversion (tail-to-tail)
                                                inverse_tandem_dup = {
                                                    'chromosome': first_chromosome,
                                                    'begin': break_first_cluster['break_first'],
                                                    'end': break_first_cluster['break_second'],
                                                    'length': break_first_cluster['break_second'] - break_first_cluster['break_first'],
                                                    'clusters': [cluster, break_first_cluster],
                                                    'orient': 'tail-totail',
                                                    'igv': f"{first_chromosome}:{break_first_cluster['break_first']}-{cluster['break_second']}",
                                                }
                                                self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'INV_DUP_{first_chromosome}')
                                                self.cluster_dict[first_chromosome][second_chromosome][break_first_cluster['i']].setdefault('sv', []).append(f'INV_DUP_{first_chromosome}')
                                                self._sv_inverted_tandem_duplications.append(inverse_tandem_dup)
                                        else:  # Signature of inverted upstream dup.

                                            sorted_sv = sorted([cluster, break_first_cluster], key=lambda x:x['break_second'])
                                            if (
                                                self.only_left_clipped_cluster(sorted_sv[0]) and
                                                self.only_right_clipped_cluster(sorted_sv[1]) and
                                                not is_breakpoint_shared(cluster['break_second'], break_first_cluster['break_second']) and
                                                self.cluster_alignments_are_bidirectional(cluster) and
                                                self.cluster_alignments_are_bidirectional(break_first_cluster)
                                            ):
                                                inverse_dup = {  # Upstream inverted insertion
                                                    'donor': second_chromosome,
                                                    'acceptor': first_chromosome,
                                                    'begin': sorted_sv[0]['break_second'],
                                                    'end': sorted_sv[1]['break_second'],
                                                    'length': sorted_sv[1]['break_second'] - sorted_sv[0]['break_second'],
                                                    'insertion_point': cluster['break_first'],
                                                    'clusters': sorted_sv,
                                                    'donor-igv': f"{second_chromosome}:{sorted_sv[0]['break_second']}-{sorted_sv[1]['break_second']}",
                                                    'acceptor-igv': f"{first_chromosome}:{cluster['break_first']}"
                                                }
                                                self.cluster_dict[first_chromosome][second_chromosome][sorted_sv[0]['i']].setdefault('sv', []).append(f'INV_DUP_{first_chromosome}_{second_chromosome}')
                                                self.cluster_dict[first_chromosome][second_chromosome][sorted_sv[1]['i']].setdefault('sv', []).append(f'INV_DUP_{first_chromosome}_{second_chromosome}')
                                                self._sv_inverted_distant_duplications.append(inverse_dup)

                            for break_second_cluster in break_second_clusters:
                                if break_second_cluster['i'] != cluster['i']: # Two clusters with shared second breakpoint_1
                                    if (self.only_left_clipped_cluster(cluster) and self.only_right_clipped_cluster(break_second_cluster)) or (self.only_right_clipped_cluster(cluster) and self.only_left_clipped_cluster(break_second_cluster)):
                                        if self.is_duplex_cluster(break_second_cluster):  # Tandem downstream dup with inversion (tail-to-tail):
                                            if first_chromosome == second_chromosome:
                                                inverse_tandem_dup = {
                                                    'chromosome': first_chromosome,
                                                    'begin': cluster['break_first'],
                                                    'end': cluster['break_second'],
                                                    'length': cluster['break_second'] - cluster['break_first'],
                                                    'clusters': [cluster, break_second_cluster],
                                                    'orient': 'tail-to-tail',
                                                    'igv': f"{first_chromosome}:{cluster['break_first']}-{cluster['break_second']}",
                                                }
                                                self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'INV_DUP_{first_chromosome}')
                                                self.cluster_dict[first_chromosome][second_chromosome][break_second_cluster['i']].setdefault('sv', []).append(f'INV_DUP_{first_chromosome}')
                                                self._sv_inverted_tandem_duplications.append(inverse_tandem_dup)
                                        else:
                                            sorted_sv = sorted([cluster, break_second_cluster], key=lambda x: x['break_first'])
                                            if (
                                                    self.only_left_clipped_cluster(sorted_sv[0]) and
                                                    self.only_right_clipped_cluster(sorted_sv[1]) and
                                                    not is_breakpoint_shared(cluster['break_first'],break_first_cluster['break_first']) and
                                                    self.cluster_alignments_are_bidirectional(cluster) and
                                                    self.cluster_alignments_are_bidirectional(break_first_cluster)
                                            ):
                                                inverse_dup = {  # Downstream inverted insertion
                                                    'donor': first_chromosome,
                                                    'acceptor': second_chromosome,
                                                    'begin': sorted_sv[0]['break_first'],
                                                    'end': sorted_sv[1]['break_first'],
                                                    'length': sorted_sv[1]['break_first'] - sorted_sv[0]['break_first'],
                                                    'insertion_point': cluster['break_second'],
                                                    'clusters': sorted_sv,
                                                    'donor-igv': f"{first_chromosome}:{sorted_sv[0]['break_first']}-{sorted_sv[1]['break_first']}",
                                                    'acceptor-igv': f"{second_chromosome}:{cluster['break_second']}"
                                                }
                                                self.cluster_dict[first_chromosome][second_chromosome][sorted_sv[0]['i']].setdefault('sv', []).append(f'DIST_INV_DUP_{first_chromosome}_{second_chromosome}')
                                                self.cluster_dict[first_chromosome][second_chromosome][sorted_sv[1]['i']].setdefault('sv', []).append(f'DIST_INV_DUP_{first_chromosome}_{second_chromosome}')
                                                self._sv_inverted_distant_duplications.append(inverse_dup)


    def call_large_chromosomal_svs(self):
        """
        Calls large chromosomal rearrangements such as reciprocal translocations, chromosome arm fusions, and ring
        chromosomes.
        :return:
        """

        for first_chromosome in self.cluster_dict:
            for second_chromosome in self.cluster_dict[first_chromosome]:
                clusters = self.cluster_dict[first_chromosome][second_chromosome]

                first_chrom_sorted_clusters = sorted(clusters, key=lambda x: x['break_first'])

                for cluster in first_chrom_sorted_clusters:
                    if len(cluster['alignments']) >= self.min_cluster_reads:

                        if first_chromosome == second_chromosome:
                            if self.is_duplex_cluster(cluster):  # Todo: Potentially add median alignment score as quality measure, here
                                break_first_clusters = self.get_shared_breakpoint_clusters(clusters, cluster, 'break_first')
                                if len(break_first_clusters) == 1:
                                    arm_fusion = {
                                        'chromosome': first_chromosome,
                                        'fusion_point': cluster['break_first'],
                                        'clusters': cluster,
                                        'length': 'Nan',
                                        'igv': f"{first_chromosome}:{cluster['break_first']}"
                                    }
                                    self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'ARM_FUS{first_chromosome}')
                                    self._sv_chromosome_arm_fusions.append(arm_fusion)

                            elif self.outward_clipping_cluster(cluster):  # Can also be a large tandem duplication
                                    length = cluster['break_second'] - cluster['break_first']  # Todo: Potentially add median alignment score as quality measure, here or at cluster creation
                                    if length > 30000000:  # Above 30 million. Duplications are typically not this size # Todo: This typically shows up with low alignment score clusters. Condition on that: If all reads have a bad score at one location, it is probably bad
                                        ring_chromosome = {
                                            'chromosome': first_chromosome,
                                            'breakpoint_1': cluster['break_first'],
                                            'breakpoint_2': cluster['break_second'],
                                            'clusters': cluster,
                                            'length': length,
                                            'igv_1': f"{first_chromosome}:{cluster['break_first']}",
                                            'igv_2': f"{first_chromosome}:{cluster['break_second']}"
                                        }
                                        self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'RING_CHROM{first_chromosome}')
                                        self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'RING_CHROM{first_chromosome}')
                                        self._sv_ring_chromosomes.append(ring_chromosome)
                            else:
                                break_first_clusters = self.get_shared_breakpoint_clusters(clusters, cluster, 'break_first')
                                for break_first_cluster in break_first_clusters:
                                    if break_first_cluster['i'] != cluster['i']:  # Two clusters with shared first breakpoint_1
                                        if self.only_left_clipped_cluster(break_first_cluster) and self.only_left_clipped_cluster(cluster) and is_breakpoint_shared(break_first_cluster['break_second'], cluster['break_second'], self.max_breakpoint_wobble):

                                            # Potential uneven chromosome arm fusion
                                            if first_chromosome == second_chromosome:
                                                arm_fusion = {
                                                    'chromosome': first_chromosome,
                                                    'fusion_point_1': cluster['break_first'],
                                                    'fusion_point_2': cluster['break_second'],
                                                    'clusters': [cluster, break_first_cluster],
                                                    'length': cluster['break_second'] - cluster['break_first'],
                                                    'igv': f"{first_chromosome}:{cluster['break_first']}-{cluster['break_second']}"
                                                }
                                                if 'sv' in cluster:
                                                    print('Ambiguous cluster - uneven chromosome arm fusion')
                                                    arm_fusion['previous sv'] = cluster['sv']
                                                self._sv_uneven_chromosome_arm_fusions.append(arm_fusion)
                                                self.cluster_dict[first_chromosome][second_chromosome][break_first_cluster['i']].setdefault('sv', []).append(f'ARM_FUS_{first_chromosome}')
                                            else:
                                                print(
                                                    'Found uneven arm fusion signature, but breakpoints are on different chromosomes')

                        else: # Translocations in this section

                            if self.bidirectional_clipping_cluster(cluster):  # Starting with the most likely: Balanced translocation. Signature will also be seen with 3:1 malsegregation with interchange trisomy
                                malsegregated_translocation_data = self.get_translocation_with_malsegregation_data(cluster, first_chromosome, second_chromosome)
                                if malsegregated_translocation_data:

                                    if self.cluster_reads_have_inversion_signature(cluster):  # Signature for inversions or pq translocations
                                        self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'PQ_MALSEG_TRANS_{first_chromosome}_{second_chromosome}')
                                        self._sv_pq_translocations_with_malsegregation.append(malsegregated_translocation_data)

                                    elif self.cluster_alignments_are_unidirectional(cluster):
                                        self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'QQ_MALSEG_TRANS_{first_chromosome}_{second_chromosome}')
                                        self._sv_qq_translocations_with_malsegregation.append(malsegregated_translocation_data)

                                else:  # We assume that it is a balanced translocation if the coverage profile does not indicate malsegregation - because the clipping pattern is unique for translocations
                                    translocation = {
                                        'chromosome_1': first_chromosome,
                                        'chromosome_2': second_chromosome,
                                        'breakpoint_1': cluster['break_first'],
                                        'breakpoint_2': cluster['break_second'],
                                        'clusters': cluster,
                                        'length': 'Nan',
                                        'igv_1': f"{first_chromosome}:{cluster['break_first']}",
                                        'igv_2': f"{second_chromosome}:{cluster['break_second']}",
                                    }

                                    if self.cluster_reads_have_inversion_signature(cluster):
                                        self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault( 'sv', []).append(f'PQ_BAL_TRANS_{first_chromosome}_{second_chromosome}')
                                        self._sv_pq_balanced_translocations.append(translocation)

                                    elif self.cluster_alignments_are_unidirectional(cluster):
                                        self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'QQ_BAL_TRANS_{first_chromosome}_{second_chromosome}')
                                        self._sv_qq_balanced_translocations.append(translocation)
                            else:
                                break_first_clusters = self.get_shared_breakpoint_clusters(clusters, cluster, 'break_first')
                                break_second_clusters = self.get_shared_breakpoint_clusters(clusters, cluster, 'break_second')

                                if len(break_first_clusters) == 1 and len(break_second_clusters) == 1: # Potential unbalanced translocation or malsegregation or balanced translocation with double deletion or malsegregation with double deletion
                                    unbalanced_translocation_data = self.get_unbalanced_translocation_data(cluster, first_chromosome, second_chromosome)  # Deletions here are invisible
                                    if unbalanced_translocation_data:

                                        if self.cluster_reads_have_inversion_signature(cluster):
                                            self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'PQ_UNBAL_TRANS_{first_chromosome}_{second_chromosome}')
                                            self._sv_pq_unbalanced_translocations.append(unbalanced_translocation_data)

                                        elif self.cluster_alignments_are_unidirectional(cluster):
                                            self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'QQ_UNBAL_TRANS_{first_chromosome}_{second_chromosome}')
                                            self._sv_qq_unbalanced_translocations.append(unbalanced_translocation_data)
                                    else:
                                        translocation_with_dual_deletion_data = self.get_balanced_translocation_with_dual_deletion_data(cluster, clusters, first_chromosome, second_chromosome)
                                        if translocation_with_dual_deletion_data:

                                            # testing for malsegregation before concluding
                                            malsegregated_translocation_with_dual_deletion_data = self.get_malsegregated_translocation_with_dual_deletion_data(cluster, clusters, first_chromosome, second_chromosome)
                                            if malsegregated_translocation_with_dual_deletion_data:
                                                sv_clusters = malsegregated_translocation_with_dual_deletion_data['clusters']
                                                if self.cluster_reads_have_inversion_signature(cluster):
                                                    self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[0]['i']].setdefault('sv', []).append(f'PQ_MALSEG_TRANS_DOUBLE_DEL_{first_chromosome}_{second_chromosome}')
                                                    self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[1]['i']].setdefault('sv', []).append(f'PQ_MALSEG_TRANS_DOUBLE_DEL_{first_chromosome}_{second_chromosome}')
                                                    self._sv_pq_translocations_with_deletion_and_malsegregation.append(malsegregated_translocation_with_dual_deletion_data)
                                                elif self.cluster_alignments_are_unidirectional(cluster):
                                                    self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[0]['i']].setdefault('sv', []).append(f'QQ_MALSEG_TRANS_DOUBLE_DEL_{first_chromosome}_{second_chromosome}')
                                                    self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[1]['i']].setdefault('sv', []).append( f'QQ_MALSEG_TRANS_DOUBLE_DEL_{first_chromosome}_{second_chromosome}')
                                                    self._sv_qq_translocations_with_deletion_and_malsegregation.append(malsegregated_translocation_with_dual_deletion_data)

                                            else:
                                                sv_clusters = translocation_with_dual_deletion_data['clusters']
                                                if self.cluster_reads_have_inversion_signature(cluster):
                                                    self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[0]['i']].setdefault('sv', []).append(f'PQ_BAL_TRANS_DOUBLE_DEL_{first_chromosome}_{second_chromosome}')
                                                    self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[1]['i']].setdefault('sv', []).append(f'PQ_BAL_TRANS_DOUBLE_DEL_{first_chromosome}_{second_chromosome}')
                                                    self._sv_pq_balanced_translocations_with_deletion.append(translocation_with_dual_deletion_data)
                                                elif self.cluster_alignments_are_unidirectional(cluster):
                                                    self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[0]['i']].setdefault('sv', []).append(f'QQ_BAL_TRANS_DOUBLE_DEL_{first_chromosome}_{second_chromosome}')
                                                    self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[1]['i']].setdefault('sv', []).append(f'QQ_BAL_TRANS_DOUBLE_DEL_{first_chromosome}_{second_chromosome}')
                                                    self._sv_qq_balanced_translocations_with_deletion.append(translocation_with_dual_deletion_data)
                                        else:
                                            translocation_with_malsegregation_data = self.get_translocation_with_malsegregation_data(cluster, first_chromosome, second_chromosome)  # Placing this at the bottom because it may test True with double deletions because the coverage measurements regions are much larger than deletions
                                            if translocation_with_malsegregation_data:
                                                if self.cluster_reads_have_inversion_signature(cluster):
                                                    self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'PQ_MALSEG_TRANS_{first_chromosome}_{second_chromosome}')
                                                    self._sv_pq_translocations_with_malsegregation.append(translocation_with_malsegregation_data)
                                                elif self.cluster_alignments_are_unidirectional(cluster):
                                                    self.cluster_dict[first_chromosome][second_chromosome][cluster['i']].setdefault('sv', []).append(f'QQ_MALSEG_TRANS_{first_chromosome}_{second_chromosome}')
                                                    self._sv_qq_translocations_with_malsegregation.append(translocation_with_malsegregation_data)

                                elif (len(break_first_clusters) == 2 and len(break_second_clusters) == 1) or (len(break_first_clusters) == 1 and len(break_second_clusters) == 2): # Potential  balanced translocation with single deletion or malsegregation with single deletion

                                    balanced_translocation_with_deletion_data = self.get_balanced_tranlocation_with_deletion_data(break_first_clusters, break_second_clusters, first_chromosome, second_chromosome)
                                    if balanced_translocation_with_deletion_data:
                                        sv_clusters = balanced_translocation_with_deletion_data['clusters']
                                        if self.cluster_reads_have_inversion_signature(cluster):
                                            self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[0]['i']].setdefault('sv', []).append(f'PQ_BAL_TRANS_DEL{first_chromosome}_{second_chromosome}')
                                            self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[1]['i']].setdefault('sv', []).append(f'PQ_BAL_TRANS_DEL{first_chromosome}_{second_chromosome}')
                                            self._sv_pq_balanced_translocations_with_deletion.append(balanced_translocation_with_deletion_data)
                                        elif self.cluster_alignments_are_unidirectional(cluster):
                                            self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[0]['i']].setdefault('sv', []).append(f'QQ_BAL_TRANS_DEL{first_chromosome}_{second_chromosome}')
                                            self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[1]['i']].setdefault('sv', []).append(f'QQ_BAL_TRANS_DEL{first_chromosome}_{second_chromosome}')
                                            self._sv_qq_balanced_translocations_with_deletion.append(balanced_translocation_with_deletion_data)

                                    else:
                                        translocation_with_malsegregation_and_deletion_data = self.get_translocation_with_malsegregation_data_and_deletion(break_first_clusters, break_second_clusters, first_chromosome, second_chromosome)
                                        if translocation_with_malsegregation_and_deletion_data:
                                            sv_clusters = translocation_with_malsegregation_and_deletion_data['clusters']
                                            if self.cluster_reads_have_inversion_signature(cluster):
                                                self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[0]['i']].setdefault('sv', []).append(f'PQ_MALSEG_TRANS_DEL{first_chromosome}_{second_chromosome}')
                                                self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[1]['i']].setdefault('sv', []).append(f'PQ_MALSEG_TRANS_DEL{first_chromosome}_{second_chromosome}')
                                                self._sv_pq_translocations_with_deletion_and_malsegregation.append(translocation_with_malsegregation_and_deletion_data)
                                            elif self.cluster_alignments_are_unidirectional(cluster):
                                                self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[0]['i']].setdefault('sv', []).append(f'QQ_MALSEG_TRANS_DEL{first_chromosome}_{second_chromosome}')
                                                self.cluster_dict[first_chromosome][second_chromosome][sv_clusters[1]['i']].setdefault('sv', []).append(f'QQ_MALSEG_TRANS_DEL{first_chromosome}_{second_chromosome}')
                                                self._sv_qq_translocations_with_deletion_and_malsegregation.append(translocation_with_malsegregation_and_deletion_data)


                                elif len(break_first_clusters) > 1 and len(break_second_clusters) > 1:
                                    self.ambiguous_breakpoints.append({
                                        'signature': 'first>1;second>1',
                                        'chromosome_1': first_chromosome,
                                        'chromosome_2': second_chromosome,
                                        'break_first': cluster['break_first'],
                                        'break_second': cluster['break_second'],
                                        'break_first_clusters': break_first_clusters,
                                        'break_second_clusters': break_second_clusters,
                                    })
                                elif len(break_first_clusters) > 2:
                                    self.ambiguous_breakpoints.append({
                                        'signature': 'first>2',
                                        'chromosome_1': first_chromosome,
                                        'chromosome_2': second_chromosome,
                                        'break_first': cluster['break_first'],
                                        'break_first_clusters': break_first_clusters,
                                    })
                                elif len(break_second_clusters) > 2:
                                    self.ambiguous_breakpoints.append({
                                        'signature': 'second>2',
                                        'chromosome_1': first_chromosome,
                                        'chromosome_2': second_chromosome,
                                        'break_second': cluster['break_second'],
                                        'break_second_clusters': break_second_clusters
                                    })


    def do_sv_scan(self):
        """
        Does an end-to-end scan-and-trace of clusters found at each chromosome. Rather than independent calling of each SV via signature recognition, this approach may be better at handling complex
        rearrangements.

        For alignments, each element has the structure:
        {
            'first': [1844, 7207, 7207, '+', 'l'],
            'second': [26056879, 26061816, 26056879, '-', 'l'],
            'read': 'cd9436f1-8176-4292-9633-669a3acfaaba',
        }
        :return:
        """

        # # delete. Structure check while developing
        # for first_chromosome in self.cluster_dict:
        #
        #     for second_chromosome in self.cluster_dict[first_chromosome]:
        #         clusters = self.cluster_dict[first_chromosome][second_chromosome]
        #
        #         first_chrom_sorted_clusters = sorted(clusters, key=lambda x: x['break_first'])
        #         first_chrom_clusters = iter(first_chrom_sorted_clusters)
        #         previous_cluster = first_chrom_clusters.__next__()
        #         for cluster in first_chrom_clusters:
        #             if is_breakpoint_shared(previous_cluster['break_first'], cluster['break_first']):
        #                 print('shared_break_first_cluster in self.cluster_dict')
        #                 previous_cluster = cluster
        #
        #         second_chrom_sorted_clusters = sorted(clusters, key=lambda x: x['break_second'])
        #         second_chrom_clusters = iter(second_chrom_sorted_clusters)
        #         previous_cluster = second_chrom_clusters.__next__()
        #         for cluster in second_chrom_clusters:
        #             if is_breakpoint_shared(previous_cluster['break_second'], cluster['break_second']):
        #                 print('shared_break_second in self.cluster_dict')
        #                 previous_cluster = cluster
        # exit()

        # Add unique id to each cluster
        print('Reindexing clusters')
        cluster_id = 0
        for first_chromosome in self.cluster_dict:

            for second_chromosome in self.cluster_dict[first_chromosome]:
                clusters = self.cluster_dict[first_chromosome][second_chromosome]

                for cluster in clusters:
                    cluster['i'] = cluster_id
                    cluster_id += 1

        # Rearrange structure to collect all clusters for each chromosome. Change break-first to refer to that chromosome for chromosomes with lower lexical order
        chromosome_clusters = {}
        for first_chromosome in self.cluster_dict:
            for second_chromosome in self.cluster_dict[first_chromosome]:
                clusters = self.cluster_dict[first_chromosome][second_chromosome]
                for cluster in clusters:
                    cluster['second_chromosome'] = second_chromosome
                    chromosome_clusters.setdefault(first_chromosome, []).append(cluster)

                    if first_chromosome != second_chromosome:
                        second_chrom_alignments = []
                        for alignment in cluster['alignments']:
                            second_chrom_alignments.append({
                                'first': alignment['second'],
                                'second': alignment['first'],
                                'read': alignment['read'],
                            })

                        break_first = cluster['break_second']
                        break_second = cluster['break_first']

                        cluster['break_first'] = break_first
                        cluster['break_second'] = break_second
                        cluster['alignments'] = second_chrom_alignments
                        cluster['second_chromosome'] = first_chromosome

                        chromosome_clusters.setdefault(second_chromosome, []).append(cluster)

        # Sort and index
        id_to_index = {}
        for chromosome in chromosome_clusters:
            clusters = chromosome_clusters[chromosome]
            sorted_clusters = sorted(clusters, key=lambda x:x['break_first'])

            for cluster_index, cluster in enumerate(sorted_clusters):
                id_to_index.setdefault(chromosome, {})[cluster['i']] = cluster_index
            chromosome_clusters[chromosome] = sorted_clusters
        for chromosome in chromosome_clusters:
            for cluster in chromosome_clusters[chromosome]:
                cluster['second_chromosome_index'] = id_to_index[cluster['second_chromosome']][cluster['i']]

        # Register shared-breakpoint clusters
        for chromosome in chromosome_clusters:
            for index in range(1, len(chromosome_clusters[chromosome])):
                if is_breakpoint_shared(chromosome_clusters[chromosome][index - 1]['break_first'], chromosome_clusters[chromosome][index]['break_first'], self.max_breakpoint_wobble):
                    chromosome_clusters[chromosome][index - 1]['shared_break_cluster'] = True
                    chromosome_clusters[chromosome][index]['shared_break_cluster'] = True
                    print(f'Shared break on {chromosome} at {chromosome_clusters[chromosome][index]["break_first"]}')
        exit()


        # Scan and assign SVs
        SVs = {}
        for chromosome in chromosome_clusters:
            print('Calling structural variants')
            for cluster in chromosome_clusters[chromosome]:

                if cluster['second_chromosome'] == chromosome:
                    # Either cluster is left or right-clipped
                    if self.only_left_clipped_alignments(cluster, 'break_first'):
                        pass
                    elif self.only_right_clipped_alignments(cluster, 'break_first'):
                        pass

                else: # Interchromosomal SV
                    # Either cluster is left or right-clipped
                    if self.only_left_clipped_alignments(cluster, 'break_first'):
                        pass
                    elif self.only_right_clipped_alignments(cluster, 'break_first'):
                        pass



    def register_unplaced(self):
        for first_chromosome in self.cluster_dict:

            for second_chromosome in self.cluster_dict[first_chromosome]:
                clusters = self.cluster_dict[first_chromosome][second_chromosome]

                for cluster in clusters:
                    if not 'sv' in cluster:
                        self.unplaced.setdefault(first_chromosome, {}).setdefault(second_chromosome, []).append(cluster)


    def sanitize_cluster_dict(self):
        """
        Sanitizes the clusters dictionary by throwing out clusters that do not contain the minimal number of reads
        an adds an index number to each cluster.
        :return:
        """
        sanitized_cluster_dict = {}
        for first_chromosome in self.cluster_dict:

            for second_chromosome in self.cluster_dict[first_chromosome]:
                clusters = self.cluster_dict[first_chromosome][second_chromosome]
                index = 0
                for cluster in clusters:
                    if len(cluster['alignments']) >= self.min_cluster_reads:
                        cluster['i'] = index
                        index += 1
                        sanitized_cluster_dict.setdefault(first_chromosome, {}).setdefault(second_chromosome, []).append(cluster)

        self.cluster_dict = sanitized_cluster_dict



    def __init__(self, cluster_dict, bam_data, min_cluster_reads=3, max_breakpoint_wobble=20, annotation_file=None):
        """
        :param cluster_dict:
        :param min_cluster_reads:
        :param annotation_file:
        """
        self.cluster_dict = cluster_dict
        self.bam_data = bam_data
        self.min_cluster_reads = min_cluster_reads
        self.max_breakpoint_wobble = max_breakpoint_wobble
        self.annotation_file = annotation_file
        self.sanitize_cluster_dict()

        # SV attributes - MUST ALL BE PREFIXED WITH _sv_
        self._sv_distant_duplications = []
        self._sv_inverted_distant_duplications = []
        self._sv_tandem_duplications = []
        self._sv_inverted_tandem_duplications = []
        self._sv_inversions = []
        self._sv_deletions = []
        
        self._sv_chromosome_arm_fusions = []
        self._sv_uneven_chromosome_arm_fusions = []
        self._sv_ring_chromosomes = []
        
        self._sv_qq_balanced_translocations = []
        self._sv_qq_balanced_translocations_with_deletion = []
        self._sv_qq_unbalanced_translocations = []
        self._sv_qq_unbalanced_translocations_with_deletion = []
        self._sv_qq_translocations_with_malsegregation = []
        self._sv_qq_translocations_with_deletion_and_malsegregation = []

        self._sv_pq_balanced_translocations = []
        self._sv_pq_balanced_translocations_with_deletion = []
        self._sv_pq_unbalanced_translocations = []
        self._sv_pq_unbalanced_translocations_with_deletion = []
        self._sv_pq_translocations_with_malsegregation = []
        self._sv_pq_translocations_with_deletion_and_malsegregation = []

        self.ambiguous_breakpoints = []
        self.unplaced = {}

        print('Inferring Del/Dup variants with no inversions')
        self.call_non_inversion_svs()
        print('Inferring Del/Dup variants with inversions')
        self.call_inversion_svs()
        print('Inferring large chromosomal structural variants')
        self.call_large_chromosomal_svs()
        self.register_unplaced()

        # self.do_sv_scan()



def infer_svs_from_clusters(cluster_dict, bam_data, output_folder, sample_suffix, min_cluster_reads=3):
    """
    Infers SVs from the pairwise clusters by calling the StuctVarFinder class.
    Secondly, it outputs a json file for each type of SV.
    Lastly, it outputs jsons for any unplaced clusters with a json file for each number of reads in the unplaced
    clusters.
    :param cluster_dict:
    :param bam_data:
    :param output_folder:
    :param sample_suffix:
    :param min_cluster_reads:
    :return:
    """

    sv_finder = StuctVarFinder(cluster_dict, bam_data, min_cluster_reads, 20)

    svs = [sv for sv in vars(sv_finder) if sv.startswith('_sv_')]

    os.makedirs(os.path.join(output_folder, 'sv_signatures'), exist_ok=True)

    sv_stats = {
        'Structural variant': [],
        'Count': []
    }
    for sv in svs:

        sv_calls = getattr(sv_finder, sv)
        sv_calls = sorted(sv_calls, key=lambda x: 0 if type(x['length']) == str else x['length'], reverse=True)
        per_read_number_svs = {}

        if len(sv_calls) > 0:
            for entry in sv_calls:
                if isinstance(entry['clusters'], dict):
                    per_read_number_svs.setdefault(len(entry['clusters']['alignments']), []).append(entry)
                else:
                    per_read_number_svs.setdefault(min(len(entry['clusters'][0]['alignments']), len(entry['clusters'][1]['alignments'])), []).append(entry)

            with open(os.path.join(output_folder, 'sv_signatures', f'{sample_suffix}{sv}.json'), 'w') as json_file:
                json_file.write(json.dumps(per_read_number_svs))

            # Save stats for pandas array
            sv_stats['Structural variant'].append(sv)
            sv_stats['Count'].append(len(sv_calls))


    # Write unplaced clusters and ambiguous breakpoints to json files
    unplaced = sv_finder.unplaced
    with open(os.path.join(output_folder, 'sv_signatures', 'unplaced_clusters.json'), 'w') as json_file:
        json_file.write(json.dumps(unplaced))

    ambiguous_breakpoints = sv_finder.ambiguous_breakpoints
    with open(os.path.join(output_folder, 'sv_signatures', 'ambiguous_breakpoints.json'), 'w') as json_file:
        json_file.write(json.dumps(ambiguous_breakpoints))

    generate_cluster_read_counts_jsons(unplaced, os.path.join(output_folder, 'unplaced_clusters_per_nr_of_supporting_reads'))

    sv_stats_df = pandas.DataFrame.from_dict(sv_stats)
    sv_stats_html_table = sv_stats_df.to_html(col_space=200, index=False)

    placement_stats = {}
    for first_chromosome in cluster_dict:

        for second_chromosome in cluster_dict[first_chromosome]:
            clusters = cluster_dict[first_chromosome][second_chromosome]
            for cluster in clusters:
                if len(cluster['alignments']) >= min_cluster_reads:
                    if 'sv' in cluster:
                        placement_stats.setdefault(len(cluster['alignments']), [0, 0])[0] += 1
                    else:
                        placement_stats.setdefault(len(cluster['alignments']), [0, 0])[1] += 1
    for number_of_reads in placement_stats:
        placement_stats[number_of_reads] = placement_stats[number_of_reads][0] / (placement_stats[number_of_reads][0] + placement_stats[number_of_reads][1])

    sorted_placement_stats = {}
    for key in sorted(placement_stats):
        sorted_placement_stats[key] = placement_stats[key]

    placement_dict = {
        'Reads per cluster': [],
        'Proportion of clusters in identified SVs': []
    }
    for number_of_reads in sorted_placement_stats:
        placement_dict['Reads per cluster'].append(number_of_reads)
        placement_dict['Proportion of clusters in identified SVs'].append(sorted_placement_stats[number_of_reads])

    data_frame = pandas.DataFrame.from_dict(placement_dict)
    placement_html_table = data_frame.to_html(col_space=200, index=False)

    return sv_stats_html_table, placement_html_table



class HTML_Element:

    def get_element(self):
        elm_string = f'<{self.elm}{self.id}{self.classes}>'
        if self.content:
            elm_string += self.content
        for child in self.children:
            elm_string += '\n' + child.get_element() + '\n'
        elm_string += f'</{self.elm}>'
        return elm_string


    def add_element(self, element, content=None, id=None, classes=[]):
        self.children.append(HTML_Element(element, content, id, classes))
        return self.children[-1]


    def __init__(self, element, content=None, id=None, classes=[]):
        self.elm = element
        self.content = content
        self.id = f' id="{id}"' if id else ''
        self.classes = f' class="{" ".join(classes)}"' if classes else ''
        self.children = []


class HTML_Report:

    def get_body_inner(self):
        return '\n'.join([element.get_element() for element in self.elements])


    def get_html(self):
        html = f"""
        <html>
            <head>
                <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            </head>
            <style>
                {self.css}
            </style>
            <body>
                {self.get_body_inner()}
            </body>
        </html>
        """
        soup = beautiful_soup(html, features="html.parser")
        return soup.prettify()

    def add_element(self, element, content=None, id=None, classes=[]):
        self.elements.append(HTML_Element(element, content, id, classes))
        return self.elements[-1]

    def write_report(self):
        with open(self.file_name, 'w') as html_file:
            html_file.write(self.get_html())

    def __init__(self, file_name):
        self.file_name = file_name if file_name.endswith('.html') else file_name + '.html'
        self.elements = []
        self.css = ''


def get_sample_name(bam_file):
    """
    Extracts sample name from a bam file if it is included in the header
    :param bam_file:
    :return: str
    """
    samtools_header = subprocess.Popen(["samtools", "view", "-H", bam_file],stdout=subprocess.PIPE)
    grep = subprocess.Popen(["grep", "^@RG"],stdin=samtools_header.stdout,stdout=subprocess.PIPE)
    grep_out, grep_err = grep.communicate()

    header_line_elements = grep_out.decode().split()
    for elm in header_line_elements:
        if elm.startswith('SM'):
            return elm.split(':')[1]


# def get_ucsc_rest_data(endpoint_function, ep_params):
#     """
#     https://api.genome.ucsc.edu/getData/track?genome=hg38;track=centromeres;chrom=chr1
#     :param endpoint_function:
#     :param ep_params:
#     :return:
#     """
#     import requests, sys
#
#     server = 'https://api.genome.ucsc.edu/'
#     request_url = server + endpoint_function
#     if ep_params:
#         request_url += '?'
#         ';'.join([f'{param[0]}={param[1]}'  for param in ep_params])
#
#     return requests.get(request_url)
#
#     request_result = requests.get(request_url)
#
#
#
# def get_centromere_positions(genome):
#     import requests, sys
#
#     server = 'https://api.genome.ucsc.edu/'
#     endpoint =


def run_with_snakemake_input(log_file):
    """
    Not currently used. If needed, make it correspond to main
    :return:
    """
    import sys
    orig_stdout = sys.stdout
    orig_stderr = sys.stderr
    log = open(log_file, 'w')
    sys.stdout = log
    sys.stderr = log

    target_file = snakemake.input['BAM']
    number_of_genes = snakemake.params['number_of_genes']
    quantify_genes(target_file, number_of_genes)

    sys.stdout = orig_stdout
    sys.stderr = orig_stderr
    log.close()

def main():
    develop = False
    # annotation_file = "Q:/resources/reference_genomes/T2T/T2T-CHM13v2.0/logic_names/ALL_CAPS/T2T-CHM13v2_annotation.gtf.gz"
    #
    # get_gene_exons_annotation_data(['SDHA'], annotation_file)

    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-b', '--bam_file', help='Path to bam file', required=False)
    parser.add_argument('-d', '--bam_dir', help='Path to directory containing bam file(s).', required=False)
    parser.add_argument('-n', '--sample_name', help='Sample name. Used as suffix for output files. Ignored if the recursive flag is set.', required=False)
    parser.add_argument('-s', '--output_folder_suffix', default=None, help='Suffix to add to the output folder', required=False)
    parser.add_argument('-z', '--zip_archive', default=False, action="store_true", help='Flag. If set a zip archive of all result directories will be created in the bam dir. If a single bam file is defined it wil be placed in the folder containing it.', required=False)
    parser.add_argument('-f', '--use_file_suffix', default=False, action="store_true", help='Flag. If set, the file name is used as or added to the output_folder_suffix. Useful when the recursive flag is set if a folder holds multiple bam files.', required=False)
    parser.add_argument('-r', '--recursive', default=False, action="store_true", help='Flag. Only efficient with the -d option. A "set and forget" option that causes recursive traversal of the directory and analysis of all encountered bam files.', required=False)
    parser.add_argument('-w', '--max_breakpoint_wobble', default=20, help='Maximally allowed breakpoint_1 wobble counted as bases upstream and downstream from a bp when two breakpoints are evaluated for equality. Default = 20.', type=int, required=False)
    parser.add_argument('-q', '--minimum_alignment_score', default=20, help='Minimum alignment_score. Default = 20. Range: 0 - 60, higher/lower values will be adjusted.', type=int, required=False)
    parser.add_argument('-m', '--minimum_cluster_reads', default=3, help='Minimum cluster reads. Default = 3', type=int, required=False)
    parser.add_argument('-c', '--number_of_cores', default=1, help='Number of cores to use for the initial and most expensive bam-file querying.', type=int, required=False)
    parser.add_argument('-C', '--available_cores', help='Shows available cores and exits', action='version', version=f'The machine has {multiprocessing.cpu_count()} logical cores available.')
    parser.add_argument('-v', '--version', help='Shows script version and exits.', action='version', version=f'Version: {__version__}')

    args = parser.parse_args()
    bam_dir = args.bam_dir
    output_folder_suffix = args.output_folder_suffix
    max_breakpoint_wobble = args.max_breakpoint_wobble
    minimum_alignment_score = max(0, min(args.minimum_alignment_score, 60))
    minimum_cluster_reads = args.minimum_cluster_reads
    recursive = args.recursive
    use_file_suffix = args.use_file_suffix
    zip_archive = args.zip_archive
    bam_files = get_bamfiles(args.bam_file) if args.bam_file else get_bamfiles(bam_dir, recursive)
    number_of_cores = args.number_of_cores

    if output_folder_suffix:
        if os.sep in output_folder_suffix:
            no_trace_error(f'"{os.sep}" is not allowed in the output_folder_suffix. Please choose a different one.')
            exit()

        output_folder_base_suffix = f'{output_folder_suffix}_'
    else:
        output_folder_base_suffix = ''

    if not bam_files:
        no_trace_error('\nNo bam file found or specified. Please provide a path to a bam file or to a folder containing it. It is ok if it is in a sub-folder of the specified folder')

    output_folders = []
    for bam_file in bam_files:
        print(f'Processing file: {bam_file}')

        if use_file_suffix:
            file_name = os.path.splitext(os.path.basename(bam_file))[0]
            output_folder_suffix = output_folder_base_suffix + f'{file_name}_'

        sample_name = get_sample_name(bam_file) or args.sample_name if not recursive else None
        sample_suffix = f'{sample_name}_' if sample_name else ''  # Used as suffix for all output files

        output_folder_name = f'{output_folder_suffix}structural_variants'
        output_folder = os.path.join(os.path.dirname(bam_file), output_folder_name)
        os.makedirs(output_folder, exist_ok=True)

        try:  # Handle ubams (pysam cannot process ubams)
            bam_data = BamData(bam_file)
        except ValueError as error:
            if 'sequence' in str(error):
                print('File has no alignments - it might be a ubam file. Skipping')
                continue

        # Check bam file depth
        if bam_data.median_depth < 5:
            print(f'Median depth of BM file is {bam_data.median_depth}. Skipping')
            continue

        output_folders.append(output_folder)

        pairwise_alignment_clusters_json = os.path.join(output_folder, f'{sample_suffix}pairwise_alignment_clusters.json')
        if develop:
            chromosome_to_chromosome_reads_json = os.path.join(output_folder, f'{sample_suffix}chromosome_to_chromosome_reads.json')

            if os.path.isfile(pairwise_alignment_clusters_json):
                with open(pairwise_alignment_clusters_json, "r") as json_file:
                    pairwise_alignment_clusters = json.loads(json_file.read())
            else:
                if os.path.isfile(chromosome_to_chromosome_reads_json):
                    print('Loading chromosome_to_chromosome_reads.json')
                    with open(chromosome_to_chromosome_reads_json, "r") as json_file:
                        chromosome_to_chromosome_reads = json.loads(json_file.read())
                else:
                    chromosome_to_chromosome_reads = get_chromosome_to_chromosome_reads(bam_data, minimum_alignment_score, cores=number_of_cores)

                    if chromosome_to_chromosome_reads:
                        with open(chromosome_to_chromosome_reads_json, 'w') as json_file:
                            json_file.write(json.dumps(chromosome_to_chromosome_reads))


                pairwise_alignment_clusters = get_chromosome_to_chromosome_clusters(chromosome_to_chromosome_reads, max_breakpoint_wobble)

                if pairwise_alignment_clusters:
                    with open(pairwise_alignment_clusters_json, 'w') as json_file:
                        json_file.write(json.dumps(pairwise_alignment_clusters))

        else:
            chromosome_to_chromosome_reads = get_chromosome_to_chromosome_reads(bam_data, minimum_alignment_score, cores=number_of_cores)

            pairwise_alignment_clusters = get_chromosome_to_chromosome_clusters(chromosome_to_chromosome_reads, max_breakpoint_wobble)

        #######################################################################################
        #                       Infer SVs from pairwise_alignment_clusters                    #
        #######################################################################################


        # pairwise_alignment_clusters is the basic data format for SV signature analysis
        if pairwise_alignment_clusters:

            read_counts_table_html = generate_cluster_read_counts_jsons(pairwise_alignment_clusters, os.path.join(output_folder, 'clusters_per_nr_of_supporting_reads'))

            sv_stats_html_table, placement_html_table = infer_svs_from_clusters(pairwise_alignment_clusters, bam_data, output_folder, sample_suffix, minimum_cluster_reads)

            # Legacy output
            generate_max_cluster_count_heatmap_and_json(pairwise_alignment_clusters, bam_data.contig_names, output_folder, sample_suffix)  # Legacy output

            #######################################################################################
            #                                 HTML report                                         #
            #######################################################################################

            html_report = HTML_Report(os.path.join(output_folder, f'{sample_suffix}sv_report.html'))
            header = html_report.add_element('header')
            if sample_name:
                header.add_element('h1', f'Structural variants report for the sample: {sample_name}')
            else:
                header.add_element('h1', 'Structural variants report')

            coverage = html_report.add_element('div', id='coverage', classes=['plot_container'])
            coverage.add_element('h2', 'Genome wide coverage')
            coverage.add_element('div', bam_data.depth_box_plot_html)

            cluster_stats = html_report.add_element('div', id='cluster_stats')
            cluster_stats.add_element('h2', 'Cluster statistics')
            read_counts_table = cluster_stats.add_element('div', classes=['cluster_stats_table'])
            read_counts_table.add_element('h3', 'Overall cluster counts')
            read_counts_table.add_element('div', read_counts_table_html)

            placement_proportions_table = cluster_stats.add_element('div', classes=['cluster_stats_table'])
            placement_proportions_table.add_element('h3', 'SV identification success ratio')
            placement_proportions_table.add_element('div', placement_html_table)

            sv_stats = html_report.add_element('div', id='sv_stats')
            sv_stats.add_element('h2', 'SV statistics')
            sv_stats.add_element('div', sv_stats_html_table)

            html_report.css = """
    div.cluster_stats_table {
        display: inline-grid;
        margin-right: 200px;
    }
                """
            html_report.write_report()

    if zip_archive:
        print('Creating results archive')
        print('bam_dir', bam_dir)
        print('output_folders', output_folders)
        archive_name = args.output_folder_suffix or 'SV'
        with zipfile.ZipFile(os.path.join(bam_dir, f'{archive_name}_results_archive.zip'), 'w', zipfile.ZIP_DEFLATED) as zip_container:
            for folder in output_folders:
                zipdir(folder, zip_container)
        print('Done.')



if 'snakemake' in locals() or 'snakemake' in globals():
    log_file = str(snakemake.log[0])
    run_with_snakemake_input(log_file)

elif __name__ == "__main__":
    main()
