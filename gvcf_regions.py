#!/usr/bin/env python

import gzip

def is_header(line):
    """Check if a line is header."""

    return line.startswith('#')

def has_END(line):
    """Check if a line has the 'END' tag."""

    return line.find('END') != -1

# FIELD index
# CHROM 0, POS 1, REF 3, QUAL 5, INFO 7, FORMAT 8, sample 9

def get_END(line):
    """Extract the END position of a line from the 'END' tag under INFO."""

    INFO = line.strip().split('\t')[7]
    INFO_fields = INFO.split(';')
    for INFO_field in INFO_fields:
        if INFO_field.split('=')[0] == 'END':
            END = int(INFO_field.split('=')[1])
            return END

def get_bed_region(line):
    """Extract the bed region of a line.

    Return:
        tuple: (line_start, line_end)"""

    fields = line.strip().split('\t')
    POS = int(fields[1])
    REF = fields[3]

    # vcf: 1 index, inclusive start, inclusive end
    # bed: 0 index, inclusive start, exlusive end

    # starting postion of the line in bed index
    if POS == 0:
        line_start = 0
    else:
        line_start = POS - 1
    # ending position of the line in bed index
    if has_END(line):
        line_end = get_END(line)
    else:
        line_end = POS -1 + len(REF)
    return (line_start, line_end)

def get_GQ(line):
    """Extract the GQ (Genotype Quality) of a line."""

    fields = line.strip().split('\t')
    FORMAT, sample = fields[8], fields[9]
    FORMAT_fields = FORMAT.split(':'); sample_fields = sample.split(':')

    GQ_index = FORMAT_fields.index('GQ')
    GQ = int(sample_fields[GQ_index])
    return GQ

def is_considered(line, ignore_phrases):
    """Check if a line should be considered. A line that contains any of the
    ignore phrases is discarded. E.g., to discard CNV and ME lines of Complete Genomics
    gvcf, set ignore_phrases = ['CNV', 'INS:ME'].

    Args:
       ignore_phrases (None or non-empty list)"""

    if ignore_phrases != None:
        for ignore_phrase in ignore_phrases:
            if line.find(ignore_phrase) != -1:
                return False
    return True

def is_called(line, min_GQ, min_QUAL, pass_phrases):
    """Check if a line is considered as 'called', i.e., if its GQ is at least
    min_GQ, if its QUAL is at least min_QUAL, and if it contains any of the
    pass phrases. E.g., for platypus gvcf, set pass_phrases = ['PASS', 'REFCALL'].

    Args:
        min_GQ (None or int)
        min_QUAL (None or float)
        pass_phrases (None or non-empty list)"""

    fields = line.strip().split('\t')
    GQ_passes = (min_GQ == None or get_GQ(line) >= min_GQ)
    QUAL_passes = (min_QUAL == None or float(fields[5]) >= min_QUAL)

    line_has_pass_phrases = True
    if pass_phrases != None:
        for pass_phrase in pass_phrases:
            if line.find(pass_phrase) != -1:
                line_has_pass_phrases = True
                break
            else:
                line_has_pass_phrases = False

    return GQ_passes and QUAL_passes and line_has_pass_phrases

def gvcf_regions(gvcf, unreported_is_called, ignore_phrases,
                        min_GQ, min_QUAL, pass_phrases):
    """Generate the called regions of a gvcf as a bed file.

    Args:
        unreported_is_called (bool): whether an unreported site is considered
            as called. E.g., in Complete Genomics gvcf, an unreported site
            is hom-ref (called), whereas in freebayes/gatk/platypus
            gvcf, it is no-call.
        ignore_phrases (None or non-empty list)
        min_GQ (None or int)
        min_QUAL (None or float)
        pass_phrases (None or non-empty list)"""

    g = gzip.open(gvcf) if gvcf.endswith(".gz") else open(gvcf)

    region_CHROM = ''
    for line in g:
        # filter out header and lines with ignore phrases
        if not is_header(line) and is_considered(line, ignore_phrases):
            fields = line.strip().split('\t')
            CHROM = fields[0]
            is_line_called = is_called(line, min_GQ, min_QUAL, pass_phrases)
            (line_start, line_end) = get_bed_region(line)

            # handle cases depending on if the previous block is called,
            # if the current line is called, and if the previous block
            # and the current line have a gap

            # when chromosome changes
            if region_CHROM != CHROM:
                # when chromosome changes from one to another,
                # write the previous region
                if region_CHROM != '':
                    # we assume a chromosome ends with a string of 'N's
                    if is_previous_block_called:
                        region_end = previous_block_end
                        print('\t'.join([region_CHROM,
                            str(region_start), str(region_end)]))
                region_CHROM = CHROM
                if is_line_called:
                    region_start = line_start
                    is_previous_block_called = True
                else:
                    # we assume a chromosome starts with a string of 'N's
                    is_previous_block_called = False
            # when line and previous block are on the same chromosome and
            # they have a gap
            elif line_start > previous_block_end:
                if is_previous_block_called:
                    if unreported_is_called:
                        if not is_line_called:
                            region_end = line_start
                            print('\t'.join([region_CHROM,
                                str(region_start), str(region_end)]))
                            is_previous_block_called = False
                    else:
                        region_end = previous_block_end
                        print('\t'.join([region_CHROM,
                            str(region_start), str(region_end)]))
                        if is_line_called:
                            region_start = line_start
                        else:
                            is_previous_block_called = False
                else:
                    if unreported_is_called:
                        region_start = previous_block_end
                        if is_line_called:
                            is_previous_block_called = True
                        else:
                            region_end = line_start
                            print('\t'.join([region_CHROM,
                                str(region_start), str(region_end)]))
                    else:
                        if is_line_called:
                            region_start = line_start
                            is_previous_block_called = True
            # when line and previous block are on the same chromosome and
            # they have no gap
            else:
                if is_previous_block_called:
                    if not is_line_called:
                        region_end = line_start
                        print('\t'.join([region_CHROM,
                            str(region_start), str(region_end)]))
                        is_previous_block_called = False
                else:
                    if is_line_called:
                        region_start = line_start
                        is_previous_block_called = True
            # update previous_block_end
            previous_block_end = line_end
    # at the end of file, write the previous region
    else:
        # we assume a chromosome ends with a string of 'N's
        if is_previous_block_called:
            region_end = previous_block_end
            print('\t'.join([region_CHROM,
                str(region_start), str(region_end)]))

    g.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Output the called regions \
        of a gvcf file to stdout in bed format.')
    parser.add_argument('gvcf', metavar='GVCF', help='input gvcf file')

    parser.add_argument('--unreported_is_called', action='store_true',
        help='use this flag to treat unreported sites as called')
    parser.add_argument('--ignore_phrases', nargs='+',
        help='list of phrases considered as discarded, e.g., CNV, ME. \
        A line that contains any of the ignore phrases is discarded.')
    parser.add_argument('--min_GQ', type=int,
        help='minimum GQ (Genotype Quality) considered as called')
    parser.add_argument('--min_QUAL', type=float,
        help='minimum QUAL considered as called')
    parser.add_argument('--pass_phrases',  nargs='+',
        help='list of phrases considered as called, e.g., PASS, REFCALL. \
        A line must contain any of the pass phrases to be considered as called.')

    # presets of unreported_is_called, ignore_phrases, min_GQ, min_QUAL, pass_phrases
    complete_genomics_preset = [True, ['CNV', 'INS:ME'], None, None, ['PASS']]
    freebayes_preset = [False, None, None, None, ['PASS']]
    gatk_preset = [False, None, 5, None, None]
    platypus_preset = [False, None, None, None, ['PASS', 'REFCALL']]

    gvcf_type_help = '''type of gvcf output.
            [unreported_is_called, ignore_phrases, min_GQ, min_QUAL, pass_phrases] presets:
            complete_genomics: %s.
            freebayes: %s.
            gatk: %s.
            platypus: %s.''' % (complete_genomics_preset,
                freebayes_preset, gatk_preset, platypus_preset)

    parser.add_argument("--gvcf_type",
        choices=['complete_genomics', 'freebayes', 'gatk', 'platypus'], help=gvcf_type_help)
    args = parser.parse_args()

    par = [args.unreported_is_called, args.ignore_phrases,
                args.min_GQ, args.min_QUAL, args.pass_phrases]

    if args.gvcf_type == 'complete_genomics':
        par = complete_genomics_preset
    elif args.gvcf_type == 'freebayes':
        par = freebayes_preset
    elif args.gvcf_type == 'gatk':
        par = gatk_preset
    elif args.gvcf_type == 'platypus':
        par = platypus_preset

    all_par = [args.gvcf] + par
    gvcf_regions(*all_par)
