'''
f_line1 = 'chr22	14257	.	CTG	C	26.01	PASS	AC=2;AF=1;AN=2;DECOMPOSED;GC=53.47;HRun=0;LEN=2;MQ=58.68;NS=1;QD=8.67;TYPE=del	GT:AO:DP:GQ:PL:QA:QR:RO	1|1:2:2:19:38,6,0:40:0:0'
f_line2 = 'chr22	14273	.	C	<*>	0	PASS	AC=0;AF=nan;AN=0;DP=3;END=14282;GC=52.48;MIN=3;MQ=71.25;NS=1	DP:GQ:MIN:QA:QR	3:99:3:0:1203'
g_line1 = '1	63735	.	CCTA	C,<NON_REF>	41.73	.	BaseQRankSum=0.736;ClippingRankSum=0.736;DP=3;MLEAC=1,0;MLEAF=0.500,0.00;MQ=41.10;MQ0=0;MQRankSum=0.736;ReadPosRankSum=0.736	GT:AD:DP:GQ:PL:SB	0/1:1,2,0:3:45:81,0,45,84,51,135:0,1,0,2'
g_line2 = '1	63739	.	C	<NON_REF>	.	.	END=63739	GT:DP:GQ:MIN_DP:PL	0/0:3:0:3:0,0,0'
p_line1 = '1	30867	.	CCTCT	C	49	PASS	BRF=0.0;FR=0.9990;HP=2;HapScore=2;MGOF=20;MMLQ=16;MQ=60.0;NF=0;NR=1;PP=49;QD=53.0;SC=TCTGTGTCTCCCTCTCTCTCT;SbPval=1.0;Source=Platypus;TC=1;TCF=0;TCR=1;TR=1;WE=30879;WS=30857	GT:GL:GOF:GQ:NR:NV	./.:-5.1,-0.0,0.0:20:3:1:1'
p_line2 = '1	30872	.	C	N	5	REFCALL	BRF=.;END=31871;FR=.;FS=.;HP=.;HapScore=.;MGOF=.;MMLQ=.;MQ=.;NF=.;NR=.;PP=.;QD=.;ReadPosRankSum=.;SC=.;START=.;SbPval=.;Size=1000;Source=.;TC=.;TCF=.;TCR=.;TR=.;WE=.;WS=.	GT:GL:GOF:GQ:NR:NV	./.:-1,-1,-1:-1:-1:2:0'
'''

def is_header(line):
    """Check if line is header."""
    
    return line.startswith('#')

def is_var(line):
    """Check if line indicates a variant or a hom-ref block. A hom-ref 
    block is characterized by the 'END' tag under INFO."""
    
    return line.find('END') == -1
    
# FIELD index
# CHROM 0, POS 1, REF 3, ALT 4, QUAL 5, FILTER 6, INFO 7, FORMAT 8, sample 9

def get_END(line):
    """Extract the END position of a hom-ref line."""
    
    INFO = line.strip().split('\t')[7]
    INFO_fields = INFO.split(';')
    for INFO_field in INFO_fields:
        if INFO_field.split('=')[0] == 'END':
            END = int(INFO_field.split('=')[1])
            return END

def get_GQ(line):
    """Extract the GQ (Genotype Quality) of a line."""
    
    fields = line.strip().split('\t')
    FORMAT, sample = fields[8], fields[9]
    FORMAT_fields = FORMAT.split(':'); sample_fields = sample.split(':')
    
    GQ_index = FORMAT_fields.index('GQ')
    GQ = int(sample_fields[GQ_index])
    return GQ
    
def pass_req(line, min_GQ, min_QUAL, FILTER_phrases):
    """Check if a line passes quality requirements."""
    
    fields = line.strip().split('\t')
    GQ = get_GQ(line)
    QUAL = float(fields[5])
    FILTER = fields[6]
    
    return (min_GQ == None or GQ >= min_GQ) and (min_QUAL == None or QUAL >= min_QUAL) \
                and (FILTER_phrases == None or FILTER in FILTER_phrases)
    
def gvcf_regions(gvcf, bed, min_GQ, min_QUAL, FILTER_phrases):
    """Generate the called regions of a gvcf in a bed file."""
    
    g = open(gvcf); b = open(bed, 'w')
    
    region_CHROM = ''
    for line in g:
        # filter header lines and lines that fail requirements
        if not is_header(line) and pass_req(line, min_GQ, min_QUAL, FILTER_phrases):
            # FIELD index, CHROM 0, POS 1, REF 3
            fields = line.strip().split('\t')
            CHROM = fields[0]
            POS = int(fields[1])
            REF = fields[3]
            
            # vcf: 1 index, inclusive start, inclusive end
            # bed: 0 index, inclusive start, exlusive end
            
            # starting postion of the line in bed index
            line_start = POS - 1            
            # ending position of the line in bed index
            if is_var(line):
                line_end = POS + len(REF)
            else:
                line_end = get_END(line)
            
            if region_CHROM != CHROM:
                # when chromosome changes non-trivally
                if region_CHROM != '':
                    region_end = previous_line_end
                    b.write('\t'.join([region_CHROM, 
                        str(region_start), str(region_end)]) + '\n')
                region_CHROM = CHROM
                region_start = line_start
            # when called regions have a gap between the previous line
            # and the current
            elif line_start > previous_line_end:
                region_end = previous_line_end
                b.write('\t'.join([region_CHROM, 
                        str(region_start), str(region_end)]) + '\n')
                region_start = line_start
                
            previous_line_end = line_end
                
    g.close(); b.close()
    
import argparse
parser = argparse.ArgumentParser(description='Get the called regions of a gvcf file.')
parser.add_argument('gvcf', metavar='GVCF', help='input gvcf file')
parser.add_argument('bed', metavar='BED', help='output bed file')

parser.add_argument('--min_GQ', type=int, help='minimum GQ (Genotype Quality) considered as pass')
parser.add_argument('--min_QUAL', type=float, help='minimum QUAL considered as pass')
parser.add_argument('--FILTER_phrases',  nargs='+', help='list of FILTER phrases considered as pass, e.g., PASS, REFCALL')

# presets of min_GQ, min_QUAL, FILTER_phrases
freebayes_preset = [5, None, ['PASS']]
gatk_preset = [5, None, None]
platypus_preset = [None, 5, ['PASS', 'REFCALL']]

gvcf_type_help = '''type of gvcf output.
        [min_GQ, min_QUAL, FILTER_phrases] presets:
        freebayes: %s.
        gatk: %s.
        platypus: %s.''' % (freebayes_preset, gatk_preset, platypus_preset)

parser.add_argument("--gvcf_type", choices=['freebayes', 'gatk', 'platypus'],
                                help=gvcf_type_help)
args = parser.parse_args()

pass_par = [args.min_GQ, args.min_QUAL, args.FILTER_phrases]

if args.gvcf_type == 'freebayes':
    pass_par = freebayes_preset
elif args.gvcf_type == 'gatk':
    pass_par = gatk_preset
elif args.gvcf_type == 'platypus':
    pass_par = platypus_preset
    
par = [args.gvcf, args.bed] + pass_par
    
print(par)
# gvcf_regions(*par)

