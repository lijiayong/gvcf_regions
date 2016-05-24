gVCF Regions
========

gVCF Regions is a command line tool that output the called regions of a gVCF file in BED format. 
It handles four main types of gVCFs (Complete Genomics, Freebayes, GATK, Platypus), with the capability 
to customize the settings of 'called regions'.

### gVCF Regions Command
```
gvcf_regions.py [-h] [--unreported_is_called]
                       [--ignore_phrases IGNORE_PHRASES [IGNORE_PHRASES ...]]
                       [--min_GQ MIN_GQ] [--min_QUAL MIN_QUAL]
                       [--pass_phrases PASS_PHRASES [PASS_PHRASES ...]]
                       [--gvcf_type {complete_genomics,freebayes,gatk,platypus}]
                       GVCF

Output the called regions of a gvcf file to stdout in bed format.

positional arguments:
  GVCF                  input gvcf file

optional arguments:
  -h, --help            show this help message and exit
  --unreported_is_called
                        use this flag to treat unreported sites as called
  --ignore_phrases IGNORE_PHRASES [IGNORE_PHRASES ...]
                        list of phrases considered as discarded, e.g., CNV,
                        ME. A line that contains any of the ignore phrases is
                        discarded.
  --min_GQ MIN_GQ       minimum GQ (Genotype Quality) considered as called
  --min_QUAL MIN_QUAL   minimum QUAL considered as called
  --pass_phrases PASS_PHRASES [PASS_PHRASES ...]
                        list of phrases considered as called, e.g., PASS,
                        REFCALL. A line must contain any of the pass phrases
                        to be considered as called.
  --gvcf_type {complete_genomics,freebayes,gatk,platypus}
                        type of gvcf output. [unreported_is_called,
                        ignore_phrases, min_GQ, min_QUAL, pass_phrases]
                        presets: complete_genomics: [True, ['CNV', 'INS:ME'],
                        None, None, ['PASS']]. freebayes: [False, None, None,
                        None, ['PASS']]. gatk: [False, None, 5, None, None].
                        platypus: [False, None, None, None, ['PASS',
                        'REFCALL']].
```
