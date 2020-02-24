#!/usr/bin/env python3

import multiprocessing

busco = 'shub://TomHarrop/assembly-utils:busco_4.0.4'
bbduk = 'shub://TomHarrop/seq-utils:bbmap_38.76'
flye = 'shub://TomHarrop/assemblers:flye_2.6-g47548b8'

raw_pb = 'data/pbraw/P01DY19168939-1_r64053_20191111_075118_1_A01.subreads.bam'

rule target:
    input:
        ('output/030_busco/'
         'busco/run_arthropoda_odb10/'
         'full_table.tsv'),

rule busco:
    input:
        fasta = 'output/020_flye/assembly.fasta',
        lineage = 'data/arthropoda_odb10'
    output:
        ('output/030_busco/'
         'busco/run_arthropoda_odb10/'
         'full_table.tsv'),
    log:
        Path(('output/logs/'
              'busco.log')).resolve()
    params:
        wd = 'output/020_busco',
        name = 'busco',
        fasta = lambda wildcards, input: Path(input.fasta).resolve(),
        lineage = lambda wildcards, input:
            Path(input.lineage).resolve()
    threads:
        multiprocessing.cpu_count()
    singularity:
        busco
    shell:
        'cd {params.wd} || exit 1 ; '
        'busco '
        '--force '
        '--in {params.fasta} '
        '--out {params.name} '
        '--lineage_dataset {params.lineage} '
        '--cpu {threads} '
        '--augustus_species honeybee1 '
        '--mode genome '
        '&> {log}'

rule flye:
    input:
        fa = 'output/010_reads/pb_raw.fasta.gz'
    output:
        'output/020_flye/assembly.fasta'
    params:
        outdir = 'output/020_flye',
        size = '500m'
    threads:
        min(128, multiprocessing.cpu_count())
    log:
        'output/logs/020_flye.log'
    singularity:
        flye
    shell:
        'flye '
        '--resume '
        '--pacbio-raw {input.fa} '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&>> {log}'


rule convert_to_fasta:
    input:
        raw_pb
    output:
        'output/010_reads/pb_raw.fasta.gz'
    singularity:
        bbduk
    log:
        'output/logs/convert_to_fasta.log'
    shell:
        'reformat.sh '
        '-Xmx10g '
        'in={input} '
        'out={output} '
        'int=f '
        'zl=9 '
        '2> {log}'
