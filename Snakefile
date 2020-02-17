#!/usr/bin/env python3

bbduk = 'shub://TomHarrop/seq-utils:bbmap_38.76'


raw_pb = 'data/pbraw/P01DY19168939-1_r64053_20191111_075118_1_A01.subreads.bam'

rule target:
    input:
        'output/010_reads/pb_raw.fasta.gz'

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
