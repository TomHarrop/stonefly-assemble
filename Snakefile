#!/usr/bin/env python3

import multiprocessing
from pathlib import Path

def resolve_path(x):
    return(Path(x).resolve().as_posix())

bbduk = 'shub://TomHarrop/seq-utils:bbmap_38.76'
busco = 'shub://TomHarrop/assembly-utils:busco_4.0.4'
flye = 'shub://TomHarrop/assemblers:flye_2.6-g47548b8'
funannotate = 'shub://TomHarrop/funannotate-singularity:funannotate_1.7.4'
te_tools = 'shub://TomHarrop/funannotate-singularity:tetools_1.1'

raw_pb = 'data/pbraw/P01DY19168939-1_r64053_20191111_075118_1_A01.subreads.bam'
cpus = min(multiprocessing.cpu_count(), 128)


rule target:
    input:
        expand(('output/030_busco/'
                '{assembly}/run_insecta_odb10/'
                'full_table.tsv'),
               assembly=['assembly', 'scaffolds']),
        expand('output/040_stats/stats.{assembly}.tsv',
               assembly=['assembly', 'scaffolds']),
        'output/060_funannotate/predict_results/ASW.mrna-transcripts.fa'


# run funannotate
rule funannotate_predict:
    input:
        'output/060_funannotate/training/funannotate_train.transcripts.gff3',
        fasta = 'output/095_repeatmasker/assembly.fa.masked',
        db = 'data/fundb_20200227',
        trinity = 'data/Trinity.fasta'
    output:
        'output/060_funannotate/predict_results/ASW.gff3',
        'output/060_funannotate/predict_results/ASW.mrna-transcripts.fa'
    params:
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        db = lambda wildcards, input: resolve_path(input.db),
        wd = resolve_path('output/060_funannotate')
    log:
        'output/logs/funannotate_predict.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        funannotate
    shell:
        'cp /genemark/gm_key_64 ${{HOME}}/.gm_key ; '
        'funannotate predict '
        '-i {params.fasta} '
        '-s stonefly '
        '--transcript_evidence {input.trinity} '
        '-o {params.wd} '
        '-d {params.db} '
        '--cpus {threads} '
        '--augustus_species stonefly '
        '--optimize_augustus '
        '--busco_seed_species fly '
        '--busco_db insecta '
        '--organism other '
        '--repeats2evm '
        '--max_intronlen 10000 '
        '&> {log}'

rule funannotate_train:
    input:
        fasta = 'output/095_repeatmasker/assembly.fa.masked',
        trinity = 'data/Trinity.fasta'
    output:
        'output/060_funannotate/training/funannotate_train.transcripts.gff3',
    params:
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        wd = resolve_path('output/060_funannotate'),
    log:
        'output/logs/funannotate_train.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        funannotate
    shell:
        'cp /genemark/gm_key_64 ${{HOME}}/.gm_key ; '
        'funannotate train '
        '--input {params.fasta} '
        '--out {params.wd} '
        '--trinity {input.trinity} '
        '--no_trimmomatic ' # disabling trimmomatic seems to disable normalising, bleuch
        '--max_intronlen 10000 '
        '--species stonefly '
        '--cpus {threads} '
        '&> {log}'


# repeat modeller / masker
rule rm_mask:
    input:
        cons = 'output/095_repeatmasker/consensi.fa',
        fasta = 'output/095_repeatmasker/assembly.fa'
    output:
        'output/095_repeatmasker/assembly.fa.masked'
    params:
        wd = resolve_path('output/095_repeatmasker'),
        lib = lambda wildcards, input: resolve_path(input.cons),
        fasta = lambda wildcards, input: resolve_path(input.fasta)
    log:
        resolve_path('output/logs/rm_mask.log')
    singularity:
        te_tools
    shell:
        'cd {params.wd} || exit 1 ; '
        'RepeatMasker '
        '-engine ncbi '
        '-pa {cpus} '
        '-lib {params.lib} '
        '-dir {params.wd} '
        '-gccalc -xsmall -gff -html '
        '{params.fasta} '
        '&> {log}'

rule rm_model:
    input:
        'output/095_repeatmasker/assembly.translation'
    output:
        'output/095_repeatmasker/families.stk',
        'output/095_repeatmasker/consensi.fa'
    params:
        wd = resolve_path('output/095_repeatmasker'),
    log:
        resolve_path('output/logs/rm_model.log')
    singularity:
        te_tools
    shell:
        'cd {params.wd} || exit 1 ; '
        'RepeatModeler '
        '-database assembly '
        '-engine ncbi '
        '-pa {cpus} '
        '-dir {params.wd} '
        # '-recoverDir {params.wd} '
        '&> {log}'

rule rm_build:
    input:
        fasta = 'output/020_flye/assembly.fasta',
    output:
        fa = 'output/095_repeatmasker/assembly.fa',
        tx = 'output/095_repeatmasker/assembly.translation'
    params:
        wd = resolve_path('output/095_repeatmasker')
    log:
        resolve_path('output/logs/rm_build.log')
    singularity:
        te_tools
    shell:
        'cp {input.fasta} {output.fa} ; '
        'cd {params.wd} || exit 1 ; '
        'BuildDatabase '
        '-name assembly '
        '-engine ncbi '
        '-dir {params.wd} '
        '&> {log} '


rule assembly_stats:
    input:
        fasta = 'output/020_flye/{assembly}.fasta',
    output:
        stats = 'output/040_stats/stats.{assembly}.tsv'
    log:
        'output/logs/assembly_stats.{assembly}log'
    threads:
        1
    singularity:
        bbduk
    shell:
        'stats.sh '
        'in={input} '
        'minscaf=1000 '
        'format=3 '
        'threads={threads} '
        '> {output} '
        '2> {log}'


rule busco:
    input:
        fasta = 'output/020_flye/{assembly}.fasta',
        lineage = 'data/insecta_odb10'
    output:
        ('output/030_busco/'
         '{assembly}/run_insecta_odb10/'
         'full_table.tsv'),
    log:
        Path(('output/logs/'
              'busco.{assembly}.log')).resolve()
    params:
        wd = 'output/030_busco',
        name = '{assembly}',
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
        '--mode genome '
        '&> {log}'

rule extract_scaffolds:
    input:
        'output/020_flye/assembly.fasta'
    output:
        'output/020_flye/scaffolds.fasta'
    singularity:
        bbduk
    log:
        'output/logs/extact_scaffolds.log'
    shell:
        'filterbyname.sh '
        'in={input} '
        'names=scaffold '
        'include=t '
        'prefix=t '
        'out={output} '
        '2> {log}'


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
