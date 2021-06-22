#!/usr/bin/env nextflow
/* 
=================================================================================================
                              LONGMETARG
=================================================================================================
Nextflow pipeline for the identifiation of antibiotic resistance genes in metagenomic long reads
-------------------------------------------------------------------------------------------------
@ Author
Louisa Wessels Perelo <louperelo@gmail.com>
#### Homepage / Documentation
 https://github.com/louperelo/longmetarg
-------------------------------------------------------------------------------------------------
*/
def helpMessage() {
    
    //log.info nfcoreHeader()
    log.info"""
    ==============================================================
                              LONGMETARG
    ==============================================================
    Nextflow pipeline for the identifiation of antibiotic resistance genes in metagenomic long reads
    --------------------------------------------------------------
   
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --input metagenome.fasta -profile conda

    Mandatory arguments:
      --input [file]                  Path to input fasta file
      -profile [str]                  Configuration profile to use. 
                                      Available: conda. Soon available: docker, singularity

    Options:
      --readtype [str]               long read type (ONT or PB) vs reference mapping in minimap2
                                     Available: map-ont (default), map-pb
      --flye [bool]                  run assembly with metaFlye before mapping
                                     default: false
      --flyeDt [str]                 specify sequencing type in metaFlye
                                     Available: nano-raw (default), nano-corr, subassemblies,
                                     pacbio-raw, pacbio-corr, pacbio-hifi
      --length [int]                 cutoff ARG alignment length (%), default: 95
      --ident [int]                  cutoff ARG cover identity (%), default: 80
      --genomeSize [int]             estimated genome size for metaFlye assembly, default: 1000000
      --taxon [str]                  taxonomy alignment tool
                                     Available: Diamond (default), Blast
      --pfThreshold [float]          manually specified threshold for probability filtering in PlasFlow
                                     default: 0.7
      --pfModels [file]              custom location of models used for prediction in PlasFlow
                                     default: $PWD/env/plasflow

    References for use with -profile conda                      
      --argDb [file]                 Path to the CARD database fasta
                                     default: $PWD/card_db/card_database_v*.fasta
      --drugClass [file]             Path to CARD aro_index file
                                     default: $PWD/aro_index.tsv
      --blastDbDir [file]            Path to folder containing Blast database
                                     default: $PWD/blast_db
      --blastDbName [str]            Name of Blast database
                                     default: ref_prok_rep_genomes
      --diamondDbDir [file]          Path to folder containing Diamond database
                                     default: $PWD/diamond_db
      --diamondDbName [str]          Name of Diamond database
                                     default: nr.dmnd

    //Other options:
      --outdir [file]                The output directory where the results will be saved
      --threads [int]                Number of parallel threads, default: 15
   
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


// Setup

/*move params to nextflow.config file
*/
/*
params.help = false
params.input = "$PWD/sample.fa"
params.readtype = 'map-ont'
params.argDb = "$PWD/card_db/card_database_v*.fasta"
params.drugClass = "$PWD/aro_index.tsv"
params.blastDbDir =  "$PWD/blast_db"
params.blastDbName = "ref_prok_rep_genomes"
params.diamondDbDir = "$PWD/diamond_db"
params.diamondDbName = "proteinDb.dmnd"
params.flye = false
params.flyeDt = "nano-raw"
params.length = 95
params.ident = 80
params.genomeSize = 1000000
params.outdir = "$PWD/out_dir"
params.taxon = "Diamond" 
params.threads = 15
params.pfThreshold = 0.7
params.pfModels = "$PWD/env/plasflow/models"
*/

in_ch = Channel.fromPath(params.input)
               .into{in_qc; in_flye; in_align; in_analysis} 


/*
 input FASTA statistical summary and plots using NanoPlot
 */
process stats_NanoPlot {
    label 'mapcard' //use to define env in nextflow.config

    input:  
    file in_qc

    publishDir "$params.outdir/nanoPlot"

    output:  
    file '*' //into qcheck                

    script:
    """
    NanoPlot --fasta $in_qc -t $params.threads

    """
}

/*
MetaFlye assembly
execute only if params.flye is set to true
*/
process assembly_MetaFlye {
    label 'metaflye' 

    input:
    file in_flye

    publishDir "$params.outdir"

    output:     
    file 'out_flye/assembly.fasta' into (out_flye, flye_analysis)
    
    when:
    params.flye == true

    script:
    /*
    minimum overlap was set to default 2000 (2kb), also used in Flye paper for Cow rumen set. 
    genome size default was set to 1 G, as in Flye paper
    */
    """
    flye --$params.flyeDt $in_flye --threads $params.threads --out-dir out_flye --meta --plasmids --genome-size $params.genomeSize 

    """
}

/*
Minimap2 alignment to CARD database
*/
process map_CARD {
    label 'mapcard' 

    input:
	file fa from in_align
	
    publishDir "$params.outdir", mode: 'copy'//, overwrite: true

    output:  
    file '*' into card_aln 

	when:
	params.flye	== false
        
    script:
	
    """
    minimap2 -ax $params.readtype $params.argDb $fa -t $params.threads -c | samtools view -S -b | samtools view -b -F 4 > aln_F4.bam

    """ 	
}

process mapFlye_CARD {
    label 'mapcard' 

    input:
    file out_flye
	
    publishDir "$params.outdir", mode: 'copy'

    output:  
    file '*' into flyecard_aln 

	when:
	params.flye	== true
        
    script:
	
    """
    minimap2 -ax $params.readtype $params.argDb $out_flye -t $params.threads | samtools view -S -b | samtools view -b -F 4 > aln_F4.bam

    """ 	
} 

/*
Extract information from minimap2 alignment output
*/
process aln_analysis {
    label 'mapcard' 
  
    input:
    file card_aln
    file in_analysis

    publishDir "$params.outdir", mode:'copy'

    output:   
    file "argHitsDf.tsv" into in_summary
    file "argHits.fasta" into(in_tax, in_pf)
	
	when:
	params.flye == false
           
    script:
    """
    read_analysis.py --bam_file $card_aln --reads_fasta $in_analysis --ident $params.ident --cov_length $params.length --folder $params.outdir --aro_idx $PWD/bin/aro_index.tsv
    
    """
}

process alnFlye_analysis {
    label 'mapcard' 
  
    input:
    file flyecard_aln
    file flye_analysis

    publishDir "$params.outdir", mode:'copy'

    output:   
    file "argHitsDf.tsv" into inflye_summary
    file "argHits.fasta" into(inflye_tax, inflye_pf)
	
	when:
	params.flye == true
             
    script:
    """
    read_analysis.py --bam_file $flyecard_aln --reads_fasta $flye_analysis --ident $params.ident --cov_length $params.length --folder $params.outdir --aro_idx $PWD/bin/aro_index.tsv
    
    """
}

/*
PlasFlow analysis of results.fasta from alignment_analysis, adding column with Plasmid/Chromosome information to results_df
*/ 
process plasflow {
    label 'plasflow'
    
    input:
	file argHits from in_pf.mix(inflye_pf)

    publishDir "$params.outdir", mode:'copy'

    output:   
    file 'result_pf.tsv' into out_pf  
             
    script:
    """
    PlasFlow.py --input $argHits --output result_pf.tsv --threshold $params.pfThreshold --models $params.pfModels
    
    """
} 

/*
BLAST or DIAMOND taxonomic analysis for result.fasta reads
*/
process taxonomy {
    if( params.taxon == "Blast")
        label 'blast' 
    if( params.taxon == "Diamond")
        label 'diamond'
 
    input:
	file argHits from in_tax.mix(inflye_tax)

    publishDir "$params.outdir", mode:'copy'

    output:   
    file 'result_tax.tsv' into (in_ssciname, out_taxonomy)
        
    script:
    if( params.taxon == "Blast")
    """
    blastn -task blastn -max_target_seqs 1 -num_threads $params.threads -outfmt "6 delim=    qseqid staxid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -db $params.blastDbDir/$params.blastDbName -query $argHits -out result_tax.tsv 
    
    """
    else if( params.taxon == "Diamond")
    """
    diamond blastx -d $params.diamondDbDir/$params.diamondDbName -q $argHits -p $params.threads -F 15 -f 6 qseqid staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore --max-target-seqs 1 -o result_tax.tsv
    
    """
    else
    throw new IllegalArgumentException("Unknown taxonomy aligner $params.taxon")   
} 

/*
Get scientific name from NCBI taxid in process taxonomy output
*/ 
process get_ssciname {
    label 'mapcard'
  
    input:
    file in_ssciname

    publishDir "$params.outdir", mode:'copy'

    output:   
    file "*" into out_ssciname
             
    script:
    """
    get_taxon.py --blast_out $in_ssciname

    """
}

/*
join results to output summary
*/ 
process summary {
    label 'mapcard'
 
    input: 
	file argHitsDf from in_summary.mix(inflye_summary)
    file out_pf
    file out_taxonomy
    file out_ssciname

    publishDir "$params.outdir", mode:'copy'

    output:   
    file 'summary.tsv' //into summary
	
    script:
    """
    #!/usr/bin/python
import pandas as pd
aln_df = pd.read_csv("$argHitsDf", sep="\\t")
pf_df = pd.read_csv("$out_pf", sep="\\t")
taxon_df = pd.read_csv("$out_taxonomy", sep="\\t", header=None, names=['contig_name', 'taxid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
ssciname_df = pd.read_csv("$out_ssciname", sep="\\t")
tax_df = pd.merge(taxon_df, ssciname_df[['taxid', 'species', 'class']], on='taxid', how='left')

aln_pf_df = pd.merge(aln_df, pf_df[['contig_name', 'label']], on='contig_name', how='left')
aln_pf_df.rename(columns={'label':'gene_location'}, inplace = True)
aln_pf_df['gene_location'] = aln_pf_df['gene_location'].apply(lambda x: x.split('.')[0])
aln_pf_tax_df = pd.merge(aln_pf_df, tax_df[['contig_name', 'species', 'class', 'pident', 'evalue', 'bitscore']], on='contig_name', how='left')
	
summary_df = aln_pf_tax_df.drop(['index', 'contig_length', 'aln_length', 'aln_nm', 'ref_start', 'ref_length', 'ref_ARO', 'ref_ID', 'ref_NCBI', 'c_tuples', 'matches', 'insertions', 'deletions'], axis=1)
	
summary_df.to_csv('summary.tsv', sep="\\t")
	
    """
} 

