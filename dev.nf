#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/wgsfastqtobam
========================================================================================
 nf-core/wgsfastqtobam Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/wgsfastqtobam
----------------------------------------------------------------------------------------
// Dev test runs

/*
    Steps in a typical fastq to bam workflow
    1) Perform QC analysis with FASTQC
    
*/

// Help message
def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/wgsfastqtobam \
    --ann genome.gtf \
    --bamlist bamlist.txt \
    --outdir = 'results'
    
    Mandatory arguments:
      --ann                         Path to the annotation file
      --roi                         Path to the region of interest file
      --bamlist                     Path to a tab seperated bamlist file

    Optional arguments:
      --outdir                      The output directory where the results will be saved
      --choi                        Path to chromosomes of interest file
      --chrom_sizes                 Path to file with chromosome sizes
      --roi_script                  Path to the Region Of Interest script to parse regions of interest from the annotation file
      -resume                       Resume the pipeline from prior runs
      -with-docker                  Launch docker containers
      --help                        Dsiplay this help message

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help) {
    // Invoke the function above to display the help message
    helpMessage()
    // Exit out and do not run
    exit 0
}


/*
========================================================================================
                         ANALYSIS STEPS
========================================================================================
 */

/*
 * STEP 1. Read in manifest file
 */

// Create the channels
// Arrays
//config_barcodes_ch = Channel.create()
//config_lanes_ch = Channel.create()
//config_reads_ch = Channel.create()
// Variables
//config_suffix_ch =Channel.create()
//config_sample_label_ch =Channel.create()
//config_sample_type_ch =Channel.create()

// Collect all fields into channel
Channel
    .fromPath("${params.sample_cfg}")
    .splitCsv(skip:1, sep: "\t")
    .map { row -> tuple(
        row[0],row[1],row[2],row[3],row[4],row[5],row[6])
    }
    .set{manifest_ch}

/*
 * STEP: Perform fastqc analysis
 */

// Collect fastq files into channel
read_pairs_fastqc = Channel.fromFilePairs(params.reads, flat: true)

// TODO: read 1 and read 2 are hard-coded, redesign to be more reproducible
process fastqc_1 {
    cache true
    container "steepale/fastqc:1.0"
    publishDir "${params.workdir}/results/fastqc_1", mode: 'copy'
    if (params.echo) {
        echo true
    }

    input:
    set pair_id, file(read1), file(read2) from read_pairs_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastQCreport

    // conditional to skip prcoess if needed
    when:
    !params.skip_fastqc_1

    script:
    """
    ### Perform fastqc on all fastq files
    fastqc \
    -t 2 \
    -q \
    ${read1} \
    ${read2}
    """
}

/*
 * STEP: Trim adapters with Trimmomatic
 */

// Collect fastq files into channel
read_pairs_trimmomatic = Channel.fromFilePairs(params.reads, flat: true)

// Perform adapter trimming with trimmomatic
process trimmomatic {
    cache true
    container "steepale/trimmomatic:1.0"
    publishDir "${params.workdir}/test", mode: 'copy'
    if (params.echo) {
        echo true
    }

    input:
    set pair_id, file(read1), file(read2) from read_pairs_trimmomatic
    // set val(sample_id), val(barcode), val(lane), val(suffix), val(sample_label), val(type) from manifest_ch

    output:
    file "*paired.fastq.gz" into trimmomatic_ch_out

    // conditional to skip prcoess if needed
    when:
    !params.skip_trimmomatic

    script:
    """
    ### Perform fastqc on all fastq files
    java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar \
    PE \
    -threads 2 \
    ${read1} \
    ${read2} \
    ${read1.baseName}_paired.fastq.gz \
    ${read1.baseName}_unpaired.fastq.gz \
    ${read2.baseName}_paired.fastq.gz \
    ${read2.baseName}_unpaired.fastq.gz \
    ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 \
    HEADCROP:9
    """
}

/*
 * STEP: Trim reads with Sickle
 */

// Mix the output channel from trimmomatic into sickle

// Perform adapter trimming with trimmomatic
process sickle {
    cache true
    container "steepale/sickle:1.0"
    publishDir "${params.workdir}/test", mode: 'copy'
    if (params.echo) {
        echo true
    }

    input:
    set pair_id, file(read1), file(read2) from read_pairs_sickle
    // set val(sample_id), val(barcode), val(lane), val(suffix), val(sample_label), val(type) from manifest_ch

    output:
    file "*paired.fastq.gz" into sickle_ch_out

    // conditional to skip prcoess if needed
    when:
    !params.skip_sickle

    script:
    """
    ### Perform fastqc on all fastq files
    java -jar /opt/sickle-0.39/sickle-0.39.jar \
    PE \
    -threads 2 \
    ${read1} \
    ${read2} \
    ${read1.baseName}_paired.fastq.gz \
    ${read1.baseName}_unpaired.fastq.gz \
    ${read2.baseName}_paired.fastq.gz \
    ${read2.baseName}_unpaired.fastq.gz \
    ILLUMINACLIP:/opt/sickle-0.39/adapters/TruSeq2-PE.fa:2:30:10 \
    HEADCROP:9
    """
}






/*
fastqc \
    -t 2 \
    -q \
    ${read1} \
    ${read2}
 */


// A function for a header to the help message
def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/wgsfastqtobam v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

