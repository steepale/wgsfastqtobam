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

(manifest_ch1, manifest_ch2, manifest_ch3, manifest_ch4) = manifest_ch.into(4)

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
    //publishDir "${params.workdir}/test", mode: 'copy'
    if (params.echo) {
        echo true
    }

    input:
    set pair_id, file(read1), file(read2) from read_pairs_trimmomatic
    // set val(sample_id), val(barcode), val(lane), val(suffix), val(sample_label), val(type) from manifest_ch

    output:
    file "*.fastq.gz" into trimmomatic_ch_out

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
    ${read1.simpleName}_paired.fastq.gz \
    ${read1.simpleName}_unpaired.fastq.gz \
    ${read2.simpleName}_paired.fastq.gz \
    ${read2.simpleName}_unpaired.fastq.gz \
    ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 \
    HEADCROP:9
    """
}

/*
 * STEP: Trim reads with Sickle
 */

// Seperate the output channel from trimmomatic into sickle
(trim_ch1, trim_ch2, trim_ch3, trim_ch4) = trimmomatic_ch_out.into(4)

// Perform read trimming with sickle
process sickle {
    cache true
    container "steepale/sickle:1.0"
    publishDir "${params.workdir}/test", mode: 'copy'
    if (params.echo) {
        echo true
    }

    input:
    file(read1_paired) from trim_ch1.flatten().filter( ~/.*R1.*_paired.fastq.gz/ )
    file(read2_paired) from trim_ch2.flatten().filter( ~/.*R2.*_paired.fastq.gz/ )
    file(read1_unpaired) from trim_ch3.flatten().filter( ~/.*R1.*_unpaired.fastq.gz/ )
    file(read2_unpaired) from trim_ch4.flatten().filter( ~/.*R2.*_unpaired.fastq.gz/ )

    output:
    file "*_sickle.fastq" into sickle_out_ch

    // conditional to skip prcoess if needed
    when:
    !params.skip_sickle

    script:
    """
    ### Perform fastqc on all fastq files
    echo "read1_paired: ${read1_paired}"
    echo "read2_paired: ${read2_paired}"
    echo "read1_unpaired: ${read1_unpaired}"
    echo "read2_unpaired: ${read2_unpaired}"
    # Trim the paired reads
    sickle pe -f ${read1_paired} \
    -r ${read2_paired} \
    -t sanger \
    -o ${read1_paired.simpleName}_sickle.fastq \
    -p ${read2_paired.simpleName}_sickle.fastq \
    -s singles_PE_sickle.fastq \
    -q 20 -l 50 -g
    """
}

/*
 * STEP: Perform mapping with bwa
 */

// Split the sickle output
(sic_ch1, sic_ch2) = sickle_out_ch.into(2)
ref_bwa_ch1 = Channel.from(file(params.genome))
ref_bwa_ch2 = Channel.fromPath("${params.genome_dir}/*{amb,ann,bwt,pac,sa}").collect()
//(ref_bwa_ch1, ref_bwa_ch2) = ref_bwa_ch.into(2)

// Perform read alignment with bwa
process bwa {
    cache true
    container "steepale/bwa:1.0"
    publishDir "${params.workdir}/test", mode: 'copy'
    if (params.echo) {
        echo true
    }

    input:
    file(read1_paired) from sic_ch1.flatten().filter( ~/.*R1.*_paired_sickle.fastq/ )
    file(read2_paired) from sic_ch2.flatten().filter( ~/.*R2.*_paired_sickle.fastq/ )
    //file genome from ref_bwa_ch1.flatten().filter( ~/galgal5.fa/ )
    file genome from ref_bwa_ch1
    file genome_all from ref_bwa_ch2.collect()
    set val(sample_id), val(barcode), val(lane), val(suffix), val(sample_label), val(type) from manifest_ch1

    output:
    file "*.sam" into bwa_out_ch

    // conditional to skip prcoess if needed
    when:
    !params.skip_sickle

    script:
    """
    echo ${genome}
    echo ${read1_paired}
    echo ${read2_paired}
    echo ${sample_id}
    echo ${genome_all}
    bwa mem \
    -t 2 \
    -T 20 \
    ${genome} \
    ${read1_paired} \
    ${read2_paired} \
    > ${sample_id}.sam
    """
}

/*
 * STEP: Add ReadGroups with Picard
 */

// Add ReadGroups with Picard
process readgroups {
    cache true
    container "steepale/gatk:3.5"
    publishDir "${params.workdir}/test", mode: 'copy'
    if (params.echo) {
        echo true
    }

    input:
    val rg_mem from params.rg_mem
    file bwa_sam from bwa_out_ch
    set val(sample_id), val(barcode), val(lane), val(suffix), val(sample_label), val(type) from manifest_ch2
    val platform from params.platform


    output:
    file "*_rg.sam" into rg_out_ch

    // conditional to skip prcoess if needed
    when:
    !params.skip_readgroups

    script:
    """
    echo ${rg_mem}
    echo ${bwa_sam}
    echo ${bwa_sam.simpleName}
    echo ${lane}
    echo ${sample_id}
    echo ${type}

    # Add readgroups with picard
    java -Xmx${rg_mem} \
    -jar /opt/picard/build/libs/picard.jar \
    AddOrReplaceReadGroups \
    INPUT=${bwa_sam} \
    OUTPUT=${bwa_sam.simpleName}_${lane}_rg.sam \
    RGID=${sample_id}_${lane} \
    RGPL= ${platform} \
    RGPU=${barcode}_${lane} \
    RGSM=${sample_id} \
    RGLB=${type}
    """
}

/*
 * STEP: Sort, Merge, and Realign Bam Files (by lane)
 */

// Input Channels
ref_smr_by_lane_ch = Channel.from(file(params.genome))
ref_smr_by_lane_ch2 = Channel.fromPath("${params.genome_dir}/*{amb,ann,bwt,pac,sa,fai,dict}").collect()

// Sort, Merge, and Realign Bam Files (by lane)
process smr_by_lane {
    cache true
    container "steepale/gatk:3.5"
    publishDir "${params.workdir}/test", mode: 'copy'
    if (params.echo) {
        echo true
    }

    input:
    val smr_by_lane_mem from params.smr_by_lane_mem
    file rg_sam from rg_out_ch
    set val(sample_id), val(barcode), val(lane), val(suffix), val(sample_label), val(type) from manifest_ch3
    file genome from ref_smr_by_lane_ch
    file genome_ann from ref_smr_by_lane_ch2

    output:
    file "*.{bam,bai,list,txt}" into rsmr_by_lane_out_ch

    // conditional to skip prcoess if needed
    when:
    !params.skip_smr_by_lane

    script:
    """
    # Sort the sam file
    java -Xmx${smr_by_lane_mem} \
    -jar /opt/picard/build/libs/picard.jar \
    SortSam \
    I=${rg_sam} \
    O=${rg_sam.baseName}.bam \
    SORT_ORDER=coordinate

    echo "Bam file sorted"

    # Mark duplicates in the bam file
	java -Xmx${smr_by_lane_mem} \
    -jar /opt/picard/build/libs/picard.jar \
    MarkDuplicates \
    I=${rg_sam.baseName}.bam \
    O=${rg_sam.baseName}_marked.bam \
    METRICS_FILE=${rg_sam.baseName}_metrics.txt \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

    echo "Bam file duplicates marked"

    # Build an index of bam file
	java -Xmx${smr_by_lane_mem} \
    -jar /opt/picard/build/libs/picard.jar \
    BuildBamIndex \
    I=${rg_sam.baseName}_marked.bam \
    O=${rg_sam.baseName}_marked.bai

    echo "Bam indexed"

    # Create targets for indel realignment
	java -Xmx${smr_by_lane_mem} \
    -jar /opt/gatk3.5/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R ${genome} \
    -I ${rg_sam.baseName}_marked.bam \
    -o ${rg_sam.baseName}_intervals.list

    echo "Targets created for indel relaignment"

    # Realign around indels
	java -Xmx${smr_by_lane_mem} \
    -jar /opt/gatk3.5/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R ${genome} \
    -I ${rg_sam.baseName}_marked.bam \
    -targetIntervals ${rg_sam.baseName}_intervals.list \
    -o ${rg_sam.baseName}_realigned.bam

    echo "Bam realigned around indels"

    """
}




/*

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

