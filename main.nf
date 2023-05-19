#!/usr/bin/env nextflow

// Use DSL2
nextflow.preview.dsl=2

// QUEST nextflow version message
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

// Variables
date = new Date().format('yyyyMMdd')

// Parameters
params.help = null
params.data = null
params.boot_n = null
params.seed = 99
params.out = "${workflow.projectDir}/BOOT_results_${date}"
params.bin_dir = "${workflow.projectDir}/bin"

// LOG AND HELP MESSAGE SETUP
if (!params.help) {
log.info '''
B O O T    I T   U P ! ! !
===============================================
'''
    log.info ""
    log.info "data                     = ${params.data}"
    log.info "boot_n                   = ${params.boot_n}"
    log.info "seed                     = ${params.seed}"
    log.info "output directory         = ${params.out}"
    log.info "bin dir                  = ${params.bin_dir}"
    log.info ""
    } else {
log.info '''
B O O T    I T   U P ! ! !
===============================================
'''
    log.info "Usage:"
    log.info "The typical command for running the pipeline is as follows:"
    log.info "nextflow run main.nf --data /projects/b1059/projects/Tim/DFE_boots-nf/input_data/06_DFE_joined_imputed_geno_10kbMax_pheno.csv --boot_n 100 --seed 99"
    log.info ""
    log.info "Required Arguments:"
    log.info "--data        String       Full path to the fitness data for RI(A)ILs"
    log.info "--boot_n      Float        The number of bootstraps to run"
    log.info ""
    log.info "Optional Arguments:"
    log.info "--seed        Float        A seed for reproducing analysis, default is 99"
    log.info "--out         String       The output directory, default is BOOT_results_<date>"
    log.info ""
    log.info "Flags:"
    log.info "--help                     Display this message"
    log.info ""
    log.info "--------------------------------------------------------"
        exit 1
    }

/*
~ ~ ~ > * WORKFLOW
*/

workflow {
    // Get inputs
    config = Channel.fromPath("${params.bin_dir}/config.boots.R")
        .combine(Channel.from("${params.data}"))
        .combine(Channel.from("${params.boot_n}"))
        .combine(Channel.from("${params.seed}"))
        .combine(Channel.from("${params.out}"))
        //.view()

    // configure boots
    configBoot(config)

    // get the configuration for boots
    boots = configBoot.out.boots_file
        .splitCsv(header:true, sep: "\t")
        .map { row ->
                [file("${params.out}/data/configBoots.rda"), row.boots, "${params.out}", file("${params.bin_dir}/runBoots.R")]
        }
    //.view()

    // run the bootstraps
    runBoots(boots)

    // send output to procBoots channel
    procBoots_ch = runBoots.out.bootsOutput.collect() | procBoots
    
}

process configBoot {
    publishDir "${params.out}/data", mode: 'copy', pattern: "*.rda"
    conda '/projects/b1059/software/conda_envs/R4.2.2'

    input:
        tuple file(config_script), val(data), val(boot_n), val(seed), val(out)

    output:
        path "configBoots.rda", emit: boots_setup_rda
        path "boots.tsv", emit: boots_file
        

    """
        # Configure boots with configBoots.R
        Rscript --vanilla ${config_script} ${data} ${boot_n} ${seed} ${out}

    """
}

process runBoots {
    //publishDir "${params.out}/data", mode: 'copy', pattern: "*.tsv"
    conda '/projects/b1059/software/conda_envs/R4.2.2'

    input:
        tuple file(data), val(boot), val(out), file(script)

    output:
        path "*.tsv", emit: bootsOutput

    """
        # Configure boots with configBoots.R
        Rscript --vanilla ${script} ${data} ${boot} ${out}

    """
}

process procBoots {
    publishDir "${params.out}/data", mode: 'copy', pattern: "*.tsv"
    
    input:
        path("*")

    output:
        path "*.tsv", emit: mergedBoots

    """
        # merge all the boots
        awk 'FNR==1 && NR!=1 {next} 1' *.tsv >> mergedBoots.tsv

    """
}
