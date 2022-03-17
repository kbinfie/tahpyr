nextflow.enable.dsl=2


process RADINITIO_MAKE_POPULATION {
    tag "make_population"

    publishDir params.publish_dir

    input:
        path genome
        path chromosome_list
        path output_dir
    output:
        path "${output_dir}/radinitio.log", emit: log
        path "${output_dir}/popmap.tsv", emit: popmap
        path "${output_dir}/msprime_vcfs/*", emit: vcfs
        path "${output_dir}/ref_loci_vars/*", emit: master_vcf
    
    script:
    """
    mkdir -p {output_dir}

    radinitio --make-population \\
              --genome ${genome} \\
              --chromosomes ${chromosome_list} \\
              --out-dir ${output_dir} \\
              --n-pops ${params.n_pops} \\
              --pop-eff-size ${params.pop_eff_size} \\
              --n-seq-indv ${params.n_seq_indiv}
    """

    stub:
    """
    mkdir -p ${output_dir}

    """
}

process RADINITIO_MAKE_LIBRARY {
    tag "make_library"

    publishDir params.publish_dir

    input:
        path genome
        path chromosome_list
        path output_dir
        path popmap
        path vcfs
        path master_vcf
    output:
        path "rad_alleles/*", emit: alleles
        path "rad_reads/*", emit: reads
        path "radinitio.log", emit: log
        path "ref_loci_var/reference_rad_loci.fa.gz", emit: ref_vcf
        path "ref_loci_var/reference_rad_loci.stats.gz", emit: stats
        path "sequence_clone_distrib.tsv", emit: clone_distrib
    
    script:
    def enzymes = params.library_type == "sdRAD" ? "--library-type sdRAD --enz ${params.enzyme_1}" : "--library-type ddRAD --enz ${params.enzyme_1} --enz2 ${params.enzyme_2}"
    """
    # recreate the population folder
    mkdir -p pop_simulation/msprime_vcfs
    mkdir -p pop_simulation/ref_loci_vars
    mv ${vcfs} pop_simulation/msprime_vcfs/
    mv ${master_vcf} pop_simulation/ref_loci_vars/
    mv ${popmap} pop_simulation/

    mkdir -p ${output_dir}
    
    radinitio --make-library \\
              --genome ${genome} \\
              --chromosomes ${chromosome_list} \\
              --out-dir ${output_dir} \\
              --make-pop-sim-dir pop_simulation \\
              ${enzymes} \\
              --insert-mean ${params.insert_mean} \\
              --insert-stdev ${params.insert_sd} \\
              --pcr-cycles ${params.pcr_cycles} \\
              --coverage ${params.coverage} \\
              --read-length ${params.read_length}
    """

    stub:
    """
    mkdir -p ${output_dir}

    """
}

process RADINITIO_TALLY {
    tag "tally"

    publishDir params.publish_dir

    input:
        path genome
        path chromosomes
        val vars
        val out_dir
    output:
        path "${out_dir}/${vars.ID}_radinitio.log", emit: log
        path "${out_dir}/ref_loci_vars/${vars.ID}_reference_rad_loci.fa.gz", emit: ref_fasta
        tuple val(vars),path("${out_dir}/ref_loci_vars/${vars.ID}_reference_rad_loci.stats.gz"), emit: stats

    
    script:
    def enzymes = vars.LIBRARY_TYPE == "sdRAD" ? "--library-type sdRAD --enz ${vars.ENZ1}" : "--library-type ddRAD --enz ${vars.ENZ1} --enz2 ${vars.ENZ2}"
    """
    mkdir -p ${out_dir}

    /usr/local/anaconda3/envs/tahpyr/bin/radinitio --tally-rad-loci \\
              --genome ${genome} \\
              --chromosomes ${chromosomes} \\
              --out-dir ${out_dir} \\
                ${enzymes}

    mv ${out_dir}/radinitio.log ${out_dir}/${vars.ID}_radinitio.log
    mv ${out_dir}/ref_loci_vars/reference_rad_loci.fa.gz ${out_dir}/ref_loci_vars/${vars.ID}_reference_rad_loci.fa.gz
    mv ${out_dir}/ref_loci_vars/reference_rad_loci.stats.gz ${out_dir}/ref_loci_vars/${vars.ID}_reference_rad_loci.stats.gz
    """

    stub:
    """
    mkdir -p ${output_dir}

    """
}

process COUNT_LOCI {
    tag "count_loci"

    input:
        tuple val(vars), path(stats)
    output:
        path "loci_count.csv", emit: loci_count
    """
    count_loci=`gzip -dc ${stats} | grep kept | wc -l`
    echo "ID,LIBRARY_TYPE,ENZ1,ENZ2,MEDIAN_INSERT_SIZE,STDEV_INSERT_SIZE,COVERAGE,PCR_CYCLES,READ_LENGTH,LOCI_COUNT" > loci_count.csv
    echo "${vars.ID},${vars.LIBRARY_TYPE},${vars.ENZ1},${vars.ENZ2},${vars.MEDIAN_INSERT_SIZE},${vars.STDEV_INSERT_SIZE},${vars.COVERAGE},${vars.PCR_CYCLES},${vars.READ_LENGTH},\${count_loci}" >> loci_count.csv
    """
}

// workflows

workflow TALLY_WF {
    take:
        genome
        chromosomes
        vars
        out_dir
    main:
        RADINITIO_TALLY(genome=genome, chromosomes=chromosomes, vars=vars, out_dir=out_dir)
        COUNT_LOCI(stats=RADINITIO_TALLY.out.stats) |
        collectFile(keepHeader: true, name: "summary_stats.csv", storeDir: params.publish_dir, sort: {it.text.split()[1].split(",")[9].toInteger() * -1})
    // emit:
    //     logs = RADINITIO_TALLY.out.log
    //     ref_fasta = RADINITIO_TALLY.out.ref_fasta
    //     stats = RADINITIO_TALLY.out.stats
}


workflow {
    ref_genome = file(params.ref_genome)
    chromosome_list = file(params.chromosome_list)
    Channel.fromPath(params.exp_design).splitCsv(header: true).take(10).set{exp_design}
    TALLY_WF(genome=ref_genome, chromosomes=chromosome_list, vars=exp_design, out_dir=params.out_dir)
}