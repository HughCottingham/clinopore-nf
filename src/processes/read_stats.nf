process RAW_ONT_STATS {
        label "short_job"
        conda "${params.clinopore_env}"

        input:
        tuple val(isolate_id), path(fastq_file, stageAs: "raw_reads.fastq.gz")

        output:
        path("${isolate_id}_read_stats.txt")

        script:
        """
        NanoPlot -t ${params.threads} --fastq raw_reads.fastq.gz --plots dot -o nanoplot_out --tsv_stats
        clean_up_nanostats.py nanoplot_out/NanoStats.txt ${isolate_id} ${isolate_id}_read_stats.txt
        """
}

process COMBINE_RAW_ONT_STATS {
        label "short_job"
        conda "${params.clinopore_env}"
        publishDir path:("${params.outdir}/stats"), mode: 'copy'

        input:
        file("*_stats.txt")

        output:
        file('ONT_raw_read_stats.txt')

        script:
        """
        for f in *stats.txt;do head -n 2 \$f|tail -n 1 >> ONT_raw_read_stats.txt;done
        (echo -e "Isolate\tReads_>Q10\tReads_>Q12\tReads_>Q15\tReads_>Q5\tReads_>Q7\thighest_Q_read_(with_length):1\thighest_Q_read_(with_length):2\thighest_Q_read_(with_length):3\thighest_Q_read_(with_length):4\thighest_Q_read_(with_length):5\tlongest_read_(with_Q):1\tlongest_read_(with_Q):2\tlongest_read_(with_Q):3\tlongest_read_(with_Q):4\tlongest_read_(with_Q):5\tmean_qual\tmean_read_length\tmedian_qual\tmedian_read_length\tn50\tnumber_of_bases\tnumber_of_reads\tread_length_stdev" && cat ONT_raw_read_stats.txt) > ONT_raw_read_stats.txt1 && mv ONT_raw_read_stats.txt1 ONT_raw_read_stats.txt
        """
}

process FILTERED_ONT_STATS {
        label "short_job"
        conda "${params.clinopore_env}"

        input:
        tuple val(isolate_id), path(fastq_file, stageAs: "filtered_reads.fastq.gz")

        output:
        path("${isolate_id}_read_stats.txt")

        script:
        """
        NanoPlot -t ${params.threads} --fastq filtered_reads.fastq.gz --plots dot -o nanoplot_out --tsv_stats
        clean_up_nanostats.py nanoplot_out/NanoStats.txt ${isolate_id} ${isolate_id}_read_stats.txt
        """
}

process COMBINE_FILTERED_ONT_STATS {
        label "short_job"
        publishDir path:("${params.outdir}/stats"), mode: 'copy'

        input:
        file("*_stats.txt")

        output:
        file('ONT_filtered_read_stats.txt')

        script:
        """
        for f in *stats.txt;do head -n 2 \$f|tail -n 1 >> ONT_filtered_read_stats.txt;done
        (echo -e "Isolate\tReads_>Q10\tReads_>Q12\tReads_>Q15\tReads_>Q5\tReads_>Q7\thighest_Q_read_(with_length):1\thighest_Q_read_(with_length):2\thighest_Q_read_(with_length):3\thighest_Q_read_(with_length):4\thighest_Q_read_(with_length):5\tlongest_read_(with_Q):1\tlongest_read_(with_Q):2\tlongest_read_(with_Q):3\tlongest_read_(with_Q):4\tlongest_read_(with_Q):5\tmean_qual\tmean_read_length\tmedian_qual\tmedian_read_length\tn50\tnumber_of_bases\tnumber_of_reads\tread_length_stdev" && cat ONT_filtered_read_stats.txt) > ONT_filtered_read_stats.txt1 && mv ONT_filtered_read_stats.txt1 ONT_filtered_read_stats.txt
        """
}

process ILLUMINA_STATS {
        label "short_job"
        conda "${params.clinopore_env}"
        input:
        tuple val(isolate_id), path(reads_fwd), path(reads_rev)

        output:
        path("*.zip")

        script:
        """
        fastqc --noextract $reads_fwd $reads_rev
        """
}

process COMBINE_ILLUMINA_STATS {
        label "short_job"
        conda "${params.clinopore_env}"
        publishDir path: {"${params.outdir}/stats"}, mode: 'copy', saveAs: {filename -> "Illumina_read_stats.txt"}

        input:
        path("*")

        output:
        path("multiqc_data/multiqc_fastqc.txt")

        script:
        """
        multiqc --force .
        """
}
