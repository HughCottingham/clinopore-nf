process POLISHING_STATS {
        label "short_job"
        conda "${params.clinopore_env}"
        input:
        tuple val(isolate_id), path(input_fasta)
        tuple val(isolate_id), path(output_fasta)
        val(polishing_step)
        output:
        path("${isolate_id}_${polishing_step}_stats.txt")
        """
        bioawk -c fastx '{ print \$name, length(\$seq) }' < $input_fasta > assembly_input
        bioawk -c fastx '{ print \$name, length(\$seq) }' < $output_fasta > assembly_output
        for f in assembly_*;do base=\$(basename \$f _stats.txt);awk '{print FILENAME"\t" \$0}' \$f> \$f"1";mv \$f"1" \$f;done
        cat assembly_* >> polishing_summary1.txt
        awk '{print \$0, "\t${isolate_id}\t${polishing_step}"}' polishing_summary1.txt > ${isolate_id}_${polishing_step}_stats.txt
        """
}

process POLCA_POLISHING_STATS {
        label "short_job"
        conda "${params.clinopore_env}"
        input:
        tuple val(isolate_id), path(input_fasta)
        tuple val(isolate_id), path(output_fasta)
        val(polishing_step)
        output:
        path("${isolate_id}_${polishing_step}_stats.txt")
        """
        bioawk -c fastx '{ print \$name, length(\$seq) }' < $input_fasta > assembly_input
        bioawk -c fastx '{ print \$name, length(\$seq) }' < $output_fasta > assembly_output
        for f in assembly_*;do base=\$(basename \$f _stats.txt);awk '{print FILENAME"\t" \$0}' \$f> \$f"1";mv \$f"1" \$f;done
        cat assembly_* >> polishing_summary1.txt
        awk '{print \$0, "\t${isolate_id}\t${polishing_step}"}' polishing_summary1.txt > ${isolate_id}_${polishing_step}_stats.txt
        """
}

process COMBINE_POLISHING_STATS {
        label "short_job"
        conda "${params.clinopore_env}"

        input:
        file("*_stats.txt")
        val(polishing_step)

        output:
        path("${polishing_step}_polishing_stats.txt"), emit: polishing_stats

        script:
        """
        cat *_stats.txt >> medaka_polishing_stats_in.txt
        combine_single_polishing.py medaka_polishing_stats_in.txt $polishing_step ${polishing_step}_polishing_stats.txt 
        """
}

process COMBINE_POLCA_POLISHING_STATS {
        label "short_job"
        conda "${params.clinopore_env}"

        input:
        file("*_stats.txt")
        val(polishing_step)

        output:
        path("${polishing_step}_polishing_stats.txt")

        script:
        """
        cat *_stats.txt >> medaka_polishing_stats_in.txt
        combine_single_polishing.py medaka_polishing_stats_in.txt $polishing_step ${polishing_step}_polishing_stats.txt 
        """
}

process COMBINE_ALL_POLISHING_STATS {
        label "short_job"
        conda "${params.clinopore_env}"
        publishDir path:("${params.outdir}/stats"), mode: 'copy'

        input:
        file("*_stats.txt")

        output:
        path("polishing_summary.txt")

        script:
        """
        cat *_stats.txt >> polishing_summary_in.txt
        multi_polishing_cleanup.py polishing_summary_in.txt polishing_summary.txt 
        """
}

process ASSEMBLY_STATS {
        label "short_job"
        conda "${params.clinopore_env}"
        input:
        tuple val(isolate_id), path(fasta_file)

        output:
        path("${isolate_id}_stats.txt")

        script:
        """
        assembly_stats.py -a $fasta_file --id $isolate_id > ${isolate_id}_stats.txt
        """
}

process COMBINE_ASSEMBLY_STATS {
        label "short_job"
        conda "${params.clinopore_env}"
        publishDir path:("${params.outdir}/stats"), mode: 'copy'

        input:
        file("*_stats.txt")

        output:
        file('assembly_stats.txt')

        script:
        """
        awk 'FNR==1 && NR!=1 { while (/^assembly/) getline; } 1 {print}' *_stats.txt > assembly_stats.txt
        """
}