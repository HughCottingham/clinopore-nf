process FILTER {
        label "short_job"
        conda "${params.clinopore_env}"
        input:
        tuple val(isolate_id), path(reads_se)

        output:
        tuple val(isolate_id), path("${isolate_id}_assembly_input.fastq.gz"), emit: assembly_input

        script:
        """
        filtlong --min_length 1000 --keep_percent 95 $reads_se |gzip > ${isolate_id}_assembly_input.fastq.gz
        """
}

process ASSEMBLE {
        label "long_job"
        conda "${params.clinopore_env}"
        publishDir path:("${params.outdir}/flye"), mode: 'copy', saveAs: {filename -> "${isolate_id}_flye.fasta"}, pattern: '*fasta'
        publishDir path:("${params.outdir}/gfa"), mode: 'copy', saveAs: {filename -> "${isolate_id}.gfa"}, pattern: '*gfa'
        input:
        tuple val(isolate_id), path(input_ch, stageAs: "assembly_input.fastq.gz")

        output:
        tuple val(isolate_id), path("${isolate_id}_flye.fasta"), emit: assembly_fasta
        path("${isolate_id}.gfa")
        
        script:
        """
        flye --nano-hq assembly_input.fastq.gz --out-dir flye_out --threads ${params.threads}
        reordering_contigs.py flye_out/assembly.fasta flye_out/${isolate_id}_inter1.fasta flye_out/${isolate_id}_inter2.fasta flye_out/assembly_graph.gfa ${isolate_id}.gfa flye_out/assembly_info.txt flye_out/${isolate_id}_assembly_info.txt flye_out/${isolate_id}_inter3.fasta ${isolate_id}_flye.fasta
        mv flye_out/assembly.fasta flye_out/${isolate_id}_old.fasta
        mv flye_out/assembly_graph.gfa flye_out/${isolate_id}_old.gfa
        mv flye_out/assembly_info.txt flye_out/${isolate_id}_assembly_info_old.txt
        """
}

process MEDAKA {
        label "medium_job"
        conda "${params.clinopore_env}"
        publishDir path:("${params.outdir}/medaka"), mode: 'copy', saveAs: {filename -> "${isolate_id}_medaka.fasta"}, pattern: '*fasta'
        publishDir path:("${params.outdir}/stats"), mode: 'copy', saveAs: {filename -> "${isolate_id}_medaka_polishing_stats.txt"}, pattern: '*stats.txt'

        input:
        tuple val(isolate_id), path(assembly), path(filtered_reads)

        output:
        tuple val(isolate_id), path("${isolate_id}_medaka.fasta"), emit: assembly_fasta

        script:
        """
        medaka_consensus -d ${assembly} -o . -i ${filtered_reads} -t ${params.threads} -m ${params.medaka_model}
        mv consensus.fasta ${isolate_id}_medaka_inter1.fasta
        contig_renaming.py ${assembly} ${isolate_id}_medaka_inter1.fasta ${isolate_id}_medaka_inter2.fasta ${isolate_id}_medaka_inter3.fasta
        seqkit sort --by-length --reverse ${isolate_id}_medaka_inter3.fasta > ${isolate_id}_medaka.fasta
        """
}

process POLYPOLISH {
        label "short_job"
        conda "${params.clinopore_env}"
        publishDir path:("${params.outdir}/polypolish"), mode: 'copy', saveAs: {filename -> "${isolate_id}_polypolish.fasta"}, pattern: '*polypolish.fasta'        

        input:
        //tuple from medaka polished
        //tuple pe reads
        tuple val(isolate_id), path(assembly), path(illumina1), path(illumina2), path(flye_assembly, stageAs: "flye_input.fasta")
        output:
        tuple val(isolate_id), path("${isolate_id}_polypolish.fasta"), emit: assembly_fasta

        script:
        """
        bwa index ${assembly}
        bwa mem -t ${params.threads} -a ${assembly} ${illumina1} > ${isolate_id}_r1.sam
        bwa mem -t ${params.threads} -a ${assembly} ${illumina2} > ${isolate_id}_r2.sam
        polypolish_insert_filter.py --in1 ${isolate_id}_r1.sam --in2 ${isolate_id}_r2.sam --out1 ${isolate_id}_filtered_r1.sam --out2 ${isolate_id}_filtered_r2.sam
        polypolish ${assembly} ${isolate_id}_filtered_r1.sam ${isolate_id}_filtered_r2.sam| sed 's/_polypolish//' > ${isolate_id}_polypolish1.fasta
        contig_renaming.py ${flye_assembly} ${isolate_id}_polypolish1.fasta ${isolate_id}_inter.fasta ${isolate_id}_polypolish2.fasta
        seqkit sort --by-length --reverse ${isolate_id}_polypolish2.fasta > ${isolate_id}_polypolish.fasta
        rm *sam
        """
}

process POLCA {
        label "short_job"
        conda "${params.polca_env}"
        publishDir path:("${params.outdir}"), mode: 'copy', saveAs: {filename -> "${isolate_id}_polca.fasta"}, pattern: '*polca.fasta'        


        input:
        tuple val(isolate_id), path(assembly, stageAs: "input_assembly.fasta"), path(illumina1), path(illumina2), path(flye_assembly, stageAs: "flye_input.fasta")
        output:
        tuple val(isolate_id), path("${isolate_id}_polca.fasta"), emit: assembly_fasta

        script:
        """
        polca.sh -a ${assembly} -r '${illumina1} ${illumina2}' -t ${params.threads} -m 4G
        mv input_assembly.fasta.PolcaCorrected.fa ${isolate_id}_polca_intermediate.fasta
        seqkit sort --by-length --reverse ${isolate_id}_polca_intermediate.fasta > ${isolate_id}_polca.fasta
        contig_renaming.py ${flye_assembly} ${isolate_id}_polca.fasta ${isolate_id}_inter.fasta ${isolate_id}_polca.fasta

        """
}
