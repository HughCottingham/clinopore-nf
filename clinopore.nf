nextflow.enable.dsl=2

log.info """\
         CLINOPORE
         ===================================
         reads        : ${params.reads}
         """
         .stripIndent()


def check_arguments(params) {
  // Check required input and outputs
  if (! params.reads) {
    exit 1, "ERROR: option 'reads' must be set in nextflow.config or on the command line using --reads"
  }
  if (! params.outdir) {
    exit 1, "ERROR: option 'outdir' must be set in nextflow.config or on the command line using --outdir"
  }
  if (! params.clinopore_env) {
    exit 1, "ERROR: option 'clinopore_env' must be set in nextflow.config or on the command line using --clinopore_env"
  }
  if (! params.polca_env) {
    exit 1, "ERROR: option 'polca_env' must be set in nextflow.config or on the command line using --polca_env"
  }
}

check_arguments(params)

def get_read_prefix_and_type(filepath) {
  // NOTE: nf requires escaping '$'
  regex_se = "^(.+?).fastq(?:.gz)?\$"
  regex_pe = "^(.+?)_R?[12](?:_001)?.fastq(?:.gz)?\$"
  //regex_both = "^(.+?)(_[12])?.fastq(?:.gz)?\$"

  java.util.regex.Matcher matcher;
  String read_type;

  if ((matcher = (filepath.getName() =~ /$regex_pe/))) {
    read_type = 'pe'
  } else if ((matcher = (filepath.getName() =~ /$regex_se/))) {
    read_type = 'se'
  } else {
    exit 1, "ERROR: did not find any readsets with the provided glob: ${filepath}"
  }
  return [read_type, matcher.group(1), filepath]
}



// Get reads and seperate into pe and se channels based on prefix
reads = Channel.fromPath(params.reads).ifEmpty {
    exit 1, "ERROR: did not find any read files with '${params.reads}'"
  }.map {
    get_read_prefix_and_type(it)
  }.branch {
    paired: it[0] == 'pe'
    single: it[0] == 'se'
}
reads_se = reads.single.map { it[1..-1] }.groupTuple()
reads_pe = reads.paired.map { it[1..-1] }.groupTuple()
// Check that we have the expected number of reads for each prefix in pe and se channels and flatten tuple
reads_pe = reads_pe.map {
  if (it[1].size() != 2) {
    exit 1, "ERROR: didn't get exactly two readsets prefixed with ${it[0]}:\n${it[1]}"
  }
  [it[0], *it[1]]
}
reads_se = reads_se.map {
  if (it[1].size() != 1) {
    exit 1, "ERROR: didn't get exactly one readset prefixed with ${it[0]}:\n${it[1]}"
  }
  [it[0], *it[1]]
}
//reads_se.view()
//reads_pe.view()

//combined_ch = reads_se
    //.combine(reads_pe, by: 0)



def check_host(workflow) {
  // Do not run on MASSIVE unless user specifies profile to use to avoid inadvertently using a local executor
  massive_hostnames = ['m3-login1', 'm3-login2']
  on_massive = massive_hostnames.contains(InetAddress.getLocalHost().getHostName())
  profile_explicit = workflow.commandLine.tokenize(' ').contains('-profile')
  if (on_massive && ! profile_explicit) {
    exit 1, "ERROR: to run on MASSIVE you must explicitly set -profile massive"
  }
}


// For an optional stage param variable, check that it is either a Boolean or String
// If it is a string and either 'true' or 'false', return the boolean equivalent
def check_boolean_option(option, name) {
  if (option.getClass() == java.lang.Boolean) {
    return option
  } else if (option.getClass() == java.lang.String) {
    if (option.toLowerCase() == 'true') {
      return true
    } else if (option.toLowerCase() == 'false') {
      return false
    }
  }
  exit 1, "ERROR: ${name} option must be true or false"
}

//run_rasusa = check_boolean_option(params.run_rasusa, 'run_rasusa')
run_filtlong = check_boolean_option(params.run_filtlong, 'run_filtlong')
run_medaka = check_boolean_option(params.run_medaka, 'run_medaka')
run_polypolish = check_boolean_option(params.run_polypolish, 'run_polypolish')
run_polca = check_boolean_option(params.run_polca, 'run_polca')



process FILTER {
        conda "${params.clinopore_env}"
        input:
        tuple val(isolate_id), path(reads_se)

        output:
        tuple val(isolate_id), path("${isolate_id}_filtered.fastq")

        script:
        """
        filtlong --min_length 1000 --keep_percent 95 $reads_se > ${isolate_id}_filtered.fastq
        """
}

process ASSEMBLE {
        conda "${params.clinopore_env}"
        publishDir path:("${params.outdir}/flye"), mode: 'copy', saveAs: {filename -> "${isolate_id}_flye.fasta"}, pattern: '*fasta'
        publishDir path:("${params.outdir}/gfa"), mode: 'copy', saveAs: {filename -> "${isolate_id}.gfa"}, pattern: '*gfa'
        input:
        tuple val(isolate_id), path(filtered_ch)

        output:
        tuple val(isolate_id), path("${isolate_id}_flye.fasta"), path("${isolate_id}.gfa")
        
        script:
        """
        flye --nano-hq ${isolate_id}_filtered.fastq --out-dir flye_out --threads 16
        reordering_contigs.py flye_out/assembly.fasta flye_out/${isolate_id}_inter1.fasta flye_out/${isolate_id}_inter2.fasta flye_out/assembly_graph.gfa ${isolate_id}.gfa flye_out/assembly_info.txt flye_out/${isolate_id}_assembly_info.txt ${isolate_id}_flye.fasta
        mv flye_out/assembly.fasta flye_out/${isolate_id}_old.fasta
        mv flye_out/assembly_graph.gfa flye_out/${isolate_id}_old.gfa
        mv flye_out/assembly_info.txt flye_out/${isolate_id}_assembly_info_old.txt
        """
}

process MEDAKA {
        conda "${params.clinopore_env}"
        publishDir path:("${params.outdir}/medaka"), mode: 'copy', saveAs: {filename -> "${isolate_id}_medaka.fasta"}, pattern: '*fasta'

        input:
        tuple val(isolate_id), path(assembly), path(gfa_file), path(filtered_reads)

        output:
        tuple val(isolate_id), path("${isolate_id}_medaka.fasta")

        script:
        """
        medaka_consensus -d ${assembly} -o . -i ${filtered_reads} -t 16 -m r941_min_sup_g507
        mv consensus.fasta ${isolate_id}_medaka_inter1.fasta
        contig_renaming.py ../../../${params.outdir}/flye/${isolate_id}_flye.fasta ${isolate_id}_medaka_inter1.fasta ${isolate_id}_medaka_inter2.fasta ${isolate_id}_medaka.fasta
        seqkit sort --by-length --reverse ${isolate_id}_medaka.fasta > ${isolate_id}_medaka.fasta

        """
}



process POLYPOLISH {
        conda "${params.clinopore_env}"
        publishDir path:("${params.outdir}/polypolish"), mode: 'copy', saveAs: {filename -> "${isolate_id}_medaka_polypolish.fasta"}, pattern: '*polypolish.fasta'        

        input:
        //tuple from medaka polished
        //tuple pe reads
        tuple val(isolate_id), path(medaka_polished_assembly), path(illumina1), path(illumina2)
        output:
        tuple val(isolate_id), path("${isolate_id}_medaka_polypolish.fasta")

        script:
        """
        bwa index ${medaka_polished_assembly}
        bwa mem -t 16 -a ${medaka_polished_assembly} ${illumina1} > ${isolate_id}_r1.sam
        bwa mem -t 16 -a ${medaka_polished_assembly} ${illumina2} > ${isolate_id}_r2.sam
        polypolish_insert_filter.py --in1 ${isolate_id}_r1.sam --in2 ${isolate_id}_r2.sam --out1 ${isolate_id}_filtered_r1.sam --out2 ${isolate_id}_filtered_r2.sam
        polypolish ${medaka_polished_assembly} ${isolate_id}_filtered_r1.sam ${isolate_id}_filtered_r2.sam| sed 's/_polypolish//' > ${isolate_id}_medaka_polypolish.fasta
        contig_renaming.py ../../../${params.outdir}/flye/${isolate_id}_flye.fasta ${isolate_id}_medaka_polypolish.fasta ${isolate_id}_inter.fasta ${isolate_id}_medaka_polypolish.fasta
        seqkit sort --by-length --reverse ${isolate_id}_medaka_polypolish.fasta > ${isolate_id}_medaka_polypolish.fasta
        """
}

process POLCA {
        conda "${params.polca_env}"
        publishDir path:("${params.outdir}"), mode: 'copy', saveAs: {filename -> "${isolate_id}_medaka_polypolish_polca.fasta"}, pattern: '*polca.fasta'        


        input:
        tuple val(isolate_id), path(polypolish_polished_assembly), path(illumina1), path(illumina2)
        output:
        tuple val(isolate_id), path("${isolate_id}_medaka_polypolish_polca.fasta")

        script:
        """
        polca.sh -a ${polypolish_polished_assembly} -r '${illumina1} ${illumina2}' -t 16 -m 4G
        mv ${isolate_id}_medaka_polypolish.fasta.PolcaCorrected.fa ${isolate_id}_medaka_polypolish_polca_intermediate.fasta
        seqkit sort --by-length --reverse ${isolate_id}_medaka_polypolish_polca_intermediate.fasta > ${isolate_id}_medaka_polypolish_polca.fasta
        contig_renaming.py ../../../${params.outdir}/flye/${isolate_id}_flye.fasta ${isolate_id}_medaka_polypolish_polca.fasta ${isolate_id}_inter.fasta ${isolate_id}_medaka_polypolish_polca.fasta

        """
}


workflow {
  if ( run_filtlong ) {
    filtered_ch=FILTER(reads_se)
  }
  else {
    filtered_ch=(reads_se)
  }
  assembled_ch=ASSEMBLE(filtered_ch)
  if ( run_medaka ) {
    flye_filtered_ch=assembled_ch
    .combine(filtered_ch,by:0)
    .view()
    medaka_polished_ch=MEDAKA(flye_filtered_ch)
  }

  if ( run_polypolish ) {
    medaka_illumina_ch=medaka_polished_ch
    .combine(reads_pe, by: 0)
    .view()
    polypolished_ch=POLYPOLISH(medaka_illumina_ch)
  }
  if ( run_polca ) {
    polypolish_illumina_ch=polypolished_ch
    .combine(reads_pe,by:0)
    .view()
    polca_polished_ch=POLCA(polypolish_illumina_ch)
  }
}
