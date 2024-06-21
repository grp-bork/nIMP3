process bwa_index {
	container "registry.git.embl.de/schudoma/align-docker:latest"
	label "small"

	input:
	tuple val(sample), path(fasta)
	val(stage)

	output:
	tuple val(sample.id), val(sample.library_source), path("index/${stage}/${sample.library_source}/${sample.id}/${sample.id}*"), emit: index

	script:
	"""
	mkdir -p index/${stage}/${sample.library_source}/${sample.id}/
	bwa index -p index/${stage}/${sample.library_source}/${sample.id}/${sample.id} ${fasta}
	"""


}


// rule bwa_index:
//     input:
//         "{fasta}"
//     output:
//         "{fasta}.amb",
//         "{fasta}.bwt",
//         "{fasta}.pac",
//         "{fasta}.sa",
//         "{fasta}.ann"
//     resources:
//         runtime = "24:00:00",
//         mem = MEMCORE
//     threads: 1
//     conda: ENVDIR + "/IMP_mapping.yaml"
// #    log: "logs/assembly_bwa_index.log"
//     message: "bwa_index: Indexing {wildcards.fasta} for bwa."
//     shell:
//         """
//         bwa index {wildcards.fasta} 
//         """