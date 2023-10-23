process bwa_index {

	input:
	tuple val(sample), path(fasta)
	val(stage)

	output:
	tuple val(sample), path("index/${stage}/${sample.id}/${sample.id}*"), emit: index

	script:
	"""
	bwa index -p index/${stage}/${sample.id}/${sample.id} ${fasta}
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