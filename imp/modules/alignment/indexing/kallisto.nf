process kallisto_index {
	container "quay.io/biocontainers/kallisto:0.50.1--hc877fd6_1"
	label "medium"
	tag "${sample.id}.${stage}"

	input:
	tuple val(sample), path(fasta)
	val(stage)

	output:
	tuple val(sample), path("kallisto/index/${stage}/${sample.library_source}/${sample.id}/${sample.id}.idx"), emit: index

	script:
	"""
	mkdir -p kallisto/index/${stage}/${sample.library_source}/${sample.id}/
	kallisto index -i index/${stage}/${sample.library_source}/${sample.id}/${sample.id}.idx ${fasta}
	"""
	
}
