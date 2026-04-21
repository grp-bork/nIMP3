process kallisto_index {
	container "quay.io/biocontainers/kallisto:0.52.0--h13ff97a_0"
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
	kallisto index -i kallisto/index/${stage}/${sample.library_source}/${sample.id}/${sample.id}.idx ${fasta}
	"""
	
}
