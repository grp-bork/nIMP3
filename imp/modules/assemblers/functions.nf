process concatenate_contigs {
	input:
		tuple val(sample), path(contigs)
		val(stage)
		val(assembler)
		
	output:
		tuple val(sample), path("assemblies/${assembler}/${stage}/${sample.library_type}/${sample.id}/${sample.id}*"), emit: contigs
	
	script:
		"""
		mkdir -p assemblies/${assembler}/${stage}/${sample.library_type}/${sample.id}/

		cat <(awk -f rename_contigs.awk -v preprocessing ${contigs[0]}) <(awk -f rename_contigs.awk -v unmapped ${contigs[1]}) > assemblies/${assembler}/${stage}/${sample.library_type}/${sample.id}/${sample.id}.final_contigs.fasta
		"""

}