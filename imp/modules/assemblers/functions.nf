process concatenate_contigs {
	input:
		tuple val(sample), path(contigs)
		val(stage)
		val(assembler)
		
	output:
		tuple val(sample), path("assemblies/${assembler}/${stage}/${sample.library_source}/${sample.id}/${sample.id}*"), emit: contigs
	
	script:
		"""
		mkdir -p assemblies/${assembler}/${stage}/${sample.library_source}/${sample.id}/

		cat <(rename_contigs.awk -v ctg_prefix=preprocessing ${contigs[0]}) <(rename_contigs.awk -v ctg_prefix=unmapped ${contigs[1]}) > assemblies/${assembler}/${stage}/${sample.library_source}/${sample.id}/${sample.id}.final_contigs.fasta
		"""

}


def filter_fastq(fastq, is_paired, library_source) {

	def metaT_p = ~/.*metaT.*/
	def metaG_p = ~/.*metaG.*/
	def se_p = ~/.*(singles|orphans|chimeras).*/

	return (
		library_source == "metaG"
			? metaG_p.matcher(fastq.name).matches()
			: metaT_p.matcher(fastq.name).matches()
		) && (
			is_paired
				? !se_p.matcher(fastq.name).matches()
				: se_p.matcher(fastq.name).matches()
		)

}