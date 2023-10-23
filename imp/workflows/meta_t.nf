include { spades } from "../modules/assemblers/spades"
include { bwa_index } from "../modules/alignment/indexing/bwa_index"
include { extract_unmapped } from "../modules/alignment/extract"
include { concatenate_contigs } from "../modules/assemblers/functions"


workflow metaT_initial_assembly {
	take:
		fastq_ch
	main:
		initial_assembly_ch = fastq_ch
			.filter { it[0].library_type == "metaT" }
			.map { sample, fastqs -> 
				meta = [:]
				meta.id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				meta.library = sample.library
				meta.library_type = sample.library_type
				
				return tuple(meta, [fastqs].flatten())
			}
			.groupTuple(by: 0, size: 2, remainder: true)
			.map { sample, fastqs -> return tuple(sample, fastqs.flatten())}

		spades(initial_assembly_ch, "initial")
		bwa_index(spades.out.contigs, "initial")

		post_assembly_check_ch = fastq_ch
			.filter { it[0].library_type == "metaT" }
			.map { sample, fastqs -> 
				sample_base_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_base_id, sample.library_type, sample, [fastqs].flatten())
			}
			.combine(bwa_index.out.index, by: [0, 1])
			.map { sample_id, libtype, sample, fastqs, index ->
				sample.index_id = sample_id
				return tuple(sample, fastqs, index) 
			}

		extract_unmapped(post_assembly_check_ch, "initial")

		unmapped_ch = extract_unmapped.out.fastqs
			.map { sample, fastqs ->
				sample.id = sample.index_id
				return tuple(sample.id, sample, fastqs)
			}
			.groupTuple(by: 0, size: 2, remainder: true)
			.map { sample_id, sample, fastqs -> 
				meta = [:]
				meta.id = sample_id
				meta.library = sample.library
				meta.library_type = sample.library_type
				return tuple(meta, fastqs.flatten())
			}
			.groupTuple(by: 0, size: 2, remainder: true, sort: true)
			.map { sample, fastqs ->
				sample.library = sample.library[0]
				sample.library_type = sample.library_type[0]
				return tuple(sample, fastqs.flatten())
			}

		emit:
			unmapped_reads = unmapped_ch
			contigs = spades.out.contigs
}


workflow metaT_assembly {
	take:
		fastq_ch
	main:
		metaT_initial_assembly(fastq_ch)
		metaT_initial_assembly.out.unmapped_reads.view()
		spades(metaT_initial_assembly.out.unmapped_reads, "unmapped")

		contigs_ch = metaT_initial_assembly.out.contigs
			.concat(spades.out.contigs)
			.groupTuple(by: 0, size: 2, remainder: true, sort: true)
		concatenate_contigs(contigs_ch, "final", "spades")

	emit:
		initial_contigs = metaT_initial_assembly.out.contigs
		unmapped_contigs = spades.out.contigs
		final_contigs = concatenate_contigs.out.contigs

}
