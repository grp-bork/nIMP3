include { extract_unmapped } from "../modules/alignment/extract"

workflow get_unmapped_reads {
	take:
		fastq_ch
	main:
		post_assembly_check_ch = fastq_ch
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
		reads = unmapped_ch
}