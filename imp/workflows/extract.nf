include { extract_unmapped } from "../modules/alignment/extract"


workflow get_unmapped_reads {
	take:
		fastq_ch
        index_ch
	main:
        // [[id:sample1.metaG, is_paired:true, library:paired, library_type:metaG, merged:true], [/scratch/schudoma/imp3_test/work/da/eea24ab29b720733d771235b1d7a15/no_host/sample1.metaG/sample1.metaG_R1.fastq.gz, /scratch/schudoma/imp3_test/work/da/eea24ab29b720733d771235b1d7a15/no_host/sample1.metaG/sample1.metaG_R2.fastq.gz]]
        // [[id:sample1.metaT, is_paired:true, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/32/89e706473152c55350201e65cd16a3/no_host/sample1.metaT/sample1.metaT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/32/89e706473152c55350201e65cd16a3/no_host/sample1.metaT/sample1.metaT_R2.fastq.gz]]
		post_assembly_check_ch = fastq_ch
			.map { sample, fastqs -> 
				sample_base_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_base_id, sample.library_type, sample, [fastqs].flatten())
			}
			.combine(index_ch, by: [0, 1])
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
				def meta = sample.clone()
				meta.id = sample_id
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