include { fastq_input } from "../../nevermore/workflows/input"


workflow metaT_input {
	take:
		fastq_ch
	main:
		fastq_input(fastq_ch, Channel.of("metaT"))
	emit:
		reads = fastq_input.out.fastqs
			.map {
				sample, files ->					
					def new_sample = sample.clone()
					new_sample.library_source = "metaT"					
				return tuple(new_sample, files)
			}
}


workflow metaG_input {
	take:
		fastq_ch
	main:
		fastq_input(fastq_ch, Channel.of("metaG"))
	emit:
		reads = fastq_input.out.fastqs
			.map {
				sample, files ->
					def new_sample = sample.clone()
					new_sample.library_source = "metaG"
				return tuple(new_sample, files)
			}
			
}

// collect all fastq files for a sample
// identify related fastqs via sample id : <prefix>.<library_source>[.singles]
// each sample has at most 2 groups of files: [2 x PE, 1 x orphan], [1 x singles]
workflow assembly_prep {
	take:
		fastq_ch
	main:
		initial_assembly_ch = fastq_ch
			.map { sample, fastqs -> 
				def new_sample = sample.clone()
				new_sample.id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				
				return tuple(new_sample.id, new_sample.library_source, [fastqs].flatten())
			}
			.groupTuple(by: [0, 1], size: 2, remainder: true)

		initial_assembly_ch.dump(pretty: true, tag: "in_assemply_prep")

		initial_assembly_ch = initial_assembly_ch
			.map { sample_id, sample_source, fastqs -> 
				def meta = [:]
				meta.id = sample_id
				meta.library_source = sample_source
				return tuple(meta, [fastqs].flatten())
			}
	emit:
		reads = initial_assembly_ch
}

