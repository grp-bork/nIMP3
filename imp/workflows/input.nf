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
					sample.library_type = "metaT"					
				return tuple(sample, files)
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
					sample.library_type = "metaG"
				return tuple(sample, files)
			}
			
}


workflow assembly_prep {
	take:
		fastq_ch
	main:
		initial_assembly_ch = fastq_ch
			.map { sample, fastqs -> 
				def meta = sample.clone()
				meta.id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				
				return tuple(meta, [fastqs].flatten())
			}
			.groupTuple(by: 0, size: 2, remainder: true)
			.map { sample, fastqs -> return tuple(sample, fastqs.flatten())}
	emit:
		reads = initial_assembly_ch
}

