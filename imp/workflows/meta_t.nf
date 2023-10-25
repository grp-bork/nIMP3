include { rnaspades } from "../modules/assemblers/spades"
include { metaT_megahit } from "../modules/assemblers/megahit"
include { bwa_index } from "../modules/alignment/indexing/bwa_index"
include { extract_unmapped } from "../modules/alignment/extract"
include { concatenate_contigs } from "../modules/assemblers/functions"

include { assembly_prep } from "./input"
include { get_unmapped_reads } from "./extract"


workflow metaT_initial_assembly {
	take:
		fastq_ch
	main:
		assembly_prep(fastq_ch)
		initial_assembly_ch = assembly_prep.out.reads
		// initial_assembly_ch = fastq_ch
		// 	.map { sample, fastqs -> 
		// 		meta = [:]
		// 		meta.id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
		// 		meta.library = sample.library
		// 		meta.library_type = sample.library_type
				
		// 		return tuple(meta, [fastqs].flatten())
		// 	}
		// 	.groupTuple(by: 0, size: 2, remainder: true)
		// 	.map { sample, fastqs -> return tuple(sample, fastqs.flatten())}

		if (params.assembler == "spades") {
			rnaspades(initial_assembly_ch, "initial")
			contigs_ch = rnaspades.out.contigs
		} else {
			metaT_megahit(initial_assembly_ch, "initial")
			contigs_ch = metaT_megahit.out.contigs
		}
		bwa_index(contigs_ch, "initial")

		// get_unmapped_reads(fastq_ch, bwa_index.out.index)

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
			// unmapped_reads = get_unmapped_reads.out.reads
			unmapped_reads = unmapped_ch
			contigs = contigs_ch
			reads = initial_assembly_ch
}


workflow metaT_assembly {
	take:
		fastq_ch
	main:
		metaT_initial_assembly(fastq_ch)
		metaT_initial_assembly.out.unmapped_reads.view()

		def assembler = ""
		if (params.assembler == "spades") {
			rnaspades(metaT_initial_assembly.out.unmapped_reads, "unmapped")
			contigs_ch = rnaspades.out.contigs
			assembler = "rnaspades"
		} else {
			metaT_megahit(metaT_initial_assembly.out.unmapped_reads, "unmapped")
			contigs_ch = metaT_megahit.out.contigs
			assembler = "metaT_megahit"
		}

		all_contigs_ch = metaT_initial_assembly.out.contigs
			.concat(contigs_ch)
			.groupTuple(by: 0, size: 2, remainder: true, sort: true)
		concatenate_contigs(all_contigs_ch, "final", assembler)

	emit:
		initial_contigs = metaT_initial_assembly.out.contigs
		unmapped_contigs = contigs_ch
		final_contigs = concatenate_contigs.out.contigs
		reads = metaT_initial_assembly.out.reads

}
