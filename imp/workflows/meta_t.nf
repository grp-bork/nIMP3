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

		if (params.assembler == "spades") {
			rnaspades(initial_assembly_ch, "initial")
			contigs_ch = rnaspades.out.contigs
		} else {
			metaT_megahit(initial_assembly_ch, "initial")
			contigs_ch = metaT_megahit.out.contigs
		}
		bwa_index(contigs_ch, "initial")

		post_assembly_check_ch = fastq_ch
			.map { sample, fastqs -> 
				sample_base_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_base_id, sample.library_source, sample, [fastqs].flatten())
			}
			.combine(bwa_index.out.index, by: [0, 1])
			.map { sample_id, libsrc, sample, fastqs, index ->
				def new_sample = sample.clone()
				new_sample.index_id = sample_id
				return tuple(new_sample, fastqs, index) 
			}

		post_assembly_check_ch.dump(pretty: true, tag: "post_assembly_check_ch")

		post_assembly_check_ch
			.filter { it[0].is_paired }
			.dump(pretty: true, tag: "post_assembly_check_ch_is_paired")
		post_assembly_check_ch.subscribe { println it[0].values() }
		

		extract_unmapped(post_assembly_check_ch, "initial")

		unmapped_ch = extract_unmapped.out.fastqs
			.map { sample, fastqs ->
				def new_sample = sample.clone()
				new_sample.id = sample.index_id
				return tuple(new_sample, [fastqs].flatten())
			}
			.groupTuple(by: 0, size: 2, remainder: true, sort: true)
		unmapped_ch.dump(pretty: true, tag: "unmapped_after_metaT_assembly_1")

		unmapped_ch = unmapped_ch
			.map { sample, fastqs ->
				def new_sample = sample.clone()
				new_sample.library = sample.library[0]
				new_sample.library_source = sample.library_source[0]
				return tuple(new_sample, [fastqs].flatten())
			}

		emit:
			unmapped_reads = unmapped_ch
			contigs = contigs_ch
			reads = initial_assembly_ch
}


workflow metaT_assembly {
	take:
		fastq_ch
	main:
		metaT_initial_assembly(fastq_ch)

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

		metaT_initial_assembly.out.contigs.dump(pretty: true, tag: "metaT_initial_assembly.out.contigs")

		all_contigs_ch = metaT_initial_assembly.out.contigs
			.concat(contigs_ch)
			.groupTuple(by: 0, size: 2, remainder: true, sort: true)

		all_contigs_ch.dump(pretty: true, tag: "all_contigs_ch")

		concatenate_contigs(all_contigs_ch, "final", assembler)

	emit:
		initial_contigs = metaT_initial_assembly.out.contigs
		unmapped_contigs = contigs_ch
		final_contigs = concatenate_contigs.out.contigs
		reads = metaT_initial_assembly.out.reads

}
