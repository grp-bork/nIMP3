include { rnaspades; metaspades } from "../modules/assemblers/spades"
include { quast } from "../modules/assemblers/quast"
include { bwa_index } from "../modules/alignment/indexing/bwa_index"
include { extract_unmapped } from "../modules/alignment/extract"

include { metaT_assembly } from "./meta_t"
include { assembly_prep } from "./input"
include { hybrid_megahit } from "../modules/assemblers/megahit"
include { get_unmapped_reads } from "./extract"
include { concatenate_contigs; filter_fastq } from "../modules/assemblers/functions"

include { kallisto_index; kallisto_quant } from "../modules/alignment/indexing/kallisto"


params.assembler = "megahit"


workflow imp_main {
	take:
		fastq_ch

	main:	
		metaT_assembly(
			fastq_ch
				.filter { it[0].library_source == "metaT" }			
		)
		
		// collect metaG fastqs per sample
		assembly_prep(
			fastq_ch
				.filter { it[0].library_source == "metaG" }
		)

		metaG_assembly_ch = assembly_prep.out.reads

		metaG_assembly_ch.dump(pretty: true, tag: "metaG_hybrid_input")

		// assign proper sample labels to metaT contigs
		metaT_contigs_ch = metaT_assembly.out.final_contigs
			.map { sample, contigs ->
				def meta = [:]
				meta.id = sample.id.replaceAll(/\.singles$/, "").replaceAll(/\.metaT/, "")
				return [ meta.id, meta, contigs ]
			}
		metaT_contigs_ch.dump(pretty: true, tag: "metaT_contigs_ch")

		// group metaT files by sample id
		hybrid_assembly_input_ch = metaT_assembly.out.reads
			.map { sample, fastqs ->
				def meta = [:]
				meta.id = sample.id.replaceAll(/\.singles$/, "").replaceAll(/\.metaT/, "")
				return [ meta, fastqs ]				
			}
		hybrid_assembly_input_ch.dump(pretty: true, tag: "metaT_hybrid_input")

		// combine the metaT and metaG reads
		hybrid_assembly_input_ch = hybrid_assembly_input_ch
			.mix(
				metaG_assembly_ch
					.map { sample, fastqs ->
						def meta = [:]
						meta.id = sample.id.replaceAll(/\.singles$/, "").replaceAll(/\.metaG/, "")
						return [ meta, fastqs ]
					}			
			)
			.map { sample, fastqs -> [ sample.id, sample, fastqs] }
			.groupTuple(by: 0, size: 4, remainder: true)
			.map { sample_id, sample, fastqs -> [ sample_id, sample[0], fastqs.flatten() ] }

		hybrid_assembly_input_ch.dump(pretty: true, tag: "all_reads_hybrid_input")

		// add the metaT contigs to the metaG/T input reads
		/*
		hybrid_assembly_input_ch = hybrid_assembly_input_ch
			.mix(
				metaT_contigs_ch
			)
			.groupTuple(by: 0, size: 2)
			.map { sample_id, sample, data -> [ sample[0], data[0], data[1] ] }
		*/

		hybrid_assembly_input_ch = hybrid_assembly_input_ch
			.join(metaT_contigs_ch.map { sample_id, sample, contigs -> [ sample_id, contigs ] }, by: 0)
			.map { sample_id, sample, fastqs, contigs -> [ sample, fastqs, contigs ] }

		hybrid_assembly_input_ch.dump(pretty: true, tag: "hybrid_assembly_input_ch")




		// perform initial hybrid assembly, label the resulting contigs as hybrid and build bwa index
		if (params.assembler == "spades") {
			metaspades(hybrid_assembly_input_ch, "initial")
			contigs_ch = metaspades.out.contigs		
		} else {
			hybrid_megahit(hybrid_assembly_input_ch, "initial")
			contigs_ch = hybrid_megahit.out.contigs
		}

		contigs_ch = contigs_ch.map {
			sample, fastqs -> 
			def new_sample = sample.clone()
			new_sample.library_source = "hybrid"
			return [ new_sample, fastqs ]
		}

		bwa_index(contigs_ch, "initial")

		bwa_index.out.index.dump(pretty: true, tag: "bwa_index.out.index")

		kallisto_index(contigs_ch, "initial")


		fastq_ch.dump(pretty: true, tag: "fastq_ch")

		// add the bwa indices to the input reads
		combined_assembly_input_index_ch = hybrid_assembly_input_ch
			.map { sample, fastqs, contigs -> [ sample.id, sample, fastqs ] }
			// .map { sample, fastqs, contigs -> [ sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "").replaceAll(/.meta[GT]$/, ""), sample, fastqs ] }

			
			// .join(bwa_index.out.index.map { sample, index -> [ sample.id, index ] }, by: 0)
			.join(bwa_index.out.index, by: 0)
			.map { sample_id, sample, fastqs, libsrc, index -> [ sample_id, sample, fastqs, index ] }
		combined_assembly_input_index_ch.dump(pretty: true, tag: "combined_assembly_input_index_ch")


		metaT_paired_unmapped_ch = combined_assembly_input_index_ch
			.map { sample_id, sample, fastqs, index ->
				def new_sample = [:]
				new_sample.id = sample.id + ".metaT"
				new_sample.library_source = "metaT"
				new_sample.is_paired = true
				new_sample.index_id = sample_id
				def wanted_fastqs = fastqs
					.findAll({ filter_fastq(it, true, "metaT") })
				return  [ new_sample, wanted_fastqs, index ]
			}
			.filter { it[1].size() > 0 }
		metaT_single_unmapped_ch = combined_assembly_input_index_ch
			.map { sample_id, sample, fastqs, index ->
				def new_sample = [:]
				new_sample.id = sample.id + ".metaT.singles"
				new_sample.library_source = "metaT"
				new_sample.is_paired = false
				new_sample.index_id = sample_id
				def wanted_fastqs = fastqs
					.findAll({ filter_fastq(it, false, "metaT") })
				return [ new_sample, wanted_fastqs, index ]
			}
			.filter { it[1].size() > 0 }
		metaG_paired_unmapped_ch = combined_assembly_input_index_ch
			.map { sample_id, sample, fastqs, index ->
				def new_sample = [:]
				new_sample.id = sample.id + ".metaG"
				new_sample.library_source = "metaG"
				new_sample.is_paired = true
				new_sample.index_id = sample_id
				def wanted_fastqs = fastqs
					.findAll({ filter_fastq(it, true, "metaG") })
				return [ new_sample, wanted_fastqs, index ]
			}
			.filter { it[1].size() > 0 }
		metaG_single_unmapped_ch = combined_assembly_input_index_ch
			.map { sample_id, sample, fastqs, index ->
				def new_sample = [:]
				new_sample.id = sample.id + ".metaG.singles"
				new_sample.library_source = "metaG"
				new_sample.is_paired = false
				new_sample.index_id = sample_id
				def wanted_fastqs = fastqs
					.findAll({ filter_fastq(it, false, "metaG") })
				return [ new_sample, wanted_fastqs, index ]
			}
			.filter { it[1].size() > 0 }

		extract_unmapped_ch = Channel.empty()
			.mix(metaT_paired_unmapped_ch)
			.mix(metaT_single_unmapped_ch)
			.mix(metaG_paired_unmapped_ch)
			.mix(metaG_single_unmapped_ch)
		
		extract_unmapped_ch.dump(pretty: true, tag: "extract_unmapped_ch")

		base_id_ch = fastq_ch
			.map { sample, fastqs -> 
				def sample_base_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "").replaceAll(/.meta[GT]$/, "")
				return [ sample_base_id, sample, [fastqs].flatten() ]
			}
		
		base_id_ch.dump(pretty: true, tag: "base_id_ch")

		index_and_fastqs_ch = bwa_index.out.index.combine(base_id_ch, by: 0)
		index_and_fastqs_ch.dump(pretty: true, tag: "index_and_fastqs_ch")

		with_index_ch = base_id_ch.combine(bwa_index.out.index)
		with_index_ch.dump(pretty: true, tag: "with_index_ch")

		joined_ch = base_id_ch.join(bwa_index.out.index, by: 0)
		joined_ch.dump(pretty: true, tag: "joined_ch")

		// base_id_ch.combine(bwa_index.out.index, by: 0).dump(pretty: true, tag: "base_id_ch")

		// post_assembly_check_ch = fastq_ch
		// 	.map { sample, fastqs -> 
		// 		sample_base_id = sample.id //
		// 		sample_base_id = sample_base_id.replaceAll(/.(orphans|singles|chimeras)$/, "").replaceAll(/.meta[GT]$/, "")
		// 		return tuple(sample_base_id, sample, [fastqs].flatten())
		// 	}
		// post_assembly_check_ch = with_index_ch
		// 	.map { sample_id, sample, fastqs, slib, index ->
		// 		sample.index_id = sample_id
		// 		return tuple(sample, fastqs, index) 
		// 	}

		// post_assembly_check_ch.dump(pretty: true, tag: "post_assembly_check_ch")
		extract_unmapped(extract_unmapped_ch, "initial")
		extract_unmapped.out.fastqs.dump(pretty: true, tag: "extract_unmapped_fastqs_ch")
		
		unmapped_ch = extract_unmapped.out.fastqs
			.map { sample, fastqs ->
				def new_sample = sample.clone()
				new_sample.id = sample.index_id
				return [ new_sample.id, new_sample, fastqs ]
			}
			.groupTuple(by: 0, size: 2, remainder: true)
			.map { sample_id, sample, fastqs -> 
				def meta = [:]
				meta.id = sample_id
				return [ meta.id, meta, fastqs.flatten() ]
			}
			.groupTuple(by: 0, size: 2, remainder: true) //, sort: true)
			.map { sample_id, sample, fastqs ->
				def new_sample = sample[0].clone()
				new_sample.library_source = "hybrid"
				return [ new_sample, fastqs.flatten() ]
			}

		unmapped_ch.dump(pretty: true, tag: "unmapped_ch")
		
		empty_file = file("${params.output_dir}/NO_INPUT")
		empty_file.text = "NOTHING TO SEE HERE."
		print empty_file


		unmapped_contigs_ch = Channel.empty()
		if (params.assembler == "megahit") {
			megahit_hybrid_unmapped(unmapped_ch, Channel.of(empty_file))
			unmapped_contigs_ch = megahit_hybrid_unmapped.out.contigs
		}
		unmapped_contigs_ch.dump(pretty: true, tag: "unmapped_contigs_ch")

		all_contigs_ch = contigs_ch
			.mix(unmapped_contigs_ch)
			.groupTuple(by: 0, size: 2, remainder: true, sort: true)

		all_contigs_ch.dump(pretty: true, tag: "all_contigs_ch")
		
		concatenate_contigs(all_contigs_ch, "final", params.assembler)
		// hybrid_megahit(unmapped_ch.combine(Channel.of(empty_file)))
		// final_assembly_ch = get_unmapped_reads.out.reads
		// 	.map { sample, fastqs -> return tuple(sample, fastqs, [empty_file])}
		// final_assembly_ch.view()

		quast(concatenate_contigs.out.contigs)
		
	
}

workflow megahit_hybrid_unmapped {
	take:
		fastq_ch
		contigs_ch
	main:
		hybrid_megahit(fastq_ch.combine(contigs_ch), "final")
	emit:
		contigs = hybrid_megahit.out.contigs
}
