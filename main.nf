#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"

include { rnaspades } from "./imp/modules/assemblers/spades"
include { bwa_index } from "./imp/modules/alignment/indexing/bwa_index"
include { extract_unmapped } from "./imp/modules/alignment/extract"

// if (params.input_dir && params.remote_input_dir) {
// 	log.info """
// 		Cannot process both --input_dir and --remote_input_dir. Please check input parameters.
// 	""".stripIndent()
// 	exit 1
// } else if (!params.input_dir && !params.remote_input_dir) {
// 	log.info """
// 		Neither --input_dir nor --remote_input_dir set.
// 	""".stripIndent()
// 	exit 1
// }

def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir

// params.ignore_dirs = ""
params.remote_input_dir = false
// params.has_assay_suffix = true

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

		rnaspades(initial_assembly_ch, "initial")
		bwa_index(rnaspades.out.contigs, "initial")

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

		extract_unmapped(metaT_post_assembly_check_ch, "initial")

		unmapped_ch = extract_unmapped.out.fastqs
			.map { sample, fastqs ->
				sample.id = sample.index_id
				return tuple(sample, fastqs)
			}
			.groupTuple(by: 0, size: 2, remainder: true)
			.map { sample, fastqs -> return tuple(sample, fastqs.flatten()) }

		emit:
			unmapped_reads = unmapped_ch
			contigs = rnaspades.out.contigs
}

workflow metaT_assembly {
	take:
		fastq_ch
	main:
		metaT_initial_assembly(fastq_ch)
		metaT_initial_assembly.out.unmapped_reads.view()
		metaT_initial_assembly.out.contigs.view()
	emit:
		contigs = metaT_initial_assembly.out.contigs

}



workflow {

	metaT_input(
		Channel.fromPath(params.metaT_input_dir + "/*", type: "dir")
	)
	metaG_input(
		Channel.fromPath(params.metaG_input_dir + "/*", type: "dir")
	)

	metaT_ch = metaT_input.out.reads		
	metaG_ch = metaG_input.out.reads		

	nevermore_main(metaT_ch.concat(metaG_ch))

	metaT_assembly(nevermore_main.out.fastqs)
	metaT_assembly.out.contigs.view()

	// nevermore_main.out.fastqs.view()
	
	// metaT_preprocessed_ch = nevermore_main.out.fastqs
	// 	.filter { it[0].library_type == "metaT" }
	// 	.map { sample, fastqs -> 
	// 		meta = [:]
	// 		meta.id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
	// 		meta.library = sample.library
	// 		meta.library_type = sample.library_type
			
	// 		return tuple(meta, [fastqs].flatten())
	// 	}
	// 	.groupTuple(by: 0, size: 2, remainder: true)
	// 	.map { sample, fastqs -> return tuple(sample, fastqs.flatten())}
	
	// // metaT_preprocessed_ch.view()

	// rnaspades(metaT_preprocessed_ch, "initial")
	// bwa_index(rnaspades.out.contigs, "initial")

	// metaT_post_assembly_check_ch = nevermore_main.out.fastqs
	// 	.filter { it[0].library_type == "metaT" }
	// 	.map { sample, fastqs -> 
	// 		sample_base_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
	// 		return tuple(sample_base_id, sample.library_type, sample, [fastqs].flatten())
	// 	}
	// 	.combine(bwa_index.out.index, by: [0, 1])
	// 	.map { sample_id, libtype, sample, fastqs, index ->
	// 		sample.index_id = sample_id
	// 		return tuple(sample, fastqs, index) }
	
	// extract_unmapped(metaT_post_assembly_check_ch, "initial")
	// extract_unmapped.out.fastqs.view()

}
