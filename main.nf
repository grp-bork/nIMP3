#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
// include { fastq_input } from "./nevermore/workflows/input"
include { metaT_input; metaG_input } from "./imp/workflows/input"

include { rnaspades } from "./imp/modules/assemblers/spades"
include { bwa_index } from "./imp/modules/alignment/indexing/bwa_index"
include { extract_unmapped } from "./imp/modules/alignment/extract"

include { metaT_assembly } from "./imp/workflows/meta_t"
include { assembly_prep } from "./imp/workflows/input"

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

params.remote_input_dir = false


workflow {
	// last working revision: 06e468bf12

	metaT_input(
		Channel.fromPath(params.metaT_input_dir + "/*", type: "dir")
	)
	metaG_input(
		Channel.fromPath(params.metaG_input_dir + "/*", type: "dir")
	)

	metaT_ch = metaT_input.out.reads		
	metaG_ch = metaG_input.out.reads		

	nevermore_main(metaT_ch.concat(metaG_ch))

	metaT_assembly(
		nevermore_main.out.fastqs
			.filter { it[0].library_type == "metaT" }			
	)

	metaT_assembly.out.reads.view()
	metaT_assembly.out.final_contigs.view()

	assembly_prep(
		nevermore_main.out.fastqs
			.filter { it[0].library_type == "metaG" }
	)

	metaG_assembly_ch = assembly_prep.out.reads
	metaG_assembly_ch.view()
}

