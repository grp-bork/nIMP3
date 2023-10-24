#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
// include { fastq_input } from "./nevermore/workflows/input"
include { metaT_input; metaG_input } from "./imp/workflows/input"

include { rnaspades; metaspades } from "./imp/modules/assemblers/spades"
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
	// [[id:sample1.metaT, library:paired, library_type:metaT], [/scratch/schudoma/imp3_test/work/0a/c32c5be6621089f9c71d87fc2fd308/no_host/sample1.metaT/sample1.metaT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/0a/c32c5be6621089f9c71d87fc2fd308/no_host/sample1.metaT/sample1.metaT_R2.fastq.gz, /scratch/schudoma/imp3_test/work/f3/2d793cfd78eda426ffc94ce4f5712a/merged/sample1.metaT.singles_R1.fastq.gz]]
	metaT_assembly.out.final_contigs.view()
	// [[id:sample1.metaT, library:paired, library_type:metaT], /scratch/schudoma/imp3_test/work/29/885cb85b918cada4ee1e07111a2434/assemblies/rnaspades/final/metaT/sample1.metaT/sample1.metaT.final_contigs.fasta]

	assembly_prep(
		nevermore_main.out.fastqs
			.filter { it[0].library_type == "metaG" }
	)

	metaG_assembly_ch = assembly_prep.out.reads
	metaG_assembly_ch.view()

	hybrid_assembly_input_ch = metaT_assembly.out.reads
		.map { sample, fastqs ->
			meta = [:]
			meta.id = sample.id.replaceAll(/\.metaT/, "") 
			// return tuple(meta.id, meta, fastqs)
			
			return tuple(meta, fastqs)
			
		}
		.concat(
			metaG_assembly_ch
				.map { sample, fastqs ->
					meta = [:]
					meta.id = sample.id.replaceAll(/\.metaG/, "")
					// return tuple(meta.id, meta, fastqs)
					// sample.id = sample.id.replaceAll(/\.metaT/, "")
					return tuple(meta, fastqs)

				}			
		)
		.groupTuple()
		.map { sample, fastqs -> return tuple(sample, fastqs.flatten()) }
		.concat(
			metaT_assembly.out.final_contigs
				.map { sample, contigs ->
					meta = [:]
					meta.id = sample.id.replaceAll(/\.metaT/, "")
					return tuple(meta, contigs)
				},			
		)
		.groupTuple()
		.map { sample, data -> return tuple(sample, data[0], data[1]) }

	hybrid_assembly_input_ch.view()

	metaspades(hybrid_assembly_input_ch)

	metaspades.out.contigs.view()



}

