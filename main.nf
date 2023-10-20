#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"

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

workflow metaT_input {
	take:
		fastq_ch
	main:
		fastq_input(fastq_ch)
	emit:
		reads = fastq_input.out.fastqs
			.map {
				sample, files ->
					new_sample = sample.clone()
					new_sample.library_type = "metaT"
					new_sample.id = new_sample.id + ".MT"
				return tuple(new_sample, files)
			}
}

workflow metaG_input {
	take:
		fastq_ch
	main:
		fastq_input(fastq_ch)
	emit:
		reads = fastq_input.out.fastqs
			.map {
				sample, files ->
					new_sample = sample.clone()
					new_sample.library_type = "metaG"
					new_sample.id = new_sample.id + ".MG"
				return tuple(new_sample, files)
			}
			
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

	nevermore_main.out.fastqs.view()
	

	// // this is for later
	// joined = metaT_ch.concat(metaG_ch).groupTuple(by:0, size: 2, remainder: true)
	// joined.view()

}
