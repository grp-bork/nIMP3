#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { metaT_input; metaG_input } from "./imp/workflows/input"
include { imp_main } from "./imp/workflows/imp"

params.remote_input_dir = false


workflow {

	metaT_input(
		// Channel.fromPath(params.meta_t_input + "/*", type: "dir", checkIfExists: true)
		Channel.fromPath("${params.meta_t_input}/**.{fastq,fastq.gz,fq,fq.gz}", checkIfExists: true)
	)
	metaG_input(
		// Channel.fromPath(params.meta_g_input + "/*", type: "dir", checkIfExists: true)
		Channel.fromPath("${params.meta_g_input}/**.{fastq,fastq.gz,fq,fq.gz}", checkIfExists: true)
	)

	metaT_ch = metaT_input.out.reads		
	metaG_ch = metaG_input.out.reads		

	nevermore_main(metaT_ch.mix(metaG_ch))

	nevermore_main.out.fastqs.dump(pretty: true, tag: "nvm_main_fastqs")

	imp_main(nevermore_main.out.fastqs)
	// 		.map { sample, files -> 
	// 			def new_sample = sample.clone()
	// 			new_sample.library_source = (new_sample.id.endsWith(".metaT") ? "metaT" : ((new_sample.id.endsWith(".metaG")) ? "metaG": null)) //library_source
	// 			return [ new_sample, files ]
	// 		}
	// )


	
}
