#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"

include { rnaspades } from "./imp/modules/assemblers/spades"

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
		fastq_input(fastq_ch, Channel.of("MT"))
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
		fastq_input(fastq_ch, Channel.of("MG"))
	emit:
		reads = fastq_input.out.fastqs
			.map {
				sample, files ->
					sample.library_type = "metaG"
				return tuple(sample, files)
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

	// nevermore_main.out.fastqs.view()
	
	// [[id:sample1.MG, is_paired:true, library:paired, library_type:metaG, merged:true], [/scratch/schudoma/imp3_test/work/b3/ce455386c39117375509edef448095/no_host/sample1.MG/sample1.MG_R1.fastq.gz, /scratch/schudoma/imp3_test/work/b3/ce455386c39117375509edef448095/no_host/sample1.MG/sample1.MG_R2.fastq.gz]]
	// [[id:sample1.MT, is_paired:true, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/fe/9d2f24508e4c6dc8bfba488ae5585e/no_host/sample1.MT/sample1.MT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/fe/9d2f24508e4c6dc8bfba488ae5585e/no_host/sample1.MT/sample1.MT_R2.fastq.gz]]

	// [[id:sample1.MT.singles, is_paired:false, library:paired, library_type:metaT, merged:true], /scratch/schudoma/imp3_test/work/bc/8c33cc37b6b62f57dea67ba440ae08/merged/sample1.MT.singles_R1.fastq.gz]
	// [[id:sample1.MG.singles, is_paired:false, library:paired, library_type:metaG, merged:true], /scratch/schudoma/imp3_test/work/73/09c4dbcb47693deddabd884990b5d8/merged/sample1.MG.singles_R1.fastq.gz]

	metaT_preprocessed_ch = nevermore_main.out.fastqs
		.filter { it[0].library_type == "metaT" }
		.map { sample, fastqs -> 
			meta = [:]
			meta.sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
			meta.library = sample.library
			meta.library_type = sample.library_type
			
			return tuple(meta, [fastqs].flatten())
		}
		.groupTuple(by: 0, size: 2, remainder: true)
		.map { sample, fastqs -> return tuple(sample, fastqs.flatten())}
	
	metaT_preprocessed_ch.view()

	// [
	// 	[sample_id:sample1.MT, library:paired, library_type:metaT],
	// 	[/scratch/schudoma/imp3_test/work/fe/9d2f24508e4c6dc8bfba488ae5585e/no_host/sample1.MT/sample1.MT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/fe/9d2f24508e4c6dc8bfba488ae5585e/no_host/sample1.MT/sample1.MT_R2.fastq.gz, /scratch/schudoma/imp3_test/work/bc/8c33cc37b6b62f57dea67ba440ae08/merged/sample1.MT.singles_R1.fastq.gz]
	// ]

	rnaspades(metaT_preprocessed_ch, "initial")

	// [
	// 	[sample_id:sample1.MT, library:paired, library_type:metaT], 
	// 	[
	// 		[/scratch/schudoma/imp3_test/work/fe/9d2f24508e4c6dc8bfba488ae5585e/no_host/sample1.MT/sample1.MT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/fe/9d2f24508e4c6dc8bfba488ae5585e/no_host/sample1.MT/sample1.MT_R2.fastq.gz], [/scratch/schudoma/imp3_test/work/bc/8c33cc37b6b62f57dea67ba440ae08/merged/sample1.MT.singles_R1.fastq.gz]]]

	// gq_input_ch = nevermore_main.out.fastqs
	// 	.map { sample, fastqs ->
	// 		sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
	// 			return tuple(sample_id, [fastqs].flatten())
	// 		}
	// 		.groupTuple()
	// 		.map { sample_id, fastqs -> return tuple(sample_id, fastqs.flatten()) }
	// 		gq_input_ch.view()


	// // this is for later
	// joined = metaT_ch.concat(metaG_ch).groupTuple(by:0, size: 2, remainder: true)
	// joined.view()

}
