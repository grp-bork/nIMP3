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
			meta.id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
			meta.library = sample.library
			meta.library_type = sample.library_type
			
			return tuple(meta, [fastqs].flatten())
		}
		.groupTuple(by: 0, size: 2, remainder: true)
		.map { sample, fastqs -> return tuple(sample, fastqs.flatten())}
	
	// metaT_preprocessed_ch.view()

	// [
	// 	[sample_id:sample1.MT, library:paired, library_type:metaT],
	// 	[/scratch/schudoma/imp3_test/work/fe/9d2f24508e4c6dc8bfba488ae5585e/no_host/sample1.MT/sample1.MT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/fe/9d2f24508e4c6dc8bfba488ae5585e/no_host/sample1.MT/sample1.MT_R2.fastq.gz, /scratch/schudoma/imp3_test/work/bc/8c33cc37b6b62f57dea67ba440ae08/merged/sample1.MT.singles_R1.fastq.gz]
	// ]


// [sample1.metaT, metaT, [id:sample1.metaT, is_paired:true, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/5a/d315af98dd81607d915fd6ced7e74a/no_host/sample1.metaT/sample1.metaT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/5a/d315af98dd81607d915fd6ced7e74a/no_host/sample1.metaT/sample1.metaT_R2.fastq.gz], [/scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.amb, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.ann, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.bwt, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.pac, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.sa]]
// [sample1.metaT, metaT, [id:sample1.metaT.singles, is_paired:false, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/b5/f5162a9919436a454a91a11ee65dd7/merged/sample1.metaT.singles_R1.fastq.gz], null]


// [sample1.metaT, metaT, [id:sample1.metaT, is_paired:true, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/5a/d315af98dd81607d915fd6ced7e74a/no_host/sample1.metaT/sample1.metaT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/5a/d315af98dd81607d915fd6ced7e74a/no_host/sample1.metaT/sample1.metaT_R2.fastq.gz]]
// [sample1.metaT, metaT, [id:sample1.metaT.singles, is_paired:false, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/b5/f5162a9919436a454a91a11ee65dd7/merged/sample1.metaT.singles_R1.fastq.gz]]


// [sample1.metaT, metaT, [/scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.amb, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.ann, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.bwt, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.pac, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.sa]]
// [sample1.metaT, metaT, [id:sample1.metaT, is_paired:true, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/5a/d315af98dd81607d915fd6ced7e74a/no_host/sample1.metaT/sample1.metaT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/5a/d315af98dd81607d915fd6ced7e74a/no_host/sample1.metaT/sample1.metaT_R2.fastq.gz], [/scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.amb, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.ann, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.bwt, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.pac, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.sa]]

// [sample1.metaT, metaT, [id:sample1.metaT.singles, is_paired:false, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/b5/f5162a9919436a454a91a11ee65dd7/merged/sample1.metaT.singles_R1.fastq.gz]]
// [sample1.metaT, metaT, [/scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.amb, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.ann, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.bwt, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.pac, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.sa]]
// [sample1.metaT, metaT, [id:sample1.metaT, is_paired:true, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/5a/d315af98dd81607d915fd6ced7e74a/no_host/sample1.metaT/sample1.metaT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/5a/d315af98dd81607d915fd6ced7e74a/no_host/sample1.metaT/sample1.metaT_R2.fastq.gz], [/scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.amb, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.ann, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.bwt, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.pac, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.sa]]
// [sample1.metaT, metaT, [id:sample1.metaT.singles, is_paired:false, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/b5/f5162a9919436a454a91a11ee65dd7/merged/sample1.metaT.singles_R1.fastq.gz], null]





	rnaspades(metaT_preprocessed_ch, "initial")
	bwa_index(rnaspades.out.contigs, "initial")

	metaT_post_assembly_check_ch = nevermore_main.out.fastqs
		.filter { it[0].library_type == "metaT" }
		.map { sample, fastqs -> 
			sample_base_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
			return tuple(sample_base_id, sample.library_type, sample, [fastqs].flatten())
		}

	// metaT_post_assembly_check_ch.view()

	// metaT_post_assembly_check_ch = metaT_post_assembly_check_ch	
		.combine(bwa_index.out.index, by: [0, 1])
		.map { sample_id, libtype, sample, fastqs, index ->
			sample.index_id = sample_id
			return tuple(sample, fastqs, index) }
	
	// bwa_index.out.index.view()
	// metaT_post_assembly_check_ch.view()

	extract_unmapped(metaT_post_assembly_check_ch, "initial")
	extract_unmapped.out.fastqs.view()

	// [sample1.metaT, metaT, [id:sample1.metaT, is_paired:true, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/5a/d315af98dd81607d915fd6ced7e74a/no_host/sample1.metaT/sample1.metaT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/5a/d315af98dd81607d915fd6ced7e74a/no_host/sample1.metaT/sample1.metaT_R2.fastq.gz], [/scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.amb, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.ann, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.bwt, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.pac, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.sa]]
	// [sample1.metaT, metaT, [id:sample1.metaT.singles, is_paired:false, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/b5/f5162a9919436a454a91a11ee65dd7/merged/sample1.metaT.singles_R1.fastq.gz], [/scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.amb, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.ann, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.bwt, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.pac, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.sa]]

	// [sample1.metaT, metaT, [id:sample1.metaT, is_paired:true, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/5a/d315af98dd81607d915fd6ced7e74a/no_host/sample1.metaT/sample1.metaT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/5a/d315af98dd81607d915fd6ced7e74a/no_host/sample1.metaT/sample1.metaT_R2.fastq.gz], [/scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.amb, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.ann, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.bwt, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.pac, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.sa]]
	// [sample1.metaT, metaT, [id:sample1.metaT.singles, is_paired:false, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/b5/f5162a9919436a454a91a11ee65dd7/merged/sample1.metaT.singles_R1.fastq.gz], [/scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.amb, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.ann, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.bwt, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.pac, /scratch/schudoma/imp3_test/work/9e/7aeb0b74556a625e5831d7f7ab44ee/index/initial/metaT/sample1.metaT/sample1.metaT.sa]]

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
