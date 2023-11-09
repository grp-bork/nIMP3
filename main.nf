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
include { hybrid_megahit } from "./imp/modules/assemblers/megahit"
include { get_unmapped_reads } from "./imp/workflows/extract"
include { concatenate_contigs } from "./imp/modules/assemblers/functions"


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

params.assembler = "megahit"


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

	// metaT_assembly.out.reads.view()
	// [[id:sample1.metaT, library:paired, library_type:metaT], [/scratch/schudoma/imp3_test/work/0a/c32c5be6621089f9c71d87fc2fd308/no_host/sample1.metaT/sample1.metaT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/0a/c32c5be6621089f9c71d87fc2fd308/no_host/sample1.metaT/sample1.metaT_R2.fastq.gz, /scratch/schudoma/imp3_test/work/f3/2d793cfd78eda426ffc94ce4f5712a/merged/sample1.metaT.singles_R1.fastq.gz]]
	// metaT_assembly.out.final_contigs.view()
	// [[id:sample1.metaT, library:paired, library_type:metaT], /scratch/schudoma/imp3_test/work/29/885cb85b918cada4ee1e07111a2434/assemblies/rnaspades/final/metaT/sample1.metaT/sample1.metaT.final_contigs.fasta]

	assembly_prep(
		nevermore_main.out.fastqs
			.filter { it[0].library_type == "metaG" }
	)

	metaG_assembly_ch = assembly_prep.out.reads
	// metaG_assembly_ch.view()

	hybrid_assembly_input_ch = metaT_assembly.out.reads
		.map { sample, fastqs ->
			def meta = [:]
			meta.id = sample.id.replaceAll(/\.metaT/, "") 
			// return tuple(meta.id, meta, fastqs)
			
			return tuple(meta, fastqs)
			
		}
		.concat(
			metaG_assembly_ch
				.map { sample, fastqs ->
					def meta = [:]
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
					def meta = [:]
					meta.id = sample.id.replaceAll(/\.metaT/, "")
					return tuple(meta, contigs)
				},			
		)
		.groupTuple()
		.map { sample, data -> return tuple(sample, data[0], data[1]) }

	hybrid_assembly_input_ch.dump(pretty: true, tag: "hybrid_assembly_input_ch")


	if (params.assembler == "spades") {
		metaspades(hybrid_assembly_input_ch, "initial")
		contigs_ch = metaspades.out.contigs		
	} else {
		hybrid_megahit(hybrid_assembly_input_ch, "initial")
		contigs_ch = hybrid_megahit.out.contigs
	}
	// contigs_ch.view()

	contigs_ch = contigs_ch.map {
		sample, fastqs -> 
		sample.library_type = "hybrid"
		return tuple(sample, fastqs)
	}

	bwa_index(contigs_ch, "initial")

	bwa_index.out.index.dump(pretty: true, tag: "bwa_index.out.index")
	nevermore_main.out.fastqs.dump(pretty: true, tag: "nevermore_main.out.fastqs")

	/*[DUMP: hybrid_assembly_input_ch] [
    {
        "id": "sample1"
    },
    [
        "/scratch/schudoma/imp3_test/work/9f/7654e40f9169728a3f845b272edcbc/no_host/sample1.metaT/sample1.metaT_R1.fastq.gz",
        "/scratch/schudoma/imp3_test/work/9f/7654e40f9169728a3f845b272edcbc/no_host/sample1.metaT/sample1.metaT_R2.fastq.gz",
        "/scratch/schudoma/imp3_test/work/ee/b9a3d3515e91e67ad2b0d28a25a827/merged/sample1.metaT.singles_R1.fastq.gz",
        "/scratch/schudoma/imp3_test/work/7e/72d69ed3806911f6d562adc2f11bab/no_host/sample1.metaG/sample1.metaG_R1.fastq.gz",
        "/scratch/schudoma/imp3_test/work/7e/72d69ed3806911f6d562adc2f11bab/no_host/sample1.metaG/sample1.metaG_R2.fastq.gz",
        "/scratch/schudoma/imp3_test/work/3b/764b7ee1663c14ed712dc019176639/merged/sample1.metaG.singles_R1.fastq.gz"
    ],
    "/scratch/schudoma/imp3_test/work/3c/7cf332ec2b49e3bfdd2bfed1f2d342/assemblies/metaT_megahit/final/metaT/sample1.metaT/sample1.metaT.final_contigs.fasta"
]
	*/
	combined_assembly_input_index_ch = hybrid_assembly_input_ch
		.map { sample, fastqs, contigs -> return tuple(sample.id, sample, fastqs) }
		.join(bwa_index.out.index, by: 0)
		.map { sample_id, sample, fastqs, libtype, index -> return tuple(sample_id, sample, fastqs, index) }
	combined_assembly_input_index_ch.dump(pretty: true, tag: "combined_assembly_input_index_ch")

	// metaT_paired_unmapped_ch = hybrid_assembly_input_ch
	metaT_paired_unmapped_ch = combined_assembly_input_index_ch
		// .map { sample, fastqs, contigs ->
		.map { sample_id, sample, fastqs, index ->
			def new_sample = [:]
			new_sample.id = sample.id + ".metaT"
			new_sample.library_type = "metaT"
			new_sample.is_paired = true
			new_sample.index_id = sample_id
			def wanted_fastqs = []
			wanted_fastqs.addAll(fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)singles(.*)") && it.name.matches("(.*)metaT(.*)") } ))
			wanted_fastqs.addAll(fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") && it.name.matches("(.*)metaT(.*)") } ))
			return tuple(new_sample, wanted_fastqs, index)
		}
	metaT_single_unmapped_ch = combined_assembly_input_index_ch
		.map { sample_id, sample, fastqs, index ->
			def new_sample = [:]
			new_sample.id = sample.id + ".metaT.singles"
			new_sample.library_type = "metaT"
			new_sample.is_paired = false
			new_sample.index_id = sample_id
			def wanted_fastqs = []
			wanted_fastqs.addAll(fastqs.findAll( { it.name.matches("(.*)singles(.*)") && it.name.matches("(.*)metaT(.*)") } ))
			return tuple(new_sample, wanted_fastqs, index)
		}
	metaG_paired_unmapped_ch = combined_assembly_input_index_ch
		.map { sample_id, sample, fastqs, index ->
			def new_sample = [:]
			new_sample.id = sample.id + ".metaG"
			new_sample.library_type = "metaG"
			new_sample.is_paired = true
			new_sample.index_id = sample_id
			def wanted_fastqs = []
			wanted_fastqs.addAll(fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)singles(.*)") && it.name.matches("(.*)metaG(.*)") } ))
			wanted_fastqs.addAll(fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") && it.name.matches("(.*)metaG(.*)") } ))
			return tuple(new_sample, wanted_fastqs, index)
		}
	metaG_single_unmapped_ch = combined_assembly_input_index_ch
		.map { sample_id, sample, fastqs, index ->
			def new_sample = [:]
			new_sample.id = sample.id + ".metaG.singles"
			new_sample.library_type = "metaG"
			new_sample. = false
			new_sample.index_id = sample_id
			def wanted_fastqs = []
			wanted_fastqs.addAll(fastqs.findAll( { it.name.matches("(.*)singles(.*)") && it.name.matches("(.*)metaG(.*)") } ))
			return tuple(new_sample, wanted_fastqs, index)
		}
	extract_unmapped_ch = Channel.empty()
		.concat(metaT_paired_unmapped_ch)
		.concat(metaT_single_unmapped_ch)
		.concat(metaG_paired_unmapped_ch)
		.concat(metaG_single_unmapped_ch)
	
	extract_unmapped_ch.dump(pretty: true, tag: "extract_unmapped_ch")


	// [sample1.metaT, [id:sample1.metaT, is_paired:true, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/32/89e706473152c55350201e65cd16a3/no_host/sample1.metaT/sample1.metaT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/32/89e706473152c55350201e65cd16a3/no_host/sample1.metaT/sample1.metaT_R2.fastq.gz]]
	// [[id:sample1.metaT, is_paired:true, library:paired, library_type:metaT, merged:true], [/scratch/schudoma/imp3_test/work/32/89e706473152c55350201e65cd16a3/no_host/sample1.metaT/sample1.metaT_R1.fastq.gz, /scratch/schudoma/imp3_test/work/32/89e706473152c55350201e65cd16a3/no_host/sample1.metaT/sample1.metaT_R2.fastq.gz]]
	// [sample1.metaG, [id:sample1.metaG, is_paired:true, library:paired, library_type:metaG, merged:true], [/scratch/schudoma/imp3_test/work/da/eea24ab29b720733d771235b1d7a15/no_host/sample1.metaG/sample1.metaG_R1.fastq.gz, /scratch/schudoma/imp3_test/work/da/eea24ab29b720733d771235b1d7a15/no_host/sample1.metaG/sample1.metaG_R2.fastq.gz]]
	// [[id:sample1.metaG, is_paired:true, library:paired, library_type:metaG, merged:true], [/scratch/schudoma/imp3_test/work/da/eea24ab29b720733d771235b1d7a15/no_host/sample1.metaG/sample1.metaG_R1.fastq.gz, /scratch/schudoma/imp3_test/work/da/eea24ab29b720733d771235b1d7a15/no_host/sample1.metaG/sample1.metaG_R2.fastq.gz]]
	// [sample1, null, [/scratch/schudoma/imp3_test/work/52/0bf03b27daa9dbeb9a72a2936ef999/index/initial/null/sample1/sample1.amb, /scratch/schudoma/imp3_test/work/52/0bf03b27daa9dbeb9a72a2936ef999/index/initial/null/sample1/sample1.ann, /scratch/schudoma/imp3_test/work/52/0bf03b27daa9dbeb9a72a2936ef999/index/initial/null/sample1/sample1.bwt, /scratch/schudoma/imp3_test/work/52/0bf03b27daa9dbeb9a72a2936ef999/index/initial/null/sample1/sample1.pac, /scratch/schudoma/imp3_test/work/52/0bf03b27daa9dbeb9a72a2936ef999/index/initial/null/sample1/sample1.sa]]

	// get_unmapped_reads(nevermore_main.out.fastqs, bwa_index.out.index)

	base_id_ch = nevermore_main.out.fastqs
		.map { sample, fastqs -> 
			def sample_base_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "").replaceAll(/.meta[GT]$/, "")
			return tuple(sample_base_id, sample, [fastqs].flatten())
		}
	
	base_id_ch.dump(pretty: true, tag: "base_id_ch")

	index_and_fastqs_ch = bwa_index.out.index.combine(base_id_ch, by: 0)
	index_and_fastqs_ch.dump(pretty: true, tag: "index_and_fastqs_ch")

	with_index_ch = base_id_ch.combine(bwa_index.out.index)
	with_index_ch.dump(pretty: true, tag: "with_index_ch")

	joined_ch = base_id_ch.join(bwa_index.out.index, by: 0)
	joined_ch.dump(pretty: true, tag: "joined_ch")

	// base_id_ch.combine(bwa_index.out.index, by: 0).dump(pretty: true, tag: "base_id_ch")

	// post_assembly_check_ch = nevermore_main.out.fastqs
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
	extract_unmapped.out.fastqs.dump(pretty: true, tag: "extract_unmapped_ch")
	
	unmapped_ch = extract_unmapped.out.fastqs
		.map { sample, fastqs ->
			sample.id = sample.index_id
			return tuple(sample.id, sample, fastqs)
		}
		.groupTuple(by: 0, size: 2, remainder: true)
		.map { sample_id, sample, fastqs -> 
			def meta = [:]
			meta.id = sample_id
			// meta.library_type = sample.library_type
			return tuple(meta, fastqs.flatten())
		}
		.groupTuple(by: 0, size: 2, remainder: true) //, sort: true)
		.map { sample, fastqs ->
			// sample.library_type = sample.library_type[0]
			sample.library_type = "hybrid"
			return tuple(sample, fastqs.flatten())
		}

	unmapped_ch.dump(pretty: true, tag: "unmapped_ch")
	
	empty_file = file("${launchDir}/NO_INPUT")
	empty_file.text = "NOTHING TO SEE HERE."
	print empty_file


	unmapped_contigs_ch = Channel.empty()
	if (params.assembler == "megahit") {
		megahit_hybrid_unmapped(unmapped_ch, Channel.of(empty_file))
		unmapped_contigs_ch = megahit_hybrid_unmapped.out.contigs
	}
	unmapped_contigs_ch.dump(pretty: true, tag: "unmapped_contigs_ch")

	all_contigs_ch = contigs_ch
		.concat(unmapped_contigs_ch)
		.groupTuple(by: 0, size: 2, remainder: true, sort: true)

	all_contigs_ch.dump(pretty: true, tag: "all_contigs_ch")
	
	concatenate_contigs(all_contigs_ch, "final", params.assembler)
	// hybrid_megahit(unmapped_ch.combine(Channel.of(empty_file)))
	// final_assembly_ch = get_unmapped_reads.out.reads
	// 	.map { sample, fastqs -> return tuple(sample, fastqs, [empty_file])}
	// final_assembly_ch.view()
	
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