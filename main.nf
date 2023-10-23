#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
// include { fastq_input } from "./nevermore/workflows/input"
include { metaT_input; metaG_input } from "./imp/workflows/input"

include { spades } from "./imp/modules/assemblers/spades"
include { bwa_index } from "./imp/modules/alignment/indexing/bwa_index"
include { extract_unmapped } from "./imp/modules/alignment/extract"

include { metaT_assembly } from "./imp/workflows/meta_t"

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

}


//  rule metaspades_hybrid_assembly_1:
//         input:
//             'Preprocessing/mg.r1.preprocessed.fq',
//             'Preprocessing/mg.r2.preprocessed.fq',
//             'Preprocessing/mg.se.preprocessed.fq',
//             'Preprocessing/mt.r1.preprocessed.fq',
//             'Preprocessing/mt.r2.preprocessed.fq',
//             'Preprocessing/mt.se.preprocessed.fq',
//             'Assembly/intermediary/mt.metaspades_preprocessed.1.final.contigs.fa'
//         output:
//             'Assembly/intermediary/mgmt.metaspades_hybrid.1/contigs.fa',
//             'Assembly/intermediary/mgmt.metaspades_hybrid.1.fa',
//             directory('Assembly/intermediary/mgmt.metaspades_hybrid.1')
//         params:
//             contigs = "--trusted-contigs Assembly/intermediary/mt.metaspades_preprocessed.1.final.contigs.fa"
//         resources:
//             runtime = "120:00:00",
//             mem = BIGMEMCORE
//         threads: getThreads(BIGCORENO)
//         conda: ENVDIR + "/IMP_assembly.yaml"
//         log: "logs/assembly_metaspades_hybrid_assembly_1.log"
//         message: "metaspades_hybrid_assembly_1: Performing hybrid assembly 1 from preprocessed reads using METASPADES"
//         shell:
//             """
//             if [ -d "{output[2]}" ]; then
//                 rm -rf {output[2]}
//             fi
//             METASPADES_ASSEMBLY_SHELL
//         """