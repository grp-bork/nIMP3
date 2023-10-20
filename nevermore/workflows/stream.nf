#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { classify_sample_with_library_info } from "../modules/functions"


def asset_dir = "${projectDir}/nevermore/assets"


workflow gffquant_stream {

	take:
		fastq_ch
	
	main:
		fastq_ch = fastq_ch
			.map { sample, fastqs ->
				sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				// return tuple(sample_id, sample.library, fastqs)
				return tuple(sample_id, fastqs)
			}
			//.groupTuple(by: [0, 1], size: 3, remainder: true, sort: true)
			.groupTuple(by: 0, size: 4, remainder: true, sort: true)
			

	emit:
		fastq_ch 


}