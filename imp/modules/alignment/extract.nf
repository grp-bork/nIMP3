process extract_unmapped {
	label "align"

	input:
	tuple val(sample), path(fastqs), path(index)
	val(stage)

	output:
	tuple val(sample), path("unmapped/${stage}/${sample.library_type}/${sample.id}/*.fastq.gz"), emit: fastqs, optional: true

	script:

	println "SAMPLE"
	println sample
	def reads1 = ""
	def reads2 = ""
	def filter_cmd_base = "samtools view --threads {task.cpus} -u"
	def filter_cmd =  ""

	def outpath = "unmapped/${stage}/${sample.library_type}/${sample.id}"
	def extract_cmd = ""

	// def check_cmd = "if [[ -z \"\$(gzip -dc ${outpath}/${sample.id}_R1.fastq.gz | head -n 1)\" ]]; then rm -f ${outpath}/${sample.id}_R1.fastq.gz; fi"

	if (sample.is_paired == true) {
		print "IN IS_PAIRED_BLOCK"
		reads1 += "${sample.id}_R1.fastq.gz"
		reads2 += "${sample.id}_R2.fastq.gz"
		filter_cmd += "${filter_cmd_base} -f4 -F 264 alignment.bam > self_unmapped.bam 2>> error.log\n"
		filter_cmd += "${filter_cmd_base} -f8 -F 260 alignment.bam > mate_unmapped.bam 2>> error.log\n"
		filter_cmd += "${filter_cmd_base} -f12 -F 256 alignment.bam > both_unmapped.bam 2>> error.log\n"

		filter_cmd += "samtools merge --threads ${task.cpus} -u - *_unmapped.bam 2>> error.log"
		filter_cmd += " | samtools collate -@ ${task.cpus} -o unmapped.bam - 2>> error.log"

		extract_cmd += "samtools fastq -1 ${outpath}/${sample.id}_R1.fastq.gz -2 ${outpath}/${sample.id}_R2.fastq.gz unmapped.bam"
		// check_cmd = "if [[ -z \"\$(gzip -dc ${outpath}/${sample.id}_R1.fastq.gz | head -n 1)\" ]]; then rm -f ${outpath}/*.fastq.gz; fi"
	} else {
		print "IN IS_SINGLE_BLOCK"
		reads1 += "${fastqs[0]}"
		filter_cmd += "${filter_cmd_base} -f 4 alignment.bam > unmapped.bam 2>> error.log"
		extract_cmd += "samtools fastq -0 ${outpath}/${sample.id}_R1.fastq.gz unmapped.bam"
	}



	"""
	mkdir -p unmapped/${stage}/${sample.library_type}/${sample.id}/
	bwa mem -v 1 -t ${task.cpus} ${sample.index_id} ${reads1} ${reads2} 2>> error.log | \
	samtools view --threads {task.cpus} -bS - > alignment.bam 2>> error.log

	${filter_cmd}
	${extract_cmd}
	if [[ -z "\$(gzip -dc ${outpath}/${sample.id}_R1.fastq.gz | head -n 1)" ]]; then rm -vf ${outpath}/*.fastq.gz; fi

	"""


	

}

// -f 4 -F 264: unmapped + (not mate unmapped and not secondary alignment)
// -f 8 -F 260: mate unmapped + (not read unmapped and not secondary alignment)
// -f 12 -F 256: both unmapped + not secondary alignment

// -F 0x800: not supplementary alignment

// SAMTOOLS_MEM = str(round(float(BIGMEMCORE[:-1]) * 0.75 - 0.5)) + "G"

// EXTRACT_UNMAPPED_SHELL = """
//  TMP_FILE=$(mktemp --tmpdir={TMPDIR} -t "alignment_XXXXXX.bam")
//  bwa mem -v 1 -t {threads} {input[3]} {input[0]} {input[1]} 2>> {log} | \
//   samtools view --threads {threads} -bS - > $TMP_FILE 2> {log}
//  samtools merge --threads {threads} -u - \
//   <(samtools view --threads {threads} -u  -f 4 -F 264 $TMP_FILE 2>> {log}) \
//   <(samtools view --threads {threads} -u -f 8 -F 260 $TMP_FILE 2>> {log}) \
//   <(samtools view --threads {threads} -u -f 12 -F 256 $TMP_FILE 2>> {log}) 2>> {log}| \
//   samtools view --threads {threads} -bF 0x800 - 2>> {log} | \
//   samtools sort --threads {threads} -m {SAMTOOLS_MEM} -n - 2>> {log} | \
//   bamToFastq -i stdin -fq {output[0]} -fq2 {output[1]} >> {log} 2>&1

//  if [[ -s {input[2]} ]]
//  then
//     bwa mem -v 1 -t {threads} {input[3]} {input[2]} 2>> {log}| \
//      samtools view --threads {threads} -bS - 2>> {log}| \
//      samtools view --threads {threads} -uf 4 - 2>> {log}| \
//      bamToFastq -i stdin -fq {output[2]} >> {log} 2>&1
//  else
//     echo "There are no singletons reads. {input[2]} is empty. Generating empty {output[2]}" >> {log} 
//     touch {output[2]}
//  fi

//  rm -rf $TMP_FILE
// """

// rule extract_unmapped:
//     input:
//         'Preprocessing/{type}.r1.preprocessed.fq',
//         'Preprocessing/{type}.r2.preprocessed.fq',
//         'Preprocessing/{type}.se.preprocessed.fq',
//         'Assembly/intermediary/{type}.%s_preprocessed.1.fa' % IMP_ASSEMBLER,
//         'Assembly/intermediary/{type}.%s_preprocessed.1.fa.amb' % IMP_ASSEMBLER,
//         'Assembly/intermediary/{type}.%s_preprocessed.1.fa.bwt' % IMP_ASSEMBLER,
//         'Assembly/intermediary/{type}.%s_preprocessed.1.fa.pac' % IMP_ASSEMBLER,
//         'Assembly/intermediary/{type}.%s_preprocessed.1.fa.sa' % IMP_ASSEMBLER,
//         'Assembly/intermediary/{type}.%s_preprocessed.1.fa.ann' % IMP_ASSEMBLER
//     output:
//         'Assembly/intermediary/{type}.r1.unmapped.fq',
//         'Assembly/intermediary/{type}.r2.unmapped.fq',
//         'Assembly/intermediary/{type}.se.unmapped.fq'
//     resources:
//         runtime = "12:00:00",
//         mem = BIGMEMCORE
//     threads: getThreads(BIGCORENO)
//     conda: ENVDIR + "/IMP_mapping.yaml"
//     log: "logs/assembly_extract_unmapped.{type}.log"
//     message: "extract_unmapped: Extracting unmapped {wildcards.type} reads from megahit assembly."
//     shell:
//         EXTRACT_UNMAPPED_SHELL
