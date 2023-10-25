params.mink = 25
params.maxk = 100
params.stepk = 4

process metaT_megahit {
	label "megahit"

	input:
	tuple val(sample), path(fastqs)
	val(stage)

	output:
	tuple val(sample), path("assemblies/metaT_megahit/${stage}/${sample.library_type}/${sample.id}/${sample.id}.${stage}.*.fasta"), emit: contigs

	script:
	def mem_gb = task.memory.toGiga()
	def mem = task.memory.toBytes()

	def input_files = ""
	// we cannot auto-detect SE vs. PE-orphan!
	r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)singles(.*)") } )
	r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	orphans = fastqs.findAll( { it.name.matches("(.*)singles(.*)") } )

	if (r1_files.size() != 0) {
		input_files += "-1 ${r1_files.join(' ')}"
	}
	if (r2_files.size() != 0) {
		input_files += " -2 ${r2_files.join(' ')}"
	}
	if (orphans.size() != 0) {
		input_files += " -r ${orphans.join(' ')}"
	}

	def kmer_params = "--k-min ${params.mink} --k-max ${params.maxk} --k-step ${params.stepk}"
	def outdir = "assemblies/metaT_megahit/${stage}/${sample.library_type}/${sample.id}"
	
	"""
	mkdir -p ${outdir}/
	megahit -t ${task.cpus} --cpu-only -m ${mem} ${input_files} ${kmer_params} --bubble-level 0 --mem-flag 1
	cp -v megahit_out/final.contigs.fa ${outdir}/${sample.id}.${stage}.transcripts.fasta
	"""
	// rnaspades.py -t ${task.cpus} -m ${mem_gb} -o assemblies/rnaspades/${stage}/${sample.library_type}/${sample.id} ${stranded} ${kmers} ${input_files}
	// mv -v transcripts.fasta assemblies/spades/${stage}/${sample.library_type}/${sample.id}/${sample.id}.${stage}.transcripts.fasta 
	// rnaspades.py --meta -t ${task.cpus} -m ${mem_gb} -o assemblies/rnaspades/${stage}/${sample.id} ${stranded} ${kmers} ${input_files}
}

// MEGAHIT_ASSEMBLY_SHELL = """
// if [ -d "{output[0]}" ]; then
//     rm -rf {output[0]}
// fi
// MAX_MEM="$(({BIGMEMTOTAL} * 1000000000))"
// megahit -1 {input[0]} \
//  -2 {input[1]} \
//  -r {input[2]} \
//  -o {output[0]} \
//  --k-min {config[assembly][mink]} \
//  --k-max {config[assembly][maxk]} \
//  --k-step {config[assembly][step]} \
//  --bubble-level 0 \
//  -t {threads} --cpu-only \
//  -m "${{MAX_MEM}}" \
//  --mem-flag 1 > {log} 2>&1
// ln -fs $(echo {output[1]} | cut -f 3,4 -d /) {output[2]} && touch -h {output[2]}
// """

// rule megahit_assembly_from_preprocessing:
//     input:
//         'Preprocessing/{type}.r1.preprocessed.fq',
//         'Preprocessing/{type}.r2.preprocessed.fq',
//         'Preprocessing/{type}.se.preprocessed.fq'
//     output:
//         directory('Assembly/intermediary/{type}.megahit_preprocessed.1'),
//         'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa',
//         'Assembly/intermediary/{type}.megahit_preprocessed.1.fa'
//     resources:
//         runtime = "120:00:00",
//         mem = BIGMEMCORE
//     threads: getThreads(BIGCORENO)
//     conda: ENVDIR + "/IMP_assembly.yaml"
//     log: "logs/assembly_megahit_assembly_from_preprocessing.{type}.log"
//     message: "megahit_assembly_from_preprocessing: Performing {wildcards.type} assembly step 1 from preprocessed reads using MEGAHIT"
//     shell:
//         MEGAHIT_ASSEMBLY_SHELL

// rule megahit_assembly_from_unmapped:
//     input:
//         'Assembly/intermediary/{type}.r1.unmapped.fq',
//         'Assembly/intermediary/{type}.r2.unmapped.fq',
//         'Assembly/intermediary/{type}.se.unmapped.fq'
//     output:
//         directory('Assembly/intermediary/{type}.megahit_unmapped.2'),
//         'Assembly/intermediary/{type}.megahit_unmapped.2/final.contigs.fa',
//         'Assembly/intermediary/{type}.megahit_unmapped.2.fa'
//     resources:
//         runtime = "120:00:00",
//         mem = BIGMEMCORE
//     threads: getThreads(BIGCORENO)
//     conda: ENVDIR + "/IMP_assembly.yaml"
//     log: "logs/assembly_megahit_assembly_from_unmapped.{type}.log"
//     message: "megahit_assembly_from_unmapped: Performing {wildcards.type} assembly step 2 from unmapped reads using MEGAHIT"
//     shell:
//         MEGAHIT_ASSEMBLY_SHELL
