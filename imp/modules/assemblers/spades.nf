params.stranded = null
// params.kmer_steps = "25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89,93,97"
params.kmer_steps = "25 29 33 37 41 45 49 53 57 61 65 69 73 77 81 85 89 93 97"


process rnaspades {
	container "docker://quay.io/biocontainers/spades:3.12.0--h9ee0642_3"
	label "spades"

	input:
	tuple val(sample), path(fastqs)
	val(stage)

	output:
	tuple val(sample), path("assemblies/rnaspades/${stage}/${sample.library_source}/${sample.id}/${sample.id}.${stage}.*.fasta"), emit: contigs

	script:

	def stranded = (params.stranded) ? params.stranded : ""
	def kmers = "-k ${params.kmer_steps}"
	def mem_gb = task.memory.toGiga()

	def input_files = ""
	// we cannot auto-detect SE vs. PE-orphan!
	r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)singles(.*)") } )
	r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	orphans = fastqs.findAll( { it.name.matches("(.*)singles(.*)") } )

	if (r1_files.size() != 0) {
		input_files += "--pe1-1 ${r1_files.join(' ')}"
	}
	if (r2_files.size() != 0) {
		input_files += " --pe1-2 ${r2_files.join(' ')}"
	}
	if (orphans.size() != 0) {
		input_files += " --pe1-s ${orphans.join(' ')}"
	}

	// --meta and --rna seem mutually exclusive?

	"""
	rnaspades.py -t ${task.cpus} -m ${mem_gb} -o assemblies/rnaspades/${stage}/${sample.library_source}/${sample.id} ${stranded} ${kmers} ${input_files}
	mv -v assemblies/rnaspades/${stage}/${sample.library_source}/${sample.id}/transcripts.fasta assemblies/rnaspades/${stage}/${sample.library_source}/${sample.id}/${sample.id}.${stage}.transcripts.fasta 
	"""
	// mv -v transcripts.fasta assemblies/spades/${stage}/${sample.library_source}/${sample.id}/${sample.id}.${stage}.transcripts.fasta 
	// rnaspades.py --meta -t ${task.cpus} -m ${mem_gb} -o assemblies/rnaspades/${stage}/${sample.id} ${stranded} ${kmers} ${input_files}
}

// spades.py --meta \
//  --pe1-1 {input[0]} \
//  --pe1-2 {input[1]} \
//  --pe1-s {input[2]} \
//  --pe2-1 {input[3]} \
//  --pe2-2 {input[4]} \
//  --pe2-s {input[5]} \
//   {params.contigs} \
//   -t {threads} \
//   -m {BIGMEMTOTAL} \
//   -k {KMER_STEPS} \
//   {LONG_READ_ARG} \
//   -o {output[0]} > {log} 2>&1
// ln -fs  {output[1]} {output[2]}


process metaspades {
	container "docker://quay.io/biocontainers/spades:3.12.0--h9ee0642_3"
	label "spades"

	input:
	tuple val(sample), path(fastqs), path(contigs)
	val(stage)

	output:
	tuple val(sample), path("assemblies/metaspades/${stage}/${sample.id}/${sample.id}.${stage}.*.fasta"), emit: contigs

	script:

	def kmers = "-k ${params.kmer_steps}"
	def mem_gb = task.memory.toGiga()

	def input_files = ""
	
	// we cannot auto-detect SE vs. PE-orphan!
	r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)singles(.*)") } )
	r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	orphan_files = fastqs.findAll( { it.name.matches("(.*)singles(.*)") } )

	r1_mt_files = r1_files.findAll( { it.name.matches("(.*).metaT_R1.fastq.gz" ) })
	r2_mt_files = r2_files.findAll( { it.name.matches("(.*).metaT_R2.fastq.gz" ) })
	orphan_mt_files = orphan_files.findAll( { it.name.matches("(.*).metaT.singles(.*)")})

	r1_mg_files = r1_files.findAll( { it.name.matches("(.*).metaG_R1.fastq.gz" ) })
	r2_mg_files = r2_files.findAll( { it.name.matches("(.*).metaG_R2.fastq.gz" ) })
	orphan_mg_files = orphan_files.findAll( { it.name.matches("(.*).metaG.singles(.*)(.*)")})


	if (r1_mg_files.size() != 0) {
		input_files += "--pe1-1 ${r1_mg_files.join(' ')}"
	}
	if (r2_mg_files.size() != 0) {
		input_files += " --pe1-2 ${r2_mg_files.join(' ')}"
	}
	if (orphan_mg_files.size() != 0) {
		input_files += " --pe1-s ${orphan_mg_files.join(' ')}"
	}

	if (r1_mt_files.size() != 0) {
		input_files += " --pe2-1 ${r1_mt_files.join(' ')}"
	}
	if (r2_mt_files.size() != 0) {
		input_files += " --pe2-2 ${r2_mt_files.join(' ')}"
	}
	if (orphan_mt_files.size() != 0) {
		input_files += " --pe2-s ${orphan_mt_files.join(' ')}"
	}

	def contig_str = "--trusted-contigs ${contigs}"

	"""
	spades.py --meta -t ${task.cpus} -m ${mem_gb} -o assemblies/metaspades/${stage}/${sample.id} ${kmers} ${input_files} ${contig_str}
	mv -v assemblies/metaspades/${stage}/${sample.id}/contigs.fasta assemblies/metaspades/${stage}/${sample.id}/${sample.id}.${stage}.contigs.fasta 
	"""
	// mv -v transcripts.fasta assemblies/spades/${stage}/${sample.library_source}/${sample.id}/${sample.id}.${stage}.transcripts.fasta 
	// rnaspades.py --meta -t ${task.cpus} -m ${mem_gb} -o assemblies/rnaspades/${stage}/${sample.id} ${stranded} ${kmers} ${input_files}
}


"""


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


// METASPADES_ASSEMBLY_SHELL = """
// if [ -d "{output[0]}" ]; then
//     rm -rf {output[0]}
// fi
// rnaspades.py --meta \
//  --pe1-1 {input[0]} \
//  --pe1-2 {input[1]} \
//  --pe1-s {input[2]} \
//   -t {threads} \
//   -m {BIGMEMTOTAL} \
//   -k {KMER_STEPS} \
//   {STRANDED_ARG} \
//   -o {output[0]} > {log} 2>&1
// ln -fs  {output[1]} {output[2]}
// """

// kmersteps=range(config['assembly']['mink'],config['assembly']['maxk']+1,config['assembly']['step'])
// KMER_STEPS =",".join(map(str,kmersteps))
// STRANDED_ARG = ""
// if MT_STRANDED == 1:
//     STRANDED_ARG =  "-ss-rf"
// elif MT_STRANDED == 2:
//     STRANDED_ARG =  "-ss-fr" 

// rule metaspades_assembly_from_preprocessing:
//     input:
//         'Preprocessing/mt.r1.preprocessed.fq',
//         'Preprocessing/mt.r2.preprocessed.fq',
//         'Preprocessing/mt.se.preprocessed.fq'
//     output:
//         directory('Assembly/intermediary/mt.metaspades_preprocessed.1'),
//         'Assembly/intermediary/mt.metaspades_preprocessed.1/transcripts.fasta',
//         'Assembly/intermediary/mt.metaspades_preprocessed.1.fa'
//     resources:
//         runtime = "120:00:00",
//         mem = BIGMEMCORE
//     threads: getThreads(BIGCORENO)
//     conda: ENVDIR + "/IMP_assembly.yaml"
//     log: "logs/assembly_metaspades_assembly_from_preprocessing.mt.log"
//     message: "metaspades_assembly_from_preprocessing: Performing mt assembly step 1 from preprocessed reads using MetaSpades"
//     shell:
//         METASPADES_ASSEMBLY_SHELL
