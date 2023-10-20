params.stranded = null
params.kmer_steps = "25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89,93,97"


process rnaspades {
	label "spades"

	input:
	tuple val(sample), path(fastqs)
	val(stage)

	output:
	tuple val(sample), path("assemblies/rnaspades/${stage}/${sample.id}/${sample.id}.transcripts.fasta")

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


	"""
	rnaspades.py --meta -t ${task.cpus} -m ${mem_gb} -o assemblies/rnaspades/${stage}/${sample.id} ${stranded} ${kmers} ${input_files}
	"""
}

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
