params.mink = 25
params.maxk = 100
params.stepk = 4
params.kmer_steps = "25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89,93,97"


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

	def kmer_params = "--k-list ${params.kmer_steps}" //--k-min ${params.mink} --k-max ${params.maxk} --k-step ${params.stepk}"
	def outdir = "assemblies/metaT_megahit/${stage}/${sample.library_type}/${sample.id}"
	
	"""
	mkdir -p ${outdir}/
	megahit -t ${task.cpus} --cpu-only -m ${mem} ${input_files} ${kmer_params} --bubble-level 0 --mem-flag 1
	cp -v megahit_out/final.contigs.fa ${outdir}/${sample.id}.${stage}.transcripts.fasta
	"""
	
}


process hybrid_megahit {
	label "megahit"

	input:
	tuple val(sample), path(fastqs), path(contigs)
	val(stage)

	output:
	tuple val(sample), path("assemblies/hybrid_megahit/${stage}/${sample.library_type}/${sample.id}/${sample.id}.${stage}.*.fasta"), emit: contigs

	script:
	def mem_gb = task.memory.toGiga()
	def mem = task.memory.toBytes()

	def input_files = ""
	def r1_files = []
	def r2_files = []
	def orphan_files = []

	// we cannot auto-detect SE vs. PE-orphan!
	r1_files.addAll(fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)singles(.*)") && it.name.matches("(.*)metaG(.*)") } ))
	r2_files.addAll(fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") && it.name.matches("(.*)metaG(.*)") } ))
	orphan_files.addAll(fastqs.findAll( { it.name.matches("(.*)singles(.*)") && it.name.matches("(.*)metaG(.*)") } ))

	r1_files.addAll(fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)singles(.*)") && it.name.matches("(.*)metaT(.*)") } ))
	r2_files.addAll(fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") && it.name.matches("(.*)metaT(.*)") } ))
	orphan_files.addAll(fastqs.findAll( { it.name.matches("(.*)singles(.*)") && it.name.matches("(.*)metaT(.*)") } ))

	orphan_files.addAll(contigs.findAll( { it.name != "NO_INPUT" }))

	if (r1_files.size() != 0) {
		input_files += "-1 ${r1_files.join(',')}"
	}
	if (r2_files.size() != 0) {
		input_files += " -2 ${r2_files.join(',')}"
	}	
	if (orphan_files.size() != 0) {
		input_files += " -r ${orphan_files.join(',')}"
	}	

	def kmer_params = "--k-list ${params.kmer_steps}"
	def outdir = "assemblies/hybrid_megahit/${stage}/${sample.library_type}/${sample.id}"
	
	"""
	mkdir -p ${outdir}/
	megahit -t ${task.cpus} --cpu-only -m ${mem} ${input_files} ${kmer_params} --bubble-level 0 --mem-flag 1
	cp -v megahit_out/final.contigs.fa ${outdir}/${sample.id}.${stage}.contigs.fasta
	"""
	
}

