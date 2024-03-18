process stream_gffquant_genome {
	label "gffquant"
	tag "gffquant.${sample}"

	input:
		tuple val(sample), path(fastqs), path(assembly_ref), path(gff), path(emapper_annotation)
	output:
		tuple val(sample), path("profiles/${sample}/*.txt.gz"), emit: results
		tuple val(sample), path("logs/${sample}.log")
		tuple val(sample), path("profiles/${sample}/*.pd.txt"), emit: profiles
		tuple val(sample), path("profiles/${sample}/*.coverage.txt"), emit: coverage_profiles

	script:
			def gq_output = "-o profiles/${sample}/${sample}"

			def gq_params = "-m small_genome --ambig_mode ${params.gq_ambig_mode} --db ${gff},${emapper_annotation}"
			
			gq_params += (params.gq_min_seqlen) ? (" --min_seqlen " + params.gq_min_seqlen) : ""
			gq_params += (params.gq_min_identity) ? (" --min_identity " + params.gq_min_identity) : ""
			gq_params += (params.gq_restrict_metrics) ? " --restrict_metrics ${params.gq_restrict_metrics}" : ""
			gq_params += (params.gq_keep_alignments) ? " --keep_alignment_file ${sample}.sam" : ""
			gq_params += (params.gq_unmarked_orphans) ? " --unmarked_orphans" : ""
			gq_params += (params.gq_with_coverage) ? " --with_coverage" : ""


			gq_params += " -t ${task.cpus}"

			def input_files = ""
			// we cannot auto-detect SE vs. PE-orphan!
			if (params.gq_single_end_library) {
				//input_files += "--singles \$(find . -maxdepth 1 -type l -name '*_R1.fastq.gz')"	
				input_files += "--fastq-singles ${fastqs}"
			} else {
				r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
				r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
				orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

				if (r1_files.size() != 0) {
					input_files += "--fastq-r1 ${r1_files.join(' ')}"
				}
				if (r2_files.size() != 0) {
					input_files += " --fastq-r2 ${r2_files.join(' ')}"
				}
				if (orphans.size() != 0) {
					input_files += " --fastq-orphans ${orphans.join(' ')}"
				}

				// input_files += "--fastq-r1 \$(find . -maxdepth 1 -type l -name '*_R1.fastq.gz' | grep -v singles)"
				// input_files += " --fastq-r2 \$(find . -maxdepth 1 -type l -name '*_R2.fastq.gz')"
				// input_files += " --fastq-orphans \$(find . -maxdepth 1 -type l -name '*singles*.fastq.gz')"
			}
	
			def gq_cmd = "gffquant ${gq_output} ${gq_params} --reference ${sample}.mmi --aligner ${params.gq_aligner} ${input_files}"

			"""
			set -e -o pipefail
			mkdir -p logs/ tmp/ profiles/

			minimap2 -x sr -d ${sample}.mmi ${assembly_ref}

			${gq_cmd} &> logs/${sample}.log
			"""

}
