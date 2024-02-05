#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"

include { minimap2_index } from "./nevermore/modules/align/index"
include { stream_gffquant_genome } from "./nevermore/modules/profilers/gffquant_genome"

include { collate_feature_counts } from "./nevermore/modules/profilers/gffquant"

params.gq_collate_columns = "uniq_scaled,combined_scaled"


if (params.input_dir && params.remote_input_dir) {
	log.info """
		Cannot process both --input_dir and --remote_input_dir. Please check input parameters.
	""".stripIndent()
	exit 1
} else if (!params.input_dir && !params.remote_input_dir) {
	log.info """
		Neither --input_dir nor --remote_input_dir set.
	""".stripIndent()
	exit 1
}

def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir

params.ignore_dirs = ""



workflow {

	fastq_input(
		Channel.fromPath(input_dir + "/*", type: "dir")
			.filter { !params.ignore_dirs.split(",").contains(it.name) },
		Channel.of(null)
	)

	fastq_ch = fastq_input.out.fastqs
	
	// /g/scb2/bork/data/MAGs/annotations/Larkin_2021_GO-SHIP_marine/psa_megahit/eggnog-mapper/SAMN15782773/SAMN15782773.emapper.annotations.gz
	emapper_ch = Channel.fromPath(params.annotation_dir + "/" + params.study + "/psa_megahit/eggnog-mapper/**.emapper.annotations.gz")
		.map { file ->
			return tuple(file.name.replaceAll(/\.emapper.annotations.gz$/, ""), file)
		}
	emapper_ch.dump(pretty: true, tag: "emapper_ch")
	// /g/scb2/bork/data/MAGs/annotations/Larkin_2021_GO-SHIP_marine/psa_megahit/prodigal/SAMN15782773.psa_megahit.prodigal.gff.gz
	gff_ch = Channel.fromPath(params.annotation_dir + "/" + params.study + "/psa_megahit/prodigal/**.psa_megahit.prodigal.gff.gz")
		.map { file ->
			return tuple(file.name.replaceAll(/\.psa_megahit\.prodigal\.gff\.gz$/, ""), file)
		}
	gff_ch.dump(pretty: true, tag: "gff_ch")
	// /g/scb2/bork/data/MAGs/Larkin_2021_GO-SHIP_marine/psa_megahit/assemblies/SAMN15782773-assembled.fa.gz
	assembly_ch = Channel.fromPath(params.assembly_dir + "/" + params.study + "/psa_megahit/assemblies/**-assembled.fa.gz")
		.map { file ->
			return tuple(file.name.replaceAll(/-assembled\.fa\.gz$/, ""), file)
		}
	assembly_ch.dump(pretty: true, tag: "assembly_ch")

	// minimap2_index(assembly_ch)

	gq_resources_ch = assembly_ch
		.join(gff_ch, by: 0)
		.join(emapper_ch, by: 0)
	gq_resources_ch.dump(pretty: true, tag: "gq_resources_ch")

	nevermore_main(fastq_ch)

	gq_input_ch = nevermore_main.out.fastqs
		.map { sample, fastqs ->
			sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
			return tuple(sample_id, [fastqs].flatten())
		}
		.groupTuple()
		.map { sample_id, fastqs -> return tuple(sample_id, [fastqs].flatten()) }
		.join(gq_resources_ch, by: 0)
	
	stream_gffquant_genome(gq_input_ch)
	feature_count_ch = stream_gffquant_genome.out.results
	
	feature_count_ch = feature_count_ch
		.map { sample, files -> return files }
		.flatten()
		.filter { !it.name.endsWith("Counter.txt.gz") }
		.filter { params.collate_gene_counts || !it.name.endsWith("gene_counts.txt.gz") }
		.map { file -> 
			def category = file.name
				.replaceAll(/\.txt\.gz$/, "")
				.replaceAll(/.+\./, "")
			return tuple(category, file)
		}
		.groupTuple(sort: true)
		.combine(
			Channel.from(params.gq_collate_columns.split(","))
		)

	collate_feature_counts(feature_count_ch)
	

}
