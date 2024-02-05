

process minimap2_index {
	input:
		tuple val(genome_id), path(fasta)
	output:
		tuple val(genome_id), path("indexes/minimap2/${genome_id}/${genome_id}.mmi"), emit: index
	
	script:
	"""
	mkdir -p indexes/minimap2/${genome_id}/

	minimap2 -x sr -d indexes/minimap2/${genome_id}/${genome_id}.mmi ${fasta}
	"""
}