#/usr/bin/env bash

_oat_completions(){
	if [ "${#COMP_WORDS[@]}" != "2" ]; then
		compopt -o default
		COMPREPLY=()
		return 0
	fi


	COMPREPLY=($(compgen -W "align anm_parameterize backbone_flexibility bond_analysis centroid clustering config contact_map db_to_force decimate deviations distance duplex_angle_plotter duplex_finder file_info forces2pairs generate_force mean minify multidimensional_scaling_mean output_bonds oxDNA_PDB pca persistence_length plot_energy subset_trajectory superimpose" "${COMP_WORDS[COMP_CWORD]}"))
}

complete -F _oat_completions oat
