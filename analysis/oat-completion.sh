# You can add these completions by either adding
# source /path/to/oxDNA/analysis/oat-completions.sh 
# to your rc script
# Or by copying this script to your autocompletes folder (~/.bash_completion/ in bash)

_oat_completions(){
	local opts
    opts="align anm_parameterize backbone_flexibility bond_analysis centroid clustering config contact_map db2forces decimate deviations distance duplex_angle_plotter duplex_finder file_info forces2pairs generate_forces mean minify multidimensional_scaling_mean output_bonds oxDNA_PDB pairs2db pca persistence_length plot_energy subset_trajectory superimpose"

    if [ -n "$ZSH_VERSION" ]; then
        # For zsh, use compadd directly
        compadd -- ${=opts} -M "r:|=*"
    else
        # For bash, use compgen and COMPREPLY
        local cur="${COMP_WORDS[COMP_CWORD]}"
        COMPREPLY=($(compgen -W "$opts" -- "$cur"))
    fi
}

if [ -n "$BASH" ]; then
	complete -F _oat_completions oat
fi
if [ -n "$ZSH_NAME" ]; then
	compdef _oat_completions oat
fi
