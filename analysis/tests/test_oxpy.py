import oxpy

with oxpy.Context():
    inp = oxpy.InputFile()
    inp.init_from_filename('input_rna')
    inp["trajectory_file"] = 'minitraj.dat'

    backend = oxpy.analysis.AnalysisBackend(inp)

    while backend.read_next_configuration():
        particles = backend.config_info().particles()
        print(particles[0].pos, particles[0].backbone_site())