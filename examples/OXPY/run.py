import numpy as np
import oxpy

with oxpy.Context():
    # init the manager with the given input file
    manager = oxpy.OxpyManager("input")

    # run 1k steps
    manager.run(1000)

    # run 10k steps more
    manager.run(10000)
    
    # do some computation with the current configuration
    particles = manager.config_info().particles()
    
    # compute the average position of the particles' backbones
    avg_pos = np.average(list(map(lambda p: p.backbone_site(), particles)), axis=0)
    print("Average final position:", avg_pos)
    
    # and the interaction energy between the first two particles
    print("Interaction energy between particle 0 and particle 1:", manager.config_info().interaction.pair_interaction(particles[0], particles[1]))
    