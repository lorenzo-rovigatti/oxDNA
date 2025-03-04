//
// Created by josh on 3/4/25.
//

#ifndef ESTIMATE_TM_PY_RASPBERRYPARTICLE_H
#define ESTIMATE_TM_PY_RASPBERRYPARTICLE_H

/**
 * extremely minimal patchy particle class (most particle data is stored in RaspberryInteraction class
 */
class RaspberryParticle : public BaseParticle{
public:
    RaspberryParticle(std::vector<LR_vector> geom) : _base_geometry(geom){
        int_centers.resize(_base_geometry.size());
    };
    void set_positions(){
        for (int i = 0; i < _base_geometry.size(); i++){
            int_centers[i] = orientation * _base_geometry[i];
        }
    }

    virtual bool is_rigid_body() {
        return true;
    }
protected:
    const std::vector<LR_vector> _base_geometry;
};


#endif //ESTIMATE_TM_PY_RASPBERRYPARTICLE_H
