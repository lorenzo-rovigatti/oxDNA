/*
 * meta_utils.cpp
 *
 *  Created on: Apr 3, 2022
 *      Author: lorenzo
 */

#include "meta_utils.h"
#include "../../Boxes/BaseBox.h"
#include "../../Utilities/ConfigInfo.h"
#include "../../Particles/DNANucleotide.h"
#include "../../Interactions/DNAInteraction.h"

#include <fast_double_parser/fast_double_parser.h>

namespace meta {

LR_vector particle_list_com(const std::vector<BaseParticle *> &list, BaseBox *box_ptr) {
	auto sum_abs_pos = [&](LR_vector sum, BaseParticle *p) { return sum + box_ptr->get_abs_pos(p); };
	return std::accumulate(list.begin(), list.end(), LR_vector(), sum_abs_pos) / (number) list.size();
}

std::tuple<std::vector<int>, std::vector<BaseParticle *>> get_particle_lists(input_file &inp, std::string key, std::vector<BaseParticle *> &particles, std::string description) {
	std::string p_string;
	getInputString(&inp, key.c_str(), p_string, 1);

	std::vector<int> idx_vector;
	std::vector<BaseParticle *> p_vector;
	idx_vector = Utils::get_particles_from_string(particles, p_string, description);
	for(auto p_idx : idx_vector) {
		p_vector.push_back(particles[p_idx]);
	}

	return std::make_tuple(idx_vector, p_vector);
}

number lexical_cast(const std::string &source) {
	double result;

	if(fast_double_parser::parse_number(source.c_str(), &result) == nullptr) {
		throw oxDNAException("Cannot convert '%s' to a number", source.c_str());
	}

	return result;
}

std::vector<number> split_to_numbers(const std::string &str, const std::string &delims) {
	std::vector<number> output;

	const char *ptr = str.c_str();
	while(ptr) {
		auto base = ptr;
		ptr = std::strpbrk(ptr, delims.c_str());
		if(ptr) {
			// this check makes sure that no empty strings are added to the output
			if(ptr - base) {
				output.emplace_back(lexical_cast(std::string(base, ptr - base)));
			}
			ptr++;
		}
		else {
			std::string remainder(base);
			if(remainder.size() > 0) {
				output.emplace_back(lexical_cast(remainder));
			}
		}
	}

	return output;
}

CoordSettings::CoordSettings() {
	coord_mode = CoordMode::HB_ENERGY;
	hb_energy_cutoff = -0.2;
	hb_transition_width = 0.1;
	d0 = 0.4;
	r0 = 0.5;
	n = 6;
}

void CoordSettings::get_settings(input_file &inp) {
	// Parse coordination mode
    std::string coord_type("hb_cutoff");
    getInputString(&inp, "coordination_type", coord_type, 0);
    
    // we need to parse all the parameters first, since they may be needed for multiple coordination modes
    // (e.g. if coord_type == "mixed" we need both the HB energy parameters and the switching function parameters)
    if(coord_type == "hb_cutoff" || coord_type == "mixed") {
        coord_mode = CoordMode::HB_ENERGY;
        // Get HB parameters
        getInputNumber(&inp, "hb_energy_cutoff", &hb_energy_cutoff, 0);
        getInputNumber(&inp, "hb_transition_width", &hb_transition_width, 0);
        OX_LOG(Logger::LOG_INFO, "Coordination: hb_energy_cutoff = %.2f, hb_transition_width = %.2f", hb_energy_cutoff, hb_transition_width);
    }
    if(coord_type == "switching_function" || coord_type == "mixed") {
        coord_mode = CoordMode::SWITCHING_FUNCTION;

        getInputNumber(&inp, "d0", &d0, 0);
        getInputNumber(&inp, "r0", &r0, 0);
        getInputInt(&inp, "n", &n, 0);

        if(n % 2 != 0) {
            throw oxDNAException("LTCoordination: exponent n must be an even integer");
        }
        OX_LOG(Logger::LOG_INFO, "Coordination: switching function with d0 = %.2f, r0 = %.2f, n = %d", d0, r0, n);
    } 

    if(coord_type == "hb_cutoff") {
        OX_LOG(Logger::LOG_INFO, "Coordination: Using HB energy cutoff coordination");
    }
    else if(coord_type == "switching_function") {
        OX_LOG(Logger::LOG_INFO, "Coordination: Using switching function coordination");
    }
    else if(coord_type == "mixed") {
        coord_mode = CoordMode::MIXED;
        getInputNumber(&inp, "mixed_weight", &mixed_weight, 1);
        if(mixed_weight < 0.0 || mixed_weight > 1.0) {
            throw oxDNAException("Coordination: mixed_weight must be between 0 and 1");
        }
        OX_LOG(Logger::LOG_INFO, "Coordination: Using mixed coordination with mixed_weight = %.2f", mixed_weight);
    }
    else {
        throw oxDNAException("Coordination: unknown coordination_type '%s' (valid options: 'switching_function', 'hb_cutoff', 'mixed')", coord_type.c_str());
    }
}

LR_vector distance(std::pair<BaseParticle *, BaseParticle *> &pair) {
    LR_vector p_base = pair.first->pos + pair.first->int_centers[DNANucleotide::BASE];
    LR_vector q_base = pair.second->pos + pair.second->int_centers[DNANucleotide::BASE];

    return CONFIG_INFO->box->min_image(p_base, q_base);
}

number coordination(CoordSettings &settings, std::vector<std::pair<BaseParticle *, BaseParticle *>> &all_pairs) {
    number coordination = 0.0;
    for(auto &pair : all_pairs) {
        coordination += get_pair_contribution(settings, pair);
    }
    return coordination;
}

number get_pair_contribution(CoordSettings &settings, std::pair<BaseParticle*, BaseParticle*> &pair) {
    // Avoid dynamic loopkup by caching the hydrogen bonding function pointer. This is important since this function is called many times during the force computation. 
        static BaseInteraction::energy_function HB_function = CONFIG_INFO->interaction->get_interaction_function(DNAInteraction::HYDROGEN_BONDING);

    switch(settings.coord_mode) {
        case CoordSettings::CoordMode::HB_ENERGY: {
            LR_vector r = CONFIG_INFO->box->min_image(pair.first->pos, pair.second->pos);
            CONFIG_INFO->interaction->set_computed_r(r);
            number hb_energy = HB_function(pair.first, pair.second, false, false);
            return smooth_hb_contribution(settings.hb_energy_cutoff, settings.hb_transition_width, hb_energy);
        }
        case CoordSettings::CoordMode::SWITCHING_FUNCTION: {
            number r_mod = distance(pair).module();
            return 1.0 / (1.0 + std::pow((r_mod - settings.d0) / settings.r0, settings.n));
        }
        case CoordSettings::CoordMode::MIXED: {
            // we compute both contributions and then mix them with the specified weight
            LR_vector r = CONFIG_INFO->box->min_image(pair.first->pos, pair.second->pos);
            CONFIG_INFO->interaction->set_computed_r(r);
            number hb_energy = HB_function(pair.first, pair.second, false, false);
            number hb_contribution = smooth_hb_contribution(settings.hb_energy_cutoff, settings.hb_transition_width, hb_energy);

            number r_mod = distance(pair).module();
            number switching_contribution = 1.0 / (1.0 + std::pow((r_mod - settings.d0) / settings.r0, settings.n));

            return settings.mixed_weight * hb_contribution + (1.0 - settings.mixed_weight) * switching_contribution;
        }
    }

    return 0.0; // should never be reached
}

number smooth_hb_contribution(number hb_energy_cutoff, number hb_transition_width, number hb_energy) {
    // Smooth step function that transitions around hb_energy_cutoff
    // Using tanh for smoothness: transitions from 1 (favorable HB, low energy) 
    // to 0 (unfavorable, high energy)
    number x = (hb_energy_cutoff - hb_energy) / hb_transition_width;
    
    // Clamp to avoid numerical overflow with tanh
    if(x > 10.0) return 1.0;
    if(x < -10.0) return 0.0;

    return 0.5 * (1.0 + tanh(x));
}

}
