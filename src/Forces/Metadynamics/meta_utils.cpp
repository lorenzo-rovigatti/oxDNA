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

    number contrib = 0.0;
    switch(settings.coord_mode) {
        case CoordSettings::CoordMode::HB_ENERGY: {
            LR_vector force, torque;
            number hb_energy = hb_interaction(pair.first, pair.second, force, torque, false);
            contrib = smooth_hb_contribution(settings.hb_energy_cutoff, settings.hb_transition_width, hb_energy);
            break;
        }
        case CoordSettings::CoordMode::SWITCHING_FUNCTION: {
            LR_vector r = CONFIG_INFO->box->min_image(pair.first->pos + pair.first->int_centers[DNANucleotide::BASE], pair.second->pos + pair.second->int_centers[DNANucleotide::BASE]);
            number r_mod = r.module();
            contrib =  1.0 / (1.0 + std::pow((r_mod - settings.d0) / settings.r0, settings.n));
            break;
        }
        case CoordSettings::CoordMode::MIXED: {
            LR_vector force, torque;
            number hb_energy = hb_interaction(pair.first, pair.second, force, torque, false);
            number hb_contribution = smooth_hb_contribution(settings.hb_energy_cutoff, settings.hb_transition_width, hb_energy);

            LR_vector r = CONFIG_INFO->box->min_image(pair.first->pos + pair.first->int_centers[DNANucleotide::BASE], pair.second->pos + pair.second->int_centers[DNANucleotide::BASE]);
            number r_mod = r.module();
            number switching_contribution = 1.0 / (1.0 + std::pow((r_mod - settings.d0) / settings.r0, settings.n));

            contrib = settings.mixed_weight * hb_contribution + (1.0 - settings.mixed_weight) * switching_contribution;
            break;
        }
    }

    return contrib;
}

std::pair<LR_vector, LR_vector> get_pair_force_torque_contribution(CoordSettings &settings, std::pair<BaseParticle*, BaseParticle*> &pair, BaseParticle *current_particle) {
    // Avoid dynamic loopkup by caching the hydrogen bonding function pointer. This is important since this function is called many times during the force computation. 
    static BaseInteraction::energy_function HB_function = CONFIG_INFO->interaction->get_interaction_function(DNAInteraction::HYDROGEN_BONDING);

    LR_vector force, torque;
    BaseParticle *other_particle = (current_particle == pair.first) ? pair.second : pair.first;

    switch(settings.coord_mode) {
        case CoordSettings::CoordMode::HB_ENERGY: {
            number hb_energy = hb_interaction(current_particle, other_particle, force, torque, true);
            number dcoord_dhb_energy = der_smooth_hb_contribution(settings.hb_energy_cutoff, settings.hb_transition_width, hb_energy);
            force *= dcoord_dhb_energy;
            torque *= dcoord_dhb_energy;
            break;
        }
        case CoordSettings::CoordMode::SWITCHING_FUNCTION: {
            LR_vector r_vec = CONFIG_INFO->box->min_image(current_particle->pos + current_particle->int_centers[DNANucleotide::BASE], other_particle->pos + other_particle->int_centers[DNANucleotide::BASE]);
            number r_mod = r_vec.module();

            number x = (r_mod - settings.d0) / settings.r0;
            number xn = std::pow(x, settings.n);
            number dcoord_dr = (settings.n / settings.r0) * std::pow(x, settings.n - 1) / (SQR(1.0 + xn));

            force = r_vec / r_mod * dcoord_dr;
            torque = current_particle->int_centers[DNANucleotide::BASE].cross(force);
            break;
        }
        case CoordSettings::CoordMode::MIXED: {
            number hb_energy = hb_interaction(current_particle, other_particle, force, torque, true);
            number dcoord_dhb_energy = der_smooth_hb_contribution(settings.hb_energy_cutoff, settings.hb_transition_width, hb_energy);
            LR_vector hb_force_contribution = force * dcoord_dhb_energy;
            LR_vector hb_torque_contribution = torque * dcoord_dhb_energy;

            LR_vector r_vec = CONFIG_INFO->box->min_image(current_particle->pos + current_particle->int_centers[DNANucleotide::BASE], other_particle->pos + other_particle->int_centers[DNANucleotide::BASE]);
            number r_mod = r_vec.module();

            number x = (r_mod - settings.d0) / settings.r0;
            number xn = std::pow(x, settings.n);
            number dcoord_dr = (settings.n / settings.r0) * std::pow(x, settings.n - 1) / (SQR(1.0 + xn));

            LR_vector switching_force_contribution = r_vec / r_mod * dcoord_dr;
            LR_vector switching_torque_contribution = current_particle->int_centers[DNANucleotide::BASE].cross(switching_force_contribution);

            force = settings.mixed_weight * hb_force_contribution + (1.0 - settings.mixed_weight) * switching_force_contribution;
            torque = settings.mixed_weight * hb_torque_contribution + (1.0 - settings.mixed_weight) * switching_torque_contribution;
            break;
        }
    }

    return std::make_pair(force, torque);
}

// here we got ridden of the smoothing function that is present in the original oxDNA code
number _f1(number r) {
    static number shift = HYDR_EPS_OXDNA2 * SQR(1.0 - exp(-(HYDR_RC - HYDR_R0) * HYDR_A));

	number val = 0.0;
	if(r < HYDR_RCHIGH) {
        number tmp = 1.0 - exp(-(r - HYDR_R0) * HYDR_A);
        val = HYDR_EPS_OXDNA2 * SQR(tmp) - shift;
	}

	return val;
}

number _f1D(number r) {
	number val = 0.0;
	if(r < HYDR_RCHIGH) {
        number tmp = exp(-(r - HYDR_R0) * HYDR_A);
        val = 2.0 * HYDR_EPS_OXDNA2 * (1 - tmp) * tmp * HYDR_A;
	}

	return val;
}

number _f4(number t, number t0, number a) {
	number val = 0.0;
	t -= t0;
	if(t < 0.0) {
		t *= -1.0;
    }

    // here we compute the cut-off "by hand", since the oxDNA cut-off takes 
    // into account also the smoothing function, which we don't use.
	if(t < 1.0 / std::sqrt(a)) {
        val = 1.0 - a * SQR(t);
	}

	return val;
}

number _f4Dsin(number t, number t0, number a) {
	number val = 0.0;
	number m = 1.0;
	number tt0 = t - t0;
	// this function is a parabola centered in t0. If t < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	if(tt0 < 0.0) {
		tt0 *= -1.0;
		m = -1.0;
	}

	if(tt0 < 1.0 / std::sqrt(a)) {
        number sint = sin(t);
        val = m * 2.0 * a * tt0 / sint;
	}

	return val;
}

number hb_interaction(BaseParticle *p, BaseParticle *q, LR_vector &force, LR_vector &torque, bool compute_force_torque) {
    force = torque = LR_vector();
    // true if p and q are Watson-Crick-like pairs
	bool is_pair = (q->btype + p->btype == 3);

    LR_vector r = CONFIG_INFO->box->min_image(p->pos, q->pos);
	LR_vector rhydro = r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
	number rhydromod = rhydro.module();
	number energy = (number) 0.f;
	if(is_pair && HYDR_RCLOW < rhydromod && rhydromod < HYDR_RCHIGH) {
		// vector, versor and magnitude of the base-base separation
		LR_vector rhydrodir = rhydro / rhydromod;

		// particle axes according to Allen's paper
		LR_vector &a1 = p->orientationT.v1;
		LR_vector &a3 = p->orientationT.v3;
		LR_vector &b1 = q->orientationT.v1;
		LR_vector &b3 = q->orientationT.v3;

		// angles involved in the HB interaction
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rhydrodir;
		number cost3 =  a1 * rhydrodir;

		number cost4 =  a3 * b3;
		number cost7 = -b3 * rhydrodir;
		number cost8 =  a3 * rhydrodir;

        number t1 = Utils::safe_acos(cost1);
        number t2 = Utils::safe_acos(cost2);
        number t3 = Utils::safe_acos(cost3);
        number t4 = Utils::safe_acos(cost4);
        number t7 = Utils::safe_acos(cost7);
        number t8 = Utils::safe_acos(cost8);

		// functions called at their relevant arguments
		number f1 = _f1(rhydromod);
		number f4t1 = _f4(t1, HYDR_THETA1_T0, HYDR_THETA1_A);
		number f4t2 = _f4(t2, HYDR_THETA2_T0, HYDR_THETA2_A);
		number f4t3 = _f4(t3, HYDR_THETA3_T0, HYDR_THETA3_A);

		number f4t4 = _f4(t4, HYDR_THETA4_T0, HYDR_THETA4_A);
		number f4t7 = _f4(t7, HYDR_THETA7_T0, HYDR_THETA7_A);
		number f4t8 = _f4(t8, HYDR_THETA8_T0, HYDR_THETA8_A);

		energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		// makes sense, since the above functions may return 0. exactly
		if(compute_force_torque && energy != 0.) {
			// derivatives called at the relevant arguments
			number f1D = _f1D(rhydromod);
			number f4t1Dsin =  _f4Dsin(t1, HYDR_THETA1_T0, HYDR_THETA1_A);
			number f4t2Dsin =  _f4Dsin(t2, HYDR_THETA2_T0, HYDR_THETA2_A);
			number f4t3Dsin = -_f4Dsin(t3, HYDR_THETA3_T0, HYDR_THETA3_A);

			number f4t4Dsin = -_f4Dsin(t4, HYDR_THETA4_T0, HYDR_THETA4_A);
			number f4t7Dsin =  _f4Dsin(t7, HYDR_THETA7_T0, HYDR_THETA7_A);
			number f4t8Dsin = -_f4Dsin(t8, HYDR_THETA8_T0, HYDR_THETA8_A);

			// RADIAL PART
			force = -rhydrodir * (f1D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);

			// THETA4; t4 = LRACOS (a3 * b3);
			LR_vector dir = a3.cross(b3);
			number torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8;

			torque += dir * torquemod;

			// THETA1; t1 = LRACOS (-a1 * b1);
			dir = a1.cross(b1);
			torquemod = -f1 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			torque += dir * torquemod;

			// THETA2; t2 = LRACOS (-b1 * rhydrodir);
			number fact = f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;
			force += (b1 + rhydrodir * cost2) * (fact / rhydromod);
            // t2 generates a torque on q, which is not needed here, so we don't compute it

			// THETA3; t3 = LRACOS (a1 * rhydrodir);
			force += (a1 - rhydrodir * cost3) * (f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8 / rhydromod);

			LR_vector t3dir = rhydrodir.cross(a1);
			torquemod = f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torque += t3dir * torquemod;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
			force += (b3 + rhydrodir * cost7) * (f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8 / rhydromod);
            // t7 generates a torque on q, which is not needed here, so we don't compute it

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			force += (a3 - rhydrodir * cost8) * (f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin / rhydromod);

			LR_vector t8dir = rhydrodir.cross(a3);
			torquemod = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torque += t8dir * torquemod;

			torque += p->int_centers[DNANucleotide::BASE].cross(force);
		}
	}

	return energy;
}

number der_smooth_hb_contribution(number hb_energy_cutoff, number hb_transition_width, number hb_energy) {
    number x = (hb_energy_cutoff - hb_energy) / hb_transition_width;
    
    // Clamp to avoid numerical overflow with tanh
    if(x > 10.0) return 0.0;
    if(x < -10.0) return 0.0;

    // derivative of the smooth step function with respect to the HB energy
    number tanh_x = std::tanh(x);
    return -0.5 * (1.0 - SQR(tanh_x)) / hb_transition_width;
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
