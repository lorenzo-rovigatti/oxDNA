/*
 * FFS_MD_CUDAMixedBackend.cu
 *
 *  Created on: 18/apr/2013
 *      Author: Ben Snodin
 */

#include "FFS_MD_CUDAMixedBackend.h"
#include "../../Managers/SimManager.h"
#include "FFS_CUDA_MD.cuh"
#include "../../Interactions/DNAInteraction.h"
#include <sstream>

// adapted from FFS_MD_CPUBackend.cpp parsed_expression::parse_expression()
bool master_condition_element::parse_master_element(const char *expression, const char *condition_fname) {
	string s(expression);
	//remove spaces
	s.erase(remove_if(s.begin(), s.end(), static_cast<int(*)(int)>( isspace )), s.end());

	int condition_type = -1;
	int pos = 0;
	unsigned int i;
	for (i = 0; i < s.length() - 1; i++) {
		if (s[i] == '<') {
			if (s[i + 1] == '=') {
				condition_type = 1; // LEQ
				pos = i + 2;
				break;
			} else {
				condition_type = 2; // LT
				pos = i + 1;
				break;
			}
		} else if (s[i] == '>') {
			if (s[i + 1] == '=') {
				condition_type = 0; // GEQ
				pos = i + 2;
				break;
			} else {
				condition_type = 3; // GT
				pos = i + 1;
				break;
			}

		} else if (s[i] == '=' && s[i + 1] == '=') {
			condition_type = 4; // EQ
			pos = i + 2;
			break;
		}

	}
	if (condition_type == -1) {
		return false;
	}
	type = condition_type;

	stringstream str;
	str << s.substr(pos);
	str >> value;
	
	string parname = s.substr(0, i);
	strcpy(conditions_name, parname.c_str());

	return true;

}

// adapted from FFS_MD_CPUBackend.cpp parsed_condition::parse_condition()
bool master_condition::parse_master_condition(const char *expression, const char *condition_fname){
	string str(expression);

	//remove { and }
	str.erase(std::remove(str.begin(), str.end(), '{'), str.end());
	str.erase(std::remove(str.begin(), str.end(), '}'), str.end());

	//remove spaces
	str.erase(remove_if(str.begin(), str.end(), static_cast<int(*)(int)>( isblank )), str.end());



	std::stringstream input(str.c_str());
	std::string str_expression;

	while(std::getline(input,str_expression,'\n'))
	{
	    master_condition_element mce;

	    if(str_expression.length() > 1)
	    {
		    if (! mce.parse_master_element(str_expression.c_str(), condition_fname))
	      {
		      return false;
	      }
	      else
	      {
	     	elements.push_back(mce);
	      }
	    }

	}

	return true;
}

FFS_MD_CUDAMixedBackend::FFS_MD_CUDAMixedBackend() : CUDAMixedBackend(){
	_use_master_conditions = false;
	_gen_flux = false;
	_gen_flux_save_every = -1;
	_gen_flux_desired_cc = -1;
	_gen_flux_cross_count = -1;
	_gen_flux_saved_cross_count = -1;
	_gen_flux_debug = -1;
	_flux_direction = -1;
	_n_hb_pairs = -1;
	_n_hb_regions = -1;
	_n_dist_pairs = -1;
	_n_dist_regions = -1;
	_d_hb_pairs1 = NULL;
	_d_hb_pairs2 = NULL;
	_d_dist_pairs1 = NULL;
	_d_dist_pairs2 = NULL;
	_d_op_dists = NULL;
	_d_hb_energies = NULL;
	_d_dist_region_lens = NULL;
	_d_dist_region_rows = NULL;
	_d_hb_region_lens = NULL;
	_d_hb_region_rows = NULL;
	_d_hb_cutoffs = NULL;
	_h_region_is_nearhb = NULL;
	_unexpected_master_reached = false;
}


FFS_MD_CUDAMixedBackend::~FFS_MD_CUDAMixedBackend(){
	if (_gen_flux){
		if (_use_master_conditions){
			_free_master_conditions(_master_forward_conditions);
			_free_master_conditions(_master_backward_conditions);
		}
		else{
			_free_simple_conditions(_sc_fwd);
			_free_simple_conditions(_sc_bwd);
		}
	}
	else{
		if (_use_master_conditions){
			_free_master_conditions(_master_conditions);
		}
		else{
			_free_simple_conditions(_sc);
		}
	}
	CUDA_SAFE_CALL( cudaFree(_d_hb_pairs1) );
	CUDA_SAFE_CALL( cudaFree(_d_hb_pairs2) );
	CUDA_SAFE_CALL( cudaFree(_d_dist_pairs1) );
	CUDA_SAFE_CALL( cudaFree(_d_dist_pairs2) );
	CUDA_SAFE_CALL( cudaFree(_d_op_dists) );
	CUDA_SAFE_CALL( cudaFree(_d_hb_energies) );
	CUDA_SAFE_CALL( cudaFree(_d_dist_region_lens) );
	CUDA_SAFE_CALL( cudaFree(_d_dist_region_rows) );
	CUDA_SAFE_CALL( cudaFree(_d_hb_region_lens) );
	CUDA_SAFE_CALL( cudaFree(_d_hb_region_rows) );
	CUDA_SAFE_CALL( cudaFree(_d_hb_cutoffs) );
	CUDA_SAFE_CALL( cudaFree(_d_region_is_nearhb) );
	free(_h_region_is_nearhb);
	free(_h_hb_energies);
	free(_h_op_dists);
	free(_h_nearhb_states);
}


void FFS_MD_CUDAMixedBackend::_free_master_conditions(vector<master_condition> master_conditions){
	for (vector<master_condition>::iterator master_condition = (master_conditions).begin() ; master_condition != (master_conditions).end() ; master_condition++){
		for (vector<master_condition_element>::iterator mce = (*master_condition).elements.begin() ; mce != (*master_condition).elements.end() ; mce++){
			_free_simple_conditions((*mce).simple_conditions);
		}
	}
}

void FFS_MD_CUDAMixedBackend::_free_simple_conditions(SimpleConditions sc){
	free(sc.h_nearhb_cond_lens);
	free(sc.h_nearhb_cond_rows);
	free(sc.h_hb_cond_lens);
	free(sc.h_hb_cond_rows);
	free(sc.h_dist_cond_lens);
	free(sc.h_dist_cond_rows);
	free(sc.h_dist_cond_mags);
	free(sc.h_dist_cond_types);
	free(sc.h_hb_cond_mags);
	free(sc.h_hb_cond_types);
	free(sc.hb_cond_set);
	free(sc.h_nearhb_cond_mags);
	free(sc.h_nearhb_cond_types);
	free(sc.nearhb_cond_set);
	free(sc.dist_cond_set);
	free(sc.h_ffs_stop);
	CUDA_SAFE_CALL( cudaFree(sc.d_dist_cond_lens) );
	CUDA_SAFE_CALL( cudaFree(sc.d_dist_cond_rows) );
	CUDA_SAFE_CALL( cudaFree(sc.d_hb_cond_lens) );
	CUDA_SAFE_CALL( cudaFree(sc.d_hb_cond_rows) );
	CUDA_SAFE_CALL( cudaFree(sc.d_hb_cond_mags) );
	CUDA_SAFE_CALL( cudaFree(sc.d_hb_cond_types) );
	CUDA_SAFE_CALL( cudaFree(sc.d_nearhb_cond_lens) );
	CUDA_SAFE_CALL( cudaFree(sc.d_nearhb_cond_rows) );
	CUDA_SAFE_CALL( cudaFree(sc.d_nearhb_cond_mags) );
	CUDA_SAFE_CALL( cudaFree(sc.d_nearhb_cond_types) );
	CUDA_SAFE_CALL( cudaFree(sc.d_dist_cond_mags) );
	CUDA_SAFE_CALL( cudaFree(sc.d_dist_cond_types) );
	CUDA_SAFE_CALL( cudaFree(sc.d_ffs_stop) );
}

void FFS_MD_CUDAMixedBackend::sim_step(llint curr_step){
	// valgrind complains if _check_stop() happens before sim_step() has been called for the first time
	// so, I put sim_step before _check_stop, which has the disadvantage that the system will (I think) slightly advance
	// even if it should immediately stop (this is the same behaviour as with the CPU version)
	CUDAMixedBackend::sim_step(curr_step);
	if (_check_stop()) {
		SimManager::stop = true;
		OX_LOG(Logger::LOG_INFO,
				"Reached stop conditions, stopping in step %lld", curr_step);
		char tmp[1024];
		sprintf_names_and_values(tmp);
		OX_LOG(Logger::LOG_INFO,
				"FFS final values: %s", tmp);
	}
}

void FFS_MD_CUDAMixedBackend::get_settings(input_file &inp){
	CUDAMixedBackend::get_settings(inp);

	getInputString(&inp, "order_parameters_file", _order_parameters_file, 1);
	getInputString(&inp, "ffs_file", _ffs_file, 1);
	getInputInt(&inp, "print_energy_every", &_print_energy_every, 1);

	int tmpi;
	_gen_flux = false;
	if (getInputInt(&inp, "ffs_generate_flux", &tmpi, 0) == KEY_FOUND) {
		if (tmpi > 0){
			_gen_flux = true;
			getInputInt(&inp, "gen_flux_save_every", &_gen_flux_save_every, 1);
			getInputInt(&inp, "gen_flux_total_crossings", &_gen_flux_desired_cc, 1);
			int tmpi2;
			if (getInputInt(&inp, "gen_flux_debug", &tmpi2, 0) == KEY_FOUND){
				_gen_flux_debug = (int)tmpi2;
				if (_gen_flux_debug){
				OX_LOG(Logger::LOG_INFO, "saving configurations on backwards crossing for debugging");
				}
			}
			else{
				_gen_flux_debug = 0;
			}
			getInputString(&inp, "gen_flux_conf_prefix", _conf_prefix, 1);
			OX_LOG(Logger::LOG_INFO, "doing ffs flux generation, saving every %d to get a total of %d saved crossings", _gen_flux_save_every, _gen_flux_desired_cc);
			
			_check_initial_state = 0;
			if (getInputInt(&inp, "check_initial_state", &tmpi2, 0) == KEY_FOUND){
				_check_initial_state = (int)tmpi2;
				if (_check_initial_state){
					OX_LOG(Logger::LOG_INFO, "Checking state of initial configuration; if it satisfies the forward conditions at t=0, the code will complain and exit");
				}
			}

			_die_on_unexpected_master = 0;
			if (getInputInt(&inp, "die_on_unexpected_master", &tmpi2, 0) == KEY_FOUND) {
				if (tmpi2 > 0){
					_die_on_unexpected_master = (int)tmpi2;
					if (_die_on_unexpected_master){
						getInputString(&inp, "unexpected_master_prefix", _unexpected_master_prefix, 1);
						OX_LOG(Logger::LOG_INFO, "Code will exit and save the configuration to a file with prefix %s if a master condition other than master_forward1 or master_backward1 is satisfied", _unexpected_master_prefix);
					}
				}
			}
		}
	}
	
	// initialize the ObservableOutput that prints custom configurations
	// print_every = 0 makes it not to print automatically
	// the value of this variable will be overwritten each time _prepare_configuration is called
	// so we can put whatever we want
	char customconf_file[256] = "custom.dat";
	std::string fake = Utils::sformat("{\n\tname = %s\n\tprint_every = 0\n\tonly_last = 1\n}\n", customconf_file);
	_obs_output_custom_conf = new ObservableOutput<float>(fake, inp);
	_obs_output_custom_conf->add_observable("type = configuration");
	this->_obs_outputs.push_back(_obs_output_custom_conf);
}

void FFS_MD_CUDAMixedBackend::init(){
	CUDAMixedBackend::init();
	_op.init_from_file(_order_parameters_file, this->_particles, this->_N);

	// initialise the order parameter region arrays (variable width, pseudo-2D arrays)
	int *h_dist_region_lens;
	int *h_dist_region_rows;
	int *h_hb_region_lens;
	int *h_hb_region_rows;
	_n_dist_pairs = _op.get_dist_pairs_total_count();
	_n_hb_pairs = _op.get_hb_pairs_total_count();
	_n_dist_regions = _op.get_distance_parameters_count();
	_n_hb_regions = _op.get_hb_parameters_count();
	h_dist_region_lens = (int *)malloc(_n_dist_regions * sizeof(int));
	h_hb_region_lens = (int *)malloc(_n_hb_regions * sizeof(int));
	_op.get_dist_pairs_count(h_dist_region_lens);
	_op.get_hb_pairs_count(h_hb_region_lens);
	h_dist_region_rows = _get_2D_rows(_n_dist_regions, h_dist_region_lens);
	h_hb_region_rows = _get_2D_rows(_n_hb_regions, h_hb_region_lens);

	int *h_hb_pairs1;
	int *h_hb_pairs2;
	int *h_dist_pairs1;
	int *h_dist_pairs2;
	float *h_hb_cutoffs;

	h_hb_pairs1 = (int *)malloc(_n_hb_pairs * sizeof(int));
	h_hb_pairs2 = (int *)malloc(_n_hb_pairs * sizeof(int));
	h_dist_pairs1 = (int *)malloc(_n_dist_pairs * sizeof(int));
	h_dist_pairs2 = (int *)malloc(_n_dist_pairs * sizeof(int));
	h_hb_cutoffs = (float *)malloc(_n_hb_regions * sizeof(float));
	_h_hb_energies = (float *)malloc(_n_hb_pairs * sizeof(float));
	_h_op_dists = (float *)malloc(_n_dist_pairs * sizeof(float));
	_h_nearhb_states = (bool *)malloc(_n_hb_pairs * sizeof(bool));
	_h_region_is_nearhb = (bool *)malloc(_n_hb_regions * sizeof(bool));

	// array of hb cutoff values for each hb order parameter
	for (int ii=0 ; ii < _n_hb_regions ; ii++ ){
		h_hb_cutoffs[ii] = _op.get_hb_cutoff(ii);
	}

	// these arrays contain the particle indices for the particles involved in each order parameter
	_op.get_hb_pairs(h_hb_pairs1, h_hb_pairs2);
	_op.get_dist_pairs(h_dist_pairs1, h_dist_pairs2);

	// check for any 'nearly hydrogen bonded' order parameters
	for (int ii=0 ; ii < _n_hb_regions ; ii++){
		if ((int)h_hb_cutoffs[ii] == 64){
			OX_LOG(Logger::LOG_INFO, "Using near hbond for bond order parameter number %d", ii+1);
			_h_region_is_nearhb[ii] = true;
		}
		else{
			_h_region_is_nearhb[ii] = false;
		}
	}

	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_hb_region_rows, _n_hb_regions * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_hb_region_lens, _n_hb_regions * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_dist_region_rows, _n_dist_regions * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_dist_region_lens, _n_dist_regions * sizeof(int)) );

	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_hb_pairs1, _n_hb_pairs * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_hb_pairs2, _n_hb_pairs * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<float>(&_d_hb_energies, _n_hb_pairs * sizeof(float)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<bool>(&_d_nearhb_states, _n_hb_pairs * sizeof(bool)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_dist_pairs1, _n_dist_pairs * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_dist_pairs2, _n_dist_pairs * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<float>(&_d_op_dists, _n_dist_pairs * sizeof(float)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<float>(&_d_hb_cutoffs, _n_hb_regions * sizeof(float)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<bool>(&_d_region_is_nearhb, _n_hb_regions * sizeof(bool)) );

	CUDA_SAFE_CALL( cudaMemcpy(_d_dist_region_lens, h_dist_region_lens, _n_dist_regions * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_dist_region_rows, h_dist_region_rows, _n_dist_regions * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_hb_region_lens, h_hb_region_lens, _n_hb_regions * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_hb_region_rows, h_hb_region_rows, _n_hb_regions * sizeof(int), cudaMemcpyHostToDevice) );

	CUDA_SAFE_CALL( cudaMemcpy(_d_hb_pairs1, h_hb_pairs1, _n_hb_pairs * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_hb_pairs2, h_hb_pairs2, _n_hb_pairs * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_dist_pairs1, h_dist_pairs1, _n_dist_pairs * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_dist_pairs2, h_dist_pairs2, _n_dist_pairs * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_hb_cutoffs, h_hb_cutoffs, _n_hb_regions * sizeof(float), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_region_is_nearhb, _h_region_is_nearhb, _n_hb_regions * sizeof(bool), cudaMemcpyHostToDevice) );

	_init_ffs_kernel_config(&_ffs_hb_precalc_kernel_cfg, _n_hb_pairs);
	_init_ffs_kernel_config(&_ffs_dist_precalc_kernel_cfg, _n_dist_pairs);
	_init_ffs_kernel_config(&_ffs_hb_eval_kernel_cfg, _n_hb_regions);
	_init_ffs_kernel_config(&_ffs_dist_eval_kernel_cfg, _n_dist_regions);

	_init_ffs_from_file(_ffs_file);

	if (_gen_flux){
		_flux_direction = FORWARD;
		_gen_flux_cross_count = 0;
		_gen_flux_saved_cross_count = 0;
		if (!_use_master_conditions){
			_sc_fwd = _get_simple_conditions(_forward_conditions, "forward");
			_sc_bwd = _get_simple_conditions(_backward_conditions, "backward");
		}
	}
	else{
		if (!_use_master_conditions){
			_sc = _get_simple_conditions(_conditions, "condition");
		}
	}

	free(h_dist_region_lens);
	free(h_dist_region_rows);
	free(h_dist_pairs1);
	free(h_dist_pairs2);
	free(h_hb_region_lens);
	free(h_hb_region_rows);
	free(h_hb_pairs1);
	free(h_hb_pairs2);
	free(h_hb_cutoffs);

	cudaThreadSynchronize();
}

SimpleConditions FFS_MD_CUDAMixedBackend::_get_simple_conditions(std::vector<parsed_condition> conditions, const char type[256]){
	// For the stopping conditions: convert the nice C++ vector of structs into a struct of arrays. The arrays are then copied to the device.
	
	SimpleConditions sc;

	strcpy(sc.type, type);
	// _dist_cond_len: total number of distance conditions
	// _h_dist_cond_lens: array of distance condition length for each distance region (order parameter).
	sc.dist_cond_len = 0;
	sc.h_dist_cond_lens = (int *)malloc(_n_dist_regions * sizeof(int));
	for (int ii=0 ; ii<_n_dist_regions; ii++)
		sc.h_dist_cond_lens[ii] = 0;
	sc.hb_cond_len = 0;
	sc.h_hb_cond_lens = (int *)malloc(_n_hb_regions * sizeof(int));
	for (int ii=0 ; ii<_n_hb_regions; ii++)
		sc.h_hb_cond_lens[ii] = 0;
	sc.nearhb_cond_len = 0;
	sc.h_nearhb_cond_lens = (int *)malloc(_n_hb_regions * sizeof(int));
	for (int ii=0 ; ii<_n_hb_regions; ii++)
		sc.h_nearhb_cond_lens[ii] = 0;


	// it might be possible to do this in a better way
	// count up the number of conditions for each region/(order parameter name)

	sc.stop_set_count = 0;
	int par_ind;
	for (vector<parsed_condition>::iterator jj = conditions.begin() ; jj != conditions.end() ; jj++){
		sc.stop_set_count += 1; // count up the number of sets of stopping conditions
		for (vector<parsed_expression>::iterator i = (*jj).all_expressions.begin() ; i != (*jj).all_expressions.end() ; i++){
			if ((*i).expression_type == 0){
				par_ind = (*i).parameter_index;
				if (_h_region_is_nearhb[par_ind]){
					sc.h_nearhb_cond_lens[par_ind]++;
					sc.nearhb_cond_len++;
				}
				else{
					sc.h_hb_cond_lens[par_ind]++;
					sc.hb_cond_len++;
				}
			}
			else if ((*i).expression_type == 1){
				par_ind = (*i).parameter_index;
				sc.h_dist_cond_lens[par_ind]++;
				sc.dist_cond_len++;
			}
			else{
				throw oxDNAException("unknown expression_type %d, dying now\n", (*i).expression_type);
			}
		}
	}
	sc.stop_length = sc.dist_cond_len + sc.hb_cond_len + sc.nearhb_cond_len;
	sc.h_ffs_stop = (bool *)malloc(sc.stop_length * sizeof(bool));

	sc.h_dist_cond_mags = (float *)malloc(sc.dist_cond_len * sizeof(float));
	sc.h_dist_cond_types = (int *)malloc(sc.dist_cond_len * sizeof(int));
	sc.h_hb_cond_mags = (float *)malloc(sc.hb_cond_len * sizeof(float));
	sc.h_hb_cond_types = (int *)malloc(sc.hb_cond_len * sizeof(int));
	sc.h_nearhb_cond_mags = (float *)malloc(sc.nearhb_cond_len * sizeof(float));
	sc.h_nearhb_cond_types = (int *)malloc(sc.nearhb_cond_len * sizeof(int));

	sc.h_dist_cond_rows = _get_2D_rows(_n_dist_regions, sc.h_dist_cond_lens);
	sc.h_hb_cond_rows = _get_2D_rows(_n_hb_regions, sc.h_hb_cond_lens);
	sc.h_nearhb_cond_rows = _get_2D_rows(_n_hb_regions, sc.h_nearhb_cond_lens);
	int *dist_inds = (int *)malloc(_n_dist_regions * sizeof(int));
	for (int ii=0 ; ii<_n_dist_regions; ii++){
		dist_inds[ii] = 0;
	}
	int *hb_inds = (int *)malloc(_n_hb_regions * sizeof(int));
	for (int ii=0 ; ii<_n_hb_regions; ii++){
		hb_inds[ii] = 0;
	}
	int *nearhb_inds = (int *)malloc(_n_hb_regions * sizeof(int));
	for (int ii=0 ; ii<_n_hb_regions; ii++){
		nearhb_inds[ii] = 0;
	}
	sc.hb_cond_set = (int *)malloc(sc.hb_cond_len*sizeof(int));
	sc.nearhb_cond_set = (int *)malloc(sc.nearhb_cond_len*sizeof(int));
	sc.dist_cond_set = (int *)malloc(sc.dist_cond_len*sizeof(int));
	// get the data I need from the parsed_expression structs held in the conditions vector
	int set_index = 0;
	for (vector<parsed_condition>::iterator ii = conditions.begin() ; ii != conditions.end() ; ii++){
		for (vector<parsed_expression>::iterator jj = (*ii).all_expressions.begin() ; jj != (*ii).all_expressions.end() ; jj++){
			if ((*jj).expression_type == 0){
				// hydrogen bond expression
				par_ind = (*jj).parameter_index;
				if (_h_region_is_nearhb[par_ind]){
					sc.h_nearhb_cond_mags[sc.h_nearhb_cond_rows[par_ind]+nearhb_inds[par_ind]] = (*jj).value_to_compare_with_hb;
					sc.h_nearhb_cond_types[sc.h_nearhb_cond_rows[par_ind]+nearhb_inds[par_ind]] = (*jj).compare_type;
					sc.nearhb_cond_set[sc.h_nearhb_cond_rows[par_ind]+nearhb_inds[par_ind]] = set_index;
					nearhb_inds[par_ind]++;
				}
				else{
					sc.h_hb_cond_mags[sc.h_hb_cond_rows[par_ind]+hb_inds[par_ind]] = (*jj).value_to_compare_with_hb;
					sc.h_hb_cond_types[sc.h_hb_cond_rows[par_ind]+hb_inds[par_ind]] = (*jj).compare_type;
					sc.hb_cond_set[sc.h_hb_cond_rows[par_ind]+hb_inds[par_ind]] = set_index;
					hb_inds[par_ind]++;
				}
			}
			else if ((*jj).expression_type == 1){
				// mindist expression
				par_ind = (*jj).parameter_index;
				sc.h_dist_cond_mags[sc.h_dist_cond_rows[par_ind]+dist_inds[par_ind]] = (*jj).value_to_compare_with_dist;
				sc.h_dist_cond_types[sc.h_dist_cond_rows[par_ind]+dist_inds[par_ind]] = (*jj).compare_type;
				sc.dist_cond_set[sc.h_dist_cond_rows[par_ind]+dist_inds[par_ind]] = set_index;
				dist_inds[par_ind]++;
			}
		}
		set_index += 1;
	}
	free(nearhb_inds);
	free(hb_inds);
	free(dist_inds);

	// rows and lens arrays help with the variable-width '2d' array implementation
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<float>(&sc.d_dist_cond_mags, sc.dist_cond_len * sizeof(float)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&sc.d_dist_cond_types, sc.dist_cond_len * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&sc.d_dist_cond_lens, _n_dist_regions * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&sc.d_dist_cond_rows, _n_dist_regions * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<float>(&sc.d_hb_cond_mags, sc.hb_cond_len * sizeof(float)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&sc.d_hb_cond_types, sc.hb_cond_len * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&sc.d_hb_cond_lens, _n_hb_regions * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&sc.d_hb_cond_rows, _n_hb_regions * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<float>(&sc.d_nearhb_cond_mags, sc.nearhb_cond_len * sizeof(float)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&sc.d_nearhb_cond_types, sc.nearhb_cond_len * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&sc.d_nearhb_cond_lens, _n_hb_regions * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&sc.d_nearhb_cond_rows, _n_hb_regions * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<bool>(&sc.d_ffs_stop, sc.stop_length * sizeof(bool)) );

	CUDA_SAFE_CALL( cudaMemcpy(sc.d_hb_cond_mags, sc.h_hb_cond_mags, sc.hb_cond_len * sizeof(float), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(sc.d_hb_cond_types, sc.h_hb_cond_types, sc.hb_cond_len * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(sc.d_hb_cond_lens, sc.h_hb_cond_lens, _n_hb_regions * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(sc.d_hb_cond_rows, sc.h_hb_cond_rows, _n_hb_regions * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(sc.d_nearhb_cond_mags, sc.h_nearhb_cond_mags, sc.nearhb_cond_len * sizeof(float), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(sc.d_nearhb_cond_types, sc.h_nearhb_cond_types, sc.nearhb_cond_len * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(sc.d_nearhb_cond_lens, sc.h_nearhb_cond_lens, _n_hb_regions * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(sc.d_nearhb_cond_rows, sc.h_nearhb_cond_rows, _n_hb_regions * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(sc.d_dist_cond_mags, sc.h_dist_cond_mags, sc.dist_cond_len * sizeof(float), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(sc.d_dist_cond_types, sc.h_dist_cond_types, sc.dist_cond_len * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(sc.d_dist_cond_lens, sc.h_dist_cond_lens, _n_dist_regions * sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(sc.d_dist_cond_rows, sc.h_dist_cond_rows, _n_dist_regions * sizeof(int), cudaMemcpyHostToDevice) );

	return sc;
}

void FFS_MD_CUDAMixedBackend::_init_CUDA_MD_symbols() {
	// this line tells the mixed backend to initialize the constants needed by its kernels (which are in CUDA_mixed.cuh), so that sim_step will work properly
	CUDAMixedBackend::_init_CUDA_MD_symbols();
}

bool FFS_MD_CUDAMixedBackend::_read_conditions(const char *fname, const char *condition_set_type, std::vector<parsed_condition> *conditions) {
	bool return_val = true;
	FILE * fin = 0;
	fin = fopen(fname, "r");
	if (fin == NULL)
		throw oxDNAException("Cannot open %s", fname);

	input_file input;

	loadInput(&input, fin);

	char pairname[512];
	string strexpr;
	int pairid = 1;
	sprintf(pairname, "%s%d", condition_set_type, pairid);
	while (getInputString(&input, pairname, strexpr, 0) == KEY_FOUND) {
		parsed_condition newcondition;
		if (strcmp(condition_set_type, "master") == 0){
			// special case for the master 'conditions' since we don't expect an entry in the order parameter file
			// corresponding to each master condition element
			newcondition.parse_condition(strexpr.c_str(), &_op);
			for (std::vector<parsed_expression>::iterator expression=newcondition.all_expressions.begin() ; expression != newcondition.all_expressions.end() ; expression++){
				if ((*expression).compare_type == -1)
					throw oxDNAException("Failed to parse condition %s, please check the file format and parameter names ", pairname);
			}
		}
		else {
			if (!newcondition.parse_condition(strexpr.c_str(), &_op))
				throw oxDNAException("Failed to parse condition %s, please check the file format and parameter names ", pairname);
		}
		(*conditions).push_back(newcondition);
		pairid++;
		sprintf(pairname, "%s%d", condition_set_type, pairid);
	}

	if ((*conditions).size() == 0) throw oxDNAException("Found no conditions while parsing conditions of type %s, dying now\n", condition_set_type);
	
	cleanInputFile (&input);
	fclose(fin);
	return return_val;
}

// overloaded for master conditions
bool FFS_MD_CUDAMixedBackend::_read_conditions(const char *fname, const char *condition_set_type, std::vector<master_condition> *conditions) {
	bool return_val = true;
	FILE * fin = 0;
	fin = fopen(fname, "r");
	if (fin == NULL)
		throw oxDNAException("Cannot open %s", fname);

	input_file input;
	loadInput(&input, fin);

	char pairname[512];
	string strexpr;
	int pairid = 1;
	sprintf(pairname, "%s%d", condition_set_type, pairid);
	while (getInputString(&input, pairname, strexpr, 0) == KEY_FOUND) {
		master_condition newcondition;
		if (!newcondition.parse_master_condition(strexpr.c_str(), fname)){
			throw oxDNAException("Failed to parse condition %s, please check the file format and parameter names ", pairname);
		}
		strcpy(newcondition.name, pairname);
		(*conditions).push_back(newcondition);
		pairid++;
		sprintf(pairname, "%s%d", condition_set_type, pairid);
	}

	if ((*conditions).size() == 0) {
		if (_use_master_conditions){
			if ((*conditions).size() == 0) throw oxDNAException("Found no conditions while parsing conditions of type %s, dying now\n", condition_set_type);
		}
		else{
			// since use_master_conditions isn't set, it's ok for us to find no entries
			return_val = false;
		}
	}

	cleanInputFile (&input);
	fclose(fin);
	return return_val;
}

// prepare the simple_conditions for each master condition element for each master condition
void FFS_MD_CUDAMixedBackend::_master_conditions_prepare_simple_conditions(vector<master_condition> *master_conditions, const char *fname) {
	for (vector<master_condition>::iterator master_condition = (*master_conditions).begin() ; master_condition != (*master_conditions).end() ; master_condition++){
		for (vector<master_condition_element>::iterator mce = (*master_condition).elements.begin() ; mce != (*master_condition).elements.end() ; mce++){
			std::vector<parsed_condition> element_conditions;
			_read_conditions(fname, (*mce).conditions_name, &element_conditions);
			(*mce).simple_conditions = _get_simple_conditions(element_conditions, (*mce).conditions_name);
		}
	}
}

void FFS_MD_CUDAMixedBackend::_init_ffs_from_file(const char *fname) {
	if (_gen_flux){
		if (_read_conditions(fname, "master_forward", &_master_forward_conditions)){
			_use_master_conditions = true;
			OX_LOG(Logger::LOG_INFO, "Using master conditions for ffs simulation");
			_read_conditions(fname, "master_backward", &_master_backward_conditions);
			_master_conditions_prepare_simple_conditions(&_master_forward_conditions, fname);
			_master_conditions_prepare_simple_conditions(&_master_backward_conditions, fname);
		}
		else{
			_read_conditions(fname, "forward", &_forward_conditions);
			_read_conditions(fname, "backward", &_backward_conditions);
		}
	}
	else{
		if (_read_conditions(fname, "master", &_master_conditions)){
			_use_master_conditions = true;
			OX_LOG(Logger::LOG_INFO, "Using master conditions for ffs simulation");
			_master_conditions_prepare_simple_conditions(&_master_conditions, fname);
		}
		else{
			_read_conditions(fname, "condition", &_conditions);
		}
	}
}

void FFS_MD_CUDAMixedBackend::_init_ffs_kernel_config(CUDA_kernel_cfg *kernel_cfg, int total_threads) {
	(*kernel_cfg).threads_per_block = 2*this->_device_prop.warpSize;
	(*kernel_cfg).blocks.x = total_threads / (*kernel_cfg).threads_per_block + ((total_threads % (*kernel_cfg).threads_per_block == 0) ? 0 : 1);
	if((*kernel_cfg).blocks.x == 0) (*kernel_cfg).blocks.x = 1;
	(*kernel_cfg).blocks.y = (*kernel_cfg).blocks.z = 1;
}

void FFS_MD_CUDAMixedBackend::_eval_order_parameter_states(){
	this->_cuda_interaction->_hb_op_precalc(this->_d_poss, this->_d_orientations, _d_hb_pairs1, _d_hb_pairs2, _d_hb_energies, _n_hb_pairs, _d_region_is_nearhb, _ffs_hb_precalc_kernel_cfg);

	this->_cuda_interaction->_near_hb_op_precalc(this->_d_poss, this->_d_orientations, _d_hb_pairs1, _d_hb_pairs2, _d_nearhb_states, _n_hb_pairs, _d_region_is_nearhb,  _ffs_hb_precalc_kernel_cfg);

	this->_cuda_interaction->_dist_op_precalc(this->_d_poss, this->_d_orientations, _d_dist_pairs1, _d_dist_pairs2, _d_op_dists, _n_dist_pairs, _ffs_dist_precalc_kernel_cfg);
	cudaThreadSynchronize();
}

void FFS_MD_CUDAMixedBackend::_eval_stop_conditions(SimpleConditions sc){
	// check the values of the order parameters against the stopping conditions given -- do it on the GPU to reduce device-host data transfer
	// the result is stored in sc.h_ffs_stop
	hb_op_eval<float, float4><<<_ffs_hb_eval_kernel_cfg.blocks, _ffs_hb_eval_kernel_cfg.threads_per_block>>>(_d_hb_energies, _d_hb_region_lens, _d_hb_region_rows, sc.d_hb_cond_lens, sc.d_hb_cond_rows, sc.d_hb_cond_mags, sc.d_hb_cond_types, sc.d_ffs_stop, _n_hb_regions, _d_hb_cutoffs);
	CUT_CHECK_ERROR("hb_op_eval error");

//	nearhb_op_eval<float, float4><<<_ffs_nearhb_eval_kernel_cfg.blocks, _ffs_nearhb_eval_kernel_cfg.threads_per_block>>>(_d_nearhb_states, _d_nearhb_region_lens, _d_nearhb_region_rows, sc.d_nearhb_cond_lens, sc.d_nearhb_cond_rows, sc.d_nearhb_cond_mags, sc.d_nearhb_cond_types, sc.d_ffs_stop, _n_nearhb_regions, sc.hb_cond_len + sc.dist_cond_len);
	nearhb_op_eval<float, float4><<<_ffs_hb_eval_kernel_cfg.blocks, _ffs_hb_eval_kernel_cfg.threads_per_block>>>(_d_nearhb_states, _d_hb_region_lens, _d_hb_region_rows, sc.d_nearhb_cond_lens, sc.d_nearhb_cond_rows, sc.d_nearhb_cond_mags, sc.d_nearhb_cond_types, sc.d_ffs_stop, _n_hb_regions, sc.hb_cond_len + sc.dist_cond_len);
	CUT_CHECK_ERROR("nearhb_op_eval error");

	dist_op_eval<float, float4><<<_ffs_dist_eval_kernel_cfg.blocks, _ffs_dist_eval_kernel_cfg.threads_per_block>>>(_d_op_dists, _d_dist_region_lens, _d_dist_region_rows, sc.d_dist_cond_lens, sc.d_dist_cond_rows, sc.d_dist_cond_mags, sc.d_dist_cond_types, sc.d_ffs_stop, _n_dist_regions, sc.hb_cond_len);
	CUT_CHECK_ERROR("dist_op_eval error");
	cudaThreadSynchronize();

	CUDA_SAFE_CALL( cudaMemcpy(sc.h_ffs_stop, sc.d_ffs_stop, sc.stop_length * sizeof(bool), cudaMemcpyDeviceToHost) );
}

int FFS_MD_CUDAMixedBackend::_test_crossing(SimpleConditions sc, bool suppress_logging){
	// returns the number of conditions that are satisfied
	// this could be done on the device to reduce the memory transfer per step - is it worth it?
	int success_count = 0;
	// there are (potentially) many sets of stopping conditions; stop if a particular set's conditions are all true
	int this_set;
	bool *set_stop = (bool *)malloc(sc.stop_set_count*sizeof(bool));
	// loop for hydrogen bond conditions
	for (int ii=0; ii<sc.stop_set_count; ii++)
		set_stop[ii] = true;
	for (int ii=0; ii<sc.hb_cond_len; ii++){
		this_set = sc.hb_cond_set[ii];
		if (!sc.h_ffs_stop[ii])
			set_stop[this_set] = false;
	}
	// loop for distance conditions
	for (int ii=0; ii<sc.dist_cond_len; ii++){
		this_set = sc.dist_cond_set[ii];
		if (!sc.h_ffs_stop[ii+sc.hb_cond_len])
			set_stop[this_set] = false;
	}
	// loop for nearly hydrogen bonded conditions
	for (int ii=0; ii<sc.nearhb_cond_len; ii++){
		this_set = sc.nearhb_cond_set[ii];
		if (!sc.h_ffs_stop[ii+sc.hb_cond_len+sc.dist_cond_len])
			set_stop[this_set] = false;
	}
	// for every interface
	for (int ii=0; ii<sc.stop_set_count; ii++){
		if (set_stop[ii]){
			success_count++;
			if (!suppress_logging) OX_LOG(Logger::LOG_INFO, "Reached %s%d", sc.type, ii+1);
		}
	}
	free(set_stop);
	return success_count;
}

bool FFS_MD_CUDAMixedBackend::_test_master_condition(master_condition master_condition){
	// for each element of a particular master condition, check the state of all of the element's conditions
	// if the right number of each element's conditions are satisfied, return true, otherwise return false
	SimpleConditions simple_conditions;
	int current_value;
	int this_value;
	int this_type;
	bool stop_sim = true;
	for (vector<master_condition_element>::iterator master_element = master_condition.elements.begin() ; master_element != master_condition.elements.end() ; master_element++){
		simple_conditions = (*master_element).simple_conditions;
		// it might speed things up to do all the eval_stop_conditions together and do the sync and copy after they're all done
		_eval_stop_conditions(simple_conditions);
		current_value = _test_crossing(simple_conditions, true);
		this_value = (*master_element).value;
		this_type = (*master_element).type;
		switch (this_type) {
		case 0: //GEQ
			if (!(current_value >= this_value)){
				stop_sim = false;
			}
			break;
		case 1: //LEQ
			if (!(current_value <= this_value)){
				stop_sim = false;
			}
			break;
		case 2: //LT
			if (!(current_value < this_value)){
				stop_sim = false;
			}
			break;
		case 3: //GT
			if (!(current_value > this_value)){
				stop_sim = false;
			}
			break;
		case 4: //EQ
			if (!(current_value == this_value)){
				stop_sim = false;
			}
			break;
		}
	}
	return stop_sim;
}

bool FFS_MD_CUDAMixedBackend::_test_master_conditions(vector<master_condition> master_conditions)
{
	for (vector<master_condition>::iterator master_condition = master_conditions.begin() ; master_condition != master_conditions.end() ; master_condition++){
		if (_test_master_condition(*master_condition)){
			OX_LOG(Logger::LOG_INFO, "Reached %s", (*master_condition).name);
			_log_master_state(*master_condition);
			// check for unexpected master condition in flux generation
			if (_gen_flux){
				int master_type;
				if (sscanf((*master_condition).name, "master_forward%d", &master_type)){
					if (master_type > 1){
						_unexpected_master_reached = true;
						strcpy(_unexpected_master_name, (*master_condition).name);
					}
				}
				else if (sscanf((*master_condition).name, "master_backward%d", &master_type)){
					if (master_type > 1){
						_unexpected_master_reached = true;
						strcpy(_unexpected_master_name, (*master_condition).name);
					}
				}
				else{
					throw oxDNAException("unknown master condition name '%s', dying\n", (*master_condition).name); // should never happen
				}
			}
			return true;
		}
	}
	return false;
}
void FFS_MD_CUDAMixedBackend::_log_master_state(master_condition master_condition){
	// intended for use after a master condition has been satisfied; loop through all elements and their conditions
	// and report which conditions have been reached
	SimpleConditions simple_conditions;
	for (vector<master_condition_element>::iterator master_element = master_condition.elements.begin() ; master_element != master_condition.elements.end() ; master_element++){
		simple_conditions = (*master_element).simple_conditions;
		OX_LOG(Logger::LOG_INFO, "States for master condition element %s:", simple_conditions.type);
		_eval_stop_conditions(simple_conditions);
		_test_crossing(simple_conditions);
	}
}

void FFS_MD_CUDAMixedBackend::_handle_unexpected_master(){
	// call when an unexpected master condition is reached: print a configuration and print a message to screen
	char conf_str[256];
	sprintf(conf_str, "%s_%s_N%d.dat", _unexpected_master_prefix, _unexpected_master_name, this->_gen_flux_saved_cross_count);
	_prepare_configuration(conf_str);
	_obs_output_custom_conf->print_output(_curr_step);
	OX_LOG(Logger::LOG_INFO, "saved configuration %s at step %d", conf_str, this->_curr_step);
	_gen_flux_saved_cross_count += 1;
	OX_LOG(Logger::LOG_INFO, "Crossed interface with master condition integer other than 1 and die_on_unexpected_master is set; exiting now");
}



bool FFS_MD_CUDAMixedBackend::_check_stop(){
	_eval_order_parameter_states();

	bool stop_sim = false;
	bool crossing = false;

	if (_gen_flux){
		/// flux generation
		if (_flux_direction == FORWARD){
			if (_use_master_conditions){
				crossing = _test_master_conditions(_master_forward_conditions);
			}
			else{
				_eval_stop_conditions(_sc_fwd);
				crossing = _test_crossing(_sc_fwd);
			}
			if (crossing){
				// check whether we began the simulation with the system already in the astate (and the relevant flag is set in the input file)
				if (this->_curr_step == 0 && _check_initial_state){
					OX_LOG(Logger::LOG_INFO, "Began simulation out of the astate, but the check_initial_state option is set; terminating simulation");
					return true;
				}
				
				// check whether a master condition other than master_forward1 and master_backward1 was reached (and the relevant flag is set in the input file)
				if (_unexpected_master_reached && _die_on_unexpected_master){
					_handle_unexpected_master();
					return true;
				}

				if ((_gen_flux_cross_count % _gen_flux_save_every) == 0){
					char conf_str[256];
					sprintf(conf_str, "%s_N%d.dat", _conf_prefix, this->_gen_flux_saved_cross_count);
					_prepare_configuration(conf_str);
					_obs_output_custom_conf->print_output(_curr_step);
					OX_LOG(Logger::LOG_INFO, "saved configuration %s at step %d", conf_str, this->_curr_step);
					_gen_flux_saved_cross_count += 1;

					if (_gen_flux_saved_cross_count == _gen_flux_desired_cc) stop_sim = true;
				}
				else{
					OX_LOG(Logger::LOG_INFO, "Forwards condition reached at step %d", this->_curr_step);
				}
				_gen_flux_cross_count += 1;
				_flux_direction = BACKWARD;
			}
		}
		else if (_flux_direction == BACKWARD){
			if (_use_master_conditions){
				crossing = _test_master_conditions(_master_backward_conditions);
			}
			else{
				_eval_stop_conditions(_sc_bwd);
				crossing = _test_crossing(_sc_bwd);
			}
			if (crossing){
				// check whether a master condition other than master_forward1 and master_backward1 was reached (and the relevant flag is set in the input file)
				if (_unexpected_master_reached && _die_on_unexpected_master){
					_handle_unexpected_master();
					return true;
				}

				if (_gen_flux_debug) {
					char conf_str[256];
					sprintf(conf_str, "debug_backwards_N%d.dat", this->_gen_flux_cross_count);
					_prepare_configuration(conf_str);
					_obs_output_custom_conf->print_output(_curr_step);
					OX_LOG(Logger::LOG_INFO, "saved debug configuration %s at step %d", conf_str, this->_curr_step);
				}
				else{
					OX_LOG(Logger::LOG_INFO, "Backwards condition reached at step %d", this->_curr_step);
				}
				_flux_direction = FORWARD;
			}
		}
		else{
			throw oxDNAException("unknown flux direction, dying\n"); // should never happen
		}
	}
	else{
		/// shooting simulation
		if (_use_master_conditions){
			if (_test_master_conditions(_master_conditions)) stop_sim = true;
		}
		else{
			_eval_stop_conditions(_sc);
			if (_test_crossing(_sc)) stop_sim = true;
		}
	}

	return stop_sim;
}

void FFS_MD_CUDAMixedBackend::_prepare_configuration(char *conf_str){
	// set config name and set the correct energy
	_obs_output_custom_conf->change_output_file(conf_str);
	// these 3 lines to get the correct energy in the configuration file
	this->_gpu_to_host_particles();
	this->_U = GpuUtils::sum_4th_comp<float, float4>(this->_h_forces, this->_N);
	this->_K = GpuUtils::sum_4th_comp<float, float4>(this->_h_vels, this->_N) + GpuUtils::sum_4th_comp<float, float4>(this->_h_Ls, 	this->_N);
}

int *FFS_MD_CUDAMixedBackend::_get_2D_rows(int rows_len, int *lens){
	// create list of indices for the beginning of each row for a variable-width 2D array A using an array storing the lengths of each row of A
	int *rows = (int*)malloc(rows_len * sizeof(int));
	int counter=0;
	for (int ii=0 ; ii<rows_len; ii++){
		rows[ii] = counter;
		counter += lens[ii];
	}
	return rows;
}

//Used to print the order parameters whenever the energy is printed, to match CPU behaviour
char * FFS_MD_CUDAMixedBackend::get_op_state_str(void) { 
	CUDA_SAFE_CALL(cudaMemcpy(_h_hb_energies, _d_hb_energies, _n_hb_pairs * sizeof(float), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL(cudaMemcpy(_h_op_dists, _d_op_dists, _n_dist_pairs * sizeof(float), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL(cudaMemcpy(_h_nearhb_states, _d_nearhb_states, _n_hb_pairs * sizeof(bool), cudaMemcpyDeviceToHost) );

	// work out the lengths and row indices for the hbond and distance order parameter arrays
	int * h_hb_region_lens = (int *)malloc(_n_hb_regions * sizeof(int)); //These could probably be redefined as class variables, rather than being redeclared here and in init. 
	int * h_dist_region_lens = (int *)malloc(_n_dist_regions * sizeof(int));
	_op.get_hb_pairs_count(h_hb_region_lens);
	_op.get_dist_pairs_count(h_dist_region_lens);
	int * h_hb_region_rows = _get_2D_rows(_n_hb_regions, h_hb_region_lens);
	int * h_dist_region_rows = _get_2D_rows(_n_dist_regions, h_dist_region_lens);

	char * aux;
	aux = (char *) _state_str;
	for (int i = 0; i < _n_hb_regions; i++) {
		if (_op.get_hb_cutoff(i)!=64) {
			int tot_hb = 0;
			for (int j=0 ; j < h_hb_region_lens[i] ; j++){
				if (_h_hb_energies[h_hb_region_rows[i] + j] < _op.get_hb_cutoff(i)) tot_hb += 1;
			}
			sprintf(aux, "%2d ", tot_hb);
		} else {
			int tot_near_hb = 0;
			for (int j=0 ; j < h_hb_region_lens[i] ; j++) tot_near_hb+=_h_nearhb_states[h_hb_region_rows[i] + j];
			sprintf(aux, "%2d ", tot_near_hb);
			printf("Tot near is %2d \n", tot_near_hb);
		}
		aux = (char *) _state_str + strlen(_state_str);
	}
	for (int i = 0; i < _n_dist_regions; i++) {
		float min_dist = _h_op_dists[h_dist_region_rows[i]];
		for (int j=1 ; j < h_dist_region_lens[i] ; j++){
			if (_h_op_dists[h_dist_region_rows[i] + j] < min_dist) {min_dist=_h_op_dists[h_dist_region_rows[i] + j];}
		}
		sprintf(aux, "%2f ", min_dist);
		aux = (char *) _state_str + strlen(_state_str);
	}

	free(h_hb_region_lens);
	free(h_dist_region_lens);
	free(h_hb_region_rows);
	free(h_dist_region_rows);
	return _state_str;
}

//Used to print the order parameters at the end of the simulation, to match CPU behaviour
void FFS_MD_CUDAMixedBackend::sprintf_names_and_values(char *str) { 
	// copy order parameter states from GPU
	CUDA_SAFE_CALL(cudaMemcpy(_h_hb_energies, _d_hb_energies, _n_hb_pairs * sizeof(float), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL(cudaMemcpy(_h_op_dists, _d_op_dists, _n_dist_pairs * sizeof(float), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL(cudaMemcpy(_h_nearhb_states, _d_nearhb_states, _n_hb_pairs * sizeof(bool), cudaMemcpyDeviceToHost) );

	// work out the lengths and row indices for the hbond and distance order parameter arrays
	int * h_hb_region_lens = (int *)malloc(_n_hb_regions * sizeof(int));
	int * h_dist_region_lens = (int *)malloc(_n_dist_regions * sizeof(int));
	_op.get_hb_pairs_count(h_hb_region_lens);
	_op.get_dist_pairs_count(h_dist_region_lens);
	int * h_hb_region_rows = _get_2D_rows(_n_hb_regions, h_hb_region_lens);
	int * h_dist_region_rows = _get_2D_rows(_n_dist_regions, h_dist_region_lens);

	char * aux;
	aux = (char *) str;
	for (int i = 0; i < _n_hb_regions; i++) {
		if (_op.get_hb_cutoff(i)!=64) {
			int tot_hb = 0;
			for (int j=0 ; j < h_hb_region_lens[i] ; j++){
				if (_h_hb_energies[h_hb_region_rows[i] + j] < _op.get_hb_cutoff(i)) tot_hb += 1;
			}
			sprintf(aux, "%s: %d; ", _op.get_name_from_hb_id(i), tot_hb);
		} else {
			int tot_near_hb = 0;
			for (int j=0 ; j < h_hb_region_lens[i] ; j++) tot_near_hb+=_h_nearhb_states[h_hb_region_rows[i] + j];
			sprintf(aux, "%s: %d; ", _op.get_name_from_hb_id(i), tot_near_hb);
		}
		aux = (char *) str + strlen(str);
	}
	for (int i = 0; i < _n_dist_regions; i++) {
		float min_dist = _h_op_dists[h_dist_region_rows[i]];
		for (int j=1 ; j < h_dist_region_lens[i] ; j++){
			if (_h_op_dists[h_dist_region_rows[i] + j] < min_dist) {min_dist=_h_op_dists[h_dist_region_rows[i] + j];}
		}
		sprintf(aux, "%s: %12.9f; ", _op.get_name_from_distance_id(i), min_dist);
		aux = (char *) str + strlen(str);
	}

	free(h_hb_region_lens);
	free(h_dist_region_lens);
	free(h_hb_region_rows);
	free(h_dist_region_rows);
}

void FFS_MD_CUDAMixedBackend::print_observables(llint curr_step) {
	if (_curr_step % _print_energy_every == _print_energy_every-1) this->_backend_info  = get_op_state_str();
	CUDAMixedBackend::print_observables(curr_step);
}
