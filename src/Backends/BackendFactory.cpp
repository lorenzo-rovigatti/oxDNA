/*
 * BackendFactory.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: rovigatti
 */

#include "BackendFactory.h"

#include "MD_CPUBackend.h"
#include "MC_CPUBackend.h"
#include "MC_CPUBackend2.h"
#include "VMMC_CPUBackend.h"
#include "FFS_MD_CPUBackend.h"
#include "FFS_MD_CPUBackend.h"
#include "MinBackend.h"
#include "FIREBackend.h"
#ifndef NOCUDA
#include "../CUDA/Backends/MD_CUDABackend.h"
#ifndef CUDA_DOUBLE_PRECISION
#include "../CUDA/Backends/MD_CUDAMixedBackend.h"
#include "../CUDA/Backends/FFS_MD_CUDAMixedBackend.h"
#endif
#endif

#ifdef HAVE_MPI
#include "PT_VMMC_CPUBackend.h"
#endif

std::shared_ptr<SimBackend> BackendFactory::make_backend(input_file &inp) {
	std::string backend_opt, backend_prec, sim_type;
	getInputString(&inp, "backend", backend_opt, 1);

	int precision_state = getInputString(&inp, "backend_precision", backend_prec, 0);

	if(precision_state == KEY_FOUND && backend_opt == "CPU") {
		OX_LOG(Logger::LOG_NOTHING, "");
		OX_LOG(Logger::LOG_WARNING, "The 'backend_precision' option cannot be set by input file when running on CPU\n");
	}

	if(getInputString(&inp, "sim_type", sim_type, 0) == KEY_NOT_FOUND) {
		OX_LOG(Logger::LOG_INFO, "Simulation type not specified, using MD");
		sim_type = "MD";
	}
	else {
		OX_LOG(Logger::LOG_INFO, "Simulation type: %s", sim_type.c_str());
	}

	SimBackend *new_backend = nullptr;
	if(sim_type == "MD") {
		if(backend_opt == "CPU") {
			new_backend = new MD_CPUBackend();
		}
#ifndef NOCUDA
		else if(backend_opt == "CUDA") {
#ifndef CUDA_DOUBLE_PRECISION
			if(precision_state == KEY_NOT_FOUND) {
				backend_prec = "mixed";
			}
			if(backend_prec == "mixed") {
				new_backend = new CUDAMixedBackend();
			}
			else if(backend_prec == "float") {
				new_backend = new MD_CUDABackend();
			}
			else {
				throw oxDNAException("Backend precision '%s' is not allowed, as the code has been compiled with 'float' and 'mixed' support only", backend_prec.c_str());
			}
#else
			if(precision_state == KEY_NOT_FOUND) {
				backend_prec = "double";
			}
			if(backend_prec == "double") {
				new_backend = new MD_CUDABackend();
			}
			else {
				throw oxDNAException("Backend precision '%s' is not allowed, as the code has been compiled with 'double' support only", backend_prec.c_str());
			}
#endif
			OX_LOG(Logger::LOG_INFO, "CUDA backend precision: %s", backend_prec.c_str());
		}
#endif
		else {
			throw oxDNAException("Backend '%s' not supported", backend_opt.c_str());
		}
	}
	else if(sim_type == "MC") {
		if(backend_opt == "CPU") {
			new_backend = new MC_CPUBackend();
		}
		else {
			throw oxDNAException("Backend '%s' not supported", backend_opt.c_str());
		}
	}
	else if(sim_type == "MC2") {
		if(backend_opt == "CPU") {
			new_backend = new MC_CPUBackend2();
		}
		else {
			throw oxDNAException("Backend '%s' not supported", backend_opt.c_str());
		}
	}
	else if(sim_type == "VMMC") {
			if(backend_opt == "CPU") {
				new_backend = new VMMC_CPUBackend();
			}
			else {
				throw oxDNAException("Backend '%s' not supported", backend_opt.c_str());
			}
	}
#ifdef HAVE_MPI
	else if(sim_type == "PT_VMMC") {
			if(backend_opt == "CPU") {
				new_backend = new PT_VMMC_CPUBackend();
			}
			else {
				throw oxDNAException("Backend '%s' not supported", backend_opt.c_str());
			}

	}
#endif
	else if(sim_type == "min") {
			if(backend_opt == "CPU") {
				new_backend = new MinBackend();
			}
			else {
				throw oxDNAException("Backend '%s' not supported with sim_type = %s", backend_opt.c_str(), sim_type.c_str());
			}
	}
	else if(sim_type == "FIRE") {
			if(backend_opt == "CPU") {
				new_backend = new FIREBackend();
			}
			else {
				throw oxDNAException("Backend '%s' not supported with sim_type = %s", backend_opt.c_str(), sim_type.c_str());
			}
	}
	else if(sim_type == "FFS_MD") {
		if(backend_opt == "CPU") {
			new_backend = new FFS_MD_CPUBackend();
		}
#ifndef NOCUDA
#ifndef CUDA_DOUBLE_PRECISION
		else if(backend_opt == "CUDA") {
			if(backend_prec == "mixed") {
				new_backend = new FFS_MD_CUDAMixedBackend();
			}
			else {
				throw oxDNAException("Backend precision '%s' for FFS simulations with CUDA is not supported", backend_prec.c_str());
			}
		}
#endif
#endif
		else {
			throw oxDNAException("Backend '%s' not supported", backend_opt.c_str());
		}
	}
	else {
		throw oxDNAException("Simulation type '%s' not supported", sim_type.c_str());
	}

	return std::shared_ptr<SimBackend>(new_backend);
}
