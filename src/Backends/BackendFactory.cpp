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
	char backend_opt[256], backend_prec[256], sim_type[256];
	getInputString(&inp, "backend", backend_opt, 1);

	if(getInputString(&inp, "backend_precision", backend_prec, 0) == KEY_NOT_FOUND) {
		OX_LOG(Logger::LOG_INFO, "Backend precision not specified, using double");
		sprintf(backend_prec, "%s", "double");
	}
	else if(!strcmp(backend_opt, "CPU")) {
		OX_LOG(Logger::LOG_NOTHING, "");
		OX_LOG(Logger::LOG_WARNING, "The 'backend_precision' option cannot be set by input file when running on CPU\n");
	}
	else {
		OX_LOG(Logger::LOG_INFO, "Backend precision: %s", backend_prec);
	}

	if(getInputString(&inp, "sim_type", sim_type, 0) == KEY_NOT_FOUND) {
		OX_LOG(Logger::LOG_INFO, "Simulation type not specified, using MD");
		sprintf(sim_type, "%s", "MD");
	}
	else OX_LOG(Logger::LOG_INFO, "Simulation type: %s", sim_type);

	SimBackend *new_backend = NULL;
	if(!strcmp(sim_type, "MD")) {
		if(!strcmp(backend_opt, "CPU")) {
			new_backend = new MD_CPUBackend();
		}
#ifndef NOCUDA
		else if(!strcmp(backend_opt, "CUDA")) {
#ifndef CUDA_DOUBLE_PRECISION
			if(!strcmp(backend_prec, "mixed")) new_backend = new CUDAMixedBackend();
			else if(!strcmp(backend_prec, "float")) new_backend = new MD_CUDABackend();
			else throw oxDNAException("Backend precision '%s' is not allowed, as the code has been compiled with 'float' and 'mixed' support only", backend_prec);
#else
			if(!strcmp(backend_prec, "double")) new_backend = new MD_CUDABackend();
			else throw oxDNAException("Backend precision '%s' is not allowed, as the code has been compiled with 'double' support only", backend_prec);
#endif
		}
#endif
		else throw oxDNAException("Backend '%s' not supported", backend_opt);
	}
	else if(!strcmp(sim_type, "MC")) {
		if(!strcmp(backend_opt, "CPU")) {
			new_backend = new MC_CPUBackend();
		}
		else throw oxDNAException("Backend '%s' not supported", backend_opt);
	}
	else if(!strcmp(sim_type, "MC2")) {
		if(!strcmp(backend_opt, "CPU")) {
			new_backend = new MC_CPUBackend2();
		}
		else throw oxDNAException("Backend '%s' not supported", backend_opt);
	}
	else if(!strcmp (sim_type, "VMMC")) {
			if(!strcmp(backend_opt, "CPU")) {
				new_backend = new VMMC_CPUBackend();
			}
			else throw oxDNAException("Backend '%s' not supported", backend_opt);

	}
#ifdef HAVE_MPI
	else if(!strcmp (sim_type, "PT_VMMC")) {
			if(!strcmp(backend_opt, "CPU")) {
				if(!strcmp(backend_prec, "double")) new_backend = new PT_VMMC_CPUBackend<double>();
				else if(!strcmp(backend_prec, "float")) {
					new_backend = new PT_VMMC_CPUBackend<float>();
				}
				else throw oxDNAException("Backend precision '%s' is not supported", backend_prec);
			}
			else throw oxDNAException("Backend '%s' not supported", backend_opt);

	}
#endif
	else if(!strcmp (sim_type, "min")) {
			if(!strcmp(backend_opt, "CPU")) {
				new_backend = new MinBackend();
			}
			else throw oxDNAException("Backend '%s' not supported with sim_type = %s", backend_opt, sim_type);
	}
	else if(!strcmp (sim_type, "FIRE")) {
			if(!strcmp(backend_opt, "CPU")) {
				new_backend = new FIREBackend();
			}
			else throw oxDNAException("Backend '%s' not supported with sim_type = %s", backend_opt, sim_type);
	}
	else if(!strcmp (sim_type, "FFS_MD")) {
		if(!strcmp(backend_opt, "CPU")) {
			new_backend = new FFS_MD_CPUBackend();
		}
#ifndef NOCUDA
		else if(!strcmp(backend_opt, "CUDA")) {
			if(!strcmp(backend_prec, "mixed")) new_backend = new FFS_MD_CUDAMixedBackend();
			else throw oxDNAException("Backend precision '%s' for FFS simulations with CUDA is not supported", backend_prec);
		}
#endif
		else throw oxDNAException("Backend '%s' not supported", backend_opt);
	}
	else throw oxDNAException("Simulation type '%s' not supported", sim_type);

	return std::shared_ptr<SimBackend>(new_backend);
}
