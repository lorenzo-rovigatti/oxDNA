/*
 * BackendFactory.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: rovigatti
 */

#include "BackendFactory.h"

#include "MD_CPUBackend.h"
#include "FFS_MD_CPUBackend.h"
#include "MC_CPUBackend.h"
#include "MC_CPUBackend2.h"
#include "FFS_MD_CPUBackend.h"
#include "VMMC_CPUBackend.h"
#ifndef NOCUDA
#include "../CUDA/Backends/MD_CUDABackend.h"
#include "../CUDA/Backends/MD_CUDAMixedBackend.h"
#include "../CUDA/Backends/FFS_MD_CUDAMixedBackend.h"
#endif

#ifdef HAVE_MPI
#include "MD_MPIBackend.h"
#include "PT_VMMC_CPUBackend.h"
#endif

BackendFactory::BackendFactory() {

}

BackendFactory::~BackendFactory() {

}

ISimBackend *BackendFactory::make_backend(input_file &inp) {
	char backend_opt[256], backend_prec[256], sim_type[256];
	getInputString(&inp, "backend", backend_opt, 1);

	if(getInputString(&inp, "backend_precision", backend_prec, 0) == KEY_NOT_FOUND) {
		OX_LOG(Logger::LOG_INFO, "Backend precision not speficied, using double");
		sprintf(backend_prec, "%s", "double");
	}
	else OX_LOG(Logger::LOG_INFO, "Backend precision: %s", backend_prec);

	if(getInputString(&inp, "sim_type", sim_type, 0) == KEY_NOT_FOUND) {
		OX_LOG(Logger::LOG_INFO, "Simulation type not specified, using MD");
		sprintf(sim_type, "%s", "MD");
	}
	else OX_LOG(Logger::LOG_INFO, "Simulation type: %s", sim_type);

	ISimBackend *new_backend = NULL;
	if(!strcmp(sim_type, "MD")) {
		if(!strcmp(backend_opt, "CPU")) {
#ifndef HAVE_MPI
			if(!strcmp(backend_prec, "double")) new_backend = new MD_CPUBackend<double>();
			else if(!strcmp(backend_prec, "float")) new_backend = new MD_CPUBackend<float>();
			else throw oxDNAException("Backend precision '%s' is not supported", backend_prec);
#endif
#ifdef HAVE_MPI
			if(!strcmp(backend_prec, "double")) new_backend = new MD_MPIBackend<double>();
			else if(!strcmp(backend_prec, "float")) new_backend = new MD_MPIBackend<float>();
			else throw oxDNAException("Backend precision '%s' is not supported", backend_prec);
			int procsize;
			MPI_Comm_size (MPI_COMM_WORLD, &procsize);
			OX_LOG(Logger::LOG_INFO, "MPI Simulation with %d nodes", procsize);
#endif

		}
#ifndef NOCUDA
		else if(!strcmp(backend_opt, "CUDA")) {
			// LR_double4 is defined in CUDAUtils.h
			if(!strcmp(backend_prec, "double"))  new_backend = new MD_CUDABackend<double, LR_double4>();
			else if(!strcmp(backend_prec, "float")) new_backend = new MD_CUDABackend<float, float4>();
			else if(!strcmp(backend_prec, "mixed")) new_backend = new CUDAMixedBackend();
			else throw oxDNAException("Backend precision '%s' is not supported", backend_prec);
		}
#endif
		else throw oxDNAException("Backend '%s' not supported", backend_opt);
	}
	else if(!strcmp(sim_type, "MC")) {
		if(!strcmp(backend_opt, "CPU")) {
			if(!strcmp(backend_prec, "double")) new_backend = new MC_CPUBackend<double>();
			else if(!strcmp(backend_prec, "float")) new_backend = new MC_CPUBackend<float>();
			else throw oxDNAException("Backend precision '%s' is not supported", backend_prec);
		}
		else throw oxDNAException("Backend '%s' not supported", backend_opt);
	}
	else if(!strcmp(sim_type, "MC2")) {
		if(!strcmp(backend_opt, "CPU")) {
			if(!strcmp(backend_prec, "double")) new_backend = new MC_CPUBackend2<double>();
			else if(!strcmp(backend_prec, "float")) new_backend = new MC_CPUBackend2<float>();
			else throw oxDNAException("Backend precision '%s' is not supported", backend_prec);
		}
		else throw oxDNAException("Backend '%s' not supported", backend_opt);
	}
	else if(!strcmp (sim_type, "VMMC")) {
			if(!strcmp(backend_opt, "CPU")) {
				if(!strcmp(backend_prec, "double")) new_backend = new VMMC_CPUBackend<double>();
				else if(!strcmp(backend_prec, "float")) {
					new_backend = new VMMC_CPUBackend<float>();
				}
				else throw oxDNAException("Backend precision '%s' is not supported", backend_prec);
			}
			else throw oxDNAException("Backend '%s' not supported", backend_opt);

	}
	else if(!strcmp (sim_type, "FFS_MD")) {
		if(!strcmp(backend_opt, "CPU")) {
			if(!strcmp(backend_prec, "double")) new_backend = new FFS_MD_CPUBackend<double>();
			else if(!strcmp(backend_prec, "float")) new_backend = new FFS_MD_CPUBackend<float>();
			else throw oxDNAException("Backend precision '%s' is not supported", backend_prec);
		}
#ifndef NOCUDA
		else if(!strcmp(backend_opt, "CUDA")) {
			if(!strcmp(backend_prec, "mixed")) new_backend = new FFS_MD_CUDAMixedBackend();
			else throw oxDNAException("Backend precision '%s' is not supported", backend_prec);
		}
#endif
		else throw oxDNAException("Backend '%s' not supported", backend_opt);
	}
#ifdef HAVE_MPI
	else if(!strcmp (sim_type, "PT_VMMC")) {
		if(!strcmp(backend_opt, "CPU")) {
			if(!strcmp(backend_prec, "double")) new_backend = new PT_VMMC_CPUBackend<double>();
			else if(!strcmp(backend_prec, "float")) new_backend = new PT_VMMC_CPUBackend<float>();
			else throw oxDNAException("Backend precision '%s' is not supported", backend_prec);
		}
		else throw oxDNAException("Backend '%s' not supported", backend_opt);
	}
#endif
	else throw oxDNAException("Simulation type '%s' not supported", sim_type);

	return new_backend;
}
