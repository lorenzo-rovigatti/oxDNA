/*
 * SignalManager.cpp
 *
 *  Created on: 12/giu/2015
 *      Author: lorenzo
 */

#include "SignalManager.h"

#ifdef SIGNAL
#include "../defs.h"
#include "oxDNAException.h"

void SignalManager::segfault_handler(int sig, siginfo_t *info, void *secret) {
	if(sig != SIGSEGV) throw oxDNAException("segfault_handler should handle segmentation faults only, got %d", sig);

	void *trace[16];
	char **messages = (char **) NULL;
	int i, trace_size = 0;
	ucontext_t *uc = (ucontext_t *) secret;

	OX_LOG(Logger::LOG_ERROR, "Segmentation fault identified, faulty address is %p, from %p", info->si_addr, uc->uc_mcontext.gregs[MY_REG_RIP]);

	trace_size = backtrace(trace, 16);
	/* overwrite sigaction with caller's address */
	trace[1] = (void *) uc->uc_mcontext.gregs[MY_REG_RIP];

	messages = backtrace_symbols(trace, trace_size);
	/* skip first stack frame (points here) */

	OX_LOG(Logger::LOG_INFO, "Execution path:");
	for (i = 1; i < trace_size; ++i) OX_LOG(Logger::LOG_NOTHING, "\t[bt] %s", messages[i]);

	exit(0);
}

void SignalManager::manage_segfault() {
	struct sigaction sa;

	sa.sa_sigaction = segfault_handler;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = SA_RESTART | SA_SIGINFO;

	sigaction(SIGSEGV, &sa, NULL);
}

#else
void SignalManager::segfault_handler(int sig, siginfo_t *info, void *secret) {}
void SignalManager::manage_segfault() {}

#endif
