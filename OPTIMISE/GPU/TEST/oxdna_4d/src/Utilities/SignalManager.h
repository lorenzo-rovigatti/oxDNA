/*
 * SignalManager.h
 *
 *  Created on: 12/giu/2015
 *      Author: lorenzo
 */

#ifndef SIGNALMANAGER_H_
#define SIGNALMANAGER_H_

// see http://www.linuxjournal.com/files/linuxjournal.com/linuxjournal/articles/063/6391/6391l3.html
#include <csignal>
#include <execinfo.h>

#ifdef SIGNAL
/* get REG_EIP from ucontext.h */
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <ucontext.h>

#if __WORDSIZE == 64
#define MY_REG_RIP REG_RIP
#else
#define MY_REG_RIP REG_EIP
#endif
#endif



namespace SignalManager {

void segfault_handler(int sig, siginfo_t *info, void *secret);
void manage_segfault();


}

#endif /* SIGNALMANAGER_H_ */
