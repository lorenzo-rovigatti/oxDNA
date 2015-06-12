/*
 * SignalManager.h
 *
 *  Created on: 12/giu/2015
 *      Author: lorenzo
 */

#ifndef SIGNALMANAGER_H_
#define SIGNALMANAGER_H_

// see http://www.linuxjournal.com/files/linuxjournal.com/linuxjournal/articles/063/6391/6391l3.html
#include <signal.h>
#include <execinfo.h>
/* get REG_EIP from ucontext.h */
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <ucontext.h>

namespace SignalManager {

void segfault_handler(int sig, siginfo_t *info, void *secret);
void manage_segfault();

}

#endif /* SIGNALMANAGER_H_ */
