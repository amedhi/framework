/*---------------------------------------------------------------------------
* Scheduler: A class that handles user jobs.  
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-08-17 12:44:04
* Last Modified by:   amedhi
* Last Modified time: 2015-09-30 11:59:07
*----------------------------------------------------------------------------*/
// File: scheduler.cc 
// Implementation of the Scheduler class.

#include <iostream>
#include "scheduler.h"

namespace scheduler {

int start(int argc, const char *argv[], AbstractTask& theTask)
{

  Scheduler* theScheduler;
  theScheduler = new MasterScheduler(argc, argv);
  int res = theScheduler->run(theTask);

	return res;	
}

int Scheduler::run(AbstractTask& theTask) 
{
  // Task task;
  // valid = input.read_params(0);
  // task.init_task_param()
  return 0;
}

} // end namespace scheduler
