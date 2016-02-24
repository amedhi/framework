/*---------------------------------------------------------------------------
* master_scheduler: Scheduler running on master
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 19:51:04
* Last Modified by:   amedhi
* Last Modified time: 2015-09-28 22:50:02
*----------------------------------------------------------------------------*/
#include "scheduler.h"

namespace scheduler {

MasterScheduler::MasterScheduler(int argc, const char *argv[])
  : Scheduler(), cmdarg(argc, argv), input(), task_size{0}
{
}

int MasterScheduler::run(AbstractTask& theTask) 
{
  if (cmdarg.not_valid()) return 0;
  std::cout << " starting..." << std::endl;

  // job parameters
  read_parameters(cmdarg.filename());
  if (!valid) return -1;

  for (unsigned task=0; task<task_size; ++task) {
    set_task_parameters(task);
    theTask.start(parms);
    //params << pstore(task_id);
  }

  return 0;
}

void MasterScheduler::read_parameters(const std::string& filename) 
{
  valid = input.read_params(filename);
  if (valid) {
    input.init_task_param(parms);
    task_size = input.task_size();
  }
}

void MasterScheduler::set_task_parameters(const unsigned& task) 
{
  input.set_task_param(parms, task);
  //parms.show(task);
}


} // end namespace scheduler
