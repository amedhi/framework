/*---------------------------------------------------------------------------
* Scheduler: A class that handles user jobs.  
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-08-17 13:33:19
* Last Modified by:   amedhi
* Last Modified time: 2015-09-28 22:48:00
*----------------------------------------------------------------------------*/
// File: scheduler.h 
// Definition of the Scheduler class.

#ifndef SCHEDULER_SCHEDULER_H
#define SCHEDULER_SCHEDULER_H

#include <iostream>
#include "task.h"
#include "cmdargs.h"
#include "inputparams.h"

namespace scheduler {

int start(int argc, const char *argv[], AbstractTask& theTask);

class Scheduler 
{
public:
  Scheduler(): simmaster(0) {};
  ~Scheduler() {};
  virtual int run(AbstractTask& theTask);

protected:
  input::Parameters parms;
  bool valid;

private:
  unsigned simmaster;
  //TaskParams parms;
};

class MasterScheduler : public Scheduler
{
public:
  MasterScheduler(int argc, const char *argv[]);
  MasterScheduler() = delete;
  ~MasterScheduler() {};
  void read_parameters(const std::string& filename);
  void set_task_parameters(const unsigned& task);
  int run(AbstractTask& theTask) override;

private:
  CmdArg cmdarg;
  input::JobInput input;
  unsigned int task_size;
};

} // end namespace input

#endif
