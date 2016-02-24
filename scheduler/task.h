/*---------------------------------------------------------------------------
* task.h: Class representing a task to be run.
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 10:57:03
* Last Modified by:   amedhi
* Last Modified time: 2015-09-28 22:50:08
*----------------------------------------------------------------------------*/
#include <iostream>
#include "inputparams.h"

#ifndef SCHEDULER_TASK_H
#define SCHEDULER_TASK_H

namespace scheduler {

class AbstractTask
{
public:
  virtual ~AbstractTask() {};
  virtual void start(input::Parameters& p) = 0; // start all runs
  virtual void run() = 0; // run for some time (in seconds)
  virtual void halt() = 0; // halt all runs, simulation is finished        
};

class Task : public AbstractTask
{
public:
  Task();
  virtual ~Task();
  
  //virtual void construct(); // needs to be called to finish construction
  void start(input::Parameters& p) override; // start simulation
  void run() override;// run a few steps and return control
  virtual void dostep()=0; // do a step
  //void finish(); // mark as finished
  // bool started() const { return started_;}
  void halt() override;
  
protected:
  bool finished_;

private:
  bool started_; // is the task running?
};


} // end namespace input

#endif
