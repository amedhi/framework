/*---------------------------------------------------------------------------
* task.cc: Implementation of the task classes.  
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 11:01:37
* Last Modified by:   amedhi
* Last Modified time: 2015-09-28 22:02:47
*----------------------------------------------------------------------------*/
#include "task.h"

namespace scheduler {

Task::Task() : finished_(false), started_(false)
{
}

Task::~Task()
{
}

void Task::start(input::Parameters& p)
{
  started_ = true;
}

void Task::run()
{
  //if(started() && !finished_)
    dostep();
}

void Task::halt()
{
  started_=false;
}

/*
void Task::finish()
{
  finished_=true;
}

// halt all active runs
*/


} // end namespace input
