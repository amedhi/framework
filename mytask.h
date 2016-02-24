/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 12:08:30
* Last Modified by:   amedhi
* Last Modified time: 2016-02-22 13:37:29
*----------------------------------------------------------------------------*/
#include <iostream>
#include "scheduler/task.h"
#include "lattice/lattice.h"

/*---------------MYTask object----------------------*/
class MYTask : public scheduler::Task
{
public:
  MYTask() {};
  ~MYTask() {};
  void start(input::Parameters& parms); 
  void run(); 
  void dostep(); 
  void halt(); 
  void run_test_code(const input::Parameters& parms); 
private:
  int nsites;
} mytask;

void MYTask::start(input::Parameters& parms)
{

  int i = parms.set_value("nsite", 0);
  double beta = parms.set_value("beta", 0.0);
  std::cout << "i=" << i << std::endl;
  std::cout << "b=" << beta << std::endl;


  //double u = parms.set_value("U", 1.0);
  std::cout << "task " << parms.task_id()+1 << " of " << parms.task_size() << std::endl;

  std::cout << "\nChecking latticelibrary" << std::endl;
  lattice::Lattice lattice(parms); 



}

void MYTask::run()
{
  //double u = parms.set_value("U", 1.0);
}

void MYTask::dostep()
{
}

void MYTask::halt()
{
}


