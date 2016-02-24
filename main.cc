#include <iostream>
#include "scheduler/scheduler.h"
#include "mytask.h"

int main(int argc, const char *argv[])
{
  try {
	 return scheduler::start(argc, argv, mytask);
  }
  catch (std::exception& exn) {
    std::cout << exn.what() << std::endl;
    return -1;
  }
  catch (...) {
    std::cout << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
}
