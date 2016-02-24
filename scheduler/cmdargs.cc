/*---------------------------------------------------------------------------
* DMRG Project: DMRG using Matrix Product States 
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-08-17 13:33:19
* Last Modified by:   amedhi
* Last Modified time: 2015-09-28 20:25:50
*----------------------------------------------------------------------------*/
// File: cmdopts.cpp 
// Class definitions for handling command line arguments

#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include "cmdargs.h"
#include "optionparser.h"

namespace scheduler {

// constructor
CmdArg::CmdArg(int argc, const char *argv[]): progname("program"), inputfile("")
{
	if (argc > 0) {
		boost::filesystem::path p(argv[0]);
		progname = p.filename().string();
		--argc; ++argv; // skip program name for the rest
	}
	// option stats
  option::Stats stats(usage, argc, argv);
  options = new option::Option[stats.options_max];
  buffer = new option::Option[stats.buffer_max];
  // parse the arguments
  option::Parser parse(usage, argc, argv, options, buffer);
  option_count = parse.optionsCount();
  nonOption_count = parse.nonOptionsCount();
  if (nonOption_count) inputfile = parse.nonOption(0);
  valid = process_options();
  if (!valid) inputfile.clear();
}

// default desctructor
CmdArg::CmdArg(): progname("program"), inputfile("input.parm")
{
  options = new option::Option[1];
  buffer = new option::Option[1];
  option_count = 0;
}

// desctructor
CmdArg::~CmdArg()
{
	delete[] options;
	delete[] buffer;
}

// member functions
bool CmdArg::process_options(void) const
{
	// handle a few of the options
	if (options[unknown]) {
		std::cout << "Unknown option: " << options[unknown].name << "\n\n";
		option::printUsage(std::cout, usage);
		return false;
	}

	if (options[help]) {
		option::printUsage(std::cout, usage);
		return false;
	}

  if (options[info]) {
    std::cout << "Build info:" << std::endl;
    return false;
  }

  //option::Option* opt = options[file];
  //if (options[file]) inputfile = options[file].arg;

	return true;
}


} // end namespace
