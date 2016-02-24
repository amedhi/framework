/*---------------------------------------------------------------------------
* InputParameter: Classes for processing input parameters.
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-08-17 13:33:19
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-01-25 18:43:46
*----------------------------------------------------------------------------*/
// File: inputparams.h 

#ifndef INPUT_PARAMETERS_H
#define INPUT_PARAMETERS_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>

namespace input {

enum class value_type {boo, num, str, nan};

class Parameters;  // forward declaration

class JobInput
{
public:
  JobInput() {n_tasks=n_params=0; valid=false;} 
  JobInput(const std::string& inputfile); 
  ~JobInput() {} 
  bool read_params(const std::string& inputfile);
  bool not_valid(void) const {return !valid;};
  unsigned int task_size(void) {return n_tasks;}
  void get_task_param(const unsigned& task_id); 
  void init_task_param(Parameters& p);
  void set_task_param(Parameters& p, const unsigned& task_id); 

private:
  struct parameter {std::string name; value_type type; unsigned size;};
  unsigned int n_params;
  unsigned int n_tasks;
  bool valid;
  std::string infile;
  std::vector<parameter> param_list;
  std::map<std::string, std::vector<bool> > boo_params;
  std::map<std::string, std::vector<double> > num_params;
  std::map<std::string, std::vector<std::string> > str_params;

  unsigned int parse(const std::string& inputfile);

  // bad_input exception
  class bad_input: public std::runtime_error
  {
  public:
    bad_input(const std::string& msg, const int& ln=-1);
    std::string message(void) const; 
  private:
    int lnum;
  };
};


class Parameters 
{
public:
  Parameters() {};
  ~Parameters() {};
  int set_value(const std::string& pname, const int& defval) const;  
  int set_value(const std::string& pname, const int& defval, int& info) const;
  int set_value(const std::string& pname, const int& defval, int& info, bool& is_const) const;
  double set_value(const std::string& pname, const double& defval) const;  
  double set_value(const std::string& pname, const double& defval, int& info) const;  
  double set_value(const std::string& pname, const double& defval, int& info, bool& is_const) const;  
  std::string set_value(const std::string& pname, const std::string& defval) const;  
  std::string set_value(const std::string& pname, const std::string& defval, int& info) const;  
  bool is_constant(const std::string& pname) const;  
  unsigned task_id(void) const { return this_task; }
  unsigned task_size(void) const { return n_tasks; }
  void show(const unsigned&) const;
  friend void JobInput::init_task_param(Parameters& p);
  friend void JobInput::set_task_param(Parameters& p, const unsigned& task_id);

private:
  struct pval {bool is_const; value_type type; bool bool_val; double num_val; std::string str_val;};
  std::map<std::string, pval> params;
  unsigned n_params;
  unsigned this_task;
  unsigned n_tasks;
  mutable std::map<std::string, pval>::const_iterator it;
  void warn_not_found(const std::string& pname) const;
  void warn_type_mismatch(const std::string& pname, const std::string& type) const;
};



} // end namespace input

#endif
