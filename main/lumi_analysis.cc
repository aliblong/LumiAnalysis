#include <iostream>
#include <string>

#include "analysis.h"
#include "error.h"

using std::string;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char *argv[]) {
  // Accepts a command line argument for the parameters json file path.
  // Otherwise, defaults to params/parameters.json
  string params_filepath;
  if (argc > 2) {
    cerr << "ERROR: invalid number of command line args. This program takes"
         << "one optional argument for the filepath of the parameters file."
         << endl;
    return 1;
  } else if (argc == 2) {
    params_filepath = argv[1];
  } else {
    params_filepath = "params/param_trees/default.json";
  }

  Analysis analysis(params_filepath);
  LOG_IF_ERR( analysis.RunAnalysis() )
}
