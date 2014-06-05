#include <iostream>
#include <string>

#include "analysis.h"

using std::string;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char *argv[]) {
  string params_filepath;
  if (argc > 2) {
    cerr << "ERROR: invalid number of command line args. This program takes"
         << "one optional argument for the filepath of the parameters file."
         << endl;
  } else if (argc == 2) {
    params_filepath = argv[1];
  } else {
    params_filepath = "params/parameters.json";
  }

  Analysis analysis(params_filepath);
  int err = analysis.RunAnalysis();
  if (!err) {
    cout << "Execution successful" << endl;
  } else {
    cout << "Execution failed" << endl;
  }
}
