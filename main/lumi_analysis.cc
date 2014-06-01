#include <iostream>
#include <string>

#include "analysis.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  Analysis analysis("params/parameters.json");
  int err = analysis.RunAnalysis();
  if (!err) {
    cout << "Execution successful" << endl;
  } else {
    cout << "Execution failed" << endl;
  }
}
