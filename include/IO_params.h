#ifndef IO_params_H_
#define IO_params_H_

#include <array>
#include <cstring> 
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>

namespace cluster {

struct LatticeDataContainer { // Just the thing that holds all variables
  // lattice parameter
  int L[4];
  int V;
  // action parameter
  std::string formulation;
  double kappa;
  double lambda;
  // metropolis parameter
  int metropolis_local_hits;
  int metropolis_global_hits;
  double metropolis_delta;
  // cluster parameter
  int cluster_hits;
  double cluster_min_size;
  // run parameter
  int start_measure;
  int total_measure;
  int measure_every_X_updates;
  std::string outpath;
};
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class IO_params {

private:

  inline LatticeDataContainer read_infile(int argc, char** argv) {

    int opt = -1;
    int reader = 0;
    char infilename[200];
    char readin[256];
    FILE* infile = NULL;

    LatticeDataContainer data;

    // search for command line option and put filename in "infilename"
    for(int i = 0; i < argc; ++i) {
      if(std::strcmp(argv[i], "-i") == 0) {
        opt = i+1;
        break;
      }
    }
    if(opt < 0) {
      mdp << "No input file specified, trying infile.in" << endl;
      sprintf(infilename, "infile.in");
    } else {
      sprintf(infilename, "%s", argv[opt]);
      mdp << "Trying input file " << infilename << endl;
    }
    // open file for reading
    if ((infile = fopen(infilename, "r")) == NULL ) {
      std::cerr << "Could not open file " << infilename << endl;
      std::cerr << "Aborting..." << endl;
      exit(-10);
    }
    // lattice size
    reader += fscanf(infile, "L = %d %d %d %d \n", &data.L[0], &data.L[1], 
                                                   &data.L[2], &data.L[3]);
    data.V = data.L[0]*data.L[1]*data.L[2]*data.L[3];
    // kappa and lambda
    reader += fscanf(infile, "formulation = %255s\n", readin);
    data.formulation.assign(readin);
    reader += fscanf(infile, "kappa = %lf\n", &data.kappa);
    reader += fscanf(infile, "lambda = %lf\n", &data.lambda);

    if(data.formulation == "continuum"){
      data.lambda = 4.*data.kappa*data.kappa*data.lambda;
      mdp << "Parameters lambda and kappa are changed to lattice versions: \n"
          << "\tlambda = " << data.lambda << " kappa = " << data.kappa << endl; 
    }
    // metropolis 
    reader += fscanf(infile, "metropolis_local_hits = %d\n", 
                             &data.metropolis_local_hits);
    reader += fscanf(infile, "metropolis_global_hits = %d\n", 
                             &data.metropolis_global_hits);
    reader += fscanf(infile, "metropolis_delta = %lf\n", &data.metropolis_delta);
    // cluster
    reader += fscanf(infile, "cluster_hits = %d\n", &data.cluster_hits);
    reader += fscanf(infile, "cluster_min_size = %lf\n", &data.cluster_min_size);
    // configs
    reader += fscanf(infile, "start_measure = %d\n", &data.start_measure);
    reader += fscanf(infile, "total_measure = %d\n", &data.total_measure);
    reader += fscanf(infile, "measure_every_X_updates = %d\n", 
                             &data.measure_every_X_updates);
    reader += fscanf(infile, "outpath = %255s\n", readin);
    data.outpath.assign(readin);

    // close input file
    fclose(infile);

    return data;
  };

public:
  const LatticeDataContainer data; 
  
  IO_params(int argc, char** argv) : data(read_infile(argc, argv)) {};

}; // end of class definition

} // end of namespace

#endif // IO_params
