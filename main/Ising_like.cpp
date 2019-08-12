#include <array>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>

#include "mdp.h"

#include "Ising_params.h" 

// the random number generator
mdp_random_generator random1;
const int comp = 2;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
typedef enum cluster_state_t {
  CLUSTER_UNCHECKED=0,
  CLUSTER_FLIP 
} cluster_state_t;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline double compute_magnetisation(mdp_field<std::array<double,comp> >& phi, 
                                    mdp_site& x){

  double m[comp];
  std::fill_n(m, comp, 0.0);
  forallsites(x){
    for(int i = 0; i < comp ; i++) m[i] += phi(x)[i];
  }
  double result = 0.0;
  for(int i = 0; i < comp ; i++) result += m[i]*m[i];
  return sqrt(result);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline std::array<double,comp> create_phi_update(const double delta){
  
  double update[comp];
  for(int i = 0; i < comp ; i++) update[i] = random1.plain()*2. - 1.)*delta;
  return update;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*
inline double metropolis_hits(mdp_field<std::array<double, 2> >& phi ,mdp_site& x,
                              const double kappa, const double delta, 
                              const size_t nb_of_hits){
    for(int parity = EVEN, parity <= ODD, parity++){
        forallsitesofparity(x,parity){
          // running over all components
          
        
      }    
    phi.update(parity); // communicate boundaries
  }
                              
}
  
*/                            

int main(int argc, char** argv){

  mdp.open_wormholes(argc,argv);

  cluster::IO_params params(argc, argv); // reading infile

  // lattice parameters
  const int dim = params.data.dimension ;
  int L[dim] ;
  std::fill_n(L, dim, params.data.L); 
  const int T = params.data.T; // useless at this point
  const int V = pow(params.data.L,dim);
  //const int comp = params.data.comp;
  
  // setup the lattice and filds
  mdp_lattice hypercube(dim,L); // declare lattice
  mdp_field<std::array<double, comp> > phi(hypercube); // declare phi field
  mdp_site x(hypercube); // declare lattice lookuptable
  
  mdp << "\ndimension = " << dim << endl;
  mdp << "L = " << L[0] << "\t" << L[1] << endl;
  mdp << "V = " << V <<endl;
  mdp << "componenets = " << comp << endl;
  mdp << "kappa = " << params.data.kappa << endl;

  random1.initialize(1227);

  /*
  // random start configuration
  forallsites(x)
    phi(x) = create_phi_update(1.,comp); 

  // compute magnetisation on start config
  //rotate_phi_field(phi, x, double(V));
  //double M = compute_magnetisation(phi, x);
  //mdp.add(M);
  //mdp << "\n\n\tmagnetization at start = " << M/V << endl;

  
  std::string mag_file = params.data.outpath + 
                         "/mag.T" + std::to_string(params.data.L[0]) +
                         "X" + std::to_string(params.data.L[1]) +
                         "Y" + std::to_string(params.data.L[2]) +
                         "Z" + std::to_string(params.data.L[3]) +
                         "kap" + std::to_string(params.data.kappa) + 
                         "lam" + std::to_string(params.data.lambda) + 
                         ".rep_" + std::to_string(params.data.replica) +
                         ".dat";
  FILE *f_mag = fopen(mag_file.c_str(), "w"); 
  if (f_mag == NULL) {
      printf("Error opening file!\n");
      exit(1);
  }
  std::string cluster_size_file = params.data.outpath + 
                         "/clusterSize.T" + std::to_string(params.data.L[0]) +
                         "X" + std::to_string(params.data.L[1]) +
                         "Y" + std::to_string(params.data.L[2]) +
                         "Z" + std::to_string(params.data.L[3]) +
                         "kap" + std::to_string(params.data.kappa) + 
                         "lam" + std::to_string(params.data.lambda) + 
                         ".rep_" + std::to_string(params.data.replica) +
                         ".dat";
  FILE *f_cluster = fopen(cluster_size_file.c_str(), "w"); 
  if (f_cluster == NULL) {
      printf("Error opening file!\n");
      exit(1);
  }
  */
  
  mdp.close_wormholes();
  return 0;
}

