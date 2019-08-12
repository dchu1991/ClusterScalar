#include <array>
#include <cmath>
#include <ctime>
#include <vector>

#include "mdp.h"

#include "Ising_params.h" 

// the random number generator
mdp_random_generator random1;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
typedef enum cluster_state_t {
  CLUSTER_UNCHECKED=0,
  CLUSTER_FLIP 
} cluster_state_t;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline double compute_magnetisation(mdp_field<std::array<double, 2> >& phi, 
                                    mdp_site& x){

  double m[2] = {0.0,0.0};
  forallsites(x){
    for(int i = 0; i < 2 ; i++) m[i] += phi(x)[i];
  }

  return sqrt(m[0]*m[0]+m[1]*m[1]);  
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline double compute_staggermagnetisation(
                                    mdp_field<std::array<double, 2> >& phi, 
                                    mdp_site& x){

  double m[2] = {0.0,0.0};
  forallsites(x){
    for(int i = 0; i < 2 ; i++) m[i] += phi(x)[i]*pow(-1, x.parity() );
  }

  return sqrt(m[0]*m[0]+m[1]*m[1]);  
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline double compute_totalEnergy(mdp_field<std::array<double, 2> >& phi, 
                                  mdp_site& x, const int dim){

  double E = 0.0;
  double m[2];
  forallsitesofparity(x,ODD){
  m[0] = 0.0;
  m[1] = 0.0;
    for(int dir = 0; dir < dim ; dir++) {
      for(int i = 0; i < 2 ; i++) m[i] += phi(x+dir)[i]+phi(x-dir)[i];
    }
    E += m[0]*phi(x)[0] + m[1]*phi(x)[1];
  }

  return E;  
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline double compute_TCorrelation(mdp_field<std::array<double, 2> >& phi, 
                                   mdp_site& x, const int dim){
  
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline std::array<double, 2> create_phi_update(const double delta){
  
  double random_angle = (random1.plain()*2. - 1.)*delta;
  return {cos(random_angle),sin(random_angle)};
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline double metropolis_hits(mdp_field<std::array<double, 2> >& phi ,mdp_site& x,
                              const double kappa, 
                              const size_t nb_of_hits, const size_t dim){
                              
    double acc = 0.0;
    for(int parity=EVEN; parity<=ODD; parity++){
        forallsitesofparity(x,parity){
          //multihits
          for(size_t hit = 0; hit < nb_of_hits; hit++){
            auto newphi = create_phi_update(M_PI);
            auto dS = 0.0;
            // components
            for(size_t comp = 0; comp<2; comp++){
              auto neighbour_sum = 0.0;
              for (size_t dir = 0; dir < dim; dir++)
                neighbour_sum += phi(x-dir)[comp] + phi(x+dir)[comp];
              dS += kappa*neighbour_sum*(phi(x)[comp]-newphi[comp]);    
          }// loop over components          
          //monte-carlo step
          if(random1.plain() < exp(-dS)){
            phi(x) = newphi;
            acc++;
          }         
       }//multi hits ends here
    }//parity ends here
    phi.update(parity); // communicate boundaries
  }
return  acc/(nb_of_hits);       
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline void check_neighbour(const size_t y,
                            const double scalar_x,
                            const std::array<double, 2>& phi_y,
                            const std::array<double, 2>& r,
                            size_t& cluster_size,
                            std::vector<cluster_state_t>& checked_points,
                            std::vector<int>& look,
                            size_t& look_size ){
                            
  if(checked_points.at(y) == CLUSTER_UNCHECKED){
    double scalar_y = phi_y[0]*r[0] + phi_y[1]*r[1];
    double dS = scalar_x*scalar_y;
    if(dS < 0 && 1.0-exp(dS) > random1.plain()){  
      look[look_size++] = y;
      checked_points.at(y) = CLUSTER_FLIP;
      cluster_size++;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline int cluster_update(mdp_field<std::array<double,2> >& phi, mdp_site& x,
                             const double kappa, const int dim,
                             std::vector<int>& look_1,
                             std::vector<int>& look_2, const double min_size){
  //define a rotation plane
  std::array<double,2> r = create_phi_update(M_PI);
  //look-up table to check which lattice points will be flipped
  std::vector<cluster_state_t> 
              checked_points(x.lattice().nvol, CLUSTER_UNCHECKED);
  
  //while loop to at least some percentage of cluster is marked
  size_t cluster_size = 0;
  while(double(cluster_size)/x.lattice().nvol <= min_size){
    //choose a random start point and check it's not overlapping 
    size_t start_x = size_t(random1.plain()*x.lattice().nvol);
    while(checked_points.at(start_x) == CLUSTER_FLIP)
      start_x = size_t(random1.plain()*x.lattice().nvol);
    checked_points.at(start_x) = CLUSTER_FLIP;
    look_1[0] = start_x;
    cluster_size++;
    
    size_t look_1_size = 1, look_2_size = 0; 
    // run over both lookuptables until there are no more points to update -----
    while(look_1[0] != -1){ 
      // run over first lookuptable and building up second lookuptable
      for(const auto& x_look : look_1){
        if(x_look == -1) break;
        double scalar_x = -2.0*kappa*(phi(x_look)[0]*r[0] + phi(x_look)[1]*r[1]);
        for(size_t dir = 0; dir < dim; dir++){ 
          // negative direction
          auto y = x.lattice().dw[x_look][dir];
          check_neighbour(y, scalar_x, phi(y), r, cluster_size,
                          checked_points, look_2, look_2_size);
          // positive direction
          y = x.lattice().up[x_look][dir];
          check_neighbour(y, scalar_x, phi(y), r, cluster_size,
                          checked_points, look_2, look_2_size);
        }
      }
      // run over second lookuptable and building up first lookuptable
      for(size_t i = 0; i <= look_1_size; i++)
        look_1[i] = -1;
      look_1_size = 0;
      for(const auto& x_look : look_2){
        if(x_look == -1) break;
        double scalar_x = -2.0*kappa*(phi(x_look)[0]*r[0] + phi(x_look)[1]*r[1]);
        for(size_t dir = 0; dir < dim; dir++){ 
          // negative direction
          auto y = x.lattice().dw[x_look][dir];
          check_neighbour(y, scalar_x, phi(y), r, cluster_size,
                          checked_points, look_2, look_2_size);
          // positive direction
          y = x.lattice().up[x_look][dir];
          check_neighbour(y, scalar_x, phi(y), r, cluster_size,
                          checked_points, look_2, look_2_size);
        }
      }
      for(size_t i = 0; i <= look_2_size; i++) 
        look_2[i] = -1;
      look_2_size = 0;
    } // while loop to build cluster
  }//while loop to ensure minimal cluster ends here
  
  //perform flip
  forallsites(x){
    if(checked_points.at(x.idx) == CLUSTER_FLIP){
      double scalar = -2.*(phi(x)[0]*r[0] + phi(x)[1]*r[1]);
      for(int dir = 0; dir < dim; dir++) 
        phi(x)[dir] += scalar*r[dir];      
    }
  }
return cluster_size;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){

  mdp.open_wormholes(argc,argv);

  cluster::IO_params params(argc, argv); // reading infile

  // lattice parameters
  const int dim = params.data.dimension ;
  if(dim > 9) {
    mdp << "dimension of lattice cannot be greater than 10" << endl;
    exit(1);
  }
  int L[dim] = {params.data.T};
  for(int i = 1; i < dim ; i++) L[i] = params.data.L;
  const int V = pow(params.data.L,dim-1)*params.data.T;
  //const int comp = params.data.comp;
  
  // setup the lattice and fields
  mdp_lattice hypercube(dim,L); // declare lattice
  mdp_field<std::array<double, 2> > phi(hypercube); // declare phi field
  mdp_site x(hypercube); // declare lattice lookuptable
  
  random1.initialize(params.data.seed);
  
  // random start configuration
  forallsites(x)
    phi(x) = create_phi_update(M_PI); 

  // compute magnetisation on start config
  double M = compute_magnetisation(phi, x);
  double sM = compute_staggermagnetisation(phi,x);
  double E = compute_totalEnergy(phi,x,dim);
  mdp.add(M);
  mdp << "\tmagnetization at start = " << M/V << endl;
  mdp << "stagger magnetization at start = " << sM/V << endl;
  mdp << "Total Energy at start = " << -1*E/V << endl;
  
  
  // outfiles inits
  std::string file_ending = ".T" + std::to_string(params.data.T) +
                         ".L" + std::to_string(params.data.L) +
                         ".dim" + std::to_string(params.data.dimension) +
                         ".kap" + std::to_string(params.data.kappa) +  
                         ".rep_" + std::to_string(params.data.replica) ;
  std::string mag_file = params.data.outpath + "XYmag" + file_ending +
                         ".dat";
  FILE *f_mag = fopen(mag_file.c_str(), "w"); 
  if (f_mag == NULL) {
      printf("Error opening file!\n");
      exit(1);
  } 


  for(int ii = 0; ii < params.data.start_measure+params.data.total_measure; ii++) {
    
    clock_t begin = clock(); // start time for one update step
    
    std::vector<int> look_1(V, -1), look_2(V, -1); // lookuptables for the cluster

    //metropolis update
    double acc = 0.0;
    for(int global_metro_hits = 0; 
        global_metro_hits < params.data.metropolis_global_hits; 
        global_metro_hits++)
      acc += metropolis_hits(phi, x, params.data.kappa,  
                             params.data.metropolis_local_hits, dim);
    acc /= params.data.metropolis_global_hits;
    
    clock_t mid = clock(); //end time for metropolis update step
    
    
    //cluster update
    double cluster_size = 0.0 ;
    size_t tmp ;
    for(size_t nb = 0; nb < params.data.cluster_hits; nb++){
      tmp = cluster_update(phi, x, params.data.kappa, dim, look_1, look_2,
                           params.data.cluster_min_size);
      cluster_size += tmp;
      //mdp << "cluster of total size = " << tmp << " formed\n" ;
      //fprintf(f_cluster,"%d\t",tmp);
      }
    cluster_size /= params.data.cluster_hits;
    
    clock_t end = clock(); // end time for one update step
    //fprintf(f_cluster,"\n");
    //fflush(f_cluster);

    
    if(ii > params.data.start_measure &&
       ii%params.data.measure_every_X_updates == 0){
       
      std::string Conf_file = params.data.outpath + "/Confs/XYModelConf" 
      + file_ending + ".conf" + std::to_string(ii) + ".dat";
      phi.save(Conf_file.c_str());
       
      mdp_field<std::array<double, 2> > phi_rot(phi); // copy field
      M = compute_magnetisation(phi_rot, x);
      sM = compute_staggermagnetisation(phi_rot, x);
      E = compute_totalEnergy(phi_rot,x,dim);
      mdp.add(M); // adding magnetisation and acceptance rate in parallel
      mdp.add(acc);
      mdp << ii << "\tmag = " << M/V 
          << "\t stagger mag = " << sM/V 
          << "\t energy = " << -1*E/V
          << endl;
      mdp << "  \tacc. rate = " << acc/V 
          << "  \tcluster size = " << 100.*cluster_size/V 
          << "\ttime metro = " << double(mid - begin) / CLOCKS_PER_SEC 
          << "\ttime clust = " << double(end - mid) / CLOCKS_PER_SEC 
          << "\n" << endl;
      fprintf(f_mag, "%.14lf\n", M/V);
      fflush(f_mag);
      
    }// end of measurements
  } // end of updates
   
  //end everything
  fclose(f_mag);
  mdp.close_wormholes();
  return 0;
}

