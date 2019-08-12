#include <array>
#include <cmath>
#include <ctime>
#include <vector>

#include "mdp.h"

#include "IO_params.h" 

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
inline double compute_magnetisation(mdp_field<std::array<double, 4> >& phi, 
                                    mdp_site& x){

  double m0 = 0.0, m1 = 0.0, m2 = 0.0, m3 = 0.0;
  forallsites(x){
    m0 += phi(x)[0];
    m1 += phi(x)[1];
    m2 += phi(x)[2];
    m3 += phi(x)[3];
  }
  return sqrt(m0*m0 + m1*m1 + m2*m2 + m3*m3);

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline void get_phi_field_direction(mdp_field<std::array<double, 4> >& phi, 
                                    mdp_site& x, std::array<double, 4>& dir, 
                                    const double V){

  dir = {{0.0, 0.0, 0.0, 0.0}};
  forallsites(x) {
    dir[0] += phi(x)[0];
    dir[1] += phi(x)[1];
    dir[2] += phi(x)[2];
    dir[3] += phi(x)[3];
  }
  for (int i = 0; i < 4; ++i)
    dir[i] /= V;

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline double get_angle (const double x, const double y) {

  double d = sqrt (x * x + y * y);
  if (d < 1E-10)
    return 0.0;
  double w = asin (x / d);
  if (y < 0)
    w = M_PI - w;
  return w;

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline void rotate_phi_field_component(mdp_field<std::array<double, 4> >& phi, 
                                       mdp_site& x, const int ind1, 
                                       const int ind2, const double w){
  double c = cos (w);
  double s = sin (w);

  forallsites(x) {
    double y = phi(x)[ind1];
    double z = phi(x)[ind2];
    phi(x)[ind1] =  c*y + s*z;
    phi(x)[ind2] = -s*y + c*z;
  }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void rotate_phi_field (mdp_field<std::array<double, 4> >& phi, mdp_site& x,
                       const double V) {

  double angle;
  std::array<double, 4> dir;

  get_phi_field_direction (phi, x, dir, V);
  angle = get_angle (dir[1], dir[0]);
  rotate_phi_field_component (phi, x, 1, 0, -angle);

  get_phi_field_direction (phi, x, dir, V);
  angle = get_angle (dir[2], dir[0]);
  rotate_phi_field_component (phi, x, 2, 0, -angle);

  get_phi_field_direction (phi, x, dir, V);
  angle = get_angle (dir[3], dir[0]);
  rotate_phi_field_component (phi, x, 3, 0, -angle);

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline std::array<double, 4> create_phi_update(const double delta){

  return {{(random1.plain()*2. - 1.)*delta,
           (random1.plain()*2. - 1.)*delta,
           (random1.plain()*2. - 1.)*delta,
           (random1.plain()*2. - 1.)*delta,
         }};

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
double metropolis_update(mdp_field<std::array<double, 4> >& phi, mdp_site& x,
                         const double kappa, const double lambda, 
                         const double delta, const size_t nb_of_hits){

  double acc = .0;
  for(int parity=EVEN; parity<=ODD; parity++) {
    forallsitesofparity(x, parity) {
      // computing phi^2 on x
      auto phiSqr = phi(x)[0]*phi(x)[0] + phi(x)[1]*phi(x)[1] + 
                    phi(x)[2]*phi(x)[2] + phi(x)[3]*phi(x)[3];
      // running over the four components, comp, of the phi field - Each 
      // component is updated individually with multiple hits
      for(size_t comp = 0; comp < 4; comp++){
        auto& Phi = phi(x)[comp]; // this reference gives a speedup
        // compute the neighbour sum
        auto neighbourSum = 0.0;
        for(size_t dir = 0; dir < 4; dir++) // dir = direction
          neighbourSum += phi(x-dir)[comp] + phi(x+dir)[comp];
        // doing the multihit

        for(size_t hit = 0; hit < nb_of_hits; hit++){
          auto deltaPhi = (random1.plain()*2. - 1.)*delta;
          auto deltaPhiPhi = deltaPhi * Phi;
          auto deltaPhideltaPhi = deltaPhi * deltaPhi;
          // change of action
          auto dS = -2.*kappa*deltaPhi*neighbourSum + 
                     2.*deltaPhiPhi*(1. - 2.*lambda*(1. - phiSqr - deltaPhideltaPhi)) +
                     deltaPhideltaPhi*(1. - 2.*lambda*(1. - phiSqr)) +
                     lambda*(4.*deltaPhiPhi*deltaPhiPhi + deltaPhideltaPhi*deltaPhideltaPhi);
          // Monate Carlo accept reject step -------------------------------------
          if(random1.plain() < exp(-dS)) {
            phiSqr -= Phi*Phi;
            Phi += deltaPhi;
            phiSqr += Phi*Phi;
            acc++; 
          }
        } // multi hit ends here
      } // loop over components ends here
    } // loop over parity ends here
    phi.update(parity); // communicate boundaries
  }

  return acc/(4*nb_of_hits); // the 4 accounts for updating the component indiv.

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline void check_neighbour(const size_t y, 
                            const double scalar_x,
                            const std::array<double, 4>& phi_y,
                            const std::array<double, 4>& r, 
                            size_t& cluster_size,
                            std::vector<cluster_state_t>& checked_points,
                            std::vector<int>& look, 
                            size_t& look_size){

  if(checked_points.at(y) == CLUSTER_UNCHECKED){
    double scalar_y = phi_y[0]*r[0] + phi_y[1]*r[1] + 
                      phi_y[2]*r[2] + phi_y[3]*r[3];
    double dS = scalar_x * scalar_y;
    if((dS < 0.0) && (1.-exp(dS)) > random1.plain()){
      look[look_size++] = y; // y will be used as a starting point in next iter.
      checked_points.at(y) = CLUSTER_FLIP;
      cluster_size++;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline double cluster_update(mdp_field<std::array<double, 4> >& phi, 
                      mdp_site& x, 
                      std::vector<int>& look_1, 
                      std::vector<int>& look_2,
                      const double kappa, 
                      const double min_size){

  // lookuptable to check which lattice points will be flipped
  std::vector<cluster_state_t> 
              checked_points(x.lattice().nvol, CLUSTER_UNCHECKED);

  // vector which defines rotation plane ---------------------------------------
  std::array<double, 4> r = 
                         {{random1.plain()*2.-1., random1.plain()*2.-1., 
                           random1.plain()*2.-1., random1.plain()*2.-1.}};
  double len = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);
  r[0]/=len; r[1]/=len; r[2]/=len; r[3]/=len; // normalisation

  // while-loop: until at least some percentage of the lattice is updated ------
  size_t cluster_size = 0;
  while(double(cluster_size)/x.lattice().nvol <= min_size){

    // Choose a random START POINT for the cluster: 0 <= xx < volume and check 
    // if the point is already part of another cluster - if so another start 
    // point is choosen
    size_t xx = size_t(random1.plain()*x.lattice().nvol);
    while(checked_points.at(xx) == CLUSTER_FLIP)
      xx = size_t(random1.plain()*x.lattice().nvol);
    checked_points.at(xx) = CLUSTER_FLIP;
    look_1[0] = xx;
    cluster_size++; 

    size_t look_1_size = 1, look_2_size = 0; 
    // run over both lookuptables until there are no more points to update -----
    while(look_1[0] != -1){ 
      // run over first lookuptable and building up second lookuptable
      for(const auto& x_look : look_1){
        if(x_look == -1) break;
        double scalar_x = -4.*kappa * (phi(x_look)[0]*r[0] + phi(x_look)[1]*r[1] + 
                                       phi(x_look)[2]*r[2] + phi(x_look)[3]*r[3]);
        for(size_t dir = 0; dir < 4; dir++){ 
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
        double scalar_x = -4.*kappa * (phi(x_look)[0]*r[0] + phi(x_look)[1]*r[1] + 
                                       phi(x_look)[2]*r[2] + phi(x_look)[3]*r[3]);
        for(size_t dir = 0; dir < 4; dir++){ 
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
    } // while loop to build the cluster ends here
  } // while loop to ensure minimal total cluster size ends here

  // perform the phi flip ------------------------------------------------------
  forallsites(x)
    if(checked_points.at(x.idx) == CLUSTER_FLIP){
      double scalar = -2.*(phi(x)[0]*r[0] + phi(x)[1]*r[1] + 
                           phi(x)[2]*r[2] + phi(x)[3]*r[3]);
      for(int dir = 0; dir < 4; dir++)
        phi(x)[dir] += scalar*r[dir];
    }

  return cluster_size;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {

  mdp.open_wormholes(argc,argv);

  cluster::IO_params params(argc, argv); // reading infile

  // lattice parameters
  int L[]={params.data.L[0], params.data.L[1],
           params.data.L[2], params.data.L[3]} ;
  const int V = params.data.V;
  const double kappa = 0.13137;
  const double lambda = 0.01035;
  const double delta = 4.7; // update parameter for new phi 
  const size_t nb_of_hits = 10;


  // setup the lattice and filds
  mdp_lattice hypercube(4,L); // declare lattice
  mdp_field<std::array<double, 4> > phi(hypercube); // declare phi field
  mdp_site x(hypercube); // declare lattice lookuptable

  random1.initialize(1227);

  // random start configuration
  forallsites(x)
    phi(x) = create_phi_update(1.); 

  // compute magnetisation on start config
  rotate_phi_field(phi, x, double(V));
  double M = compute_magnetisation(phi, x);
  mdp.add(M);
  mdp << "\n\n\tmagnetization at start = " << M/V << endl;

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
  /*
  if(params.data.restart != 0){
    std::string conf_file = params.data.outpath + 
                            "/T" + std::to_string(params.data.L[0]) +
                            ".X" + std::to_string(params.data.L[1]) +
                            ".Y" + std::to_string(params.data.L[2]) +
                            ".Z" + std::to_string(params.data.L[3]) +
                            ".kap" + std::to_string(params.data.kappa) + 
                            ".lam" + std::to_string(params.data.lambda)+
                            ".conf" + std::to_string(params.data.restart);
    phi.load(conf_file.c_str());
    std::string rnd_state_filename = params.data.outpath + 
                            "/T" + std::to_string(params.data.L[0]) +
                            ".X" + std::to_string(params.data.L[1]) +
                            ".Y" + std::to_string(params.data.L[2]) +
                            ".Z" + std::to_string(params.data.L[3]) +
                            ".kap" + std::to_string(params.data.kappa) + 
                            ".lam" + std::to_string(params.data.lambda)+
                            ".conf" + std::to_string(params.data.restart) + 
                            ".random_generator_state";
    random1.read_state(rnd_state_filename);
  }*/

  std::vector<int> look_1(V, -1), look_2(V, -1); // lookuptables for the cluster
  // The update ----------------------------------------------------------------
  for(int ii = 0; ii < params.data.start_measure+params.data.total_measure; ii++) {

    clock_t begin = clock(); // start time for one update step
    // metropolis update
    double acc = 0.0;
    for(int global_metro_hits = 0; 
        global_metro_hits < params.data.metropolis_global_hits; 
        global_metro_hits++)
      acc += metropolis_update(phi, x, params.data.kappa, params.data.lambda, 
                               params.data.metropolis_delta, 
                               params.data.metropolis_local_hits);
    acc /= params.data.metropolis_global_hits;

    clock_t mid = clock(); // start time for one update step

    // cluster update
    double cluster_size = 0.0 ;
    size_t tmp ;
    for(size_t nb = 0; nb < params.data.cluster_hits; nb++){
      tmp = cluster_update(phi, x, look_1, look_2, params.data.kappa, 
                                     params.data.cluster_min_size);
      cluster_size += tmp;
      //mdp << "cluster of total size = " << tmp << " formed\n" ;
      fprintf(f_cluster,"%d\t",tmp);
      }
    cluster_size /= params.data.cluster_hits;
    clock_t end = clock(); // end time for one update step
    fprintf(f_cluster,"\n");
    fflush(f_cluster);
    //mdp << "\n" << cluster_size << "\n";
    

    // compute magnetisation every ZZZ configuration
    if(ii > params.data.start_measure &&
       ii%params.data.measure_every_X_updates == 0){
      mdp_field<std::array<double, 4> > phi_rot(phi); // copy field
      rotate_phi_field(phi_rot, x, double(V)); 
      M = compute_magnetisation(phi_rot, x);
      mdp.add(M); // adding magnetisation and acceptance rate in parallel
      mdp.add(acc);
      mdp << ii << "\tmag after rot = " << M/V;
      mdp << "  \tacc. rate = " << acc/V 
          << "  \tcluster size = " << 100.*cluster_size/V 
          << "\ttime metro = " << double(mid - begin) / CLOCKS_PER_SEC 
          << "\ttime clust = " << double(end - mid) / CLOCKS_PER_SEC 
          << endl;
      fprintf(f_mag, "%.14lf\n", M/V);
      fflush(f_mag);
    }
    
    /*if(params.data.save_config == "yes" && ii > params.data.start_measure &&
       ii%params.data.save_config_every_X_updates == 0){
      std::string conf_file = params.data.outpath + 
                              "/T" + std::to_string(params.data.L[0]) +
                              ".X" + std::to_string(params.data.L[1]) +
                              ".Y" + std::to_string(params.data.L[2]) +
                              ".Z" + std::to_string(params.data.L[3]) +
                              ".kap" + std::to_string(params.data.kappa) + 
                              ".lam" + std::to_string(params.data.lambda)+
                              ".conf" + std::to_string(ii);
      phi.save(conf_file.c_str());
      std::string rnd_state_filename = params.data.outpath + 
                              "/T" + std::to_string(params.data.L[0]) +
                              ".X" + std::to_string(params.data.L[1]) +
                              ".Y" + std::to_string(params.data.L[2]) +
                              ".Z" + std::to_string(params.data.L[3]) +
                              ".kap" + std::to_string(params.data.kappa) + 
                              ".lam" + std::to_string(params.data.lambda)+
                              ".conf" + std::to_string(ii) + 
                              ".random_generator_state";
      random1.
      //random1.write_state(rnd_state_filename);
    }*/
    
  }

  // end everything
  fclose(f_mag);
  fclose(f_cluster);
  mdp.close_wormholes();
  return 0;
}
