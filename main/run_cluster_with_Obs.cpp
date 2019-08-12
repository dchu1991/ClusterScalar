#include <array>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>

#include <fftw3.h>

#include "mdp.h"

#include "IO_params.h" 

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
////////////////////////  Functions for Propagators  ///////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline void get_phi_field_unit_vec(mdp_field<std::array<double, 4> >& phi, 
                                   mdp_site& x, std::array<double, 4>& dir){

  double inv_length;
  dir = {{0.0, 0.0, 0.0, 0.0}};
  forallsites(x) {
    dir[0] += phi(x)[0];
    dir[1] += phi(x)[1];
    dir[2] += phi(x)[2];
    dir[3] += phi(x)[3];
  }
  inv_length = 1/sqrt( dir[0]*dir[0] + dir[1]*dir[1] +
		                   dir[2]*dir[2] + dir[3]*dir[3] );
  for (int i = 0; i < 4; i++)
    dir[i] *= inv_length; // such that ||dir|| = 1

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline void Rescale(mdp_field<std::array<double, 4> >& phi_rescale, mdp_site& x,
                    const double Scale_factor){
  forallsites(x){
    for (int i = 0; i < 4; i++ )  
      phi_rescale(x)[i] *= sqrt(Scale_factor);
  }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline int Flag(const double momSqr, const int SlotCnt, 
         const std::vector<double>& DifferentMomenta){

  for (int I = 0; I < SlotCnt; I++){
    if ( fabs( momSqr - DifferentMomenta[I] ) < 1E-9 ){
      return I;
    }
  }
  return -1;

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// holy shit......


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline void Not_Projection(mdp_field<std::array<double, 4> >& phi, mdp_site& x,
		                    fftw_complex* const output){
		                   
  std::array<double,4> dir;
  get_phi_field_unit_vec(phi, x, dir);
  
  forallsites(x){
    for(int i = 0; i < 4 ; i++){
      output[4*x.global_index()+i][0] = phi(x)[i];
      output[4*x.global_index()+i][1] = 0.0;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline double PropComponent(const int comp, fftw_complex const * const output,
                   const double V, const int i){
  double tmp = output[4*comp+i][0] * output[4*comp+i][0] +
               output[4*comp+i][1] * output[4*comp+i][1];
  return tmp/V;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline double GetPropComponent(fftw_complex const * const output,
                    const std::vector<double>& sinPSqr, 
                    const std::vector<double>& DifferentMomenta, 
                    const int comp, const int V, const int j){
  int counter = 0;
  double tmp = 0.0;
  for (int i = 0; i < V; i++){
    if ( fabs(DifferentMomenta[comp]-sinPSqr[i]) < 1E-9 ){
      tmp += PropComponent(i, output, V, j);
      counter++;
    }
  }
  return tmp/counter;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline void Compute_UnprojectedCorrelatorObs(mdp_field<std::array<double,4> >& phi,
                               mdp_site& x, 
                               std::vector<double>& CorrelatorObs,
                               const double factor){

  forallsites(x){
    CorrelatorObs[4*x(0)+0] += phi(x)[0];
    CorrelatorObs[4*x(0)+1] += phi(x)[1];
    CorrelatorObs[4*x(0)+2] += phi(x)[2];
    CorrelatorObs[4*x(0)+3] += phi(x)[3];
  }
  
  for(std::vector<int>::size_type i = 0; i != CorrelatorObs.size(); i++)    
    CorrelatorObs[i] *= factor;
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

  return {{(mdp_random.plain()*2. - 1.)*delta,
           (mdp_random.plain()*2. - 1.)*delta,
           (mdp_random.plain()*2. - 1.)*delta,
           (mdp_random.plain()*2. - 1.)*delta,
         }};

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
double metropolis_update(mdp_field<std::array<double, 4> >& phi, mdp_site& x,
                         const double kappa, const double lambda, 
                         const double delta, const size_t nb_of_hits){

  double acc = .0;
  for(int parity=EVEN; parity<=ODD; parity++) {
    forallsitesofparity(x,parity) {
      // computing phi^2 on x
      auto phiSqr = phi(x)[0]*phi(x)[0] + phi(x)[1]*phi(x)[1] + 
                    phi(x)[2]*phi(x)[2] + phi(x)[3]*phi(x)[3];
      // running over the four components, comp, of the phi field - Each 
      // component is updated individually with multiple hits
      for(size_t comp = 0; comp < 4; comp++){
        auto& Phi = phi(x)[comp]; // just a copy for simplicity
        // compute the neighbour sum
        auto neighbourSum = 0.0;
        for(size_t dir = 0; dir < 4; dir++) // dir = direction
          neighbourSum += phi(x-dir)[comp] + phi(x+dir)[comp];
        // doing the multihit
        for(size_t hit = 0; hit < nb_of_hits; hit++){
          auto deltaPhi = (mdp_random.plain()*2. - 1.)*delta;
          auto deltaPhiPhi = deltaPhi * Phi;
          auto deltaPhideltaPhi = deltaPhi * deltaPhi;
          // change of action
          auto dS = -2.*kappa*deltaPhi*neighbourSum + 
                     2.*deltaPhiPhi*(1. - 2.*lambda*(1. - phiSqr - deltaPhi*deltaPhi)) +
                     deltaPhideltaPhi*(1. - 2.*lambda*(1. - phiSqr)) +
                     lambda*(4.*deltaPhiPhi*deltaPhiPhi + deltaPhideltaPhi*deltaPhideltaPhi);
          // Monate Carlo accept reject step -------------------------------------
          if(mdp_random.plain() < exp(-dS)) {
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
inline void check_neighbour(const size_t x_look, const size_t y, 
                            const double kappa, 
                            mdp_field<std::array<double, 4> >& phi,
                            const std::array<double, 4>& r, size_t& cluster_size,
                            std::vector<cluster_state_t>& checked_points,
                            std::vector<size_t>& look){

  if(checked_points.at(y) == CLUSTER_UNCHECKED){
    double scalar_x = phi(x_look)[0]*r[0] + phi(x_look)[1]*r[1] + 
                      phi(x_look)[2]*r[2] + phi(x_look)[3]*r[3];
    double scalar_y = phi(y)[0]*r[0] + phi(y)[1]*r[1] + 
                      phi(y)[2]*r[2] + phi(y)[3]*r[3];
    double dS = -4.*kappa * scalar_x * scalar_y;
    if((dS < 0.0) && (1.-exp(dS)) > mdp_random.plain()){
      look.emplace_back(y); // y will be used as a starting point in next iter.
      checked_points.at(y) = CLUSTER_FLIP;
      cluster_size++;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
double cluster_update(mdp_field<std::array<double, 4> >& phi, mdp_site& x, 
                      const double kappa, const double min_size){

  // lookuptable to check which lattice points will be flipped
  std::vector<cluster_state_t> 
              checked_points(x.lattice().nvol, CLUSTER_UNCHECKED);

  // vector which defines rotation plane ---------------------------------------
  std::array<double, 4> r = 
                         {{mdp_random.plain()*2.-1., mdp_random.plain()*2.-1., 
                           mdp_random.plain()*2.-1., mdp_random.plain()*2.-1.}};
  double len = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);
  r[0]/=len; r[1]/=len; r[2]/=len; r[3]/=len; // normalisation

  // while-loop: until at least some percentage of the lattice is updated ------
  size_t cluster_size = 0;
  while(double(cluster_size)/x.lattice().nvol <= min_size){

    // lookuptables to build the cluster
    std::vector<size_t> look_1(0, 0), look_2(0, 0);

    // Choose a random START POINT for the cluster: 0 <= xx < volume and check 
    // if the point is already part of another cluster - if so another start 
    // point is choosen
    size_t xx = size_t(mdp_random.plain()*x.lattice().nvol);
    while(checked_points.at(xx) == CLUSTER_FLIP)
      xx = size_t(mdp_random.plain()*x.lattice().nvol);
    checked_points.at(xx) = CLUSTER_FLIP;
    look_1.emplace_back(xx);
    cluster_size++; 
 
    // run over both lookuptables until there are no more points to update -----
    while(look_1.size()){ 
      // run over first lookuptable and building up second lookuptable
      look_2.resize(0);
      for(const auto& x_look : look_1){ 
        for(size_t dir = 0; dir < 4; dir++){ 
          // negative direction
          auto y = x.lattice().dw[x_look][dir];
          check_neighbour(x_look, y, kappa, phi, r, cluster_size,
                          checked_points, look_2);
          // positive direction
          y = x.lattice().up[x_look][dir];
          check_neighbour(x_look, y, kappa, phi, r, cluster_size,
                          checked_points, look_2);
        }
      }
      // run over second lookuptable and building up first lookuptable
      look_1.resize(0);
      for(const auto& x_look : look_2){ 
        for(size_t dir = 0; dir < 4; dir++){ 
          // negative direction
          auto y = x.lattice().dw[x_look][dir];
          check_neighbour(x_look, y, kappa, phi, r, cluster_size,
                          checked_points, look_1);
          // positive direction
          y = x.lattice().up[x_look][dir];
          check_neighbour(x_look, y, kappa, phi, r, cluster_size,
                          checked_points, look_1);
        }
      }
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
  const double rescale_factor = double(L[0])/double(V);
  const double kappa = 0.13137;
  const double lambda = 0.01035;
  const double delta = 4.7; // update parameter for new phi 
  const size_t nb_of_hits = 10;

  // setup the lattice and filds
  mdp_lattice hypercube(4,L); // declare lattice
  mdp_field<std::array<double, 4> > phi(hypercube); // declare phi field
  mdp_site x(hypercube); // declare lattice lookuptable

  // initialise the random number generator
  mdp_random.initialize(params.data.seed);

  // random start configuration
  forallsites(x)
    phi(x) = create_phi_update(1.); 
    
  // compute magnetisation on start config
  rotate_phi_field(phi, x, double(V));
  double M = compute_magnetisation(phi, x);
  mdp.add(M);
  mdp << "\n\n\tmagnetization at start = " << M/V << endl;

  // creating output file names and files *************************************
  std::string file_ending = ".X" + std::to_string(params.data.L[1]) +
                            ".Y" + std::to_string(params.data.L[2]) +
                            ".Z" + std::to_string(params.data.L[3]) +
                            ".kap" + std::to_string(params.data.kappa) + 
                            ".lam" + std::to_string(params.data.lambda) + 
                            ".rep_" + std::to_string(params.data.replica);
  std::string mag_file = params.data.outpath + "/mag.T" + 
                         std::to_string(params.data.L[0]) + file_ending + ".dat";
  std::string HiggsGoldstoneProp_file = params.data.outpath + "/HiggsPropagator.T" + 
                               std::to_string(params.data.L[0]) + file_ending + ".dat";
  std::string UnprojectedProp_file = params.data.outpath + "/UnProjectedPropagator.T" + 
                                   std::to_string(params.data.L[0]) + file_ending + ".dat";
  std::string HiggsGoldstoneCorr_file = params.data.outpath + "/HiggsCorrelator.T" + 
                               std::to_string(params.data.L[0]) + file_ending + ".dat";
  std::string UnprojectedCorr_file = params.data.outpath + "/UnProjectedCorrelator.T" + 
                                   std::to_string(params.data.L[0]) + file_ending + ".dat";
  std::string AvgPhiDir_file = params.data.outpath + "/AvgPhiDir.T" + 
                         std::to_string(params.data.L[0]) + file_ending + ".dat";     
                           
  FILE *f_mag = fopen(mag_file.c_str(), "w"); 
  FILE *f_HiggsGoldstoneProp = fopen(HiggsGoldstoneProp_file.c_str(), "wb"); 
  FILE *f_UnprojectedProp = fopen(UnprojectedProp_file.c_str(), "wb");
  FILE *f_HiggsGoldstoneCorr = fopen(HiggsGoldstoneCorr_file.c_str(), "wb"); 
  FILE *f_UnprojectedCorr = fopen(UnprojectedCorr_file.c_str(), "wb");
  FILE *f_AvgPhiDir = fopen(AvgPhiDir_file.c_str(), "w");
  

  if (f_mag == NULL || f_HiggsGoldstoneProp == NULL 
   ||f_UnprojectedProp == NULL || f_HiggsGoldstoneCorr == NULL ||
   f_UnprojectedCorr == NULL ||
   f_AvgPhiDir == NULL) {
      cout << params.data.outpath << endl;
      printf("Error opening data files\n");
      exit(1);
  }
  
  // Propagator initiation ****************************************************
  // ini FFT by creating a plan at first
  fftw_complex* output = new fftw_complex[4*V]; 
 
  int n[4],inembed[4],onembed[4];
  for (int j = 0; j < 4; j++){
    n[j] = params.data.L[j];
    inembed[j] = params.data.L[j];
    onembed[j] = params.data.L[j];
  }
  
  fftw_plan Plan = fftw_plan_many_dft(4, n, 4, &(output[0]), inembed, 4, 1,
	   		                              &(output[0]), onembed, 4, 1,
	  			                            FFTW_FORWARD, FFTW_MEASURE);
	  			                            
	fftw_complex* Unprojected_output = new fftw_complex[4*V]; 
   
  fftw_plan Unprojected_Plan = fftw_plan_many_dft(4, n, 4, &(Unprojected_output[0]), 
                                      inembed, 4, 1,
	   		                              &(Unprojected_output[0]), onembed, 4, 1,
	  			                            FFTW_FORWARD, FFTW_MEASURE);
  
  // create a list of \sum sin^2(P/2)
  std::vector<double> sinPSqr(V);
  std::vector<double> DifferentMomenta(V);
  std::array<double,4> p;
  int ctr = 0;
  int SlotCnt = 0;
  for (int x0 = 0; x0 < params.data.L[0]; x0++){
  p[0] = x0*M_PI/L[0]; // half-momentum
    for (int x1 = 0; x1 < params.data.L[1]; x1++){
    p[1] = x1 * M_PI/L[1];
      for (int x2 = 0; x2 < params.data.L[2]; x2++){
      p[2] = x2 * M_PI/L[2];
        for (int x3 = 0; x3 < params.data.L[3]; x3++){
        p[3] = x3 * M_PI/L[3];
        
        sinPSqr[ctr] = 4.0 * ( sin(p[0])*sin(p[0]) + sin(p[1])*sin(p[1]) +
                               sin(p[2])*sin(p[2]) + sin(p[3])*sin(p[3]) );
        
        int flag = Flag(sinPSqr[ctr], SlotCnt, DifferentMomenta);
        
        if (flag < 0){
          DifferentMomenta[SlotCnt] = sinPSqr[ctr];
          SlotCnt++;
        }
         
        ctr++;
        }
      }
    }
  }
  DifferentMomenta.resize(SlotCnt);
  printf("\n\n\tThere are %d distinct momenta in the end.\n",SlotCnt);
  std::sort(DifferentMomenta.begin(), DifferentMomenta.end());
  
  // initlize correlator observables
  std::vector<double> UnprojectedCorrelator(4*params.data.L[0]);
  std::vector<double> HiggsGoldstoneCorrelator(4*params.data.L[0]);
  std::vector<double> HiggsCorrelator(params.data.L[0]);
  
  // DIR
  std::array<double,4> AvgPhiDir;

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

    // cluster update
    double cluster_size = 0.0;
    for(size_t nb = 0; nb < params.data.cluster_hits; nb++)
      cluster_size += cluster_update(phi, x, params.data.kappa, 
                                     params.data.cluster_min_size);
    cluster_size /= params.data.cluster_hits;

    // compute observables every ZZZ configuration
    if(ii > params.data.start_measure &&
       ii%params.data.measure_every_X_updates == 0){
      
      
      mdp_field<std::array<double, 4> > phi_rot(phi); // copy field
      rotate_phi_field(phi_rot, x, double(V)); 
      M = compute_magnetisation(phi_rot, x);
      mdp.add(M); // adding magnetisation and acceptance rate in parallel
      mdp.add(acc);
      fprintf(f_mag, "%.14lf\n", M/V);
      fflush(f_mag);      
      
      // save average vector field...
      get_phi_field_unit_vec(phi, x, AvgPhiDir);
      fprintf(f_AvgPhiDir, "%.14lf %.14lf %.14lf %.14lf\n", 
      AvgPhiDir[0], AvgPhiDir[1], AvgPhiDir[2], AvgPhiDir[3]);
      fflush(f_AvgPhiDir);
      

    	///// Propagator working zone
    	// get re-scaled field.
    	mdp_field< std::array<double, 4> > phi_rescale(phi);
      Rescale(phi_rescale, x, 2*params.data.kappa);
    	mdp_field< std::array<double, 4> > phi_rescale_rot(phi_rot);
      Rescale(phi_rescale_rot, x, 2*params.data.kappa);
    	
    	// get projected modes
    	Not_Projection(phi_rescale_rot, x, output);
    	
    	// get not-projected modes
    	Not_Projection(phi_rescale, x, Unprojected_output);
            
    	// execute plan
    	fftw_execute(Plan);
    	fftw_execute(Unprojected_Plan);
    	
    	const int keep_components = 100;
    	std::array<double,4*keep_components> HiggsGoldstonePropOut;
    	std::array<double,4*keep_components> PropagatorOut;
    	
    	// computing components
    	for (int j = 0; j < keep_components; j++){         
        for (int k = 0; k < 4; k ++){
            HiggsGoldstonePropOut[4*j+k] = GetPropComponent(output,
            sinPSqr, DifferentMomenta, j, V, k);
            PropagatorOut[4*j+k] = GetPropComponent(Unprojected_output,
            sinPSqr, DifferentMomenta, j, V, k);
          }
    	}
      // writing data to file
//      for (int l = 0; l < keep_components-1; l++){
//        fprintf(f_Higgs, "%.14lf ", HiggsPropOut[l]);
//        fprintf(f_Goldstone, "%.14lf ", GoldstonePropOut[l]);
//      }
//      fprintf(f_Higgs, "%.14lf \n", HiggsPropOut[keep_components-1]);
//      fprintf(f_Goldstone, "%.14lf \n", GoldstonePropOut[keep_components-1]);

      fwrite(&HiggsGoldstonePropOut[0], sizeof(double), 4*keep_components, f_HiggsGoldstoneProp);
      fflush(f_HiggsGoldstoneProp);
      fwrite(&PropagatorOut[0], sizeof(double), 4*keep_components, f_UnprojectedProp);
      fflush(f_UnprojectedProp);
      // end of propagator
      
      // correlators

      Compute_UnprojectedCorrelatorObs(phi_rescale,x,UnprojectedCorrelator,rescale_factor);
      Compute_UnprojectedCorrelatorObs(phi_rescale_rot,x,HiggsGoldstoneCorrelator,rescale_factor);
      //Compute_HiggsCorrelatorObs(phi_rescale,x,HiggsCorrelator,rescale_factor);
      
      fwrite(&HiggsGoldstoneCorrelator[0], sizeof(double), 4*params.data.L[0], f_HiggsGoldstoneCorr);
      fflush(f_HiggsGoldstoneCorr);
      fwrite(&UnprojectedCorrelator[0], sizeof(double), 4*params.data.L[0], f_UnprojectedCorr);
      fflush(f_UnprojectedCorr);
      //end of correlators
      
      /*
      mdp << "Correlator Components : " ;
      for(std::vector<int>::size_type i = 0; i != HiggsCorrelator.size(); i++)
        mdp << " " << HiggsCorrelator[i] << "\t" << HiggsGoldstoneCorrelator[4*i] << endl;
      mdp << endl;
      */
      
      clock_t end = clock(); // end time for one update step
      mdp << "\n" << ii << "\tmag after rot = " << M/V;
      mdp << "  \tacc. rate = " << acc/V 
          << "  \tcluster size = " << 100.*cluster_size/V 
          << "  \ttime for 1 update= " << double(end - begin) / CLOCKS_PER_SEC 
          << endl;
      fflush(stdout);	
    }// end of cumputing observables
  }// end of the update
  
  // end everything
  fftw_destroy_plan(Plan);
  fftw_destroy_plan(Unprojected_Plan);

  fclose(f_HiggsGoldstoneProp);
  fclose(f_UnprojectedProp);
  fclose(f_HiggsGoldstoneCorr); 
  fclose(f_UnprojectedCorr);
  fclose(f_AvgPhiDir);

  fclose(f_mag);

  mdp.close_wormholes();
  return 0;
}
