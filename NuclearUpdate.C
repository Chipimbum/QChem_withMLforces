#include "qchem.h"
#include "aimdman.h"
#include "BSetMgr.hh"
#include "convfac.h"
#include "rem_values.h"
#include "functionals.h"
#include "aimdman.h"
#include "rem_values.h"
#include "spherical_hardwall.h"
#include <libqints/qchem/aobasis.h>
using libqints::qchem::aobasis;
/******************************************************************
 *                                                                *
 *  Update nuclear coordinates and velocities per one MD step.    *
 *                                                                *
 *  JMH (3/05)                                                    *
 *                                                                *
 ******************************************************************/

//RPS(11/09): Thermostatting added -- now with Nose-Hoover! -JMH

/*BRL FSSH

  I think it's poor practice to add the thermostat from within nuclear update.
  MD applications may want to know the
  nuclear accelerations at a given timestep (including the thermostat forces),
  but as the code is currently written the accelerations
  are only finalized within NuclearUpdate.  I think it is preferable therefore
  to use AddThermostat in BOMDstep.
  For (A)FSSH we need the accelerations and I will write overloaded
  NuclearUpdate and NuclearUpdatePart2 functions, but I think
  the regular functions should also be changed as well.

*/

double NuclearUpdate(double dt, Nukes &nukes, double *NucGrad,
                     double *PrevNucGrad)
{
  // Updates nuclear coordinates by one integration step, returns the
  // associated kinetic energy.  For velocity Verlet & Beeman integrators,
  // this is only Step 1 of a two stage process.
  // PrevNucGrad is needed only for Beeman integration (o/w can pass NULL)

  static int NStep = 0;
  static double *OldNucCarts;

  int IntMethod = rem_read(REM_AIMD_INTEGRATION);
  int IThermo   = rem_read(REM_AIMD_THERMOSTAT);
  int NAtoms = rem_read(REM_NATOMS);
  int   NBas = bSetMgr.crntShlsStats(STAT_NBASIS);

  //ERM added: ----------------------------------------------------------------
  //Here I read from the rem, if I find "ML_GRAD 1" then I will read the forces.
  int predGRAD = rem_read(REM_ML_GRAD); //ERM added 
  //end ERM added--------------------------------------------------------------

  double *NucCarts = nukes.Ptr2Carts();
  double *NucVeloc = nukes.Ptr2Veloc();
  double *NucMass  = nukes.Ptr2Mass();

  arma::mat Carts_old(NucCarts, 3, NAtoms, true, true);
  arma::mat Carts_new(NucCarts, 3, NAtoms, false, true);
  arma::mat Veloc_new(NucVeloc, 3, NAtoms, false, true);

  if (NStep == 0){
     if (IntMethod == BEEMAN) {
        IntMethod = VVERLET;
        printf("Beeman's algorithm initiated with a velocity Verlet step\n");
     }
     else if (IntMethod == VERLET)
        OldNucCarts = QAllocDouble(3*NAtoms);
  }


  double NuclearKE = 0.0;
  double dt2 = dt*dt, dt_2 = 0.5*dt, dt2_2 = 0.5*dt*dt;

 //ERM added --------------------------------------------------------------------------------------------------------------
 //Here I specify file where the forces from sGDML are written at each step and from where I read.
    FILE *MLF_File;
    //MLF_File = fopen("/home/scr/fanirm/h20_dyns/FORCE_QUERY-2-good/integral/read_predF_from_file/6March/MLforces.dat", "r");
    MLF_File = fopen("/home/SCR/qcscr/test/AIMD/last_MLForces.dat", "r");
    if (NULL == MLF_File) {
        printf("Warning: The file containing ML Forces from sGDML can't be opened/ doesn't exist.\n");
    }
    double MLFor[9]; //read file into an array of double

    if (predGRAD == 1){    
    printf("QChem is using Forces predicted externally using sGMDL software. \n\n"); 
    for (int i = 0; i < 3*NAtoms; i++)
    {
        fscanf(MLF_File, "%lf,", &MLFor[i]);
        // ********** KEY CHANGE *************** 
	NucGrad[i] = MLFor[i]; //Replace the calculated NucGrad by the ML Force
        // *************************************
    }
    
    }
    else {
    printf ("\n QChem is NOT using Machine Learning Forces. To activate the feature, add ML_GRAD 1 to your input. \n");
    }


//end ERM added ----------------------------------------------------------------------------------------------------------------

  if (IntMethod == EULER){
     AddThermostat(NucVeloc,NucGrad,NucMass);
     for (int i = 0; i < NAtoms; i++)
        for (int j = 0; j < 3; j++){
           int I = 3*i + j;
	   double v_old = NucVeloc[I];
           NucVeloc[I] -= dt*NucGrad[I]/NucMass[i];
	   NucCarts[I] += dt*NucVeloc[I];
           //NucCarts[I] += dt * v_old; 
           NuclearKE += 0.5 * NucMass[i] * NucVeloc[I] * NucVeloc[I];
        }

   //ERM added: only a tmp printing while debugging ----------
   printf ("NuclearKE in the NuclearUpdate.C step AFTER EULER: %lf \n", NuclearKE); //ERM added
   printf("-------------------------\n");
   printf("NucGrad [i] AFTER  EULER\n");
   printf("-------------------------\n");

   for (int i = 0; i < 3*NAtoms; i++) {
            printf(" %f \n",  NucGrad[i]);
        }

   printf("-------------------------\n");
   printf("NucVeloc [i] AFTER  EULER\n");
   printf("-------------------------\n");

   for (int i = 0; i < 3*NAtoms; i++) {
            printf(" %f \n",  NucVeloc[i]);

        }
  //end ERM added --------------------------------------------

  }else if (IntMethod == LEAPFROG) {
     AddThermostat(NucVeloc,NucGrad,NucMass);
     if (NStep == 0){ // initialize with Euler step
        for (int i = 0; i < NAtoms; i++)
           for (int j = 0; j < 3; j++)
              NucVeloc[3*i+j] += dt_2*NucGrad[3*i+j]/NucMass[i];
     }


     for (int i = 0; i < NAtoms; i++)
        for (int j = 0; j < 3; j++){
           int I = 3*i + j;
           double v_old = NucVeloc[I];
           NucVeloc[I] -= dt*NucGrad[I]/NucMass[i];
           NucCarts[I] += dt*NucVeloc[I];
           double v_t = 0.5*(v_old + NucVeloc[I]);
           NuclearKE += 0.5*NucMass[i] * v_t * v_t;
        }

  }else if (IntMethod == VERLET) {
     AddThermostat(NucVeloc,NucGrad,NucMass);
     if (NStep == 0) {  // get r(t - dt)
        VRcopy(OldNucCarts,NucCarts,3*NAtoms);
        for (int i = 0; i < NAtoms; i++)
           for (int j = 0; j < 3; j++){
              int I = 3*i+j;
              OldNucCarts[I] -= dt*NucVeloc[I] - dt2_2*NucGrad[I]/NucMass[i];
           }
     }

     double rt,rtmdt,vt;
     for (int i = 0; i < NAtoms; i++)
        for (int j = 0; j < 3; j++){
           int I = 3*i + j;
           rt = NucCarts[I];
           rtmdt = OldNucCarts[I];
           NucCarts[I] = 2.0*rt - rtmdt - dt2*NucGrad[I]/NucMass[i];
           OldNucCarts[I] = rt;
           vt = NucCarts[I] - rtmdt;
           NucVeloc[I] = vt;
           NuclearKE += vt*vt*NucMass[i];
        }
     NuclearKE *= 0.125/dt2;
     VRscale(NucVeloc,3*NAtoms,0.5/dt);



  }else if (IntMethod == VVERLET) {

     if (IThermo == NOSE_HOOVER){ // Nose-Hoover
        double T = (double)rem_read(REM_AIMD_TEMP);
        nukes.PropagateNHC(dt,T);
     }else if (IThermo > 0) // all other thermostats
        AddThermostat(NucVeloc,NucGrad,NucMass);
     for (int i = 0; i < NAtoms; ++i) {
        for (int j = 0; j < 3; ++j) {
           int I = 3*i + j;
           NucCarts[I] += dt*NucVeloc[I] - dt2_2*NucGrad[I]/NucMass[i];

           // half-step velocities and KE
           NucVeloc[I] -= dt_2*NucGrad[I]/NucMass[i];
           NuclearKE += 0.5*NucVeloc[I]*NucVeloc[I]*NucMass[i];
        }
     }
  //BR -- Step-1 for Rattle Algorithm
  }else if (IntMethod == RATTLE){
     //Rattle Algorithm - H.C. Andersen, J. Comp. Phys. 52, 24--34 (1983)
     int maxit = rem_read(REM_RATTLE_MAXIT);
     double thresh = 1.0/pow(10.0,double(rem_read(REM_RATTLE_THRESH)));
     int NBonds = rem_read(REM_RATTLE_NBONDS);
     bool converged;
     int iter;

     //read the bond index and bondlengths from the file
     int *bond_index = QAllocINTEGER(NBonds*2);
     double *bond_length = QAllocDouble(NBonds);
     FileMan(FM_READ,FILE_RATTLE_BONDINDEX,FM_INT,NBonds*2,0,FM_BEG,bond_index);
     FileMan(FM_READ,FILE_RATTLE_BONDLENGTHS,FM_DP,NBonds,0,FM_BEG,bond_length);
    
     if (IThermo == NOSE_HOOVER){ // Nose-Hoover
        double T = (double)rem_read(REM_AIMD_TEMP);
        nukes.PropagateNHC(dt,T);
     }else if (IThermo > 0) // all other thermostats
        AddThermostat(NucVeloc,NucGrad,NucMass);

     //q[t] = v[t] + 0.5 * dt * F(t)/m(t)
     double *q = QAllocDouble(3*NAtoms);
     for (int i = 0; i < NAtoms; ++i) {
        for (int j = 0; j < 3; ++j) {
           int I = 3*i + j;
           q[I] = NucVeloc[I] - dt_2*NucGrad[I]/NucMass[i];
        }
     }
  
    //Beginning of the first iteration
    converged = false; iter = 1;
    while ( !converged && iter <= maxit ){
         
        double d,s,rij,sr,s2,g,tmp;
        int I,J,m,n;

        converged = true;
        for(int i=0; i < NBonds; i++){
           d = bond_length[i];  
           m = bond_index[2*i];
           n = bond_index[2*i+1]; 
           sr=s2=0.0;
           for (int k=0; k < 3;k++){
              I = 3*m + k;
              J = 3*n + k;
              rij = NucCarts[I]- NucCarts[J];
              s = rij + dt*(q[I] - q[J]);
              s2 += s*s;
              sr += s*rij;
           }
           if (fabs(s2-d*d) > thresh){
              converged=false;
              tmp = (0.5*(s2-d*d)*NucMass[m]*NucMass[n]);
              g = tmp/(dt*sr*(NucMass[m]+NucMass[n]));
              for (int k=0; k < 3; k++){
                 I = 3*m + k;
                 J = 3*n + k;
                 rij = g*(NucCarts[I] - NucCarts[J]);
                 q[I] -= rij/NucMass[m];
                 q[J] += rij/NucMass[n];
              }
           }
        }
        iter++;
   }

   if(iter>maxit) QCrash("Max number of iterations reached for RATTLE (Step 1) w/o convergence"); 

   #ifdef DEVELOPMENT
   cout << "Number of iterations for RATTLE (Step 1) = " << iter << endl;
   #endif

   for (int i = 0; i < NAtoms; ++i) {
      for (int j = 0; j < 3; ++j) {
         int I = 3*i + j;
         NucCarts[I] += dt*q[I];

         // half-step velocities and KE
         NucVeloc[I] = q[I];
         NuclearKE += 0.5*NucVeloc[I]*NucVeloc[I]*NucMass[i];
      }
   }
  QFree(q);
  QFree(bond_index); 
  QFree(bond_length);
  }else if (IntMethod == BEEMAN){
     AddThermostat(NucVeloc,NucGrad,NucMass);
     double dt_6 = dt/6.0, dt2_6 = dt2/6.0;
     for (int i = 0; i < NAtoms; i++)
        for (int j = 0; j < 3; j++){
           int I = 3*i + j;
           NucCarts[I] += dt*NucVeloc[I]
               - dt2_6*(4.0*NucGrad[I] - PrevNucGrad[I])/NucMass[i];

           // half-step velocities and KE
           NucVeloc[I] -= dt_6 * (5.0*NucGrad[I] - PrevNucGrad[I])/NucMass[i];
           NuclearKE += 0.5*NucVeloc[I]*NucVeloc[I]*NucMass[i];
        }

  }else QCrash("Illegal integration method");


  if(rem_read(REM_BC_SPHERICAL_RADIUS) > 0)
    apply_spherical_hardwall(Carts_old, Carts_new, Veloc_new);

  if(rem_read(REM_EWALD_ON)==1){
     double* NucCartsTemp = QAllocDouble(3*NAtoms);
     VRcopy(NucCartsTemp,NucCarts,3*NAtoms);
     VRload(NucCarts,3*NAtoms,0.0);
     EnforceBox(NucCarts,NucCartsTemp);
     QFree(NucCartsTemp);
  }

  NStep++;
  return NuclearKE;
}





double NuclearUpdatePart2(double dt, Nukes &nukes, double *NucGrad)
{
  static int NStep = 0;
  int IntMethod = rem_read(REM_AIMD_INTEGRATION);
  int IThermo = rem_read(REM_AIMD_THERMOSTAT);
  int   NBas = 0;
  if(rem_read(REM_AIMD_CHILD) == 1) NBas = aobasis.b1.get_nbsf();
  else NBas = bSetMgr.crntShlsStats(STAT_NBASIS);
  int NAtoms = rem_read(REM_NATOMS);

  double *NucCarts = nukes.Ptr2Carts();
  double *NucVeloc = nukes.Ptr2Veloc();
  double *NucMass  = nukes.Ptr2Mass();

  arma::mat Carts_old(NucCarts, 3, NAtoms, true, true);
  arma::mat Carts_new(NucCarts, 3, NAtoms, false, true);
  arma::mat Veloc_new(NucVeloc, 3, NAtoms, false, true);

  if (NStep == 0 && IntMethod == BEEMAN) IntMethod = VVERLET;

  double NuclearKE = 0.0;
  double dt_2 = 0.5*dt;
  double kT_AU_to_Kelvin = ConvFac(HARTREES_TO_JOULES)/ConvFac(BOLTZMANN_CONSTANT_IN_J_K);

  if (IntMethod == VVERLET){
     if (IThermo > 0 && IThermo != NOSE_HOOVER)
        AddThermostat(NucVeloc,NucGrad,NucMass);

     for (int i = 0; i < NAtoms; i++)
        for (int j = 0; j < 3; j++){
           int I = 3*i + j;
           NucVeloc[I] -= dt_2*NucGrad[I]/NucMass[i];
           NuclearKE += 0.5*NucMass[i]*NucVeloc[I]*NucVeloc[I];
        }

     if (IThermo == NOSE_HOOVER){
        double T = (double)rem_read(REM_AIMD_TEMP);
        nukes.PropagateNHC(dt,T);
        double h = nukes.NHC_Conserved_Qty();
        printf(" Value of NHC conserved quantity = %16.6e\n",h);
        double kT = 2.0*nukes.KE()/(3.0*(double)NAtoms);
        T = kT*kT_AU_to_Kelvin;
        printf(" Instantaneous temperature = %.1f K\n",T);
        //nukes.PrintNHC();
     }
  }
  else if (IntMethod == BEEMAN){
     AddThermostat(NucVeloc,NucGrad,NucMass);
     double dt_3 = dt/3.0;
     for (int i = 0; i < NAtoms; i++)
        for (int j = 0; j < 3; j++){
           int I = 3*i + j;
           NucVeloc[I] -= dt_3*NucGrad[I]/NucMass[i];
           NuclearKE += 0.5*NucMass[i]*NucVeloc[I]*NucVeloc[I];
        }

  }
  else if(IntMethod == RATTLE){

     int maxit = rem_read(REM_RATTLE_MAXIT);
     double thresh = 1.0/pow(10.0,double(rem_read(REM_RATTLE_THRESH)));
     int NBonds = rem_read(REM_RATTLE_NBONDS);
     bool converged;
     int iter;

     int *bond_index = QAllocINTEGER(NBonds*2);
     double *bond_length = QAllocDouble(NBonds);
     FileMan(FM_READ,FILE_RATTLE_BONDINDEX,FM_INT,NBonds*2,0,FM_BEG,bond_index);
     FileMan(FM_READ,FILE_RATTLE_BONDLENGTHS,FM_DP,NBonds,0,FM_BEG,bond_length);

     if (IThermo > 0 && IThermo != NOSE_HOOVER)
        AddThermostat(NucVeloc,NucGrad,NucMass);

     for (int i = 0; i < NAtoms; i++){
         for (int j = 0; j < 3; j++){
            int I = 3*i + j;
            NucVeloc[I] -= dt_2*NucGrad[I]/NucMass[i];
         }
     }

     converged = false; iter = 1;
     while ( !converged && iter <= maxit ){
          double d,rij,vij,rv,g;
          int I,J,m,n;

          converged = true;
          for(int i=0; i < NBonds; i++){

             m = bond_index[2*i];
             n = bond_index[2*i+1];
             
             rv=0.0;
             for (int k=0; k < 3;k++){
                 I = 3*m + k;
                 J = 3*n + k;
                 rij = NucCarts[I]- NucCarts[J];
                 vij = NucVeloc[I] - NucVeloc[J];
                 rv += rij*vij;
             }
             if (fabs(rv) > thresh){
                converged=false;
                d = bond_length[i];
                g = (rv*NucMass[m]*NucMass[n])/(d*d*(NucMass[m]+NucMass[n]));
                for (int k=0; k < 3; k++){
                  I = 3*m + k;
                  J = 3*n + k;
                  rij = g*(NucCarts[I] - NucCarts[J]);
                  NucVeloc[I] -= rij/NucMass[m];
                  NucVeloc[J] += rij/NucMass[n];
                }
             }
         }

         iter++;
    }

    if(iter>maxit) QCrash("Max number of iterations reached for RATTLE (Step 2) w/o convergence.");

    #ifdef DEVELOPMENT
    cout << "Number of iterations for RATTLE (Step 2) = " << iter << endl;
    #endif

    for (int i = 0; i < NAtoms; ++i) {
        for (int j = 0; j < 3; ++j) {
           int I = 3*i + j;
           NuclearKE += 0.5*NucVeloc[I]*NucVeloc[I]*NucMass[i];
        }
    }

    if (IThermo == NOSE_HOOVER){
        double T = (double)rem_read(REM_AIMD_TEMP);
        nukes.PropagateNHC(dt,T);
        double h = nukes.NHC_Conserved_Qty();
        printf(" Value of NHC conserved quantity = %16.6e\n",h);
        double kT = 2.0*nukes.KE()/(3.0*(double)NAtoms);
        T = kT*kT_AU_to_Kelvin;
        printf(" Instantaneous temperature = %.1f K\n",T);
   }
  QFree(bond_index); 
  QFree(bond_length);
  }
  else QCrash("Why are you here??");

  if(rem_read(REM_BC_SPHERICAL_RADIUS) > 0)
    apply_spherical_hardwall(Carts_old, Carts_new, Veloc_new);

  NStep++;
  return NuclearKE;
}
