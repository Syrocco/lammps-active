#include "fix_underdamped_active_2D.h"
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "math_extra.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "group.h"


using namespace LAMMPS_NS;
using namespace FixConst;


/* ---------------------------------------------------------------------- */

FixUnderdampedActive2D::FixUnderdampedActive2D(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
  
{
  if (narg != 8) error->all(FLERR,"Illegal fix active_3d command");


  T = utils::numeric(FLERR,arg[3],false,lmp);
  gammaT = utils::numeric(FLERR,arg[4],false,lmp);
  Dr = utils::numeric(FLERR,arg[5],false,lmp);
  vo = utils::numeric(FLERR,arg[6],false,lmp);
  seed = utils::numeric(FLERR,arg[7],false,lmp);
  if (seed <= 0) error->all(FLERR,"Illegal fix active 3d command");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);
}
/* ---------------------------------------------------------------------- */

FixUnderdampedActive2D::~FixUnderdampedActive2D()
{

  delete random;

}

/* ---------------------------------------------------------------------- */


int FixUnderdampedActive2D::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixUnderdampedActive2D::init()
{
  double *angle = atom->q;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  for (int i = 0; i < nlocal; i++){       
      // Initialise Active vector with random direction at beginining of simulation
      angle[i]  = (random->uniform() - 0.5)*2*M_PI;       
  }


}

/* ----------------------------------------------------------------------
   set current t_target and t_sqrt
------------------------------------------------------------------------- */


void FixUnderdampedActive2D::initial_integrate(int vflag)
{
 
  // function to update x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *angle = atom->q;
  double *rmass = atom->rmass;
  double mag;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int step = update->ntimestep;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA
  dt = update->dt;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      
      angle[i] += sqrt(2*dt*Dr)*random->gaussian();

      v[i][0] += (-gammaT*v[i][0]*dt + dt*f[i][0] + sqrt(2*dt*T*gammaT)*random->gaussian() + gammaT*vo*cos(angle[i])*dt)/rmass[i];
      v[i][1] += (-gammaT*v[i][1]*dt + dt*f[i][1] + sqrt(2*dt*T*gammaT)*random->gaussian() + gammaT*vo*sin(angle[i])*dt)/rmass[i];
      x[i][0] +=  v[i][0]*dt;
      x[i][1] +=  v[i][1]*dt; 
                                                                        
      }
      
}

/* ---------------------------------------------------------------------- */
