#ifdef FIX_CLASS

FixStyle(underdampedActive2D,FixUnderdampedActive2D)

#else

#ifndef LMP_FIX_underdampedActive2D_H
#define LMP_FIX_underdampedActive2D_H

#include "fix.h"

namespace LAMMPS_NS {

class FixUnderdampedActive2D : public Fix {
 public:
  FixUnderdampedActive2D(class LAMMPS *, int, char **);
  virtual ~FixUnderdampedActive2D();
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);


 private: 
 double dt, sqrtdt;

 protected:
  class RanMars *random;
  int seed;
  double T, Dr, vo, gammaT;
  char *id_temp;

};

}

#endif
#endif