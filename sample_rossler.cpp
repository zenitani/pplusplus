#include <particle.h>
#include <stdio.h>

// ************* Rossler attractor ************************************
// O. E. Rossler, Phys. Lett. A, 57A, 397-398 (1976)
// doi:10.1016/0375-9601(76)90101-8
const double a = 0.2;
const double b = 0.2;
const double c = 5.7;
// const double a = 0.1;
// const double b = 0.1;
// const double c = 14.0;
// ************* Rossler attractor ************************************

vector3 F( const vector3& _r, const vector3& _v,
           const  double& _t, const  double& _q ){
  return vector3( - _v.y - _v.z,
                  _v.x + a * _v.y,
                  (_r.x*_v.z + _v.x*_r.z) - c*_v.z );
}
int main()
{
  particle p;
  p.sett(0);
  p.setm(1); // unused
  p.setq(1); // unused
  p.setr( 0.0, -6.78, 0.0 );
  p.setv( - p.r.y - p.r.z,
          p.r.x + a*p.r.y,
          b + p.r.x*p.r.z - c*p.r.z);

  // marching in time
  for( int i=0;p.gett()<=200.001; i++ ){
    printf( "%f %f %f %f %f %f\n",
            p.r.x, p.r.y, p.r.z, p.v.x, p.v.y, p.v.z );
    p.rk6(0.02); // nonrelativistic motion
  }

  return 0;
}
