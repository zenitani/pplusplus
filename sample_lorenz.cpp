#include <particle.h>
#include <stdio.h>

// ************* Lorenz attractor ************************************
// E. N. Lorenz, J. Atmos. Sci., 20, 130-141 (1963).
// doi:10.1175/1520-0469(1963)020<0130:DNF>2.0.CO;2
const double sigma = 10.0;
//const double a = sqrt(0.5);
const double b = 8.0/3.0;
const double r = 28.0;
// ************* Lorenz attractor ************************************

vector3 F( const vector3& _r, const vector3& _v,
           const  double& _t, const  double& _q ){
  return vector3( sigma*( -_v.x + _v.y ),
                  - (_r.x*_v.z + _v.x*_r.z) + r*_v.x - _v.y,
                  + (_r.x*_v.y + _v.x*_r.y) - b*_v.z );
}
int main()
{
  particle p;
  p.sett(0);
  p.setm(1); // unused
  p.setq(1); // unused
  p.setr(   0.0, 1.0, 0.0 );
  p.setv( sigma*( -p.r.x+p.r.y ),
          - p.r.x*p.r.z + r*p.r.x - p.r.y,
          + p.r.x*p.r.z - b*p.r.z );

  // marching in time
  for( int i=0;p.gett()<=60.0; i++ ){
    printf( "%f %f %f %f %f %f\n",
            p.r.x, p.r.y, p.r.z, p.v.x, p.v.y, p.v.z );
    p.rk6(0.01); // nonrelativistic motion
  }

  return 0;
}
