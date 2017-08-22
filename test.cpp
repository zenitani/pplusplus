#include <particle.h>
#include <stdio.h>

vector3 E( const vector3& _r ){
  vector3 _E( 0.0, 0.5, 0.0 );
  return _E;
}
vector3 B( const vector3& _r ){
  vector3 _B( 0.0, 0.0, 1.0 );
  return _B;
}
vector3 F( const vector3& _r, const vector3& _v,
           const  double& _t, const  double& _q ){
  return _q * ( E(_r) + _v * B(_r) );
}
int main()
{
  particle p;
  p.sett(0);
  p.setm(1);
  p.setq(1);
  //  p.setq(-1);
  p.setr(0.0,0.0,0.0);
  p.setv(0.3,0.0,0.1);

  // rewinding
  for( int i=0;p.gett()>-10; i++ ){
//     printf( "%f %f %f %f %f %f\n",
// 	    p.r.x, p.r.y, p.r.z,
// 	    p.v.x, p.v.y, p.v.z );
    p.rk4(-0.2); // nonrelativistic motion
  }
  // marching in time
  for( int i=0;p.gett()<50; i++ ){
    printf( "%f %f %f %f %f %f\n",
	    p.r.x, p.r.y, p.r.z,
	    p.v.x, p.v.y, p.v.z );
    p.rk4(0.2); // nonrelativistic motion
  }

  return 0;
}
