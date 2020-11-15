// LoRPro (Long-Range Projectile motion simulator)
// To modeling the action of Coriolis Effect on long-rage projectiles motion//

// Citation:
// E. D. Guarin and N. Méndez-Hincapié. Modeling Coriolis Effect on long-range projectiles motion. Revista de Enseñanza de la Física, 28(1):73–82, 2016

#include<iostream>
#include<fstream>         
#include<math.h>
#include<stdlib.h>

using namespace std;

// This code works whit the parameters for the "Schwerer Gustav" missile at the Earth

const double g = 9.8; //Gravity in m/s^2
const double w = 7.272205e-5;//Angular speed of the Earth in rad/seg.
const double Dt = 0.001; //Temporal step for the numerical methods
const double R = 6.370e6;//Average radius of the Earth in meters.

//These parameters can be changed to observe this phenomenon in other planets for example.

class Missile{
    //Here, the variables and functions responsible for carrying out the calculations are defined
 public:
  double m,x,y,z,vx,vy,vz,fx,fy,fz,fco,L,t,wx,wy,wz;//Definition of parameters
  void start(double L0,double m0,double t0,double x0,double y0,double z0,double vx0,double vy0,double vz0);//Initial conditions definition
  void force(void);//Function that calculates the forces over the bullet
  void force2(void);//Function that calculates the forces over the bullet without considering the Coriolis and centrifugal forces
  void midpoint(double dt);//Nunerical method
  //Definition of the function to calculate position, velocity, and time for the movement: 
  double Getx(void){return x;};
  double Gety(void){return y;};
  double Getz(void){return z;};
  double Getvx(void){return vx;};
  double Getvy(void){return vy;};
  double Getvz(void){return vz;};
  double Gett(void){return t;};
  double Getf(void){return fz;};
};

void Missile::start(double L0,double m0,double t0, double x0,double y0,double z0,double vx0, double vy0,double vz0) {
  L=L0;//L -> Latitude angle where missile is located by the user
  m=m0;//Missile mass (4800 Kg.)
  //Initial conditions:
  x=x0;
  y=y0;
  z=z0;
  vx=vx0;
  vy=vy0;
  vz=vz0;
  t=t0;
}

//Coriolis Force -> F=-2m (w x v)
//w -> angular velocity of the planet
//Centrifugal force -> F=m*w x (w x r), r -> circumference radius described by the missile due to Earth´s rotation

void Missile::force(void) {//It calculates the total force over the missile considering Coriolis and centrifugal forces.
    double rx,ry,rz;
    double c=1.29*0.8*0.2*0.5; //Friction constant
    double v=sqrt((vx*vx)+(vy*vy)+(vz*vz));//Total velocity of the bullet
    wz=w*sin(L);
    wy=w*cos(L);
    wx=0;
    rx=x;
    ry=y;
    rz=R+z;
    // Cartesian components of the total force:
    fx=(-2*m*(wy*vz-wz*vy))+(m*(-wy*wy*rx-wz*wz*rx))-(c*v*vx);

    fy=(-2*m*(wz*vx))+(m*(wz*wy*rz-wz*wz*ry))-(c*v*vy);

    fz=(-2*m*(-wy*vx))-(m*g)+(m*(-wy*wy*rz+wz*wy*ry))-(c*v*vz);
    
    //Total Coriolis force:
    fco=sqrt(pow(2*m*(wy*vz-wz*vy),2)+pow(2*m*(wz*vx-wx*vz),2)+pow(2*m*(wx*vy-wy*vx),2));
}

void Missile::force2(void) {//It calculates the total force over the missile without considering Coriolis and centrifugal forces.
    double rx, ry, rz;
    double c=1.29*0.8*0.2*0.5;//Friction constant
    double v=sqrt((vx*vx)+(vy*vy)+(vz*vz));//Net velocity of the bullet
    // Cartesian components of the total force:
    fx=-(c*v*vx);
    fy=-(c*v*vy);
    fz=-(m*g)-(c*v*vz);
}

void Missile::midpoint(double dt) {//Using midpoint method:
  x=x+(dt*vx)+(0.5*dt*dt*(fx/m));
  y=y+(dt*vy)+(0.5*dt*dt*(fy/m));
  z=z+(dt*vz)+(0.5*dt*dt*(fz/m));
  vx=vx+(dt*(fx/m));
  vy=vy+(dt*(fy/m));
  vz=vz+(dt*(fz/m));
  t=t+dt;
}
int main () {
    //missile1.txt save the data for the situation with coriolis and centrifugal forces
    //missile2.txt save the data for the situation without coriolis and centrifugal forces
  ofstream f("missile1.txt"),h("missile2.txt"),c("command");
  Missile Gustav1;
  Missile Gustav2;
  double a,b,vel,ori,vx,vy,vz,delta,reach,t0,x0,y0,z0,vx0,vy0,vz0, Dtt, f0;
  cout<<"Please, enter the following initial conditions (use only numbers):"<<endl;
  cout<<" "<<endl;
  cout<<"Elevation angle for the shooting in degrees (min=0, max=90)= ";
  cin>>b;
  b=b*M_PI/180;
  cout<<" "<<endl;
  cout<<"Azimuthal angle in degrees (N=90, E=0, S=-90, W=180) = ";
  cin>>ori;
  ori=ori*M_PI/180;
  cout<<" "<<endl;
  cout<<"Initial velocity (m/s)= ";
  cin>>vel;
  cout<<" "<<endl;
  cout<<"Latitude (degrees)= ";
  cin>>a;
  a=a*M_PI/180;
  vx=vel*cos(b)*cos(ori);
  vy=vel*cos(b)*sin(ori);
  vz=vel*sin(b);
  cout<<" "<<endl;
  cout<<"calculating..."<<endl;

  //Here, the code calculates the data with the action of Coriolis and centrifugal force
  Gustav1.start(a,4800,0,0,0,2,vx,vy,vz);
  //Remember that the data in brackets are (latitude, mass, initial time, x position, y position, z position, x velocity, y velocity, z velocity)
  f<<Gustav1.Gett()<<"\t"<<Gustav1.Getx()<<"\t"<<Gustav1.Gety()<<"\t"<<Gustav1.Getz()<<endl;
  //Using the functions defined before to obtain the data for the movement:
  while (Gustav1.Getz()>=0) {
        x0=Gustav1.Getx();
        y0=Gustav1.Gety();
        z0=Gustav1.Getz();
        vx0=Gustav1.Getvx();
        vy0=Gustav1.Getvy();
        vz0=Gustav1.Getvz();
        t0=Gustav1.Gett();
        f0=Gustav1.Getf();
        Gustav1.force();
        Gustav1.midpoint(Dt);
        if(Gustav1.Getz()>0){
            f<<Gustav1.Gett()<<"\t"<<Gustav1.Getx()<<"\t"<<Gustav1.Gety()<<"\t"<<Gustav1.Getz()<<endl;
        }
        else {
            Dtt=(-vz0-sqrt(vz0*vz0-2*f0*z0/4800))/(f0/4800);
            Gustav1.start(a,4800,t0,x0,y0,z0,vx0,vy0,vz0);
            Gustav1.force();
            Gustav1.midpoint(Dtt);
            f<<Gustav1.Gett()<<"\t"<<Gustav1.Getx()<<"\t"<<Gustav1.Gety()<<"\t"<<Gustav1.Getz()<<endl;
            break;
        }

  }

  //Here, the code calculates the data without the action of Coriolis and centrifugal forces
  Gustav2.start(a,4800,0,0,0,2,vx,vy,vz);
  h<<Gustav2.Gett()<<"\t"<<Gustav2.Getx()<<"\t"<<Gustav2.Gety()<<"\t"<<Gustav2.Getz()<<endl;
  while (Gustav2.Getz() >=0) {
        x0=Gustav2.Getx();
        y0=Gustav2.Gety();
        z0=Gustav2.Getz();
        vx0=Gustav2.Getvx();
        vy0=Gustav2.Getvy();
        vz0=Gustav2.Getvz();
        t0=Gustav2.Gett();
        f0=Gustav2.Getf();
        Gustav2.force2();
        Gustav2.midpoint(Dt);
        if(Gustav2.Getz()>=0)
        h<<Gustav2.Gett()<<"\t"<<Gustav2.Getx()<<"\t"<<Gustav2.Gety()<<"\t"<<Gustav2.Getz()<<endl;
        else {
            Dtt=(-vz0-sqrt(vz0*vz0-2*f0*z0/4800))/(f0/4800);
            Gustav2.start(a,4800,t0,x0,y0,z0,vx0,vy0,vz0);
            Gustav2.force();
            Gustav2.midpoint(Dtt);
            h<<Gustav2.Gett()<<"\t"<<Gustav2.Getx()<<"\t"<<Gustav2.Gety()<<"\t"<<Gustav2.Getz()<<endl;
            break;
        }
  }

  //delta gives the horizontal deviation of the bullet with respect to the trajectory without the action of Coriolis and centrifugal forces
  delta=sqrt(pow(Gustav1.Getx()-Gustav2.Getx(),2)+pow(Gustav1.Gety()-Gustav2.Gety(),2));
  //Horizontal range of the missile:
  reach=sqrt(pow(Gustav1.Getx(),2)+pow(Gustav1.Gety(),2));
  f.close();
  h.close();
  cout.precision(7);

  c<<"set xlabel 'X [m]' \n set ylabel 'Y [m]' \n set zlabel 'Z [m]' \n  splot 'missile1.txt' u 2:3:4 w l \n pause -1 \n";

  c<<"set xlabel 'Deviation [m]' \n set ylabel 'Height [m]' \n plot 'missile1.txt' u 2:4 w l \n pause mouse \n ";

  c<<"set xlabel 'Horizontal range [m]' \n set ylabel 'Height [m]' \n unset key \n plot 'missile1.txt' u 3:4 w l \n pause mouse \n ";

  c<<"set xlabel 'Deviation [m]' \n set ylabel 'Horizontal range [m]' \n plot 'missile1.txt' u 2:3 w l, 'missile2.txt' u 2:3 w l \n pause mouse \n ";

  cout<<" "<<endl;
  cout<<" "<<endl;

    //Pop up window with additional information
  if(a>0)cout<<"The launch took place in the northern hemisphere."<<endl;
  if(a<0)cout<<"The launch took place in the southern hemisphere."<<endl;
  if(a==0)cout<<"The launch took place in the Ecuador line"<<endl;
  cout<<" "<<endl;
  cout<<"The horizontal range was "<<reach<<" meters"<<endl;
  cout<<" "<<endl;
  cout<<"The horizontal deviation was "<<delta<<" meters"<<endl;
  cout<<" "<<endl;
  cout<<"The flight time was: "<<Gustav1.Gett()<<" seconds"<<endl;
  cout<<" "<<endl;
  cout<<"WARNING: -x axis -> West, +x axis -> East"<<endl;
  cout<<"             -y axis-> South, +y axis -> North"<<endl;
  cout<<" "<<endl;
  
  c.close();
  system("pgnuplot<command");
  system("pause");
  return 0;
}
