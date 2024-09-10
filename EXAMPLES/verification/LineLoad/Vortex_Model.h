#ifndef VORTEX_MODEL_H
#define VORTEX_MODEL_H

#include <vector>
#include <limits>
#include <math.h>
#include <stdio.h>
#include <iostream>

class VortexModel {
public:

VortexModel(){ }

void define_tor(double Um_in, double Sm_in, double rm_in, double z_in){
  Um = Um_in;
  S = Sm_in;
  rm = rm_in;
  zm = z_in;
  
};

void position(double xw, double yw, double zw, double &vxw, double &vyw, double &vzw ,std::string type){
 double vxf, vyf, vzf, rhof;
 if (type =="Baker"){    
   get_fluid_velocity_and_density(vxf,vyf,vzf,rhof,xw,yw,zw);
 }

else if (type =="Rankine"){    
   get_Rankine(vxf,vyf,vzf,rhof,xw,yw,zw);
 }

else if (type =="Bjerknes"){    
   get_Bjerknes(vxf,vyf,vzf,rhof,xw,yw,zw);
 }

 else if (type =="Vatistas"){    
   get_Vatistas(vxf,vyf,vzf,rhof,xw,yw,zw);
 }

else if (type == "Fujita") {
  get_Fujita(vxf,vyf,vzf,rhof,xw,yw,zw);
}

else if (type == "Burgers_rott") {
  get_Burgers_rott(vxf,vyf,vzf,rhof,xw,yw,zw);
}

 else if (type == "Sullivan") {
  get_Sullivan(vxf,vyf,vzf,rhof,xw,yw,zw);
 }

 else {
  std::cout<< "Enter the correct name this is automatically running Baker Model";
  get_fluid_velocity_and_density(vxf,vyf,vzf,rhof,xw,yw,zw);

 }
}

void get_fluid_velocity_and_density(double& vxf, double& vyf, double& vzf, double& rhof,
			  	      double xp, double yp, double zp) {
    // Define reference values and constants for use in dimensionless evaluations
    double K = S*(2.0/std::log(2.0));
    double gamma = 2.0;
    double delta = zm/rm;

    // Define the center of the vortex
    double xc = 10.0; // [m]
    double yc = 0.0; // [m]
    double zc = 0.0; // [m]

    // Compute normalized radial and height coordinates
    double rp = std::sqrt((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc));
    double inv_rp = 1.0/(rp+std::numeric_limits<double>::min());
    double cosp = (xp-xc)*inv_rp;
    double sinp = (yp-yc)*inv_rp;
    double rbar = rp/rm;
    double zbar = (zp-zc)/zm;

    // Compute normalized radial, tangential, and vertical velocity of the vortex
    double one_rbar2 = 1.0 + rbar*rbar;
    double one_zbar2 = 1.0 + zbar*zbar;
    double log_one_zbar2 = std::log(one_zbar2);
    double Ubar = -4.0*rbar*zbar/(one_rbar2*one_zbar2);
    double Vbar = K*std::pow(rbar,gamma-1.0)*std::pow(log_one_zbar2,0.5*gamma)/std::pow(one_rbar2,0.5*gamma);
    double Wbar = 4.0*delta*log_one_zbar2/(one_rbar2*one_rbar2);
    double Pbar = 8.0*std::pow(rbar,2)*zbar/std::pow(one_rbar2*one_zbar2,2)-4.15*S*S*std::pow(log_one_zbar2,2)/one_rbar2 -
                  4*log_one_zbar2*one_zbar2/std::pow(one_rbar2*one_zbar2,2);
                  

    // Compute the x,y,z components of the fluid velocity
    vxf = Um*(Ubar*cosp-Vbar*sinp);
    vyf = Um*(Ubar*sinp+Vbar*cosp);
    vzf = Um*Wbar;

    // Assume the density is constant
    rhof = 1.293; // [kg/m^3] fluid density (of air at STP)
      
  } // get_fluid_velocity_and_density()

  //Ranking Model 
  //Rankine, W.J.M., 1882. A Manual of Applied Physics, tenth ed. Charles Griff and Co.































void get_Rankine (double& vxf, double& vyf, double& vzf, double& rhof,
			  	      double xp, double yp, double zp) {
                  vzf =0;
                  double vtf;
                  rhof = 1.293; // [kg/m^3] fluid density (of air at STP)
                  double Um = 100.0;
                  double rc = 0.005;
                  double rm =0.01;
                  double E = 1 ; //decay index
                  // Define the center of the vortex
                  double xc = 10.0; // [m]
                  double yc = 0.0; // [m]
                  double zc = 0.0; // [m]

                  // Compute normalized radial and height coordinates
                  double rp = std::sqrt((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc));
                  double thita = std::atan2((yp-yc),(xp-xc));
                  double rbar = rp/rc;
                  if(rp<rm){
                    if(rp<rc){
                      vtf = Um*rbar;
                      double delP =rhof*std::pow(Um,2)/2*std::pow(rbar,2);
                    }
                    else{
                      vtf = std::pow(1/rbar,E);
                      double delP =rhof*std::pow(Um,2)/2+rhof*std::pow(Um,2)/(2*E)*(1-std::pow(1/rbar,2*E));
                    }
                  }
                  else{
                    vtf = 0;
                  }
                  vxf = vtf*std::cos(thita);
                  vyf = vtf*std::sin(thita);
                }


   //Bjerknes Model
   // https://journals.ametsoc.org/view/journals/mwre/49/1/1520-0493_1921_49_1_tmottz_2_0_co_2.xml
  void get_Bjerknes (double& vxf, double& vyf, double& vzf, double& rhof,
			  	      double xp, double yp, double zp) {
                  vzf =0;
                  double vtf;
                  double Um = 100.0;
                  double rc =0.01;
                  double E = 1 ; //decay index
                  // Define the center of the vortex
                  double xc = 10.0; // [m]
                  double yc = 0.0; // [m]
                  double zc = 0.0; // [m]
                  double rhof = 1.293; 
                  // Compute normalized radial and height coordinates
                  double rp = std::sqrt((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc));
                  double rbar = rp/rc;
                  double thita = std::atan2((yp-yc),(xp-xc));
                  vtf = 2*Um *rbar/(1+std::pow(rbar,2));
                  vxf = vtf*std::cos(thita);
                  vyf = vtf*std::sin(thita);

                  double del_P = 2*rhof*std::pow(Um,2)*rbar/(1+std::pow(rbar,2));
                }

  void get_Vatistas(double& vxf, double& vyf, double& vzf, double& rhof,
			  	      double xp, double yp, double zp) {
                  double Um = 100.0;
                  double rc =0.01;
                  double beta = 1 ; //decay index
                  // Define the center of the vortex
                  double ve = 1.5*std::pow(10,-5);//eddy viscosity m^2/s
                  double xc = 10.0; // [m]
                  double yc = 0.0; // [m]
                  double zc = 0.0; // [m]
                  double rhof = 1.293;
                  double beta = 2; // Power law index should be more than 1
                  // Compute normalized radial and height coordinates
                  double rp = std::sqrt((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc));
                  double thita = std::atan2((yp-yc),(xp-xc));
                  double rbar = rp/rc;
                  double V = Um*rbar/std::pow(1+std::pow(rbar,2*beta),1/beta);
                  double U =  -2*(beta+1)*(ve/rc)*std::pow(rbar,2*beta-1)/(1+std::pow(rbar,2*beta));
                  double W = 4*beta*(beta+1)*(ve/rc)*(zp/rc)*std::pow(rbar,2*(beta-1))/std::pow((1+std::pow(rbar,2*beta)),2);
                  vxf = U*std::cos(thita)-V*std::sin(thita);
                  vyf = U*std::sin(thita)+V*std::cos(thita);
                  vzf = W;
                }

//Fujita Model 
//https://scholar.google.com/scholar_lookup?title=Workbook%20of%20Tornadoes%20and%20High%20Winds%20for%20Engineering%20Applications&author=T.T.%20Fujita&publication_year=1978
  void get_Fujita(double& vxf, double& vyf, double& vzf, double& rhof,
			  	      double xp, double yp, double zp) {
                  double Um = 100.0;
                  double ro = 0.1;
                  double n = 0.9-0.7*exp(-0.0005*ro);
                  double rn = n*ro; 
                  double Hi = 0.55*(1-n*n)*ro;
                  double xc = 10.0; // [m]
                  double yc = 0.0; // [m]
                  double zc = 0.0; // [m]
                  double rhof = 1.293;
                  double ko = 1/6;
                  double k = 0.003;
                  double j =0.03;
                  double Fr, Fz;
                  double Tana, Tanao;
                  double Am = {0.75}, Bm ={0.0217};
                  double rp = std::sqrt((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc));
                  double thita = std::atan2((yp-yc),(xp-xc));
                  double V,U,W,Wo, delP;
                  if(zp<Hi){
                    Fz = std::pow(zc/Hi,ko);
                    Tanao = -Am*(1-std::pow(zp/Hi,1.5)); 
                    Wo = 0.059*Am*(16*std::pow(zp/Hi,7/6)-7*(pow(zp/Hi,8/3))) ;
                                       
                  }
                  else{
                    Fz = exp(-k*((zp/Hi)-1));
                    Tanao = Bm *(1-exp(-k*(zp/Hi-1)));
                    Wo = 0.55*Bm/k*(2*exp(-k*((zp/Hi)-1))-exp(-2*k*(zp/Hi-1)));
                  }

                  if (rp<rn){
                    Fr = rp/ro;
                    Tana = 0;
                    W =0;
                  }
                  else if(rn<=rp<ro){
                    Fr = rp/ro;
                    Tana = Tanao/(1-n*n)*(1-std::pow(ro/rp*n,2));
                    W = Um*Wo;
                    delP = rhof*std::pow(Fz*Um*rp/ro,2)/2;
                  }
                  else{
                    Fr = ro/rp;
                    Tana = Tanao;
                    W = 0;
                    delP = rhof*std::pow(Fz*Um,2)-rhof*std::pow(Fz*Um*rp/ro,2)/2;;
                  }       
                  U = Um*Fr*Fz;
                  V = U*Tana;
                  vxf = U*std::cos(thita)-V*std::sin(thita);
                  vyf = U*std::sin(thita)+V*std::cos(thita);
                  vzf = W;
                }

// Burgers model 
//https://www.sciencedirect.com/science/article/pii/S0065215608701005
  void get_Burgers_rott(double& vxf, double& vyf, double& vzf, double& rhof,
			  	      double xp, double yp, double zp) {
                  double U,V,W;
                  rhof = 1.293; // [kg/m^3] fluid density (of air at STP)
                  double Um = 100.0;
                  // Define the center of the vortex
                  double xc = 10.0; // [m]
                  double yc = 0.0; // [m]
                  double zc = 0.0; // [m]
                   double rp = std::sqrt((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc));
                  double thita = std::atan2((yp-yc),(xp-xc));
                  double ve = 1.5*std::pow(10,-5);//eddy viscosity m^2/s
                  double Kbr1 = {1.26} , Krb2 = {0.72};
                  double a = 1.26*2*ve/(rp*rp); //Velocity Gradiant 
                  double rc = std::sqrt(2*ve*Kbr1/a);
                  V = Um*rc/(Krb2*rp)*(1-exp(-Kbr1*std::pow(rp/rc,2))); //Tangential Component 
                  U =-a*rp;
                  W = 2*a*zp;
                  vxf = U*std::cos(thita)-V*std::sin(thita);
                  vyf = U*std::sin(thita)+V*std::cos(thita);
                  vzf = W;
                }

  void get_Sullivan (double& vxf, double& vyf, double& vzf, double& rhof,
			  	      double xp, double yp, double zp) {
                           double U,V,W;
                  double U,V,W,H={zp};
                  rhof = 1.293; // [kg/m^3] fluid density (of air at STP)
                  double Um = 100.0;
                  // Define the center of the vortex
                  double xc = 10.0; // [m]
                  double yc = 0.0; // [m]
                  double zc = 0.0; // [m]
                   double rp = std::sqrt((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc));
                  double thita = std::atan2((yp-yc),(xp-xc));
                  double ve = 1.5*std::pow(10,-5);//eddy viscosity m^2/s
                  double Ks1 = {6.24} , Ks2 = {0.88}, Hinf = {37.9};
                  double a = 1.26*2*ve/(rp*rp); //Velocity Gradiant 
                  double rc = std::sqrt(2*ve*Ks1/a);
                  V = Um/Ks2*rc/rp*H*(Ks1*(rp/rc)*(rp/rc))/Hinf ; //Tangential Component 
                  U =-a*rp +6*ve/rp*(1-exp(-a*rp*rp/2/ve));
                  W = 2*a*zp*(1-3*exp(-a*rp*rp/2/ve));
                  vxf = U*std::cos(thita)-V*std::sin(thita);
                  vyf = U*std::sin(thita)+V*std::cos(thita);
                  vzf = W;
                }



private:
 double Um, S, rm, zm;


};

#endif