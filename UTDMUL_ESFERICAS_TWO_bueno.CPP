 // Calcula la atenuacion en funcion del numero de aristas para N wedges con el
// mismo ·ngulo interior.

#include "cab.hpp"
#include <stdio.h>
#include <Complex.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "Complex.hpp"

#define pi 3.14159265358979323846
#define ONE Complex(1.0,0.0)
#define ZERO Complex(0.0,0.0)
#define N 11     /* N-1 edificios. Hay que poner 8 para 6 edificios completos y Rx a w del ˙ltimo edificio*/
#define Ngamma 1000
//#define frecuencia 1000e6
#define frecuencia 9.37e9
#define LAMBDA (3e8/frecuencia)
//#define LAMBDA 0.55


int knife=0;
double valor=5.0;
double frec=9.37e9;
fcomplex Coef_difrac(int sft, int perf_cond, fcomplex permirel,
							double angulo_int, double bta,double btacero,
							double phi, double phiprima, double L);

//fcomplex planas(int q, double y0);
//fcomplex planas2(int q, double y0);
//fcomplex Coef_refl(int soft_hard, int perf_cond, fcomplex permirel, double alpha);

main()

{
  FILE *fsal;


  fcomplex AnUTD[N],AnUTDaux[N],AnUTD2[N],A,B,C,D,E,F,G,H,Em[N];  
  fcomplex AnUTDT[N],pp;
 // fcomplex Ancte,peso,aux3,Ef,Edir,Edif,Et,E1,E1_dir,E2,E2_dir,E3,Eini,var1,var2,var3,Edif_1,Edif_2,Edif_3,Edif_4,Edif_5,Edif_6,Edif_7,Edif_8,Edif_9,Edif_10,Edif_11,Edif_12,Edif_13,Edif_14,Edif_15,Edif_16,Edif_17,Edif_18,Edif_19,Edif_20,Edif_21,Edif_22,Edif_23,Edif_24,Edif_25,Edif_26,Edif_27,Edif_28,Edif_29,Edif_30,Eref_1,Eref_2,Edif1,Edif2,Eref,Edif3,Edif4,Eref1,Eref2,Edif_doble1,Edif_doble2,Edif_doble3,Edif_doble4,Edif_doble5,Edif_doble6,Edif_triple1,Edif_triple2,Edif_triple3,Edif_triple4,Edif_cuadruple;
  fcomplex Ancte,peso,aux3,Ef,Et,Edif,var1,var2,E1,E2,E3,Eini,Edif_doble1,Edif_doble2,Edif_doble3,Edif_triple1,Edif_triple2,Edif_triple3,Edif_triple4,Edif_triple5,Edif_cuadruple,Edif_cuadruple1,Edif_cuadruple2,Edif_cuadruple3,Edif_cuadruple4,Edif_cuadruple5,Edif_cuadruple6,Edif_1,Edif_2,Edif_3,Edif_4,Edif_5,Edif_6,Edif_7,Edif_8,Edif_9,Edif_10,Edif_11,Edif_12,Edif_13,Edif_14,Edif_15,Edif_16,Edif_17,Edif_18;
  fcomplex CoefD,CoefD2,CoefD3,CoefD4,SumoDirecto,Ddd;
  double R[2*N+1],R2[2*N+1];
  double suma,alpha,alfa,alphan1,alpham2,alpha2,alfa0,alfafin,v,w,phit,phi,phiprima,phi1,phi1prima,phi2,phi2prima,phi_2,phiprima_2,Lw2,Lw,Lw3,Lw4,L1,L2,aux,aux2,s2,betha,alphaprima,mu,muprima,bethaprima,rho,rhoprima,omega,tau,tauprima,theta;
  double v1r[orden2],v1i[orden2],v2r[orden2],v2i[orden2],dd,d2,d3;
  float y0,fase,AA,ss,RR,Dist,E0,d,alt,alt2,Rt,z,R0,R1,R22,R3,R4,R5,R6,R7,R0prima,R1prima,R3prima,R4prima,R5prima,A1,t,Rr,h,l;
  int i,n,m,p,q,j,r,rmax,Ei;
  int SOFT=0;

 // fsal=fopen("settled_field_HARD_v28_w22_922.dat","w");

  fsal=fopen("prueba.dat","w");
 // fsal=fopen("prueba.dat","w");
/////////////////////////////////  w=22;   //Separacion entre edificios en longitudes de onda
/////////////////////////////////  v=28;    //Anchura de edificios en longitudes de onda
///////////////////////  w=0.4/LAMBDA;
//////////////////////  v=0.1/LAMBDA;
//  w=41.33;
//  v=4.133;
  phit=45.0*pi/180;	// Semiangulo del wedge
  phi=1.5*pi;
  phi_2=pi;
//    phi=1.5*pi-0.000001*phit;
//    phi_2=1.5*pi-0.999999*phit;
//    phi=1.5*pi+0.0012*phit;
//  phi_2=pi+0.001*phit;

//  dd=10.0*LAMBDA;
//  dd=105.25*LAMBDA;
///////////////////////////////////////////////////  dd=0.162;
///////////////////////////////////y0=0.0075;
//  dd=20.666*LAMBDA;
////////////////////////////////////////////////dd=30.0;
  Ei=1;

  //alpha=1.2*pi/180; //·ngulo de incidencia

//for (w=20.0;w<201.0;w=w+1.0) {

// for (v=20.0;v<201.0;v=v+1.0) {



//sin obstaculos
//E0=0.002555;
//for (Dist=0.05;Dist<2.5;Dist=Dist+0.05){
//      Ef=RCmul(E0/Dist,Complex(cos(2*pi*Dist/LAMBDA),-sin(2*pi*Dist/LAMBDA)));
//	fprintf(fsal,"%f %f\n", Dist, 20*log10(Cabs(Ef)));
//}

/*
//con 1 obstaculo
E0=0.002555;
dd=0.6;
alt=0.032;
v=0.03;
	
	//RAYO DIRECTO 
	for (w=0.01;w<1.87;w=w+0.05){
	Rt=dd+v+w;
	Edir=RCmul(E0/Rt,Complex(cos(2*pi*Rt/LAMBDA),-sin(2*pi*Rt/LAMBDA)));
	
	//DIFRACCION DOBLE 
	//E(0)
	R0=sqrt(pow(alt,2)+pow(dd,2));
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	alpha=atan(alt/dd);
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2 + alpha;
	Lw= ((R0/LAMBDA)*(v/LAMBDA))/((R0/LAMBDA)+(v/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((v/LAMBDA)*((R0/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA)));          
	
	//Edif doble total 
	R1=sqrt(pow(alt,2)+pow((dd+v),2));
	alphaprima= atan(alt/(w));
	z=sqrt(pow(alt,2)+pow(w,2));
	betha= atan(alt/(dd+v));
	phi=pi-alphaprima;
	phit=45.0*pi/180;
	phiprima=0;         
	Lw=(((R0+v)/LAMBDA)*(z/LAMBDA))/(((R0+v)/LAMBDA)+(z/LAMBDA));       
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v)/LAMBDA)/((z/LAMBDA)*(((R0+v)/LAMBDA)+(z/LAMBDA))));       
	var2=RCmul(A1,var1);
	Edif=Cmul(var2,Complex(cos(2*pi*z/LAMBDA),-sin(2*pi*z/LAMBDA)));

	//DIFRACCIONES SIMPLES
	//Edif1 
	t=sqrt(pow(alt,2)+pow((w+v),2));
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));
	
	mu=atan(alt/(w+v));
	phi=3*pi/2-mu;
	phit=45.0*pi/180;
	phiprima=pi/2 + alpha;
	Lw= ((R0/LAMBDA)*(t/LAMBDA))/((R0/LAMBDA)+(t/LAMBDA));

	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((t/LAMBDA)*((R0/LAMBDA)+(t/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif1=Cmul(var2,Complex(cos(2*pi*t/LAMBDA),-sin(2*pi*t/LAMBDA)));

	//Edif2 
	Eini=RCmul(E0/R1,Complex(cos(2*pi*R1/LAMBDA),-sin(2*pi*R1/LAMBDA)));
	
	phi=pi-alphaprima;
	phit=45.0*pi/180;
	phiprima=betha;
	Lw= ((R1/LAMBDA)*(z/LAMBDA))/((R1/LAMBDA)+(z/LAMBDA));

	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R1/LAMBDA)/((z/LAMBDA)*((R1/LAMBDA)+(z/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif2=Cmul(var2,Complex(cos(2*pi*z/LAMBDA),-sin(2*pi*z/LAMBDA)));
	

	//RAYOS REFLEJADOS 
	if(w>=0.57 && w<=0.63){
			Rr=sqrt(pow(alt,2)+pow(((dd+v+w)/2),2));
			Eref=RCmul(E0/(2*Rr),Complex(cos(2*pi*2*Rr/LAMBDA),-sin(2*pi*2*Rr/LAMBDA)));	
	}

	//E total
	Edif_1=Cadd(Edif1,Edif1);
	Edif_2=Cadd(Edif2,Edif2);
	Edif_3=Cadd(Edif_1,Edif_2);
	Edif_4=Cadd(Edir,Edif_3);
	Edif_5=Cadd(Edif,Edif);
	Edif_6=Cadd(Edif_4,Edif_5);
	Edif_7=Cadd(Eref,Eref);
	if(w>=0.57 && w<=0.63){
			Et=Cadd(Edif_6,Edif_7);		
		
	}
	else { Et=Edif_6; }
	
	fprintf(fsal,"%f %f\n", Rt, 20*log10(Cabs(Et)));
	//fprintf(fsal,"%f\n", 20*log10(Cabs(Et)));
	}


//2 obstaculos
E0=0.002555;
dd=0.6;
d2=0.8;
alt=0.032;
v=0.03;
d3=d2-dd-v;


	//RAYO DIRECTO
	for (w=0.01;w<1.67;w=w+0.05){
	Rt=d2+v+w;
	Edir=RCmul(E0/Rt,Complex(cos(2*pi*Rt/LAMBDA),-sin(2*pi*Rt/LAMBDA)));
	
	//DIFRACCIONES SIMPLES 
	//Edif1
	R0=sqrt(pow(alt,2)+pow(dd,2));
	t=sqrt(pow(alt,2)+pow((d2-dd+v+w),2));
	alpha=atan(alt/dd);
	mu=atan(alt/(d2-dd+v+w));

	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	phi=3*pi/2-mu;
	phit=45.0*pi/180;
	phiprima=pi/2 + alpha;
	Lw= ((R0/LAMBDA)*(t/LAMBDA))/((R0/LAMBDA)+(t/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((t/LAMBDA)*((R0/LAMBDA)+(t/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif1=Cmul(var2,Complex(cos(2*pi*t/LAMBDA),-sin(2*pi*t/LAMBDA)));

	//Edif2
	R1=sqrt(pow(alt,2)+pow((dd+v),2));
	z=sqrt(pow(alt,2)+pow((d2-dd+w),2));
	betha=atan(alt/(dd+v));
	rho=atan(alt/(d2-dd+w));

	Eini=RCmul(E0/R1,Complex(cos(2*pi*R1/LAMBDA),-sin(2*pi*R1/LAMBDA)));
	
	phi=pi-rho;
	phit=45.0*pi/180;
	phiprima=betha;
	Lw= ((R1/LAMBDA)*(z/LAMBDA))/((R1/LAMBDA)+(z/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R1/LAMBDA)/((z/LAMBDA)*((R1/LAMBDA)+(z/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif2=Cmul(var2,Complex(cos(2*pi*z/LAMBDA),-sin(2*pi*z/LAMBDA)));


	//Edif3
	R0prima=sqrt(pow(alt,2)+pow(d2,2));
	h=sqrt(pow(alt,2)+pow((v+w),2));
	alphaprima=atan(alt/d2);
	muprima=atan(alt/(v+w));

	Eini=RCmul(E0/R0prima,Complex(cos(2*pi*R0prima/LAMBDA),-sin(2*pi*R0prima/LAMBDA)));
	
	phi=3*pi/2-muprima;
	phit=45.0*pi/180;
	phiprima=pi/2+alphaprima;
	Lw= ((R0prima/LAMBDA)*(h/LAMBDA))/((R0prima/LAMBDA)+(h/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0prima/LAMBDA)/((h/LAMBDA)*((R0prima/LAMBDA)+(h/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif3=Cmul(var2,Complex(cos(2*pi*h/LAMBDA),-sin(2*pi*h/LAMBDA)));


	//Edif4
	R1prima=sqrt(pow(alt,2)+pow((d2+v),2));
	l=sqrt(pow(alt,2)+pow(w,2));
	bethaprima=atan(alt/(d2+v));
	rhoprima=atan(alt/w);

	Eini=RCmul(E0/R1prima,Complex(cos(2*pi*R1prima/LAMBDA),-sin(2*pi*R1prima/LAMBDA)));

	phi=pi-rhoprima;
	phit=45.0*pi/180;
	phiprima=bethaprima;
	Lw= ((R1prima/LAMBDA)*(l/LAMBDA))/((R1prima/LAMBDA)+(l/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R1prima/LAMBDA)/((l/LAMBDA)*((R1prima/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif4=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));


	//DIFRACCIONES DOBLES
	//Edif_doble1
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2 + alpha;
	Lw= ((R0/LAMBDA)*(v/LAMBDA))/((R0/LAMBDA)+(v/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((v/LAMBDA)*((R0/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA)));   
	
	//Edif doble total 
	phi=pi-rho;
	phit=45.0*pi/180;
	phiprima=0;           
	Lw=(((R0+v)/LAMBDA)*(z/LAMBDA))/(((R0+v)/LAMBDA)+(z/LAMBDA));          
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v)/LAMBDA)/((z/LAMBDA)*(((R0+v)/LAMBDA)+(z/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_doble1=Cmul(var2,Complex(cos(2*pi*z/LAMBDA),-sin(2*pi*z/LAMBDA)));


	//Edif_doble2
	//E(0)
	Eini=RCmul(E0/R1,Complex(cos(2*pi*R1/LAMBDA),-sin(2*pi*R1/LAMBDA)));

	//E(1)
	phi=pi;
	phit=45.0*pi/180;
	phiprima=betha;
	Lw= ((R1/LAMBDA)*(d3/LAMBDA))/((R1/LAMBDA)+(d3/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R1/LAMBDA)/((d3/LAMBDA)*((R1/LAMBDA)+(d3/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*d3/LAMBDA),-sin(2*pi*d3/LAMBDA)));                   
	
	//Edif doble total 
	phi=3*pi/2-muprima;
	phit=45.0*pi/180;
	phiprima=0;              
	Lw=(((R1+d3)/LAMBDA)*(h/LAMBDA))/(((R1+d3)/LAMBDA)+(h/LAMBDA));         
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R1+d3)/LAMBDA)/((h/LAMBDA)*(((R1+d3)/LAMBDA)+(h/LAMBDA))));           
	var2=RCmul(A1,var1);
	Edif_doble2=Cmul(var2,Complex(cos(2*pi*h/LAMBDA),-sin(2*pi*h/LAMBDA)));


	//Edif_doble3
	//E(0)
	Eini=RCmul(E0/R0prima,Complex(cos(2*pi*R0prima/LAMBDA),-sin(2*pi*R0prima/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2+alphaprima;
	Lw= ((R0prima/LAMBDA)*(v/LAMBDA))/((R0prima/LAMBDA)+(v/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0prima/LAMBDA)/((v/LAMBDA)*((R0prima/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA)));             
	
	//Edif doble total 
	phi=pi-rhoprima;
	phit=45.0*pi/180;
	phiprima=0;             
	Lw=(((R0prima+v)/LAMBDA)*(l/LAMBDA))/(((R0prima+v)/LAMBDA)+(l/LAMBDA));           
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0prima+v)/LAMBDA)/((l/LAMBDA)*(((R0prima+v)/LAMBDA)+(l/LAMBDA))));                     
	var2=RCmul(A1,var1);
	Edif_doble3=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));


	//Edif_doble4
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2+alpha;
	Lw= ((R0/LAMBDA)*((v/LAMBDA)+(d3/LAMBDA)))/((R0/LAMBDA)+(v/LAMBDA)+(d3/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/(((v/LAMBDA)+(d3/LAMBDA))*((R0/LAMBDA)+(v/LAMBDA)+(d3/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*(v/LAMBDA+d3/LAMBDA)),-sin(2*pi*(v/LAMBDA+d3/LAMBDA))));   
	
	//Edif doble total 
	phi=3*pi/2-muprima;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+v+d3)/LAMBDA)*(h/LAMBDA))/(((R0+v+d3)/LAMBDA)+(h/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v+d3)/LAMBDA)/((h/LAMBDA)*(((R0+v+d3)/LAMBDA)+(h/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_doble4=Cmul(var2,Complex(cos(2*pi*h/LAMBDA),-sin(2*pi*h/LAMBDA)));


	//Edif_doble5
	//E(0)
	Eini=RCmul(E0/R1,Complex(cos(2*pi*R1/LAMBDA),-sin(2*pi*R1/LAMBDA)));

	//E(1)
	phi=pi;
	phit=45.0*pi/180;
	phiprima=betha;
	Lw= ((R1/LAMBDA)*((v/LAMBDA)+(d3/LAMBDA)))/((R1/LAMBDA)+(v/LAMBDA)+(d3/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R1/LAMBDA)/(((v/LAMBDA)+(d3/LAMBDA))*((R1/LAMBDA)+(v/LAMBDA)+(d3/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*(v/LAMBDA+d3/LAMBDA)),-sin(2*pi*(v/LAMBDA+d3/LAMBDA))));   
	
	//Edif doble total 
	phi=pi-rhoprima;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R1+d3+v)/LAMBDA)*(l/LAMBDA))/(((R1+d3+v)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R1+d3+v)/LAMBDA)/((l/LAMBDA)*(((R1+d3+v)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_doble5=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));


	//Edif_doble6
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2+alpha;
	Lw= ((R0/LAMBDA)*(2*(v/LAMBDA)+(d3/LAMBDA)))/((R0/LAMBDA)+2*(v/LAMBDA)+(d3/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((2*(v/LAMBDA)+(d3/LAMBDA))*((R0/LAMBDA)+2*(v/LAMBDA)+(d3/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*(2*(v/LAMBDA)+(d3/LAMBDA))),-sin(2*pi*(2*(v/LAMBDA)+(d3/LAMBDA)))));   
	
	//Edif doble total 
	phi=pi-rhoprima;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+2*v+d3)/LAMBDA)*(l/LAMBDA))/(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+2*v+d3)/LAMBDA)/((l/LAMBDA)*(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_doble6=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	//DIFRACCIONES TRIPLES 
	//Edif_triple1
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2 + alpha;
	Lw= ((R0/LAMBDA)*(v/LAMBDA))/((R0/LAMBDA)+(v/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((v/LAMBDA)*((R0/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	  
	//E(2)
	phi=pi;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+v)/LAMBDA)*(d3/LAMBDA))/(((R0+v)/LAMBDA)+(d3/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v)/LAMBDA)/((d3/LAMBDA)*(((R0+v)/LAMBDA)+(d3/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*d3/LAMBDA),-sin(2*pi*d3/LAMBDA))); 
	
	//Dif triple total
	phi=3*pi/2-muprima;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+d3+v)/LAMBDA)*(h/LAMBDA))/(((R0+d3+v)/LAMBDA)+(h/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R0+d3+v)/LAMBDA)/((h/LAMBDA)*(((R0+d3+v)/LAMBDA)+(h/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_triple1=Cmul(var2,Complex(cos(2*pi*h/LAMBDA),-sin(2*pi*h/LAMBDA)));

	//Edif_triple2
	//E(0)
	Eini=RCmul(E0/R1,Complex(cos(2*pi*R1/LAMBDA),-sin(2*pi*R1/LAMBDA)));

	//E(1)
	phi=pi;
	phit=45.0*pi/180;
	phiprima=betha;
	Lw= ((R1/LAMBDA)*(d3/LAMBDA))/((R1/LAMBDA)+(d3/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R1/LAMBDA)/((d3/LAMBDA)*((R1/LAMBDA)+(d3/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*d3/LAMBDA),-sin(2*pi*d3/LAMBDA)));
	
	//E(2)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R1+d3)/LAMBDA)*(v/LAMBDA))/(((R1+d3)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R1+d3)/LAMBDA)/((v/LAMBDA)*(((R1+d3)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	
	//Dif triple total 
	phi=pi-rhoprima;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R1+d3+v)/LAMBDA)*(l/LAMBDA))/(((R1+d3+v)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R1+d3+v)/LAMBDA)/((l/LAMBDA)*(((R1+d3+v)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_triple2=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	//Edif_triple3
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2 + alpha;
	Lw= ((R0/LAMBDA)*(v/LAMBDA))/((R0/LAMBDA)+(v/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((v/LAMBDA)*((R0/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	  
	//E(2)
	phi=pi;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+v)/LAMBDA)*((d3+v)/LAMBDA))/(((R0+v)/LAMBDA)+((d3+v)/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v)/LAMBDA)/(((d3+v)/LAMBDA)*(((R0+v)/LAMBDA)+((d3+v)/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*(d3+v)/LAMBDA),-sin(2*pi*(d3+v)/LAMBDA))); 
	
	//Dif triple total 
	phi=pi-rhoprima;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+2*v+d3)/LAMBDA)*(l/LAMBDA))/(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R0+2*v+d3)/LAMBDA)/((l/LAMBDA)*(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_triple3=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));


	//Edif_triple4
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2 + alpha;
	Lw= ((R0/LAMBDA)*((v+d3)/LAMBDA))/((R0/LAMBDA)+((v+d3)/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/(((v+d3)/LAMBDA)*((R0/LAMBDA)+((v+d3)/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*(v+d3)/LAMBDA),-sin(2*pi*(v+d3)/LAMBDA))); 
	  
	//E(2)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+v+d3)/LAMBDA)*(v/LAMBDA))/(((R0+v+d3)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v+d3)/LAMBDA)/((v/LAMBDA)*(((R0+v+d3)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	
	//Dif triple total 
	phi=pi-rhoprima;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+2*v+d3)/LAMBDA)*(l/LAMBDA))/(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R0+2*v+d3)/LAMBDA)/((l/LAMBDA)*(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_triple4=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	//DIFRACCION CUADRUPLE 
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2 + alpha;
	Lw= ((R0/LAMBDA)*(v/LAMBDA))/((R0/LAMBDA)+(v/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((v/LAMBDA)*((R0/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	  
	//E(2)
	phi=pi;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+v)/LAMBDA)*(d3/LAMBDA))/(((R0+v)/LAMBDA)+(d3/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v)/LAMBDA)/((d3/LAMBDA)*(((R0+v)/LAMBDA)+(d3/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*d3/LAMBDA),-sin(2*pi*d3/LAMBDA))); 
	
	//E(3)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+d3+v)/LAMBDA)*(v/LAMBDA))/(((R0+d3+v)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R0+d3+v)/LAMBDA)/((v/LAMBDA)*(((R0+d3+v)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E3=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA)));
	
	//Dif cuadruple total 
	phi=pi-rhoprima;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+2*v+d3)/LAMBDA)*(l/LAMBDA))/(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E3,CoefD);
	A1=sqrt(((R0+2*v+d3)/LAMBDA)/((l/LAMBDA)*(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_cuadruple=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	//RAYOS REFLEJADOS 
	//Eref1
	if(w>=0.37 && w<=0.43){
			Rr=sqrt(pow(alt,2)+pow(((d2+v+w)/2),2));
			Eref1=RCmul(E0/(2*Rr),Complex(cos(2*pi*2*Rr/LAMBDA),-sin(2*pi*2*Rr/LAMBDA)));	
	}

	//Eref2
	if(w>=0.77 && w<=0.83){
			Rr=sqrt(pow(alt,2)+pow(((d2+v+w)/2),2));
			Eref2=RCmul(E0/(2*Rr),Complex(cos(2*pi*2*Rr/LAMBDA),-sin(2*pi*2*Rr/LAMBDA)));	
	}
	

	//E total
   	Edif_1=Cadd(Edif1,Edif1);
	Edif_2=Cadd(Edif2,Edif2);
	Edif_3=Cadd(Edif3,Edif3);
	Edif_4=Cadd(Edif4,Edif4);

	Edif_5=Cadd(Edif_1,Edif_2);
	Edif_6=Cadd(Edif_3,Edif_4);
	Edif_7=Cadd(Edif_5,Edif_6);
	Edif_8=Cadd(Edir,Edif_7);

	Edif_9=Cadd(Edif_doble1,Edif_doble1);
	Edif_10=Cadd(Edif_doble2,Edif_doble2);
	Edif_11=Cadd(Edif_doble3,Edif_doble3);

	Edif_12=Cadd(Edif_8,Edif_9);
	Edif_13=Cadd(Edif_10,Edif_11);
	Edif_14=Cadd(Edif_12,Edif_13);

	Edif_15=Cadd(Edif_doble4,Edif_doble4);
	Edif_16=Cadd(Edif_doble5,Edif_doble5);
	Edif_17=Cadd(Edif_doble6,Edif_doble6);

	Edif_18=Cadd(Edif_15,Edif_16);
	Edif_19=Cadd(Edif_17,Edif_18);
	Edif_20=Cadd(Edif_14,Edif_19);

	Edif_21=Cadd(Edif_triple1,Edif_triple1);
	Edif_22=Cadd(Edif_triple2,Edif_triple2);
	Edif_23=Cadd(Edif_21,Edif_22);
	Edif_24=Cadd(Edif_20,Edif_23);

	Edif_25=Cadd(Edif_triple3,Edif_triple3);
	Edif_26=Cadd(Edif_triple4,Edif_triple4);
	Edif_27=Cadd(Edif_25,Edif_26);
	//Edif_28=Cadd(Edif_24,Edif_27);
	Edif_28=Cadd(Edif_14,Edif_23);
	Edif_29=Cadd(Edif_cuadruple,Edif_cuadruple);
	//Edif_30=Cadd(Edif_28,Edif_29);
	Edif_30=Cadd(Edif_28,Edif_29);

	Eref_1=Cadd(Eref1,Eref1);
	Eref_2=Cadd(Eref2,Eref2);

	if(w>=0.37 && w<=0.43){
		Et=Cadd(Edif_30,Eref_1);
	}

	else if(w>=0.77 && w<=0.83){
		Et=Cadd(Edif_30,Eref_2);
	
	}
	else {Et=Edif_30;}
	
	fprintf(fsal,"%f %f\n", Rt, 20*log10(Cabs(Et)));
	
	}

//1 obst·culo en medio 
E0=0.002555;
dd=0.42;
alt=0.0375;
v=0.03;

	for (w=0.01;w<1.87;w=w+0.05){
	Rt=dd+v+w;
	R0=sqrt(pow(alt,2)+pow(dd,2));
	alpha=atan(alt/dd);
	l=sqrt(pow(alt,2)+pow(w,2));
	rhoprima=atan(alt/w);
	
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2 - alpha;
	Lw= ((R0/LAMBDA)*(v/LAMBDA))/((R0/LAMBDA)+(v/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((v/LAMBDA)*((R0/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA)));   
	
	//Edif doble total 
	phi=pi + rhoprima;
	phit=45.0*pi/180;
	phiprima=0;           
	Lw=(((R0+v)/LAMBDA)*(l/LAMBDA))/(((R0+v)/LAMBDA)+(l/LAMBDA));          
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v)/LAMBDA)/((l/LAMBDA)*(((R0+v)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));
	
	Et=Cadd(Edif,Edif);
	fprintf(fsal,"%f %f\n", Rt, 20*log10(Cabs(Et)));
	}



//2 obst·culos en medio 
E0=0.002555;
dd=0.42;
d2=0.92;
alt=0.0375;
v=0.03;
d3=d2-dd-v;

	for (w=0.01;w<1.67;w=w+0.05){
	Rt=dd+2*v+d3+w;
	R0=sqrt(pow(alt,2)+pow(dd,2));
	alpha=atan(alt/dd);
	l=sqrt(pow(alt,2)+pow(w,2));
	rhoprima=atan(alt/w);

	//DIFRACCI”N DOBLE
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2-alpha;
	Lw= ((R0/LAMBDA)*(2*(v/LAMBDA)+(d3/LAMBDA)))/((R0/LAMBDA)+2*(v/LAMBDA)+(d3/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((2*(v/LAMBDA)+(d3/LAMBDA))*((R0/LAMBDA)+2*(v/LAMBDA)+(d3/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*(2*(v/LAMBDA)+(d3/LAMBDA))),-sin(2*pi*(2*(v/LAMBDA)+(d3/LAMBDA)))));   
	
	//Edif doble total 
	phi=pi+rhoprima;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+2*v+d3)/LAMBDA)*(l/LAMBDA))/(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+2*v+d3)/LAMBDA)/((l/LAMBDA)*(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_doble1=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	//DIFRACCIONES TRIPLES
	//Edif_triple1
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2 - alpha;
	Lw= ((R0/LAMBDA)*(v/LAMBDA))/((R0/LAMBDA)+(v/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((v/LAMBDA)*((R0/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	  
	//E(2)
	phi=pi;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+v)/LAMBDA)*((d3+v)/LAMBDA))/(((R0+v)/LAMBDA)+((d3+v)/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v)/LAMBDA)/(((d3+v)/LAMBDA)*(((R0+v)/LAMBDA)+((d3+v)/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*(d3+v)/LAMBDA),-sin(2*pi*(d3+v)/LAMBDA))); 
	
	//Dif triple total 
	phi=pi+rhoprima;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+2*v+d3)/LAMBDA)*(l/LAMBDA))/(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R0+2*v+d3)/LAMBDA)/((l/LAMBDA)*(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_triple1=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));


	//Edif_triple2
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2 - alpha;
	Lw= ((R0/LAMBDA)*((v+d3)/LAMBDA))/((R0/LAMBDA)+((v+d3)/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/(((v+d3)/LAMBDA)*((R0/LAMBDA)+((v+d3)/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*(v+d3)/LAMBDA),-sin(2*pi*(v+d3)/LAMBDA))); 
	  
	//E(2)
	phi=pi;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+v+d3)/LAMBDA)*(v/LAMBDA))/(((R0+v+d3)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v+d3)/LAMBDA)/((v/LAMBDA)*(((R0+v+d3)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	
	//Dif triple total 
	phi=pi+rhoprima;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+2*v+d3)/LAMBDA)*(l/LAMBDA))/(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R0+2*v+d3)/LAMBDA)/((l/LAMBDA)*(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_triple2=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));


	//DIFRACCI”N CUADRUPLE 
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2 - alpha;
	Lw= ((R0/LAMBDA)*(v/LAMBDA))/((R0/LAMBDA)+(v/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((v/LAMBDA)*((R0/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	  
	//E(2)
	phi=pi;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+v)/LAMBDA)*(d3/LAMBDA))/(((R0+v)/LAMBDA)+(d3/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v)/LAMBDA)/((d3/LAMBDA)*(((R0+v)/LAMBDA)+(d3/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*d3/LAMBDA),-sin(2*pi*d3/LAMBDA))); 
	
	//E(3)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+d3+v)/LAMBDA)*(v/LAMBDA))/(((R0+d3+v)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R0+d3+v)/LAMBDA)/((v/LAMBDA)*(((R0+d3+v)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E3=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA)));
	
	//Dif cuadruple total 
	phi=pi+rhoprima;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+2*v+d3)/LAMBDA)*(l/LAMBDA))/(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E3,CoefD);
	A1=sqrt(((R0+2*v+d3)/LAMBDA)/((l/LAMBDA)*(((R0+2*v+d3)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_cuadruple=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	Edif_1=Cadd(Edif_doble1,Edif_doble1);
	Edif_2=Cadd(Edif_triple1,Edif_triple1);
	Edif_3=Cadd(Edif_triple2,Edif_triple2);
	Edif_4=Cadd(Edif_cuadruple,Edif_cuadruple);

	Edif_5=Cadd(Edif_1,Edif_2);
	Edif_6=Cadd(Edif_3,Edif_4);
	Et=Cadd(Edif_5,Edif_6);
	


	fprintf(fsal,"%f %f\n", Rt, 20*log10(Cabs(Et)));
	
	}

*/

//obst·culos combinados 
E0=0.002555;
dd=0.6;
d3=0.8;
v=0.03;
d2=d3-dd-v;
alt=0.032;
alt2=0.0375;

	for (w=0.01;w<1.67;w=w+0.05){
	Rt=d3+v+w;
	R22=sqrt(pow(d3,2)+pow(alt2,2));
	R0=sqrt(pow(alt,2)+pow(dd,2));
	R1=sqrt(pow(alt,2)+pow((dd+v),2));
	alpha=atan(alt2/d3);
	betha=atan(alt/(dd+v));
	bethaprima=atan(alt/dd);
	l=sqrt(pow(w,2)+pow(alt2,2));
	rho=atan(alt2/w);

	//DIFRACCI”N DOBLE
	//Edif_doble1
	//E(0)
	Eini=RCmul(E0/R22,Complex(cos(2*pi*R22/LAMBDA),-sin(2*pi*R22/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2-alpha;
	Lw= ((R22/LAMBDA)*(v/LAMBDA))/((R22/LAMBDA)+(v/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R22/LAMBDA)/((v/LAMBDA)*((R22/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*(v/LAMBDA)),-sin(2*pi*(v/LAMBDA))));   
	
	//Edif doble total 
	phi=pi+rho;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R22+v)/LAMBDA)*(l/LAMBDA))/(((R22+v)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R22+v)/LAMBDA)/((l/LAMBDA)*(((R22+v)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_doble1=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));


	//DIFRACCIONES TRIPLES
	//Edif_triple1
	R3=sqrt(pow((alt2-alt),2)+pow(d2,2));
	tau=atan((alt2-alt)/d2);
	//E(0)
	Eini=RCmul(E0/R1,Complex(cos(2*pi*R1/LAMBDA),-sin(2*pi*R1/LAMBDA)));

	//E(1)
	phi=pi+tau;
	phit=45.0*pi/180;
	phiprima=betha;
	Lw= ((R1/LAMBDA)*(R3/LAMBDA))/((R1/LAMBDA)+(R3/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R1/LAMBDA)/((R3/LAMBDA)*((R1/LAMBDA)+(R3/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*R3/LAMBDA),-sin(2*pi*R3/LAMBDA))); 
	  
	//E(2)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2-tau;
	Lw=(((R1+R3)/LAMBDA)*(v/LAMBDA))/(((R1+R3)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R1+R3)/LAMBDA)/((v/LAMBDA)*(((R1+R3)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	
	//Dif triple total 
	phi=pi+rho;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R1+R3+v)/LAMBDA)*(l/LAMBDA))/(((R1+R3+v)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R1+R3+v)/LAMBDA)/((l/LAMBDA)*(((R1+R3+v)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_triple1=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	//Edif_triple2
	R4=sqrt(pow((alt+alt2),2)+pow((d2+v),2));
	omega=atan((alt+alt2)/(d2+v));
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2-omega;
	phit=45.0*pi/180;
	phiprima=pi/2+bethaprima;
	Lw= ((R0/LAMBDA)*(R4/LAMBDA))/((R0/LAMBDA)+(R4/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((R4/LAMBDA)*((R0/LAMBDA)+(R4/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*R4/LAMBDA),-sin(2*pi*R4/LAMBDA))); 
	  
	//E(2)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2-omega;
	Lw=(((R0+R4)/LAMBDA)*(v/LAMBDA))/(((R0+R4)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+R4)/LAMBDA)/((v/LAMBDA)*(((R0+R4)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	
	//Dif triple total 
	phi=pi+rho;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+R4+v)/LAMBDA)*(l/LAMBDA))/(((R0+R4+v)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R0+R4+v)/LAMBDA)/((l/LAMBDA)*(((R0+R4+v)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_triple2=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	//Edif_triple3
	R5=sqrt(pow((alt+alt2),2)+pow(d2,2));
	theta=atan((alt+alt2)/d2);
	//E(0)
	Eini=RCmul(E0/R1,Complex(cos(2*pi*R1/LAMBDA),-sin(2*pi*R1/LAMBDA)));

	//E(1)
	phi=pi-theta;
	phit=45.0*pi/180;
	phiprima=betha;
	Lw= ((R1/LAMBDA)*(R5/LAMBDA))/((R1/LAMBDA)+(R5/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R1/LAMBDA)/((R5/LAMBDA)*((R1/LAMBDA)+(R5/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*R5/LAMBDA),-sin(2*pi*R5/LAMBDA))); 
	  
	//E(2)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2-theta;
	Lw=(((R1+R5)/LAMBDA)*(v/LAMBDA))/(((R1+R5)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R1+R5)/LAMBDA)/((v/LAMBDA)*(((R1+R5)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	
	//Dif triple total 
	phi=pi+rho;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R1+R5+v)/LAMBDA)*(l/LAMBDA))/(((R1+R5+v)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R1+R5+v)/LAMBDA)/((l/LAMBDA)*(((R1+R5+v)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_triple3=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));


	//DIFRACCIONES CUADRUPLES
	//Edif_cuadruple1
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2+bethaprima;
	Lw= ((R0/LAMBDA)*(v/LAMBDA))/((R0/LAMBDA)+(v/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((v/LAMBDA)*((R0/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	  
	//E(2)
	phi=pi+tau;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+v)/LAMBDA)*(R3/LAMBDA))/(((R0+v)/LAMBDA)+(R3/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v)/LAMBDA)/((R3/LAMBDA)*(((R0+v)/LAMBDA)+(R3/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*R3/LAMBDA),-sin(2*pi*R3/LAMBDA))); 
	
	//E(3)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2-tau;
	Lw=(((R0+R3+v)/LAMBDA)*(v/LAMBDA))/(((R0+R3+v)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R0+R3+v)/LAMBDA)/((v/LAMBDA)*(((R0+R3+v)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E3=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA)));
	
	//Dif cuadruple total 
	phi=pi+rho;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+2*v+R3)/LAMBDA)*(l/LAMBDA))/(((R0+2*v+R3)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E3,CoefD);
	A1=sqrt(((R0+2*v+R3)/LAMBDA)/((l/LAMBDA)*(((R0+2*v+R3)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_cuadruple1=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	//Edif_cuadruple2
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));
	//fprintf(fsal,"%f %f\n", Rt,(Cabs(Eini)));

	//E(1)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2+bethaprima;
	Lw= ((R0/LAMBDA)*(v/LAMBDA))/((R0/LAMBDA)+(v/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((v/LAMBDA)*((R0/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA))); 
	  
	//E(2)
	phi=pi-theta;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+v)/LAMBDA)*(R5/LAMBDA))/(((R0+v)/LAMBDA)+(R5/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+v)/LAMBDA)/((R5/LAMBDA)*(((R0+v)/LAMBDA)+(R5/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*R5/LAMBDA),-sin(2*pi*R5/LAMBDA))); 
	
	//E(3)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2-theta;
	Lw=(((R0+R5+v)/LAMBDA)*(v/LAMBDA))/(((R0+R5+v)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R0+R5+v)/LAMBDA)/((v/LAMBDA)*(((R0+R5+v)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E3=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA)));
	
	//Dif cuadruple total 
	phi=pi+rho;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+2*v+R5)/LAMBDA)*(l/LAMBDA))/(((R0+2*v+R5)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E3,CoefD);
	A1=sqrt(((R0+2*v+R5)/LAMBDA)/((l/LAMBDA)*(((R0+2*v+R5)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_cuadruple2=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	//Edif_cuadruple3
	R6=2*alt;
	//E(0)
	Eini=RCmul(E0/R1,Complex(cos(2*pi*R1/LAMBDA),-sin(2*pi*R1/LAMBDA)));

	//E(1)
	phi=pi/2;
	phit=45.0*pi/180;
	phiprima=betha;
	Lw= ((R1/LAMBDA)*(R6/LAMBDA))/((R1/LAMBDA)+(R6/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R1/LAMBDA)/((R6/LAMBDA)*((R1/LAMBDA)+(R6/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*R6/LAMBDA),-sin(2*pi*R6/LAMBDA))); 
	  
	//E(2)
	phi=pi+tau;
	phit=45.0*pi/180;
	phiprima=pi/2;
	Lw=(((R1+R6)/LAMBDA)*(R3/LAMBDA))/(((R1+R6)/LAMBDA)+(R3/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R1+R6)/LAMBDA)/((R3/LAMBDA)*(((R1+R6)/LAMBDA)+(R3/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*R3/LAMBDA),-sin(2*pi*R3/LAMBDA))); 
	
	//E(3)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2+tau;
	Lw=(((R1+R6+R3)/LAMBDA)*(v/LAMBDA))/(((R1+R6+R3)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R1+R6+R3)/LAMBDA)/((v/LAMBDA)*(((R1+R6+R3)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E3=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA)));
	
	//Dif cuadruple total 
	phi=pi+rho;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R1+R6+R3+v)/LAMBDA)*(l/LAMBDA))/(((R1+R6+R3+v)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E3,CoefD);
	A1=sqrt(((R1+R6+R3+v)/LAMBDA)/((l/LAMBDA)*(((R1+R6+R3+v)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_cuadruple3=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	//Edif_cuadruple4
	//E(0)
	Eini=RCmul(E0/R1,Complex(cos(2*pi*R1/LAMBDA),-sin(2*pi*R1/LAMBDA)));

	//E(1)
	phi=pi/2;
	phit=45.0*pi/180;
	phiprima=betha;
	Lw= ((R1/LAMBDA)*(R6/LAMBDA))/((R1/LAMBDA)+(R6/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R1/LAMBDA)/((R6/LAMBDA)*((R1/LAMBDA)+(R6/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*R6/LAMBDA),-sin(2*pi*R6/LAMBDA))); 
	  
	//E(2)
	phi=pi/2+theta;
	phit=45.0*pi/180;
	phiprima=pi/2;
	Lw=(((R1+R6)/LAMBDA)*(R5/LAMBDA))/(((R1+R6)/LAMBDA)+(R5/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R1+R6)/LAMBDA)/((R5/LAMBDA)*(((R1+R6)/LAMBDA)+(R5/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*R5/LAMBDA),-sin(2*pi*R5/LAMBDA))); 
	
	//E(3)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2-theta;
	Lw=(((R1+R6+R5)/LAMBDA)*(v/LAMBDA))/(((R1+R6+R5)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R1+R6+R5)/LAMBDA)/((v/LAMBDA)*(((R1+R6+R5)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E3=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA)));
	
	//Dif cuadruple total 
	phi=pi+rho;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R1+R6+R5+v)/LAMBDA)*(l/LAMBDA))/(((R1+R6+R5+v)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E3,CoefD);
	A1=sqrt(((R1+R6+R5+v)/LAMBDA)/((l/LAMBDA)*(((R1+R6+R5+v)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_cuadruple4=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	//Edif_cuadruple5
	R7=sqrt(pow((2*alt),2)+pow(v,2));
	mu=atan((2*alt)/v);
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2-mu;
	phit=45.0*pi/180;
	phiprima=pi/2+bethaprima;
	Lw= ((R0/LAMBDA)*(R7/LAMBDA))/((R0/LAMBDA)+(R7/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((R7/LAMBDA)*((R0/LAMBDA)+(R7/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*R7/LAMBDA),-sin(2*pi*R7/LAMBDA))); 
	  
	//E(2)
	phi=3*pi/2+tau;
	phit=45.0*pi/180;
	phiprima=pi/2-mu;
	Lw=(((R0+R7)/LAMBDA)*(R3/LAMBDA))/(((R0+R7)/LAMBDA)+(R3/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+R7)/LAMBDA)/((R3/LAMBDA)*(((R0+R7)/LAMBDA)+(R3/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*R3/LAMBDA),-sin(2*pi*R3/LAMBDA))); 
	
	//E(3)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2+tau;
	Lw=(((R0+R7+R3)/LAMBDA)*(v/LAMBDA))/(((R0+R7+R3)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R0+R7+R3)/LAMBDA)/((v/LAMBDA)*(((R0+R7+R3)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E3=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA)));
	
	//Dif cuadruple total 
	phi=pi+rho;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+R7+R3+v)/LAMBDA)*(l/LAMBDA))/(((R0+R7+R3+v)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E3,CoefD);
	A1=sqrt(((R0+R7+R3+v)/LAMBDA)/((l/LAMBDA)*(((R0+R7+R3+v)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_cuadruple5=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));
	
	//Edif_cuadruple6
	//E(0)
	Eini=RCmul(E0/R0,Complex(cos(2*pi*R0/LAMBDA),-sin(2*pi*R0/LAMBDA)));

	//E(1)
	phi=3*pi/2-mu;
	phit=45.0*pi/180;
	phiprima=pi/2+bethaprima;
	Lw= ((R0/LAMBDA)*(R7/LAMBDA))/((R0/LAMBDA)+(R7/LAMBDA));

	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(Eini,CoefD);
	A1=sqrt((R0/LAMBDA)/((R7/LAMBDA)*((R0/LAMBDA)+(R7/LAMBDA))));
	var2=RCmul(A1,var1);
	E1=Cmul(var2,Complex(cos(2*pi*R7/LAMBDA),-sin(2*pi*R7/LAMBDA))); 
	  
	//E(2)
	phi=pi-theta;
	phit=45.0*pi/180;
	phiprima=pi/2-mu;
	Lw=(((R0+R7)/LAMBDA)*(R5/LAMBDA))/(((R0+R7)/LAMBDA)+(R5/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E1,CoefD);
	A1=sqrt(((R0+R7)/LAMBDA)/((R5/LAMBDA)*(((R0+R7)/LAMBDA)+(R5/LAMBDA))));
	var2=RCmul(A1,var1);
	E2=Cmul(var2,Complex(cos(2*pi*R5/LAMBDA),-sin(2*pi*R5/LAMBDA))); 
	
	//E(3)
	phi=3*pi/2;
	phit=45.0*pi/180;
	phiprima=pi/2-theta;
	Lw=(((R0+R7+R5)/LAMBDA)*(v/LAMBDA))/(((R0+R7+R5)/LAMBDA)+(v/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E2,CoefD);
	A1=sqrt(((R0+R7+R5)/LAMBDA)/((v/LAMBDA)*(((R0+R7+R5)/LAMBDA)+(v/LAMBDA))));
	var2=RCmul(A1,var1);
	E3=Cmul(var2,Complex(cos(2*pi*v/LAMBDA),-sin(2*pi*v/LAMBDA)));
	
	//Dif cuadruple total 
	phi=pi+rho;
	phit=45.0*pi/180;
	phiprima=0;
	Lw=(((R0+R7+R5+v)/LAMBDA)*(l/LAMBDA))/(((R0+R7+R5+v)/LAMBDA)+(l/LAMBDA));
	
	CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
	var1=Cmul(E3,CoefD);
	A1=sqrt(((R0+R7+R5+v)/LAMBDA)/((l/LAMBDA)*(((R0+R7+R5+v)/LAMBDA)+(l/LAMBDA))));
	var2=RCmul(A1,var1);
	Edif_cuadruple6=Cmul(var2,Complex(cos(2*pi*l/LAMBDA),-sin(2*pi*l/LAMBDA)));

	Edif_1=Cadd(Edif_doble1,Edif_doble1);

	Edif_2=Cadd(Edif_triple1,Edif_triple1);
	Edif_3=Cadd(Edif_1,Edif_2);
	Edif_4=Cadd(Edif_triple2,Edif_triple2);
	Edif_5=Cadd(Edif_3,Edif_4);
	Edif_6=Cadd(Edif_triple3,Edif_triple3);
	Edif_7=Cadd(Edif_5,Edif_6);

	Edif_8=Cadd(Edif_cuadruple1,Edif_cuadruple1);
	Edif_9=Cadd(Edif_7,Edif_8);
	Edif_10=Cadd(Edif_cuadruple2,Edif_cuadruple2);
	Edif_11=Cadd(Edif_9,Edif_10);
	Edif_12=Cadd(Edif_cuadruple3,Edif_cuadruple3);
	Edif_13=Cadd(Edif_11,Edif_12);
	Edif_14=Cadd(Edif_cuadruple4,Edif_cuadruple4);
	Edif_15=Cadd(Edif_13,Edif_14);
	Edif_16=Cadd(Edif_cuadruple5,Edif_cuadruple5);
	Edif_17=Cadd(Edif_15,Edif_16);
	Edif_18=Cadd(Edif_cuadruple6,Edif_cuadruple6);

	Et=Cadd(Edif_17,Edif_18);
	
	fprintf(fsal,"%f %f\n", Rt, 20*log10(Cabs(Et)));
}


 for (y0=1.0;y0<2.0;y0=y0+2.0){
/////////////////////////////////for (dd=0.10;dd<0.90;dd=dd+0.01){
///////////////////////////////// for (alfa=11.309;alfa<12.0;alfa=alfa+1.0){

/////////////////////////////////     alpha=alfa*pi/180;
     alpha=atan(y0/dd);
/////////////////////////////////     y0=dd*tan(alpha);
     for (i=1;i<N;i++)
     AnUTD[i]=ZERO;
     for (i=1;i<N;i++)
     AnUTDT[i]=ZERO;
     for (i=0;i<=2*N;i++){
       if ((i%2)==0){
           R[i]=sqrt(pow(y0,2)+pow(dd+(i/2)*(v*LAMBDA+w*LAMBDA),2));
           }
      else{
           R[i]=sqrt(pow(y0,2)+pow(dd+((i-1)/2)*(v*LAMBDA+w*LAMBDA)+v*LAMBDA,2));
           R2[i]=sqrt(pow(y0,2)+pow(dd+((i-1)/2)*(v*LAMBDA+w*LAMBDA)+w*LAMBDA,2));
           }
     		 }
 //      for (j=0;j<=2*N;j++){
 //     if ((j%2)==0){
//           R2[j]=sqrt(pow(y0,2)+pow(dd+((j-1)/2)*(v*LAMBDA+w*LAMBDA)+w*LAMBDA,2));
    //       R2[j]=0.0;
 //          }
 //     else{
   //        RR[i]=sqrt(pow(y0,2)+pow(dd+((i-1)/2)*(v*LAMBDA+w*LAMBDA)+w*LAMBDA,2));
  //         BB[j]=0;
  //         }
  //   		 }
 //    for (i=1;i<N;i++)
 //    Em[i]=ZERO;

     for (n=1;n<N;n++) {

 // Complex(cos(2*pi*v),-sin(2*pi*v))
         for (m=0;m<=n-1;m++) {
         
         phiprima=0.5*pi+alpha;
 //        phiprima=0.5*pi+alpha-0.012*phit;
         for (i=m;i<N;i++)
         Em[i]=ZERO;
               if (m==0) {
               Em[m]=RCmul(Ei/R[0],Complex(cos(2*pi*R[0]/LAMBDA),-sin(2*pi*R[0]/LAMBDA)));} else {

               for (q=0;q<=m-1;q++){

                  Lw3=(R[0]/LAMBDA*(m-q)*(v+w))/(R[0]/LAMBDA+(m-q)*(v+w));
                  CoefD3=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw3);
                  E=RCmul(R[0]/R[2*(m-q)],Complex(cos(2*pi*(R[2*m]/LAMBDA-R[2*q]/LAMBDA)),-sin(2*pi*(R[2*m]/LAMBDA-R[2*q]/LAMBDA))));
                  F=Cmul(RCmul(sqrt((R[0]/LAMBDA)/((m-q)*(v+w)*(R[0]/LAMBDA+(m-q)*(v+w)))),CoefD3),Complex(cos(2*pi*(m-q)*(v+w)),-sin(2*pi*(m-q)*(v+w))));
                      if (alpha<0.0) {
                  Em[m]=Cadd(Em[m],Cmul(Em[q],Cadd(ZERO,F)));
                                      } else {
                  Em[m]=Cadd(Em[m],Cmul(Em[q],Cadd(E,F)));
                                              }
                                    } // del for de q

               for (r=1;r<=m;r++){

               //   if (alpha<0.0) {phiprima_2=0.01;}  // Para angulos negativos te metes en los edificios!!!!
               //   else {phiprima_2=alpha;}
                  phiprima_2=alpha;
               //   phiprima_2=0.5*pi+alpha-0.988*phit;
                  Lw4=(R[0]/LAMBDA*((m-r)*(v+w)+w))/(R[0]/LAMBDA+((m-r)*(v+w)+w));
                  CoefD4=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi_2,phiprima_2,Lw4);
                  G=RCmul(R[0]/R2[2*(m-r)+1],Complex(cos(2*pi*(R[2*m]/LAMBDA-R[2*r-1]/LAMBDA)),-sin(2*pi*(R[2*m]/LAMBDA-R[2*r-1]/LAMBDA))));
                  H=Cmul(RCmul(sqrt((R[0]/LAMBDA)/(((m-r)*(v+w)+w)*(R[0]/LAMBDA+(m-r)*(v+w)+w))),CoefD4),Complex(cos(2*pi*((m-r)*(v+w)+w)),-sin(2*pi*((m-r)*(v+w)+w))));
                      if (alpha<0.0) {
                  Em[m]=Cadd(Em[m],Cmul(AnUTDT[r],Cadd(ZERO,H)));
                                     } else {
                  Em[m]=Cadd(Em[m],Cmul(AnUTDT[r],Cadd(G,H)));
                                             }
                                    } // del for de r

         		Em[m]=RCmul(1.0/(2*m),Em[m]);
/////////////////////////////////////////////////////////////               if ((n==2) && (m==n-1)) {         //Hay que poner 7 para 6 edificios completos y Rx a w del ˙ltimo edificio
///////////////////            if (m==n-1)  {
/////////////////////////////////////////////////////////////            fprintf(fsal,"%d %f %f %f\n",n,dd,R[2*n-1],-20*log10((Cabs(Em[m])*R[2*n-1])));}
////////////////////////////////            fprintf(fsal,"%d %f %f %f\n",n,dd,R[2*n-1],-20*log10((Cabs(Em[m]))));}
/////	      {fprintf(fsal,"%d %f %f %f\n",n,y0,Em[m].r,Em[m].i);}
//               if (m==(n-1)) {
//               fprintf(fsal,"%d %f %f %f\n",n,y0,Em[m].r,Em[m].i);}
                      } // del else del if (m==0)

         Lw=(R[0]/LAMBDA*((n-m)*(v+w)-w))/(R[0]/LAMBDA+((n-m)*(v+w)-w));
         CoefD=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
         A=RCmul(R[0]/R[2*(n-m)-1],Complex(cos(2*pi*(R[2*n-1]/LAMBDA-R[2*m]/LAMBDA)),-sin(2*pi*(R[2*n-1]/LAMBDA-R[2*m]/LAMBDA))));
         B=Cmul(RCmul(sqrt((R[0]/LAMBDA)/(((n-m)*(v+w)-w)*(R[0]/LAMBDA+(n-m)*(v+w)-w))),CoefD),Complex(cos(2*pi*((n-m)*(v+w)-w)),-sin(2*pi*((n-m)*(v+w)-w))));
              if (alpha<0.0){
         AnUTD[n]=Cadd(AnUTD[n],Cmul(Em[m],Cadd(ZERO,B)));
                           }  else {
         AnUTD[n]=Cadd(AnUTD[n],Cmul(Em[m],Cadd(A,B)));
                           }
                             }     // del for de m
         if (n!=1){
         for (p=1;p<=n-1;p++){

         phiprima_2=alpha;
     //    phiprima_2=0.5*pi+alpha-0.988*phit;
         Lw2=(R[0]/LAMBDA*(n-p)*(v+w))/(R[0]/LAMBDA+(n-p)*(v+w));
         CoefD2=Coef_difrac(SOFT,1,Complex(6.3,-0.35),2.0*phit,2.0*pi,pi/2.0,phi_2,phiprima_2,Lw2);
         C=RCmul(R[0]/R[2*(n-p)],Complex(cos(2*pi*(R[2*n-1]/LAMBDA-R[2*p-1]/LAMBDA)),-sin(2*pi*(R[2*n-1]/LAMBDA-R[2*p-1]/LAMBDA))));
         D=Cmul(RCmul(sqrt((R[0]/LAMBDA)/((n-p)*(v+w)*(R[0]/LAMBDA+(n-p)*(v+w)))),CoefD2),Complex(cos(2*pi*((n-p)*(v+w))),-sin(2*pi*((n-p)*(v+w)))));
               if (alpha<0.0) {
         AnUTD[n]=Cadd(AnUTD[n],Cmul(AnUTDT[p],Cadd(ZERO,D)));
                              } else {
         AnUTD[n]=Cadd(AnUTD[n],Cmul(AnUTDT[p],Cadd(C,D)));
                              }
                            }      // del for de p
                   }  //del if (n!=1)

      AnUTD[n]=RCmul(1.0/(2*n-1),AnUTD[n]);

     AnUTDT[n]=AnUTD[n];

//     if (n==N-1) {
////////////////////////////////////////////////////////////////////////     fprintf(fsal,"%d %f %f %f\n",n,y0,alpha*180/pi,-20*log10(Cabs(AnUTD[n])*R[2*n-1]));//}
//	  fprintf(fsal,"%d %f %f %f %f\n",n,y0,w,v,Cabs(AnUTD[n]));
//     fprintf(fsal,"%d %f %f\n",n,y0,Cabs(Cmul(AnUTD[n],Complex(cos(2*pi*(N-1)*(w+v)*cos(alpha)),-sin(2*pi*(N-1)*(w+v)*cos(alpha))))));}
     }    //del for de n
  //   } // del for de v
  //   } // del for de w

//     fprintf(fsal,"%d %f %f\n",n,alpha,Cabs(AnUTD[n-1])*R[2*(n-1)-1]);
//     phiprima=atan(y0/175);
//     phi=1.5*pi-83.157*pi/180;
//     Lw=(sqrt(y0*y0+175*175)/LAMBDA*(100.717/LAMBDA))/(sqrt(y0*y0+175*175)/LAMBDA+(100.717/LAMBDA));
//     Lw=100.717/LAMBDA;
//     CoefD=Coef_difrac(SOFT,1,Complex(5.5,-0.46),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
//     AA=sqrt((sqrt(y0*y0+175*175)/LAMBDA)/((100.717/LAMBDA)*((100.717/LAMBDA)+(sqrt(y0*y0+175*175)/LAMBDA))));
//     AA=1.0/sqrt(100.717/LAMBDA);
//     AnUTDaux[n]=Cmul(AnUTDT[n-1],RCmul(AA,Cmul(CoefD,Complex(cos(2*pi*(100.717/LAMBDA)),-sin(2*pi*(100.717/LAMBDA))))));
//     if  (y0>21.1) {
//           RR=sqrt((y0+12)*(y0+12)+(275)*(275));
//           AnUTDaux[n]=Cadd(AnUTDaux[n],RCmul(Cabs(AnUTDT[n-1])*(sqrt(y0*y0+175*175))/RR,Complex(cos(2*pi*(RR/LAMBDA)),-sin(2*pi*(RR/LAMBDA)))));
  //         fprintf(fsal,"%d %f %f %f\n",n,phiprima,y0,20*log10(Cabs(AnUTDaux[n])*sqrt((y0+12)*(y0+12)+275*275)));
//     					 }
  //     fprintf(fsal,"%d %f %f %f\n",n,alpha*180/pi,y0,20*log10(Cabs(AnUTDaux[n])*sqrt((y0+12)*(y0+12)+275*275)));
  //   fprintf(fsal,"%d %f %f %f\n",n,phiprima,y0,20*log10(Cabs(AnUTDaux[n])*sqrt((y0+12)*(y0+12)+275*275)));
//   fprintf(fsal,"%d %f %f %f\n",n,phiprima,y0,20*log10(Cabs(AnUTDaux[n])*sqrt((y0+12)*(y0+12)+275*275)/(Cabs(AnUTDT[n-1])*(sqrt(y0*y0+175*175)))));

      } //del for de y0 Û de alpha

//     fprintf(fsal,"%d %f %f\n",n,alfa,Cabs(AnUTD[n-1])*R[2*(n-1)-1]);
//     phiprima=1.636*pi/180;
//     phi=1.5*pi-83.157*pi/180;                                                 // Htx=5ù
//     Lw=(175.071/LAMBDA*100.717/LAMBDA)/(175.071/LAMBDA+100.717/LAMBDA);
//     CoefD=Coef_difrac(SOFT,1,Complex(5.5,-0.46),2.0*phit,2.0*pi,pi/2.0,phi,phiprima,Lw);
//     AA=sqrt(175.071/(100.717*(100.717+175.071)));
//     AnUTD[n]=Cmul(AnUTD[n-1],RCmul(AA,Cmul(CoefD,Complex(cos(2*pi*100.717),-sin(2*pi*100.717)))));
//     fprintf(fsal,"%d %f %f %f\n",n,alfa,Cabs(AnUTD[n]),Cabs(AnUTD[n])*275.524);


    fclose(fsal);


	cout<<("\n\n** FINAL DEL PROGRAMA **\n");
	int c;
   cin>>c;

   exit(1);
}


