#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include "marumi.h" 

#define Pi           3.141592
#define PARA_NUM     5
#define Rearth 6378.4
#define ERROR_UNDEF  -999.9
//#define Vf -1.0 //Fall speed of snow [m/s]
#define SDEV_THRESHOLD 1.7
#define SDEV_RANGE  4
#define RHO0 1.2754 //#kg/m^3  

int main(int argc, char* argv[]){
  
  int k,j,i,ncid,varid,bins;
  size_t tot_numray2,max_numbin2,sweep2;
  int tot_numray,max_numbin,tot_numbin,sweep,numray;

  FILE *wfp;

  char infile[200],parav[100],paras[100];
  float *vels,*sdvels, *vels1,*sdvels1,*refs,*refs1;
  float *temps,*rhods;
  double *temps1,*rhods1;
  float *az,*el,*range,az_offset;
  int *sweep_start,*sweep_end;
  float *fixed_angle,*time2;
  float *time;
  double r_alt;

  float SDVEL_THRE;

  int n,m,l,pivot;
  double max;
  double Smm[PARA_NUM][PARA_NUM];
  double Pnm[PARA_NUM];
  double Gmm1[PARA_NUM][PARA_NUM+1];
  double Km[PARA_NUM];
  double tmp1[1][PARA_NUM+1];
  float *K1,*K2,*K3,*K4,*K5,*U,*V,*WSPD,*WDRC,*HGT,*DIV,*DEF,*RADIUS;
  int data_num;

  double Az,El;
  double vel;
  float sdvel;
  double r;
  double x,y,z;
  float undefV,undefV2;
  float sumz,avez,sumel,aveel,sumxy,avexy,xy;

  double u0,v0;
  float ws,wd,u0d,v0d,divg,defm;
  double G00,G11;
  int flag;

  int nc_error;
  int range_dimid,sweep_dimid,mheight_dimid;
  int dimid2s[2],dimidrange[1],dimidsweep[1],dimidmheight[1];
  int varidk1,varidk2,varidk3,varidk4,varidk5,varidws,varidwd;
  int varidrange,varidfixed_angle,varidhgt,variddiv,varidt,varidrad;
  int variddef,varidu,varidv,varidvf;
  float undef;
  float sign;
  float Vf,Vfi,rhod;

  int a,aa,aaa,counta,count,num_h;
  float h,hc;
  //  float hmax = 10, dh = 0.05; //km
  float hmax = 13, dh = 0.05; //km // changed 2020/12/26
  float sumu,sumv,sums,sumd,sumdiv,sumdef;
  float aveu,avev,aves,aved,avediv,avedef;
  float sdevu,sdevv,sdevs,sdevd,sdevdiv,sdevdef;
  float *Ufil,*Vfil,*WSPDfil,*WDRCfil,*DIVfil,*DEFfil;
  float *Ufil2,*Vfil2,*WSPDfil2,*WDRCfil2,*DIVfil2,*DEFfil2;
  float *AU,*AV,*AWSPD,*AWDRC,*ADIV,*ADEF,*AHGT;
  float *SU,*SV,*SWSPD,*SWDRC,*SDIV,*SDEF;
  int varidu2,varidv2,varidws2,varidwd2,variddiv2,variddef2;
  int varidu3,varidv3,varidws3,varidwd3,variddiv3,variddef3,varidhgt3;
  int varidu4,varidv4,varidws4,varidwd4,variddiv4,variddef4;
  num_h = (int)(hmax / dh + 0.1);

  struct MARUMI Maru; 
  
  if(argc >= 4 ){
    strcpy(infile,argv[1]);
    printf("infile = %s \n" ,infile);
    sign = atof(argv[2]); 
    //    Vf = atof(argv[3]);
    //    printf("Fall Speed: %f \n",Vf);
    strcpy(parav,argv[3]);
    strcpy(paras,argv[4]);
    SDVEL_THRE = atof(argv[5]); 
    printf("SDVEL_THRE = %f \n",SDVEL_THRE);
    az_offset =  atof(argv[6]);
    printf("Az_OFFSETE = %f \n",az_offset);
    /*
    Z_MAX = atof(argv[3]); 
    DZ = atof(argv[4]); 
    printf("Z_MAX = %f (km)\n",Z_MAX);
    printf("DZ = %f (km)\n",DZ);*/
    }else{
    printf("error ARGV \n");
    exit(1);
  }

  //  Z_NUM = (int)(Z_MAX/DZ);

  /* bin_num ray_num */
  k = nc_open(infile,0,&ncid);
  k = nc_inq_dimid(ncid,"time",&varid);  
  k = nc_inq_dimlen(ncid,varid,&tot_numray2);
  k = nc_inq_dimid(ncid,"range",&varid);  
  k = nc_inq_dimlen(ncid,varid,&max_numbin2);
  k = nc_inq_dimid(ncid,"sweep",&varid);  
  k = nc_inq_dimlen(ncid,varid,&sweep2);
  tot_numray = tot_numray2;
  max_numbin = max_numbin2;
  tot_numbin = tot_numray*max_numbin;
  sweep = sweep2;
  printf("num_ray = %d; num_bin = %d(%d); num_sweep = %d \n",tot_numray,max_numbin,tot_numbin,sweep);

  vels   = (void *)malloc(sizeof(float)*tot_numray*max_numbin);
  sdvels = (void *)malloc(sizeof(float)*tot_numray*max_numbin);
  refs = (void *)malloc(sizeof(float)*tot_numray*max_numbin);
  temps = (void *)malloc(sizeof(float)*tot_numray*max_numbin);
  rhods = (void *)malloc(sizeof(float)*tot_numray*max_numbin);
  vels1   = (void *)malloc(sizeof(float)*tot_numray*max_numbin);
  sdvels1 = (void *)malloc(sizeof(float)*tot_numray*max_numbin);
  refs1 = (void *)malloc(sizeof(float)*tot_numray*max_numbin);
  temps1 = (void *)malloc(sizeof(double)*tot_numray*max_numbin);
  rhods1 = (void *)malloc(sizeof(double)*tot_numray*max_numbin);
  az    = (float*)malloc(sizeof(float)*tot_numray);
  el    = (float*)malloc(sizeof(float)*tot_numray); 
  time  = (float*)malloc(sizeof(float)*tot_numray); 
  range = (float*)malloc(sizeof(float)*max_numbin); 
  sweep_start = (int*)malloc(sizeof(int)*sweep); 
  sweep_end   = (int*)malloc(sizeof(int)*sweep); 
  fixed_angle = (void *)malloc(sizeof(float)*sweep);

  K1  = (void *)malloc(sizeof(float)*sweep*max_numbin);
  K2  = (void *)malloc(sizeof(float)*sweep*max_numbin);
  K3  = (void *)malloc(sizeof(float)*sweep*max_numbin);
  K4  = (void *)malloc(sizeof(float)*sweep*max_numbin);
  K5  = (void *)malloc(sizeof(float)*sweep*max_numbin);
  U   = (void *)malloc(sizeof(float)*sweep*max_numbin);
  V   = (void *)malloc(sizeof(float)*sweep*max_numbin);
  WSPD= (void *)malloc(sizeof(float)*sweep*max_numbin);
  WDRC= (void *)malloc(sizeof(float)*sweep*max_numbin);
  HGT = (void *)malloc(sizeof(float)*sweep*max_numbin);
  DIV = (void *)malloc(sizeof(float)*sweep*max_numbin);
  DEF = (void *)malloc(sizeof(float)*sweep*max_numbin);
  time2 = (void *)malloc(sizeof(float)*sweep);
  RADIUS = (void *)malloc(sizeof(float)*sweep*max_numbin);
  Ufil   = (void *)malloc(sizeof(float)*sweep*max_numbin);
  Vfil   = (void *)malloc(sizeof(float)*sweep*max_numbin);
  WSPDfil= (void *)malloc(sizeof(float)*sweep*max_numbin);
  WDRCfil= (void *)malloc(sizeof(float)*sweep*max_numbin);
  DIVfil = (void *)malloc(sizeof(float)*sweep*max_numbin);
  DEFfil = (void *)malloc(sizeof(float)*sweep*max_numbin);
  Ufil2   = (void *)malloc(sizeof(float)*sweep*max_numbin);
  Vfil2   = (void *)malloc(sizeof(float)*sweep*max_numbin);
  WSPDfil2= (void *)malloc(sizeof(float)*sweep*max_numbin);
  WDRCfil2= (void *)malloc(sizeof(float)*sweep*max_numbin);
  DIVfil2 = (void *)malloc(sizeof(float)*sweep*max_numbin);
  DEFfil2 = (void *)malloc(sizeof(float)*sweep*max_numbin);
  AU   = (void *)malloc(sizeof(float)*num_h);
  AV   = (void *)malloc(sizeof(float)*num_h);
  AWSPD= (void *)malloc(sizeof(float)*num_h);
  AWDRC= (void *)malloc(sizeof(float)*num_h);
  ADIV = (void *)malloc(sizeof(float)*num_h);
  ADEF = (void *)malloc(sizeof(float)*num_h);
  AHGT = (void *)malloc(sizeof(float)*num_h);
  SU   = (void *)malloc(sizeof(float)*num_h);
  SV   = (void *)malloc(sizeof(float)*num_h);
  SWSPD= (void *)malloc(sizeof(float)*num_h);
  SWDRC= (void *)malloc(sizeof(float)*num_h);
  SDIV = (void *)malloc(sizeof(float)*num_h);
  SDEF = (void *)malloc(sizeof(float)*num_h);

  /* Radar altitude */
  //k = nc_inq_varid(ncid,"altitude",&varid);
  //k = nc_get_var_double(ncid,varid,&r_alt);
  r_alt = 49;
  /* Azimuth */
  k = nc_inq_varid(ncid,"azimuth",&varid);
  k = nc_get_var_float(ncid,varid,az);
  /* Elevation */
  //k = nc_inq_varid(ncid,"fixed_angle",&varid);
  //k = nc_get_var_float(ncid,varid,fixed_angle);
  k = nc_inq_varid(ncid,"elevation",&varid);
  k = nc_get_var_float(ncid,varid,el);
  /* Range from Radar */
  k = nc_inq_varid(ncid,"range",&varid);
  k = nc_get_var_float(ncid,varid,range);
  /* Time */
  k = nc_inq_varid(ncid,"time",&varid);
  k = nc_get_var_float(ncid,varid,time);
  /* sweep start ray index */
  //k = nc_inq_varid(ncid,"sweep_start_ray_index",&varid);
  //k = nc_get_var_int(ncid,varid,sweep_start);
  /* sweep end ray index */
  //k = nc_inq_varid(ncid,"sweep_end_ray_index",&varid);
  //k = nc_get_var_int(ncid,varid,sweep_end);

  /* Read Valiable */
  k = nc_inq_varid(ncid,parav,&varid); 
  k = nc_get_var_float(ncid,varid,vels1);

  k = nc_inq_varid(ncid,paras,&varid); //assuming SNR
  k = nc_get_var_float(ncid,varid,sdvels1);

  k = nc_inq_varid(ncid,"reflectivity",&varid); 
  k = nc_get_var_float(ncid,varid,refs1);
  k = nc_inq_varid(ncid,"temperature",&varid); 
  k = nc_get_var_double(ncid,varid,temps1);
  k = nc_inq_varid(ncid,"dry_air_density",&varid); 
  k = nc_get_var_double(ncid,varid,rhods1);

  nc_close(ncid);

  sweep_start[0] = 0;
  sweep_end[0] = tot_numray - 1;

  for(k=0;k<tot_numray;k++){
    az[k] = az[k]-az_offset;
    if(az[k]<0){az[k]=360+az[k];}
    for(i=0;i<max_numbin;i++){
      //bins = k+i*tot_numray;
      bins = i+k*max_numbin;
      vels[i+k*max_numbin]= vels1[bins];
      sdvels[i+k*max_numbin]= sdvels1[bins];
      refs[i+k*max_numbin]= refs1[bins];
      temps[i+k*max_numbin]= (float)temps1[bins];
      rhods[i+k*max_numbin]= (float)rhods1[bins];
     if(vels1[bins]>9999){      vels[i+k*max_numbin]=-9999;}
     if(sdvels1[bins]>9999){    sdvels[i+k*max_numbin]=-9999;}
     if(refs1[bins]>9999){    refs[i+k*max_numbin]=-9999;}
    }
  }

  wfp = fopen("vad_uv.txt","w");


  for(k=0;k<sweep;k++){

    numray = sweep_end[k] - sweep_start[k] +1;
    time2[k] = (float)time[(sweep_start[k])];

    n=0;
    
    for(i=0;i<max_numbin;i++){
      r = range[i]*0.001;
      sumel = 0;
      sumz = 0;
      sumxy = 0;
      //      Vn  =(void *)malloc(sizeof(double)*numray); 
      for(l=0;l<PARA_NUM;l++){
	for(m=0;m<PARA_NUM;m++){
	  Smm[m][l] = 0.0;
	  Pnm[m] = 0.0;
	}}

      for(j=sweep_start[k];j<=sweep_end[k];j++){
	if(az[j] < 0.0){continue;}
	if(az[j] > 360.0){continue;}
	Az = az[j] / 180.0 * Pi;
	El = el[j] / 180.0 * Pi;
	
	marumi(Az, El, r, &Maru);
	x = Maru.dx_km;
	y = Maru.dy_km;
	z = Maru.dz_km + r_alt * 0.001;
	xy = sqrt(x*x + y*y);
	sumel = sumel + El;
	sumz = sumz + z;
	sumxy = sumxy + xy;
	//if(z < z1 || z >= z2){continue;}

	if(vels[j*max_numbin+i] <= -100 || vels[j*max_numbin+i]>100){continue;}

	vel = vels[j*max_numbin+i] * sign * (-1);
	sdvel = sdvels[j*max_numbin+i];

	if(sdvel < SDVEL_THRE){continue;}
   //== Fall speed based on modified Giangrande et al. (2013)      
   Vfi=0.0;
   //- Case 1: convective, T>0C                              
   if(temps[j*max_numbin+i]>=0.0 && refs[j*max_numbin+i]>40.0){
     Vfi = 2.65 *  pow( pow(10.0,refs[j*max_numbin+i]/10.0), 0.098);
   }
   //-Caase 2: convective, T<0, graupel relationship            
   if(temps[j*max_numbin+i]<0.0 && refs[j*max_numbin+i]>33.0 && refs[j*max_numbin+i]<40.0){
     Vfi = 2.2 + pow( pow(10.0,(refs[j*max_numbin+i]-33.0)/10.0), 0.5);
   }
   if(temps[j*max_numbin+i]<0.0 && refs[j*max_numbin+i]>40.0){
     Vfi = 2.65 *  pow( pow(10.0,refs[j*max_numbin+i]/10.0), 0.098);
   }
   //- Case 3: stratiform, T>0C                               
   if(temps[j*max_numbin+i]>=0.0 && refs[j*max_numbin+i]<40.0){
     Vfi = 3.15 *  pow( pow(10.0,refs[j*max_numbin+i]/10.0), 0.098);
   }
   //- Case 4: snow, T<0C                                               
   if(temps[j*max_numbin+i]<0.0 && refs[j*max_numbin+i]<33.0){
     Vfi = 0.37 *  pow( pow(10.0,refs[j*max_numbin+i]/10.0), 0.19);
   }

   // adjust for air density         
   rhod=RHO0;
   if(rhods[j*max_numbin+i]>0 && rhods[j*max_numbin+i]<2){
     rhod=rhods[j*max_numbin+i];
   }
   Vfi=Vfi*pow(RHO0/rhod,0.4);
   vel = vel-(Vfi*sin(El));
   Vf = Vfi;

    Pnm[0] = Pnm[0] + vel;
    
    Pnm[1] = Pnm[1] + (vel * sin(Az));
    
    Pnm[2] = Pnm[2] + (vel * cos(Az));
    
    Pnm[3] = Pnm[3] + (vel * sin(2*Az));
    
    Pnm[4] = Pnm[4] + (vel * cos(2*Az));
    
/*-- Gmm --*/
	Smm[0][0]=Smm[0][0] + 1;
	Smm[1][1]=Smm[1][1] + pow(sin(Az),2);
	Smm[2][2]=Smm[2][2] + pow(cos(Az),2);
	Smm[3][3]=Smm[3][3] + pow(sin(2*Az),2);
	Smm[4][4]=Smm[4][4] + pow(cos(2*Az),2);
	Smm[0][1]=Smm[0][1] + sin(Az);
	Smm[1][0]=Smm[0][1];
	Smm[0][2]=Smm[0][2] + cos(Az);
	Smm[2][0]=Smm[0][2];
	Smm[0][3]=Smm[0][3] + sin(2*Az);
	Smm[3][0]=Smm[0][3];
	Smm[0][4]=Smm[0][4] + cos(2*Az);
	Smm[4][0]=Smm[0][4];
	Smm[1][2]=Smm[1][2] + (sin(Az)*cos(Az));
	Smm[2][1]=Smm[1][2];
	Smm[1][3]=Smm[1][3] + (sin(Az)*sin(2*Az));
	Smm[3][1]=Smm[1][3];
	Smm[1][4]=Smm[1][4] + (sin(Az)*cos(2*Az));
	Smm[4][1]=Smm[1][4];
	Smm[2][3]=Smm[2][3] + (cos(Az)*sin(2*Az));
	Smm[3][2]=Smm[2][3];
	Smm[2][4]=Smm[2][4] + (cos(Az)*cos(2*Az));
	Smm[4][2]=Smm[2][4];
	Smm[3][4]=Smm[3][4] + (sin(2*Az)*cos(2*Az));
	Smm[4][3]=Smm[3][4];



		    n++;
		    if(n>=numray){break;}


      } // loop j

      aveel = sumel / (float)numray;
      avez = sumz / (float)numray;
      avexy = sumxy / (float)numray;
      data_num = n;

	for(m=0;m<PARA_NUM;m++) {
	  Km[m]=0;
	}
	//      if(data_num <= PARA_NUM){
      if(data_num < 100){
	for(m=0;m<PARA_NUM;m++) {
	  Km[m]=ERROR_UNDEF;
	}
      }else{





/*====================*/
/*====== Gauss elimination =======*/
/*====================*/
		  
/*-- matrix Gm(m+1) --*/
		  for(m=0;m<PARA_NUM;m++) {
		      for(n=0;n<(PARA_NUM+1);n++) {
			  if(n==PARA_NUM){
			      Gmm1[m][n]=Pnm[m];
			  }else{
			      Gmm1[m][n]=Smm[m][n];
			  }
		      }
		  }
/*--------------------------*/		  
/*-- Gauss Jordan ----------*/
		  
/*-- pivot  --*/
		  flag=0;
		  for(m=0;m<PARA_NUM;m++) {

		    max = 0;
		    pivot = m;
		      for(l=m;l<PARA_NUM;l++){
			/*chose max raw in m*/
			  if(max < fabs(Gmm1[l][m])){
			    max = fabs(Gmm1[l][m]);
			    pivot = l;
			  }
		      }

		      if(pivot != m){
			 for(n=0;n<(PARA_NUM+1);n++){
			   tmp1[0][n]=Gmm1[m][n];
			   Gmm1[m][n]=Gmm1[pivot][n];
			   Gmm1[pivot][n]=tmp1[0][n];
			 }
		      }
		  }


		  flag=0;
		  for(m=0;m<PARA_NUM;m++) {
		    G00 = Gmm1[m][m];  // keep diagonal matrix element
		    Gmm1[m][m] = 1;    // diagonal matrix element = 1
/*-- flag  --*/
		      if(G00 == 0.0){
			  flag=1;
			  break;
		      }

		    for(n=m+1;n<(PARA_NUM+1);n++){
		      Gmm1[m][n] /= G00;
		    }

		    for(l=m+1;l<PARA_NUM;l++){
		      G11=Gmm1[l][m];
		      for(n=m+1;n<(PARA_NUM+1);n++){
			Gmm1[l][n] -= G11*Gmm1[m][n];
		      }

		      Gmm1[l][m] = 0; 
		    }
		  }


		  //		  printf("\n %d %d flag=%d data_num=%d\n",i,k,flag,data_num);
		  if(flag == 1){
		      for(m=0;m<PARA_NUM;m++){
			  Km[m]=ERROR_UNDEF;
		      }
		  }else{
		    if(flag == 0){
		      /*calculate Km*/
		      for(m=PARA_NUM-1;m>=0;m--) {
			Km[m] = Gmm1[m][PARA_NUM];
			for(n=PARA_NUM-1;n>m;n--){
			  Km[m] -= Gmm1[m][n]*Km[n];
			}
			//			printf("%f ",Km[m]);
		      }
		      //		      printf("\n");
		    }
		  }
    
      } // if data_num <= PARA_NUM



      K1[k*max_numbin+i] = Km[0];
      K2[k*max_numbin+i] = Km[1];
      K3[k*max_numbin+i] = Km[2];
      K4[k*max_numbin+i] = Km[3];
      K5[k*max_numbin+i] = Km[4];
      //K6[k*max_numbin+i] = 0;
      //K6[k*max_numbin+i] = Km[5];
      
   
/*------------------------------------------------*/		  
 /*-- Calculation of wind speed and direction --*/
      if(Km[1] != ERROR_UNDEF && Km[2] != ERROR_UNDEF){
	u0 = Km[1] / cos(aveel) * (-1);
	v0 = Km[2] / cos(aveel) * (-1);
	ws = sqrt(u0 * u0 + v0 * v0);
	ws = sqrt(Km[1]*Km[1] + Km[2]*Km[2])/cos(aveel);
	u0d = u0 * -1;
	v0d = v0 * -1;
	if(v0d >= 0.0 && u0d >=0){
	  wd = acos(fabs(v0d)/ws) * 180 / Pi;
	}
	if(v0d < 0.0 && u0d >=0){
	  //		      wd = acos(v0d/ws) * 180 / Pi;
	  wd = 180 - acos(fabs(v0)/ws) * 180 / Pi;
	}
	if(v0d < 0.0 && u0d <0){
	  //		      wd = 360 - acos(v0d/ws) * 180 / Pi;
	  wd = 180 + acos(fabs(v0)/ws) * 180 / Pi;
	}
	if(v0d >= 0.0 && u0d <0){
	  wd = 360 - acos(fabs(v0d)/ws) * 180 / Pi;
	}
	/*
	if(Km[2] >=0){
	  wd = atan(Km[1]/Km[2]);
	}else{
	  wd = atan(Km[1]/Km[2]) + Pi;
	}
	wd = wd * 180.0/Pi;
	*/
	fprintf(wfp,"%f %f %f %f %f %f %d \n",avez,u0,v0,ws,wd,aveel,data_num);
      }else{
	u0 = ERROR_UNDEF;
	v0 = ERROR_UNDEF;
	ws = ERROR_UNDEF;
	wd = ERROR_UNDEF;
      }
/*------------------------------------------------*/		  
 /*-- Calculation of divergence etc --*/
      if(Km[0] != ERROR_UNDEF){
	divg = 2*(Vf*tan(aveel) - Km[0]/cos(aveel))/range[i];
      }else{
	divg = ERROR_UNDEF;
      }

      if(Km[3] != ERROR_UNDEF && Km[4] != ERROR_UNDEF){
	defm = 2 * sqrt(Km[3]*Km[3]+Km[4]*Km[4]) / cos(aveel) / range[i];
      }else{
	defm = ERROR_UNDEF;
      }

      U[k*max_numbin+i] = u0;
      V[k*max_numbin+i] = v0;
      WSPD[k*max_numbin+i] = ws;
      WDRC[k*max_numbin+i] = wd;
      DIV[k*max_numbin+i] = divg;
      DEF[k*max_numbin+i] = defm;

      HGT[k*max_numbin+i] = avez;
      RADIUS[k*max_numbin+i] = avexy;

    } //loop i


    /*Variance of Doppler velocity along ray*/
    for(i=0;i<max_numbin;i++){
	aa = i - SDEV_RANGE;
	if(aa < 0){ aa = 0;}
	aaa = i + SDEV_RANGE;
	if(aaa >= max_numbin){ aaa = max_numbin - 1;}
	sumu=0;sumv=0;sumdiv=0;sumdef=0;
	counta = 0;
	for(a=aa;a<aaa;a++){
	  if(U[k*max_numbin+a] > ERROR_UNDEF){
	    sumu = sumu + U[k*max_numbin+a];
	    sumv = sumv + V[k*max_numbin+a];
	    sumdiv = sumdiv + DIV[k*max_numbin+a];
	    sumdef = sumdef + DEF[k*max_numbin+a];
	    counta ++;
	  }
	}
	if(counta > 2){
	  aveu = sumu / counta;
	  avev = sumv / counta;
	  avediv = sumdiv / counta;
	  avedef = sumdef / counta;
	sumu=0;sumv=0;sumdiv=0;sumdef=0;
	  for(a=aa;a<aaa;a++){
	    if(U[k*max_numbin+a] > ERROR_UNDEF){
	      sumu = sumu + pow((U[k*max_numbin+a] - aveu),2);
	      sumv = sumv + pow((V[k*max_numbin+a] - avev),2);
	      sumdiv = sumdiv + pow((DIV[k*max_numbin+a] - avediv),2);
	      sumdef = sumdef + pow((DEF[k*max_numbin+a] - avedef),2);
	    }
	  }
	  sdevu = sqrt(sumu / counta);
	}else{
	  aveu = ERROR_UNDEF;
	  avev = ERROR_UNDEF;
	  avediv = ERROR_UNDEF;
	  avedef = ERROR_UNDEF;
	  sdevu = ERROR_UNDEF;
	}
	if(sdevu > SDEV_THRESHOLD || sdevu ==ERROR_UNDEF || fabs(U[k*max_numbin+i])>100 || fabs(V[k*max_numbin+i])>100 || sdevu <=0 ){
	  Ufil2[k*max_numbin+i] = ERROR_UNDEF;
	  Vfil2[k*max_numbin+i] = ERROR_UNDEF;
	  WSPDfil2[k*max_numbin+i] = ERROR_UNDEF;
	  WDRCfil2[k*max_numbin+i] = ERROR_UNDEF;
	  DIVfil2[k*max_numbin+i] = ERROR_UNDEF;
	  DEFfil2[k*max_numbin+i] = ERROR_UNDEF;
	}else{
	  Ufil2[k*max_numbin+i] = U[k*max_numbin+i];
	  Vfil2[k*max_numbin+i] = V[k*max_numbin+i];
	  WSPDfil2[k*max_numbin+i] = WSPD[k*max_numbin+i];
	  WDRCfil2[k*max_numbin+i] = WDRC[k*max_numbin+i];
	  DIVfil2[k*max_numbin+i] = DIV[k*max_numbin+i];
	  DEFfil2[k*max_numbin+i] = DEF[k*max_numbin+i];
	  /*
	  Ufil[k*max_numbin+i] = aveu;
	  Vfil[k*max_numbin+i] = avev;
	  DIVfil[k*max_numbin+i] = avediv;
	  DEFfil[k*max_numbin+i] = avedef;
	  WSPDfil[k*max_numbin+i] = sqrt(aveu * aveu + avev * avev);
	  */
	}
    }// loop i 2


    for(i=0;i<max_numbin;i++){

	aa = i - SDEV_RANGE;	if(aa < 0){ aa = 0;}
	aaa = i + SDEV_RANGE;	if(aaa >= max_numbin){ aaa = max_numbin - 1;}
	sumu=0;sumv=0;sumdiv=0;sumdef=0;
	counta = 0;
	for(a=aa;a<aaa;a++){
	  if(Ufil2[k*max_numbin+a] > ERROR_UNDEF && Vfil2[k*max_numbin+a] > ERROR_UNDEF){
	    sumu = sumu + Ufil2[k*max_numbin+a];
	    sumv = sumv + Vfil2[k*max_numbin+a];
	    sumdiv = sumdiv + DIVfil2[k*max_numbin+a];
	    sumdef = sumdef + DEFfil2[k*max_numbin+a];
	    counta ++;
	  }
	}
	if(counta > 0){
	  aveu = sumu / counta;
	  avev = sumv / counta;
	  avediv = sumdiv / counta;
	  avedef = sumdef / counta;
	  Ufil[k*max_numbin+i] = aveu;
	  Vfil[k*max_numbin+i] = avev;
	  DIVfil[k*max_numbin+i] = avediv;
	  DEFfil[k*max_numbin+i] = avedef;
	  WSPDfil[k*max_numbin+i] = sqrt(aveu * aveu + avev * avev);

	u0d = aveu * -1;
	v0d = avev * -1;
	if(v0d >= 0.0 && u0d >=0){
	  WDRCfil[k*max_numbin+i] = acos(fabs(v0d)/WSPDfil[k*max_numbin+i]) * 180 / Pi;
	}
	if(v0d < 0.0 && u0d >=0){
	  WDRCfil[k*max_numbin+i] = 180 - acos(fabs(avev)/WSPDfil[k*max_numbin+i]) * 180 / Pi;
	}
	if(v0d < 0.0 && u0d <0){
	  WDRCfil[k*max_numbin+i] = 180 + acos(fabs(avev)/WSPDfil[k*max_numbin+i]) * 180 / Pi;
	}
	if(v0d >= 0.0 && u0d <0){
	   WDRCfil[k*max_numbin+i] = 360 - acos(fabs(v0d)/WSPDfil[k*max_numbin+i]) * 180 / Pi;
	}
	}else{
	  Ufil[k*max_numbin+i] = ERROR_UNDEF;
	  Vfil[k*max_numbin+i] = ERROR_UNDEF;
	  WSPDfil[k*max_numbin+i] = ERROR_UNDEF;
	  WDRCfil[k*max_numbin+i] = ERROR_UNDEF;
	  DIVfil[k*max_numbin+i] = ERROR_UNDEF;
	  DEFfil[k*max_numbin+i] = ERROR_UNDEF;
	}

  } // loop i 3


		  
  } // loop k  


fclose(wfp);


 j=0;
 for(h=0;h<hmax;h=h+dh){
   hc = h+dh/2;
   sumu=0;sumv=0;sums=0;sumd=0;sumdiv=0;sumdef=0;
   count=0;
   //for(k=1;k<(sweep-1);k++){
   for(k=0;k<sweep;k++){
   for(i=0;i<max_numbin;i++){
     if(Ufil[k*max_numbin+i]>ERROR_UNDEF){
     if(HGT[k*max_numbin+i]>=h && HGT[k*max_numbin+i]<(h+dh)){
       sumu = sumu+Ufil[k*max_numbin+i];
       sumv = sumv+Vfil[k*max_numbin+i];
       sums = sums+WSPDfil[k*max_numbin+i];
       sumd = sumd+WDRCfil[k*max_numbin+i];
       sumdiv = sumdiv+DIVfil[k*max_numbin+i];
       sumdef = sumdef+DEFfil[k*max_numbin+i];
       count ++;
     }
     }
   }
  }
  if(count >0){
    aveu = sumu/count;
    avev = sumv/count;
    aves = sums/count;
    aved = sumd/count;
    avediv = sumdiv/count;
    avedef = sumdef/count;
    sumu=0;sumv=0;sums=0;sumd=0;sumdiv=0;sumdef=0;
    count=0;
    for(k=0;k<sweep;k++){
      for(i=0;i<max_numbin;i++){
	if(Ufil[k*max_numbin+i]>ERROR_UNDEF){
	  if(HGT[k*max_numbin+i]>=h && HGT[k*max_numbin+i]<(h+dh)){
	    sumu = sumu + pow((Ufil[k*max_numbin+i] - aveu),2);
	    sumv = sumv + pow((Vfil[k*max_numbin+i] - avev),2);
	    sums = sums + pow((WSPDfil[k*max_numbin+i] - aves),2);
	    sumd = sumd + pow((WDRCfil[k*max_numbin+i] - aved),2);
	    sumdiv = sumdiv + pow((DIVfil[k*max_numbin+i] - avediv),2);
	    sumdef = sumdef + pow((DEFfil[k*max_numbin+i] - avedef),2);
	    count ++;
	  }
	}
	sdevu = sqrt(sumu / count);
	sdevv = sqrt(sumv / count);
	sdevs = sqrt(sums / count);
	sdevd = sqrt(sumd / count);
	sdevdiv = sqrt(sumdiv / count);
	sdevdef = sqrt(sumdef / count);
      }
    }
  }else{
    aveu=ERROR_UNDEF;avev=ERROR_UNDEF;aves=ERROR_UNDEF;aved=ERROR_UNDEF;
    avediv=ERROR_UNDEF;avedef=ERROR_UNDEF;
    sdevu=ERROR_UNDEF;sdevv=ERROR_UNDEF;sdevs=ERROR_UNDEF;sdevd=ERROR_UNDEF;
    sdevdiv=ERROR_UNDEF;sdevdef=ERROR_UNDEF;
  }

  AU[j] = aveu; AV[j] = avev; AWSPD[j] = aves; AWDRC[j] = aved; ADIV[j] = avediv; ADEF[j] = avedef;

  SU[j] = sdevu; SV[j] = sdevv; SWSPD[j] = sdevs; SWDRC[j] = sdevd; SDIV[j] = sdevdiv; SDEF[j] = sdevdef;

  AHGT[j] = hc;
  j=j+1;
 }


 fixed_angle[0]=aveel;



  /*== NETCDF FILE OPEN for WRITE ==*/
  nc_error = nc_create("vad.nc",NC_CLOBBER,&ncid);
  if(nc_error != NC_NOERR){
    printf("OPEN ERROR vad.nc\n");
    exit(1);
  }
  undef = ERROR_UNDEF;

  /*== GLOBAL ATTRIBUTION ==*/
  k=nc_put_att_text(ncid,NC_GLOBAL,"Convention",40,"XSAPR Radar PPI VAD");

  /*== DIMENSION REGISTRATION ==*/
  nc_def_dim(ncid, "range",max_numbin, &range_dimid);
  nc_def_dim(ncid, "sweep",sweep, &sweep_dimid);
  nc_def_dim(ncid, "mheight",num_h, &mheight_dimid);
  dimid2s[1]=range_dimid;
  dimid2s[0]=sweep_dimid;
  dimidrange[0]=range_dimid;
  dimidsweep[0]=sweep_dimid;
  dimidmheight[0]=mheight_dimid;

  nc_def_var(ncid, "FallSpeed",  NC_FLOAT, 0, 0, &varidvf);
  nc_def_var(ncid, "range",          NC_FLOAT, 1, dimidrange, &varidrange);
  nc_def_var(ncid, "fixed_angle", NC_FLOAT, 1, dimidsweep, &varidfixed_angle);
  nc_def_var(ncid, "time", NC_FLOAT, 1, dimidsweep, &varidt);
 
  nc_def_var(ncid, "HEIGHT", NC_FLOAT, 2, dimid2s, &varidhgt);
  nc_put_att_float(ncid, varidhgt, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "RADIUS", NC_FLOAT, 2, dimid2s, &varidrad);
  nc_put_att_float(ncid, varidrad, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "U0", NC_FLOAT, 2, dimid2s, &varidu);
  nc_put_att_float(ncid, varidu, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "V0", NC_FLOAT, 2, dimid2s, &varidv);
  nc_put_att_float(ncid, varidv, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "WindSpeed", NC_FLOAT, 2, dimid2s, &varidws);
  nc_put_att_float(ncid, varidws, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "WindDirection", NC_FLOAT, 2, dimid2s, &varidwd);
  nc_put_att_float(ncid, varidwd, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "Divergence", NC_FLOAT, 2, dimid2s, &variddiv);
  nc_put_att_float(ncid, variddiv, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "Deformation", NC_FLOAT, 2, dimid2s, &variddef);
  nc_put_att_float(ncid, variddef, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "RADIUS", NC_FLOAT, 2, dimid2s, &varidrad);
  nc_put_att_float(ncid, varidrad, "_FillValue", NC_FLOAT, 1, &undef);


  nc_def_var(ncid, "U0_Fil", NC_FLOAT, 2, dimid2s, &varidu2);
  nc_put_att_float(ncid, varidu2, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "V0_Fil", NC_FLOAT, 2, dimid2s, &varidv2);
  nc_put_att_float(ncid, varidv2, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "WindSpeed_Fil", NC_FLOAT, 2, dimid2s, &varidws2);
  nc_put_att_float(ncid, varidws2, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "WindDirection_Fil", NC_FLOAT, 2, dimid2s, &varidwd2);
  nc_put_att_float(ncid, varidwd2, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "Divergence_Fil", NC_FLOAT, 2, dimid2s, &variddiv2);
  nc_put_att_float(ncid, variddiv2, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "Deformation_Fil", NC_FLOAT, 2, dimid2s, &variddef2);
  nc_put_att_float(ncid, variddef2, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "K1", NC_FLOAT, 2, dimid2s, &varidk1);
  nc_put_att_float(ncid, varidk1, "_FillValue", NC_FLOAT, 1, &undef);
  nc_def_var(ncid, "K2", NC_FLOAT, 2, dimid2s, &varidk2);
  nc_put_att_float(ncid, varidk2, "_FillValue", NC_FLOAT, 1, &undef);
  nc_def_var(ncid, "K3", NC_FLOAT, 2, dimid2s, &varidk3);
  nc_put_att_float(ncid, varidk3, "_FillValue", NC_FLOAT, 1, &undef);
  nc_def_var(ncid, "K4", NC_FLOAT, 2, dimid2s, &varidk4);
  nc_put_att_float(ncid, varidk4, "_FillValue", NC_FLOAT, 1, &undef);
  nc_def_var(ncid, "K5", NC_FLOAT, 2, dimid2s, &varidk5);
  nc_put_att_float(ncid, varidk5, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "MEAN_HEIGHT", NC_FLOAT, 1, dimidmheight, &varidhgt3);
  nc_put_att_float(ncid, varidhgt3, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "MEAN_U0", NC_FLOAT, 1, dimidmheight, &varidu3);
  nc_put_att_float(ncid, varidu3, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "MEAN_V0", NC_FLOAT, 1, dimidmheight, &varidv3);
  nc_put_att_float(ncid, varidv3, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "MEAN_WindSpeed", NC_FLOAT, 1, dimidmheight, &varidws3);
  nc_put_att_float(ncid, varidws3, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "MEAN_WindDirection", NC_FLOAT, 1, dimidmheight, &varidwd3);
  nc_put_att_float(ncid, varidwd3, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "MEAN_Divergence", NC_FLOAT, 1, dimidmheight, &variddiv3);
  nc_put_att_float(ncid, variddiv3, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "MEAN_Deformation", NC_FLOAT, 1, dimidmheight, &variddef3);
  nc_put_att_float(ncid, variddef3, "_FillValue", NC_FLOAT, 1, &undef);


  nc_def_var(ncid, "SD_U0", NC_FLOAT, 1, dimidmheight, &varidu4);
  nc_put_att_float(ncid, varidu4, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "SD_V0", NC_FLOAT, 1, dimidmheight, &varidv4);
  nc_put_att_float(ncid, varidv4, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "SD_WindSpeed", NC_FLOAT, 1, dimidmheight, &varidws4);
  nc_put_att_float(ncid, varidws4, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "SD_WindDirection", NC_FLOAT, 1, dimidmheight, &varidwd4);
  nc_put_att_float(ncid, varidwd4, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "SD_Divergence", NC_FLOAT, 1, dimidmheight, &variddiv4);
  nc_put_att_float(ncid, variddiv4, "_FillValue", NC_FLOAT, 1, &undef);

  nc_def_var(ncid, "SD_Deformation", NC_FLOAT, 1, dimidmheight, &variddef4);
  nc_put_att_float(ncid, variddef4, "_FillValue", NC_FLOAT, 1, &undef);

  k=nc_enddef(ncid);


  /*== PUT ==*/
  nc_put_var_float(ncid, varidvf, &Vf);
  nc_put_var_float(ncid, varidrange, range);
  nc_put_var_float(ncid, varidfixed_angle, fixed_angle);
  nc_put_var_float(ncid, varidt, time2);

  nc_put_var_float(ncid, varidhgt, HGT);
  nc_put_var_float(ncid, varidrad, RADIUS);
  nc_put_var_float(ncid, varidu, U);
  nc_put_var_float(ncid, varidv, V);
  nc_put_var_float(ncid, varidws, WSPD);
  nc_put_var_float(ncid, varidwd, WDRC);
  nc_put_var_float(ncid, variddiv, DIV);
  nc_put_var_float(ncid, variddef, DEF);
  nc_put_var_float(ncid, varidu2, Ufil);
  nc_put_var_float(ncid, varidv2, Vfil);
  nc_put_var_float(ncid, varidws2, WSPDfil);
  nc_put_var_float(ncid, varidwd2, WDRCfil);
  nc_put_var_float(ncid, variddiv2, DIVfil);
  nc_put_var_float(ncid, variddef2, DEFfil);
  nc_put_var_float(ncid, varidk1, K1);
  nc_put_var_float(ncid, varidk2, K2);
  nc_put_var_float(ncid, varidk3, K3);
  nc_put_var_float(ncid, varidk4, K4);
  nc_put_var_float(ncid, varidk5, K5);
  nc_put_var_float(ncid, varidhgt3, AHGT);
  nc_put_var_float(ncid, varidu3, AU);
  nc_put_var_float(ncid, varidv3, AV);
  nc_put_var_float(ncid, varidws3, AWSPD);
  nc_put_var_float(ncid, varidwd3, AWDRC);
  nc_put_var_float(ncid, variddiv3, ADIV);
  nc_put_var_float(ncid, variddef3, ADEF);
  nc_put_var_float(ncid, varidu4, SU);
  nc_put_var_float(ncid, varidv4, SV);
  nc_put_var_float(ncid, varidws4, SWSPD);
  nc_put_var_float(ncid, varidwd4, SWDRC);
  nc_put_var_float(ncid, variddiv4, SDIV);
  nc_put_var_float(ncid, variddef4, SDEF);

  nc_close(ncid);
  
}
