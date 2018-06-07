/* 3D elastic modeling (mu=0) with 2Nth order staggered-grid FD */

#include <rsf.h>
#include <math.h>

static int nb,nz,nx,ny,nzpad,nxpad,nypad,nt,N;
static float dz,dx,dy,dt,fm;

/*< expand domain of 'a' to 'b': source(a)-->destination(b) >*/
void expand3d(float ***b,float ***a){
	int iz,ix,i3;
	for(i3=0;i3<ny;i3++){
		for(ix=0;ix<nx;ix++){
			for(iz=0;iz<nz;iz++){
			b[nb+i3][nb+ix][nb+iz]=a[i3][ix][iz];
			}
		}
	}
	for(i3=0;i3<nypad;i3++){
		for(ix=0; ix<nxpad;ix++){
			for(iz=0;iz<nb;iz++){
				b[i3][ix][iz]=b[i3][ix][nb];
				b[i3][ix][nzpad-iz-1]=b[i3][ix][nzpad-nb-1];
			}
		}
	}


	for(i3=0;i3<nypad;i3++){
		for(ix=0;ix<nb;ix++){
			for (iz=0; iz<nzpad;iz++){
				b[i3][ix][iz]=b[i3][nb][iz];
				b[i3][nxpad-ix-1][iz]=b[i3][nxpad-nb-1][iz];
			}
		}
	}

	for(i3=0;i3<nb;i3++){
		for(ix=0;ix<nxpad;ix++){
			for(iz=0;iz<nzpad;iz++){
				b[i3][ix][iz]=b[nb][ix][iz];
				b[nypad-i3-1][ix][iz]=b[nypad-nb-1][ix][iz];
			}
		}
	}
}

/*< window 'b' to 'a': source(b)-->destination(a) >*/
void window3d(float ***a, float ***b){
	
	int i3,ix,iz;
	for(i3=0;i3<ny;i3++){
		for(ix=0;ix<nx;ix++){
			for(iz=0;iz<nz;iz++){
				a[i3][ix][iz]=b[nb+i3][nb+ix][nb+iz];
			}
		}
	}
}

/*< apply absorbing boundary condition >*/
void apply_sponge(float ***uu, float *spo)
{
	int iz,ix,iy,ib,ibz,ibx,iby;
	float w;

	for(ib=0;ib<nb;ib++){
		w=spo[ib];
		ibz=nzpad-ib-1;
		for(iy=0;iy<nypad;iy++){
			for(ix=0;ix<nxpad;ix++){
				uu[iy][ix][ib]*=w; /* z min */
				uu[iy][ix][ibz]*=w; /* z max */
			}
		}
		ibx=nxpad-ib-1;
		for(iy=0;iy<nypad;iy++){
			for(iz=0;iz<nzpad;iz++){
				uu[iy][ib][iz]*=w; /* x min */
				uu[iy][ibx][iz]*=w; /* x max */
			}
		}
		iby=nypad-ib-1;
		for(ix=0;ix<nxpad;ix++) {
			for(iz=0;iz<nzpad;iz++){
				uu[ib ][ix][iz]*=w; /* y min */
				uu[iby][ix][iz]*=w; /* y max */
			}
		}

	}
}

/* compute the coefficient */

float* compute_coef(int N){
	
	int n,m;
	float tempc;
	float *cc=calloc(N+1,sizeof(float));
	for(n=1;n<=N;n++){
		tempc=1;
		for(m=1;m<=N;m++){
			if(m!=n){
				tempc*=(2*(float)m-1)*(2*(float)m-1)/((2*(float)n-1)*(2*(float)n-1)-(2*(float)m-1)*(2*(float)m-1));
			}
		}
		if(tempc<0) tempc=0-tempc;
		cc[n]=(pow((float)(-1),(float)(n+1))/(2*(float)n-1))*tempc;
	}
	return cc;
}

/*< forward modeling step >*/
void step_forward3d(float ***pxx, float ***pzz, float ***pyy, float ***vz, float ***vx, float ***vy, float ***k, float ***rho, float *coef){

	int iz,ix,iy,iop;
	float tmp1,tmp2,tmp3,tmpx;
	for(ix=N;ix<nxpad-N-1;ix++){
		for(iz=N;iz<nzpad-N;iz++){
			for(iy=N;iy<nypad-N;iy++){
				tmpx=0;
				for(iop=1;iop<N+1;iop++){
					tmpx=tmpx+(pxx[iy][ix+iop][iz]-pxx[iy][ix-iop+1][iz])*coef[iop];
				}
				vx[iy][ix][iz]=vx[iy][ix][iz]+tmpx*dt/rho[iy][ix][iz]*(1.0/dx);
			}
		}
	}
	for(ix=N;ix<nxpad-N;ix++){
		for(iz=N;iz<nzpad-N-1;iz++){
			for(iy=N;iy<nypad-N;iy++){
				tmpx=0;
				for(iop=1;iop<N+1;iop++){
					tmpx=tmpx+(pzz[iy][ix][iz+iop]-pzz[iy][ix][iz-iop+1])*coef[iop];
				}
				vz[iy][ix][iz]=vz[iy][ix][iz]+tmpx*dt/rho[iy][ix][iz]*(1.0/dz);
			}
		}
	}
	for(ix=N;ix<nxpad-N;ix++){
		for(iz=N;iz<nzpad-N;iz++){
			for(iy=N;iy<nypad-N-1;iy++){
				tmpx=0;
				for(iop=1;iop<N+1;iop++){
					tmpx=tmpx+(pyy[iy+iop][ix][iz]-pyy[iy-iop+1][ix][iz])*coef[iop];
				}
				vy[iy][ix][iz]=vy[iy][ix][iz]+tmpx*dt/rho[iy][ix][iz]*(1.0/dy);
			}
		}
	}
	for(ix=N;ix<nxpad-N;ix++){
		for(iz=N;iz<nzpad-N;iz++){
			for(iy=N;iy<nypad-N;iy++){
				tmp1=0;
				tmp2=0;
				tmp3=0;
				for(iop=1;iop<N+1;iop++){
					tmp1=tmp1+(vx[iy][ix+iop-1][iz]-vx[iy][ix-iop][iz])*coef[iop];
				}
				for(iop=1;iop<N+1;iop++){
					tmp2=tmp2+(vz[iy][ix][iz+iop-1]-vz[iy][ix][iz-iop])*coef[iop];
				}
				for(iop=1;iop<N+1;iop++){
					tmp3=tmp3+(vy[iy+iop-1][ix][iz]-vy[iy-iop][ix][iz])*coef[iop];
				}
				pxx[iy][ix][iz]=pxx[iy][ix][iz]+(tmp1*(1.0/dx)+tmp2*(1.0/dz)+tmp3*(1.0/dy))*dt*k[iy][ix][iz];
				pzz[iy][ix][iz]=pzz[iy][ix][iz]+(tmp1*(1.0/dx)+tmp2*(1.0/dz)+tmp3*(1.0/dy))*dt*k[iy][ix][iz];
				pyy[iy][ix][iz]=pyy[iy][ix][iz]+(tmp1*(1.0/dx)+tmp2*(1.0/dz)+tmp3*(1.0/dy))*dt*k[iy][ix][iz];
			}
		}
	}
}

int main(int argc, char* argv[]){

	int it,ib,sx,sz,sy,ix,iz,iy;
	float tmp,strike,dip,slip;
	float *wlt,*bndr,*coef;
	float ***temp,***vel,***rho,***vx,***vz,***vy,***pxx,***pzz,***pyy,***k;
	sf_file Fvel,Frho,Fwfl;
	
	sf_init(argc,argv);   /* initial madagascar */
	
	Fvel=sf_input("in");  /* velocity model     */
	Frho=sf_input("rho"); /* density            */
	Fwfl=sf_output("out"); /* wavefield snapshot */

	/* input the parameters and data */
	if(!sf_histint(Fvel,"n1",&nz))          sf_error("No nz in velocity model");
	if(!sf_histint(Fvel,"n2",&nx))          sf_error("No nx in velocity model");
	if(!sf_histint(Fvel,"n3",&ny))          sf_error("No ny in velocity model");
	if(!sf_histfloat(Fvel,"d1",&dz))        sf_error("No dz in velocity model");
	if(!sf_histfloat(Fvel,"d2",&dx))        sf_error("No dx in velocity model");
	if(!sf_histfloat(Fvel,"d3",&dy))        sf_error("No dy in velocity model");
        if(!sf_getint("nb",&nb))                nb=100;
	if(!sf_getint("nt",&nt))                sf_error("No nt found");
	if(!sf_getfloat("dt",&dt))              sf_error("No dt found");
	if(!sf_getfloat("strike",&strike))      sf_error("No strike found");
	if(!sf_getfloat("dip",&dip))            sf_error("No dip found");
	if(!sf_getfloat("slip",&slip))          sf_error("No slip found");
	if(!sf_getfloat("fm",&fm))              fm=20.0;
	if(!sf_getint("N",&N))                N=2;

	/* initial the output data */
	sf_putint(Fwfl,"n1",nz);
	sf_putint(Fwfl,"n2",nx);
	sf_putint(Fwfl,"n3",ny);
	sf_putint(Fwfl,"n4",nt/30);

	nzpad=nz+2*nb;
	nxpad=nx+2*nb;
	nypad=ny+2*nb;
	sx=nxpad/2;
	sz=nzpad/2;
	sy=nypad/2;



	/* allocate variables */
	wlt  =sf_floatalloc(nt);
	bndr =sf_floatalloc(nb);
	temp =sf_floatalloc3(nz,nx,ny);
	vel  =sf_floatalloc3(nzpad,nxpad,nypad);
	rho  =sf_floatalloc3(nzpad,nxpad,nypad);
	pxx  =sf_floatalloc3(nzpad,nxpad,nypad);
	pzz  =sf_floatalloc3(nzpad,nxpad,nypad);
	pyy  =sf_floatalloc3(nzpad,nxpad,nypad);
	vz   =sf_floatalloc3(nzpad,nxpad,nypad);
	vx   =sf_floatalloc3(nzpad,nxpad,nypad);
	vy   =sf_floatalloc3(nzpad,nxpad,nypad);
	k    =sf_floatalloc3(nzpad,nxpad,nypad);

	/* Ricky wavelet */
	for(it=0;it<nt;it++){

	tmp=SF_PI*fm*(it*dt-1.2/fm);
	tmp*=tmp;
	wlt[it]=(1.0-2.0*tmp)*expf(-tmp);

	}

	/* Boundary */
	for(ib=0;ib<nb;ib++){

	tmp=0.015*(nb-ib);
	bndr[ib]=expf(-tmp*tmp);

	}

	/* initialization */
	sf_floatread(temp[0][0],nz*nx*ny,Fvel);
	expand3d(vel, temp);
	sf_floatread(temp[0][0],nz*nx*ny,Frho);
	expand3d(rho, temp);
	memset(pxx[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(pzz[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(pyy[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(vx[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(vz[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(vy[0][0],0,nzpad*nxpad*nypad*sizeof(float));

	/* k=v^2*rho */
	for(ix=0;ix<nxpad;ix++){
		for(iz=0;iz<nzpad;iz++){
			for(iy=0;iy<nypad;iy++){
			k[iy][ix][iz]=rho[iy][ix][iz]*vel[iy][ix][iz]*vel[iy][ix][iz];
			}
		}
	}
	/* calculate the strike,dip,and slip */
	strike=strike*SF_PI/180;
	dip=dip*SF_PI/180;
	slip=slip*SF_PI/180;
	
	/* calculate the coefficent */
	coef=compute_coef(N);

	/* staggered grid finite difference */
	for(it=0; it<nt; it++)
	{
		if(it%30==0){
			window3d(temp,pxx);
			sf_floatwrite(temp[0][0],nz*nx*ny,Fwfl);
		}
		pxx[sy][sx][sz]+=wlt[it]*cos(slip)*sin(dip)*sin(2*strike)-wlt[it]*sin(slip)*sin(2*dip)*sin(strike)*sin(strike);
		pyy[sy][sx][sz]+=-wlt[it]*cos(slip)*sin(dip)*sin(2*strike)-wlt[it]*sin(slip)*sin(2*dip)*sin(strike)*cos(strike);
		pzz[sy][sx][sz]+=wlt[it]*sin(slip)*sin(2*dip);
		step_forward3d(pxx,pzz,pyy,vz, vx, vy,k, rho,coef);
		apply_sponge(pxx, bndr);
		apply_sponge(pzz, bndr);
		apply_sponge(pyy, bndr);
		apply_sponge(vx, bndr);
		apply_sponge(vz, bndr);
		apply_sponge(vy, bndr);

	}
	
	/* free variables */
	free(wlt);
	free(bndr);
	free(**rho);  free(*rho);  free(rho);
	free(**temp); free(*temp); free(temp);
	free(**vel);  free(*vel);  free(vel);
	free(**pxx);  free(*pxx);  free(pxx);
	free(**pzz);  free(*pzz);  free(pzz);
	free(**pyy);  free(*pyy);  free(pyy);
	free(**vx);   free(*vx);   free(vx);
	free(**vz);   free(*vz);   free(vz);
	free(**vy);   free(*vy);   free(vy);
	free(**k);    free(*k);    free(k);

	
	exit(0);
}


	
	



































	
