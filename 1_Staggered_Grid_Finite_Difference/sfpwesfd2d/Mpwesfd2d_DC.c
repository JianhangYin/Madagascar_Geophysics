/* 2D acoustic modeling wiht 2Nth order staggered-grid FD */
/*
  Copyright (C) 2014 Xi'an Jiaotong University (Pengliang Yang)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>

static int nb,nz,nx,nzpad,nxpad,nt,N;
static float dz,dx,dt,fm;

/*< expand domain of 'a' to 'b': source(a)-->destination(b) >*/
void expand2d(float **b,float **a){
	int iz,ix;
	for(ix=0;ix<nx;ix++){
		for(iz=0;iz<nz;iz++){
	    	b[nb+ix][nb+iz]=a[ix][iz];
		}
    	}

    	for(ix=0;ix<nxpad;ix++){
		for(iz=0;iz<nb;iz++){
	   	b[ix][iz]=b[ix][nb];
	    	b[ix][nzpad-iz-1]=b[ix][nzpad-nb-1];
		}
    	}

    	for(ix=0;ix<nb;ix++){
		for(iz=0;iz<nzpad;iz++){
	    	b[ix][iz]=b[nb][iz];
	    	b[nxpad-ix-1][iz]=b[nxpad-nb-1][iz];
		}
    	}
}

/*< window 'b' to 'a': source(b)-->destination(a) >*/
void window2d(float **a, float **b){
	int iz,ix;
	for(ix=0;ix<nx;ix++){
		for(iz=0;iz<nz;iz++){
		a[ix][iz]=b[nb+ix][nb+iz] ;
		}
	}
}

/*< apply absorbing boundary condition >*/
void apply_sponge(float **u, float *bndr)
{
	int ix,iz;
	for(ix=0; ix<nxpad; ix++){
		for(iz=0;iz<nb;iz++){	/* top ABC */			
		u[ix][iz]=bndr[iz]*u[ix][iz];
		}
		for(iz=nz+nb;iz<nzpad;iz++){/* bottom ABC */			
		u[ix][iz]=bndr[nzpad-iz-1]*u[ix][iz];
		}
	}
	for(iz=0; iz<nzpad; iz++){
		for(ix=0;ix<nb;ix++){	/* left ABC */			
	    	u[ix][iz]=bndr[ix]*u[ix][iz];
		}	
		for(ix=nx+nb;ix<nxpad;ix++){/* right ABC */		
	   	u[ix][iz]=bndr[nxpad-ix-1]*u[ix][iz];
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
void step_forward(float **pxx, float **pzz, float **vz, float **vx, float **k, float **rho, float *coef){

	int iz,ix,iop;
	float tmp1,tmp2,tmpx;
	for(ix=N;ix<nxpad-N-1;ix++){
		for(iz=N;iz<nzpad-N;iz++){
			tmpx=0;
			for(iop=1;iop<N+1;iop++){
				tmpx=tmpx+(pxx[ix+iop][iz]-pxx[ix-iop+1][iz])*coef[iop];
			}
			vx[ix][iz]=vx[ix][iz]+tmpx*dt/rho[ix][iz]*(1.0/dx);
		}
	}
	for(ix=N;ix<nxpad-N;ix++){
		for(iz=N;iz<nzpad-N-1;iz++){
			tmpx=0;
			for(iop=1;iop<N+1;iop++){
				tmpx=tmpx+(pzz[ix][iz+iop]-pzz[ix][iz-iop+1])*coef[iop];
			}
			vz[ix][iz]=vz[ix][iz]+tmpx*dt/rho[ix][iz]*(1.0/dz);
		}
	}
	for(ix=N;ix<nxpad-N;ix++){
		for(iz=N;iz<nzpad-N;iz++){
			tmp1=0;
			tmp2=0;
			for(iop=1;iop<N+1;iop++){
				tmp1=tmp1+(vx[ix+iop-1][iz]-vx[ix-iop][iz])*coef[iop];
			}
			for(iop=1;iop<N+1;iop++){
				tmp2=tmp2+(vz[ix][iz+iop-1]-vz[ix][iz-iop])*coef[iop];
			}
			pxx[ix][iz]=pxx[ix][iz]+(tmp1*(1.0/dx)+tmp2*(1.0/dz))*dt*k[ix][iz];
			pzz[ix][iz]=pzz[ix][iz]+(tmp1*(1.0/dx)+tmp2*(1.0/dz))*dt*k[ix][iz];
		}
	}
}

int main(int argc, char* argv[]){

	int it,ib,sx,sz,ix,iz;
	float tmp;
	float *wlt,*bndr,*coef;
	float **temp,**vel,**rho,**vx,**vz,**pxx,**pzz,**k;
	sf_file Fvel,Frho,Fwfl,Fvx,Fvz;
	
	sf_init(argc,argv);   /* initial madagascar */
	
	Fvel=sf_input("in");  /* velocity model     */
	Frho=sf_input("rho"); /* density            */
	Fvx=sf_output("vx");
	Fvz=sf_output("vz");
	Fwfl=sf_output("out"); /* wavefield snapshot */

	/* input the parameters and data */
	if(!sf_histint(Fvel,"n1",&nz))          sf_error("No nz in velocity model");
	if(!sf_histint(Fvel,"n2",&nx))          sf_error("No nx in velocity model");
	if(!sf_histfloat(Fvel,"d1",&dz))        sf_error("No dz in velocity model");
	if(!sf_histfloat(Fvel,"d2",&dx))        sf_error("No dx in velocity model");
        if(!sf_getint("nb",&nb))                nb=100;
	if(!sf_getint("nt",&nt))                sf_error("No nt found");
	if(!sf_getfloat("dt",&dt))              sf_error("No dt found");
	if(!sf_getfloat("fm",&fm))              fm=20.0;
	if(!sf_getint("N",&N))                N=2;

	/* initial the output data */
	sf_putint(Fwfl,"n1",nz);
	sf_putint(Fwfl,"n2",nx);
	sf_putint(Fwfl,"n3",nt);
	sf_putint(Fvx,"n1",nz);
	sf_putint(Fvx,"n2",nx);
	sf_putint(Fvx,"n3",nt);
	sf_putint(Fvz,"n1",nz);
	sf_putint(Fvz,"n2",nx);
	sf_putint(Fvz,"n3",nt);


	nzpad=nz+2*nb;
	nxpad=nx+2*nb;
	sx=nxpad/2;
	sz=nzpad/2;


	/* allocate variables */
	wlt  =sf_floatalloc(nt);
	bndr =sf_floatalloc(nb);
	temp =sf_floatalloc2(nz,nx);
	vel  =sf_floatalloc2(nzpad,nxpad);
	rho  =sf_floatalloc2(nzpad,nxpad);
	pxx  =sf_floatalloc2(nzpad,nxpad);
	pzz  =sf_floatalloc2(nzpad,nxpad);
	vz   =sf_floatalloc2(nzpad,nxpad);
	vx   =sf_floatalloc2(nzpad,nxpad);
	k  =sf_floatalloc2(nzpad,nxpad);

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
	sf_floatread(temp[0],nz*nx,Fvel);
	expand2d(vel, temp);
	sf_floatread(temp[0],nz*nx,Frho);
	expand2d(rho, temp);
	memset(pxx[0],0,nzpad*nxpad*sizeof(float));
	memset(pzz[0],0,nzpad*nxpad*sizeof(float));
	memset(vx[0],0,nzpad*nxpad*sizeof(float));
	memset(vz[0],0,nzpad*nxpad*sizeof(float));

	/* k=v^2*rho */
	for(ix=0;ix<nxpad;ix++){
		for(iz=0;iz<nzpad;iz++){
		k[ix][iz]=rho[ix][iz]*vel[ix][iz]*vel[ix][iz];
		}
	}

	coef=compute_coef(N);
	/* staggered grid finite difference */
	for(it=0; it<nt; it++)
	{
		window2d(temp,pxx);
		sf_floatwrite(temp[0],nz*nx,Fwfl);
		window2d(temp,vx);
		sf_floatwrite(temp[0],nz*nx,Fvx);
		window2d(temp,vz);
		sf_floatwrite(temp[0],nz*nx,Fvz);
		pxx[sx][sz]+=wlt[it];
		pzz[sx][sz]+=-wlt[it];
		step_forward(pxx,pzz,vz, vx, k, rho,coef);
		apply_sponge(pxx, bndr);
		apply_sponge(pzz, bndr);
		apply_sponge(vx, bndr);
		apply_sponge(vz, bndr);

	}
	
	/* free variables */
	free(wlt);
	free(bndr);
	free(*rho);  free(rho);
	free(*temp); free(temp);
	free(*vel);  free(vel);
	free(*pxx);  free(pxx);
	free(*pzz);  free(pzz);
	free(*vx);   free(vx);
	free(*vz);   free(vz);
	free(*k);    free(k);

	
	exit(0);
}


	
	



































	
