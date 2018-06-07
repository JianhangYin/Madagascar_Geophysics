/* Diffraction identification based on Common Reflection Surface */
/* 
  Copyright (C) 2017 OU Jianhang Yin
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
#include <math.h>
# include <omp.h>

static int na;
static int w,ntr;
static float vns;

int main(int argc, char* argv[]){
	
	int nt,nm;
	float dt,dm,ot,om;
	float **dat,**vel,**alp,**rnp,**c;

	int im,it,ia,ic,itr;
	float alphatmp,K1,K2,K1tmp,K2tmp;
	float *coh;
	float ***picktmp;

	sf_file data,velo,alph,rnip,cohe;

	sf_init(argc,argv);

	data = sf_input("in");
        velo = sf_input("vel");
	alph = sf_output("alph");
	cohe = sf_output("coh");
	rnip = sf_output("out");


	/* input the parameters and data */
	if(!sf_histint(data,"n1",&nt))         sf_error("No nt found");	
	if(!sf_histint(data,"n2",&nm))         sf_error("No nm found"); 
	if(!sf_histfloat(data,"d1",&dt))       sf_error("No dt found");
	if(!sf_histfloat(data,"d2",&dm))       sf_error("No dm found");
	if(!sf_histfloat(data,"o1",&ot))       sf_error("No ot found");
	if(!sf_histfloat(data,"o2",&om))       sf_error("No om found");

	if(!sf_getint("na",&na))               sf_error("No na found");
	if(!sf_getint("ntr",&ntr))             sf_error("No ntr found");
	if(!sf_getfloat("vns",&vns))           sf_error("No vns found");

	/* initial the output data */

	sf_putint(alph,"n1",nt);
	sf_putint(alph,"n2",nm);
	sf_putfloat(alph,"d1",dt);
	sf_putfloat(alph,"d2",dm);
	sf_putfloat(alph,"o1",ot);
	sf_putfloat(alph,"o2",om);

	sf_putint(cohe,"n1",nt);
	sf_putint(cohe,"n2",nm);
	sf_putfloat(cohe,"d1",dt);
	sf_putfloat(cohe,"d2",dm);
	sf_putfloat(cohe,"o1",ot);
	sf_putfloat(cohe,"o2",om);

	sf_putint(rnip,"n1",nt);
	sf_putint(rnip,"n2",nm);
	sf_putfloat(rnip,"d1",dt);
	sf_putfloat(rnip,"d2",dm);
	sf_putfloat(rnip,"o1",ot);
	sf_putfloat(rnip,"o2",om);
	
	/* allocate variable */
	dat = sf_floatalloc2(nt,nm);
	vel = sf_floatalloc2(nt,nm);
	alp = sf_floatalloc2(nt,nm);
	rnp = sf_floatalloc2(nt,nm);
	c   = sf_floatalloc2(nt,nm);
	coh = sf_floatalloc(na);

	picktmp = sf_floatalloc3(nt,nm,2*ntr+1);

	/* initialization */
	sf_floatread(dat[0],nt*nm,data);
	sf_floatread(vel[0],nt*nm,velo);
	memset(alp[0],0,nt*nm*sizeof(float));
	memset(rnp[0],0,nt*nm*sizeof(float));
	memset(c[0],0,nt*nm*sizeof(float));
	memset(coh,0,na*sizeof(float));
	memset(picktmp[0][0],0,nt*nm*(2*ntr+1)*sizeof(float));

	/* CRS core part 2 */
	for(im=0; im<nm; im++){

		for(it=0; it<nt; it++){

			sf_warning("%d threads: %s %d of %d, %s %d of %d;",omp_get_max_threads(),"cmp",im+1,nm,"time",it+1,nt);

			/* Linear ZO stack */
			for(ia=1; ia<na+1; ia++){
				alphatmp=-SF_PI/2+ia*SF_PI/(na+1);
				/* Travel time picking */
				for(itr=0; itr<2*ntr+1; itr++){
					picktmp[itr][im][it]=it*dt+2*sin(alphatmp)*(itr-ntr)*dm/vns;
				}
				K1=0;
				K2=0;
				/* Coherent Analysis */
				for(ic=-w/2; ic<w/2+1; ic++){
					K1tmp=0;
					K2tmp=0;
					for(itr=0; itr<2*ntr+1; itr++){
						if((int)(picktmp[itr][im][it]/dt)+ic<nt && (int)(picktmp[itr][im][it]/dt)+ic>=0 && im-ntr+itr>=0 && im-ntr+itr<nm){
							K1tmp+=dat[im-ntr+itr][(int)(picktmp[itr][im][it]/dt)+ic];
							K2tmp+=dat[im-ntr+itr][(int)(picktmp[itr][im][it]/dt)+ic]*dat[im-ntr+itr][(int)(picktmp[itr][im][it]/dt)+ic];
						}
					}
					K1+=K1tmp*K1tmp;
					K2+=K2tmp;
				}	
				/* Choosing the maxiam value of coherence */	
				if(K2==0){
					coh[ia]=0;
				}else{		
					coh[ia]=K1/(K2*(2*ntr+1));
				}
				if(ia==1 || c[im][it]<coh[ia]){
					alp[im][it]=alphatmp;
					c[im][it]=coh[ia];
				}
			}
			rnp[im][it]=vel[im][it]*vel[im][it]*it*dt*cos(alp[im][it])*cos(alp[im][it])/(2*vns);

		}
	}
	
	sf_floatwrite(rnp[0],nt*nm,rnip);
	sf_floatwrite(alp[0],nt*nm,alph);
	sf_floatwrite(c[0],nt*nm,cohe);

	exit(0);
	
	
}

