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

# include <rsf.h>
# include <math.h>

static int nv;
static float dv;
static float v0;
static int w;

int main(int argc, char* argv[]){
	
	int nt,nm,nh;
	float dt,dm,dh,ot,om,oh;
	float ***dat,***sca;

	int im,it,ih,iv,ic;
	float vtmp,K1,K2,K1tmp,K2tmp;
	float ***picktmp;

	sf_file data,scan;

	sf_init(argc,argv);

	data = sf_input("in");
	scan = sf_output("out");


	/* input the parameters and data */
	if(!sf_histint(data,"n1",&nt))         sf_error("No nt found");	
	if(!sf_histint(data,"n2",&nm))         sf_error("No nm found"); 
	if(!sf_histint(data,"n3",&nh))         sf_error("No nh found");
	if(!sf_histfloat(data,"d1",&dt))       sf_error("No dt found");
	if(!sf_histfloat(data,"d2",&dm))       sf_error("No dm found");
	if(!sf_histfloat(data,"d3",&dh))       sf_error("No dh found");
	if(!sf_histfloat(data,"o1",&ot))       sf_error("No ot found");
	if(!sf_histfloat(data,"o2",&om))       sf_error("No om found");
	if(!sf_histfloat(data,"o3",&oh))       sf_error("No oh found");

	if(!sf_getint("nv",&nv))               sf_error("No nv found");

	if(!sf_getfloat("dv",&dv))             sf_error("No dv found");

	if(!sf_getfloat("v0",&v0))             sf_error("No v0 found");

	if(!sf_getint("w",&w))                 sf_error("No w found");

	/* initial the output data */
	sf_putint(scan,"n1",nt);
	sf_putint(scan,"n2",nm);
	sf_putint(scan,"n3",nv);
	sf_putfloat(scan,"d1",dt);
	sf_putfloat(scan,"d2",dm);
	sf_putfloat(scan,"d3",dv);
	sf_putfloat(scan,"o1",ot);
	sf_putfloat(scan,"o2",om);
	sf_putfloat(scan,"o3",v0);
	
	/* allocate variabl */
	dat = sf_floatalloc3(nt,nm,nh);
	sca = sf_floatalloc3(nt,nm,nv);
	picktmp = sf_floatalloc3(nt,nm,nh);

	/* initialization */
	sf_floatread(dat[0][0],nt*nm*nh,data);
	memset(sca[0][0],0,nt*nm*nv*sizeof(float));
	memset(picktmp[0][0],0,nt*nm*nh*sizeof(float));
		
	/* CRS core part 1 */
	for(im=0; im<nm; im++){
		sf_warning("%s %d of %d;","cmp",im+1,nm);

		for(it=0; it<nt; it++){

			/* Automatic CMP stack */
			for(iv=0; iv<nv; iv++){
				vtmp=v0+iv*dv;
				/* Travel time Picking (it) */
				for(ih=0; ih<nh; ih++){
					picktmp[ih][im][it]=sqrt(it*dt*it*dt+4*ih*dh*ih*dh/(vtmp*vtmp));
				}
				K1=0;
				K2=0;
				/* Coherent Analysis */
				for(ic=-w/2; ic<w/2+1; ic++){
					K1tmp=0;
					K2tmp=0;
					for(ih=0; ih<nh; ih++){
						if((int)(picktmp[ih][im][it]/dt)+ic<nt && (int)(picktmp[ih][im][it]/dt)+ic>=0){
							K1tmp+=dat[ih][im][(int)(picktmp[ih][im][it]/dt)+ic];
							K2tmp+=dat[ih][im][(int)(picktmp[ih][im][it]/dt)+ic]*dat[ih][im][(int)(picktmp[ih][im][it]/dt)+ic];
						}
					}
					K1+=K1tmp*K1tmp;
					K2+=K2tmp;
				}
				if(K2==0){
					sca[iv][im][it]=0;
				}else{
					sca[iv][im][it]=K1/(K2*nh);
				}
			}

		}
	}
	sf_floatwrite(sca[0][0],nt*nm*nv,scan);

	exit(0);
	
	
}

