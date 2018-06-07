/* Extract diffraction information in 3D cmp based on travel time picking */
/*
  Copyright (C) 2016 OU Jianhang Yin

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
# include <omp.h>

static int k;
static int tfdif;

int main(int argc, char* argv[]){

	int nt,nm,nh;
	float dt,dm,dh,ot,om,oh;
	float ***dat,***pic,***dif,**tf;

	int it,im,ih,itt;
	
	sf_file data,pick,tfpa,difr;
	sf_init(argc,argv);

	data = sf_input("in");
	pick = sf_input("pick");
	tfpa = sf_input("Tf");
	difr = sf_output("out");

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

	if(!sf_getint("k",&k))                 sf_error("No k found");

	/* initial the output data */
	sf_putint(difr,"n1",nt);
	sf_putint(difr,"n2",nm);
	sf_putint(difr,"n3",nh);
	sf_putfloat(difr,"d1",dt);
	sf_putfloat(difr,"d2",dm);
	sf_putfloat(difr,"d3",dh);
	sf_putfloat(difr,"o1",ot);
	sf_putfloat(difr,"o2",om);
	sf_putfloat(difr,"o3",oh);

	/* allocate variabl */
	dat = sf_floatalloc3(nt,nm,nh);
	pic = sf_floatalloc3(nt,nm,nh);
	dif = sf_floatalloc3(nt,nm,nh);
	tf  = sf_floatalloc2(nt,nm);

	/* initialization */
	sf_floatread(dat[0][0],nt*nm*nh,data);
	sf_floatread(pic[0][0],nt*nm*nh,pick);
	sf_floatread(tf[0],nt*nm,tfpa);
	memset(dif[0][0],0,nt*nm*nh*sizeof(float));

	/* Extract diffraction based on Tf */
	for(im=0; im<nm; im++){

		for(it=0; it<nt; it++){	

			sf_warning("%d threads: %s %d of %d, %s %d of %d;",omp_get_max_threads(),"cmp",im+1,nm,"time",it+1,nt);

			if(tf[im][it]<=1 && tf[im][it]>=tfdif){
				for(ih=0; ih<nh; ih++){
					for(itt=0; itt<nt; itt++){
						if(itt<=(int)(pic[ih][im][it]/dt)+k && itt>=(int)(pic[ih][im][it]/dt)-k){
							dif[ih][im][itt]=dat[ih][im][itt]*tf[im][it]*tf[im][it]*100;
						}
					}
				}
			}
		}
	}
	sf_floatwrite(dif[0][0],nt*nm*nh,difr);
	exit(0);
}

















