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

static float vns;

int main(int argc, char* argv[]){
	
	int nt,nm,nh;
	float dt,dm,dh,ot,om,oh;
	float ***dat,***pic,**alp,**rnp;	

	int im,it,ih;
	float tmp1,tmp2;
		
	
	sf_file data,alph,rnip,pick;	

	sf_init(argc,argv);

	data = sf_input("in");
	alph = sf_input("alph");
	rnip = sf_input("rnip");
	pick = sf_output("out");

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

	if(!sf_getfloat("vns",&vns))           sf_error("No vns found");

	/* initial the output data */
	sf_putint(pick,"n1",nt);
	sf_putint(pick,"n2",nm);
	sf_putint(pick,"n3",nh);
	sf_putfloat(pick,"d1",dt);
	sf_putfloat(pick,"d2",dm);
	sf_putfloat(pick,"d3",dh);
	sf_putfloat(pick,"o1",ot);
	sf_putfloat(pick,"o2",om);
	sf_putfloat(pick,"o3",oh);

	/* allocate variable */
	dat = sf_floatalloc3(nt,nm,nh);
	pic = sf_floatalloc3(nt,nm,nh);
	alp = sf_floatalloc2(nt,nm);
	rnp = sf_floatalloc2(nt,nm);

	/* initialization */
	sf_floatread(dat[0][0],nt*nm*nh,data);
	sf_floatread(alp[0],nt*nm,alph);
	sf_floatread(rnp[0],nt*nm,rnip);
	memset(pic[0][0],0,nt*nm*nh*sizeof(float));

	/* CRS core part 4 */
	for(im=0; im<nm; im++){
		sf_warning("%s %d of %d;","cmp",im+1,nm);
		for(it=0; it<nt; it++){
			tmp1=0;
			tmp2=0;
			if(rnp[im][it]!=0){
				for(ih=0; ih<nh; ih++){
					tmp1=it*dt;
					tmp2=2*it*dt*cos(alp[im][it])*cos(alp[im][it])*ih*dh*ih*dh/(vns*rnp[im][it]);
					pic[ih][im][it]=sqrt(tmp1*tmp1+tmp2);
				}
			}

		}
	}
	sf_floatwrite(pic[0][0],nt*nm*nh,pick);

	exit(0);
}

























