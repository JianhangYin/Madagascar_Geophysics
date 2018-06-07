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

int main(int argc, char* argv[]){
	
	int nt,nm;
	float dt,dm,ot,om;
	float **rnp,**rns,**tf;
	sf_file rnip,rnsf,Tf;
	
	int it,im;
	float tmp1,tmp2;
	
	sf_init(argc,argv);
	
	rnip = sf_input("in");
	rnsf = sf_input("rn");
	Tf   = sf_output("out");
	
	/* input the parameters and data */
	if(!sf_histint(rnip,"n1",&nt))         sf_error("No nt found");	
	if(!sf_histint(rnip,"n2",&nm))         sf_error("No nm found"); 
	if(!sf_histfloat(rnip,"d1",&dt))       sf_error("No dt found");
	if(!sf_histfloat(rnip,"d2",&dm))       sf_error("No dm found");
	if(!sf_histfloat(rnip,"o1",&ot))       sf_error("No ot found");
	if(!sf_histfloat(rnip,"o2",&om))       sf_error("No om found");

	/* initial the output data */

	sf_putint(Tf,"n1",nt);
	sf_putint(Tf,"n2",nm);
	sf_putfloat(Tf,"d1",dt);
	sf_putfloat(Tf,"d2",dm);
	sf_putfloat(Tf,"o1",ot);
	sf_putfloat(Tf,"o2",om);	

	/* allocate variable */
	rnp = sf_floatalloc2(nt,nm);
	rns = sf_floatalloc2(nt,nm);
	tf  = sf_floatalloc2(nt,nm);

	/* initialization */
	sf_floatread(rnp[0],nt*nm,rnip);
	sf_floatread(rns[0],nt*nm,rnsf);
	memset(tf[0],0,nt*nm*sizeof(float));

	tmp1=0;
	tmp2=0;

	/* Tf calculation */
	for(im=0; im<nm; im++){

		for(it=0; it<nt; it++){

			sf_warning("%s %d of %d, %s %d of %d;","cmp",im+1,nm,"time",it+1,nt);
			tmp1=abs(rns[im][it]-rnp[im][it]);
			tmp2=abs(rns[im][it]+rnp[im][it]);
			if(tmp2==0){
				tf[im][it]=0;
			}else{
				tf[im][it]=exp(-tmp1/tmp2);
			}
		}
	}
	sf_floatwrite(tf[0],nt*nm,Tf);
	exit(0);
}






























