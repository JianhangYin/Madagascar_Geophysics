/* transp CRS parameter to a 2D matrix */
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

#include <rsf.h>

static int nt;
static float dt,dc;

int main(int argc, char* argv[]){
	int it,ic,nc;
	float tvt;
        float *ti,*ci;
        float **ccc;
        sf_file tt,cc,crs;
       
        sf_init(argc,argv);   /* initial madagascar */
	
	cc=sf_input("in");    /* crs parameters     */
	tt=sf_input("pick");  /* travel time        */
        crs=sf_output("out"); /* 2D crs figure      */



	/* input the parameters and data */

        if(!sf_histint(cc,"n1",&nc))            sf_error("No nc found");
        if(!sf_getint("nt",&nt))                sf_error("No nt found"); 

        if(!sf_getfloat("dc",&dc))              sf_error("No dc found");       
	if(!sf_getfloat("dt",&dt))              sf_error("No dt found");

	/* initial the output data */

	sf_putint(crs,"n1",nc);
        sf_putint(crs,"n2",nt);
	sf_putfloat(crs,"d1",dc);
	sf_putfloat(crs,"d2",dt);


	
	/* allocate variable */
	ccc =sf_floatalloc2(nc,nt);
	ti  =sf_floatalloc(nc);
	ci  =sf_floatalloc(nc);

	/* initialization */
        sf_floatread(ti,nc,tt);
	sf_floatread(ci,nc,cc);
	memset(ccc[0],0,nt*nc*sizeof(float));

        for(ic=0; ic<nc; ic++)
	{	
		tvt=ti[ic];
		it=(int)(tvt/dt);
		ccc[it][ic]=ci[ic];
	}
	sf_floatwrite(ccc[0],nt*nc,crs);

	exit(0);
}
		
		























