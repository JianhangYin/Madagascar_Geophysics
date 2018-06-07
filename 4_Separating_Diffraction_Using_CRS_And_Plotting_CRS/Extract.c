/* Extract diffraction information in cmp based on travel time picking */
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

static int k;

int main(int argc, char* argv[]){
	int n1,n2,x,y;                                                      /* n1 (half offset) and n2 (time)       */
        float d1,d2;                                                        /* d1 (half offset) and d2 (time)       */
        float *mt,*cc;                                                      /* 1D travel time matrix                */
        float **md,**mk;                                                    /* 2D data matrix and 2D output maxtrix */
        sf_file oridata,tratime,extdata,criteri;                            /* RSF file                             */
       
        sf_init(argc,argv);                                                 /* initial madagascar                   */
	
	oridata = sf_input("in");                                           /* origin data                          */
	tratime = sf_input("pick");                                         /* travel time                          */
        criteri = sf_input("cri");                                          /* criterion                            */
        extdata = sf_output("out");                                         /* extract data                         */


	/* input the parameters and data */

        if(!sf_histint(oridata,"n1",&n1))         sf_error("No n1 found");
        if(!sf_histint(oridata,"n2",&n2))         sf_error("No n2 found"); 

        if(!sf_histfloat(oridata,"d1",&d1))       sf_error("No d1 found");
        if(!sf_histfloat(oridata,"d2",&d2))       sf_error("No d2 found");

	if(!sf_getint("k",&k))                    sf_error("No k found");

	/* initial the output data */

	sf_putint(extdata,"n1",n1);
        sf_putint(extdata,"n2",n2);
	sf_putfloat(extdata,"d1",d1);
	sf_putfloat(extdata,"d2",d2);

	/* allocate variable */
	md = sf_floatalloc2(n1,n2);
	mk = sf_floatalloc2(n1,n2);
	mt = sf_floatalloc(n1);
	cc = sf_floatalloc(1);

	/* initialization */
	sf_floatread(cc,1,criteri);
        sf_floatread(md[0],n1*n2,oridata);
	sf_floatread(mt,n1,tratime);
	memset(mk[0],0,n1*n2*sizeof(float));
	
	if(cc[0]<=1 && cc[0]>=0.7)
	{
        	for(x=0; x<n1; x++)
		{	
			for(y=0; y<n2; y++)
			{
				if(y<=(int)(mt[x]/d2)+k && y>=(int)(mt[x]/d2)-k)
					mk[y][x]=md[y][x];
			}
		}
	}
	sf_floatwrite(mk[0],n1*n2,extdata);

	exit(0);
}
		
		























