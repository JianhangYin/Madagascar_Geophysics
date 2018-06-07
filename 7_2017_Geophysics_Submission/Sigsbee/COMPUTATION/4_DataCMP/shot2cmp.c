/* Tranfer shotgather to CMP (Midpoint-halfoffset-time) */
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
# include <string.h>

static int hof;

int main(int argc, char* argv[]){

	int ns,nr,nm,nh,nt;
	float ds,dr,dm,dh,dt;
	float os,or,om,oh,ot;
	char *trace;
	off_t pos, step;
	
	int im,ih,it;

	sf_file cmpd,shot;
	sf_init(argc,argv);

	shot = sf_input("in");
	cmpd = sf_output("out");

	/* input the parameters and data */
	if(!sf_histint(shot,"n1",&nt))         sf_error("No nt found");	
	if(!sf_histint(shot,"n2",&nr))         sf_error("No nr found"); 
	if(!sf_histint(shot,"n3",&ns))         sf_error("No ns found");
	if(!sf_histfloat(shot,"d1",&dt))       sf_error("No dt found");
	if(!sf_histfloat(shot,"d2",&dr))       sf_error("No dr found");
	if(!sf_histfloat(shot,"d3",&ds))       sf_error("No ds found");
	if(!sf_histfloat(shot,"o1",&ot))       sf_error("No ot found");
	if(!sf_histfloat(shot,"o2",&or))       sf_error("No or found");
	if(!sf_histfloat(shot,"o3",&os))       sf_error("No os found");

	if(!sf_getint("hof",&hof))             sf_error("No hof found");

	/* initial the output data */

	nm=nr-2*hof;
	nh=hof+1;
	dm=dr;
	dh=ds;
	om=or+hof*dm;
	oh=0;

	sf_warning("nr=%d",nr);
	sf_warning("ns=%d",ns);

	sf_warning("nm=%d",nm);
	sf_warning("nh=%d",nh);

	sf_putint(cmpd,"n1",nt);
	sf_putint(cmpd,"n2",nm);
	sf_putint(cmpd,"n3",nh);
	sf_putfloat(cmpd,"d1",dt);
	sf_putfloat(cmpd,"d2",dm);
	sf_putfloat(cmpd,"d3",dh);
	sf_putfloat(cmpd,"o1",ot);
	sf_putfloat(cmpd,"o2",om);
	sf_putfloat(cmpd,"o3",oh);

	nt   *= sf_esize(shot);
	trace = sf_charalloc(nt);
	sf_fileflush(cmpd,shot);
	sf_setform(shot,SF_NATIVE);
	sf_setform(cmpd,SF_NATIVE);
	sf_unpipe(shot,(off_t) ns*nr*nt);
	pos   = sf_tell(shot);

	/* Main part */
	for(ih=0; ih<nh; ih++){
		for(im=0; im<nm; im++){
			sf_warning("%s %d of %d,%s %d of %d;","half_Offset",ih+1,nh,"Midpoint",im+1,nm);
			step = (off_t)(hof+im+ih)*nr+(hof+im-ih);
			step*= nt;
			sf_seek(shot,pos+step,SEEK_SET);
			sf_charread(trace,nt,shot);
			sf_charwrite(trace,nt,cmpd);
		}
	}

	exit(0);
}	






















	
