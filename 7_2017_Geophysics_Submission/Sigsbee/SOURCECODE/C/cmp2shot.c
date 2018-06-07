/* Tranfer CMP (Midpoint-halfoffset-time) to shotgather */
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
# include <string.h>

int main(int argc, char* argv[]){

	int ns,nr,nm,nh,nt;
	float ds,dr,dm,dh,dt;
	float os,or,om,oh,ot;
	char *trace, *zero;
	off_t pos, step;
	
	int is,ir;

	sf_file cmpd,shot;
	sf_init(argc,argv);

	cmpd = sf_input("in");
	shot = sf_output("out");

	/* input the parameters and data */
	if(!sf_histint(cmpd,"n1",&nt))         sf_error("No nt found");	
	if(!sf_histint(cmpd,"n2",&nm))         sf_error("No nm found"); 
	if(!sf_histint(cmpd,"n3",&nh))         sf_error("No nh found");
	if(!sf_histfloat(cmpd,"d1",&dt))       sf_error("No dt found");
	if(!sf_histfloat(cmpd,"d2",&dm))       sf_error("No dm found");
	if(!sf_histfloat(cmpd,"d3",&dh))       sf_error("No dh found");
	if(!sf_histfloat(cmpd,"o1",&ot))       sf_error("No ot found");
	if(!sf_histfloat(cmpd,"o2",&om))       sf_error("No om found");
	if(!sf_histfloat(cmpd,"o3",&oh))       sf_error("No oh found");

	/* initial the output data */
	nr=nm+2*(nh-1);
	ns=nm+2*(nh-1);

	dr=dm;
	ds=dh;

	or=om-(nh-1)*dh;
	os=om-(nh-1)*dh;

	sf_warning("nr=%d",nr);
	sf_warning("ns=%d",ns);

	sf_warning("nm=%d",nm);
	sf_warning("nh=%d",nh);

	sf_putint(shot,"n1",nt);
	sf_putint(shot,"n2",nr);
	sf_putint(shot,"n3",ns);
	sf_putfloat(shot,"d1",dt);
	sf_putfloat(shot,"d2",dr);
	sf_putfloat(shot,"d3",ds);
	sf_putfloat(shot,"o1",ot);
	sf_putfloat(shot,"o2",or);
	sf_putfloat(shot,"o3",os);

	nt   *= sf_esize(cmpd);
	trace = sf_charalloc(nt);
	zero  = sf_charalloc(nt);
	memset(zero,0,nt);

	sf_fileflush(shot,cmpd);
	sf_setform(cmpd,SF_NATIVE);
	sf_setform(shot,SF_NATIVE);
	sf_unpipe(cmpd,(off_t) nh*nm*nt);
	pos   = sf_tell(cmpd);

	/* Main part */
	for(is=0; is<ns; is++){
		for(ir=0; ir<nr; ir++){
			sf_warning("%s %d of %d,%s %d of %d;","Receiver",ir+1,nr,"Source",is+1,ns);
			if(((is+ir)/2-nh)>=0 && ((is+ir)/2-nh)<nm && abs(is-ir)/2>=0 && abs(is-ir)/2<nh){
				step = (off_t)(abs(is-ir)/2*nm+(is+ir)/2-nh);
				step*= nt;
				sf_seek(cmpd,pos+step,SEEK_SET);
				sf_charread(trace,nt,cmpd);
				sf_charwrite(trace,nt,shot);
			}else{
				sf_charwrite(zero,nt,shot);
			}
		}
	}

	exit(0);
}	






















	
