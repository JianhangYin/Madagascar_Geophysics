#include <cmath>

double **getkernelmatrix(double **a,double var,int rows,int columns,int sign){
	double kernel(double var,double *x,double *y,int columns,int sign);
	int i,j,k;
	double *x,*y;
	double **K;
	x=new double[columns];
	y=new double[columns];
	K=new double*[rows];
	for(i=0;i<rows;i++){
		K[i]=new double[rows];
	}
	for(i=0;i<rows;i++){
		for(j=0;j<rows;j++){
			for(k=0;k<columns;k++){
			   x[k]=a[i][k];
			   y[k]=a[j][k];
			}
			K[i][j]=kernel(var,x,y,columns,sign);
		}
	}
	return K;
	delete x;
	delete y;
	for(i=0;i<rows;i++){
		delete[]K[i];
	}
	delete []K;
}

double kernel(double var,double *x,double *y,int columns,int sign){
	double product(double *a,double *b,int size);
	int i;
	double t,s;
	if(sign==1){
	    for(i=0;i<columns;i++){
		    x[i]=x[i]-y[i];
		}
	    t=product(x,x,columns);
	    s=exp(-t/(2*var*var));
	}
	else if(sign==2){
		t=product(x,y,columns);
		s=t;
	}
	return s;
}

double product(double *a,double *b,int size){
	double sum=0;
	for(int i=0;i<size;i++){
		sum+=a[i]*b[i];
	}
	return sum;
}

double **modifykernelmatrix(double **K,int rows){
   int i,j,p,q;
   double s1,s2,s3=0;
   double **KL;
   KL=new double*[rows];
   for(i=0;i<rows;i++){
	  KL[i]=new double[rows];
    }
   for(p=0;p<rows;p++){
	   for(q=0;q<rows;q++){
		   s3+=K[p][q];
	    }
   }
   for(i=0;i<rows;i++){
	   for(p=0,s1=0;p<rows;p++){
		   s1+=K[i][p];
	   }
	   for(j=0;j<rows;j++){
		   for(p=0,s2=0;p<rows;p++){
			   s2+=K[p][j];
		   }
		   KL[i][j]=K[i][j]-s1/rows-s2/rows+s3/(rows*rows);
	   }
   }
   return KL;
   for(i=0;i<rows;i++){
		delete[]KL[i];
   }
   delete []KL;
}

void zhengjiao(double **v,int rows){
	double product(double *a,double *b,int size);
	double **b;
	double *xx,*yy;
	b=new double*[rows];
	for(int i=0;i<rows;i++){
		b[i]=new double[rows];
	}
	xx=new double[rows];
	yy=new double[rows];
	for(int i=0;i<rows;i++){
		b[i][0]=v[i][0];
	}
	for(int j=1;j<rows;j++){
		for(int i=0;i<rows;i++){
			for(int k=0;k<j;k++){
				for(int t=0;t<rows;t++){
					xx[t]=b[t][k];
					yy[t]=v[t][j];
				}
		    	b[i][j]=v[i][j]-(product(xx,yy,rows)/product(xx,xx,rows))*b[i][k];
			}
		}
	}
	for(int i=0;i<rows;i++){
		for(int j=0;j<rows;j++){
			xx[j]=b[j][i];
		}
		yy[i]=sqrt(product(xx,xx,rows));
	}
	for(int i=0;i<rows;i++){
		for(int j=0;j<rows;j++){
			v[i][j]=b[i][j]/yy[j];
		}
	}
	for(int i=0;i<rows;i++){
		delete []b[i];
	}
	delete []b;
	delete[]xx;
	delete[]yy;
}

void selectionsort(double *A,double **v,int rows){
	void swap(double &x,double &y);
	int maxindex;
	int i,j;
	for(int i=0;i<rows-1;i++){
		maxindex=i;
		for(int j=i+1;j<rows;j++){
			if(A[j]>A[maxindex]){
				maxindex=j;
			}
		}
		swap(A[i],A[maxindex]);
		for(int k=0;k<rows;k++){
			swap(v[k][i],v[k][maxindex]);
		}
	}
}

void swap(double &x,double &y){
	double d;
	d=x;
	x=y;
	y=d;
}
