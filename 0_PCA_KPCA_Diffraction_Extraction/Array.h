#include <cstring>

unsigned char **New_2D_Chr(long N1,long N2){
	unsigned char **P=new unsigned char *[N1];
	for(long i=0;i<N1;i++){
		P[i]=new unsigned char [N2];
		memset(P[i],0,N2);
	}
	return P;
}

unsigned char ***New_3D_Chr(long N1,long N2,long N3){
	unsigned char ***P=new unsigned char **[N1];
	for(long i=0;i<N1;i++){
		P[i]=New_2D_Chr(N2,N3);
	}
	return P;
}

float **New_2D_Flt(long N1,long N2){
	float **P=new float *[N1];
	for(long i=0;i<N1;i++){
		P[i]=new float [N2];
		memset(P[i],0,N2*4);
	}
	return P;
}

float ***New_3D_Flt(long N1,long N2,long N3){
	float ***P=new float **[N1];
	for(long i=0;i<N1;i++){
		P[i]=New_2D_Flt(N2,N3);
	}
	return P;
}

double **New_2D_Dbl(long N1,long N2){
	double **P=new double *[N1];
	for(long i=0;i<N1;i++){
		P[i]=new double [N2];
		memset(P[i],0,N2*8);
	}
	return P;
}

double ***New_3D_Dbl(long N1,long N2,long N3){
	double ***P=new double **[N1];
	for(long i=0;i<N1;i++){
		P[i]=New_2D_Dbl(N2,N3);
	}
	return P;
}

void Del_2D_Chr(unsigned char **const P,long N1,long N2){
	for(long i=0;i<N1;i++){
		delete [] P[i];
	}
	delete [] P;
}

void Del_3D_Chr(unsigned char ***const P,long N1,long N2,long N3){
	for(long i=0;i<N1;i++){
		Del_2D_Chr(P[i],N2,N3);
	}
	delete [] P;
}

void Del_2D_Flt(float **const P,long N1,long N2){
	for(long i=0;i<N1;i++){
		delete [] P[i];
	}
	delete [] P;
}

void Del_3D_Flt(float ***const P,long N1,long N2,long N3){
	for(long i=0;i<N1;i++){
		Del_2D_Flt(P[i],N2,N3);
	}
	delete [] P;
}

void Del_2D_Dbl(double **const P,long N1,long N2){
	for(long i=0;i<N1;i++){
		delete [] P[i];
	}
	delete [] P;
}

void Del_3D_Dbl(double ***const P,long N1,long N2,long N3){
	for(long i=0;i<N1;i++){
		Del_2D_Dbl(P[i],N2,N3);
	}
	delete [] P;
}

void Cpy_1D_Chr(const unsigned char *const C1,unsigned char *const C2,long N){
	for(long i=0;i<N;i++){
		C2[i]=C1[i];
	}
}

void Cpy_2D_Chr(unsigned char **const C1,unsigned char **const C2,long N1,long N2){
	for(long i=0;i<N1;i++){
		Cpy_1D_Chr(C1[i],C2[i],N2);
	}
}

void Cpy_3D_Chr(unsigned char ***const D1,unsigned char ***const D2,long N1,long N2,long N3){
	for(long i=0;i<N1;i++){
		Cpy_2D_Chr(D1[i],D2[i],N2,N3);
	}
}

void Cpy_1D_Flt(const float *const D1,float *const D2,long N){
	for(long i=0;i<N;i++){
		D2[i]=D1[i];
	}
}

void Cpy_2D_Flt(float **const D1,float **const D2,long N1,long N2){
	for(long i=0;i<N1;i++){
		Cpy_1D_Flt(D1[i],D2[i],N2);
	}
}

void Cpy_3D_Flt(float ***const D1,float ***const D2,long N1,long N2,long N3){
	for(long i=0;i<N1;i++){
		Cpy_2D_Flt(D1[i],D2[i],N2,N3);
	}
}

void Cpy_1D_Dbl(const double *const D1,double *const D2,long N){
	for(long i=0;i<N;i++){
		D2[i]=D1[i];
	}
}

void Cpy_2D_Dbl(double **const D1,double **const D2,long N1,long N2){
	for(long i=0;i<N1;i++){
		Cpy_1D_Dbl(D1[i],D2[i],N2);
	}
}

void Cpy_3D_Dbl(double ***const D1,double ***const D2,long N1,long N2,long N3){
	for(long i=0;i<N1;i++){
		Cpy_2D_Dbl(D1[i],D2[i],N2,N3);
	}
}
