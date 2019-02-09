#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <Eigen/Dense>
#include <Segy.h>
#include <Array.h>

// diffraction extraction on post-stack seismic data
// kernel principal component analysis(KPCA)
// Author @Jianhang Yin 2016
int main(){
	using namespace std;
	using namespace Eigen;
	int iRangeLin[2],iRangeTrc[2],iSizeSlipWin,iLag[3],iStepDis,iSizeSeis,iAsWin,iPCAno,sign;
	string cFILEname,cFILEnamepca,cFILETemp;
	ifstream fileidpara;
	// read the parameters from file
	fileidpara.open("parameter.txt");
	getline(fileidpara,cFILEname);
	getline(fileidpara,cFILEname);
	getline(fileidpara,cFILEnamepca);
	getline(fileidpara,cFILEnamepca);
	getline(fileidpara,cFILETemp);
	fileidpara>>iRangeLin[0];
	fileidpara>>iRangeLin[1];
	getline(fileidpara,cFILETemp);
	getline(fileidpara,cFILETemp);
	fileidpara>>iRangeTrc[0];
	fileidpara>>iRangeTrc[1];
	getline(fileidpara,cFILETemp);
	getline(fileidpara,cFILETemp);
	fileidpara>>iSizeSlipWin;
	getline(fileidpara,cFILETemp);
	getline(fileidpara,cFILETemp);
	fileidpara>>iLag[0];
	fileidpara>>iLag[1];
	fileidpara>>iLag[2];
	getline(fileidpara,cFILETemp);
	getline(fileidpara,cFILETemp);
	fileidpara>>iStepDis;
	getline(fileidpara,cFILETemp);
	getline(fileidpara,cFILETemp);
	fileidpara>>iSizeSeis;
	getline(fileidpara,cFILETemp);
	getline(fileidpara,cFILETemp);
	fileidpara>>iAsWin;
	getline(fileidpara,cFILETemp);
	getline(fileidpara,cFILETemp);
	fileidpara>>iPCAno;
	fileidpara.close();

	string cFILEnamepcaD[iPCAno];
	ifstream fileidsegy;
	ofstream fileidsegypca[iPCAno];
	stringstream ss;
	fileidsegy.open(cFILEname.c_str(),ios::binary);
	for(int nNum=0;nNum<iPCAno;nNum++){
		ss<<nNum+1;
		cFILEnamepcaD[nNum]=cFILEnamepca+ss.str();
		fileidsegypca[nNum].open(cFILEnamepcaD[nNum].c_str(),ios::binary);
		ss.str("");
	}

	unsigned char *fSeisHd=new unsigned char [3600];
	fileidsegy.read((char *)fSeisHd,3600);
	for(int nNum=0;nNum<iPCAno;nNum++){
		fileidsegypca[nNum].write((char *)fSeisHd,3600);
	}

	int iHalfLinRange=iLag[0]+iSizeSlipWin;
	int iHalfTrcRange=iLag[1]+iSizeSlipWin;
	int iHalfTimRange=iLag[2]+iSizeSlipWin;

	// End of Part A

	// Start of Part B

	int iLinCal[2];
	int nLin,nTempLinCount,nTempLinCountHd,nTempTrcCount,nTempTrcHd,nTempTimeCount;
	unsigned char ***fTrcHd=New_3D_Chr(iHalfLinRange+1,iRangeTrc[1]-iRangeTrc[0]+1,240);
	unsigned char *fTempTrc=new unsigned char [iSizeSeis*sizeof(float)];
	float *fTrc=new float [iSizeSeis];
	// Step 1: extend the seismic data to deal with the boundary problem
	double ***fCalLin=New_3D_Dbl(1+2*iHalfLinRange,iRangeTrc[1]-iRangeTrc[0]+1+2*iHalfTrcRange,iSizeSeis+2*iHalfTimRange);
	for(nLin=iRangeLin[0];nLin<=iRangeLin[0];nLin++){
		iLinCal[0]=nLin-iHalfLinRange;
		iLinCal[1]=nLin+iHalfLinRange;
		nTempLinCount=0;
		nTempLinCountHd=-1;
		for(int nRdLin=iLinCal[0];nRdLin<=iLinCal[1];nRdLin++){
			nTempTrcCount=0;
			nTempTrcHd=0;
			if(nRdLin<iRangeLin[0]||nRdLin>iRangeLin[1]){
				for(int nTempTrc=iRangeTrc[0]-iHalfTrcRange;nTempTrc<=iRangeTrc[1]+iHalfTrcRange;nTempTrc++){
					for(int nTempTime=0;nTempTime<iSizeSeis+2*iHalfTimRange;nTempTime++){
						fCalLin[nTempLinCount][nTempTrcCount][nTempTime]=-9999.9;
					}
					nTempTrcCount=nTempTrcCount+1;
				}
			}else{
				for(int nTempTrc=iRangeTrc[0]-iHalfTrcRange;nTempTrc<=iRangeTrc[1]+iHalfTrcRange;nTempTrc++){
					if(nTempTrc<iRangeTrc[0]||nTempTrc>iRangeTrc[1]){
						for(int nTempTime=0;nTempTime<iSizeSeis+2*iHalfTimRange;nTempTime++){
							fCalLin[nTempLinCount][nTempTrcCount][nTempTime]=-9999.9;
						}
						nTempTrcCount=nTempTrcCount+1;
					}else{
						if(nTempTrcHd==0){
							nTempLinCountHd=nTempLinCountHd+1;
						}
						fileidsegy.read((char *)fTrcHd[nTempLinCountHd][nTempTrcHd],240);
						fileidsegy.read((char *)fTempTrc,iSizeSeis*sizeof(float));
						Bin_Flt_IBM_BE(fTempTrc,fTrc,iSizeSeis);
						nTempTrcHd=nTempTrcHd+1;
						nTempTimeCount=0;
						for(int nTempTime=-iHalfTimRange;nTempTime<iSizeSeis+iHalfTimRange;nTempTime++){
							if(nTempTime<0){
								fCalLin[nTempLinCount][nTempTrcCount][nTempTimeCount]=fTrc[0];
								nTempTimeCount=nTempTimeCount+1;
							}else if(nTempTime>=iSizeSeis){
								fCalLin[nTempLinCount][nTempTrcCount][nTempTimeCount]=fTrc[iSizeSeis-1];
								nTempTimeCount=nTempTimeCount+1;
							}else{
								fCalLin[nTempLinCount][nTempTrcCount][nTempTimeCount]=fTrc[nTempTimeCount-iHalfTimRange];
								nTempTimeCount=nTempTimeCount+1;
							}
						}
						nTempTrcCount=nTempTrcCount+1;
					}
				}
			}
			nTempLinCount=nTempLinCount+1;
		}
	}
	// Step 2: assign values to the matrix boundary
    nTempLinCount=iHalfLinRange;
	nTempTrcCount=0;
    for(int nTempTrc=iRangeTrc[0]-iHalfTrcRange;nTempTrc<=iRangeTrc[1]+iHalfTrcRange;nTempTrc++){
		if(nTempTrc<iRangeTrc[0]){
			Cpy_1D_Dbl(fCalLin[nTempLinCount][iHalfTrcRange],fCalLin[nTempLinCount][nTempTrcCount],iSizeSeis+2*iHalfTimRange);
			nTempTrcCount=nTempTrcCount+1;
		}else if(nTempTrc>iRangeTrc[1]){
			Cpy_1D_Dbl(fCalLin[nTempLinCount][nTempTrcCount-1],fCalLin[nTempLinCount][nTempTrcCount],iSizeSeis+2*iHalfTimRange);
			nTempTrcCount=nTempTrcCount+1;
		}else{
			nTempTrcCount=nTempTrcCount+1;
		}
	}
	nTempLinCount=iHalfLinRange*2;
	for(int nRdLin=iLinCal[1];nRdLin>=iLinCal[0];nRdLin--){
		if(nRdLin!=iRangeLin[0]){
			Cpy_2D_Dbl(fCalLin[iHalfLinRange],fCalLin[nTempLinCount],iRangeTrc[1]-iRangeTrc[0]+1+2*iHalfTrcRange,iSizeSeis+2*iHalfTimRange);
		}
		nTempLinCount=nTempLinCount-1;
	}
	double max=fCalLin[0][0][0];
	double min=fCalLin[0][0][0];
	for(int k=0;k<1+2*iHalfLinRange;k++){
		for(int p=0;p<iRangeTrc[1]-iRangeTrc[0]+1+2*iHalfTrcRange;p++){
			for(int q=0;q<iSizeSeis+2*iHalfTimRange;q++){
				if( max<fCalLin[k][p][q]){
					max=fCalLin[k][p][q];
				}else if(min>fCalLin[k][p][q]){
					min=fCalLin[k][p][q];
				}
			}
		}
	}
	for(int k=0;k<1+2*iHalfLinRange;k++){
		for(int p=0;p<iRangeTrc[1]-iRangeTrc[0]+1+2*iHalfTrcRange;p++){
			for(int q=0;q<iSizeSeis+2*iHalfTimRange;q++){
				fCalLin[k][p][q]=10*(fCalLin[k][p][q]-min)/(max-min);
			}
		}
	}
	cout<<"the processing line NO. is : "<<nLin-1<<endl;
	cout<<"choose the kernel function: 1. Guess 2. Linear"<<endl;
	cin>>sign;
	int rbf_var;
    if(sign==1){
		cout<<"choose the rbf_var"<<endl;
		cin>>rbf_var;
	}

	// End of Part B

	// Start of Part C

	nTempTrcCount=0;
	int nNum,nTemp1,nTemp2,nTemp3;
	float ***fPCAcubeResult=New_3D_Flt(iPCAno,iRangeTrc[1]-iRangeTrc[0]+1,iSizeSeis);
	for(int nTempTrc=iRangeTrc[0];nTempTrc<=iRangeTrc[1];nTempTrc+=iAsWin){
		for(int nTempTim=0;nTempTim<iSizeSeis;nTempTim+=iAsWin){
			nNum=0;
			// d = the number of the small cubes contained in the large cube as the attributes
			// M = the number of the points contained in the small cube as sample numbers
			int rows=(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1);
			int columns=(2*(iLag[0]/iStepDis)+1)*(2*(iLag[1]/iStepDis)+1)*(2*(iLag[2]/iStepDis)+1);
			double **fCalCube=New_2D_Dbl((2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1),(2*(iLag[0]/iStepDis)+1)*(2*(iLag[1]/iStepDis)+1)*(2*(iLag[2]/iStepDis)+1));
			double **K=New_2D_Dbl((2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1),(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1));
			double **KL=New_2D_Dbl((2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1),(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1));
			for(int nLoopLin=0;nLoopLin<2*(iLag[0]/iStepDis)+1;nLoopLin++){
				for(int nLoopTrc=0;nLoopTrc<2*(iLag[1]/iStepDis)+1;nLoopTrc++){
					for(int nLoopTim=0;nLoopTim<2*(iLag[2]/iStepDis)+1;nLoopTim++){
						nNum=nNum+1;
						for(int nCalLin=0;nCalLin<2*iSizeSlipWin+1;nCalLin++){
							for(int nCalTrc=0;nCalTrc<2*iSizeSlipWin+1;nCalTrc++){
								for(int nCalTim=0;nCalTim<2*iSizeSlipWin+1;nCalTim++){
									nTemp1=nCalTim+nLoopTim*iStepDis+nTempTim;
									nTemp2=nCalTrc+nLoopTrc*iStepDis+nTempTrcCount;
									nTemp3=nCalLin+nLoopLin*iStepDis;
									fCalCube[nCalLin*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)+nCalTim+nCalTrc*(2*iSizeSlipWin+1)][nNum-1]=fCalLin[nTemp3][nTemp2][nTemp1];
								}
							}
						}
					}
				}
			}
			// calculate the kernel matrix k
			for(int i=0;i<rows;i++){
				for(int j=0;j<rows;j++){
					double temp;
					double *Ktemp1=new double [columns];
			        double *Ktemp2=new double [columns];
					temp=0;
					for(int k=0;k<columns;k++){
						Ktemp1[k]=fCalCube[i][k];
						Ktemp2[k]=fCalCube[j][k];
					}
					if(sign==1){
						for(int n=0;n<columns;n++){
							Ktemp1[n]=Ktemp1[n]-Ktemp2[n];
						}
						for(int m=0;m<columns;m++){
							temp+=Ktemp1[m]*Ktemp1[m];
						}
						K[i][j]=exp(-temp/(2*rbf_var*rbf_var));
					}
					else if(sign==2){
						for(int m=0;m<columns;m++){
							temp+=Ktemp1[m]*Ktemp2[m];
						}
						K[i][j]=temp;
					}
					delete Ktemp1;
			        delete Ktemp2;
				}
			}
			// calculate the revised kernel matrix KL
			double s1,s2,s3;
			s3=0;
			for(int p=0;p<rows;p++){
				for(int q=0;q<rows;q++){
					s3+=K[p][q];
				}
			}
			for(int i=0;i<rows;i++){
				s1=0;
				for(int p=0;p<rows;p++){
					s1+=K[i][p];
				}
				for(int j=0;j<rows;j++){
					s2=0;
					for(int p=0;p<rows;p++){
						s2+=K[p][j];
					}
					KL[i][j]=K[i][j]-s1/rows-s2/rows+s3/(rows*rows);
				}
			}
			MatrixXd fK_n(rows,rows);
			for(int i=0;i<rows;i++){
				for(int j=0;j<rows;j++){
					fK_n(j,i)=KL[i][j];
				}
			}
			// calculate the eigenvalue and eigenvector of KL
			EigenSolver<MatrixXd> es(fK_n);
			MatrixXcd Vc=es.eigenvectors();
			VectorXcd Dc=es.eigenvalues();
			double **V=New_2D_Dbl(rows,rows);
			double *D=new double [rows];
			for(int i=0;i<rows;i++){
				for(int j=0;j<rows;j++){
					V[i][j]=Vc(j,i).real();
				}
			}
			for(int i=0;i<rows;i++){
				D[i]=Dc[i].real();
			}
			// double **B=New_2D_Dbl(rows,rows);
			// double *XX=new double [rows];
			// double *YY=new double [rows];
			// for(int i=0;i<rows;i++){
				// B[i][0]=V[i][0];
			// }
			// for(int j=1;j<rows;j++){
				// for(int i=0;i<rows;i++){
					// for(int k=0;k<j;k++){
						// for(int t=0;t<rows;t++){
							// XX[t]=B[t][k];
							// YY[t]=V[t][j];
						// }
						// double temp1=0;
                        // for(int m=0;m<rows;m++){
							// temp1+=XX[m]*YY[m];
						// }
						// double temp2=0;
                        // for(int n=0;n<rows;n++){
							// temp2+=XX[n]*XX[n];
						// }
						// B[i][j]=V[i][j]-(temp1/temp2)*B[i][k];
					// }
				// }
			// }
			// for(int i=0;i<rows;i++){
				// for(int j=0;j<rows;j++){
					// XX[j]=B[j][i];
				// }
				// double temp3=0;
				// for(int n=0;n<rows;n++){
					// temp3+=XX[n]*XX[n];
				// }
				// YY[i]=sqrt(temp3);
			// }
	        // for(int i=0;i<rows;i++){
				// for(int j=0;j<rows;j++){
					// V[i][j]=B[i][j]/YY[j];
				// }
			// }
			// Del_2D_Dbl(B,rows,rows);
			// delete XX;
			// delete YY;
			// sort the eigenvalue and eigenvector and select the first 5
			double *Vtemp=new double [rows];
			double Dtemp;
			for(int nNumy=0;nNumy<rows-1;nNumy++){
				for(int nNumx=nNumy;nNumx<rows-1;nNumx++){
					if(D[nNumx+1]>D[nNumy]){
						Cpy_1D_Dbl(V[nNumx+1],Vtemp,rows);
						Cpy_1D_Dbl(V[nNumy],V[nNumx+1],rows);
						Cpy_1D_Dbl(Vtemp,V[nNumy],rows);
						Dtemp=D[nNumx+1];
						D[nNumx+1]=D[nNumy];
						D[nNumy]=Dtemp;
					}
				}
			}
			// map the sorted eigenvector to the revised kernel function
			double **fEigVet=New_2D_Dbl(iPCAno,rows);
			for(int nNumx=0;nNumx<iPCAno;nNumx++)
				Cpy_1D_Dbl(V[nNumx],fEigVet[nNumx],rows);
			double *fDiffCube=new double [rows];
			double ***fPCAcubeWhole=New_3D_Dbl(2*iSizeSlipWin+1,(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1),iPCAno);
			for(int nNumy=0;nNumy<2*iSizeSlipWin+1;nNumy++){
				for(int nNumz=0;nNumz<(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1);nNumz++){
					for(int nNumx=0;nNumx<rows;nNumx++){
						fDiffCube[nNumx]=KL[nNumx][nNumy*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)+nNumz];
					}
					for(int i=0;i<iPCAno;i++){
						for(int j=0;j<rows;j++){
							fPCAcubeWhole[nNumy][nNumz][i]+=fDiffCube[j]*fEigVet[i][j];
						}
					}
				}
			}
			double ***fPCAcube=New_3D_Dbl(iPCAno,2*iSizeSlipWin+1,(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1));
			for(int nNumy=0;nNumy<iPCAno;nNumy++){
				for(int nNumz=0;nNumz<2*iSizeSlipWin+1;nNumz++){
					for(int nNumx=0;nNumx<(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1);nNumx++){
						fPCAcube[nNumy][nNumz][nNumx]=fPCAcubeWhole[nNumz][nNumx][nNumy];
					}
				}
			}
		    // reconstruct the data
			double fRightDis=sqrt(iAsWin*iAsWin+iAsWin*iAsWin),fRight;
			int nPosAsRes[2],nPosAs;
			for(int nNumx=0;nNumx<iPCAno;nNumx++){
				for(int nNumz=-iAsWin;nNumz<=iAsWin;nNumz++){
					nPosAsRes[1]=nTempTrcCount+nNumz;
					for(int nNumT=-iAsWin;nNumT<=iAsWin;nNumT++){
						nPosAs=(iSizeSlipWin+nNumz)*(2*iSizeSlipWin+1)+(iSizeSlipWin+nNumT);
						nPosAsRes[0]=nTempTim+nNumT;
						if(nPosAsRes[0]>=0&&nPosAsRes[0]<iSizeSeis&&nPosAsRes[1]>=0&&nPosAsRes[1]<iRangeTrc[1]-iRangeTrc[0]+1){
							fRight=(fRightDis-sqrt(nNumz*nNumz+nNumT*nNumT))/fRightDis;
							fPCAcubeResult[nNumx][nPosAsRes[1]][nPosAsRes[0]]+=fPCAcube[nNumx][iSizeSlipWin+1][nPosAs]*fRight;
						}
					}
				}
			}
			if(iSizeSeis-nTempTim<=iAsWin){
				for(int nNumz=-iAsWin;nNumz<=iAsWin;nNumz++){
					nPosAsRes[1]=nTempTrcCount+nNumz;
					for(int nNumx=0;nNumx<iPCAno;nNumx++){
						for(int nNumT=0;nNumT<=iSizeSeis-nTempTim;nNumT++){
							if(nPosAsRes[1]>=0&&nPosAsRes[1]<iRangeTrc[1]-iRangeTrc[0]+1){
								fPCAcubeResult[nNumx][nPosAsRes[1]][nTempTim+nNumT]+=fPCAcubeResult[nNumx][nPosAsRes[1]][nTempTim];
							}
						}
					}
				}
			}
			Del_2D_Dbl(fCalCube,(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1),(2*(iLag[0]/iStepDis)+1)*(2*(iLag[1]/iStepDis)+1)*(2*(iLag[2]/iStepDis)+1));
			Del_2D_Dbl(K,(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1),(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1));
			Del_2D_Dbl(KL,(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1),(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1));
			Del_2D_Dbl(V,rows,rows);
			delete D;
			delete Vtemp;
			Del_2D_Dbl(fEigVet,iPCAno,rows);
			delete fDiffCube;
			Del_3D_Dbl(fPCAcubeWhole,2*iSizeSlipWin+1,(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1),iPCAno);
			Del_3D_Dbl(fPCAcube,iPCAno,2*iSizeSlipWin+1,(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1));
		}
		cout<<"the trace NO. is : "<<nTempTrc<<endl;
		if(iRangeTrc[1]-nTempTrc<iAsWin){
			for(int nNumz=1;nNumz<=iRangeTrc[1]-nTempTrc;nNumz++){
				for(int nNumx=0;nNumx<iPCAno;nNumx++){
					Cpy_1D_Flt(fPCAcubeResult[nNumx][nTempTrcCount],fPCAcubeResult[nNumx][nTempTrcCount+nNumz],iSizeSeis);
				}
			}
		}
		nTempTrcCount=nTempTrcCount+iAsWin;
	}
	// for(int k=0;k<iPCAno;k++){
		// for(int p=0;p<iRangeTrc[1]-iRangeTrc[0]+1;p++){
			// for(int q=0;q<iSizeSeis;q++){
				// fPCAcubeResult[k][p][q]=(fPCAcubeResult[k][p][q]*(max-min)+min);
			// }
		// }
	// }
	unsigned char *fTrcOut=new unsigned char [iSizeSeis*sizeof(float)];
	for(int nNumx=0;nNumx<iPCAno;nNumx++){
		nTempTrcCount=0;
		for(int nTempTrc=iRangeTrc[0];nTempTrc<=iRangeTrc[1];nTempTrc++){
			fileidsegypca[nNumx].write((char *)fTrcHd[0][nTempTrcCount],240);
			Flt_Bin_IBM_BE(fPCAcubeResult[nNumx][nTempTrcCount],fTrcOut,iSizeSeis);
			fileidsegypca[nNumx].write((char *)fTrcOut,iSizeSeis*sizeof(float));
			nTempTrcCount=nTempTrcCount+1;
		}
	}
	delete fTrcOut;
	Del_3D_Flt(fPCAcubeResult,iPCAno,iRangeTrc[1]-iRangeTrc[0]+1,iSizeSeis);
	nTempTrcCount=0;
	for(int nNumx=0;nNumx<iHalfLinRange;nNumx++){
		Cpy_2D_Chr(fTrcHd[nNumx+1],fTrcHd[nNumx],iRangeTrc[1]-iRangeTrc[0]+1,240);
	}

	// End of Part C

	// Start of Part A

	fileidsegy.close();
	for(int nNum=0;nNum<iPCAno;nNum++){
		fileidsegypca[nNum].close();
	}

	return 0;
}

