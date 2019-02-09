#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <Eigen/Dense>
#include <Segy.h>
#include <Array.h>

// diffraction extraction on post-stack seismic data
// using principal component analysis(PCA)
// Author @Jianhang Yin 2016
int main(){
	using namespace std;
	using namespace Eigen;
	int iRangeLin[2],iRangeTrc[2],iSizeSlipWin,iLag[3],iStepDis,iSizeSeis,iAsWin,iPCAno;
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
	cout<<"the processing line NO. is : "<<nLin-1<<endl;

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
			double ***fCalCube=New_3D_Dbl((2*(iLag[0]/iStepDis)+1)*(2*(iLag[1]/iStepDis)+1)*(2*(iLag[2]/iStepDis)+1),2*iSizeSlipWin+1,(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1));
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
									fCalCube[nNum-1][nCalLin][nCalTim+nCalTrc*(2*iSizeSlipWin+1)]=fCalLin[nTemp3][nTemp2][nTemp1];
								}
							}
						}
					}
				}
			}
			// calculate the covariance matrix fCov (d * M)
			double **fCov=New_2D_Dbl(nNum,nNum);
			for(int nNumx=0;nNumx<nNum;nNumx++){
				for(int nNumy=0;nNumy<nNum;nNumy++){
					double sumx=0,sumy=0,tempxy=0,avex,avey;
					for(int i=0;i<2*iSizeSlipWin+1;i++){
						for(int j=0;j<(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1);j++){
							sumx+=fCalCube[nNumx][i][j];
							sumy+=fCalCube[nNumy][i][j];
						}
					}
					avex=sumx/((2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1));
					avey=sumy/((2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1));
					for(int i=0;i<2*iSizeSlipWin+1;i++){
						for(int j=0;j<(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1);j++){
							tempxy+=(fCalCube[nNumx][i][j]-avex)*(fCalCube[nNumy][i][j]-avey);
						}
					}
					fCov[nNumx][nNumy]=tempxy/((2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1)-1);
				}
			}
			MatrixXd fCovx(nNum,nNum);
			for(int i=0;i<nNum;i++){
				for(int j=0;j<nNum;j++){
					fCovx(i,j)=fCov[i][j];
				}
			}
			// calculate the eigenvalue and eigenvector of fCov
			EigenSolver<MatrixXd> es(fCovx);
			MatrixXcd Vc=es.eigenvectors();
			VectorXcd Dc=es.eigenvalues();
			double **V=New_2D_Dbl(nNum,nNum);
			double *D=new double [nNum];
			for(int i=0;i<nNum;i++){
				for(int j=0;j<nNum;j++){
					V[i][j]=Vc(j,i).real();
				}
			}
			for(int i=0;i<nNum;i++){
				D[i]=Dc[i].real();
			}
			double *Vtemp=new double [nNum];
			double Dtemp;
			// sort the eigenvalue and eigenvector
			for(int nNumy=0;nNumy<nNum-1;nNumy++){
				for(int nNumx=nNumy;nNumx<nNum-1;nNumx++){
					if(D[nNumx+1]>D[nNumy]){
						Cpy_1D_Dbl(V[nNumx+1],Vtemp,nNum);
						Cpy_1D_Dbl(V[nNumy],V[nNumx+1],nNum);
						Cpy_1D_Dbl(Vtemp,V[nNumy],nNum);
						Dtemp=D[nNumx+1];
						D[nNumx+1]=D[nNumy];
						D[nNumy]=Dtemp;
					}
				}
			}
			// map the sorted eigenvector to data
			double **fEigVet=New_2D_Dbl(iPCAno,nNum);
			for(int nNumx=0;nNumx<iPCAno;nNumx++)
				Cpy_1D_Dbl(V[nNumx],fEigVet[nNumx],nNum);
			double *fDiffCube=new double [nNum];
			double ***fPCAcubeWhole=New_3D_Dbl(2*iSizeSlipWin+1,(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1),iPCAno);
			for(int nNumy=0;nNumy<2*iSizeSlipWin+1;nNumy++){
				for(int nNumz=0;nNumz<(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1);nNumz++){
					for(int nNumx=0;nNumx<nNum;nNumx++){
						fDiffCube[nNumx]=fCalCube[nNumx][nNumy][nNumz];
					}
					for(int i=0;i<iPCAno;i++){
						for(int j=0;j<nNum;j++){
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
			// reconstruct the seismic data
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
			Del_3D_Dbl(fCalCube,(2*(iLag[0]/iStepDis)+1)*(2*(iLag[1]/iStepDis)+1)*(2*(iLag[2]/iStepDis)+1),2*iSizeSlipWin+1,(2*iSizeSlipWin+1)*(2*iSizeSlipWin+1));
			Del_2D_Dbl(fCov,nNum,nNum);
			Del_2D_Dbl(V,nNum,nNum);
			delete D;
			delete Vtemp;
			Del_2D_Dbl(fEigVet,iPCAno,nNum);
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
	// output the result
	delete fTrcOut;
	Del_3D_Flt(fPCAcubeResult,iPCAno,iRangeTrc[1]-iRangeTrc[0]+1,iSizeSeis);
	nTempTrcCount=0;
	for(int nNumx=0;nNumx<iHalfLinRange;nNumx++){
		Cpy_2D_Chr(fTrcHd[nNumx+1],fTrcHd[nNumx],iRangeTrc[1]-iRangeTrc[0]+1,240);
	}

	// End of Part C

	// Start of Part D

	fileidsegy.close();
	for(int nNum=0;nNum<iPCAno;nNum++){
		fileidsegypca[nNum].close();
	}
	return 0;
}
