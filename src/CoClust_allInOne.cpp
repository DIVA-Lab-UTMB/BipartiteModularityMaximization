#define _XOPEN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <new>          // ::operator new[]
#include <Rcpp.h>
using namespace Rcpp;

#define eps 1E-4
#define alpha 1E-3
#define beta 1E-3
#define toler 1E-5

/* Structure definition*/
typedef struct{
	double modularity;
	int *community;
}comstruct;


/* Prototypes */
double smprng(void);
comstruct barber(double (*rng)(void), int run, double **a, int na, int nt);
double BBbisection(double (*rng)(void), double **b, int *g, int *s, int n);
double eigen(double **B, double *v, int n);
double abmax(double *v,int n);
double BBKLtuning(double (*rng)(void), double **B, int *s, int na, int nt);
double BBFinaltuning(double (*rng)(void), double **B, int **S, int na, int nt, int *ng);
double BBAgglomeration(double **B, int na, int nt, int *ng, int *gs, int **G);

double BBModularity(int na, int nt, int *s,double **b);
double BBAddmatrix(double **b,int na, int nt);
void trans(int **S, int ng, int *gs, int **G);
void invertrans(int **S, int n, int *ng, int *gs, int **G);
void con_to_G(int *con, int n, int *ng, int *gs, int **G);

int _NA_;
double m,Q;
static int **G;

double smprng()
{
	return unif_rand();
}

comstruct barber(double (*rng)(void), int run, double **a, int na, int nt) {
    int i, j, l,  ii, flag=0, n, ng,ngtemp,gstemp, nn, np;
	int *s,*gs,*Gitemp, **S, *comm;
	double max1,max2,dQ,Qmax=0;
	double **B, *k, *d, **b;
	comstruct result;

	result.modularity=0;

	_NA_=na;
	n=na+nt;

	//srand((long int) time(NULL)); done using GetRNGstate()

	G=(int **)malloc(sizeof(int *)*n);
	gs=(int *)malloc(sizeof(int)*n);
	s=(int *)malloc(sizeof(int)*n);
	Gitemp=(int *)malloc(sizeof(int)*n);
	S=(int **)malloc(sizeof(int *)*n);
	comm=(int *)malloc(sizeof(int)*n);
	result.community=(int *)malloc(sizeof(int)*n);


	k=(double *)malloc(sizeof(double)*na);
	d=(double *)malloc(sizeof(double)*nt);
	b=(double **)malloc(sizeof(double *)*na);
	for (i=0; i<na; i++) b[i]=(double *)malloc(sizeof(double)*nt);
	B=(double **)malloc(sizeof(double *)*n);
	for (i=0; i<n; i++) {
		B[i]=(double *)malloc(sizeof(double)*n);
		G[i]=(int *)malloc(sizeof(int)*n);
		S[i]=(int *)malloc(sizeof(int)*n);
		result.community[i]=0;
	}


	m=0;

	for (i=0; i<na; i++) {
		k[i]=0;
		for (j=0; j<nt; j++)
			k[i]=k[i]+a[i][j];
		m=m+k[i];
	}
	if (fabs(m)<1e-6) {
		Rprintf("Adjacency matrix is 0!!");
		return result;
	}

	for (i=0; i<nt; i++) {
		d[i]=0;
		for (j=0; j<na; j++)
			d[i]=d[i]+a[j][i];
	}

	for (i=0; i<na; i++)
		for (j=0; j<nt; j++) {
			b[i][j]=a[i][j]-k[i]*d[j]/m;
			if (b[i][j]!=0) flag=1;
		}

	if (flag==0) {
		Rprintf("Modularity matrix is 0!");
		return result;
	}


	for (ii=0; ii<run; ii++) {
		Q=0;
		max1=-2;
		max2=-1;

		for (i=0; i<n; i++) {
			s[i]=0;
			gs[i]=0;
			Gitemp[i]=0;
			G[0][i]=i;
		}

		for (i=1; i<n; i++)
			for (j=0; j<n; j++)
				G[i][j]=0;
		nn=np=0;
		ng=1;
		gs[0]=n;


		while (max2-max1>toler) {
		  Rcpp::checkUserInterrupt();// check interrupt in R
			max1=max2;
			ngtemp=ng;

			for (i=0; i<ngtemp; i++){
				dQ=BBbisection(rng,b,G[i],s,gs[i]);
				if (dQ>toler) {
					ng++;

					gstemp=gs[i];

					for (j=0; j<gstemp; j++) {
						Gitemp[j]=G[i][j];
						G[i][j]=0;
					}


					for (j=0; j<gstemp; j++) {
						if (s[j]==1) {
							G[i][np]=Gitemp[j];
							np++;
						}
						else{
							G[ng-1][nn]=Gitemp[j];
							nn++;
						}
					}

					gs[i]=np;
					gs[ng-1]=nn;
					np=0;
					nn=0;

					for (j=0; j<gstemp; j++){
						s[i]=0;
						Gitemp[j]=0;
					}

					Q=Q+dQ;
				}
			}

			for (i=0; i<n; i++)
				for (j=i; j<n; j++)
					S[i][j]=S[j][i]=0;

			trans(S,ng,gs,G);
			Q=Q+BBFinaltuning(rng,b,S,na,nt,&ng);
			invertrans(S,n,&ng,gs,G);
			Q=Q+BBAgglomeration(b,na,nt,&ng,gs,G);

			max2=Q;

		}
		//Rprintf("Run #%d, Modularity is %15.10f\n", ii+1, Q);
		if (Q>Qmax) {
			Qmax=Q;
			for (i=0; i<ng; i++)
				for (j=0; j<gs[i]; j++)
					comm[G[i][j]]=i;
		}

	}

	result.modularity=Qmax;
	for(i=0; i<n; i++) result.community[i]=comm[i];

	for (i=0; i<n; i++) {
		free(G[i]);
		free(B[i]);
		free(S[i]);
	}
	for (i=0; i<na; i++) free(b[i]);

	free(S);
	free(G);
	free(B);
	free(gs);
	free(comm);
	free(Gitemp);
	free(s);
	free(b);
	free(k);
	free(d);


    return result;
}


double BBbisection(double (*rng)(void), double **b, int *g, int *s, int n){
	int i,j,l, nn=0, np=0, na=0, nt=0, min, flag=0;
	double x0,x,bik,dQ, xa, xt;
	double *va, *vt, **c, **bb;

	for (i=0; i<n; i++) {
		if (g[i]<_NA_) na++;
		else nt++;
	}
	if (na<=1 || nt<=1) return 0;

	va=(double *)malloc(sizeof(double)*na);
	vt=(double *)malloc(sizeof(double)*nt);
	c=(double **)malloc(sizeof(double *)*na);
	for (i=0; i<na; i++) c[i]=(double *)malloc(sizeof(double)*nt);

	for (i=0; i<na; i++) {
		for (j=0; j<nt; j++) {
			c[i][j]=b[g[i]][g[j+na]-_NA_];
		}
	}



	if (na<=nt) {
		min=na;
		bb=(double **)malloc(sizeof(double *)*na);
		for (i=0; i<na; i++) bb[i]=(double *)malloc(sizeof(double)*na);

		for (i=0; i<na; i++)
			for (j=i; j<na; j++) {
				bb[i][j]=0;
				for (l=0; l<nt; l++)
					bb[i][j]=bb[i][j]+c[i][l]*c[j][l];
				if(bb[i][j]!=0) flag=1;
				bb[j][i]=bb[i][j];
			}

		if (flag==0) {
			for (i=0; i<na; i++) {
				free(c[i]);
				free(bb[i]);
			}
			free(va);
			free(vt);
			free(c);
			free(bb);
			return 0;
		}

		xa=eigen(bb, va, na);
		for (i=0; i<nt; i++) {
			vt[i]=0;
			for (j=0; j<na; j++)
				vt[i]=vt[i]+c[j][i]*va[j];
			vt[i]=vt[i]/xa;
		}
	}

	else {
		min=nt;
		bb=(double **)malloc(sizeof(double *)*nt);
		for (i=0; i<nt; i++) bb[i]=(double *)malloc(sizeof(double)*nt);

		for (i=0; i<nt; i++)
			for (j=i; j<nt; j++) {
				bb[i][j]=0;
				for (l=0; l<na; l++)
					bb[i][j]=bb[i][j]+c[l][i]*c[l][j];
				if(bb[i][j]!=0) flag=1;
				bb[j][i]=bb[i][j];
			}
		if (flag==0) {
			for (i=0; i<na; i++) {
				free(c[i]);
			}
			for (i=0; i<nt; i++) {
				free(bb[i]);
			}
			free(va);
			free(vt);
			free(c);
			free(bb);
			return 0;
		}

		xt=eigen(bb, vt, nt);
		for (i=0; i<na; i++) {
			va[i]=0;
			for (j=0; j<nt; j++)
				va[i]=va[i]+c[i][j]*vt[j];
			va[i]=va[i]/xt;
		}
	}

	for (i=0; i<na; i++){
		if (va[i]>0.0) s[i]=1;
		else if (va[i]<0.0) s[i]=-1;
		else if (rng()>0.5) s[i]=1;
		else s[i]=-1;
	}

	for (i=0; i<nt; i++){
		if (vt[i]>0.0) s[i+na]=1;
		else if (vt[i]<0.0) s[i+na]=-1;
		else if (rng()>0.5) s[i+na]=1;
		else s[i+na]=-1;
	}

	BBKLtuning(rng,c,s,na,nt);
	dQ=(BBModularity(na,nt,s,c)-BBAddmatrix(c,na,nt))/(2*m);

	for (i=0; i<na; i++)
		free(c[i]);
	for (i=0; i<min; i++) {
		free(bb[i]);
	}
	free(bb);
	free(c);
	free(va);
	free(vt);
	return dQ;
}


double eigen(double **B, double *v, int n) {
	int i,j,N=0;
	double w1=0,w2=10,sum;
	double *u,*ad,*temp;

	u=(double *)malloc(sizeof(double)*n);
	ad=(double *)malloc(sizeof(double)*n);
	temp=(double *)malloc(sizeof(double)*n);

	for (i=0; i<n; i++)
		ad[i]=B[i][i];//to preserve the diagonal elements of **a

	for (i=0; i<(int)(n/2); i++) {//the initial u and v are random
		u[i]=1;
		v[i]=1;
	}
	for (i=(int)(n/2); i<n; i++) {
		u[i]=2;
		v[i]=2;
	}
	for (i=0; i<n; i++) temp[i]=v[i]-u[i];
	w1=abmax(u,n);
	while (fabs(w1-w2)>eps|| fabs(abmax(temp,n))>eps) {
		N++;
		if (N>300) {
			for (i=0; i<n; i++)
				B[i][i]=B[i][i]+alpha;
			N=0;
		}
		for (i=0; i<n; i++) {
			sum=0;
			for (j=0; j<n; j++)
				sum+=B[i][j]*v[j];
			u[i]=sum;
		}
		w1=w2;
		w2=abmax(u, n);
		while (w2==0) {
			for (i=0; i<n; i++)
				v[i]=i+0.5;//change
			for (i=0; i<n; i++) {
				sum=0;
				for (j=0; j<n; j++)
					sum+=B[i][j]*v[j];
				u[i]=sum;
			}
			w2=abmax(u, n);
		}

		for (i=0; i<n; i++) {
			u[i]=u[i]/w2;
			temp[i]=u[i]-v[i];
			v[i]=u[i];
		}

	}

	for (i=0; i<n; i++)
		B[i][i]=ad[i];

	free(u);
	free(ad);
	free(temp);
	return w2;
}

double BBKLtuning(double (*rng)(void), double **B, int *s, int na, int nt){
	int i, j, k, p, r, n;
	int *al, *c, *d;
	double dQ,qq,max1,max2;
	double *W, *dW, *dq, *q;

	n=na+nt;
	max1=-2;
	max2=0;

	al=(int *)malloc(sizeof(int)*n);
	c =(int *)malloc(sizeof(int)*n);
	d =(int *)malloc(sizeof(int)*n);

	W =(double *)malloc(sizeof(double)*n);
	dW=(double *)malloc(sizeof(double)*n);
	dq=(double *)malloc(sizeof(double)*n);
	q =(double *)malloc(sizeof(double)*n);

	while (max2-max1>toler) {
		max1=max2;
		q[0]=max2;

		for (i=0; i<n; i++) {
			al[i]=-1;
			c[i]=0;//c[i]=0 means node i not moved yet
			dW[i]=0;
			W[i]=0;
		}

		for (i=0; i<na; i++) {
			for (j=0; j<nt; j++)
				W[i]=W[i]+(double)s[j+na]*B[i][j];
			dq[i]=-(double)s[i]*W[i]/m;
		}
		for (i=na; i<n; i++) {
			for (j=0; j<na; j++)
				W[i]=W[i]+(double)s[j]*B[j][i-na];
			dq[i]=-s[i]*W[i]/m;
		}

		/* to pick up the largest dQ*/
		dQ=-10;
		for (i=0; i<n; i++) {
			if (dq[i]>dQ) {
				p=1;
				d[0]=i;
				dQ=dq[i];
			}
			else if (dq[i]==dQ) {
			    d[p]=i;
				p++;
			}
		}

		if (p==1) {
			al[0]=d[0];
			c[al[0]]=1;//c[i]=1 means the i-th node has been choosen
		}
		else {
			p=(int)p*rng();
			al[0]=d[p];
			c[al[0]]=1;
		}

		/*moving nodes*/
		for (i=1; i<n; i++) {
			q[i]=q[i-1]+dQ;
			dQ=-10;
			for (j=0; j<na; j++)
				if (c[j]==0) {
					if (al[i-1]>=na)
						dW[j]=(double)s[al[i-1]]*B[j][al[i-1]-na]*2;//!!!
					else dW[j]=0;
				}
			for (j=na; j<n; j++)
				if (c[j]==0) {
					if (al[i-1]<na)
						dW[j]=(double)s[al[i-1]]*B[al[i-1]][j-na]*2; //z
					else dW[j]=0;
				}
			for (j=0; j<n; j++)
				if (c[j]==0) {
					W[j]=W[j]-dW[j];
					dq[j]=-(double)s[j]*W[j]/m;

					if (dq[j]>dQ) {
						p=1;
						d[0]=j;
						dQ=dq[j];
					}
					else if (dq[j]==dQ) {
						d[p]=j;
						p++;
					}
				}
			if (p==1) {
				al[i]=d[0];
			}
			else {
				p=(int)p*rng();
				al[i]=d[p];
			}
			c[al[i]]=1;
		}

		/*to find the largest configuration*/
		qq=-1;
		for (i=0; i<n; i++) {
			if (q[i]>qq) {
				p=1;
				d[0]=i;
				qq=q[i];
			}
			else if (q[i]==qq) {
				d[p]=i;
				p++;
			}
		}
		if (p==1)
			r=d[0];
		else {
			p=(int)p*rng();
			r=d[p];
		}//the #r-th state has the largest configuration
		for (i=0; i<r; i++)
			s[al[i]]=-s[al[i]];
		dQ=q[r]-max2;
		max2=q[r];
	}

	free(al);
	free(c);
	free(d);
	free(W);
	free(dW);
	free(dq);
    free(q);

	return max1;
}

double BBFinaltuning(double (*rng)(void), double **B, int **S, int na, int nt, int *ng){
	int i,j,k, p, r, n;
	int *con, *al, *be, *c, *dn, *dg;
	double dQ, qq, max1, max2;
	double **W, **dW, **dq, *q;

	n=na+nt;
	max1=-2;
	max2=0;

	con=(int *)malloc(sizeof(int)*n);//one dimensional configuration
	al=(int *)malloc(sizeof(int)*n);//alpha, record the node to be moved each time
	be=(int *)malloc(sizeof(int)*n);//beta, record the group each node being moved to
	c =(int *)malloc(sizeof(int)*n);
	dn=(int *)malloc(sizeof(int)*n*n);
	dg=(int *)malloc(sizeof(int)*n*n);

	W=(double **)malloc(sizeof(double *)*n);
	dW=(double **)malloc(sizeof(double *)*n);
	dq=(double **)malloc(sizeof(double *)*n);
	q =(double *)malloc(sizeof(double)*n);

	for (i=0; i<n; i++) {
		W[i]=(double *)malloc(sizeof(double)*n);
		dW[i]=(double *)malloc(sizeof(double)*n);
		dq[i]=(double *)malloc(sizeof(double)*n);
	}

	for (i=0; i<n; i++) {
		j=0;
		while (S[i][j]==0) {
			j++;
		}
		con[i]=j;
	}

	while (max2-max1>toler) {
		max1=max2;
		q[0]=max2;

		for (i=0; i<n; i++) {
			al[i]=-1;
			c[i]=0;
			for (j=0; j<n; j++) {
				dW[i][j]=0;
				W[i][j]=0;
			}
		}

		for (i=0; i<na; i++)
			for (j=0; j<*ng+1; j++)
				for (k=0; k<nt; k++)
					W[i][j]=W[i][j]+(double)B[i][k]*S[k+na][j];
		for (i=na; i<n; i++)
			for (j=0; j<*ng+1; j++)
				for (k=0; k<na; k++)
					W[i][j]=W[i][j]+(double)B[k][i-na]*S[k][j];

		for (i=0; i<n; i++)
			for (j=0; j<*ng+1; j++) {
				if (j==con[i]) dq[i][j]=0;
				else {
					dq[i][j]=(-W[i][con[i]]+W[i][j])/m;
				}
			}

		dQ=-10;
		for (i=0; i<n; i++)
			for (j=0; j<*ng+1; j++)
				if (j!=con[i]) {
					if (dq[i][j]>dQ) {
						p=1;
						dn[0]=i;
						dg[0]=j;
						dQ=dq[i][j];
					}
					else if (dq[i][j]==dQ) {
						dn[p]=i;
						dg[p]=j;
						p++;
					}
				}

		if (p==1) {
			al[0]=dn[0];
			c[al[0]]=1;
			be[0]=dg[0];
		}
		else {
			p=(int)p*rng();
			al[0]=dn[p];
			c[al[0]]=1;
			be[0]=dg[p];
		}

		if (be[0]==*ng) *ng=*ng+1;

		/*moving nodes*/
		for (i=1; i<n; i++) {
			q[i]=q[i-1]+dQ;
			dQ=-10;

			for (j=0; j<na; j++)
				if (c[j]==0)
					if (al[i-1]>=na)
						for (k=0; k<*ng+1; k++) {
							if (k==con[al[i-1]]) dW[j][k]=dW[j][k]+B[j][al[i-1]-na];
							else if (k==be[i-1]) dW[j][k]=dW[j][k]-B[j][al[i-1]-na];
						}
			for (j=na; j<n; j++)
				if (c[j]==0)
					if (al[i-1]<na)
						for (k=0; k<*ng+1; k++) {
							if (k==con[al[i-1]]) dW[j][k]=dW[j][k]+B[al[i-1]][j-na];
							else if (k==be[i-1]) dW[j][k]=dW[j][k]-B[al[i-1]][j-na];
						}

			for (j=0; j<n; j++)
				if (c[j]==0)
					for (k=0; k<*ng+1; k++) {
						if (k==con[j]) dq[j][k]=0;
						else {
							dq[j][k]=(-W[j][con[j]]+W[j][k]+dW[j][con[j]]-dW[j][k])/m;
						}
					}

			for (j=0; j<n; j++)
				if (c[j]==0)
					for (k=0; k<*ng+1; k++)
						if (k!=con[j]) {
							if (dq[j][k]>dQ) {
								p=1;
								dn[0]=j;
								dg[0]=k;
								dQ=dq[j][k];
							}
							else if (dq[j][k]==dQ) {
								dn[p]=j;
								dg[p]=k;
								p++;
							}
						}

			if (p==1) {
				al[i]=dn[0];
				be[i]=dg[0];
			}
			else {
				p=(int)p*rng();
				al[i]=dn[p];
				be[i]=dg[p];
			}
			c[al[i]]=1;

			if (be[i]==*ng) *ng=*ng+1;
		}



		qq=-1;
		for (i=0; i<n; i++) {
			if (q[i]>qq) {
				p=1;
				dn[0]=i;
				qq=q[i];
			}
			else if (q[i]==qq) {
				dn[p]=i;
				p++;
			}
		}

		if (p==1) r=dn[0];
		else {
			p=(int)p*rng();
			r=dn[p];
		}
		for (i=0; i<r; i++){
			S[al[i]][con[al[i]]]=0;
			S[al[i]][be[i]]=1;
            con[al[i]]=be[i];
		}

		dQ=q[r]-max2;
		max2=q[r];
	}

	free(al);
	free(be);
	free(c);
	free(dn);
	free(dg);
	free(con);
	free(q);

	for (i=0; i<n; i++) {
		free(W[i]);
		free(dW[i]);
		free(dq[i]);
	}
	free(W);
	free(dW);
	free(dq);

	return max1;
}


double BBAgglomeration(double **B, int na, int nt, int *ng, int *gs, int **G) {
	int i,j,k,ii,jj,R,commnum,g1,g2,n;
	double dQ,temp,max;
	int *con, **config;
	double *Qt;

	n=na+nt;
	commnum=*ng;
	Qt=(double *)malloc(sizeof(double)*n*n/2);
	con=(int *)malloc(sizeof(int)*n);
	config=(int **)malloc(sizeof(int *)*commnum);
	for (i=0; i<commnum; i++) config[i]=(int *)malloc(sizeof(int)*n);

	Qt[0]=0;

	for (i=0; i<commnum-1; i++) {
		/*record */
		for (j=0; j<*ng; j++)
			for (k=0; k<gs[j]; k++)
				config[i][G[j][k]]=j;
		for (j=0; j<n; j++) con[j]=config[i][j];

		temp=-1;
		for (j=0; j<*ng-1;j++ )
			for (k=j+1; k<*ng; k++) {
				dQ=0;
				for (ii=0; ii<gs[j]; ii++)
					for (jj=0; jj<gs[k]; jj++) {
						if (G[j][ii]<na && G[k][jj]>=na) dQ=dQ+B[G[j][ii]][G[k][jj]-na];
						else if (G[j][ii]>=na && G[k][jj]<na) dQ=dQ+B[G[k][jj]][G[j][ii]-na];
					}
				dQ=dQ/m;
				if (dQ>temp) {
					temp=dQ;
					g1=j;
					g2=k;
				}
			}

		Qt[i+1]=Qt[i]+temp;
		for (j=0; j<n; j++) {
			if (con[j]==g2)
				con[j]=g1;
			else if (con[j]>g2)
				con[j]=con[j]-1;
		}
		con_to_G(con,n,ng,gs,G);
	}

	for (i=0; i<*ng; i++)
		for (j=0; j<gs[i]; j++)
			config[commnum-1][G[i][j]]=i;

	temp=-1;
	for (i=0; i<commnum; i++)
		if (Qt[i]>=temp) {
			temp=Qt[i];
			R=i;
		}

	con_to_G(config[R],n,ng,gs,G);
	max=Qt[R];

	for (i=0; i<commnum; i++)
		free(config[i]);
	free(config);
	free(Qt);
    free(con);

	return max;
}


double abmax(double *v,int n){
	int i;
	double temp=v[0];
	for(i=1;i<n;i++)
		if(fabs(v[i])>fabs(temp) || (fabs(v[i])==fabs(temp) && v[i]>temp))//take the positive one
			temp=v[i];
	return temp;
}

double BBModularity(int na, int nt, int *s, double **b){
	int i,j;
	double q=0;

	for (i=0; i<na; i++)
		for (j=0; j<nt; j++)
			q=q+(double)s[i]*b[i][j]*s[j+na];

	return(q);
}



double BBAddmatrix(double **b, int na, int nt){
	int i,j;
	double sum=0;

	for (i=0; i<na; i++)
		for (j=0; j<nt; j++)
			sum=sum+b[i][j];

	return(sum);
}

void trans(int **S, int ng, int *gs, int **G){
	int i,j;

	for (i=0; i<ng; i++)
		for (j=0; j<gs[i]; j++)
			S[G[i][j]][i]=1;

	return;
}

void invertrans(int **S, int n, int *ng, int *gs, int **G){
	int i,j,k,temp,flag=0;

	temp=*ng;

	for (i=0; i<n; i++) {
		gs[i]=0;
		for (j=0; j<n; j++)
			G[i][j]=0;
	}//initialize gs and G

	k=0;
	for (i=0; i<*ng; i++) {
		for (j=0; j<n; j++)
			if (S[j][i]==1) {
				G[k][gs[k]]=j;
				gs[k]++;
				flag=1;
			}
		if (flag==1) k++;
		flag=0;
	}

	*ng=k;

	return;

}

void con_to_G(int *con, int n, int *ng, int *gs, int **G) {
	int i,j;

	*ng=0;
	for (i=0; i<n; i++) {
		gs[i]=0;
		for (j=0; j<n; j++)
			G[i][j]=0;
	}//initialize gs and G

	for (i=0; i<n; i++) {
		G[con[i]][gs[con[i]]]=i;
		gs[con[i]]++;
		if (con[i]>*ng-1)
			*ng=con[i]+1;
	}

	return;
}


//' C++ code to to partition a bipartite network into non-overlapping biclusters, by optimizing bipartite modularity.
//'
//' @param nr Number of rows of incidence matrix.
//' @param nc Number of columns of incidence matrix.
//' @param data Vectorized incidence matrix.
//' @param ITER Number of iterations.
//' @return MODULARITY Modularity value.
//' @return ASSIGN Partition of rows and columns.
//' @keywords internal
// [[Rcpp::export]]
List CoClust(int nr, int nc, NumericVector data, int ITER)
{
	GetRNGstate();//R
	int i,j, n, na=nr, nt=nc,run=ITER;
	double **A;
	comstruct example;

	/* Convert R double to C matrix */
	//A = malloc(na*sizeof(double*));
	//A[0] = malloc(na*nt*sizeof(double));
	A=new double*[na];//cpp
	A[0]=new double[na*nt];//cpp
	for (i=1; i<na; i++) A[i] = A[i-1] + nt;

	for (i=0; i<na; i++)
			for (j=0; j<nt; j++)
			A[i][j] = data[i*nt + j];

    n=na+nt;

    /* Compute communities */
	example=barber(smprng,run,A,na,nt);

	/* Output */
	double MODULARITY = example.modularity;
	IntegerVector ASSIGN(n);
	for (i = 0; i<n; i++) ASSIGN[i] = example.community[i] + 1;
	free(example.community); //z
	//free(A[0]); //z
	//free(A); //z
	delete[] A[0]; //cpp
	delete[] A; //cpp
	PutRNGstate();//R
	//return ASSIGN;
	List L = List::create(Named("MODULARITY") = MODULARITY , _["ASSIGN"] = ASSIGN);
	return L;
}
