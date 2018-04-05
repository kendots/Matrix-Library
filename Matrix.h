#ifndef MATRIX_H
#define MATRIX_H


#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include"Kmath.h"
#include"Vector.h"


typedef struct
{
	int r;
	int c;
	double ** v;
} matrix;

matrix mnull={0,0,NULL};

void     MatInit (matrix *);
void     MatFree (matrix *);
matrix   Matrix (int , int , ...);
int      MatZero (matrix);
void     MatGet (matrix);
void     MatPrint (matrix);
void     MatAddRow (matrix, matrix, matrix *);
void     MatAddCol (matrix, matrix, matrix *);
matrix   I (int);
matrix   MatCst (int, int, double);
matrix   MatRand (int, int, int);
void     MatCpy(matrix, matrix);
void     MatAdd (matrix, matrix, matrix *);
void     MatScalar (matrix, double, matrix *);
void     MatProduct (matrix, matrix, matrix *);
int	 MatCmp (matrix, matrix);
void  	 Transpose (matrix, matrix *);
double	 Trace (matrix);
void  	 Ref (matrix, matrix *);
double	 Det (matrix);
void   	 Inverse (matrix, matrix *);
void 	 lineq (matrix, vector, vector *);


//allocate a matrix
void MatInit (matrix *a)
{
	if (a->v)
	{
		puts("Error in MatInit: Pointer might be allocated");
		MatFree(a);
	}
	int i;
	a->v= malloc(a->r*sizeof(double *));
	for(i=0; i<a->r; i++)
		a->v[i]= calloc(a->c, sizeof(double));
}


//free all memory in the matrix
void MatFree (matrix * a)
{
	if (a->v)
	{
		int i;
		for (i=0; i<a->r; i++)
			free(a->v[i]);
		free(a->v);
	}
	else
		puts("Error in MatFree: Pointer is NULL");
}


matrix Matrix (int r, int c, ...)
{
	int i,j;
	matrix a={r,c,NULL};
	MatInit(&a);
	va_list val;
	c=c*r;
	va_start(val,c);
	c=c/r;
	for (i=0; i<r; i++)
	{
		for (j=0; j<c; j++)
			a.v[i][j]=va_arg(val, double);
	}
	va_end(val);
	return a;
}


int MatZero(matrix a)
{
	int i,j,k=0;
	for (i=0; i<a.r; i++)
		for (j=0; j<a.c; j++)
		{
			if (a.v[i][j]<eps && a.v[i][j]>-eps)
			{ 
				a.v[i][j]=0;
				k++;
			}
		}
	return k;
}


//Input a matrix from stdin
void MatGet (matrix a)
{
	int i,j;
	for (i=0; i<a.r; i++)
		for (j=0; j<a.c; j++)
			scanf("%lf",&a.v[i][j]);
}

//Print a matrix to stdout
void MatPrint (matrix a)
{
	int i,j;
	puts("");
	for (i=0; i<a.r; i++)
	{
		for (j=0; j<a.c; j++)
		{
			if (a.v[i][j]<eps && a.v[i][j]>-eps) a.v[i][j]=0;
			printf("%g \t",a.v[i][j]);
		}
		puts("");
	}
}


void MatAddRow(matrix a, matrix b, matrix * c)
{
	int i,j;

	if (c->v==NULL)
	{
		c->r=a.r+b.r;
		c->c=a.c;
		MatInit(c);
	}

	for (i=0; i<a.r; i++)
		for (j=0; j<a.c; j++)
		{
			c->v[i][j]=a.v[i][j];
			if (i<b.r)
				c->v[i+a.r][j]=b.v[i][j];
		}
}


void MatAddCol(matrix a, matrix b, matrix * c)
{
	int i,j;

	if (c->v==NULL)
	{
		c->r=a.r;
		c->c=a.c+b.c;
		MatInit(c);
	}
	for (i=0; i<a.r; i++)
	{
		for (j=0; j<a.c; j++)
		{
			c->v[i][j]=a.v[i][j];
			if (j<b.c)
				c->v[i][j+a.c]=b.v[i][j];
		}
	}
}


//a=In
matrix I (int n)
{
	int i,j;
	matrix a={n,n,NULL};
	MatInit(&a);
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			if (i==j) a.v[i][j]=1;
			else a.v[i][j]=0;
		}
	}
	return a;
}


//a[i][j]=x for every i,j
matrix MatCst (int m, int n, double x)
{
	int i,j;
	matrix a={m,n,NULL};
	MatInit(&a);
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			a.v[i][j]=x;
	return a;
}


matrix MatRand(int m, int n, int x)
{
	int i,j;
	matrix a={m,n,NULL};
	MatInit(&a);
	for (i=0; i<m; i++)
	{
		for (j=0; j<n; j++)
		{
			a.v[i][j]=rand()%x;
		}
	}
	return a;
}


void MatCpy (matrix a, matrix b)
{
	int i,j;
	for (i=0; i<a.r; i++)
		for (j=0; j<a.c; j++)
			a.v[i][j]=b.v[i][j];
}

//c=a+b
void MatAdd (matrix a, matrix b, matrix * c)
{
	int i,j;

	if (c->v==NULL)
	{
		c->r=a.r;
		c->c=a.c;
		MatInit(c);
	}

	for (i=0; i<a.r; i++)
	{
		for (j=0; j<a.c; j++)
		{
			c->v[i][j]=a.v[i][j]+b.v[i][j];
			if (c->v[i][j]<eps && c->v[i][j]>-eps) c->v[i][j]=0;
		}
	}
}


//z = x*a (x in R)
void MatScalar (matrix a, double x, matrix * z)
{
	int i,j;

	if (z->v==NULL)
	{
		z->r=a.r;
		z->c=a.c;
		MatInit(z);
	}

	for (i=0; i<a.r; i++)
		for (j=0; j<a.c; j++)
			z->v[i][j]=x*a.v[i][j];
}


//c = a*b (matrices)
void MatProduct(matrix a, matrix b, matrix * c)
{
	int i,j,k;
	matrix z={a.r, b.c,NULL};
	MatInit(&z);

	if (c->v==NULL)
	{
		c->r=a.r;
		c->c=a.c;
		MatInit(c);
	}

	for (i=0; i<a.r; i++)
	{
		for (j=0; j<b.c; j++)
		{
			z.v[i][j]=0;
			for (k=0; k<b.r; k++)
				z.v[i][j]+=a.v[i][k]*b.v[k][j];
			if (z.v[i][j]<eps && z.v[i][j]>-eps) 
				z.v[i][j]=0;
		}
	}
	MatCpy(*c,z);
	MatFree(&z);
}



int MatCmp(matrix a, matrix z)
{
	if ((a.r- z.r) || (a.c- z.c))
		return 0;
	int i,j;
	double x;
	for (i=0; i<a.r; i++)
		for (j=0; j<a.c; j++)
		{
			x=a.v[i][j]-z.v[i][j];
			if (x>eps || x<-eps)
			{
				return 0;
			}
		}
	return 1;
}


void Transpose (matrix a, matrix * z)
{
	int i,j;
	matrix b={a.c,a.r,NULL};
	MatInit(&b);

	if (z->v==NULL)
	{
		z->r=a.c;
		z->c=a.r;
		MatInit(z);
	}

	for (i=0; i<a.c; i++)
		for (j=0; j<a.r; j++)
		{
			b.v[i][j]=a.v[j][i];
		}
	MatCpy(*z,b);
	MatFree(&b);
}



//Trace of a square matrix is the sum of the elements on the diagonal
double Trace (matrix a)
{
	if (a.r-a.c)
	{
		puts("Error in Trace Function:\nThe matrix should be a square matrix");
		return 0;
	}

	int i,s=0;
	for (i=0; i<a.r; i++)
		s+=a.v[i][i];
	return s;
}


//Row echelon form
void Ref (matrix a, matrix * z)
{
	int i,j,x,y,m=a.r,n=a.c;
	double b,t;
	if (z->v==NULL)
	{
		z->r=a.r;
		z->c=a.c;
		MatInit(z);
	}
	MatCpy(*z,a);

	i=0; j=0;
	for (i=0,y=0; i<m; i++,y++)
	{
		if (z->v[y][i]<eps && z->v[y][i]>-eps)
		{
			for (x=y; x<m; x++)
				if (z->v[x][i]>eps || z->v[x][i]<-eps )
				{
					for (j=0; j<n; j++)
					{
						t=z->v[y][j];
						z->v[y][j]=z->v[x][j];
						z->v[x][j]=t;
					}
					break;
				}
			if (x==m)
			{
				y--;
				continue;
			}
		}

		for (x=y+1; x<m; x++)
		{
			b=z->v[x][i]/z->v[y][i];
			int c=b;
			if (b-c>eps)
				for (j=y; j<n; j++)
					z->v[x][j]=(z->v[y][i]*z->v[x][j])-b* (z->v[y][i]*z->v[y][j]);
			else
				for (j=y; j<n; j++)
					z->v[x][j]=z->v[x][j]-b*z->v[y][j];
		}
	}
}


//Determinant
double Det(matrix a)
{
	if (a.r- a.c)
	{
		puts("Error in the Det Function:\nThe matrix should be a square matrix");
		return 0;
	}

	int i,j,z,x=1,l;
	double s;
	matrix b={a.r-1,a.r-1,NULL};
	MatInit(&b);

	if (a.r==1)
		s = a.v[0][0];
	if (a.r==2)
		s = a.v[0][0]*a.v[1][1]-a.v[0][1]*a.v[1][0];

	if (a.r>2)
	{
		s=0;
		for (z=0; z<a.r; z++)
		{
			l=0;
			for (i=1; i<a.r; i++)
			{
				for (j=0; j<a.r; j++)
				{
					if(j==z) continue;
					b.v[i-1][l]=a.v[i][j];
					l++;
				}
				l=0;
			}
			s=s+x*a.v[0][z]*Det(b);
			x=-x;
		}
	}

	MatFree(&b);
	return s;
}


void Inverse (matrix a, matrix * z)
{
	double x=Det(a),y;
	if (x<eps && x>-eps)
	{
		puts("Error in the Inverse Function:\nThe matrix is not Invertable");
		return;
	}

	if (z->v==NULL)
	{
		z->r=a.r;
		z->c=a.c;
		MatInit(z);
	}

	if (a.r<2)
	{
		z->v[0][0]=1.0/x;
		if (z->v[0][0]<eps && z->v[0][0]>-eps) z->v[0][0]=0;
	}

	else
	{
		int i,j,k,l,t,r,n=a.r;
		matrix b={n-1, n-1,NULL},c={n,n,NULL};
		MatInit(&b);
		MatInit(&c);


		for (i=0; i<n; i++)
		{
			for (j=0; j<n; j++)
			{
				for (k=0, t=0; k<n; k++)
				{
					if (k==i) continue;
					for (l=0, r=0; l<n; l++)
					{
						if (l==j) continue;
						b.v[t][r]=a.v[k][l];
						r++;
					}
					t++;
				}

				if ((i+j)%2) y=-x;
				else y=x;
				c.v[j][i]=Det(b)/y;
				if (c.v[j][i]<eps && c.v[j][i]>-eps) c.v[j][i]=0;
			}
		}

		MatCpy(*z,c);
		MatFree(&c);
		MatFree(&b);
	}
}


//finds the solution vector x of the equation a*x=b if possible
void lineq(matrix a, vector b, vector * x)
{
	int n=a.r;
	double z=Det(a);

	if (ab(z)<eps)
	{
		puts("Can't be solved, Matrix is singular");
		return ;
	}

	if (x->t==NULL)
	{
		x->n=a.r;
		VectInit(x);
	}

	int i,j,k;
	matrix c={n,n,NULL};
	MatInit(&c);
	for (k=0; k<n; k++)
	{
		for (i=0; i<n; i++)
			for (j=0; j<n; j++)
			{
				if (j==k) c.v[i][j]=b.t[i];
				else 	  c.v[i][j]=a.v[i][j];

			}
		x->t[k]=Det(c)/z;
		if (ab(x->t[k])<eps) x->t[k]=0;
	}
	MatFree(&c);
}

#endif
