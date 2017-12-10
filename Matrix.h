#include<stdlib.h>
#include<time.h>

typedef struct {
int row;
int col;
double ** v;
} matrix;


//free all memory in the matrix
void MatFree (matrix a){
int i;
for (i=0; i<a.row; i++)
free(a.v[i]);
free(a.v);
}


//allocate a matrix
void MatInit (matrix * a){
int i;
if (a->v!=NULL){
puts("It was full");
MatFree(*a);
}

a->v=(double **) malloc(a->row*sizeof(double *));
for(i=0; i<a->row; i++)
a->v[i]=(double *) malloc(a->col*sizeof(double));
}


//Input a matrix from the user
void MatGet (matrix a){
int i,j;
for (i=0; i<a.row; i++)
for (j=0; j<a.col; j++)
scanf("%lf",&a.v[i][j]);
}

//Print a matrix to the user
void MatShow (matrix a){
int i,j;
puts("");
for (i=0; i<a.row; i++){
for (j=0; j<a.col; j++)
printf("%g \t",a.v[i][j]);
puts("");
}}


void AddRow(matrix a, matrix b, matrix c){
if ((a.col-b.col) || (a.col -c.col)){
puts("Error in the AddRow Function:\nThe coloumns of the matrices should be equal");
return;
}

int i,j;
for (i=0; i<a.row; i++)
for (j=0; j<a.col; j++){
c.v[i][j]=a.v[i][j];
if (i<b.row)
c.v[i+a.row][j]=b.v[i][j];
}}


void AddCol(matrix a, matrix b, matrix c){
if ((a.row-b.row) || (a.row - c.row)){
puts("Error in the AddCol Function:\nThe rows of the matrices should be equal");
return;
}

int i,j;
for (i=0; i<a.row; i++)
for (j=0; j<a.col; j++){
c.v[i][j]=a.v[i][j];
if (j<b.col)
c.v[i][j+a.col]=b.v[i][j];
}}


//a=In
void I (matrix a){
if (a.row- a.col){
puts("Error in the I Function:\nThe matrix should be a square matrix");
return;
}

int i,j;
for (i=0; i<a.row; i++)
for (j=0; j<a.col; j++){
if (i==j) a.v[i][j]=1;
else a.v[i][j]=0;
}}


//a[i][j]=x for every i,j
void Cst (matrix a, double x){
int i,j;
for (i=0; i<a.row; i++)
for (j=0; j<a.col; j++)
a.v[i][j]=x;
}


void RandFill(matrix a, int x){
int i,j;
srand(time(NULL));
for (i=0; i<a.row; i++)
for (j=0; j<a.col; j++){
a.v[i][j]=rand()%x;
}}


//c=a+b
void Add (matrix a, matrix b, matrix c){
int i,j;
if ((a.row-b.row) || (a.col-b.col) || (a.row- c.row) || (a.col - c.col)){
puts("Error in the Add Function:\nThe matrices should have the same sizes");
return;
}

for (i=0; i<a.row; i++)
for (j=0; j<a.col; j++)
c.v[i][j]=a.v[i][j]+b.v[i][j];
}


//z = x*a (x in R)
void Scalar (matrix a, matrix z, double x){
if ((a.row-z.row) || (a.col- z.col)){
puts("Error in the Scalar Function:\nThe matrices should have the same sizes");
return;
}

int i,j;
for (i=0; i<a.row; i++)
for (j=0; j<a.col; j++)
z.v[i][j]=x*a.v[i][j];
}


//c = a*b (matrices)
void Product(matrix a, matrix b, matrix c){
int i,j,k;
if ((a.col-b.row) ||(a.row -c.row) || (b.col - c.col)){
puts("Error in the Product Function:\nThe inner sizes should be equal");
return;
}

for (i=0; i<a.row; i++)
for (j=0; j<b.col; j++){
c.v[i][j]=0;
for (k=0; k<b.row; k++)
c.v[i][j]+=a.v[i][k]*b.v[k][j];
}}


int Equal(matrix a, matrix z){
if ((a.row - z.row) || (a.col - z.col))
return 0;
int i,j;
for (i=0; i<a.row; i++)
for (j=0; j<a.col; j++)
if (a.v[i][j]-z.v[i][j]) return 0;
return 1;
}


void Transpose (matrix a, matrix z){
if ((a.row- z.row) || (a.col-z.col)){
puts("Error in the Transpose Function:\nThe matrices should have the same sizes");
return;
}

int i,j;
for (i=0; i<z.row; i++)
for (j=0; j<z.col; j++)
z.v[i][j]=a.v[j][i];
}



//Trace of a square matrix is the sum of the elements on the diagonal
double Trace (matrix a){
if (a.row-a.col){
puts("Error in Trace Function:\nThe matrix should be a square matrix");
return 0;
}

int i,s=0;
for (i=0; i<a.row; i++)
s+=a.v[i][i];
return s;
}


//Row echelon form
void Ref (matrix a, matrix z){
if ((a.row-z.row) || (a.col- z.col)){
puts("Error in the Ref Function:\nThe matrices should have the same sizes");
return;
}

int i,j,x,y,m=a.row,n=a.col;
double b,t;

for (i=0; i<m; i++)
for (j=0; j<n; j++)
z.v[i][j]=a.v[i][j];

i=0; j=0;
for (i=0,y=0; i<m; i++,y++){
if (z.v[y][i]<1e-4 && z.v[y][i]>-1e-4){
for (x=y; x<m; x++)
if (z.v[x][i]){
for (j=0; j<n; j++){
t=z.v[y][j];
z.v[y][j]=z.v[x][j];
z.v[x][j]=t;
}
break;
}
if (x==m){
y--;
continue;
}}
        
for (x=y+1; x<m; x++){
b=z.v[x][i]/z.v[y][i];
int c=b;
if (b-c)
for (j=y; j<n; j++)
z.v[x][j]=(z.v[y][i]*z.v[x][j])-b* (z.v[y][i]*z.v[y][j]);
else
for (j=y; j<n; j++)
z.v[x][j]=z.v[x][j]-b*z.v[y][j];
}}}


//Determinant
double Det(matrix a){
if (a.row - a.col){
puts("Error in the Det Function:\nThe matrix should be a square matrix");
return 0;
}

int i,j,k,z,x=1,l;
double s;
matrix b={a.row-1,a.row-1};
MatInit(&b);

if (a.row==1)
s = a.v[0][0];
if (a.row==2)
s = a.v[0][0]*a.v[1][1]-a.v[0][1]*a.v[1][0];

if (a.row>2){
s=0;
for (z=0; z<a.row; z++){
k=0; l=0;
for (i=1; i<a.row; i++, k++){
for (j=0; j<a.row; j++){
if(j==z) continue;
b.v[k][l]=a.v[i][j];
l++;
}
l=0;
}
s=s+x*a.v[0][z]*Det(b);
x=-x;
}}
MatFree(b);
return s;
}


void Inverse (matrix a, matrix z){
if ((a.row-z.row) || (a.col- z.col) || (z.row-z.col)){
puts("Error in the Inverse Function:\nThe matrices should have the same sizes and be square matrices");
return;
}

double x=Det(a),y;
printf("Det = %g\n",x);
if (!x){
puts("Error in the Inverse Function:\nThe matrix is not Invertable");
return;
}

if (a.row<2) z.v[0][0]=1/x;

else{
int i,j,k,l,t,r,n=a.row;
matrix b={n-1, n-1};
MatInit(&b);


for (i=0; i<n; i++){
for (j=0; j<n; j++){
for (k=0, t=0; k<n; k++){
if (k==i) continue;
for (l=0, r=0; l<n; l++){
if (l==j) continue;
b.v[t][r]=a.v[k][l];
r++;
}
t++;
}
if ((i+j)%2 )y=-x;
else y=x;
z.v[j][i]=Det(b)/y;
}}
MatFree(b);
}}
