#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <omp.h>
#include <unistd.h>
#include <arcContinuation.h>
#include "nodo-sarea.h"

extern void arcContinuation(/* input */
                  int dim, double *var0, double *b0, double init, double end, 
                  int action, double *preal, long numpreal, int *pint, long numpint,
		     /* output */
                  double *var1, double *b1, double *taptr, double *tbptr, int * infoptr);
extern void vectorcp(int dim, double *optr, double *dptr);
extern void rest(int dim, val_type *var, time_type t, val_type *res, 
                 double *preal, long prelem, int *pint, long pielem);

extern void jacobian(int dim, val_type *var, time_type t, val_type *jac,
                     double *preal, long prelem, int *pint, long pielem);

extern int irakurri_hazia (FILE * fitx, double *haziaptr,int alkop);

#ifdef LOGFILE
extern FILE *loga;
#endif
/**************************************************************
typedef struct elem
{
  int maila;
  double hf;  // helburu funtzioaren balioa
  double *vptr;  // x_i eta landa_i (denak bektore berean, dim elementu)
  double *betaptr;  // hiperplanoaren Beta_i
  //int *kontlist;  // hau maila bakoitzean zenbatgarren zeroa den jakiteko erabil daiteke
  struct elem *hptr;
} elem;

typedef struct
{
  elem *lehenaptr; // hedatu beharreko hurrengoa hau izango da
  elem *azkenaptr;  // lortutako azkena honen ondoren sartuko dut
  int dim;
  int landakop;
  int *pint;
  long pintkop;
  double *preal;
  long prealkop;
  double tol;
  int arkukozerokopmax;
  double hfjasotzeko;
  double hfjarraitzeko;
} elementu_zerrenda;
********************************************************/

//double hfjasotzeko;  // helburu funtzioak ze baliotik gora eduki behar duen hurrengo mailakoen zerrendan jasotzeko
//double hfjarraitzeko;  // helburu funtzioa balio honetatik behara pasatzen bada arkuaren jarripena eten.
//int arkuko;   // zenbat zero hartuko ditugun kontuan arkuaren alde bakoitzean.


void hasieratu_landak(double *vptr, elementu_zerrenda *zptr)
{
  int dim,i,j;
  double *jac;
  double t;
  int xkop;
  double *gradientea;

xkop = zptr->dim -zptr->landakop; //s31-> 16 baina s33-> 17
  #ifdef SOLVER_GSL
    /* GSL */
    //M>=N izan behar du!!!!
  //  nik M<N daukat, beraz, iraulitako matrizea deskonposa dezaket... 
  //  eta horrela lortuko nukeen A = U S Q^T ren ordez A^T = Q S U^T izango dira (?)
  #else
    char jobu,jobvt;
    int dimm, dimn, info, lda, ldu, ldvt, dimwork;
    double *S,*U,*VT,*work;
    /* LAPACK */
    dimm = xkop; // 16 edo 17
    dimn= zptr->landakop;
    lda = xkop;
    ldu = xkop;
    ldvt = zptr->landakop;
    dimwork = 10000;
    jobu = 'A';
    jobvt= 'A';
    if ((S = (double *)malloc ((ldu+1)*sizeof(double)))==NULL)  // Vector dim 16 edo 17
      {
      #ifdef LOGFILE
      fprintf(loga,"malloc problem! landetan ipiv\n");
      #endif
      return;
      }
    if ((VT = (double *)malloc (ldvt*ldvt*sizeof(double)))==NULL)  // Vector dim 17x17
      {
      #ifdef LOGFILE
      fprintf(loga,"malloc problem! landetan2 ipiv\n");
      #endif
      return;
      }
    if ((U = (double *)malloc ((ldu+1)*(ldu+1)*sizeof(double)))==NULL)  // Vector dim 16x16 edo 17x17
      {
      #ifdef LOGFILE
      fprintf(loga,"malloc problem! landetan3 ipiv\n");
      #endif
      return;
      }
    if ((work = (double *)malloc (dimwork*sizeof(double)))==NULL)  // Vector dim 10000
      {
      #ifdef LOGFILE
      fprintf(loga,"malloc problem! landetan4 ipiv\n");
      #endif
      return;
      }
  #endif
  jac = (double *)malloc(sizeof(double) * zptr->dim * zptr->dim);
  gradientea = (double *)malloc (sizeof(double) * zptr->landakop * xkop);
	// jakobiarra kalkulatuko dut, kontuan izan iraulitako jakobiarra kalkulatzen dugula
  jacobian(zptr->dim, vptr, 0.0, jac, zptr->preal, zptr->prealkop, zptr->pint, zptr->pintkop);
	//jakobiarretik goiko eskubiko kuadranteko 16 lerro hartu behar ditut (lehenengoa ez!)

  /* jakobiarreko indizeak horrela daude:


  jak = (df1/x1 ....df1/dxn
         ...
         dfn/dx1    dfn/dxn)

   s31 denean aldagaia 16 dira eta landak 17
   
  jak = (0   33  66 ... (16.lerroa) 495           528  561     azkena 1056
                                        ************************************       
         1   34  67                 496 *(gureak) 529  562            1057 *
         2   35                     497 *         530  563            1058 *
         ...                            *         531  564            1059 *
         ...                            *         ...  ...             ... *
         ...                            *(gureak) 544  577            1072 *
                                        ************************************
         ...                                      ...  ...             ...
         32  65  98                 511           560  593            1088 )
  */
  
  /*
    s33 denean aldagaiak 17 dira eta landak ere 17 dira
    
                                   17 zutabe
  jak = 0  ...  (17-1)*34   xkop*34  ...  (17+17-1)*34
                          **************************
        1                 * 17*34+1                *
                          *                        *      xkop lerro
        ...               *  ...                   *
       17                 * 17*34+17               *
                          **************************
       18
       ...
       33                   611                1155 ) 
                           
  */
  for (i=0; i<xkop; i++)  //aldagai adina lerro: s33 -> 17 baina s31 -> 16
    {
    for (j=0; j<zptr->landakop; j++)  // 17 zutabe
      {
      gradientea[j*xkop+i] = jac[(j+xkop)*zptr->dim + i+1];  // gradientea ere iraulita daukat.
      //printf("M[%d] = %lf ", j*16+i, gradientea[j*16+i]);
      }
    //printf("\n");
    }
  
  /*   A = U * SIGMA * transpose(V)

  subroutine dgesvd ( 	
		character  	JOBU,    ---> 'N'  (no interesting)
		character  	JOBVT,   ---> 'A'  (all in V^T)
		integer  	M,       ---> 16  (16 row)
		integer  	N,       ---> 17  (17 column)
		double *(M, N)  A,       ---> gradientea
		integer  	LDA,     ---> M
		double *( 16 )  S,       ---> diag[16]
		double *(M,M)   U,       ---> U[16*16]
		integer  	LDU,     ---> 16
		double *(N,N)  	VT,      ---> Q[17*17]
		integer  	LDVT,    ---> 17
		double *( >5*N) WORK,    ---> work[100]
		integer  	LWORK,   ---> 100
		integer  	INFO     ---> return value (if 0 OK)
	) 	
  */
  // niri landa = V.(0,0, ... 0, 1)^T kalkulatzea interesatzen zait.
  dgesvd_(&jobu, &jobvt,&dimm, &dimn, gradientea, &lda, S, U, &ldu, VT, &ldvt,work, &dimwork, &info);   
  /*
  Funtzioak V iraulia itzultzen du, VT da niri
  interesatzen zaidana, eta bere azken zutabean dauzkat landa balioak!!!!
  */
/*
printf("dgesvd info (0 behar luke) = %d\n",info);
for (i=0; i<17; i++)  //17 lerro
    {
    for (j=0; j<17; j++)  // 17 zutabe
      {
      printf("Q[%d] = %lg ", j*17+i, VT[j*17+i]);
      }
    printf("\n");
    }
*/
/*
 for (i = 0; i<17; i++)
    {
    printf("singular value (%d) = %.8lf ",i, gradientea[16*i+i]);
    }
*/
  for (i = 0; i<zptr->landakop; i++)
    {
    vptr[zptr->dim-zptr->landakop+i] = VT[i*17+16];   
    //printf("landa(%d) = %.8lg ",i, vptr[zptr->dim-zptr->landakop+i]);
    }
  //printf("\nlortu ditut!! orain, memoriak liberatzera....\n");
  free(jac); 
  //printf("jac lixto. orain gradientea liberatzera....\n");
  free(gradientea);
  #ifdef SOLVER_GSL

  #else
  free(U);
  free(VT);
  free(S);
  free(work);
  #endif
}

void esferaratu_elementua(double *eptr, int osagaikop, double radioa)
{
  int i;
  double binorma_s, batnorma_s;
  
  radioa = 3;
  /* esferaren baldintza:
     \sum_j^s x_i^2 = radioa^2   <---> \sum_j^{dim-landakop-1}2*x_i^2  + x_{dim-landakop}^2 = radioa^2
     
     esferara eramateko aldaketa:
     esf_i = radioa* x_i / ||x||
       non ||x|| = \sum_j^s x_i^2 
  */
  for(i=0,binorma_s=0.0,batnorma_s=0.0;i<osagaikop-1;i++)
      {
      binorma_s += 2*eptr[i]*eptr[i];//printf("batura akumulatua = %lf ",binorma_s);
      batnorma_s += 2*fabs(eptr[i]);
      }

  binorma_s += eptr[osagaikop-1]*eptr[osagaikop-1]; // s inparea denez azkena behin bakarrik
  batnorma_s += fabs(eptr[osagaikop-1]);                  // s inparea denez azkena behin bakarrik
  printf("binormaren karratua = %lf ",binorma_s);
  //binorma_s = sqrt(binorma_s);
  //printf("haziaren batnorma eta binorma= %lf eta %lf dira\n",batnorma_s,binorma_s);
 // for (i=0;i<osagaikop;i++) eptr[i] = radioa * eptr[i]/binorma_s;
}


/* lortu_elementua
   fitxategitik edo taula batetik haziaren datuak lortu ondoren 
   hazi horri dagokion landa bektorea kalkulatzen du.
   Horrez gain, Beta bektorea ere hasieratzen du.
* /
void lortu_elementua(double *haziaptr, elementu_zerrenda *zptr)
{
  int i;
  elem *eptr;
  double *vptr;
  double *betaptr;

  eptr = (elem *)malloc(sizeof(elem));
  vptr = (double *)malloc(zptr->dim*sizeof(double));
  betaptr = (double *)malloc(zptr->dim*sizeof(double));
  eptr->maila=0;
  eptr->vptr = vptr;
  eptr->betaptr = betaptr;
  //eptr->kontlist = (int *)0;   // ez dauka kontlist taularik (0 mailako elementua da)
  eptr->hptr = (struct elem *)0;  // ez dauka hurrengo elementurik! hau lehenengoa da
  zptr->lehenaptr = eptr;
  //printf("hasieran landak zerorekin hasieratuko ditut... lehenengo %d haziarenak eta %d landa\n",zptr->dim - zptr->landakop, zptr->landakop);
  for(i=0;i<(zptr->dim - zptr->landakop);i++) vptr[i]=haziaptr[i];
  for ( ;i < zptr->dim; i++) vptr[i]= 0.0;
  // haziak konsistentzia baldintza betetzen badu, hau da, ez balego esferan esferara pasa behar da
  //esferaratu_elementua(vptr,zptr->dim - zptr->landakop, 3.0);
  
  hasieratu_landak(&vptr[0],zptr);  // rest funtzioaren jakobiarrean dauzkat 
	 // f_iren gradientearen elementuak.
         // rest_2-tik 17ra bitartean \lambda_i \gradientea f_i(x) daukagu, beraz, jakobiarrean 
         // 2.lerrotik 17. lerrora bitartean lerro bakoitzeko azken 17 elementuak dira f_i bakoitzaren
         // x_j-rekiko deribatuak
  //  Orain Beta balioak zehaztu behar ditut...
  for (i = 0; i<zptr->dim; i++) betaptr[i] = 0.0;
  betaptr[zptr->dim - zptr->landakop] = 1.0;
}

**********/


void of_lortu(int alkop, double *xptr, double *hfptr)
{
int i;
double hf;

for (hf=0, i=0; i<alkop; i++)  // hf = sum_1^s x_i
    hf +=xptr[i];
hf -= xptr[alkop-1]/2.0;        // hf = sum_1^{s-1} x_i + x_s/2
*hfptr = hf*hf;
}





void irakurri_elementua(FILE *solfitx, elementu_zerrenda *zptr, int *mailaptr, int *ekintzaptr,int *gmugaptr, int *bmugaptr)
{
  int i;
  elem *eptr;
  double *vptr;
  double *betaptr;
  double nondikzer[4];  // elementuaren maila, 
                        // egin beharreko ekintza (-1 behera, 0 gora, n behera/gora buelta kopurua)
                        // goiko maila,
                        // beheko maila

  eptr = (elem *)malloc(sizeof(elem));
  vptr = (double *)malloc(zptr->dim*sizeof(double));
  betaptr = (double *)malloc(zptr->dim*sizeof(double));
  eptr->vptr = vptr;
  eptr->betaptr = betaptr;
  //eptr->kontlist = (int *)0;   // ez dauka kontlist taularik (0 mailako elementua da)
  eptr->hptr = (struct elem *)0;  // ez dauka hurrengo elementurik! hau lehenengoa da
  irakurri_hazia(solfitx,&(vptr[0]),zptr->dim); // fitxategitik irakurri hasierako elementua
  i=irakurri_hazia(solfitx,&(nondikzer[0]),4); // fitxategitik irakurri hasierako elementuaren maila,
                                               // elementuarekin egin beharreko ekintza (double bezala)
                                               // goiko muga ( ez badago ez da ezer pasatzen...)
                                               // beheko muga ( ez badago ez da ezer pasatzen...)
  if (i > 3)
      {
      i= nondikzer[3]; // beheko muga da hau
      *bmugaptr = i;
      i= nondikzer[2]; // goiko muga da hau
      *gmugaptr = i;
      }
    else
      {
      *gmugaptr = 11;
      *bmugaptr = 0;
      }
  i = nondikzer[0]; // maila int bezala behar dut, konbertsioa egin dezala!!!
  *mailaptr = i;
  eptr->maila=i;
  zptr->grafoaptr[i].lehenaptr = eptr;
  zptr->grafoaptr[i].azkenaptr = eptr;
  zptr->grafoaptr[i].goraptr = eptr;
  zptr->grafoaptr[i].beheraptr = eptr;
  i = nondikzer[1];  // ekintza: -1 behera, 0 gora, n behera eta gero gora n aldiz
  *ekintzaptr = i;
  of_lortu(zptr->dim-zptr->landakop,&(vptr[0]),&(eptr->hf));
  //  Orain Beta balioak zehaztu behar ditut...
  for (i = 0; i<zptr->dim; i++) betaptr[i] = 0.0;
  betaptr[zptr->dim - zptr->landakop] = 1.0;
  
}



void arku_jarraipena_egin(int dim, double *xptr, double *betaptr, double tinit,double tend, int ekintza, elementu_zerrenda *zptr, 
            double *xzeroptr, double *betazeroptr, double *newt0ptr, double *newtendptr, int *infoptr)
{
int birsaiatu;
double *lagptr;
double newtend;
int info;

info = 0;
birsaiatu = 0;
newtend = tinit;
//printf("jarraipenera... \n");
arcContinuation(dim, xptr, betaptr, tinit, tend, ekintza, 
                zptr->preal, zptr->prealkop, zptr->pint, zptr->pintkop,
		xzeroptr, betazeroptr, newt0ptr, &newtend, &info);
//printf("jarraipenetik\n");
if (((info == 10)/*||(info == 40)*/) && (newtend < tend)) birsaiatu = 1;
        else birsaiatu = 0;

while (birsaiatu == 1)
    {
    lagptr=xptr;
    xptr=xzeroptr;
    xzeroptr=lagptr;
    lagptr=betaptr;
    betaptr=betazeroptr;
    betazeroptr=lagptr;
    tinit = newtend;
    arcContinuation(dim, xptr, betaptr, tinit, tend, ekintza, 
                     zptr->preal, zptr->prealkop, zptr->pint, zptr->pintkop,
		     xzeroptr, betazeroptr, newt0ptr, &newtend, &info);
    //printf("    birjarraipena... info = %d, tinit = %lf\n",info,tinit);
    if (((info == 10) /*||(*infoptr == 40)*/) && (*newtendptr < tend)) birsaiatu = 1;
        else birsaiatu = 0;
    }
*newtendptr = newtend;
*infoptr = info;
//printf("    jarraipena egina: info = %d\n",info);
}


double norma(int dim,double *lptr)
{
int i;
double emaitza;

emaitza = 0.;
for (i=0; i<dim; i++)
    emaitza += lptr[i]*lptr[i];
emaitza = sqrt(emaitza);
return(emaitza);
}

void copy_aldagaien_karratua(int dim, double *xzeroptr,double *lehenaptr)
{
int i;

for (i=0; i<dim;i++) lehenaptr[i]=xzeroptr[i]*xzeroptr[i];
}


int karratuaren_berdina_da(int dim, double *xzeroptr,double *lehenaptr,double tol)
{
double n1;
int i;
double *lptr;

lptr = (double *)malloc(sizeof(double)*dim);

for (i=0; i<dim; i++) 
    lptr[i]=xzeroptr[i]*xzeroptr[i] - lehenaptr[i];
n1 = norma(dim,lptr);
free(lptr);
if(n1<tol) return(1);
    else return(0);
}


int berdinak_dira(int dim, double *xzeroptr,double *besteaptr,double tol)
{
double n1;
int i;
double diff;

for (i=0,diff=0.0; i<dim; i++) 
    {
    n1 = xzeroptr[i] - besteaptr[i];
    diff += (n1*n1);
    }
n1 = sqrt(diff);
if(n1<tol) 
        {
        //printf("errepikatutako bat!!!!(diff = %lg)\n",n1);
        return(1);
        }
    else return(0);
}


extern double kalkulatu_distantzia(double *vptr, double *sol);


// return: if errepikatuta dago --> 0
//         if  hf txarra --> -1
//         if berria da  -->  1 
int berria_eta_hedagarria_da(double *xzeroptr,double hf, int zeroaren_maila, elementu_zerrenda *zptr)
{
int berria_da;
elem *elagptr;
double dist;

berria_da = 1;
//hemen elagptr elemantuak maila berriko lenenengo elemantua dauka (baldin balego)
elagptr = zptr->grafoaptr[zeroaren_maila].lehenaptr;
for( ; (elagptr != 0) && berria_da; elagptr = elagptr->hptr)
    {
    if(berdinak_dira(zptr-> dim-zptr->landakop, xzeroptr,elagptr->vptr,zptr->tol)) 
        {
        berria_da = 0;
        //printf("errepikapen bat!!! \n");
        }
     /*   else
       {
       dist = kalkulatu_distantzia(xzeroptr,elagptr->vptr);
       if (dist < 0.001) printf("berria dela dio, baina, zeroaren eta aurretik jasotakoaren arteko distantzia = %.12lf \n",
                dist);
       } */
    }
    
if (berria_da) 
    {
    if (hf < zptr->hfjasotzeko) return(-1);
        else return(1);
    }
return (0);
}



void gehitu_zerrendan(elementu_zerrenda *zptr, int zeroaren_maila, double *xzeroptr,double *betazeroptr, double off,int zerokop)
{

elem *bptr;

bptr = (elem*)malloc(sizeof(elem)); // zerrendako elementu berria
bptr->maila= zeroaren_maila;   // maila 
bptr->hf = off;  // helburu funtzioaren balioa
bptr->vptr = xzeroptr;  // x_i eta landa_i (denak bektore berean, dim elementu)
bptr->betaptr = betazeroptr;  // hiperplanoaren Beta_i
//bptr->kontlist = (int *)malloc(sizeof(int)*(maila+1));
//for (i=0; i<maila; i++) bptr->kontlist[i] = zptr->lehenaptr->kontlist[i];
//bptr->kontlist[maila]=zeropkop;

bptr->hptr=(struct elem *)0;
if (zptr->grafoaptr[zeroaren_maila].lehenaptr ==0)
    {  // mailako lehenengo elemntua lortu dut
    zptr->grafoaptr[zeroaren_maila].lehenaptr=bptr;
    zptr->grafoaptr[zeroaren_maila].azkenaptr=bptr;
    }
  else // azkenaren ondoren sartu behar dut
    {
    zptr->grafoaptr[zeroaren_maila].azkenaptr->hptr = bptr;
    zptr->grafoaptr[zeroaren_maila].azkenaptr = bptr;
    }
// aurretik gora edo behera hedatu beharrekorik balego erakusle horiek ez dira aldatu behar
// baina ez badago gora edo behera hedatu beharrekorik, berri hau izango da hedatu beharrekoa
if (zptr->grafoaptr[zeroaren_maila].goraptr ==0)
    zptr->grafoaptr[zeroaren_maila].goraptr=bptr;
if (zptr->grafoaptr[zeroaren_maila].beheraptr ==0)
    zptr->grafoaptr[zeroaren_maila].beheraptr=bptr;
}



/* hedatu_elementua
   elementu jakin bat hedatzen du, zerrendako lehena. Gainera, badaki elementua hedatu egin behar dela.
   Lehenengo kandidatua aurkitu bezain pronto jaso egiten du, bide itxia den ala ez kontrolatzeko.
   */
void hedatu_elementua(int berriaren_maila, elem * eptr, elementu_zerrenda *zptr)
{
int i,j;
int itxia;
int zerokop;
int jarraitu;
int hedagarria;
int ekintza;
int direkzioa;
elem lehena; // bidea itxia den ala ez kontrolatzeko
double *xptr,*xzeroptr;
double *betaptr,*betazeroptr; 
double *lagptr;
double *lehenaber2ptr;
double hfzeroarena;
double tinit, tend, newt0, newtend;
int info;
  
// printf("%d mailako elementua hedatzen, malloc aurretik\n",eptr->maila);
itxia = 0;   // bidea, oraingoz, ez da bide itxia.
zerokop = 0;  // oraindik ez dut zerorik aurkitu...
xptr = malloc(sizeof(double)*zptr->dim);
  //printf("elementua hedatzen, vectorcp aurretik\n");
vectorcp(zptr->dim,eptr->vptr,xptr);
betaptr = malloc(sizeof(double)*zptr->dim);
  //printf("elementua hedatzen, 2. vectorcp aurretik\n");
vectorcp(zptr->dim, eptr->betaptr,betaptr);
xzeroptr = malloc(sizeof(double)*zptr->dim);
betazeroptr = malloc(sizeof(double)*zptr->dim);
  //printf("elementua hedatzen, malloc ondoren\n");
ekintza = 2; // zeroak bilatzeko eskatuko diot arcContinuation funtzioari.
tinit = 0.0;
tend = 40.0;
for (jarraitu = 1, direkzioa = 1, hedagarria = 1; 
     jarraitu; 
     jarraitu = (!itxia) && (hedagarria != 0) && (info == 111) && (zptr->hfjarraitzeko < hfzeroarena) && 
                ((zerokop * direkzioa) < zptr->arkukozerokopmax),
     //xptr berria azken zeroa da, eta zaharra hurrengoa lortzeko erabiliko dut.
     // betaptr-rekin berdin jokatuko dut
     lagptr=xptr,xptr=xzeroptr, xzeroptr=lagptr,lagptr=betaptr, betaptr=betazeroptr, betazeroptr=lagptr
    )
  {
  arku_jarraipena_egin(zptr->dim, xptr,betaptr,tinit,tend,ekintza,zptr, xzeroptr, betazeroptr, &newt0,&newtend,&info);
  of_lortu(zptr->dim - zptr->landakop, xzeroptr, &hfzeroarena);  
  //printf("jarraipenaren ondorengo of balioa %lf\n",hfzeroarena);
  if (info == 111)
    {
    zerokop +=direkzioa;
    //printf("zero bat (%d)!!!\n",zerokop);
    if (zerokop == 1)
      {
      lehenaber2ptr = (double *)malloc(zptr->dim*sizeof(double));
      copy_aldagaien_karratua(zptr->dim, xzeroptr,lehenaber2ptr);
      }
     else
      {
      itxia = karratuaren_berdina_da(zptr->dim,xzeroptr,lehenaber2ptr,zptr->tol);
      }
    hedagarria = berria_eta_hedagarria_da(xzeroptr,hfzeroarena,berriaren_maila,zptr); // -1 baldin hf txarra, 0 aurretik badago, 1 berria da
    if (hedagarria == 1 ) 
          {
          gehitu_zerrendan(zptr, berriaren_maila,xzeroptr,betazeroptr, hfzeroarena,zerokop); 
		//jaso egin dudanez, memoria hori ezin dut liberatu!!!
          lagptr = xzeroptr;
          xzeroptr = malloc(sizeof(double)*zptr->dim);
	  vectorcp(zptr->dim,lagptr,xzeroptr);
          lagptr = betazeroptr;
          betazeroptr = malloc(sizeof(double)*zptr->dim);
	  vectorcp(zptr->dim,lagptr,betazeroptr);
          //printf("zerrendan gehitu dut(%lf)\n",hfzeroarena);
          }
        //else printf("hedagarritasuna = %d\n",hedagarria); //zerokop -= direkzioa;
    /*
    if (!itxia) printf( "ez da bide itxia\n"); 
    if (info == 111) printf("bat aurkitu dudanez gehiago bilatzeko aukera daukat\n");
    printf("jarraitzeko hf = %lf: ",zptr->hfjarraitzeko);
    if(zptr->hfjarraitzeko < hfzeroarena) printf("zeroarena handiagoa (%lf), berz, segi bila! \n", hfzeroarena);
    printf("zerokopmax = %d, eta hau %d-garrena da\n",zptr->arkukozerokopmax, zerokop * direkzioa);
    */
    }
   else 
     {
     //printf("ez da 111: %d +\n",info);
     if (info == 20) /* arku jarraipen honek berak kurba itxia dela ikusi du */
         itxia =1;
     }
  }
// Orain beste aldera...
//if ( itxia || (hedagarria == 0)) printf("beste alderantz ez dut egin behar!!\n");
zerokop = 0;
vectorcp(zptr->dim, eptr->vptr,xptr);
vectorcp(zptr->dim, eptr->betaptr, betaptr);  // hasierako balioak hartu ditut
for (i=0; i<zptr->dim; i++) betaptr[i] = -betaptr[i];  // betari ikurra aldatu diot
for (jarraitu = !itxia && (hedagarria != 0), direkzioa = -1; 
     jarraitu; 
     jarraitu = (!itxia) && (hedagarria != 0) && (info == 111) && (zptr->hfjarraitzeko < hfzeroarena) && 
                ((zerokop * direkzioa) < zptr->arkukozerokopmax),
     //xptr berria azken zeroa da, eta zaharra hurrengoa lortzeko erabiliko dut.
     // betaptr-rekin berdin jokatuko dut
     lagptr=xptr,xptr=xzeroptr, xzeroptr=lagptr,lagptr=betaptr, betaptr=betazeroptr, betazeroptr=lagptr
    )
  {
  arku_jarraipena_egin(zptr->dim, xptr,betaptr,tinit,tend,ekintza,zptr, xzeroptr, betazeroptr, &newt0,&newtend,&info);
  of_lortu(zptr->dim -zptr->landakop, xzeroptr, &hfzeroarena); 
  //printf("jarraipenaren ondorengo of balioa %lf\n",hfzeroarena);
  if (info == 111)
    {
    zerokop +=direkzioa;
    //printf("zero bat (%d)!!!\n",zerokop);
    if ((zerokop*direkzioa) == 1)
      {
      lehenaber2ptr = (double *)malloc(zptr->dim*sizeof(double));
      copy_aldagaien_karratua(zptr->dim, xzeroptr,lehenaber2ptr);
      }
     else
      {
      itxia = karratuaren_berdina_da(zptr->dim,xzeroptr,lehenaber2ptr,zptr->tol);
      }
    hedagarria = berria_eta_hedagarria_da(xzeroptr,hfzeroarena,berriaren_maila,zptr); // -1 baldin hf txarra, 0 aurretik badago, 1 berria da
    if (hedagarria == 1 )  
          {
          gehitu_zerrendan(zptr,berriaren_maila, xzeroptr,betazeroptr, hfzeroarena,zerokop); 
		//jaso egin dudanez, memoria hori ezin dut liberatu!!!
          lagptr = xzeroptr;
          xzeroptr = malloc(sizeof(double)*zptr->dim);
	  vectorcp(zptr->dim,lagptr,xzeroptr);
          lagptr = betazeroptr;
          betazeroptr = malloc(sizeof(double)*zptr->dim);
	  vectorcp(zptr->dim,lagptr,betazeroptr);
          //printf("zerrendan gehitu dut(%lf)\n",hfzeroarena);
          }
       // else  zerokop -= direkzioa;
    }
   else
     {
     //printf("ez da 111: %d -\n",info);
     if (info == 20) /* arku jarraipen honek berak kurba itxia dela ikusi du */
         itxia =1;
     }
  }
free(xptr);
free(betaptr);
free(xzeroptr);
free(betazeroptr);
//printf("%d mailako elementua hedatuta\n",eptr->maila);
}


/* 
  mailari dagozkion parametroak zehazteko...
*/
void mailari_dagozkion_parametroak_zehaztu(int maila, elementu_zerrenda *zptr,int noranzkoa)
  {
  /* alfaren, betaren eta gammaren indizeak 1-etik hasten dira. 
     10 alfa daude, 11 beta eta 11 gamma.
     noranzkoa beherantz bada:
       goiko maila 11 da: oinarizko 5 baldintza + beste 11 betetzen ditu hasieran. 
         Gainera, 16. baldintza ere betetzen du, baina hori askatu eta bere ordez 
         \landa_16 zero noiz egiten den bilatu nahi dugu
       gamma_{maila+1} = 1 da beti, baina beste beta guztiak 0 dira.
       alfa_1 etik hasi eta alfa_maila bitartean 1 balioa dagokie, besteei 0. 
         Kontuan izan, hasieran maila 10 denez, alfa denak 1 direla. 
       beta guztiak 0 dira beti.
    noranzkoa gora bada:
       beta_{maila+1} = 1 da beti, baina beste beta guztiak 0 dira.
       alfa_1 etik hasi eta alfa_maila bitartean 1 balioa dagokie, besteei 0. 
          Kontuan izan, hasieran maila 0 denez, alfa denak 0 direla. 
       gamma guztiak 0 dira beti.
  */
  int i,azkenalfa;

  if (noranzkoa == 1) azkenalfa = maila;
      else azkenalfa = maila-1; // baldintza denak betetzen baditu maila 11 dauka!
                                // beherantz jo behar du, horretarako 10 alfak 
                                // 1 balioarekin jarri eta azken gamma ere =1
  for (i=0; i<azkenalfa; i++)
    {
    zptr->preal[i] = 1;  // alfa_i = 1 ,  i = 0..maila-1
    }
  for (;i<zptr->prealkop;i++) zptr->preal[i]=0.0;
  if (noranzkoa == -1)  // beherantz
      zptr->preal[10+11+maila-1]= 1;  // gamma_{i+1} = 1
    else  zptr->preal[10+maila]= 1;  // beta{i+1} = 1
 /* printf("parametro erreal kpurua = %ld: ",zptr->prealkop);
  for (i=0;i<10;i++) printf ("a(%d) = %lf ",i,zptr->preal[i]);
  for (i=0;i<11;i++) printf ("b(%d) = %lf ",i,zptr->preal[10+i]);
  for (i=0;i<11;i++) printf ("g(%d) = %lf ",i,zptr->preal[21+i]);*/
  }



void liberatu_elementua (elem *eptr)
{
//printf("aldagaiak liberatzera...\n");
free(eptr->vptr);
//printf("beta liberatzera...\n");
free(eptr->betaptr);
//printf("elementua bera liberatzera...\n");
free(eptr);

//printf("elementua liberatuta\n");
}

/* behera_hedatu_maila:
   maila jakin bateko kandidatuak hedatzen ditu banan bana.
  */
void hedatu_maila(int maila, elementu_zerrenda *zptr, int noranzkoa)
  {
  int i;
  elem * eptr;

  /* */
  if (noranzkoa == -1) eptr = zptr->grafoaptr[maila].beheraptr;
      else eptr = zptr->grafoaptr[maila].goraptr;
  for (i=0;
       (eptr !=0);
       eptr=eptr->hptr,i++);
  printf("%d mailan %d direkzioan hedatu beharreko elementu kopurua = %d da\n",maila,noranzkoa,i);
  /* */
  mailari_dagozkion_parametroak_zehaztu(maila, zptr,noranzkoa);
//probatu ea 0 itzultzen duen....
/* * /
  double *lptr;

  lptr = (double *)malloc(zptr->dim * sizeof (double));
  rest(zptr->dim, zptr->grafoaptr[maila].lehenaptr->vptr, 0.0, lptr, zptr->preal, zptr->prealkop, zptr->pint, zptr->pintkop);
  for (i=0; i< zptr->dim; i++) printf("Rest[%d] = %lg\n",i,lptr[i]);
  free(lptr);
 /* */
  if (noranzkoa == -1) eptr = zptr->grafoaptr[maila].beheraptr;
    else  eptr = zptr->grafoaptr[maila].goraptr;
  for ( ;
       ((eptr !=0) && (eptr->maila == maila));
       eptr = eptr->hptr  //tratatutako elementua saltatu
      )
    {
    //printf ("elementua hedatzera\n");
    hedatu_elementua(maila+noranzkoa,eptr,zptr);  // zein mailatako elementua lortuko duen esaten diot
    if (noranzkoa == -1) zptr->grafoaptr[maila].beheraptr = eptr->hptr;
        else zptr->grafoaptr[maila].goraptr = eptr->hptr;
    printf(".");fflush(stdout);
    }
  }

double kalkulatu_distantzia(double *vptr, double *sol)
{
double d,binorma_s;
int i;

 for(i=0,binorma_s=0.0;i<17;i++)
      {
      d = vptr[i]-sol[i];
      binorma_s += 2*d*d;//printf("batura akumulatua = %lf ",binorma_s);
      }
  d = vptr[17]-sol[17];
  binorma_s += d*d; // s inparea denez azkena behin bakarrik

return binorma_s;
}

void jaso_emaitzak_distantziarekin(elementu_zerrenda *zptr,char *base_dir, int multzoa, int threadzbkia,int hazia,double *sol)
{
elem *lptr;
int i,maila;
char fitxizena[80];
FILE *irteera;
int kop;
double dist;
for (maila = 0; maila < zptr->mailakop+1; maila++)
    {
    sprintf(fitxizena,"%setxea%dm%d-%d.txt",base_dir,maila,multzoa,threadzbkia);
    irteera=fopen(fitxizena,"w");
    //printf("%s fitxategian idaztera\n",fitxizena);
    for (kop = 0,lptr = zptr->grafoaptr[maila].lehenaptr;lptr != (elem *)0; lptr= lptr->hptr, kop++)
        {
        dist = kalkulatu_distantzia(lptr->vptr,sol);
        fprintf(irteera," soluzio bat: helburu funtzioa = %.12lf,dist^2 = %.12lf\n{",lptr->hf,dist);
        for (i=0; i<zptr->dim; fprintf(irteera,", "),i++)
            fprintf(irteera,"%.18lg",lptr->vptr[i]);
        fprintf(irteera,"}\n\n");
        }
    fclose(irteera);
    printf("%d nodo %d mailan\n",kop,maila);
    }
}


void fitxategitik_parametroak_hartu(elementu_zerrenda *zptr)
{
char fitxategia[80];
int i,res;
FILE *fitx;
struct stat stata;
char word[81];

strcpy(fitxategia,"./params.conf");
if (stat(fitxategia,&stata) == 0) 
    {
    fitx = fopen(fitxategia ,"r");
    for (i = fscanf(fitx,"%80s",word);
     	 i != EOF;
     	 i = fscanf(fitx,"%80s",word))
         /* no more than 80 chars, it stops in a space or EOL */
    	{
    	if (i == 1)  /* it has read the word */
      	   {
           if ((word[0] == '\0') || (word[0] == '#'))
           		res=fscanf(fitx,"%*[^\n]\n"); /* comment line */
            else 
               {
               if (strcmp(word,"tol")==0) 	//tolerance to decide when 2 roots are equal
                 {
                 res=fscanf(fitx,"%lf%*[^\n]\n",&(zptr->tol));
                 }
                else
              	 {
	       	 if (strcmp(word,"MinOF")==0)   // if OF of the root is less than 
						// MinOF the arc continuation proccess will stop
		       	res=fscanf(fitx,"%lf%*[^\n]\n",&(zptr->hfjarraitzeko)); 
                  else
		   {
		   if (strcmp(word,"MinExpandOF")==0)  // if OF is > than this value the root 
						     // will be added to the list
		      res=fscanf(fitx,"%lf%*[^\n]\n",&(zptr->hfjasotzeko)); 
                    else
                      {
                      if (strcmp(word,"numthreads")==0) /* the number of threads to be launched */
                          res=fscanf(fitx,"%d%*[^\n]\n",&(zptr->numthreads)); 
                        else
                          {
                          if (strcmp(word,"numseeds")==0) /* the number of seeds to be spreaded */
                              res=fscanf(fitx,"%d%*[^\n]\n",&(zptr->numseeds)); 
                            else
                              {
                              if (strcmp(word,"seedsongroup")==0) /* the number of seeds to be read and stored each time */
                                  res=fscanf(fitx,"%d%*[^\n]\n",&(zptr->seedsongroup));
                                else
                                  {
                                  if (strcmp(word,"arkukozerokopmax")==0) /* the number of seeds to be read and stored each time */
                                       res=fscanf(fitx,"%d%*[^\n]\n",&(zptr->arkukozerokopmax));
                                  } /* if (strcmp(word,"seedsongroup")==0)  */
                              } /* if (strcmp(word,"numseeds")==0)  */
                          } /* if (strcmp(word,"numthreads")==0)  */
                      } /* if (strcmp(word,"MinExpandOF")==0)  */
                   }  /* if (strcmp(word,"MinOF")==0)  */
                 }  /*if (strcmp(word,"tol")==0)   */
               }  /* if (word[0] == '\0' */
            } /* if (i==1) */
         } /* for*/
    fclose(fitx);
    }
  else /* configuration file no present! */	
	{ 
	}
}



void lortu_hazia(double *hazia,int zbkia,double *hazitaula, int alkop)
{
int i;
/* * /
hazia[0] =  0.25362259185648633;
hazia[1] = 0.25362259185648633;
hazia[2] = 1.1958841435506322;
hazia[3] = 0.25362259185648633;
hazia[4] = 0.25362259185648633;
hazia[5] = -1.465580112050182;
hazia[6] = 0.25362259185648633;
hazia[7] = 0.25362259185648633;
hazia[8] = 0.25362259185648633;
hazia[9] = -1.0406905558749862;
hazia[10] = 1.1958841435506324;
hazia[11] = 0.25362259185648633;
hazia[12] = 0.25362259185648633;
hazia[13] = 0.25362259185648633;
hazia[14] = -1.0406905558749855;
hazia[15] = 1.537842195022964;
/* */

for (i = 0; i <alkop; i++) hazia[i] = hazitaula[zbkia*alkop+i];
/* */

}

void hustu_zerrenda(elementu_zerrenda *zptr)
{
elem * eptr;
int maila;

for (maila=0; maila< zptr->mailakop+1;maila++)
  {
  for (
       eptr = zptr->grafoaptr[maila].lehenaptr;
       (eptr !=0);
       zptr->grafoaptr[maila].lehenaptr = eptr->hptr,liberatu_elementua(eptr), eptr = zptr->grafoaptr[maila].lehenaptr  // elementua kendu, eta liberatu memoria
      );
  }
}

// irakurri_hazia: return value: irakurritako zenbaki kopurua (alkop) or EOF

int irakurri_hazia (FILE * fitx, double *haziaptr,int alkop)
{
int i,j;
char kaka[10];
char *kakaptr;

kakaptr = &(kaka[0]);
j=fscanf(fitx, "%[{\n]%lf",kakaptr,haziaptr); // lehenengoa irakurri dut!
if (j<2)
    {
    printf("Ez dut haziaren hasiera ongi irakurri!!:%d elementu, %s da kaka eta %lf zbkia\n",j,kakaptr, haziaptr[0]);
    return(EOF);
    }
//printf("(%d)%s%lf",j,kakaptr,haziaptr[0]);
for (i=1; (i<alkop)&&(j==2)&&(j!=EOF); i++)
    {
    j=fscanf(fitx,"%[,\n]%lf",kakaptr,&(haziaptr[i]));
    //printf("(%d)%s%lf",j,kakaptr,haziaptr[i]);
    }
if ((j==2)&&(j!=EOF)) // alkop zenbaki irakurri ditu
    {
    j=fscanf(fitx, "%[},\n ]",kakaptr);
    //printf("(%d)%s",j,kakaptr);
    if (j!=1) j = EOF;
    }
return(i);
}



// irakurri_haziak:  return values: irakurritako hazi kopurua
//                                  infoptr for End Of File
//                                  taulaptr for seeds 
int irakurri_haziak(FILE *fitx, int multzokokop,int hedatzekokop,int alkop, int tratatuta, int *infoptr, double * taulaptr)
{
int i,j;

for (j=0, i = 1;
     (i != EOF)&&(j<multzokokop) && (j+tratatuta < hedatzekokop);
     j++) 
    {
    i = irakurri_hazia(fitx,&(taulaptr[j*alkop]),alkop);
    /*printf("%d hazia irakurrita\n",j)*/;
    }
if (i == EOF) j--;
//printf("%d hazi irakurrita\n",j);
*infoptr = i;
return(j);
}   

void grafoa_hasieratu(elementu_zerrenda *zptr, int mailakop)
{
int i;

// maila bakoitzeko mailako_nodoak moduko elementu bat edukiko dut
zptr->grafoaptr = (mailako_nodoak *)malloc((mailakop+1)* sizeof(mailako_nodoak));
zptr->mailakop = mailakop;
for (i = 0; i< mailakop+1; i++)
    {
    zptr->grafoaptr[i].maila = i;
    zptr->grafoaptr[i].lehenaptr = 0;  // ez dago elementurik
    zptr->grafoaptr[i].goraptr = 0;   // gora hedatu beharrekorik ez dago
    zptr->grafoaptr[i].beheraptr = 0; // behera hedatu beharrekorik ez dago
    zptr->grafoaptr[i].azkenaptr = 0; 
    }
}


void markotik_parametroak_hartu(elementu_zerrenda *nondikptr, elementu_zerrenda *zerrendaptr)
{
zerrendaptr->arkukozerokopmax =  nondikptr->arkukozerokopmax;
zerrendaptr->landakop =  nondikptr->landakop;
zerrendaptr->dim =  nondikptr->dim;
zerrendaptr->pintkop =  nondikptr->pintkop; 
zerrendaptr->prealkop =  nondikptr->prealkop;  
zerrendaptr->pint =  nondikptr->pint; 
zerrendaptr->hfjasotzeko =  nondikptr->hfjasotzeko;
zerrendaptr->hfjarraitzeko =  nondikptr->hfjarraitzeko;
zerrendaptr->tol =  nondikptr->tol;
}

int fcallkop;
int jcallkop;

int main(int argc, char *argv[])
{
int i,j,maila,mailamax;
int elmaila, gbekintza, gmuga, bmuga;
double *haziaptr;
int threadkop,nirezbkia,hazizbkia,hazikop;
int fileinfo;
int tratatuta;
int alkop;
elementu_zerrenda *zerrendaptr;
elementu_zerrenda markoa;
int multzokokop,multzozbkia;
int azkenetakoa;
char hazifitx_izena[80];
char base_dir[80];
FILE *hazienfitx;
struct stat stata;
double hasierakoa[18];
int gorabehera, gorabeherakop;


fcallkop = 0;
jcallkop=0;
setenv("OMP_NUM_THREADS","1",1); // lapack liburutegiari paralelizazioa ez erabiltzeko agindua!

if (argc >1) sprintf(base_dir,"banatuta%s/",argv[1]);
    else sprintf(base_dir,"./");
sprintf(hazifitx_izena,"%sabiapuntuak.txt",base_dir);
if (stat(hazifitx_izena,&stata) == 0) 
    {
    hazienfitx = fopen(hazifitx_izena,"r");
    }
  else
    {
    printf("ez dut hazien fitxategia aurkitu... Aio!\n");
    exit(0);
    }
    
mailamax = 11;
markoa.arkukozerokopmax = 10;
markoa.landakop = 17;
markoa.dim = 35;
alkop = markoa.dim - markoa.landakop;
markoa.pintkop =0; 
markoa.prealkop = 32;  // 10 alfa + 11 beta + 11 gamma
markoa.pint= (int *)0; 
markoa.hfjasotzeko = 2.1;
markoa.hfjarraitzeko = 0.2;
markoa.tol = 0.000001;
markoa.numthreads = 32;
markoa.numseeds = 10;
markoa.seedsongroup= 1000;
/* hazi bakarra */
fitxategitik_parametroak_hartu(&markoa);
multzokokop =markoa.seedsongroup; // multzokokop hazi irakurriko ditu
grafoa_hasieratu(&markoa, mailamax);
irakurri_elementua(hazienfitx,&markoa,&elmaila,&gbekintza,&gmuga, &bmuga); // soluzioak mailako baldintza denak betetzen ditu. 
	// 0 maila hazien maila da (bost baldintza betetzen dituztenena)
	// 11 maila soluzioen maila da (16 baldintzak betetzen dituztenena)
        // gbekintza: -1 behera bakarrik
        //             0 mailatik hasi eta gora
        //             n mailatik hasi eta behera eta gora n aldiz

for (i = 0; i < alkop ; i ++) hasierakoa[i] = markoa.grafoaptr[elmaila].lehenaptr->vptr[i];

zerrendaptr = &markoa;
zerrendaptr->preal = (double *)malloc(markoa.prealkop*sizeof(double));
/* baldintzak betetzen dituen ikusteko: */
  double *lptr;

  lptr = (double *)malloc(markoa.dim * sizeof (double));
  mailari_dagozkion_parametroak_zehaztu(elmaila, &markoa,-1); 
  rest(markoa.dim, markoa.grafoaptr[elmaila].lehenaptr->vptr, 0.0, lptr, markoa.preal, markoa.prealkop, markoa.pint, markoa.pintkop);
  for (i=0; i< markoa.dim; i++) printf("Rest[%d] = %lg\n",i,lptr[i]);
  printf("elementuaren maila %d da eta %d ekintza egin behar dut\n",elmaila,gbekintza);
  free(lptr);
 /* */
  /* * /   
  #pragma omp parallel private(haziaptr,nirezbkia,hazizbkia,zerrendaptr,maila) num_threads(markoa.numthreads)
    {
    threadkop= omp_get_num_threads();
    nirezbkia= omp_get_thread_num();
    printf("%d thread daude!!! %d\n",threadkop,nirezbkia);
    /* */
    threadkop = 1;
    nirezbkia = 0;
    hazizbkia  = 0;
    /* */
  
 // #pragma omp barrier
      /**************************** /
      // ea elementua ondo sortuta dagoen...
      mailari_dagozkion_parametroak_zehaztu(10, zerrendaptr); 
      double *lptr;
      lptr = (double *)malloc(zerrendaptr->dim * sizeof (double));
      rest(zerrendaptr->dim, zerrendaptr->lehenaptr->vptr, 0.0, lptr, zerrendaptr->preal, zerrendaptr->prealkop, zerrendaptr->pint, zerrendaptr->pintkop);
      for (i=0; i< zerrendaptr->dim; i++) printf("Rest[%d] = %lg\n",i,lptr[i]);
      if (hazizbkia == hazikop-1) 
            {
            printf("haziaren zenbakia: %d da \n",multzozbkia * multzokokop + hazizbkia);
            for (i=0; i<(zerrendaptr->dim-zerrendaptr->landakop); i++) printf("hazia[%d] = %lg\n",i,zerrendaptr->lehenaptr->vptr[i]);
            }
      free(lptr);
      // honaino da proba.
      /******************************/
      /***********/
    if (gbekintza > 0) gorabeherakop=gbekintza;
        else gorabeherakop = 1;
    for(gorabehera = 0;gorabehera < gorabeherakop;gorabehera++)
      {
      if (gbekintza == 0)  // ez du beheranzkorik egin behar
          maila = bmuga;
        else 
          {
          if (gbekintza == -1) maila = elmaila; // behera egin behar du, baina mailatik hasita eta behin bakarrik
            else maila = gmuga; // behera eta gora egin behar du gbekintzan esaten den adina aldiz
          }
      for (; (maila > bmuga) /*&& (zerrendaptr->grafoaptr[maila].lehenaptr != 0)*/; maila --)
        {
        //printf("%d maila hedatzera: %d hazia (%d threada): \n",maila, hazizbkia, nirezbkia); 
        hedatu_maila(maila,zerrendaptr,-1); // beherantz hedatu nahi dut
        
        /* */
        elem * eptr;
        for (i=0, eptr = zerrendaptr->grafoaptr[maila-1].lehenaptr;
            (eptr !=0);
            eptr=eptr->hptr,i++);
        //printf("    %d maila hedatuta: %d hazia (%d threada). Zerrendan %d elementu sartu dira \n",maila, hazizbkia, nirezbkia,i); 
        /* */
        }
     // Orain gorantz hedatuko ditut dauden guztiak...
     if (gbekintza == 0)  // ez du beheranzkorik egin eta gbmailatik hasita gora joan behar du
          maila = elmaila;
        else
          {
          if (gbekintza == -1) maila = gmuga;  // ez du goraznko biderik egin behar, behera egin du elmailatik hasita
            else maila = bmuga;
          }
     for (; (maila <gmuga) /* && (zerrendaptr->grafoaptr[maila].lehenaptr != 0)*/; maila ++)
        {
        //printf("%d maila hedatzera: %d hazia (%d threada): \n",maila, hazizbkia, nirezbkia); 
        hedatu_maila(maila,zerrendaptr,1); // gorantz hedatu nahi dut
        
        /* */
        elem * eptr;
        for (i=0, eptr = zerrendaptr->grafoaptr[maila+1].lehenaptr;
            (eptr !=0);
            eptr=eptr->hptr,i++);
        //printf("    %d maila hedatuta: %d hazia (%d threada). Zerrendan %d elementu sartu dira \n",maila, hazizbkia, nirezbkia,i); 
        /* */
        }
        printf("buelta osoa eman dut (%d)\n",gorabehera);
      } // buelta ematea bukatu du
      //printf("%d haziaren tratamendua bukatu da (%d threadak egina)\n",hazizbkia,nirezbkia);
      
        jaso_emaitzak_distantziarekin(zerrendaptr,base_dir,gbekintza,nirezbkia,hazizbkia,&(hasierakoa[0]));
        hustu_zerrenda(zerrendaptr);
    free(zerrendaptr->grafoaptr);
    free(zerrendaptr->preal);
    //free(zerrendaptr); zerrendaptr aldagaiak markoa erakusten du
    azkenetakoa = hazizbkia-threadkop;
    /* * /
    }  //paraleloko atala bukatu da
    /* */
    tratatuta += hazikop;
printf("funtzio dei kopurua = %d, jakobiarrari dei kopurua = %d\n",fcallkop, jcallkop);
fclose(hazienfitx);
azkenetakoa ++;
}



/* programa nagusia:
  hazia hartu, zerrenda sortu eta hazia hedatzeaz arduratzen da.
  defekturzko balioak dauzka hedapenerako
*/

/*
int main(int argc, char *argv[])
{
int i,j,maila,mailamax;
double hazia[16];


mailamax = 11;
elementu_zerrenda zerrenda;
double *lptr;

//lortu_hazia(fitxategia,&(hazia[0]) );

//{0.253623, 0.253623, 1.19588, 0.253623, 0.253623, -1.46558, 0.253623,
// 0.253623, 0.253623, -1.04069, 1.19588, 0.253623, 0.253623, 0.253623,
//-1.04069, 1.53784}

hazia[0] = 0.25362259185648633;
hazia[1] = 0.25362259185648633;
hazia[2] = 1.1958841435506322;
hazia[3] = 0.25362259185648633;
hazia[4] = 0.25362259185648633;
hazia[5] = -1.465580112050182;
hazia[6] = 0.25362259185648633;
hazia[7] = 0.25362259185648633;
hazia[8] = 0.25362259185648633;
hazia[9] = -1.0406905558749862;
hazia[10] = 1.1958841435506324;
hazia[11] = 0.25362259185648633;
hazia[12] = 0.25362259185648633;
hazia[13] = 0.25362259185648633;
hazia[14] = -1.0406905558749855;
hazia[15] = 1.537842195022964;
zerrenda.landakop = 17;
zerrenda.dim = 33;
zerrenda.pintkop =0; 
zerrenda.prealkop = 32;  // 10 alfa + 11 beta + 11 gamma
zerrenda.preal = (double *)malloc(zerrenda.prealkop*sizeof(double));
zerrenda.pint= (int *)0; 
zerrenda.hfjasotzeko = 0.5;
zerrenda.hfjarraitzeko = 0.2;
zerrenda.tol = 0.000000001;
fitxategitik_parametroak_hartu(&zerrenda);
   printf("parametroak: hfjarraitzeko %lg, hfjasotzeko %lf tol %lg\nElementuaren landak lortzera \n",
          zerrenda.hfjarraitzeko, zerrenda.hfjasotzeko, zerrenda.tol); 
lortu_elementua(&(hazia[0]), &zerrenda);
zerrenda.arkukozerokopmax = 10;
zerrenda.azkenaptr = zerrenda.lehenaptr;
for (maila = 0; maila <  mailamax; maila ++)
    {
    printf("maila hedatzera (%d): \n",maila); 
    hedatu_maila(maila,&zerrenda);
    }
if (zerrenda.lehenaptr != 0) jaso_emaitzak(&zerrenda,1,1,1);
printf("bukatu da\n");
}
*/



