#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int *iptr;
    int kop;
    } JintTau;
    

typedef struct {
    double *dptr;
    int kop;
    } JdTau;

typedef struct {
    int *iptr;
    int kop;
    int kop2;
    } JintTau2;
   

typedef struct {
    double *dptr;
    int kop;
    int kop2;
    } JdTau2; 

typedef struct {
    int *iptr;
    int kop;
    int kop2;
    int kop3;
    } JintTau3;
    

typedef struct {
    double *dptr;
    int kop;
    int kop2;
    int kop3;
    } JdTau3;

typedef struct {
    int maxdegree;//::Int64
    JintTau firstti;//::Array{Int64,1}
    JintTau tp;//::Array{Int64,1}
    JintTau degree;//::Array{Int64,1}
    JintTau dec1;//::Array{Int64,1}
    JintTau dec2;//::Array{Int64,1}
} TreeSet;


void pushD(JdTau *jt,double val)
{
double *bptr;
bptr = (double *)realloc(jt->dptr,(jt->kop+1)*sizeof(double));
bptr[jt->kop]=val;
jt->kop++;
jt->dptr=bptr;
}

void pushInt(JintTau *jt,int val)
{
int *bptr;
bptr = (int *)realloc(jt->iptr,(jt->kop+1)*sizeof(int));
bptr[jt->kop]=val;
jt->kop++;
jt->iptr=bptr;
}

int itest(int m, int l, int *degree)
    {
    if ((m==0) || (degree[m-1] < degree[l-1]))
        return 0;
      else
        if (degree[m-1]>degree[l-1])
            return 1;
          else
            return(m > l);
    }


int iseven(int i)
{
if ((i %2) == 0) return 1;
    else return 0;
}   
   
int isodd(int i)
{
if ((i %2) == 0) return 0;
    else return 1;
}   


TreeSet *ConstructInfTrees(int pmax,int (*test)(int))
{
int i,n,p,q,m,l;
TreeSet *solptr;
    JintTau firstti;// = zeros(Int64,pmax+1)
    JintTau tp;// = Array{Int64}(undef,0)
    JintTau degree;// = Array{Int64}(undef,0)
    JintTau dec1;// = Array{Int64}(undef,0)
    JintTau dec2;// = Array{Int64}(undef,0)
    /* hasieraketak*/
    firstti.iptr = (int *) malloc((pmax+1)*sizeof(int));
    for (i=0;i<(pmax+1);i++) firstti.iptr[i]=0;
    firstti.kop = pmax+1;
    tp.iptr = (int *)0;
    tp.kop = 0;
    degree.iptr = (int *)0;
    degree.kop = 0;
    dec1.iptr = (int *)0;
    dec1.kop = 0;
    dec2.iptr = (int *)0;
    dec2.kop = 0;
 
solptr= (TreeSet *) malloc(sizeof(TreeSet));   
    n = 0;
    firstti.iptr[0] = 1;
    
    if (test(1))
        {
        pushInt(&tp,1);
        pushInt(&degree,1);
        pushInt(&dec1,1);
        pushInt(&dec2,0);
        n +=1;
        }
    for (p=2; p<=pmax; p++) // in 2:pmax
        {
        firstti.iptr[p-1] = n+1;
        if (test(p))
            {
            n += 1;
            pushInt(&tp,p);
            pushInt(&degree,1);
            pushInt(&dec1,n);
            pushInt(&dec2,0);
            }
        for (q=1; q <p; q++) // in 1:p-1
            {
            for (m=firstti.iptr[q-1]; m<firstti.iptr[q]; m++) //in firstti[q]:firstti[q+1]-1, 
                {
                for ( l = firstti.iptr[p-q-1]; l < firstti.iptr[p-q]; l++) // l in firstti[p-q]:firstti[p-q+1]-1
                    {
                    if (itest(m,l,degree.iptr) &&  !itest(dec2.iptr[m-1],l,degree.iptr))
                        {
                        n += 1;
                        pushInt(&tp,tp.iptr[m-1]);
                        pushInt(&degree,degree.iptr[m-1]+degree.iptr[l-1]);
                        pushInt(&dec1,m);
                        pushInt(&dec2,l);
                        }
                    }
                }
            }   
        }
    firstti.iptr[pmax] = n+1;
    
    // soluzioaren egitura betetzeko
    solptr->maxdegree=pmax;
    solptr->firstti.kop = firstti.kop;
    solptr->firstti.iptr = firstti.iptr;
    solptr->tp.kop = tp.kop;
    solptr->tp.iptr = tp.iptr;
    solptr->degree.kop = degree.kop;
    solptr->degree.iptr = degree.iptr;
    solptr->dec1.kop = dec1.kop;
    solptr->dec1.iptr = dec1.iptr;
    solptr->dec2.kop = dec2.kop;
    solptr->dec2.iptr = dec2.iptr;
    return solptr;//TreeSet(pmax,firstti, tp, degree, dec1, dec2)
}
  

typedef struct {
    JintTau symind;//::Vector{Int64}
    JintTau treeindices;//::Vector{Int64}
    TreeSet *treesetptr;//::TreeSet
    JdTau2 coeffpow;//::Array{ctype,2}
    JdTau2 phi;//::Array{ctype,2}
    JdTau2 phip;//::Array{ctype,2}
    JdTau2 hcoeffpow;//::Array{ctype,2}
    JdTau2 hphi;//::Array{ctype,2}
    JdTau2 hphip;//::Array{ctype,2}
    JdTau3 Jcoeffpow;//::Array{ctype,3}
    JdTau3 Jphi;//::Array{ctype,3}
    JdTau3 Jphip;//::Array{ctype,3}
    JdTau3 Jhcoeffpow;//::Array{ctype,3}
    JdTau3 Jhphi;//::Array{ctype,3}
    JdTau3 Jhphip;//::Array{ctype,3}
    JdTau Jaux;//::Array{ctype,1}
    JdTau Jhaux;//::Array{ctype,1}
} OrderConditionsCache;



OrderConditionsCache *OrderConditionsSymInit(int maxorder, JdTau *coeffs)
{
OrderConditionsCache *solptr;
int m,s,n, i, j, maxord,ntrees;
solptr = (OrderConditionsCache *)malloc(sizeof(OrderConditionsCache));

    m = coeffs->kop;//length(coeffs)
    s = 2*m-1;
    //ctype = eltype(coeffs)
    maxord = iseven(maxorder) ? maxorder-1 : maxorder;

    solptr->symind.kop = 2*m-1;
    solptr->symind.iptr=(int *)malloc((2*m-1)*sizeof(int));
    for (i=0; i<m; i++) solptr->symind.iptr[i]= i+1;
    for (j=m-1; j>0;i++, j--) solptr->symind.iptr[i]= j;
    solptr->treesetptr = ConstructInfTrees(maxord,isodd);
    //@unpack firstti = treeset
    ntrees = solptr->treesetptr->firstti.iptr[maxord]-1;//firstti[maxorder+1]-1 ANDER: maxord ala maxorder??
    solptr->treeindices.kop = 1;
    solptr->treeindices.iptr=(int *)malloc(sizeof(int));
    solptr->treeindices.iptr[0] = 1;
    for (n=3; n<=maxord; n+=2)//n in 3:2:maxord
        {
        //push!(treeindices, (firstti[n]:firstti[n+1]-1)...)
        for (i = solptr->treesetptr->firstti.iptr[n-1];i<=solptr->treesetptr->firstti.iptr[n]-1; i++)
            {
            pushInt(&(solptr->treeindices),i);
            }
        }
     
    solptr->coeffpow.kop=m;
    solptr->coeffpow.kop2=maxord;
    solptr->coeffpow.dptr = (double *)malloc((m*maxord)*sizeof(double)); // Array{ctype}(undef,m,maxorder)
    solptr->phi.kop= s+1;
    solptr->phi.kop2= ntrees;
    solptr->phi.dptr= (double *)malloc(((s+1)*ntrees)*sizeof(double)); //Array{ctype}(undef,s+1,ntrees)
    solptr->phip.kop= s;
    solptr->phip.kop2= ntrees;
    solptr->phip.dptr= (double *)malloc((s*ntrees)*sizeof(double)); // Array{ctype}(undef,s,ntrees)
    solptr->hcoeffpow.kop= m;
    solptr->hcoeffpow.kop2= maxord;
    solptr->hcoeffpow.dptr= (double *)malloc((m*maxord)*sizeof(double)); //hcoeffpow = Array{ctype}(undef,m,maxorder)
    solptr->hphi.kop= s+1;
    solptr->hphi.kop2= ntrees;
    solptr->hphi.dptr= (double *)malloc(((s+1)*ntrees)*sizeof(double)); // hphi = Array{ctype}(undef,s+1,ntrees)
    solptr->hphip.kop= s;
    solptr->hphip.kop2= ntrees;
    solptr->hphip.dptr= (double *)malloc((s*ntrees)*sizeof(double)); // hphip = Array{ctype}(undef,s,ntrees)
    solptr->Jcoeffpow.kop= m;
    solptr->Jcoeffpow.kop2= maxord;
    solptr->Jcoeffpow.kop3= m+s;
    solptr->Jcoeffpow.dptr= (double *)malloc((m*maxord*(m+s))*sizeof(double)); 
    for (i=0;i<m*maxord*(m+s);i++) solptr->Jcoeffpow.dptr[i] =0;  // Jcoeffpow = zeros(ctype,m,maxorder,m+s)
    solptr->Jphi.kop= s+1;
    solptr->Jphi.kop2= ntrees;
    solptr->Jphi.kop3= m+s;
    solptr->Jphi.dptr= (double *)malloc(((s+1)*ntrees*(m+s))*sizeof(double));
    for (i=0;i<(s+1)*ntrees*(m+s);i++) solptr->Jphi.dptr[i] =0; //   Jphi = zeros(ctype,s+1,ntrees,m+s)
    solptr->Jphip.kop= s;
    solptr->Jphip.kop2= ntrees;
    solptr->Jphip.kop3= m+s;
    solptr->Jphip.dptr= (double *)malloc((s*ntrees*(m+s))*sizeof(double)); 
    for (i=0;i<s*ntrees*(m+s);i++) solptr->Jphip.dptr[i] =0;//    Jphip = zeros(ctype,s,ntrees,m+s)
    solptr->Jhcoeffpow.kop= m;
    solptr->Jhcoeffpow.kop2= maxord;
    solptr->Jhcoeffpow.kop3= m+s;
    solptr->Jhcoeffpow.dptr= (double *)malloc((m*maxord*(m+s))*sizeof(double));
    for (i=0;i<m*maxord*(m+s);i++) solptr->Jhcoeffpow.dptr[i] =0; //    Jhcoeffpow = zeros(ctype,m,maxorder,m+s)
    solptr->Jhphi.kop= s+1;
    solptr->Jhphi.kop2= ntrees;
    solptr->Jhphi.kop3= m+s;
    solptr->Jhphi.dptr= (double *)malloc(((s+1)*ntrees*(m+s))*sizeof(double)); 
    for (i=0;i<(s+1)*ntrees*(m+s);i++) solptr->Jhphi.dptr[i] =0;//    Jhphi = zeros(ctype,s+1,ntrees,m+s)
    solptr->Jhphip.kop= s;
    solptr->Jhphip.kop2= ntrees;
    solptr->Jhphip.kop3= m+s;
    solptr->Jhphip.dptr= (double *)malloc((s*ntrees*(m+s))*sizeof(double)); 
    for (i=0;i<s*ntrees*(m+s);i++) solptr->Jhphip.dptr[i] =0;//     Jhphip = zeros(ctype,s,ntrees,m+s)
    solptr->Jaux.kop= m+s;
    solptr->Jaux.dptr= (double *)malloc((m+s)*sizeof(double)); 
    for (i=0;i<(m+s);i++) solptr->Jaux.dptr[i] =0;//     Jaux = zeros(ctype,m+s)
    solptr->Jhaux.kop= m+s;
    solptr->Jhaux.dptr= (double *)malloc((m+s)*sizeof(double)); 
    for (i=0;i<(m+s);i++) solptr->Jhaux.dptr[i] =0;//    Jhaux = zeros(ctype,m+s)
    return(solptr); //OrderConditionsCache(symind, treeindices, treeset, coeffpow, 
                    //              phi, phip, hcoeffpow, 
                    //              hphi, hphip, Jcoeffpow, Jphi, Jphip, Jhcoeffpow, Jhphi, Jhphip, Jaux, Jhaux)
}



/*
function OrderConditionsSymLambda(maxorder::Int64, cache::OrderConditionsCache, coeffs::AbstractVector, lambda::AbstractVector)
    ncoeffs = length(coeffs)
    s = 2*ncoeffs-1
    @unpack symind, treeindices, treeset, coeffpow, phi, phip = cache
    @unpack firstti, tp, dec1, dec2 = treeset
    noc = length(treeindices)
    length(lambda)<noc ? error("length(lambda) < ", noc) : nothing
    for i in 1:ncoeffs
        coeffpow[i,1] = coeffs[i]
        coeffpow[i,2] = coeffpow[i,1]^2
    end
    for k in 3:2:maxorder
        for i in 1:ncoeffs
            coeffpow[i,k] = coeffpow[i,k-2]*coeffpow[i,2]
        end
    end     
    ntrees = firstti[maxorder+1]-1
    for n in 1:firstti[maxorder+1]-1
        m = dec1[n]
        l = dec2[n]
        k = tp[n]
        if l==0
            for j in 1:s
                phip[j,n] = coeffpow[symind[j],k]/2
            end
        else
            for j in 1:s
                phip[j,n] = phip[j,m]*phi[j,l]
            end
        end
        phi[1,n] = phip[1,n]
        for j in 2:s
            phi[j,n] = phi[j-1,n] + phip[j-1,n] + phip[j,n]
        end
    end
    for n in treeindices
        phi[s+1,n] = phi[s,n] + phip[s,n] 
    end
    phi[s+1,1] = sum(coeffpow[i,2] for i=1:ncoeffs-1) + coeffpow[ncoeffs,2]/2 - 9
    oclambda = (sum(coeffs[i] for i=1:ncoeffs-1) + coeffs[ncoeffs]/2)^2
    for i in eachindex(treeindices)
        n = treeindices[i]    
        oclambda += phi[s+1,n]*lambda[i]  
    end
    return oclambda
end
*/

double OrderConditionsSymLambdaGrad(JdTau GL/*::AbstractVector*/, int maxorder, OrderConditionsCache *cache, JdTau coeffs/*::AbstractVector*/, JdTau lambda/*::AbstractVector*/)
{
int i,j,k, ncoeffs,s,noc,ntrees,m,l,n;
double aux,oclambda,haux;
JintTau *symind,*treeindices;//::Vector{Int64}
TreeSet *treesetptr;//::TreeSet
JdTau2 *coeffpow, *phi, *phip, *hcoeffpow, *hphi, *hphip;//::Array{ctype,2}
JintTau *firstti, *tp, *dec1, *dec2;


if (iseven(maxorder)) maxorder --; // hau nik jarritakoa da.
    ncoeffs = coeffs.kop; //length(coeffs)
    s = 2*ncoeffs-1;
    //    @unpack symind, treeindices, treeset, 
    //        coeffpow, phi, phip, hcoeffpow, hphi, hphip = cache
    symind = &(cache->symind);
    treeindices = &(cache->treeindices);
    treesetptr = cache->treesetptr;
    coeffpow = &(cache->coeffpow);
    phi = &(cache->phi);
    phip = &(cache->phip);
    hcoeffpow = &(cache->hcoeffpow);
    hphi = &(cache->hphi);
    hphip = &(cache->hphip);
    //    @unpack firstti, tp, dec1, dec2 = treeset
    firstti = &(treesetptr->firstti);
    tp = &(treesetptr->tp);
    dec1 = &(treesetptr->dec1);
    dec2 = &(treesetptr->dec2);

//printf("unpakeatuta!\n");
    noc = treeindices->kop; //length(treeindices)
    if (lambda.kop < noc)  //length(lambda)<noc ? error("length(lambda) < ", noc) : nothing
        {
        printf("error: length(lambda) < %d\n",noc);
        exit(0);
        }
    if (GL.kop < (noc+ncoeffs)) //length(GL)<(noc+ncoeffs) ? error("length(GL) < ", noc+coeffs) : nothing
        {
        printf("error: length(GL) < %d\n",noc+ncoeffs);
        exit(0);
        }
    for (i=0;i<ncoeffs;i++) // in 1:ncoeffs
        {
        coeffpow->dptr[i*coeffpow->kop2+0] = coeffs.dptr[i];//coeffpow[i,1] = coeffs[i]
        coeffpow->dptr[i*coeffpow->kop2+1] = coeffs.dptr[i]*coeffs.dptr[i]; //coeffpow[i,2] = coeffpow[i,1]^2
        }
    for (k=2; k<maxorder;k+=2) // in 3:2:maxorder
        {
        for (i=0; i<ncoeffs; i++) //in 1:ncoeffs
            {
            coeffpow->dptr[i*coeffpow->kop2+k] = 
                   coeffpow->dptr[i*coeffpow->kop2+(k-2)]*
                   coeffpow->dptr[i*coeffpow->kop2+1];  //coeffpow[i,k] = coeffpow[i,k-2]*coeffpow[i,2]
            }
        }
    
    ntrees = firstti->iptr[maxorder]-1;//firstti[maxorder+1]-1
    for (n=0; n < ntrees; n++) //in 1:firstti[maxorder+1]-1
        {
        m = dec1->iptr[n]-1; // hauek zenbakiak dira, -1 eginaz indize bihurtu ditut!!
        l = dec2->iptr[n]-1; // hauek zenbakiak dira, -1 eginaz indize bihurtu ditut!!
        k = tp->iptr[n]-1;   // hauek zenbakiak dira, -1 eginaz indize bihurtu ditut!!
        if (l == -1) //l==0  bat kendu diot!
            {
            for (j=0;j<s;j++) //in 1:s
                {
                phip->dptr[phip->kop2*j+n] = coeffpow->dptr[(symind->iptr[j]-1)*coeffpow->kop2+k]/2.0; //phip[j,n] = coeffpow[symind[j],k]/2
                }
            }
        else
            {
            for (j=0;j<s;j++) // j in 1:s
                {
                phip->dptr[phip->kop2*j+n] = phip->dptr[phip->kop2*j+m]*phi->dptr[phi->kop2*j+l]; //phip[j,n] = phip[j,m]*phi[j,l]
                
                /*if (n==2) 
                    printf("%d, %d, %d  --> phip[j,n]= %lg * %lg\n",j,n,l, phip->dptr[phip->kop2*j+n],phi->dptr[phi->kop2*j+l] );
                */
                }
            }
        phi->dptr[n] = phip->dptr[n]; //phi[1,n] = phip[1,n]
        for (j=1;j<s;j++) //j in 2:s
            {
            phi->dptr[phi->kop2*j+n] = phi->dptr[phi->kop2*(j-1)+n] +
                                      phip->dptr[phip->kop2*(j-1)+n] +
                                      phip->dptr[phip->kop2*j+n];  //phi[j,n] = phi[j-1,n] + phip[j-1,n] + phip[j,n]
            }
        }

    //  ocnorm2 = zero(coeffs[1])    aldagai hau ez da erabiltzen!!!
    for (i=1;i<noc;i++)  //i in 2:noc
        {
        n = treeindices->iptr[i]-1; // treeindices[i]
        phi->dptr[phi->kop2*s+n] =  phi->dptr[phi->kop2*(s-1)+n] + 
                                     phip->dptr[phip->kop2*(s-1)+n];   // phi[s+1,n] = phi[s,n] + phip[s,n] 
        }

    //phi[s+1,1] = sum(coeffpow[i,2] for i=1:ncoeffs-1) + coeffpow[ncoeffs,2]/2 - 9
    phi->dptr[phi->kop2*s] = -9.0;
    for (i = 0; i< ncoeffs-1; i++)
        phi->dptr[phi->kop2*s] +=coeffpow->dptr[coeffpow->kop2*i+1];
    phi->dptr[phi->kop2*s] += (coeffpow->dptr[coeffpow->kop2*(ncoeffs-1)+1]/2.0);
        
    // aux = sum(coeffpow[i,1] for i=1:ncoeffs-1) + coeffpow[ncoeffs,1]/2  
    aux = 0.0;
    for (i = 0; i< ncoeffs-1; i++)
        aux +=coeffpow->dptr[coeffpow->kop2*i];
    aux += (coeffpow->dptr[coeffpow->kop2*(ncoeffs-1)]/2.0);
                                
    oclambda = aux*aux;
    for (i=0;i<noc;i++) //i in 1:noc
        {
        n = treeindices->iptr[i]-1; //treeindices[i]    
        oclambda +=  phi->dptr[phi->kop2*s+n] * lambda.dptr[i];   // oclambda += phi[s+1,n]*lambda[i]  
        }
    //printf("oclambda kalkulatuta! = %.12lg\n",oclambda);
    // Computation of the gradient with reverse mode algorithmic differentiation
    for (i=0; i < hcoeffpow->kop*hcoeffpow->kop2; i++) hcoeffpow->dptr[i] = 0.0;  //hcoeffpow .= 0
    for (i=0; i < hphi->kop*hphi->kop2; i++) hphi->dptr[i] = 0.0;  // hphi .= 0
    for (i=0; i < hphip->kop*hphip->kop2; i++) hphip->dptr[i] = 0.0;  // hphip .= 0
    for (i=0; i < GL.kop; i++) GL.dptr[i] = 0.0;  // GL .= 0
    haux = 0.0;
    for (i=0;i<noc;i++)  //i in 1:noc
        {
        n = treeindices->iptr[i] - 1; // treeindices[i] 
        GL.dptr[ncoeffs+i] += phi->dptr[phi->kop2*s+n];  //GL[ncoeffs+i] += phi[s+1,n]
        hphi->dptr[hphi->kop2*s+n] += lambda.dptr[i]; //hphi[s+1,n] += lambda[i]   
        }

    // Gradient of aux^2                                           
    haux += 2*aux;
    for (i=0;i<(ncoeffs-1); i++) //i in 1:ncoeffs-1 
        hcoeffpow->dptr[hcoeffpow->kop2*i] += haux;  //hcoeffpow[i,1] += haux
    hcoeffpow->dptr[hcoeffpow->kop2*(ncoeffs-1)] += haux/2.0;  // hcoeffpow[ncoeffs,1] += haux/2   # end of gradient of aux^2
    for (i=0;i<ncoeffs-1; i++) // i in 1:ncoeffs-1
        hcoeffpow->dptr[hcoeffpow->kop2*i+1] +=  hphi->dptr[hphi->kop2*s];  // hcoeffpow[i,2] += hphi[s+1,1]
    hcoeffpow->dptr[hcoeffpow->kop2*(ncoeffs-1)+1] +=  hphi->dptr[hphi->kop2*s]/2.0; //  hcoeffpow[ncoeffs,2] += hphi[s+1,1]/2
    for (i=1; i<noc;i++) //  i in 2:noc
        {
        n = treeindices->iptr[i]-1;  // treeindices[i]
        hphi->dptr[hphi->kop2*(s-1)+n] += hphi->dptr[hphi->kop2*s+n];  // hphi[s,n] += hphi[s+1,n]
        hphip->dptr[hphip->kop2*(s-1)+n] += hphi->dptr[hphi->kop2*s+n];  // hphip[s,n] += hphi[s+1,n]
        }

    for (n = firstti->iptr[maxorder]-2; n>=0; n--) //  n in firstti[maxorder+1]-1:-1:1
        {
        m = dec1->iptr[n]-1; // dec1[n]
        l = dec2->iptr[n]-1; // dec2[n]
        k = tp->iptr[n]-1; // tp[n]
        for (j=s-1; j>0; j--) // j in s:-1:2
            {
            hphi->dptr[hphi->kop2*(j-1)+n] +=  hphi->dptr[hphi->kop2*j+n]; // hphi[j-1,n] += hphi[j,n]
            hphip->dptr[hphip->kop2*(j-1)+n] +=  hphi->dptr[hphi->kop2*j+n]; // hphip[j-1,n] += hphi[j,n]
            hphip->dptr[hphip->kop2*j+n] +=  hphi->dptr[hphi->kop2*j+n]; // hphip[j,n] += hphi[j,n]
            }
        hphip->dptr[n] +=  hphi->dptr[n]; // hphip[1,n] += hphi[1,n]
        if (l==-1)  // l==0
            {
            for (j=0;j<s; j++) // j in 1:s
                {
                hcoeffpow->dptr[hcoeffpow->kop2*(symind->iptr[j]-1)+k] += 
                              hphip->dptr[hphip->kop2*j+n]/2.0; // hcoeffpow[symind[j],k] += hphip[j,n]/2
                }
            }
        else
            {
            for (j=s-1; j>=0;j--)  // j in s:-1:1
                {
                hphip->dptr[hphip->kop2*j+m] += hphip->dptr[hphip->kop2*j+n] *
                                              phi->dptr[phi->kop2*j+l];  // hphip[j,m] += hphip[j,n]*phi[j,l]
                hphi->dptr[hphi->kop2*j+l] += hphip->dptr[hphip->kop2*j+n] *
                                              phip->dptr[phip->kop2*j+m]; // hphi[j,l] += hphip[j,n]*phip[j,m]
                }
            }
        }
    for (k=maxorder-1; k>1; k -= 2)  // k in maxorder:-2:3
        {
        for (i=0; i<ncoeffs; i++) // i in 1:ncoeffs
            {
            hcoeffpow->dptr[hcoeffpow->kop2*i+k-2] += hcoeffpow->dptr[hcoeffpow->kop2*i+k] *
                                                    coeffpow->dptr[coeffpow->kop2*i+1]; // hcoeffpow[i,k-2] += hcoeffpow[i,k]*coeffpow[i,2]
            hcoeffpow->dptr[hcoeffpow->kop2*i+1] += hcoeffpow->dptr[hcoeffpow->kop2*i+k] *
                                                    coeffpow->dptr[coeffpow->kop2*i+k-2]; // hcoeffpow[i,2] += hcoeffpow[i,k]*coeffpow[i,k-2]
            }
        }

    for (i=0; i<ncoeffs; i++)  //  i in 1:ncoeffs
        {
        hcoeffpow->dptr[hcoeffpow->kop2*i] += 2.0*hcoeffpow->dptr[hcoeffpow->kop2*i+1] *
                                                    coeffpow->dptr[coeffpow->kop2*i];  // hcoeffpow[i,1] += 2*coeffpow[i,1]*hcoeffpow[i,2]
        GL.dptr[i] = hcoeffpow->dptr[hcoeffpow->kop2*i];  // GL[i] = hcoeffpow[i,1]
        }
    return (oclambda);
}


double OrderConditionsSymLambdaHess(JdTau2 H/*::AbstractMatrix*/, JdTau GL/*::AbstractVector*/,int maxorder/*::Int64*/, OrderConditionsCache *cache/*::OrderConditionsCache*/, JdTau coeffs/*::AbstractVector*/, JdTau lambda/*::AbstractVector*/)
{
int nereind,i,j,k,l,m,s,ncoeffs,noc,N,n,ntrees; 
int j23,j3, jp23,jp3, jh23,jh3,jc23,jc3, jhc23,jhc3, jhp23,jhp3;
double aux,haux,oclambda;
JintTau *symind,*treeindices;     //::Vector{Int64}
TreeSet *treesetptr;             //::TreeSet
JdTau2 *coeffpow, *phi, *phip, *hcoeffpow, *hphi, *hphip;   //::Array{ctype,2}
JintTau *firstti, *tp, *dec1, *dec2;    //::Vector{Int64}
JdTau3 *Jcoeffpow, *Jphi, *Jphip, *Jhcoeffpow, *Jhphi, *Jhphip;//::Array{ctype,3}
JdTau *Jaux, *Jhaux;//::Array{ctype,1}

if (iseven(maxorder)) maxorder --; // hau nik jarritakoa da.
    ncoeffs = coeffs.kop;  //  ncoeffs = length(coeffs)
    s = 2*ncoeffs-1;
    
    //@unpack symind, treeindices, treeset, 
    //        coeffpow, phi, phip, hcoeffpow, hphi, hphip,
    //        Jcoeffpow, Jphi, Jphip, Jhcoeffpow, Jhphi, Jhphip, Jaux, Jhaux = cache

    symind = &(cache->symind);
    treeindices = &(cache->treeindices);
    treesetptr = cache->treesetptr;
    coeffpow = &(cache->coeffpow);
    phi = &(cache->phi);
    phip = &(cache->phip);
    hcoeffpow = &(cache->hcoeffpow);
    hphi = &(cache->hphi);
    hphip = &(cache->hphip);
    Jcoeffpow  =  &(cache->Jcoeffpow  );
    Jphi  =  &(cache->Jphi  );
    Jphip  =  &(cache->Jphip  );
    Jhcoeffpow  =  &(cache->Jhcoeffpow  );
    Jhphi  =  &(cache->Jhphi  );
    Jhphip  =  &(cache->Jhphip  );
    Jaux  =  &(cache->Jaux  );
    Jhaux  =  &(cache->Jhaux  );

    //    @unpack firstti, tp, dec1, dec2 = treeset
    firstti = &(treesetptr->firstti);
    tp = &(treesetptr->tp);
    dec1 = &(treesetptr->dec1);
    dec2 = &(treesetptr->dec2);

 for (i=0;i<(Jcoeffpow->kop*Jcoeffpow->kop2*Jcoeffpow->kop3);i++) Jcoeffpow->dptr[i] =0;  // Jcoeffpow = zeros(ctype,m,maxorder,m+s)
 for (i=0;i<(Jphi->kop*Jphi->kop2*Jphi->kop3);i++) Jphi->dptr[i] =0; //   Jphi = zeros(ctype,s+1,ntrees,m+s)
 for (i=0;i<(Jphip->kop*Jphip->kop2*Jphip->kop3);i++) Jphip->dptr[i] =0;//    Jphip = zeros(ctype,s,ntrees,m+s)
 for (i=0;i<(Jhcoeffpow->kop*Jhcoeffpow->kop2*Jhcoeffpow->kop3);i++) Jhcoeffpow->dptr[i] =0; //    Jhcoeffpow = zeros(ctype,m,maxorder,m+s)
 for (i=0;i<(Jhphi->kop*Jhphi->kop2*Jhphi->kop3);i++) Jhphi->dptr[i] =0;//    Jhphi = zeros(ctype,s+1,ntrees,m+s)
 for (i=0;i<(Jhphip->kop*Jhphip->kop2*Jhphip->kop3);i++) Jhphip->dptr[i] =0;//     Jhphip = zeros(ctype,s,ntrees,m+s)
 for (i=0;i<Jaux->kop;i++) Jaux->dptr[i] =0;//     Jaux = zeros(ctype,m+s)
 for (i=0;i<Jhaux->kop;i++) Jhaux->dptr[i] =0;//    Jhaux = zeros(ctype,m+s)



    
    //printf("unpakeatuta!\n");
    
    jc3= Jcoeffpow->kop3;
    jc23 = Jcoeffpow->kop2 * jc3; 
    
    j3= Jphi->kop3;
    j23 = Jphi->kop2 * j3; 
      
    jp3= Jphip->kop3;
    jp23 = Jphip->kop2 * jp3;
    
    jhc3= Jhcoeffpow->kop3;
    jhc23 = Jhcoeffpow->kop2 * jhc3;
    
    jh3= Jhphi->kop3;
    jh23 = Jhphi->kop2 * jh3;
    
    jhp3= Jhphip->kop3;
    jhp23 = Jhphip->kop2 * jhp3;

    noc = treeindices->kop;  //   noc = length(treeindices)
    N = ncoeffs+noc;
    if (lambda.kop < noc)  //length(lambda)<noc ? error("length(lambda) < ", noc) : nothing
        {
        printf("error: length(lambda) < %d\n",noc);
        exit(0);
        }
    if (GL.kop < (noc+ncoeffs)) //length(GL)<(noc+ncoeffs) ? error("length(GL) < ", noc+coeffs) : nothing
        {
        printf("error: length(GL) < %d\n",noc+ncoeffs);
        exit(0);
        }
   
    for (i=0; i<ncoeffs; i++) // i in 1:ncoeffs
        {
        coeffpow->dptr[i*coeffpow->kop2+0] = coeffs.dptr[i];        //coeffpow[i,1] = coeffs[i]
        coeffpow->dptr[i*coeffpow->kop2+1] = coeffs.dptr[i]*coeffs.dptr[i]; //coeffpow[i,2] = coeffpow[i,1]^2
        Jcoeffpow->dptr[jc23*i + i]= 1;       //Jcoeffpow[i,1,i] = 1
        Jcoeffpow->dptr[jc23*i + jc3 + i]= 2.0*coeffpow->dptr[i*coeffpow->kop2] * 
                              Jcoeffpow->dptr[jc23*i + i]; // Jcoeffpow[i,2,i] = 2*coeffpow[i,1]*Jcoeffpow[i,1,i]
        }
    
    for (k=2; k<maxorder; k+=2)   //  k in 3:2:maxorder
        {
        for  (i=0; i<ncoeffs; i++) // i in 1:ncoeffs
            {
            coeffpow->dptr[i*coeffpow->kop2+k] = 
                   coeffpow->dptr[i*coeffpow->kop2+(k-2)]*
                   coeffpow->dptr[i*coeffpow->kop2+1]; // coeffpow[i,k] = coeffpow[i,k-2]*coeffpow[i,2]
            //  Jcoeffpow[i,k,i] = Jcoeffpow[i,k-2,i]*coeffpow[i,2] + coeffpow[i,k-2]*Jcoeffpow[i,2,i]
            Jcoeffpow->dptr[jc23*i + jc3*k + i]=  
                   Jcoeffpow->dptr[jc23*i + jc3*(k-2) + i] * coeffpow->dptr[coeffpow->kop2*i+1] +  
                   coeffpow->dptr[coeffpow->kop2*i+k-2] * Jcoeffpow->dptr[jc23*i + jc3 + i];  
            }
        }
     //printf("koeffpowak kalkulatua!\n");

    ntrees = firstti->iptr[maxorder]-1;   //   firstti[maxorder+1]-1
    for (n=0;n<ntrees; n++)   //  n in 1:ntrees
        {
        m = dec1->iptr[n]-1; // dec1[n]
        l = dec2->iptr[n]-1; // dec2[n]
        k = tp->iptr[n]-1;   // tp[n]

        if (l==-1) //l==0
            {  
            for (j=0; j<s; j++)   //   j in 1:s
                {
                i = symind->iptr[j]-1;  //  symind[j]
                phip->dptr[phip->kop2*j+n] = coeffpow->dptr[coeffpow->kop2*i+k]/2.0; // phip[j,n] = coeffpow[i,k]/2 
                /*if ((j==0)&&(i==0))
                    printf("Jphip[%d,%d,%d] = Jcoeffpow[i,k,i]/2: %.12lg = %.12lg\n",j,n,i, Jphip->dptr[jp23*j+jp3*n+i], Jcoeffpow->dptr[jc23*i+jc3*k+i]/2.0);*/
                
                Jphip->dptr[jp23*j+jp3*n+i] = Jcoeffpow->dptr[jc23*i+jc3*k+i]/2.0; // Jphip[j,n,i] = Jcoeffpow[i,k,i]/2
                }
            }
        else
            { 
            for (j=0; j<s;j++)  // j in 1:s
                {
                                 // phip[j,n] = phip[j,m]*phi[j,l]
                phip->dptr[phip->kop2*j+n] = phip->dptr[phip->kop2*j+m] * phi->dptr[phi->kop2*j+l];
                for (i=0;i<ncoeffs; i++)  // i in 1:ncoeffs
                    {
                    /*if ((j==0)&&(i==0))
                        printf("Jphip[%d,%d,%d] = Jphip[j,m,i]*phi[j,l]+phip[j,m]*Jphi[j,l,i]: %.12lg = %.12lg * %.12lg + %.12lg * %.12lg \n",j,n,i, Jphip->dptr[jp23*j+jp3*n+i], Jphip->dptr[jp23*j+jp3*m+i],  phi->dptr[phi->kop2*j+l],
                            phip->dptr[phip->kop2*j+m], Jphi->dptr[j23*j+j3*l+i]);*/
                    // Jphip[j,n,i] = Jphip[j,m,i]*phi[j,l] + phip[j,m]*Jphi[j,l,i]
                    Jphip->dptr[jp23*j+jp3*n+i] = Jphip->dptr[jp23*j+jp3*m+i] * phi->dptr[phi->kop2*j+l] +
                            phip->dptr[phip->kop2*j+m] * Jphi->dptr[j23*j+j3*l+i];
                    }
                }
            }
        phi->dptr[n] = phip->dptr[n];  // phi[1,n] = phip[1,n]

        for (i=0;i<ncoeffs; i++)  //  i in 1:ncoeffs
             Jphi->dptr[j3*n+i] = Jphip->dptr[jp3*n+i];  //  Jphi[1,n,i] = Jphip[1,n,i]
        
        for (j=1; j<s; j++)// j in 2:s
            {
            // phi[j,n] = phi[j-1,n] + phip[j-1,n] + phip[j,n]
            phi->dptr[phi->kop2*j+n] = phi->dptr[phi->kop2*(j-1)+n] + phip->dptr[phip->kop2*(j-1)+n] + phip->dptr[phip->kop2*j+n];
            for (i=0; i<ncoeffs; i++)  // i in 1:ncoeffs
                {
                //  Jphi[j,n,i] = Jphi[j-1,n,i] + Jphip[j-1,n,i] + Jphip[j,n,i]
                Jphi->dptr[j23*j+j3*n+i] = Jphi->dptr[j23*(j-1)+j3*n+i] + Jphip->dptr[jp23*(j-1)+jp3*n+i] + 
                                            Jphip->dptr[jp23*j+jp3*n+i];
                }
            }
        }
    // printf("hasierako phi, phip eta jphi  kalkulatua!\n");

    for (j=1; j<noc; j++)  //  j in 2:noc
        {
        n = treeindices->iptr[j]-1; // treeindices[j]
        // phi[s+1,n] = phi[s,n] + phip[s,n] 
        phi->dptr[phi->kop2*s+n] = phi->dptr[phi->kop2*(s-1)+n] + phip->dptr[phip->kop2*(s-1)+n];
        for (i=0; i<ncoeffs; i++)  // i in 1:ncoeffs
            {
            // Jphi[s+1,n,i] = Jphi[s,n,i] + Jphip[s,n,i]    
            Jphi->dptr[j23*s+j3*n+i] = Jphi->dptr[j23*(s-1)+j3*n+i] + Jphip->dptr[jp23*(s-1)+jp3*n+i];   
            } 
        }

  //phi[s+1,1] = sum(coeffpow[i,2] for i=1:ncoeffs-1) + coeffpow[ncoeffs,2]/2 - 9
    phi->dptr[phi->kop2*s] = -9.0;
    for (i = 0; i< ncoeffs-1; i++)
        phi->dptr[phi->kop2*s] +=coeffpow->dptr[coeffpow->kop2*i+1];
    phi->dptr[phi->kop2*s] += (coeffpow->dptr[coeffpow->kop2*(ncoeffs-1)+1]/2.0);
    
    for (i=0; i<(ncoeffs-1); i++)  //  i in 1:ncoeffs-1
        {
        // Jphi[s+1,1,i] = Jcoeffpow[i,2,i]
        Jphi->dptr[j23*s+i] = Jcoeffpow->dptr[jc23*i+jc3+i];
        }
    // Jphi[s+1,1,ncoeffs] = Jcoeffpow[ncoeffs,2,ncoeffs]/2
    Jphi->dptr[j23*s+ncoeffs-1] = Jcoeffpow->dptr[jc23*(ncoeffs-1)+jc3+ncoeffs-1]/2.0;

 // aux = sum(coeffpow[i,1] for i=1:ncoeffs-1) + coeffpow[ncoeffs,1]/2  
    aux = 0.0;
    for (i = 0; i< ncoeffs-1; i++)
        aux +=coeffpow->dptr[coeffpow->kop2*i];
    aux += (coeffpow->dptr[coeffpow->kop2*(ncoeffs-1)]/2.0); 

    for (j=0; j<ncoeffs; j++)//  j in 1:ncoeffs  
        {                                                          
        // Jaux[j] = sum(Jcoeffpow[i,1,j] for i=1:ncoeffs-1) + Jcoeffpow[ncoeffs,1,j]/2 
        Jaux->dptr[j] = Jcoeffpow->dptr[jc23*(ncoeffs-1)+j]/2.0;
        for (i = 0; i< ncoeffs-1; i++)
            {
            Jaux->dptr[j] += Jcoeffpow->dptr[jc23*i+j];
            }
        }                                                                    
    oclambda = aux*aux; // aux^2
    for (i=0; i<treeindices->kop; i++) // i in eachindex(treeindices)
        {
        n = treeindices->iptr[i] -1; // treeindices[i]    
        oclambda +=  phi->dptr[phi->kop2*s+n]*lambda.dptr[i];  // phi[s+1,n]*lambda[i]  
        }
    /*printf("oclambda kalkulatuta! = %lg\n",oclambda);
for (i=0;i<24;i++) printf("Jphip[0,%d,0]= %.12lg\n",i, Jphip->dptr[jp3*i]);*/
    // # Computation of the gradient with reverse mode algorithmic differentiation
    //hcoeffpow .= 0
    for (i=0; i < hcoeffpow->kop*hcoeffpow->kop2; i++) hcoeffpow->dptr[i] = 0.0; 
//printf("jphip[0,%d,0]= %.12lg\n",0, Jphip->dptr[0]);

    // hphi .= 0 
    for (i=0; i < hphi->kop*hphi->kop2; i++) hphi->dptr[i] = 0.0; 
//printf("jphip[0,%d,0]= %.12lg\n",0, Jphip->dptr[0]);

    // hphip .= 0
    for (i=0; i < hphip->kop*hphip->kop2; i++) hphip->dptr[i] = 0.0; 
//printf("jphip[0,%d,0]= %.12lg\n",0, Jphip->dptr[0]);

    // GL .= 0
    for (i=0; i < GL.kop; i++) GL.dptr[i] = 0.0;
//printf("jphip[0,%d,0]= %.12lg\n",0, Jphip->dptr[0]);

    // Jhcoeffpow .= 0
    for (i=0; i < Jhcoeffpow->kop*Jhcoeffpow->kop2*Jhcoeffpow->kop3; i++) Jhcoeffpow->dptr[i] = 0.0;
//printf("jphip[0,%d,0]= %.12lg\n",0, Jphip->dptr[0]);

    // Jhphi .= 0
    for (i=0; i < Jhphi->kop*Jhphi->kop2*Jhphi->kop3; i++) Jhphi->dptr[i] = 0.0;
//printf("jphip[0,%d,0]= %.12lg\n",0, Jphip->dptr[0]);

    // Jhphip .= 0
    for (i=0; i < Jhphip->kop*Jhphip->kop2*Jhphip->kop3; i++) Jhphip->dptr[i] = 0.0;
//printf("jphip[0,%d,0]= %.12lg\n",0, Jphip->dptr[0]);

    // H .= 0
    for (i=0; i < H.kop*H.kop2; i++) H.dptr[i] = 0.0; 
    haux = 0;
//for (i=0;i<24;i++) printf("jphip[0,%d,0]= %.12lg\n",i, Jphip->dptr[jp3*i]);
    for (j=0; j<noc; j++)  // j in 1:noc
        {
        n = treeindices->iptr[j]-1; // treeindices[j] 
        // GL[ncoeffs+j] += phi[s+1,n]
        GL.dptr[ncoeffs+j] += phi->dptr[phi->kop2*s+n];
        // hphi[s+1,n] += lambda[j]  
        hphi->dptr[hphi->kop2*s+n] += lambda.dptr[j];  
        // Jhphi[s+1,n,ncoeffs+j] += 1  
        Jhphi->dptr[jh23*s+jh3*n+ncoeffs+j] += 1; 
        for (i=0; i<ncoeffs; i++)  // i in 1:ncoeffs
            {
            // H[ncoeffs+j,i] += Jphi[s+1,n,i]
            H.dptr[H.kop2*(ncoeffs+j)+i] += Jphi->dptr[j23*s+j3*n+i];
            }
        } 
    //printf("Hasieraketak eta hasierako H kalkulatua!\n");
 
    // # Gradient of aux^2                                           
    haux += 2*aux;
    for (j=0; j<ncoeffs; j++) // j in 1:ncoeffs
        {
        // Jhaux[j] += 2*Jaux[j]
        Jhaux->dptr[j] += 2.0*Jaux->dptr[j];
        }
    for (i=0; i<ncoeffs-1;i++) // i in 1:ncoeffs-1 
        {
        hcoeffpow->dptr[hcoeffpow->kop2*i] += haux; // hcoeffpow[i,1] += haux
        for (j=0; j<N; j++) // j in 1:N
            {
            // Jhcoeffpow[i,1,j] += Jhaux[j]
            Jhcoeffpow->dptr[jhc23*i+j] += Jhaux->dptr[j];
            }
        }
    hcoeffpow->dptr[hcoeffpow->kop2*(ncoeffs-1)] += haux/2.0;  // hcoeffpow[ncoeffs,1] += haux/2  
    for (j=0; j<N; j++) // j in 1:N
            {
            // Jhcoeffpow[ncoeffs,1,j] += Jhaux[j]/2
            Jhcoeffpow->dptr[jhc23*(ncoeffs-1)+j] += Jhaux->dptr[j]/2.0;
            }                                                    
    // # end of gradient of aux^2
/* 
printf("aux^2-ren gradientea kalkulatua!\n");
printf("aux ondoren %.12lg, eta %.12lg\n",Jhcoeffpow->dptr[0],Jhcoeffpow->dptr[1]);
*/
    for (i=0; i<ncoeffs-1; i++)  // i in 1:ncoeffs-1
        {
        hcoeffpow->dptr[hcoeffpow->kop2*i+1] += hphi->dptr[hphi->kop2*s]; // hcoeffpow[i,2] += hphi[s+1,1]
        for (j=0; j<N; j++) //  j in 1:N
            {
            // Jhcoeffpow[i,2,j] += Jhphi[s+1,1,j]
            Jhcoeffpow->dptr[jhc23*i+jhc3+j] += Jhphi->dptr[jh23*s+j];
            }
        }

    // hcoeffpow[ncoeffs,2] += hphi[s+1,1]/2
    hcoeffpow->dptr[hcoeffpow->kop2*(ncoeffs-1)+1] += hphi->dptr[hphi->kop2*s]/2.0;
    for (j=0; j<N; j++) // j in 1:N
        {
        //  Jhcoeffpow[ncoeffs,2,j] += Jhphi[s+1,1,j]/2
        Jhcoeffpow->dptr[jhc23*(ncoeffs-1)+jhc3+j] += Jhphi->dptr[jh23*s+j]/2.0;
        }

    for (j=1; j<noc;j++)  //  j in 2:noc
        {
        n = treeindices->iptr[j]-1; //  treeindices[j]
        hphi->dptr[hphi->kop2*(s-1)+n] += hphi->dptr[hphi->kop2*s+n];     //  hphi[s,n] += hphi[s+1,n]
        hphip->dptr[hphip->kop2*(s-1)+n] += hphi->dptr[hphi->kop2*s+n];   //  hphip[s,n] += hphi[s+1,n]
        for (i=0; i<N; i++) // i in 1:N
            {
            //   Jhphi[s,n,i] += Jhphi[s+1,n,i]
            Jhphi->dptr[jh23*(s-1)+jh3*n+i] += Jhphi->dptr[jh23*s+jh3*n+i];
            //   Jhphip[s,n,i] += Jhphi[s+1,n,i]
            Jhphip->dptr[jhp23*(s-1)+jhp3*n+i] += Jhphi->dptr[jh23*s+jh3*n+i];
            }
        }
//printf("mlk-ra sartzera\n");

    for (n=firstti->iptr[maxorder]-2; n>=0; n--)  // n in firstti[maxorder+1]-1:-1:1
        {
         //   printf("n = %d\n",n );
        m = dec1->iptr[n]-1;   // dec1[n]
        l = dec2->iptr[n]-1;  // dec2[n]
        k = tp->iptr[n]-1;    // tp[n]
        for (j=s-1; j>0; j--)  //  j in s:-1:2
            {
            hphi->dptr[hphi->kop2*(j-1)+n] += hphi->dptr[hphi->kop2*j+n];     // hphi[j-1,n] += hphi[j,n]
            hphip->dptr[hphip->kop2*(j-1)+n] += hphi->dptr[hphi->kop2*j+n];   // hphip[j-1,n] += hphi[j,n]
            hphip->dptr[hphip->kop2*j+n] += hphi->dptr[hphi->kop2*j+n];   // hphip[j,n] += hphi[j,n]
            for (i=0; i<N; i++) // i in 1:N
                {
/* if ((j==1)&&(n==0)&&(i==0)) 
      printf("j== 1 eta Jhphi[j-1,n,i] += Jhphi[j,n,i]: %.12lg += %.12lg\n",
            Jhphi->dptr[jh23*(j-1)+jh3*n+i], Jhphi->dptr[jh23*j+jh3*n+i]);
*/
                // Jhphi[j-1,n,i] += Jhphi[j,n,i]
                Jhphi->dptr[jh23*(j-1)+jh3*n+i] += Jhphi->dptr[jh23*j+jh3*n+i];
                // Jhphip[j-1,n,i] += Jhphi[j,n,i]
                Jhphip->dptr[jhp23*(j-1)+jhp3*n+i] += Jhphi->dptr[jh23*j+jh3*n+i];
                // Jhphip[j,n,i] += Jhphi[j,n,i]
                Jhphip->dptr[jhp23*j+jhp3*n+i] += Jhphi->dptr[jh23*j+jh3*n+i];
                }
            }
            
        hphip->dptr[n] += hphi->dptr[n]; //  hphip[1,n] += hphi[1,n]
        for (i=0; i<N; i++) // i in 1:N
            {
/*
if ((n==0)&&(i==0)) printf("mlk barruan Jhphip[1,1,1]:%.12lg + (jhphi[1,1,1]) %.12lg \n",Jhphip->dptr[0],Jhphi->dptr[0]);
*/

            //  Jhphip[1,n,i] += Jhphi[1,n,i]
            Jhphip->dptr[jhp3*n+i] += Jhphi->dptr[jh3*n+i];
            }
            
        if (l==-1)   //  l==0
            {
            for (j=0; j<s; j++) // j in 1:s
                {
                //  hcoeffpow[symind[j],k] += hphip[j,n]/2
                nereind= symind->iptr[j]-1;
                hcoeffpow->dptr[hcoeffpow->kop2*nereind+k] +=  hphip->dptr[hphip->kop2*j+n]/2.0;
                for (i=0; i<N; i++) // i in 1:N
                    {
/*
if ((nereind==0)&&(k==0)&&(i==0)) printf("mlk barruan %.12lg + (jhphip[%d,%d,%d]) %.12lg \n",Jhcoeffpow->dptr[0],j,n,i,Jhphip->dptr[jhp23*j+jhp3*n+i]/2.0);
*/

                    //   Jhcoeffpow[symind[j],k,i] += Jhphip[j,n,i]/2
                    Jhcoeffpow->dptr[jhc23*nereind+jhc3*k+i] += Jhphip->dptr[jhp23*j+jhp3*n+i]/2.0;
                    }
                }
            }
        else
            {
            for (j=s-1; j>=0; j--)  // j in s:-1:1
                {
                // hphip[j,m] += hphip[j,n]*phi[j,l]
                
                hphip->dptr[hphip->kop2*j+m] += hphip->dptr[hphip->kop2*j+n] *
                                              phi->dptr[phi->kop2*j+l]; 
                // hphi[j,l] += hphip[j,n]*phip[j,m] 
                hphi->dptr[hphi->kop2*j+l] += hphip->dptr[hphip->kop2*j+n] *
                                              phip->dptr[phip->kop2*j+m]; 
                for (i=0; i<N; i++) // i in 1:N
                    {
/*
if ((j==0)&&(l==0)&&(i==0))
                        printf("Jhphi[1,1,1] += Jhphip[j,n,i]*phip[j,m] + hphip[j,n]*Jphip[%d,%d,%d]: %.12lg + %.12lg * %.12lg + %.12lg * %12lg\n ",j,m,i, Jhphi->dptr[jh23*j+jh3*l+i], Jhphip->dptr[jhp23*j+jhp3*n+i], phip->dptr[phip->kop2*j+m],  hphip->dptr[hphip->kop2*j+n], Jphip->dptr[jp23*j+jp3*m+i]);
*/
                    // Jhphip[j,m,i] += Jhphip[j,n,i] * phi[j,l]  +  hphip[j,n] * Jphi[j,l,i]
                     Jhphip->dptr[jhp23*j+jhp3*m+i] += 
                              Jhphip->dptr[jhp23*j+jhp3*n+i] *  phi->dptr[phi->kop2*j+l] +
                              hphip->dptr[hphip->kop2*j+n] * Jphi->dptr[j23*j+j3*l+i];
                    // Jhphi[j,l,i] += Jhphip[j,n,i] * phip[j,m]  +  hphip[j,n] * Jphip[j,m,i]
                    Jhphi->dptr[jh23*j+jh3*l+i] += 
                              Jhphip->dptr[jhp23*j+jhp3*n+i] * phip->dptr[phip->kop2*j+m] + 
                              hphip->dptr[hphip->kop2*j+n] * Jphip->dptr[jp23*j+jp3*m+i]; 
                    }
                }
            }
        }
/*
printf("mlk ondoren\n");
printf("mlk ondoren %.12lg, eta %.12lg\n",Jhcoeffpow->dptr[0],Jhcoeffpow->dptr[1]);
*/
    for (k=maxorder-1; k>1; k-=2)   //  k in maxorder:-2:3
        {
        for (j=0; j<ncoeffs; j++)  //  j in 1:ncoeffs
            {
            // hcoeffpow[j,k-2] += hcoeffpow[j,k]*coeffpow[j,2]
            hcoeffpow->dptr[hcoeffpow->kop2*j+k-2] += hcoeffpow->dptr[hcoeffpow->kop2*j+k] *
                                     coeffpow->dptr[coeffpow->kop2*j+1];
            // hcoeffpow[j,2] += hcoeffpow[j,k]*coeffpow[j,k-2]
            hcoeffpow->dptr[hcoeffpow->kop2*j+1] += hcoeffpow->dptr[hcoeffpow->kop2*j+k] *
                                     coeffpow->dptr[coeffpow->kop2*j+k-2];
            for (i=0; i<N; i++) // i in 1:N
                {
                // Jhcoeffpow[j,k-2,i] += Jhcoeffpow[j,k,i]*coeffpow[j,2] + hcoeffpow[j,k]*Jcoeffpow[j,2,i]
                Jhcoeffpow->dptr[jhc23*j+jhc3*(k-2)+i] +=
                                Jhcoeffpow->dptr[jhc23*j+jhc3*k+i] * coeffpow->dptr[coeffpow->kop2*j+1] +
                                hcoeffpow->dptr[hcoeffpow->kop2*j+k] * Jcoeffpow->dptr[jc23*j+jc3+i];
                // Jhcoeffpow[j,2,i] += Jhcoeffpow[j,k,i]*coeffpow[j,k-2] + hcoeffpow[j,k]*Jcoeffpow[j,k-2,i]
                Jhcoeffpow->dptr[jhc23*j+jhc3+i] += 
                                Jhcoeffpow->dptr[jhc23*j+jhc3*k+i] *coeffpow->dptr[coeffpow->kop2*j+k-2] +
                                hcoeffpow->dptr[hcoeffpow->kop2*j+k] * Jcoeffpow->dptr[jc23*j+jc3*(k-2)+i];
                }
            }
        } 
/*
printf("azkenaren ondoren %.12lg, eta %.12lg\n",Jhcoeffpow->dptr[0],Jhcoeffpow->dptr[1]);
*/
    for (j=0; j<ncoeffs; j++)  // j in 1:ncoeffs
        {
        // hcoeffpow[j,1] += 2*coeffpow[j,1]*hcoeffpow[j,2]
        hcoeffpow->dptr[hcoeffpow->kop2*j] += 2.0 * coeffpow->dptr[coeffpow->kop2*j] * hcoeffpow->dptr[hcoeffpow->kop2*j+1];
        // GL[j] = hcoeffpow[j,1]
        GL.dptr[j] = hcoeffpow->dptr[hcoeffpow->kop2*j];
        for (i=0; i<N; i++)  //  i in 1:N 
            {
/*
if (j==0)
{
    if ((i==0)||(i==1))
        printf("%d eta %d balioentzat %.12lg + %.12lg  * %.12lg + %.12lg * %.12lg \n",j,i, Jhcoeffpow->dptr[jhc23*j+i], Jcoeffpow->dptr[jc23*j+i], hcoeffpow->dptr[hcoeffpow->kop2*j+1], coeffpow->dptr[coeffpow->kop2*j], Jhcoeffpow->dptr[jhc23*j+jhc3+i]);
}
*/
            // Jhcoeffpow[j,1,i] += 2*Jcoeffpow[j,1,i]*hcoeffpow[j,2] + 2*coeffpow[j,1]*Jhcoeffpow[j,2,i]
            Jhcoeffpow->dptr[jhc23*j+i] += 
                      2.0 * Jcoeffpow->dptr[jc23*j+i] * hcoeffpow->dptr[hcoeffpow->kop2*j+1] +
                      2.0 * coeffpow->dptr[coeffpow->kop2*j] * Jhcoeffpow->dptr[jhc23*j+jhc3+i];
            // H[j,i] = Jhcoeffpow[j,1,i]
            H.dptr[H.kop2*j+i] = Jhcoeffpow->dptr[jhc23*j+i];
/*
if ((j==0) && ((i==0)||(i==1)))
    printf("eragiketa egin ondoren %lg\n",H.dptr[i]);
*/
            }
        }
    return (oclambda);
}




int main(int argc, char *argv[])
{
int i,j,maxorder,nx,noc;
OrderConditionsCache *cache;
JdTau sofro35,lambda0,GL;
JdTau2 H;
double emaitza;

TreeSet *zptr;

//zptr= ConstructInfTrees(9,isodd);

maxorder = 10;
nx = 18;

sofro35.kop=nx;
sofro35.dptr=(double *)malloc(nx*sizeof(double));
sofro35.dptr[0] = 0.078795722521686419263907679337684;
sofro35.dptr[1] = 0.31309610341510852776481247192647;
sofro35.dptr[2] = 0.027918383235078066109520273275299;
sofro35.dptr[3] = -0.22959284159390709415121339679655;
sofro35.dptr[4] = 0.13096206107716486317465685927961;
sofro35.dptr[5] = -0.26973340565451071434460973222411;
sofro35.dptr[6] = 0.074973343155891435666137105641410;
sofro35.dptr[7] = 0.11199342399981020488957508073640;
sofro35.dptr[8] = 0.36613344954622675119314812353150;
sofro35.dptr[9] = -0.39910563013603589787862981058340;
sofro35.dptr[10] = 0.10308739852747107731580277001372;
sofro35.dptr[11] = 0.41143087395589023782070411897608;
sofro35.dptr[12] = -0.0048663605831352617621956593099771;
sofro35.dptr[13] = -0.39203335370863990644808193642610;
sofro35.dptr[14] = 0.051942502962449647037182904015976;
sofro35.dptr[15] = 0.050665090759924496335874344156866;
sofro35.dptr[16] = 0.049674370639729879054568800279461;
sofro35.dptr[17] = 0.049317735759594537917680008339338;

cache = OrderConditionsSymInit(maxorder,&sofro35);
printf("hasieraketak eginda... \n");
noc = cache->treeindices.kop; //length(cache.treeindices)
for (i=0;i < cache->treeindices.kop; i++) 
    printf("treeindices[%d]=%d\n",i,cache->treeindices.iptr[i]);
printf("baldintza kopurua = %d\n",noc );
lambda0.kop = noc;
lambda0.dptr=(double *)malloc((noc)*sizeof(double));
for (i=0;i<noc;i++) lambda0.dptr[i]=1.0; // =ones(ctype,noc)

GL.kop = nx+noc; //Array{ctype}(undef,nx+noc)
GL.dptr=(double *)malloc((nx+noc)*sizeof(double));
H.kop = (nx+noc);
H.kop2 = (nx+noc); //Array{ctype}(undef,nx+noc)
H.dptr=(double *)malloc((nx+noc)*(nx+noc)*sizeof(double));

//emaitza = OrderConditionsSymLambdaGrad(GL, maxorder,cache, sofro35, lambda0);
emaitza = OrderConditionsSymLambdaHess(H,GL, maxorder,cache, sofro35, lambda0);


zptr = cache->treesetptr;
printf("emaitza = %.12g, gradu maximoa = %d\n",emaitza, zptr->maxdegree);
/*
for (i=0;i < cache->coeffpow.kop; i++) 
    {
    printf("coeffpow[%d,j] ",i+1);
    for (j=0;j<cache->coeffpow.kop2; j++) printf("  %.8lg, ",cache->coeffpow.dptr[i*cache->coeffpow.kop2+j]);
    printf("\n");
    }
    
for (i=0;i < cache->phip.kop; i++) 
    {
    printf("phip[%d,j] ",i+1);
    for (j=0;j<cache->phip.kop2; j++) printf("  %.8lg, ",cache->phip.dptr[i*cache->phip.kop2+j]);
    printf("\n");
    }
    
for (i=0;i < cache->phi.kop; i++) 
    {
    printf("phi[%d,j] ",i+1);
    for (j=0;j<cache->phi.kop2; j++) printf("  %.8lg, ",cache->phi.dptr[i*cache->phi.kop2+j]);
    printf("\n");
    }
 */
/*
for (i=0;i < cache->coeffpow.kop; i++) printf("coeffpow[%d,3]=%.8lf\n",i+1,cache->coeffpow.dptr[i*cache->coeffpow.kop2+2]);
for (i=0;i < cache->coeffpow.kop; i++) printf("coeffpow[%d,5]=%.8lf\n",i+1,cache->coeffpow.dptr[i*cache->coeffpow.kop2+4]);
printf("hcoeffpow\n");
for (i=0;i < cache->hcoeffpow.kop; i++) printf("hcoeffpow[%d,2]=%.8lf, [%d,3]=%.8lf \n",i+1,cache->hcoeffpow.dptr[i*cache->hcoeffpow.kop2+1], i+1,cache->hcoeffpow.dptr[i*cache->hcoeffpow.kop2+2] );
for (i=0;i < cache->hcoeffpow.kop; i++) printf("hcoeffpow[%d,4]=%.8lf, [%d,5]=%.8lf \n",i+1,cache->hcoeffpow.dptr[i*cache->hcoeffpow.kop2+3], i+1,cache->hcoeffpow.dptr[i*cache->hcoeffpow.kop2+4] );
for (i=0;i < cache->hcoeffpow.kop; i++) printf("hcoeffpow[%d,6]=%.8lf, [%d,7]=%.8lf \n",i+1,cache->hcoeffpow.dptr[i*cache->hcoeffpow.kop2+5], i+1,cache->hcoeffpow.dptr[i*cache->hcoeffpow.kop2+6] );
for (i=0;i < cache->hcoeffpow.kop; i++) printf("hcoeffpow[%d,8]=%.8lf, [%d,9]=%.8lf \n",i+1,cache->hcoeffpow.dptr[i*cache->hcoeffpow.kop2+7], i+1,cache->hcoeffpow.dptr[i*cache->hcoeffpow.kop2+8] );
*/
//for (i=0;i < cache->hcoeffpow.kop; i++) printf("hcoeffpow[%d,8]=%.8lf\n",i+1,cache->hcoeffpow.dptr[i*cache->hcoeffpow.kop2+7]);
//for (i=0;i < zptr->tp.kop; i++) printf("tp[%d]=%d\n",i,zptr->tp.iptr[i]);
//for (i=0;i < zptr->degree.kop; i++) printf("degree[%d]=%d\n",i,zptr->degree.iptr[i]);
//for (i=0;i < zptr->dec1.kop; i++) printf("dec1[%d]=%d eta dec2[%d] = %d\n", i+1, zptr->dec1.iptr[i], i+1, zptr->dec2.iptr[i]);
printf("GLren elementu kopurua = %d\n",nx+noc);
for (i=0;i < GL.kop; i++) printf("GL[%d]=%.6lg\n",i,GL.dptr[i]);
printf("Hko kop = %d, eta kop2 = %d\n",H.kop,H.kop2 );
/* * /
for (i=0;i < H.kop; i++) 
    {
    printf("Jakobiarreko %d lerroa: H[%d,:]= \n",i, i);
    for (j=0;j<H.kop2; j++) printf("%.6lg,",H.dptr[H.kop2*i+j]);
    printf("\n");
    }
  /*  */
}






