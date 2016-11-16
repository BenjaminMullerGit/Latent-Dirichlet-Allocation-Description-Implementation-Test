library("coda")
library("MASS")
library("MCMCpack")
#install.packages("Matrix")
library("Matrix")
library("rootSolve")




#Nd,d -> Ndd
algo<-function(k,Ndd,d,epsilon,matrice,alpha,beta)
{ 
  pas=10
  #initialisations des param?tres variationnels
  #phipast=(1/k)*matrix(1:1,nrow=k,ncol=Nd[d])
  phi=(1/k)*matrix(1:1,nrow=k,ncol=Ndd)
  phipast=phi
  gammat=alpha+Ndd/k
  gammatpast=gammat
  #gammatpast=c()
  #gammatpast=alpha+Nd[d]/k
  #gammat=c()
  
  #S pour la normalisations 
  #S=c()
  #S[1]=1
  
  P=10
  G=10
  
  while(epsilon<pas)
  {
    for (n in 1:Ndd)
    {
      #calcul de phi ? n fix?
      for (i in 1:k)
      {
        indice_mot=matrice[d,n]
        phi[i,n]=beta[i,indice_mot]*exp(digamma(gammat[i]))
      }
      #normalisation du phi ? n fix?  
      #S[n]=colSums(phi)[n]  
      #phi[,n]=phi[,n]/S[n]
      #print(phi[,n])
      phi[,n]=phi[,n]/sum(phi[,n])
     # print(phi[,n])    
    }   
    #calcul de gammat
    #for (i in 1:k)
    
    #gammat[i]=alpha[i]+sum(phi[i,n],n<-1:Nd[d])
    gammat=alpha+rowSums(phi)
    
    P=max(abs(phi-phipast))
    G=max(abs(gammat-gammatpast))
    pas=max(P,G)
    
    phipast=phi
    gammatpast=gammat
  }
  return(list(val1=phi,val2=gammat)) 
}



#Test algo 
#vd=algo(Nd,2,0.1,matrice,alpha,beta)#[1]
#algo(10,0.1,matrice,alpha,beta)[2]

#fonction qui renvoit la taille d'un doc et qui prend en entr?e un doc (=ligne du corpus)
taille_doc<-function(matrice)
{ 
  Nd<-c()
  M=nrow(matrice)
  L=ncol(matrice)
  for (d in 1:M)
  {
    i=1
    
    while (matrice[d,i]!=0 & i<L)
    {
      i=i+1  
    }
    Nd[d]=i-1 
  }
  return(Nd)
}
#Test matrice 
#Nd=taille_doc(data_corpus_1000)

#R?cup?ration des phi et gamma pour chaque document

parcours_corpus<-function(k,matrice,epsilon,alpha,beta)
{
  #M nombre de document dans le corpus 
  M=nrow(matrice)
  Nd=taille_doc(matrice)
  algo_d_prime1=0
  #Nd vecteur des N pour chaque document
 
  phiglobal<-array(0:0, dim=c(k,max(Nd),M))
  gammatglobal=matrix(1:1,nrow=k,ncol=M) 
  for (d in 1:M)
  { 
    algo_d=algo(k,Nd[d],d,epsilon,matrice,alpha,beta)
    algo_d_prime1=algo_d$val1
    for (j in 1:ncol(algo_d_prime1))
    {
      phiglobal[,j,d]<-algo_d_prime1[,j]
    }
    
    gammatglobal[,d]=algo_d$val2
  }
  return(list(pc1=phiglobal,pc2=gammatglobal))
}

#Test Parcours Corpus
#phiglobal=parcours_corpus(matrice,0.1,alpha,beta)$pc1

#maximisation des beta pour phi et gamma fixes
maxbeta<-function(k,M,Nd,phiglobal,matrice2,V)
{
  beta=matrix(0:0,nrow=k,ncol=V)
  #tu es sur
  phi_int<-array(0:0, dim=c(k,max(Nd),M))
  for (i in 1:k)
  {
    for (j in 1:V)
    {
      for (d in 1:M)
      {
        for (n in 1:Nd[d])
        {
          phi_int[i,n,d]=symbole_kronecker(j,matrice2[d,n])*phiglobal[i,n,d]
        }
      }    
      beta[i,j]=sum(phi_int[i,,])
    }
    beta[i,]=beta[i,]/sum(beta[i,])
  }
  return(beta)
}

#Test Calcul Beta
#mb=maxbeta(M,Nd,phiglobal)

#maximisation des alpha pour phi et gamma fix?s: M?thode de Newton-Raphson

#introduction du symbole de kronecker
symbole_kronecker<- function(a,b)  
{
  if (a == b)
  {
    return(1)
  }
  else
  {
    return(0)
  }
}
#Test Kronecker
#symbole_kronecker(1,1000)

#calcul de la hessienne
calcul_hess<-function(k,alpha,M)
{
  hess=matrix(1:1,nrow=k,ncol=k)
  
  for (i in 1:k)
  {
    for (j in 1:k)
    {
      #formule ? justifier et v?rifier 
      #hess[i,j]=M*symbole_kronecker(i,j)*trigamma(alpha[i])-trigamma(sum(alpha))
      hess[i,j]=-M*(symbole_kronecker(i,j)*trigamma(alpha[i])-trigamma(sum(alpha)))
    }
  }
  return(hess)
}

#test
#alpha=c(rep(1,k))
#M=10
#calcul_hess(alpha,M)

#calcul du gradient
calcul_grad<-function(k,alpha,M,gammatglobal)
{
  grad=c()
  for (i in 1:k)
  {
    #grad[i]=M*(digamma(sum(alpha))-digamma(alpha[i]))+sum(digamma(gammatglobal[,i])-digamma(sum(gammatglobal))) 
    grad[i]=M*(digamma(sum(alpha))-digamma(alpha[i]))+sum(digamma(gammatglobal[i,]))-somme_doc(gammatglobal,M,k) 
  }  
  return(grad)
}

somme_doc<-function(gammat_global,M,k)
{
  somme=0
  for (d in 1:M)
  {
    somme=somme+digamma(sum(gammat_global[,d]))
  }
  return(somme)
}

#test
#gammatglobal=parcours_corpus(matrice,0.1,alpha,beta)$pc2
#print(gammatglobal)
#calcul_grad(alpha,M,gammatglobal)

#hessienneM=solve(calcul_hess(alpha,M))%*%calcul_grad(alpha,M,gammatglobal)

#Newton-Raphson

vect_true<-function(v){
  resu<-TRUE
  i=1
  n=length(v)
  while (resu==TRUE & i<(n+1)){
    if (v[i]<=0){resu<-FALSE}
    i=i+1
  }
  return(resu)
}


Newton_Raphson<-function(M,gammatglobal,epsilon,alpha,k)
{
  print("int normal")
  alphaold=c(rep(1,k))
  pas=10
  #pb de positivit? ? v?rifier 
  while (pas>epsilon & vect_true(alpha))
  {
    print("Hessieneok")
    print(calcul_hess(k,alpha,M))
    print("gradient")
    print(calcul_grad(k,alpha,M,gammatglobal))
    hess_1xgrad=solve(calcul_hess(k,alpha,M))%*%calcul_grad(k,alpha,M,gammatglobal)
    alpha=alpha-hess_1xgrad
    pas=max(abs(hess_1xgrad))
    
    if (vect_true(alpha))
    {alphaold=alpha}
    
  }
  #print("?tape calcul du alpha r?alis?")
  return(alphaold)
}