#set.seed(18)
library("coda")
library("MASS")
library("MCMCpack")
install.packages("Matrix")
library("Matrix")
library("rootSolve")



print(matrice)
nrow(matrice_corpus_test_25_75_50)

#On g?n?re une matrice de taille k*V, avec k le nombre de topics et V le nombre
#de mots dans le dictionnaire

#initialisations des dimensions et des param?tres 
k=2
V=max(data_corpus_1000)#Nombre de mot dans le dictionnaire propre au corpus
M=nrow(data_corpus_1000)#Nombre de doc dans le Corpus
eta=c(rep(1,V)) 
alpha=c(rep(10,k)) 

beta=rdirichlet(k,rep(1,V))

#initialisations des dimensions et des param?tres pour les donn?es simul?es 
k=2
V=max(matrice_corpus_test_1_3)#Nombre de mot dans le dictionnaire propre au corpus
M=nrow(matrice_corpus_test_1_3)#Nombre de doc dans le Corpus
eta=c(rep(1,V)) 
alpha=c(rep(10,k)) 

beta=rdirichlet(k,rep(1,V))

#Calcul du Phi et du Gamma pour chaque doc (prend en entr? le corpus,k, le num?ro du doc d, alpha et beta)
#passer le doc en argument 
algo<-function(k,Nd,d,epsilon,matrice,alpha,beta)
{ 
  pas=10
  #initialisations des param?tres variationnels
  phipast=(1/k)*matrix(1:1,nrow=k,ncol=Nd[d])
  phi=(1/k)*matrix(1:1,nrow=k,ncol=Nd[d])
  gammatpast=c()
  gammat=c()
  
  #print(beta)
  #S pour la normalisations 
  S=c()
  S[1]=1
  
  P=10
  G=10
  
  #initialisations des gammas
  
  
  gammatpast=alpha+Nd[d]/k
  gammat=alpha+Nd[d]/k
  
  
  while(epsilon<pas)
  {
    for (n in 1:Nd[d])
    {
      #calcul de phi ? n fix?
      for (i in 1:k)
      {
        #remplacer le i dans le beta par lindice du mot 
        indice_mot=matrice[d,n]
        phi[i,n]=beta[i,indice_mot]*exp(digamma(gammat[i]))
      }
      #normalisation du phi ? n fix? 
      
      S[n]=colSums(phi)[n]
      
      phi[,n]=phi[,n]/S[n]
      
    }
    
    #calcul de gammat
    for (i in 1:k)
    {
      gammat[i]=alpha[i]+sum(phi[i,n],n<-1:Nd[d])
    }
    
    P=max(abs(phi-phipast))
    G=max(abs(gammat-gammatpast))
    pas=max(P,G)
    
    phipast=phi
    gammatpast=gammat
  }
  return(list(val1=phi,val2=gammat)) 
}

#Test algo 
#algo(10,Nd,1000,0.1,data_corpus_1000,alpha,beta)


 
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
Nd=taille_doc(data_corpus_1000)


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
    algo_d=algo(k,Nd,d,epsilon,matrice,alpha,beta)
    algo_d_prime1=algo_d$val1
    for (j in 1:ncol(algo_d_prime1))
      {
        phiglobal[,j,d]<-algo_d_prime1[,j]
      }
    
    gammatglobal[,d]=algo_d$val2
  }
  #print("?tape calcul des phi et gammat corpus r?alis?")
  return(list(pc1=phiglobal,pc2=gammatglobal))
}

#Test Parcours Corpus
#phiglobal=parcours_corpus(10,data_corpus_1000,0.1,alpha,beta)$pc1
#gammatglobal=parcours_corpus(10,data_corpus_1000,0.1,alpha,beta)$pc2


#maximisation des beta pour phi et gamma fixes
maxbeta<-function(k,M,Nd,phiglobal,matrice2)
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
          #print("phi_int du topic i du doc d du mot n")
        #  print(phi_int[i,n,d])
          #print("kronoker")
          #print(symbole_kronecker(j,matrice2[d,n]))
          }
        }    
      beta[i,j]=sum(phi_int[i,,])
    }
    beta[i,]=beta[i,]/sum(beta[i,])
   # print("beta du doc topic i")
  #  print(beta[i,])
  }
#  print("?tape calcul du beta r?alis?")
  return(beta)
}

#Test Calcul Beta
#mb=maxbeta(k,M,Nd,phiglobal,data_corpus_1000)

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
      hess[i,j]=M*symbole_kronecker(i,j)*trigamma(alpha[i])-trigamma(sum(alpha))
    }
  }
  return(hess)
}

#test
#alpha=c(rep(1,k))
#M=10
#calcul_hess(alpha,M)

#calcul du gradient

somme_doc<-function(gammat_global,M,k)
{
  somme=0
  for (d in 1:M)
  {
    somme=somme+digamma(sum(gammat_global[,d]))
  }
  return(somme)
}
calcul_grad<-function(k,alpha,M,gammatglobal)
{
  grad=c()
  for (i in 1:k)
  {
    grad[i]=M*(digamma(sum(alpha))-digamma(alpha[i]))+sum(digamma(gammatglobal[i,]))-somme_doc(gammatglobal,M,k)
    #grad[i]=M*(digamma(sum(alpha))-digamma(alpha[i]))+sum(digamma(gammatglobal[,i])-digamma(sum(gammatglobal))) 
  }  
  return(grad)
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
  alphaold=c(rep(1,k))
  pas=10
  
  #pb de positivit? ? v?rifier 
  while (pas>epsilon & vect_true(alpha))
  {
    alpha=alpha-solve(calcul_hess(k,alpha,M))%*%calcul_grad(k,alpha,M,gammatglobal)
    pas=max(alpha-alphaold)
    
    if (vect_true(alpha))
    {alphaold=alpha}
    
  }
  print("?tape calcul du alpha r?alis?")
  return(alphaold)
}
#Test Newton Raphson
#Newton_Raphson(M,gammatglobal,0.1,alpha,k)

#Boucle globale, calcul de l'inf?rence:

Inference<-function(matrice2,epsilon,k)
{
  #? ajouter : initialisation alpha beta dedans  
  M=nrow(matrice2)
  Nd=taille_doc(matrice2)
  pas=1
  beta_old=rdirichlet(k,rep(1,V))
  alpha_old=c(rep(2,k)) 
  n=0
  vrais=list()
  while(pas>epsilon)
  {
    #print("etape1 gamma et phi")
    #calcul des nouveaus gamma et phi en partant des alpha et beta et des gamma et phi pr?cedent
    gammatglobal=parcours_corpus(k,matrice2,epsilon,alpha,beta)$pc2#peut-on optimiser ?a ?? 
    phiglobal=parcours_corpus(k,matrice2,epsilon,alpha,beta)$pc1
    vrais<-c(log_likelihood(k,matrice2,phiglobal,gammatglobal,alpha,beta,Nd,V,M),vrais)
    #calcul des nouveau alpha et beta en partant des alpha et beta pr?cedent 
    beta=maxbeta(k,M,Nd,phiglobal,matrice2)#corriger le beta !!!! 
    #print("etape3 Newton ")
    alpha=Newton_Raphson(M,gammatglobal,epsilon,alpha,k)
    vrais<-c(log_likelihood(k,matrice2,phiglobal,gammatglobal,alpha,beta,Nd,V,M),vrais)
    
    #print("etape4 convergence alpha beta")
    pas=max(abs(beta-beta_old),abs(alpha-alpha_old))
    #print(pas)
    beta_old=beta
    alpha_old=alpha
    n=n+1
    
  }
  return(list(value1=alpha,value2=beta,value3=phiglobal,value4=gammatglobal,value5=n,value6=vrais))
  
}



#Efficacité de l'algo
#Comparaison des Beta
Comparaison_Beta<-function(beta_0,beta_fin)
  {
  return(beta_0-beta_fin)
  }
Comparaison_Beta(beta_0_25_75_Test_2,beta_test2_25_75_50)


