source("/Users/benjaminmuller/Desktop/Ensae/2A/Stat app/Fonction_intermédiaire.R")
source("/Users/benjaminmuller/Desktop/Ensae/2A/Stat app/LDA_loglikelyhood.R")


Inference<-function(matrice2,epsilon,k,V)
{
  
  M=nrow(matrice2)
  Nd=taille_doc(matrice2)
  #aa=rdirichlet(2,rep(1,6))
  beta_old=rdirichlet(k,rep(1,V))
  alpha_old=c(rep(1,k)) 
  
  #beta=rdirichlet(k,rep(1,V))
  #alpha=c(rep(1,k)) 
  alpha=c(rep(1,k)) 
  Initialized=Initialiazing(alpha,beta,M,epsilon,k,V,matrice2)
  alpha=Initialized$val1
  beta=Initialized$val2
  vraisemblance_init=Initialized$val3
    

  vraisemblance<-list()
  
    
  pas=1
  n=0#compteur du nombre d'itération 
  while(pas>epsilon)
  {
    #calcul des nouveaus gamma et phi en partant des alpha et beta et des gamma et phi pr?cedent
    phi_gamma=parcours_corpus(k,matrice2,epsilon,alpha,beta)
    gammatglobal=phi_gamma$pc2
    phiglobal=phi_gamma$pc1
    
    vraisemblance[2*n+1]=log_likelihood(k,matrice2,phiglobal,gammatglobal,alpha,beta,Nd,V,M)
    
    #calcul des nouveau alpha et beta en partant des alpha et beta precedent 
    beta=maxbeta(k,M,Nd,phiglobal,matrice2,V)
    print("alpha1")
    print(alpha)
    alpha=Newton_Raphson(M,gammatglobal,epsilon,alpha,k)
    print("alpha2")
    print(alpha)
    vraisemblance[2*n]=log_likelihood(k,matrice2,phiglobal,gammatglobal,alpha,beta,Nd,V,M)
    

    pas=max(abs(beta-beta_old),abs(alpha-alpha_old))
    beta_old=beta
    alpha_old=alpha
    n=n+1
    
  }
  return(list(alpha,beta,phiglobal,gammatglobal,n,vraisemblance,vraisemblance_init))
  
}
aa=runif(min=0,max=50,5)
b=floor(aa)
matrix(1:1,nrow=50,ncol=3)[b,]


Initialiazing<-function(alpha,beta,M,epsilon,k,V,matrice2)
{
  beta_old=rdirichlet(k,rep(1,V))
  alpha_old=c(rep(1,k)) 
  n=0
  vraisemblance_init=c()
  M_reduit=floor(M/10)
  indice_doc=runif(min=0,max=M,M_reduit)
  matrice_reduite=matrice2[indice_doc,]
  pas=1
  Nd_reduit=taille_doc(matrice_reduite)
  
  while(pas>epsilon)
  {
    #calcul des nouveaus gamma et phi en partant des alpha et beta et des gamma et phi pr?cedent
    phi_gamma=parcours_corpus(k,matrice_reduite,epsilon,alpha,beta)
    gammatglobal=phi_gamma$pc2
    phiglobal=phi_gamma$pc1
    
    vraisemblance_init[2*n+1]=log_likelihood(k,matrice_reduite,phiglobal,gammatglobal,alpha,beta,Nd_reduit,V,M_reduit)
    
    
    #calcul des nouveau alpha et beta en partant des alpha et beta precedent 
    beta=maxbeta(k,M_reduit,Nd_reduit,phiglobal,matrice_reduite,V)
    print(alpha)
    alpha=Newton_Raphson(M_reduit,gammatglobal,epsilon,alpha,k)
    print(alpha)
    vraisemblance_init[2*n]=log_likelihood(k,matrice_reduite,phiglobal,gammatglobal,alpha,beta,Nd_reduit,V,M_reduit)
    
    
    pas=max(abs(beta-beta_old),abs(alpha-alpha_old))
    beta_old=beta
    alpha_old=alpha
    n=n+1
    
  }
  print("initialized")
  print(alpha)
  print(beta)
return(list(val1=alpha,val2=beta,val3=vraisemblance_init))
  
}
