
#calcul de la vraisemblance, divise en 5 parties: (pas termin?)


#calcul de la partie 1:

part1<-function(alpha,gammatglobal,Nd,d)
{
  res_part1= log(gamma(sum(alpha)))-sum(log(gamma(alpha)))+sum((alpha-1)*(digamma(gammatglobal[,d])-digamma(sum(gammatglobal[,d]))))
  return(res_part1)
}
#part1(alpha_test,gammat_test,Nd,456)



#calcul de la partie 2:

part2<-function(k,phiglobal,gammatglobal,Nd,d)
{
  a=c()
  for (n in 1:Nd[d])
  {
    a[n]=sum(phiglobal[,n,d]*(sapply(gammatglobal[,d],digamma)[]-digamma(sum(gammatglobal[,d])))) 
  }
  res_part2=sum(a)
  return(res_part2)
}

#part2(2,phi_test,gammat_test,Nd,2)

#calcul de la partie 3:

part3<-function(phiglobal,matrice2,beta,Nd,d,V)
{
  b=matrix(0:0,nrow=V,ncol=Nd[d])
  for (j in 1:V)
  {
    for (n in 1:Nd[d])
    {
      b[j,n]=sum(phiglobal[,n,d]*symbole_kronecker(j,matrice2[d,n])*log(beta)) 
    }
    
  }
  res_part3=sum(rowSums(b))
  return(res_part3)
}

#part3(phi_test,data_corpus_1000,beta_test,Nd,2,V)


#calcul de la partie 4:
b=matrix(1:1,nrow=4,ncol=5)
b[,1]-sum(b[,1])
part4<-function(gammatglobal,Nd,d)
{
  res_part4=lgamma(sum(gammatglobal[,d]))-sum(lgamma(gammatglobal[,d]))+sum((gammatglobal[,d]-1)*(digamma(gammatglobal[,d])-digamma(sum(gammatglobal[,d]))))
  return(res_part4)
}

#part4(gammat_test,Nd,2)
#gammatglobal=gammat_test

#calcul de la partie 5:

part5<-function(k,phiglobal,Nd,d)
{
  phiglobal2=matrix(0:0,nrow=k,ncol=Nd[d])
  phiglobal2<-phiglobal[,1:Nd[d],d]
  res_part5=sum(rowSums(phiglobal2*apply(phiglobal2,c(1,2),log)))
  return(res_part5)
}


#part5(k,phi_test,Nd,150)
#phiglobal=phi_test

#calcul de la vraisemblance totale

log_likelihood<-function(k,matrice,phiglobal,gammatglobal,alpha,beta,Nd,V,M)
{
  res_log_likelihood=c()
  for (d in 1:M)
  {
    res_log_likelihood[d]=part1(alpha,gammatglobal,Nd,d)+part2(k,phiglobal,gammatglobal,Nd,d)+part3(phiglobal,matrice,beta,Nd,d,V)-part4(gammatglobal,Nd,d)-part5(k,phiglobal,Nd,d)
  }
  res_total=sum(res_log_likelihood)
  return(res_total)
}

log_likelihood(k,matrice_corpus_test_25_75_50,phi_test,gammat_test,alpha_test,beta_test,Nd,V,M)

