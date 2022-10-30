#############################################################################################################
###The function yields number of generations n a minor outbreak (R<1) would involved before extinction   ####
###assuming the R and k were constant across generations, there is >99% of chance                        ####
###a cluster will go extinct within n generations                                                             ####
#############################################################################################################

#paramters:
#R0: reproduction number
#k: dispersion paramter
#--#
time_to_extinction=function(R0,k){
  gx=function(R0,k){
    (1+(1-g)*R0/k)^(-k)
  }
  g=0
  generation=1
  while(g<=1){
    g=gx(R0,k)
    if (g>0.99) {
      return(generation)
      break
    }
    generation=generation+1
  }
}

#--functions to generate Fig. S3--#
R0 = seq(0.5,0.95,length.out=700)
k = 10^seq(-2,1,0.01)
P = matrix(NA,nrow=length(R0),ncol=length(k))
#-#
for(i in 1:length(R0)) {
  for(j in 1:length(k)) {
    P[i,j] = time_to_extinction(R0[i],k[j])
  }
}
#-#
l.label = c(1,2,5,9,15,21,27,33,39)
#-#
image(R0,k,P,xlab=bquote("Reproduction number" ~ italic("R")),
      ylab=bquote("Dispersion parameter" ~ italic("k")),yaxt="n",
      log="y",col=hcl.colors(length(l.label)-1, "Heat 2", rev = TRUE),ylim = c(min(k),max(k)))
contour(R0,k,P,levels=l.label,labcex=1.6,add=TRUE,size=0.9,lty="dashed")
axis(2,at=c(0.01,0.1,1,10),labels=c("0.01","0.1","1","10") ,las=2)
#-#


