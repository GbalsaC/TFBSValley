setwd("B:/Documents/R")
signal<-as.numeric(read.csv("aa.csv", header=F)[1,])
### Requires multiplot.R, Rcpp, ggplot

RamseyR <- function(signal){
        index=1;
        ##Packages, sourceCpp has rebuild False
        library(ggplot2);
        library("dplyr");
        library(Rcpp);
        library("reshape2");
        sourceCpp("LemireR.cpp");
        source('~/R/functions/multiplot.R');
        
        ## Default Values and Parameters
        Window=100;
        Offset = 10;
        Depth = 0.9;
        sigma=6;
        sigma2=exp((sigma-1)*log(2)/2);
        X=qnorm(1e-3,0,sigma2);
        window2=dnorm(round(X):-round(X),0,sigma2);
        wSize=length(window2);
        ones=vector(mode="double",length=wSize);
        for(i in 1:length(ones)){ones[i]<-1;}
        vector=ones/wSize;
        
        ##Cloning variables for modification individually 
        ##is memory intensive, and should not duplicate memory usage
        M = convolve(vector, rev(as.numeric(signal)), type="open")[ceiling(wSize/2):ceiling(length(signal)-1+(wSize/2))];
        newWin = Window-Offset;
        
        #**** Ends valeyscoring 1st set
        #**** Starts minmaxfilt1 1st set        
        
        M = lemireR(a=M, w=newWin+1, offset=10, Depth=0.8); #Rewrite on M for better memory usage
        #lemireR saves both minval and maxval, apparently only maxval is necesarry
        # but it's possible V(pos) is a version of minval with less window
        
        length <- dim(M)[1];
        M<-cbind(c(index:(index-1+length)),M);
        colnames(M)<-c('Index','MD', 'Vpos', 'Vlog', 'Begin', 'End');
        #### Using Melt format just for Plotting
        M2<-melt(as.data.frame(M), id="Index", measure=c("MD", "Vpos", "Vlog"));
        ggplot(data=M2, aes(x=Index, y=value, colour=variable))+geom_line();
        ggsave(file="ramseyR.png");
        
        ### Using a for loop to get index begin and end
        count=0;
        start=0;
        for(i in 1:length){
                if(M[i,'Vlog']>0)
                {
                        if(start==0){start=M[i,'Index'];}
                        M[i,'Begin']=start;
                        count=count+1;        
                        
                }else
                {
                        if(start != 0)
                        {
                                ### Then it's the first time entering 0
                                #### All elements from start to actual of the End column
                                #### should save this index value
                                M[start:(M[i,'Index']-1),'End']=M[i,'Index'];
                                start=0;
                                count=0;
                        }
                        ##Do nothing else
                }

        }
        ## Converting to Data Frame allows factor checking
        M<-as.data.frame(M);
        M
}