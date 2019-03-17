function(amat, etype="Directed"){
outmat<-matrix("E",1,4)
outmat[1,]<-c("Source","Target","Type","Weight")
nrows = dim(amat)[1]

items<-row.names(amat)

for(i1 in c(1:nrows)){
  sname<-items[i1]   #Name of source node
  targs<-items[amat[i1,]!=0]   #Names of target nodes with non-zero weights
  nitems<-length(targs)        #total number of target nodes for this source

   if(nitems > 0){             #Make matrix of out-going edges and add to main output matrix
      tmpmat<-matrix("E", nitems, 4)
      tmpmat[,1]<-sname
      tmpmat[,2]<-targs
      tmpmat[,3]<-etype
      tmpmat[,4]<-amat[i1,amat[i1,]!=0]
      outmat<-rbind(outmat,tmpmat)
      }

#  if(length(targs)>0){
#    for(i2 in c(1:nitems)){
#      wt<-amat[items==sname, items==targs[i2]]
#      newrow<-c(sname, targs[i2], "Undirected",wt)
#      outmat<-rbind(outmat, newrow)
#      }
#   }
  }
  outmat<-as.data.frame(outmat)   #Convert to data frame
  names(outmat)<-outmat[1,]
  outmat<-outmat[2:(dim(outmat)[1]-1),]
  outmat[,4]<-as.numeric(outmat[,4])
}

