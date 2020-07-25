
# Required packages need to be installed prior to run this code

library(e1071)
library(seqinr)
library(foreign)
library(caret)
library(Biostrings)
library(ranger)

# Input sequence given by the user in FASTA format

pdb=read.fasta(file.choose())

q=length(pdb)
z3=list()

for(p in 1:q)
{

z=getSequence(pdb[[p]])

# This stage will remove the his tag if present in the given input sequence

d=paste(z, collapse='' )
x=AAString(d)
g=countPattern('hhhhhh',x)
if(g >= 1)
{
f=matchPattern('hhhhhh',x)
z1=as.matrix(z)
z1=z1[-(start(f):end(f)),]
z1=as.matrix(z1)
}
if(g == 0)
{
z1=as.matrix(z)
}




z2=t(z1)

z2=toupper(z2)

for(i in 1:length(z2))
{
if(z2[1,i]=="F")
{
z2[1,i]<-"FF"
}
else if(z2[1,i]=="T")
{
z2[1,i]<-"TT"
}
}
z3[[p]]<-z2

}



# This stage will replace the normalized value of amino acid property descriptors with corresponding residues in given sequence

z4=list()

for(t in 1:q)
{
b<-as.data.frame(z3[[t]])
r=nrow(b)
c=ncol(b)
b1=read.csv("xxx/AAP SLAA(49P).csv") #Need to replace xxx with the path of file name 'AAP SLAA(49P).csv'
z=matrix(data=NA,nrow=nrow(b1),ncol=ncol(b))
for(k in 1:nrow(b1))
{
for(j in 1:c)
{
e1<-b[1,j]
e<-as.character(e1)
p=b1[k,e]
z[k,j]<-p
}
}

z4[[t]]<-z

}





y=matrix(data=NA,nrow=q,ncol=49)
for(w in 1:q)
{
c=as.data.frame(z4[[w]])
d=t(c)
for(i in 1:49)
{
e=mean(d[,i])
y[w,i]<-e
}
}

y=as.data.frame(y)


# Training the dataset we have collected using Random Forest classifier model and then it will predict the given sequence type


iris=read.csv("yyy/R1&others.csv") #Need to replace yyy with the path of file name 'R1&others.csv'
set.seed(123)
rfFit1 <- train(V50~ ., data=iris, method="ranger")
predictions <- predict(rfFit1, y[,1:49])
predictions


if(predictions=='RR')
{
iris=read.csv("yyy/R2&others.csv") #Need to replace yyy with the path of file name 'R2&others.csv'
set.seed(123)
rfFit1 <- train(V50~ ., data=iris, method="ranger")
predictions <- predict(rfFit1, y[,1:49])
predictions
}

if(predictions=='RR')
{
iris=read.csv("yyy/R3&others.csv") #Need to replace yyy with the path of file name 'R3&others.csv'
set.seed(123)
rfFit1 <- train(V50~ ., data=iris, method="ranger")
predictions <- predict(rfFit1, y[,1:49])
predictions
}

if(predictions=='RR')
{
iris=read.csv("yyy/R4&others.csv") #Need to replace yyy with the path of file name 'R4&others.csv'
set.seed(123)
rfFit1 <- train(V50~ ., data=iris, method="ranger")
predictions <- predict(rfFit1, y[,1:49])
predictions
}

if(predictions=='RR')
{
iris=read.csv("yyy/R5&others.csv") #Need to replace yyy with the path of file name 'R5&others.csv'
set.seed(123)
rfFit1 <- train(V50~ ., data=iris, method="ranger")
predictions <- predict(rfFit1, y[,1:49])
predictions
}

if(predictions=='RR')
{
iris=read.csv("yyy/R6&others.csv") #Need to replace yyy with the path of file name 'R6&others.csv'
set.seed(123)
rfFit1 <- train(V50~ ., data=iris, method="ranger")
predictions <- predict(rfFit1, y[,1:49])
predictions
}


if(predictions=='R1')
{
y='The given protein sequence is predicted as ARM repeat type solenoid protein'
}


if(predictions=='R2')
{
y='The given protein sequence is predicted as HEAT repeat type solenoid protein'
}


if(predictions=='R3')
{
y='The given protein sequence is predicted as TPR repeat type solenoid protein'
}


if(predictions=='R4')
{
y='The given protein sequence is predicted as BETA repeat type solenoid protein'
}


if(predictions=='R5')
{
y='The given protein sequence is predicted as ANK repeat type solenoid protein'
}


if(predictions=='R6')
{
y='The given protein sequence is predicted as LRR repeat type solenoid protein'
}


if(predictions=='RR')
{
y='The given protein sequence is predicted as non-solenoid protein / any other type of repeat containing solenoid protein '
}

# Corresponding results will be printed here

print(y)
