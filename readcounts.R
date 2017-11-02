readcounts = read.table("ggg.readcounts",col.names=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11"))

nsamples = nrow(readcounts)
temp = c(0,1)
par(xaxs="i")
par(cex.axis=1.0)

bitmap("ggg.readcounts.jpg",type="jpeg",width=4,height=0.8*(nsamples+2.0),res=600,pointsize=10)

par(mar=c(0,0,0,0))
plot(temp,temp,type="n",ylab="",xlab="",main="",ylim=c(0.18,nsamples+2.0+0.18),xlim=c(-0.1,1.1),axes=FALSE)

#All reads
rect(0,c(1:nsamples)-0.8,readcounts$V2/readcounts$V2,c(1:nsamples)-0.2)

#Trimmed reads
rect(0,c(1:nsamples)-0.8,readcounts$V3/readcounts$V2,c(1:nsamples)-0.2,col=colors()[19])

base = 0
#1:$line 2:$total 3:$clipped 4:$rRNA_fwd 5:$rRNA_rev 6:$ncRNA_fwd 7:$ncRNA_rev
#      8:$mRNA_fwd 9:$mRNA_rev 10:$vRNA_fwd 11:$vRNA_rev

#vRNA fwd+rev
fwdrev = readcounts$V10+readcounts$V11
rect(base,c(1:nsamples)-0.8,base+fwdrev/readcounts$V2,c(1:nsamples)-0.2,col=colors()[552])
revfrac = readcounts$V11/fwdrev
rect(base,c(1:nsamples)-0.8,base+fwdrev/readcounts$V2,c(1:nsamples)-0.8+0.6*revfrac,col=colors()[315])
base = base + fwdrev/readcounts$V2

#mRNA fwd+rev
fwdrev = readcounts$V8+readcounts$V9
rect(base,c(1:nsamples)-0.8,base+fwdrev/readcounts$V2,c(1:nsamples)-0.2,col=colors()[494])
revfrac = readcounts$V9/fwdrev
rect(base,c(1:nsamples)-0.8,base+fwdrev/readcounts$V2,c(1:nsamples)-0.8+0.6*revfrac,col=colors()[315])
base = base + fwdrev/readcounts$V2

#ncRNA_fwd
fwdrev = readcounts$V6+readcounts$V7
rect(base,c(1:nsamples)-0.8,base+fwdrev/readcounts$V2,c(1:nsamples)-0.2,col=colors()[79])
revfrac = readcounts$V7/fwdrev
rect(base,c(1:nsamples)-0.8,base+fwdrev/readcounts$V2,c(1:nsamples)-0.8+0.6*revfrac,col=colors()[315])
base = base + fwdrev/readcounts$V2

#rRNA_fwd
fwdrev = readcounts$V4+readcounts$V5
rect(base,c(1:nsamples)-0.8,base+fwdrev/readcounts$V2,c(1:nsamples)-0.2,col=colors()[616])
revfrac = readcounts$V5/fwdrev
rect(base,c(1:nsamples)-0.8,base+fwdrev/readcounts$V2,c(1:nsamples)-0.8+0.6*revfrac,col=colors()[315])
base = base + fwdrev/readcounts$V2

text(0.5,c(1:nsamples)-0.1,labels=readcounts$V1,adj=0.5)
text(0.5,nsamples+2.1,labels="ggg",adj=0.5,cex=1.4)

#Legend
x1 = 0.05

rect(0.0+x1,nsamples+1.8,0.1+x1,nsamples+1.5,col="white")
rect(0.0+x1,nsamples+1.4,0.1+x1,nsamples+1.1,col=colors()[19])
rect(0.0+x1,nsamples+1.0,0.1+x1,nsamples+0.7,col=colors()[616])
rect(0.0+x1,nsamples+0.6,0.1+x1,nsamples+0.3,col=colors()[315])
rect(0.60+x1,nsamples+1.8,0.70+x1,nsamples+1.5,col=colors()[79])
rect(0.60+x1,nsamples+1.4,0.70+x1,nsamples+1.1,col=colors()[494])
rect(0.60+x1,nsamples+1.0,0.70+x1,nsamples+0.7,col=colors()[552])

text(0.12+x1,nsamples+1.65,adj=0.0,labels="untrimmed reads")
text(0.12+x1,nsamples+1.25,adj=0.0,labels="other trimmed reads")
text(0.12+x1,nsamples+0.85,adj=0.0,labels="rRNA")
text(0.12+x1,nsamples+0.45,adj=0.0,labels="reverse sense matches")
text(0.72+x1,nsamples+1.65,adj=0.0,labels="ncRNA")
text(0.72+x1,nsamples+1.25,adj=0.0,labels="mRNA")
text(0.72+x1,nsamples+0.85,adj=0.0,labels="vRNA")

dev.off()
