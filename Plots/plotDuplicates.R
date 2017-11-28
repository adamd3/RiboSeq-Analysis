plotDuplicates = function(TABLE, JPG){
  tab = read.csv(TABLE, sep="\t", header=FALSE)
  tab$V1 = unlist(lapply(tab$V1, function (x) gsub("^\\s+", "", x)))
  tab$freq = as.numeric(unlist(lapply(tab$V1, function (x) unlist(strsplit(x, split=" ")[[1]][1]))))
  tab$ndups = as.numeric(unlist(lapply(tab$V1, function (x) unlist(strsplit(x, split=" ")[[1]][2]))))
  tab = tab[-1]
  tab = tab[order(tab$ndups),]
  r = 1:max(tab$ndups)
  r = r[!r %in% tab$ndups]
  r = cbind(rep(0, length(r)), r)
  colnames(r) = c("freq", "ndups")
  tab = rbind(tab, r)
  tab = tab[order(tab$ndups),]
  m = ceiling(log10(max(tab$freq)))
  jpeg(JPG)
  barplot(log10(tab$freq + 1), col="blue", xlim=c(0, nrow(tab)), ylim=c(0, m), axes=F,
          xlab="Number of Copies", ylab="Frequency")
  labs = seq(0, m, 1)
  axis(2, at=labs, labels=unlist(lapply(labs, function(x) bquote(10^.(x)))))
  axis(1, at=seq(0, nrow(tab)+10, 10))
  dev.off()
}
args = commandArgs(trailingOnly=TRUE)
TABLE = args[1]
JPG = args[2]
plotDuplicates(TABLE, JPG)