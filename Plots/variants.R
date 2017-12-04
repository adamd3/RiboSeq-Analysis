plotVariants <- function(PILEUP, GTF){
  filterRead <- function(string, reference){
    # fills in the "bases" column of the table with the appropriate reference base
    string = gsub("[$|^|~|*]", "", string)
    string = as.vector(unlist(gsub('[.]', as.character(reference), string)))
    string = as.vector(unlist(gsub('[,]', tolower(as.character(reference)), string)))
    
    # splits this string into characetrs
    res = strsplit(string, "")
    
    # counts how many of each base are present and returns this table
    return (table(res))
  }
  
  # read the pileup table
  tab = read.csv(PILEUP, sep="\t", header=F)
  colnames(tab) = c("Reference", "Position", "Ref_Base", "Read_Count", "Bases", "Quality")
  
  # count the number of each base at each position
  bases = mapply(filterRead, tab$Bases, tab$Ref_Base)
  
  options = c("A", "C", "T", "G", "a", "c", "t", "g")
  rtab = data.frame()
  rtab2 = data.frame()
  i = 1
  # fill in a table with the counts of all bases at all positions 
  # and a table with the proportions of all bases at all positions
  for (b in bases){
    for (o in options){
      rtab[i, o] = unlist(b[o])
      rtab2[i, o] = unlist(b[o]) / sum(b)
    }
    i = i + 1
  }
  rtab2[is.na(rtab2)] = 0
  rtab[is.na(rtab)] = 0
  i = 1
  
  # calculates the proportion of reference bases at each position
  xs = c()
  for (base in tab$Ref_Base){
    row = rtab2[i,]
    if (sum(row) != 0){
      xs = c(xs, as.numeric(row[base]))
    }
    else {
      xs = c(xs, -1)
    }
    i = i + 1
  }
  tab$xs = xs
  
  pdf(gsub(".tsv", ".jpg", PILEUP))
  
  # reads a gtf showing where the genes are in the virus
  gtf = read.table(GTF, sep="\t", header=F, comment.char="", quote="")
  nbases = nrow(tab)
  nchunks = 4
  plot.new()
  par(mfrow=c(nchunks*3, 1), mar=c(1, 1, 1, 1), bty='n', xpd=F)
  n_seg = ceiling(nbases / nchunks)
  j = 0
  for (i in 1:nchunks){
    sub_xs = xs[j:(j + n_seg)]
    sub_readcounts = tab$Read_Count[j: (j + n_seg)]
    sub_pos = tab$Position[j: (j + n_seg)]
    plot(c(), ylim=c(0, nrow(gtf) + 2), xlim=c(j, (j+n_seg)), axes=F)
    k = 1
    for (g in 1:nrow(gtf)){
      row = gtf[g,]
      nam = row[]
      start = row$V4
      stop = row$V5
      startpos = sub_pos[1]
      stoppos = sub_pos[length(sub_pos)-1]
      
      if (start < startpos & stop > startpos){
        start = startpos
      }
      if (start < stoppos & stop > stoppos){
        stop = stoppos
      }
      s = unlist(lapply(row$V9, function(x) strsplit(as.character(x), split=";")))
      f = unlist(lapply(s, function(x) grepl("gene_name", x)))
      genename = unlist(strsplit(unlist(s)[f], split=" "))[3]
      genename = gsub("\"", "", genename)
      lines(c(start, stop), c(k, k), lw=2, col='lightgrey',ylim=c(0, nrow(gtf) + 1),
            xlim=c(j, (j+n_seg)))
      if (start >= startpos & start <= stoppos){
        text(start, k, genename)
      }
      k = k + 1
    }
    axis(1, pos=0, at=c(sub_pos[1], sub_pos[length(sub_pos)-1]), lwd.ticks = 0,
         labels=F)
    sub_xs = sub_xs[!is.na(sub_xs)]
    sub_readcounts = sub_readcounts[!is.na(sub_readcounts)]
    sub_pos = sub_pos[!is.na(sub_pos)]
    plot(sub_pos, sub_xs, ylim=c(-0.1, 1), pch=19, cex=0.5, xlim=c(j, (j+n_seg)),
         col="hotpink", axes=F)
    subsub_pos = sub_pos[sub_xs == 0]
    for (p in subsub_pos){
      lines(c(p, p), c(0,1), type='l', col='blue')
    }
    subsub_pos = sub_pos[sub_xs == -1]
    for (p in subsub_pos){
      lines(c(p, p), c(0,1), type='l', col='grey')
    }
    axis(1, pos=0, at=c(sub_pos[1], sub_pos[length(sub_pos)-1]), lwd.ticks = 0,
         labels=F)
    axis(2, pos=sub_pos[1], at=c(0, 1))
    x = sub_pos
    y = log2(sub_readcounts + 1)
    
    plot(x, y, type="l", xlim=c(j, (j+n_seg)), ylim=c(0, ceiling(max(log2(tab$Read_Count + 1)))),
         axes=F)
    polygon(c(min(x), x, max(x)), c(0, y, 0),  col = "grey")
    axis(1, pos=0, at=c(sub_pos[1], sub_pos[length(sub_pos)-1]))
    m = ceiling(max(log2(tab$Read_Count + 1)))
    axis(2, pos=sub_pos[1], at=seq(0, m, 4), labels=unlist(lapply(seq(0, m, 4),
                                                                  function(x) 2^x)), las=1)
    j = j + n_seg
  }
  dev.off()
  }
args = commandArgs(trailingOnly=TRUE)
PILEUP = args[1]
GTF = args[2]
plotVariants(PILEUP, GTF)
