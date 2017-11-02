bulk = read.table('../K562_CNV_from_ENCODE_liftOver_to_hg38.bed', stringsAsFactors = FALSE)
bulk[bulk$V4=='normal', ]
bulk$interpret = 0
for (i in 1:nrow(bulk)){
  if(bulk$V4[i]=='amp'){bulk$interpret[i]=4}
  if(bulk$V4[i]=='het.del'){bulk$interpret[i]=1}
  if(bulk$V4[i]=='homo.del'){bulk$interpret[i]=0}
  if(bulk$V4[i]=='normal'){bulk$interpret[i]=2}
}
(bulk[bulk$V1=='chr22', 'interpret'])
