args <- commandArgs(TRUE)
transcripts <- read.table(args[1])

aa <- transcripts[order(transcripts$V1, -abs(transcripts$V3) ), ]

write.table( aa[ !duplicated(aa$V1), ] ,file=args[2],row.names=F,col.names=F,quote=F)
