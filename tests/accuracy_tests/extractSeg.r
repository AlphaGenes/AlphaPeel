

subset = floor(seq(1, 1000, length.out = 200))
subset = c(1, subset + 1)

seg = read.table("outputs/peeling.seg")

write.table(seg[,subset], "outputs/seg.subset.txt", row.names=F, col.names=F, quote=F)