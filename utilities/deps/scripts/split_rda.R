library(seAMLess)
args = commandArgs(trailingOnly = TRUE)
rda_file = args[1]
expr_file = args[2]
meta_file = args[3]

load(rda_file)

exprs = Biobase::exprs(scRef)
write.csv(exprs, expr_file, row.names=TRUE, quote=TRUE)

meta = Biobase::pData(scRef)
write.csv(meta, meta_file, row.names = TRUE, quote = TRUE)
