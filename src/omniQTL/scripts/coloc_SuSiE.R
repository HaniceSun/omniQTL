library(coloc)

colocSuSiE = function(ss1, ss2, ld){
    out_file = gsub('.ld$', '_coloc.txt', ld)
    D1 = list()
    D2 = list()
    df1 = read.table(ss1, header=T)
    df2 = read.table(ss2, header=T)
    df_ld = read.table(ld, header=F)
    rownames(df_ld) = df1$snp
    colnames(df_ld) = df2$snp

    D1 = as.list(df1)
    D2 = as.list(df2)
    D1$LD = as.matrix(df_ld)
    D2$LD = as.matrix(df_ld)
    D1$N = D1$N[1]
    D2$N = D2$N[1]
    D1$type = D1$type[1]
    D2$type = D2$type[1]

    if(nrow(df_ld) > 0){ 
        res = try(
                  {
                      s1=runsusie(D1)
                      s2=runsusie(D2)
                      coloc.susie(s1,s2)
                  }, silent=T)

        if(!is.na(res[1])){
            cols = 'nsnps\thit1\thit2\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tidx1\tidx2'
            writeLines(cols, out_file)
            write.table(res$summary, out_file, row.names=F, col.names=F, sep='\t', quote=F, append=T)
        }
    }
}

args = commandArgs(trailingOnly=TRUE)
ss1 = args[1]
ss2 = args[2]
ld = args[3]
colocSuSiE(ss1, ss2, ld)
