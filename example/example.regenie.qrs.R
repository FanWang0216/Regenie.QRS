# example of Regenie.QRS GWAS
library(data.table)
library(dplyr)
library(quantreg)
library(seqminer)
options(scipen = 999)

## ACAT (Cauchy combination)
## reference => https://pubmed.ncbi.nlm.nih.gov/33012899/
## source code => https://github.com/yaowuliu/ACAT
cauchy.meta <- function(pvals) {
  # Check input
  pvals = pvals[!is.na(pvals)]
  if (length(pvals) == 0) {
    return(NA)
  }
  # pvals[pvals == 0] = 2.2E-308
  # Convert to Cauchy
  cauchy = 1 / (pvals * pi)
  cauchy[pvals >= 1e-15] = tanpi(0.5 - pvals[pvals >= 1e-15])
  stats = mean(cauchy)
  p = pcauchy(q = stats, lower.tail = F)
  return(p)
}
## QRank
## reference => https://pubmed.ncbi.nlm.nih.gov/28334222/
## source code => https://CRAN.R-project.org/package=QRank
## Quantile regression
## reference => https://doi.org/10.1201/9781315120256
## source code => https://cran.r-project.org/package=quantreg

QRank_proj <- function(gene, snp, tau = c(0.25, 0.5, 0.75)){
  ltau = length(tau)
  x = as.matrix(snp)
  y = as.matrix(gene)
  N=length(y)
  VN = matrix(0, nrow = ltau, ncol = ltau)
  for(i in 1:ltau) {
    for(j in 1:ltau) {
      VN[i,j] = min(tau[i], tau[j]) - tau[i] * tau[j]
    }
  }
  
  SN = NULL
  for(i in 1:ltau) {
    qreg=list()
    qreg$coefficients=quantile(y, probs = tau[i])
    qreg$residual = y - rep(1,N)%*% as.matrix(qreg$coefficients)
    qreg$dual = 1 * (qreg$residual > 0)
    ranks = qreg$dual - (1 - tau[i])
    Sn = as.matrix(t(x) %*% (ranks))
    SN = c(SN,Sn)
  }
  
  VN2 = matrix(outer(VN,t(x)%*% x,"*"), nrow = ltau)
  pvalue1 = pchisq(SN^2/diag(VN2), 1, lower.tail = F)
  names(pvalue1) = tau
  chi.1=SN^2/diag(VN2)
  
  e = solve(chol(VN2))
  SN2 = t(e) %*% SN
  pvalue = pchisq(sum(SN2^2), ltau, lower.tail = F)
  composite.chi=sum(SN2^2)
  
  result = list(composite.pvalue = pvalue, quantile.specific.pvalue = pvalue1, tau = tau, composite.chi=composite.chi, chi.quantile.specific=chi.1,quantile.specific.beta=SN,quantile.specific.w=(1/diag(VN2)))
  return(result)}

####################################################################################################

## quantile levels
## 0 < values < 1
qntl = (1:9) / 10

set.seed(seed = 256)
data.path <- "/path/for/data/example/"
data.path='/Users/macbook/work/error-in-variable/UKBB_quantile/example/'

y.file=fread(paste0('',data.path,'normalized_pheno.txt'))
covariate=fread(paste0('',data.path,'covariate.txt'))
X_covar=as.matrix(covariate[,-c(1,2)])
y=as.numeric(unlist(y.file[,3]))
y_proj <- lm(y ~ X_covar)$residuals
y_proj.sd=scale(y_proj,center=F,scale=T)

###Run Regenie 
Regenie.path <- "/path/for/regenie"
cmd <- paste(
  paste0(Regenie.path, "regenie"),
  "--step 1",
  paste0("--bed ", data.path, "example"),
  paste0("--covarFile ", data.path, "covariate.txt"),
  paste0("--phenoFile ", data.path, "normalized_pheno.txt"),
  "--bsize 1000",
  "--lowmem",
  paste0("--out ", data.path, "fit_bin_out"),
  sep = " \\\n"
)
cat("Running command:\n", cmd, "\n")
system(cmd)

##Extract polygenic predictions from the Regenie output
y_loco.df=as.matrix(fread(paste0('',data.path,'fit_bin_out_1.loco')))[,-1]
N=length(y)
bim_file=fread(paste0('',data.path,'example.bim'))
M=dim(bim_file)[1]
#check if the colnames ids from regenie output matches with the y ids
#y_ids=y.file[,1]
#reordered_indices <- match(paste(y_ids, y_ids, sep = "_"), colnames(y_loco.df))

#load the genotype matrix
names.ukbdata=paste0('',data.path,'example')
UKBdata <- names.ukbdata
geno.chr<- readPlinkToMatrixByIndex(UKBdata,1:N, 1:M)
standardize_genotype <- function(SNP.set) {
  AF<-apply(SNP.set,2,mean)/2
  mean_value <- 2 * AF
  sd_value <- sqrt(2 * AF * (1 - AF))
  standardized_matrix <- scale(SNP.set, center = mean_value, scale = sd_value)
  return(standardized_matrix)
}
geno.chr.scale=standardize_genotype(geno.chr)

##Regenie.QRS test
is.effect.estimated = F
df.p2 <- data.frame()
for (i in 1:10) {
  print(i)
  # Inputs
  snp <- geno.chr.scale[, i]
  chr <- as.numeric(bim_file[i, 1])
  y_loco_rgen <- y_loco.df[chr, ]
  y_resid_rgen <- y_proj.sd - y_loco_rgen
  
  # Residualize SNP
  G_proj <- scale(lm(snp ~ X_covar)$residuals, center = FALSE, scale = TRUE)

  # Run Regenie.QRS test
  p_qr <- QRank_proj(gene = y_resid_rgen, snp = G_proj, tau = qntl)$quantile.specific.pvalue
  
  beta_qr=double()
  if (is.effect.estimated){
  # Get quantile-specific effect sizes
  beta_qr <- sapply(qntl, function(tau) {
    qreg <- summary(rq(y_resid_rgen ~ G_proj, tau = tau), se = "iid")
    as.numeric(qreg$coefficients[2, 1])
  })
  }
  row_df <- data.frame(t(p_qr), t(beta_qr))
  df.p2 <- rbind(df.p2, row_df)
}
df.p2 = df.p2 %>% mutate(`P QR` = apply(X = df.p2[,c(1:length(qntl))], MARGIN = 1, FUN = cauchy.meta), .before = 1)

colnames(df.p2) = c("P_QR",
                    paste0("P_QR", qntl),
                    paste0("BETA_QR", qntl))%>% head(n = 1 + length(qntl) + length(beta_qr))
snp_ids <- as.character(bim_file[[2]])  
df.p2 <- df.p2 %>%
  mutate(ID = snp_ids[seq_len(nrow(df.p2))]) %>%
  select(ID, everything())

## Regenie.QRS summary statistics
## row: variant
## column: p-value # and beta
fwrite(x = df.p2, file = paste0('',data.path,'/example.sumstat.tsv'), quote = F, sep = "\t", row.names = F, col.names = T)

