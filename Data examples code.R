
###function to estimate the weight

estimate.weights.geom<-function(y,x, vaf)
{  #y=outcome, x=observed mutation, vaf=variant allele frequencies of observed mutations
  # unconditional mean for theta when theta is estimated conditionally on subjects with observed mutations
  
  
  uncondmean<-function(theta){
    #should be all thetas for observed mutations in this group
    1/(mean(1/theta,na.rm=T))
    
  }
  
  #unconditonal mean for cases and controls
  m.theta0.uncond <-uncondmean(vaf[x == 1 & y == 0])
  m.theta1.uncond <- uncondmean(vaf[x == 1 & y == 1])
  
  wt <- rep(NA,length(x))
  
  wt[x == 0 & y == 0] <- m.theta0.uncond
  wt[x== 0 & y == 1] <- m.theta1.uncond
  wt[x == 1] <- 1
  
  return(wt)
  
  
}


##bootstrap p-value

boot.p<-function(y,x,vaf, nsim=1000,fixedNcases=TRUE)
{#y=outcome, x=observed mutation, vaf=variant allele frequencies of observed mutations
  #fixedNcases is TRUE when exact number of cases and controls has to be preserved
 
  
  wt<-estimate.weights.geom(y,x,pmin(1,2*vaf))
  fit<-glm(y~x,family=binomial(),weight=wt)
  V_sand <- sandwich::vcovHC(fit, type = "HC0")
  
  beta_hat <- coef(fit)["x"]
  se_hat     <- sqrt(V_sand["x", "x"])
  se   <- sqrt(vcov(fit)["x", "x"])
  
  ci_beta <- beta_hat + c(-1, 1) * qnorm(0.975) * se
  
  T_obs <- beta_hat / se_hat
  dat<-as.data.frame(cbind(y,x,vaf))
  ref<-ts<-rep(NA,nsim)
  for (i in 1:nsim)
  {
    if (!fixedNcases)
    {ind <- sample(seq_len(nrow(dat)), replace = TRUE)
    dat_b <- dat[ind, , drop = FALSE]}
    else
    {ind0 <- which(dat$y == 0)
    ind1 <- which(dat$y == 1)
    
    boot_ind <- c(
      sample(ind0, size = length(ind0), replace = TRUE),
      sample(ind1, size = length(ind1), replace = TRUE)
    )
    
    dat_b <- dat[boot_ind, , drop = FALSE]
    }
    
    wt2<-estimate.weights.geom(dat_b$y,dat_b$x,pmin(2*dat_b$vaf))
    fit2<-glm(y~x,data=dat_b,family=binomial(),weight=wt2)
    ts[i]<-coef(fit2)["x"]
    ref[i]<-(coef(fit2)["x"]-beta_hat)/sqrt(sandwich::vcovHC(fit2, type = "HC0")["x", "x"])
  }

  ci_beta.b <- quantile(
    ts,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  c("OR" = exp(beta_hat),
    "95% CI" = exp(ci_beta),
    "boot 95%CI" = exp(ci_beta.b),
    "boot p" = mean(abs(T_obs) <= abs(ref)))
}

#TCGA kidney example:
#download files from https://www.cbioportal.org/study/summary?id=kirc_tcga_pub
mut<-read.table("kidney rcc tcga pancancer 2018/data_mutations.txt",sep="\t",header=T,quote="")

mut<-mut[!(mut$Variant_Classification %in% c("Intron","Silent")),]
mut$sample<-substr(mut$Tumor_Sample_Barcode,1,12)
mut$vaf<-mut$t_alt_count/(mut$t_alt_count+mut$t_ref_count)

tcga<-read.table("kidney rcc tcga pancancer 2018/data_clinical_patient.txt",sep="\t",header=T,quote="")
#tcga<-read.table("cbioportal data/kirc_tcga_pub/data_clinical_patient.txt",sep="\t",header=T,quote="")
tcga<-tcga[tcga$PATIENT_ID %in% mut$sample,]

tcga$stage2<-"1,2"
tcga$stage2[tcga$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE III",  "STAGE IV")]<-"3,4"
tcga$stage22<-as.numeric(as.factor(tcga$stage2))-1



mutg<-mut[mut$Hugo_Symbol=="BAP1",]
tcga$BAP1.vaf<-mutg$vaf[match(tcga$PATIENT_ID,mutg$sample)]
tcga$BAP1<-as.numeric(tcga$PATIENT_ID %in% mutg$sample)

mutg<-mut[mut$Hugo_Symbol=="SETD2",]
tcga$SETD2.vaf<-mutg$vaf[match(tcga$PATIENT_ID,mutg$sample)]
tcga$SETD2<-as.numeric(tcga$PATIENT_ID %in% mutg$sample)

mutg<-mut[mut$Hugo_Symbol=="KDM5C",]
tcga$KDM5C.vaf<-mutg$vaf[match(tcga$PATIENT_ID,mutg$sample)]
tcga$KDM5C<-as.numeric(tcga$PATIENT_ID %in% mutg$sample)


mutg<-mut[mut$Hugo_Symbol=="PBRM1",]
tcga$PBRM1.vaf<-mutg$vaf[match(tcga$PATIENT_ID,mutg$sample)]
tcga$PBRM1<-as.numeric(tcga$PATIENT_ID %in% mutg$sample)



#adjusted p-values

wt<-estimate.weights.geom(tcga$stage22,tcga$PBRM1,pmin(1,2*tcga$PBRM1.vaf))
summary(glm(stage22~PBRM1,data=tcga, family=binomial(),weight=wt))
#
boot.p(tcga$stage22,tcga$PBRM1,tcga$PBRM1.vaf,1000)

wt<-estimate.weights.geom(tcga$stage22,tcga$BAP1,pmin(1,2*tcga$BAP1.vaf))
summary(glm(stage22~BAP1,data=tcga, family=binomial(),weight=wt))
#
boot.p(tcga$stage22,tcga$BAP1,tcga$BAP1.vaf,100000)

wt<-estimate.weights.geom(tcga$stage22,tcga$SETD2,pmin(2*tcga$SETD2.vaf))
summary(glm(stage22~SETD2,data=tcga, family=binomial(),weight=wt))
#
boot.p(tcga$stage22,tcga$SETD2,tcga$SETD2.vaf,10000)



#figure for kidney data

library(ggplot2)
library(patchwork)

make_vaf_plot <- function(data, vaf_var, gene_label) {
  
  pval <- wilcox.test(
    data[[vaf_var]] ~ data$stage2,
    na.action = na.omit
  )$p.value
  
  ggplot(data, aes(x = stage2, y = .data[[vaf_var]])) +
    geom_boxplot(
      width = 0.6,
      linewidth = 1.2,
      colour = "black",
      fill = "grey85",
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.10,
      height = 0,
      size = 2.2,
      alpha = 0.85,
      colour = "black"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = 0.92,
      label = paste0("p = ", format.pval(pval, digits = 2, eps = 0.001)),
      size = 5,
      fontface = "bold",
      colour = "black"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "AJCC stage",
      y = paste0(gene_label, " VAF")
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.border = element_rect(
        colour = "black",
        linewidth = 1.1
      ),
      axis.line = element_line(
        colour = "black",
        linewidth = 0.8
      ),
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.8
      ),
      axis.text = element_text(
        colour = "black",
        size = 12
      ),
      axis.title = element_text(
        colour = "black",
        size = 14,
        face = "bold"
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
}

p1 <- make_vaf_plot(tcga, "SETD2.vaf", "SETD2")
p2 <- make_vaf_plot(tcga, "PBRM1.vaf", "PBRM1")
p3 <- make_vaf_plot(tcga, "BAP1.vaf", "BAP1")

p1 | p2 | p3









##########################################
#circle data
#downloaded from https://gitlab.com/sanjanalab/circle/-/tree/main/DB?ref_type=heads

patient_mutations_meta_filtered <- readRDS("circle-main/circle-main/DB/patient_mutations_meta_filtered.rds")
patient_meta <- readRDS("circle-main/DB/patient_meta.rds")
patient_meta$resp2<-as.character(patient_meta$Response2)
patient_meta$resp2[patient_meta$resp2 %in% c("CR","PR")]<-"CR/PR"
patient_meta$resp2<-as.factor(patient_meta$resp2)
a<-mcols(patient_mutations_meta_filtered)


a$VAF<-pmin(1,a$VAF)
a$ID<-as.character(a$ID)
cirmat<-cirmat2<-matrix(0, ncol=6,nrow=length(unique(a$ID)))
rownames(cirmat)<-rownames(cirmat2)<-patient_meta$ID
colnames(cirmat)<-colnames(cirmat2)<-c("BCLAF1", "TP53", "KRAS",  "BRAF","NRAS","PPP6C")
for (i in colnames(cirmat))
{cirmat[unique(a$ID[a$gene_name==i]),i]<-1
cirmat2[,i]<-a$VAF[a$gene_name==i][match(rownames(cirmat),a$ID[a$gene_name==i])]

}


#p-value from paper
for (g in c("BCLAF1", "TP53", "KRAS",  "BRAF"))
{
  
  print(g)
  
  wt<-estimate.weights.geom(as.numeric(as.factor(patient_meta$resp2))-1,cirmat[,g],pmin(1,2*cirmat2[,g]))
  print(summary(glm(patient_meta$resp2~cirmat[,g], family=binomial(),weight=wt)))
  print(exp(summary(glm(patient_meta$resp2~cirmat[,g], family=binomial(),weight=wt))$coef[2,1]))
  print(boot.p( as.numeric(as.factor(patient_meta$resp2))-1,as.numeric(cirmat[,g]),as.numeric(cirmat2[,g]),1000))
}



#figure


library(ggplot2)
library(patchwork)

make_circle_vaf_plot <- function(vaf, response, gene_label) {
  
  d <- data.frame(
    g    = vaf,
    resp = factor(response)
  )
  
  pval <- wilcox.test(g ~ resp, data = d, na.action = na.omit)$p.value
  
  ggplot(d, aes(x = resp, y = g)) +
    geom_boxplot(
      width = 0.60,
      linewidth = 1.2,
      colour = "black",
      fill = "grey85",
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.10,
      height = 0,
      size = 2.2,
      alpha = 0.85,
      colour = "black"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = 0.92,
      label = paste0(
        "p = ",
        format.pval(pval, digits = 2, eps = 0.001)
      ),
      size = 5,
      fontface = "bold",
      colour = "black"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "IO response",
      y = paste0(gene_label, " VAF")
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.border = element_rect(
        colour = "black",
        linewidth = 1.1
      ),
      axis.line = element_line(
        colour = "black",
        linewidth = 0.8
      ),
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.8
      ),
      axis.text = element_text(
        colour = "black",
        size = 12
      ),
      axis.title = element_text(
        colour = "black",
        size = 14,
        face = "bold"
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
}

p1 <- make_circle_vaf_plot(
  cirmat2[, "BCLAF1"],
  patient_meta$resp2,
  "BCLAF1"
)

p2 <- make_circle_vaf_plot(
  cirmat2[, "TP53"],
  patient_meta$resp2,
  "TP53"
)

p3 <- make_circle_vaf_plot(
  cirmat2[, "KRAS"],
  patient_meta$resp2,
  "KRAS"
)

p4 <- make_circle_vaf_plot(
  cirmat2[, "BRAF"],
  patient_meta$resp2,
  "BRAF"
)

(p1 | p2) / (p3 | p4)
