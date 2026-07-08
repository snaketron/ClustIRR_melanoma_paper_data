Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/home/simo/miniconda3/bin/", sep=.Platform$path.sep))
library(ClustIRR)
dir.create("/mnt/nfs/simo/AIRR/TCR_PT_SLN/results/")
dir.create("/mnt/nfs/simo/AIRR/TCR_PT_SLN/results/clustirr_90")

# 1. prepare data
t <- get(load(file = "/mnt/nfs/simo/AIRR/TCR_PT_SLN/results/TCRb.RData"))
d <- t[, c("cdr3_amino_acid", "sample_name", "templates")]
colnames(d) <- c("CDR3b", "sample", "clone_size")
m <- t
rm(t)
gc();gc();gc()


# 2. perform ClustIRR analysis
c <- clustirr(s = d, meta = m, control = list(blast_cores = 40, blast_gmi = 0.9))
save(c, file = "/mnt/nfs/simo/AIRR/TCR_PT_SLN/results/clustirr_90/c.RData")
cat("DONE!\n")


# 3. detect_communities
dc <- detect_communities(graph = c$graph, algorithm = "leiden", metric = "average", 
                         resolution = 1, iterations = 1000, chains = "CDR3b")
save(dc, file = "/mnt/nfs/simo/AIRR/TCR_PT_SLN/results/clustirr_90/dc.RData", compress = TRUE)
cat("DONE!\n")
graph <- dc$graph
save(graph, file = "/mnt/nfs/simo/AIRR/TCR_PT_SLN/results/clustirr_90/graph.RData", compress = TRUE)
cat("DONE!\n")


# 4. get concise summary from dc 
dc_light <- list(com = dc$community_occupancy_matrix, ns = dc$node_summary, cs = dc$community_summary)
save(dc_light, file = "/mnt/nfs/simo/AIRR/TCR_PT_SLN/results/clustirr_90/dc_light.RData", compress = TRUE)


# 5. perform differential community occupancy for each patient
dir.create("/mnt/nfs/simo/AIRR/TCR_PT_SLN/results/clustirr_90/dco/", recursive = T)
dc_light <- get(load(file = "results/clustirr_90//dc_light.RData"))
com <- dc_light$com
pts <- unique(gsub(pattern = "PT|SLN", replacement = '', x = unique(colnames(com)), ignore.case = FALSE))

lapply(X = pts, com, FUN = function(x, com) {
  cs <- com[, c(paste0(x, "PT"), paste0(x, "SLN"))]
  rownames(cs) <- 1:nrow(cs)
  j <- which(cs[,1]==0 & cs[,2]==0)
  if(length(j)!=0) {
    cs <- cs[-j,]
  }
  
  if(!file.exists(paste0("/mnt/nfs/simo/AIRR/TCR_PT_SLN/results/clustirr_90/dco/ps_", x, ".RData"))) {
    d <- ClustIRR::dco(community_occupancy_matrix = cs,
                       mcmc_control = list(mcmc_chains = 4,
                                           mcmc_cores = 4,
                                           mcmc_warmup = 500,
                                           mcmc_iter = 1500))
    
    save(d, file = paste0("/mnt/nfs/simo/AIRR/TCR_PT_SLN/results/clustirr_90/dco/", x, ".RData"))
    ps <- d$posterior_summary
    save(ps, file = paste0("/mnt/nfs/simo/AIRR/TCR_PT_SLN/results/clustirr_90/dco/ps_", x, ".RData"))
    save(cs, file = paste0("/mnt/nfs/simo/AIRR/TCR_PT_SLN/results/clustirr_90/dco/cs_", x, ".RData"))
  }
})

