setwd("~/Dropbox/2222/stat_499/project/R")

library(glasso)
library(plotrix)

col_name = read.table('labels_aal.txt')  
col_name = as.character(col_name[,2])

Names = list.files(pattern = "ID")
ADHD_data = vector('list', length(Names))
for (i in seq_along(Names))
{
  ADHD_data[[i]] = read.table(Names[i], header = T, col.names = col_name)
}


ID_ADHD_data = gsub('.txt', '', Names)
ID_ADHD_data = gsub('ID', '', ID_ADHD_data)
ID_ADHD_data = gsub('^00', '', ID_ADHD_data)
ID_ADHD_data = as.numeric(ID_ADHD_data)

Obs_info = read.table('NYU_phenotypic_LK.txt', header = T)

sum(ID_ADHD_data == sort(Obs_info$ID)) ## get 222 which is good

Obs_info = Obs_info[order(Obs_info$ID), ]

set.seed(0808)

m = 10 
Lam = seq(1e-4, 1e1, length.out = m)
output_222 = vector('list', 222)
for (i in 1:1)
{

    output_87 = vector('list', 87)
  for (j in 1:87)
  { 
    output_lambda = vector('list', m)
    matrix = ADHD_data[[i]][j:(j+85), ]
    S2 = var(matrix)
    cold_glasso = glasso(S2, rho = 0.1)
    for (k in seq_along(Lam))
    {
      result_glasso = glasso(S2, rho = Lam[k], 
                             w.init=cold_glasso$w, 
                             wi.init=cold_glasso$wi)
      output_lambda[[k]] = list(w = result_glasso$w,
                                wi = result_glasso$wi,
                                log_lik   = result_glasso$loglik)
    }
    output_87[[j]] = output_lambda
  }
  output_222[[i]] = output_87
}

P = vector('list', 222)
for (i in 1:1)  
{
  PL = vector('list', m)
  for (k in seq_along(Lam))
  {
    P_matrix = matrix(0, nrow = 116, ncol = 116)
    for (j in 1:87)
    {
      P_matrix = P_matrix + (output_222[[i]][[j]][[k]]$wi != 0) 
    }
    P_matrix = P_matrix/87
    PL[[k]] = P_matrix
  }
  P[[i]] = PL
}
save(output_222, P, file = 'out.rda')
color2D.matplot(P[[1]][[3]], c(1,0),c(1,0.3),c(1,0.55),show.legend=TRUE,
                 ylab = '', xlab = '',
                 main = '')
