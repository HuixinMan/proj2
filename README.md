# proj2
## Practical 2: SEIR model with social structure
# 随机生成个体的“家庭编号”向量 h，以模拟家庭聚集结构
set.seed(123)
hmax <- 5
n <- 1000
# 随机生成家庭规模（每个家庭的人数），总人数不超过n
sizes <- sample(1:hmax, n, replace = TRUE)          
sizes <- sizes[cumsum(sizes) <= n]                  
sizes[length(sizes)] <- sizes[length(sizes)] + (n - sum(sizes))  
h <- rep(seq_along(sizes), sizes)                   
h <- sample(h)                                      
length(h);

#生成固定社交网络（排除家庭内链接）

get.net <- function(beta, h, nc = 15) {
  n <- length(beta)
  bbar <- mean(beta)
  const <- nc / (bbar^2 * (n - 1))
  links <- vector("list", n)
  for (i in 1:(n - 1)) {
    cand_j <- (i + 1):n
    cand_j <- cand_j[h[cand_j] != h[i]]
    if (length(cand_j) == 0) next
    p <- const * beta[i] * beta[cand_j]
    p[p > 1] <- 1
    chosen <- cand_j[runif(length(cand_j)) < p]
    if (length(chosen)) {
      links[[i]] <- c(links[[i]], chosen)
      for (j in chosen) links[[j]] <- c(links[[j]], i)
    }
  }
  for (i in seq_len(n)) if (length(links[[i]]) > 1) links[[i]] <- sort(unique(links[[i]]))
  links
}

## -------------------------------------------------
## 3. 带社会结构的 SEIR 模型
## -------------------------------------------------
nseir <- function(beta, h, alink,
                  alpha = c(.1, .01, .01),
                  delta = .2, gamma = .4, nc = 15,
                  nt = 100, pinf = .005) {
  n <- length(beta)
  x <- rep(0,n)
  x[sample.int(n, max(1, round(n * pinf)))] <- 2  
  S <- E <- I <- R <- rep(0,nt)
  S[1] <- sum(x == 0); I[1] <- sum(x == 2)
  hh <- split(seq_len(n), h)
  bbar <- mean(beta)
  c_r <- alpha[3] * nc / (bbar^2 * (n - 1))
  
  for (t in 2:nt) {
    idxI <- which(x == 2)
    if (length(idxI)) x[idxI[runif(length(idxI)) < delta]] <- 3
    idxE <- which(x == 1)
    if (length(idxE)) x[idxE[runif(length(idxE)) < gamma]] <- 2
    idxI <- which(x == 2); idxS <- which(x == 0)
    if (length(idxI) && length(idxS)) {
      toE <- rep(FALSE, n)
      if (alpha[1] > 0) {
        for (i in idxI) {
          mem <- hh[[h[i]]]
          mem <- mem[mem != i]
          mem <- mem[x[mem] == 0 & !toE[mem]]
          if (length(mem))
            toE[mem[runif(length(mem)) < alpha[1]]] <- TRUE
        }
      }
      if (alpha[2] > 0) {
        for (i in idxI) {
          nbr <- alink[[i]]
          nbr <- nbr[x[nbr] == 0 & !toE[nbr]]
          if (length(nbr))
            toE[nbr[runif(length(nbr)) < alpha[2]]] <- TRUE
        }
      }
      if (alpha[3] > 0) {
        BI <- sum(beta[idxI])
        if (BI > 0) {
          Sfree <- idxS[!toE[idxS]]
          lam <- c_r * beta[Sfree] * BI
          p <- 1 - exp(-lam)
          x[Sfree[runif(length(Sfree)) < p]] <- 1
        }
      }
      
      
      x[toE & x == 0] <- 1
    }
    S[t] <- sum(x == 0); E[t] <- sum(x == 1)
    I[t] <- sum(x == 2); R[t] <- sum(x == 3)
  }
  list(S = S, E = E, I = I, R = R, t = 1:nt)
}

## -------------------------------------------------
## 4. 绘图函数（包含 S、E、I、R 四条曲线）
## -------------------------------------------------
plot_seir <- function(epi, main = "") {
  plot(epi$S,ylim=c(0,max(epi$S)),xlab="day",ylab="N",main=main)
  points(epi$E,col=4);points(epi$I,col=2);points(epi$R,col=3);
  legend("topright", c("S","E","I","R"),
         lty = 1, col = c(1,4,2,3), bty = "n", cex = .8,pch=1)
}




## 5. 模拟与比较（四种情景）
## -------------------------------------------------
beta0 <- runif(n)
alink <- get.net(beta0, h, 15)

epi1 <- nseir(beta0, h, alink, alpha = c(.1,.01,.01))
epi2 <- nseir(beta0, h, alink, alpha = c(0,0,.04))
beta_bar <- rep(mean(beta0), n)
alink_bar <- get.net(beta_bar, h, 15)
epi3 <- nseir(beta_bar, h, alink_bar, alpha = c(.1,.01,.01))
epi4 <- nseir(beta_bar, h, alink_bar, alpha = c(0,0,.04))

cat("生成比较图表...\n")
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot_seir(epi1, "Full model\nRandom beta + Social structure")
plot_seir(epi2, "Random mixing only\nαh=αc=0, αr=0.04")
plot_seir(epi3, "Full structure\nConstant beta")
plot_seir(epi4, "Random mixing + Constant beta\nαh=αc=0, αr=0.04")

cat("\n=== 结果分析和评论 ===\n")

# 计算最终感染比例和峰值感染者数量
final_infected <- function(epi) {
  last <- length(epi$R)
  (epi$R[last] + epi$I[last]) / n
}

peak_infected <- function(epi) {
  max(epi$I) / n
}

inf1 <- final_infected(epi1); peak1 <- peak_infected(epi1)
inf2 <- final_infected(epi2); peak2 <- peak_infected(epi2)
inf3 <- final_infected(epi3); peak3 <- peak_infected(epi3)
inf4 <- final_infected(epi4); peak4 <- peak_infected(epi4)

cat("最终感染比例和峰值感染者比例:\n")
cat(sprintf("情景1 (完整模型): 最终%.1f%%, 峰值%.1f%%\n", inf1 * 100, peak1 * 100))
cat(sprintf("情景2 (仅随机混合): 最终%.1f%%, 峰值%.1f%%\n", inf2 * 100, peak2 * 100))
cat(sprintf("情景3 (完整结构+常数beta): 最终%.1f%%, 峰值%.1f%%\n", inf3 * 100, peak3 * 100))
cat(sprintf("情景4 (常数beta+仅随机混合): 最终%.1f%%, 峰值%.1f%%\n", inf4 * 100, peak4 * 100))

cat("\n家庭和网络结构相对于随机混合的明显影响:\n")
# 比较情景1和情景2：社会结构的影响
if (peak1 < peak2) {
  cat("1. 社会结构减缓了传播速度: 完整模型(情景1)的峰值感染者比例(", round(peak1 * 100, 1), 
      "%)低于随机混合(情景2)的", round(peak2 * 100, 1), "%，表明家庭和网络结构延缓了疫情高峰的到来。\n")
} else {
  cat("1. 社会结构加速了局部传播: 在某些参数下，社会结构可能导致更快的局部传播。\n")
}

# 比较情景1和情景3：beta变异性的影响
if (inf1 > inf3) {
  cat("2. beta变异性增加了传播风险: 随机beta(情景1)的最终感染规模(", round(inf1 * 100, 1),
      "%)高于常数beta(情景3)的", round(inf3 * 100, 1), "%，表明社交活跃度的差异促进了超级传播事件。\n")
} else {
  cat("2. beta变异性对最终规模影响有限: 两种情景的最终感染规模相似。\n")
}

# 比较情景2和情景4：随机混合下beta变异性的影响
if (inf2 != inf4) {
  cat("3. 在随机混合中，beta变异性", ifelse(inf2 > inf4, "增加", "减少"), 
      "了最终感染规模，从", round(inf4 * 100, 1), "%变为", round(inf2 * 100, 1), "%。\n")
} else {
  cat("3. 在随机混合中，beta变异性对最终感染规模影响较小。\n")
}
