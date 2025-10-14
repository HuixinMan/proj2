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




#Establish an NSEIR model with social structure
nseir <- function(beta, h, alink,              
                  alpha = c(.1, .01, .01),     
                  delta = .2, gamma = .4,       
                  nc = 15, nt = 100, pinf = .005){
# Obtain the total population n
  n <- length(beta)                            
# Establish the individual state vector x and initialize it
x <- rep(0, n)            
# Select the initial infected person randomly (with "2" state)
x[sample.int(n, max(1, round(n * pinf)))] <- 2
  
  # Build the daily counters and initialize it: S(Susceptibility), E(Exposure), I(Infection), R(Recovery)
  S <- E <- I <- R <- rep(0, nt) 
  # Total number of susceptible individuals on Day 1
  S[1] <- sum(x == 0)  
  # Total number of infected people on Day 1
  I[1] <- sum(x == 2)                          
  # Pre-build a family group list (raising efficiency)
  hh <- split(seq_len(n), h)
  # Calculate the average value of beta
  bbar <- mean(beta)   
  # Calculate the normalization constant of random propagation
  c_r <- alpha[3] * nc / (bbar^2 * (n - 1))   
  
  # Main cycle: From day 2 to day nt
  for (t in 2:nt) {                            
    # Find the index of all current infected individuals (2 status)
    idxI <- which(x == 2)   
    # If there are infected individuals
    if (length(idxI)) 
      # Convert the infected individuals to recovered patients with probability delta (3 status)
      x[idxI[runif(length(idxI)) < delta]] <- 3  
    # Find the index of all current exposed individuals (1 status)
    idxE <- which(x == 1)    
    #If there are exposed individuals
    if (length(idxE))         
      # Convert the exposed individuals to infected individuals using probability gamma (2 status)
      x[idxE[runif(length(idxE)) < gamma]] <- 2  
    # Update the index of Infected individuals
    idxI <- which(x == 2)      
    # Find the index of all current susceptible individuals (status =0)
    idxS <- which(x == 0)                    
    # If there are both infected people and susceptible individuals
    if (length(idxI) && length(idxS)) {  
      # Initialize the tag vector(Which susceptible Individuals Will be Infected (FALSE= No)) 
      toE <- rep(FALSE, n)                    
      
# Family spread  
      #If family communication starts
      if (alpha[1] > 0) {    
        # Traverse each infected person i
        for (i in idxI) {    
          # Get the list of members of the i's family
          mem <- hh[[h[i]]]   
          # Exclude i-self
          mem <- mem[mem != i]           
          # Screening: Those who are susceptible and have not been marked by other transmission routes
          mem <- mem[x[mem] == 0 & !toE[mem]]  
          # If there are eligible family members
          if (length(mem))  
          # Mark it as being infected with the probability alpha[1]
            toE[mem[runif(length(mem)) < alpha[1]]] <- TRUE  
        }
      }
      
      # Social network spread
      # If social network starts
      if (alpha[2] > 0) {  
        # Traverse each infected person i
        for (i in idxI) {       
          # Get i's list of social network friends
          nbr <- alink[[i]]   
          # Screening: Susceptible individuals & Not marked by other transmission routes
          nbr <- nbr[x[nbr] == 0 & !toE[nbr]] 
          # If have qualified network friends
          if (length(nbr))       
            # Mark it as being infected with the probability alpha[2]
            toE[nbr[runif(length(nbr)) < alpha[2]]] <- TRUE  
        }
      }
      
      # Random spread
      #If random spread starts
      if (alpha[3] > 0) {
        # Calculate the total social activity of all infected individuals
        BI <- sum(beta[idxI])      
        # If there are infected people
        if (BI > 0) {       
          # Identify susceptible individuals who have not been marked by other pathways
          Sfree <- idxS[!toE[idxS]]   
          # Calculate the infection risk level of each susceptible person j
          lam <- c_r * beta[Sfree] * BI    
          # Convert the risk level to the daily infection probability
          p <- 1 - exp(-lam)
          # Randomly infect some susceptibles based on probability
          x[Sfree [runif(length(Sfree)) < p]] <- 1  
        }
      }
      
      # Convert individuals who have been marked and are still susceptible to exposure
      x[toE & x == 0] <- 1                    
    }
    
    #Record the number of people in each status on that day
    # Total susceptibles on day t
    S[t] <- sum(x == 0)
    # Total exposed on day t
    E[t] <- sum(x == 1)   
    # Total infected on day t
    I[t] <- sum(x == 2) 
    # Total recovered on day t
    R[t] <- sum(x == 3)                        
  }
  
  # return result
  list(S = S, E = E, I = I, R = R, t = 1:nt)
}

# Create a combined plot of all four state variables
# Uses matplot to display S, E, I, R curves simultaneously
# main: plot title
plot_seir <- function(epi, main = "") {
  plot(epi$S,ylim=c(0,max(epi$S)),xlab="day",ylab="N",main=main)
  points(epi$E,col=4);points(epi$I,col=2);points(epi$R,col=3);
  
# Add legend to identify each curve
legend("topright", c("S","E","I","R"),
       lty = 1, col = c(1,4,2,3), bty = "n", cex = .8,pch=1)
}

# Generate random beta values from U(0,1) distribution
beta0 <- runif(n)

# Create social network using random beta values
alink <- get.net(beta0, h, 15)

# Scenario 1: Full model with default parameters
# Random beta + household structure + social network + random mixing
epi1 <- nseir(beta0, h, alink, alpha = c(.1, .01, .01))

# Scenario 2: Remove household and network structure, keep only random mixing
# αh = αc = 0, αr = 0.04 
epi2 <- nseir(beta0, h, alink, alpha = c(0, 0, .04))

# Create constant beta vector (all elements equal to mean of random beta)
beta_bar <- rep(mean(beta0), n)

# Generate new social network using constant beta values
alink_bar <- get.net(beta_bar, h, 15)

# Scenario 3: Full model with constant beta
# Constant beta + household structure + social network + random mixing
epi3 <- nseir(beta_bar, h, alink_bar, alpha = c(.1, .01, .01))

# Scenario 4: Constant beta with only random mixing
# Constant beta + αh = αc = 0, αr = 0.04
epi4 <- nseir(beta_bar, h, alink_bar, alpha = c(0, 0, .04))


# Set up 2x2 plot layout with adjusted margins
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Plot all four scenarios for comparison
plot_seir(epi1, "Full model\nRandom beta + Social structure")
plot_seir(epi2, "Random mixing only\nαh=αc=0, αr=0.04")
plot_seir(epi3, "Full structure\nConstant beta")
plot_seir(epi4, "Random mixing + Constant beta\nαh=αc=0, αr=0.04")

# Results analysis and commentary
# Function to calculate final infected proportion (R + I at end)
final_infected <- function(epi) {
  last <- length(epi$R)
  (epi$R[last] + epi$I[last]) / n
}

# Function to calculate peak infected proportion (maximum I during simulation)
peak_infected <- function(epi) {
  max(epi$I) / n
}

# Calculate metrics for all four scenarios
inf1 <- final_infected(epi1); peak1 <- peak_infected(epi1)
inf2 <- final_infected(epi2); peak2 <- peak_infected(epi2)
inf3 <- final_infected(epi3); peak3 <- peak_infected(epi3)
inf4 <- final_infected(epi4); peak4 <- peak_infected(epi4)

cat("Final infected proportion and peak infected proportion:\n")
cat(sprintf("Scenario 1 (Full model): Final %.1f%%, Peak %.1f%%\n", inf1 * 100, peak1 * 100))
cat(sprintf("Scenario 2 (Random mixing only): Final %.1f%%, Peak %.1f%%\n", inf2 * 100, peak2 * 100))
cat(sprintf("Scenario 3 (Full structure + constant beta): Final %.1f%%, Peak %.1f%%\n", inf3 * 100, peak3 * 100))
cat(sprintf("Scenario 4 (Random mixing + constant beta): Final %.1f%%, Peak %.1f%%\n", inf4 * 100, peak4 * 100))


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
