# Project 2: SEIR model with social structure
# Group member: S2809410 Huixin Man, S2793485 Ruoxi Wang, S2829716 Leying Wang
# The contribution ratio of Huixin, Ruoxi and Leying are 32%, 34% and 34%.
# Github Link:

# Randomly generate family IDs (h) to simulate household structure
set.seed(123)
hmax <- 5
n <- 1000
# Generate family sizes ~ U(1, hmax) ensuring total population ≤ n
sizes <- sample(1:hmax, n, replace = TRUE)          
sizes <- sizes[cumsum(sizes) <= n]                  
sizes[length(sizes)] <- sizes[length(sizes)] + (n - sum(sizes))  
# Assign a family number to each individual based on family sizes
h <- rep(seq_along(sizes), sizes)      
# Disorder the sequence of family numbers, so that individuals are randomly distributed in the data.
h <- sample(h)                                      
length(h);


# get.nut function: return a list, the ith element of which is a vector of the indices of the regular (non-household) contacts of person i.
get.net <- function(beta, h, nc = 15) {
  n <- length(beta)
  bbar <- mean(beta)
  const <- nc / (bbar^2 * (n - 1))
  links <- vector("list", n)
  for (i in 1:(n - 1)) {
    # Potential contacts for person i, exclude individuals from the same family
    cand_j <- (i + 1):n
    cand_j <- cand_j[h[cand_j] != h[i]]
    if (length(cand_j) == 0) next
    # Link probabilities proportional to contact propensities
    p <- const * beta[i] * beta[cand_j]
    # Cap probabilities at 1
    p[p > 1] <- 1
    # Randomly form links based on probability p
    chosen <- cand_j[runif(length(cand_j)) < p]
    ## Add symmetric connections
    if (length(chosen)) {
      links[[i]] <- c(links[[i]], chosen)
      for (j in chosen) links[[j]] <- c(links[[j]], i)
    }
  }
  #Remove duplicates and sort connections
  for (i in seq_len(n)) if (length(links[[i]]) > 1) links[[i]] <- sort(unique(links[[i]]))
  links
}
# Establish an NSEIR model with social structure
nseir <- function(beta, h, alink,              
                  alpha = c(.1, .01, .01),     
                  delta = .2, gamma = .4,       
                  nc = 15, nt = 100, pinf = .005){
#Obtain the total population n
  n <- length(beta)                            
#Establish the individual state vector x and initialize it
x <- rep(0, n)            
#Select the initial infected person randomly (with "2" state)
x[sample.int(n, max(1, round(n * pinf)))] <- 2
  
  #Build the daily counters and initialize it: S(Susceptibility), E(Exposure), I(Infection), R(Recovery)
  S <- E <- I <- R <- rep(0, nt) 
  #Total number of susceptible individuals on Day 1
  S[1] <- sum(x == 0)  
  #Total number of infected people on Day 1
  I[1] <- sum(x == 2)                          
  #Pre-build a family group list (raising efficiency)
  hh <- split(seq_len(n), h)
  #Calculate the average value of beta
  bbar <- mean(beta)   
  #Calculate the normalization constant of random propagation
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
    
  # Record the number of people in each status on that day
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
plot_seir(epi1, "Full model\nDefault parameters")
plot_seir(epi2, "Random mixing\nαh=αc=0, αr=0.04")
plot_seir(epi3, "Full model\nConstant β")
plot_seir(epi4, "Random mixing + Constant β\nαh=αc=0, αr=0.04")

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

# Output the result
cat("Final infected proportion and peak infected proportion:\n")
cat(sprintf("Scenario 1 (Full model): Final %.1f%%, Peak %.1f%%\n", inf1 * 100, peak1 * 100))
cat(sprintf("Scenario 2 (Random mixing only): Final %.1f%%, Peak %.1f%%\n", inf2 * 100, peak2 * 100))
cat(sprintf("Scenario 3 (Full model + constant β): Final %.1f%%, Peak %.1f%%\n", inf3 * 100, peak3 * 100))
cat(sprintf("Scenario 4 (Random mixing + constant β): Final %.1f%%, Peak %.1f%%\n", inf4 * 100, peak4 * 100))

# Comparing Scenarios 1 and 2, the peaks of E and I occur later in Scenario 1, 
# and the peak proportion of infectious individuals is smaller than in Scenario 2. 
# This indicates that the presence of fixed household and network structure slightly slows down the epidemic spread.

# Comparing Scenarios 1 and 3, the peak and final proportion of infectious individuals in Scenario 1 
# are both lower than in Scenario 3, suggesting that variability in β weakens the epidemic peak 
# and substantially reduces the overall epidemic size.

# Comparing Scenarios 2 and 4, the final proportion infected in Scenario 2 is lower than in Scenario 4, 
# consistent with the finding that β variability reduces the total epidemic size, 
# while the peak infection proportions are similar in both scenarios.

# Above all, fixed social structure leads to a slower and flatter epidemic curve, 
# whereas individual variability in β results in a weaker and more limited epidemic spread.



