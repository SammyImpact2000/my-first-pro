##Function 1 
translate_dna <- function(dna_sequence){
# create a codon table  
  codon_table <- setNames(
    unlist(strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIIIIIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")),
    c(outer(c("T", "C", "A", "G"), c("T", "C", "A", "G"), paste0)))
# Ensure sequence length is a multiple of 3
    dna_sequence <- substr(dna_sequence, 1, (nchar(dna_sequence) %/% 3) * 3)
     # Split sequence into codons
    codons <- substring(dna_sequence, seq(1, nchar(dna_sequence), 3), seq(3, nchar(dna_sequence), 3))
   # Translate codons to amino acids
  protein_sequence <- paste(codon_table[codons], collapse = "")
 return(protein_sequence) }

 # Example DNA sequence
dna_seq <- "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
protein_seq <- translate_dna(dna_seq)

# Print the translated protein sequence
cat("Protein Sequence:", protein_seq, "\n")
#  Function 2
  library(dplyr)

simulate_logistic_growth <- function(K = 1000, r = 0.5, time_steps = 100) {
  lag_phase <- sample(5:20, 1); exp_phase <- sample(20:50, 1)
  time <- 1:time_steps; pop <- numeric(time_steps); pop[1] <- 1
  
  for (t in 2:time_steps) {
    pop[t] <- if (t <= lag_phase) pop[t-1] else if (t <= lag_phase + exp_phase) 
      pop[t-1] * exp(r) else pop[t-1] + r * pop[t-1] * (1 - pop[t-1]/K)
    pop[t] <- min(pop[t], K)
  }
  data.frame(Time = time, Population = pop)
}

growth_curves <- bind_rows(lapply(1:100, function(i) mutate(simulate_logistic_growth(), Curve_ID = i)))

#Funtion 3

time_to_80_percent_K <- function(df, K = 1000) {
  df %>% filter(Population >= 0.8 * K) %>% group_by(Curve_ID) %>% summarise(Time = min(Time), .groups = 'drop')
}

time_80_percent <- time_to_80_percent_K(growth_curves)
head(time_80_percent)
#Function 4
# Function to calculate Hamming Distance
hamming_distance <- function(str1, str2) {
  # Make both strings the same length by padding the shorter one
  len1 <- nchar(str1)
  len2 <- nchar(str2)
  
  if (len1 < len2) {
    str1 <- paste0(str1, strrep(" ", len2 - len1))
  } else if (len2 < len1) {
    str2 <- paste0(str2, strrep(" ", len1 - len2))
  }
  
  # Compute Hamming distance
  distance <- sum(substr(str1, 1:nchar(str1), 1:nchar(str1)) != 
                  substr(str2, 1:nchar(str2), 1:nchar(str2)))
  
  return(distance)
}

# Your Slack username and Twitter/X handle
slack_name <- "Jesuranti Samuel"
twitter_handle <- "Jesuranti137555"

# Calculate Hamming Distance
distance <- hamming_distance(slack_name, twitter_handle)
cat("Hamming Distance:", distance, "\n")
