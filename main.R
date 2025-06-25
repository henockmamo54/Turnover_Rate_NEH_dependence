


library(dplyr)
library(minpack.lm)
library(nlme)
library(MASS)
library(ggplot2)
library(dplyr)

source("F:/professoe_turnover/New_reading/allCSVfiles/Publication/Git/code/Helper.R")

folder_path<- "SampleData" # <- please use your dataset


############################## log Files##########
warning_ria_log_file <- "F:/professoe_turnover/New_reading/allCSVfiles/Publication/Git/code/log/warning_ria_error.txt"
error_log_file <- "F:/professoe_turnover/New_reading/allCSVfiles/Publication/Git/code/log/warning_curve_fit_errors.txt"
####################################################
############ Helper Functions ##########



quant_files <- list.files(folder_path, pattern = "\\.Quant\\.csv$", full.names = TRUE)

# Initialize list to store merged results
all_merged <- list()

for (quant_file in quant_files) {
  tryCatch({
    # Generate corresponding RateConst filename
    rate_file <- sub("\\.Quant\\.csv$", ".RateConst.csv", quant_file)
    
    # Check if both files exist
    if (!file.exists(rate_file)) {
      message(paste("Skipping", basename(quant_file), "- corresponding RateConst file not found"))
      next
    }
    
    # Read files with appropriate skip values
    data_quant <- read_csv_robust2(quant_file, skip = 3)
    data_rate <- read_csv_robust2(rate_file, skip = 0)
    
    # Skip if either file failed to read
    if (nrow(data_quant) == 0 || nrow(data_rate) == 0) {
      message(paste("Skipping", basename(quant_file), "- empty data detected"))
      next
    }
    
    # Standardize column names
    names(data_rate)[names(data_rate) == 'Peptides'] <- 'Peptide'
    
    # Filter rate data
    # data_rate_filtered <- data_rate %>%
    #   filter(as.numeric(Rsquared) >= rsquared_threshold)
    
    # Merge data frames
    merged <- inner_join(data_quant, data_rate, by = c("Peptide", "Charge"))
    
    # Store the merged result
    all_merged[[basename(quant_file)]] <- merged
    
  }, error = function(e) {
    message(paste("Error processing", quant_file, ":", e$message))
  })
}

# # Combine all merged data frames into one
final_merged_liverpool_ch60 <- bind_rows(all_merged, .id = "SourceFile")

#=====================================================================

###################
##### readin
merged<-final_merged_liverpool_ch60
rsquared_threshold <- 0.80
merged<-merged %>% filter(Rsquared>=rsquared_threshold)
####################

###########################################################################

results <- c()
results_k_ratio<-c()
results_1 <- c()
results_2 <- c()
results_2_2 <- c()
all_neh <- c()
all_rates <- c()
all_index <- c()
#merged<-mergeded
################# important ################# 

#############################################
# Create empty data frame to store pxt values and times
pxt_data <- data.frame(SourceFile = character(),
                       #Peptide = character(),
                       peptide_len= numeric(),
                       Time = numeric(),
                       Pxt = numeric(),
                       #Pxt_new = numeric(),
                       #pxt_perivious= numeric(),
                       NEH= numeric(),
                       #NH= numeric(),
                       peptide_name=character(),
                       peptide_charge=numeric(),
                       k_ratio=numeric(),
                       peptide_index=numeric(),
                       m0= numeric(),
                       #B_asymp= numeric(),
                       #I0_asymp= numeric(),
                       RootMeanRSS= numeric(),
                       NDP= numeric(),
                       RSS_csv= numeric(),
                       Rsquared_csv= numeric(),
                       stringsAsFactors = FALSE)

# exp_time <- c(0, 1, 2, 3, 6, 7, 9, 13, 16, 21, 24, 31)
exp_time <- c(0, 1, 2, 3, 4, 5, 6, 14, 21) # for used datset
# pw <- 0.046
pw <- 0.0399  # for utmb dataset
ph <- 1.5574E-4 

for (index in seq_len(nrow(merged))) {
  #if (merged$Peptide[index]=="LFAEAVQK"){ 
  tryCatch({
    # Correct column names
    neh <- merged$NEH[index]
    
   
    m0_v<-merged$M0[index] / 100
    
    
    peptide_SourceFile<-merged$SourceFile[index]
    peptide_name_v <- merged$Peptide[index]  # Get peptide name
    #m0<-compute_M0(peptide_name)
    peptide_len_v <-nchar(peptide_name_v)
    #NH_v<-ElementalComposition(toupper(peptide_name_v))
    peptide_charge_v <- merged$Charge[index]
    k_ratio_v<-merged$RateConstants[index]
    RootMeanRSS_v<-merged$RootMeanRSS[index]
    NDP_v<-merged$NDP[index]
    RSS_csv_v<- (RootMeanRSS_v)^2*NDP_v
    Rsquared<-merged$Rsquared[index]
    #if( peptide_name=="AADTIGYPVMIR"){
    #if( peptide_name=="AMLSTGFK"){
    #if( peptide_name=="VVAVDcGIK"){
    peptide_times <- c()
    peptide_pxts <- c()
    #peptide_pxts_new <- c()
    peptide_pxts_perivious<- c()
    # peptide_pxts_eq2 <- c()
    for (i in seq_along(exp_time)) {
      time <- exp_time[i] # Get time of experiment
      time_index <- i  #Correct time_index calculation
     
      
      # pxt <- get_pxt(time_index-1, index, neh,merged )
      pxt <- get_pxt(2*(time_index-1), index, neh,merged ) # for utmb dataset
      
     
      if (!is.na(pxt) && is.finite(pxt)) {
        
        peptide_times <- c(peptide_times, time)
        peptide_pxts <- c(peptide_pxts, pxt)
        
        #pxt_new <- pxt*(NH_v/neh)
        #peptide_pxts_new <- c(peptide_pxts_new, pxt_new)
        #peptide_pxts_perivious <- c(peptide_pxts_perivious,pxt_perivious)
        #peptide_pxts_eq2 <- c(peptide_pxts_eq2, pxt_eq2)
        
        # Store the pxt value and time in the data frame
        pxt_data <- rbind(pxt_data, data.frame(SourceFile=peptide_SourceFile,
                                               #Peptide = peptide_name,
                                               Time = time,
                                               Pxt = pxt,
                                               #Pxt_new = pxt_new,
                                               #pxt_perivious=pxt_perivious,
                                               NEH=neh,
                                               #NH=NH_v,
                                               peptide_name= peptide_name_v,
                                               peptide_len=peptide_len_v,
                                               peptide_charge=peptide_charge_v,
                                               k_ratio=k_ratio_v,
                                               peptide_index=index,
                                               m0=m0_v,
                                               #B_asymp= B_asymp,
                                               #I0_asymp=I0_asymp,
                                               RootMeanRSS= RootMeanRSS_v,
                                               NDP= NDP_v,
                                               RSS_csv= RSS_csv_v,
                                               Rsquared_csv=Rsquared
        ))
        
        
      }
      
      
    }
    
  }, error = function(e) {
    message(paste("Error at index", index, "- skipping iteration.",e))
  })
}

########################################################################################################################
########################################################################################################################
#### Compute adjusted Pₓ(t) ####

#pxt_data_adj<-adjust_pxt(pxt_data)
pxt_data_adj<-adjust_pxt_all(pxt_data)
#########################################################################

all_RSS <- c()
results_adj <- c()
results <- c()
all_RSS<- c()
# results_1_fm <- c()
results_k_ratio<-c()
all_index <- c()


results_adj_df<- data.frame(SourceFile= character(),peptide_name =  character(),
                            k=numeric(),RSS=numeric(),neh=numeric(),peptide_index=numeric(),
                            
                            rate=numeric(),time_cou=numeric(),
                            Rsquared_k=numeric(),
                            
                            peptide_charge=numeric(),
                            RootMeanRSS= numeric(),
                            NDP= numeric(),
                            RSS_csv= numeric(),
                            Rsquared_csv= numeric(),
                            
                            stringsAsFactors = FALSE)



##########################################################
# Initialize tracking array
T_done <- vector("logical", length = nrow(merged))
for (index in seq_len(nrow(merged))) {
  # Skip already processed indices#
  if (T_done[index]) next
  #peptide_pxts_adjust<-adjust_pxt(pxt_data[pxt_data$peptide_index==index,])
  peptide_pxts_adjust<-pxt_data_adj[pxt_data_adj$peptide_index==index,]$Pxt_adj
  
  peptide_pxts<-pxt_data_adj[pxt_data_adj$peptide_index==index,]$Pxt
  
  peptide_times2<-pxt_data_adj[pxt_data_adj$peptide_index==index,]$Time
  lt<-length(unique( peptide_times2))
  #B_asymp<-unique(pxt_data_adj[pxt_data_adj$peptide_index==index,]$B_asymp[1])
  neh_v<-unique(pxt_data_adj[pxt_data_adj$peptide_index==index,]$NEH)
  
  
  
  
  # B_values_adjust <- ((1 - (peptide_pxts_adjust-ph) / (neh * (1 - ph)))^neh)
  
  # B_values <- ((1 - (peptide_pxts-ph) / (neh * (1 - ph)))^neh)
  # B_values<- (1 - peptide_pxts-ph)^neh
  ############ Eq.4
  # B_values<- (1 - (peptide_pxts/(1-ph)))^neh
  # B_values_adjust <- (1 - (peptide_pxts_adjust/(1-ph)))^neh
  ############ Eq.9
  # B_values9<- (1 - (peptide_pxts/(neh*(1-ph))))^neh
  B_values_adjust9 <- (1 - (peptide_pxts_adjust/(neh_v*(1-ph))))^neh_v
  #B_values_adjust9_pxt <- (1 - (peptide_pxts/(neh*(1-ph))))^neh
  
  #B_values_adjust <- ((1 - (peptide_pxts_adjust) / (neh * (1 )))^neh)
  #peptide_pxts_perivious<-pxt_data_adj[pxt_data_adj$peptide_index==index,]$pxt_perivious
  m0<-unique(pxt_data_adj[pxt_data_adj$peptide_index==index,]$m0[1])
  peptide_SourceFile<-pxt_data_adj[pxt_data_adj$peptide_index==index,]$SourceFile[1]
  peptide_name_v<-pxt_data_adj[pxt_data_adj$peptide_index==index,]$peptide_name[1]
  
  peptide_charge_v<-unique(pxt_data_adj[pxt_data_adj$peptide_index==index,]$peptide_charge)
  RootMeanRSS_v<-unique(pxt_data_adj[pxt_data_adj$peptide_index==index,]$RootMeanRSS)
  NDP_v<-unique(pxt_data_adj[pxt_data_adj$peptide_index==index,]$NDP)
  RSS_csv_v<-unique(pxt_data_adj[pxt_data_adj$peptide_index==index,]$RSS_csv)
  Rsquared_csv_v<-unique(pxt_data_adj[pxt_data_adj$peptide_index==index,]$Rsquared_csv)
  
  tryCatch({
    
    
    # cof<-(1-pw)^neh
    # cof<-((1-pw)/(1-ph))^neh
    cof9<-(1-(pw/(neh_v*(1-ph))))^neh_v
    
    maxrate<-max(pxt_data_adj$k_ratio)
    # fit_result <-  d2ome_curve_fit_adj_pxt(cof9, peptide_times2, B_values9) 
    fit_result <-  d2ome_curve_fit_adj_pxt(cof9, peptide_times2,B_values_adjust9,maxrate)
    
    k_v <-fit_result$k
    RSS_v<-fit_result$RSS
    rsq_k<-fit_result$rsq_k
    if (is.na(k_v)) {
      message(paste("Failed to estimate k for peptide", peptide_name_v, "at index", index))
      msg <- paste("Failed to estimate k for peptide", peptide_name_v, "at index", index)
      message(msg)
      write(msg, file = error_log_file, append = TRUE)
    }
    
    else {
      rate_v<-unique(pxt_data_adj[pxt_data_adj$peptide_index==index,]$k_ratio)
      results <- c(results,as.numeric(k_v))
      all_RSS <- c(all_RSS,as.numeric(RSS_v))
      results_k_ratio <- c(results_k_ratio ,rate_v)
      all_index <- c(all_index, index)
      all_neh<-c(all_neh,neh_v)
      # results_2 <- c( results_2, k_2)
      # results_2_2 <- c( results_2_2, k_2_2)
      results_adj_df <- rbind(results_adj_df, data.frame(SourceFile=peptide_SourceFile,
                                                         peptide_name = peptide_name_v,
                                                         k = k_v,RSS=as.numeric(RSS_v),neh=neh_v,
                                                         peptide_index=index, 
                                                         # peptide_charge=peptide_charge,
                                                         rate=rate_v,time_cou=lt,
                                                         Rsquared_k=rsq_k,
                                                         peptide_charge=peptide_charge_v,
                                                         RootMeanRSS=RootMeanRSS_v,
                                                         NDP=NDP_v,
                                                         RSS_csv=RSS_csv_v,
                                                         Rsquared_csv=Rsquared_csv_v
      ))
      
    }
  }, error = function(e) {
    message(paste("Error in curve fitting for peptide", peptide_name_v, "at index", index, ":", e))
  })
  #############
  
  T_done[index]<-TRUE 
}
