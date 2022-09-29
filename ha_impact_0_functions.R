#..........................................................................................
### ++++++++ ESTIMATION OF THE IMPACT OF HUMANITARIAN ASSISTANCE ON MORTALITY +++++++++ ###
#..........................................................................................

#..........................................................................................
## --------------------------------- FUNCTIONS FOR ANALYSIS ---------------------------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Jan 2022)
                                          # francesco.checchi@lshtm.ac.uk 



#.........................................................................................
### Function that calculates the proportion of a given month's days that is covered by a survey's recall period
#.........................................................................................

f_calc_days <- function(f_surveys_cov, f_df) {
  # select survey
  s <- subset(f_df, survey_id == f_surveys_cov["survey_id"])
  tm_now <- as.integer(f_surveys_cov["tm"])
  c1 <- as.integer(s["tm_recall_start"]) - tm_now
  c2 <- as.integer(s["tm_recall_end"])- tm_now
  x1 <- 0
  
  # calculate proportion of the month's days that are covered by the survey's recall period
  if (c1 > 0) { x1 <- 0.0 }
  if (c1 == 0) { x1 <- (as.integer(s["days_in_month_start"]) - as.integer(s["day_start"]) ) / 
    as.integer(s["days_in_month_start"]) }
  if (c1 < 0 & c2 > 0) { x1 <- 1.0 }
  if (c2 == 0) { x1 <- as.integer(s["day_end"]) / as.integer(s["days_in_month_end"]) }
  if (c2 < 0) { x1 <- 0.0 }  
    
  return(x1)
  }
  

#.........................................................................................
### Function that calculates and formats median, range and number of observations for any quantity, overall and by year
#.........................................................................................

f_calc_svy <- function(f_quantity, f_df, f_digits, f_years) {
    
    # output of n elements, where n = 1 (total) + number of years
    x1 <- rep(NA, times=length(f_years + 1))
  
    # overall (all years)
    x1[1] <- paste( round(median(f_df[, paste(f_quantity)], na.rm=TRUE), digits=f_digits), " (", 
                     round(min(f_df[, paste(f_quantity)], na.rm=TRUE), digits=f_digits), " to ",
                     round(max(f_df[, paste(f_quantity)], na.rm=TRUE), digits=f_digits), ", ", 
                     length(na.omit(f_df[, paste(f_quantity)]) ), ")"
                     , sep="" )
    
    x1_pos <- 1
    # for each year...
    for (i in f_years) {
      x1_pos <- x1_pos + 1
      x1[x1_pos] <- paste(round(median(subset(f_df, year_survey == i)[, paste(f_quantity)], na.rm=TRUE), digits=f_digits), " (", 
                      round(min(subset(f_df, year_survey == i)[, paste(f_quantity)], na.rm=TRUE), digits=f_digits), " to ",
                      round(max(subset(f_df, year_survey == i)[, paste(f_quantity)], na.rm=TRUE), digits=f_digits), ", ", 
                      length(na.omit(subset(f_df, year_survey == i)[, paste(f_quantity)]) ), ")", sep="" )
      }
    return(x1)
} 


#.........................................................................................
### Function to plot histograms of variables
#.........................................................................................  

f_hist <- function(f_var, f_data, f_lims) {
    
  plot <- ggplot(f_data)
      
    # if the variable has >= 20 unique values...
      if (length(unique(na.omit(f_data[, f_var]))) >= 20) {
        plot <- plot + geom_histogram(aes(x = as.numeric(f_data[, f_var]) ), 
          color="seagreen", fill="seagreen3", alpha = 0.5 ) +
          theme_bw() + xlab(f_var) + scale_x_continuous(expand = c(0, 0), limits = f_lims )
      }
 
    # otherwise...
      if (length(unique(na.omit(f_data[, f_var]))) < 20) {
        plot <- plot + geom_histogram(aes(x = as.factor(f_data[, f_var]) ), stat="count", 
          color="seagreen", fill="seagreen3", alpha = 0.5) +
          theme_bw() + xlab(f_var)
      }
        
    print(plot)
  }


#.........................................................................................
### Function to standardise livelihood type nomenclature
#.........................................................................................   

f_liv <- function(f_ts, f_livelihood_substrings) {
  # Agriculturalists
  if (length(grep(paste(f_livelihood_substrings$agriculturalists, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[2]) )
  
  # Pastoralists (slightly different to avoid confusion with agropastoralists)
  if (length(grep(paste(f_livelihood_substrings$pastoralists, collapse="|"), f_ts )) > 0  & 
    length(grep(paste(f_livelihood_substrings$agropastoralists, collapse="|"), f_ts )) == 0
    ) 
    return( paste(names(livelihood_substrings)[3]) )
  
  # Agropastoralists       
  if (length(grep(paste(f_livelihood_substrings$agropastoralists, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[4]) )
  
  # Riverine        
  if (length(grep(paste(f_livelihood_substrings$riverine, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[5]) )
  
  # Fishing       
  if (length(grep(paste(f_livelihood_substrings$fishing, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[6]) )
  
  # Urban        
  if (length(grep(paste(f_livelihood_substrings$urban, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[7]) )
  
  # Displaced        
  if (length(grep(paste(f_livelihood_substrings$displaced, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[8]) )
  
  # Refugee       
  if (length(grep(paste(f_livelihood_substrings$refugee, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[9]) )
  
  # Otherwise return NA
  if (length(grep(paste(unlist(f_livelihood_substrings), collapse="|"), f_ts )) == 0 ) 
    return( NA )
    
}      



#.........................................................................................
### Function to execute sections of code 
  # (based on https://stackoverflow.com/questions/26245554/execute-a-set-of-lines-from-another-r-file )
#.........................................................................................

f_source_part <- function(f_script, f_start_tag, f_end_tag) {

  # Identify lines with start and end tags
  st <- grep(f_start_tag, f_script)
  en <- grep(f_end_tag, f_script)
  
  # Set up a connection
  tc <- textConnection(f_script[(st + 1):(en - 1)])
  
  # Run the script
  source(tc, echo = TRUE)
  
  # Close the connection
  close(tc)
}



  
#.........................................................................................
### ENDS
#.........................................................................................
