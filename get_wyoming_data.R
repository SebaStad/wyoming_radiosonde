# Packages ----------------------------------------------------------------
# can be installed with

# install.packages(c("glue", "tidyverse", "zoo"))

library(glue)
library(tidyverse)
library(zoo)

# Functions ---------------------------------------------------------------

interpol_orig <- function(.data, max_HGHT = 3000, use_all = T){
  
  lin_inter <- data.frame(HGHT = seq(max_HGHT))
  
  min_hght <- min(.data$HGHT, na.rm=T)
  
  lin_inter %>%
    dplyr::full_join(y= .data, by = "HGHT") %>%
    dplyr::filter(HGHT >= min_hght &
                    HGHT<= .data$HGHT[HGHT>max_HGHT][1]) %>%
    dplyr::arrange(HGHT) %>%
    
    dplyr:: mutate(THTA = na.approx(THTA),
                   MIXR = na.approx(MIXR),
                   UV2 = na.approx(UV2),
                   THTV = na.approx(THTV),
                   TEMP = na.approx(TEMP+273.15),
                   DIR = atan2(sin(DRCT*(pi/180)), cos(DRCT*(pi/180)))%%(pi*2),
                   DIR = na.approx(DIR),
                   DRCT = DIR*(180/pi))%>%
    dplyr::select(HGHT, THTA, MIXR, UV2, THTV, TEMP, DIR, DRCT) %>%
    dplyr::filter(HGHT<=max_HGHT)
}


get_months_wyoming <- function(year, month, region, station,
                              write_monthly = TRUE, write_daily = TRUE, leap = FALSE,
                              read_monthly_from_file = FALSE, do_interpolation = TRUE){
 months_vec <- c("january", "february", "march", "april", "may", "june", "july",
                 "august", "september", "october", "november", "december")
 outlist <- list()
 
 for(loop_month in month){
  outlist[[months_vec[loop_month]]] <- list()
  max <- 30
  if(loop_month==2){
    max <- 28
    if(leap) max <- 29
  } else if(loop_month%in%c(1,3,5,7,8,10,12)){
    max <- 31
  }
  
  if(read_monthly_from_file){
    parse <- readr::read_file(paste0( getwd(),"/", station, "/", months_vec[loop_month], ".txt"))
    p_split <- str_split(parse,  "Observations at")
    
  } else {
    mm <- formatC(loop_month, width = 2, format = "d", flag = "0")
  
    url <- glue::glue("http://weather.uwyo.edu/cgi-bin/sounding?region={region}&TYPE=TEXT%3ALIST&YEAR={year}&MONTH={mm}&FROM=0100&TO={max}18&STNM={station}")
    temp <- tempfile()
    download.file(url, temp)
    parse <- readr::read_file(temp) %>% 
      gsub("<.*?>", "", .) 
  
    if(write_monthly){
      dir.create( paste0( getwd(),"/", station))
      write_file(parse, paste0( getwd(),"/", station, "/", months_vec[loop_month], ".txt"))
    }
    
    p_split <- str_split(parse, "Observations at")
  }
  
  for(i in 2:length(p_split[[1]])){

    if(is.na(p_split[[1]][i])) break
    
    hour <- as.numeric(substr(p_split[[1]][i],2,3))
    day  <- as.numeric(substr(p_split[[1]][i],6,7))

    
    headers <- gsub(pattern = "-", replacement = "", p_split[[1]][i]) %>% 
      read_table(., col_names = T) %>% 
      .[1,] %>% 
      unlist(.) %>% 
      str_split(., pattern = "[ \t]+") %>% 
      unlist(.)
    
    data2 <-  gsub("<.*?>", "", p_split[[1]][i]) %>% 
      read_table(., skip = 5) %>%
      separate(names(.), headers, sep = "[ \t]+") %>%
      na.omit() %>%
      mutate_if(is.character, as.numeric) %>% 
      na.omit()
    
    check_wrong <- FALSE
    check_thtv  <- FALSE
    check_drct  <- FALSE
    check_height <- FALSE
    check_mixr <- FALSE
    
    if("SKNT"%in%names(data2)){
      data2 <- data2 %>% 
        mutate(UV2 = SKNT*0.5144)
      
    } else {
      data2 <- data2 %>% 
        mutate(SKNT = 0) %>% 
        mutate(UV2 = SKNT*0.5144) %>% 
        mutate(THTA = 0)
      check_wrong <- TRUE
    }
    
    if(!"THTV"%in%names(data2)){
      data2 <- data2 %>% 
        mutate(THTV = 0)
      
      check_thtv <- TRUE
    }
    
    if(!"DRCT"%in%names(data2)){
      data2 <- data2 %>% 
        mutate(DRCT = 0)
      
      check_drct <- TRUE
    }
    
    if(data2$HGHT[1] < 0){
      check_height <- TRUE
      
    }
    
    if(!"MIXR"%in%names(data2)){
      data2 <- data2 %>% 
        mutate(MIXR=0)
      
      check_mixr <- TRUE
    }
    
    
    if(write_daily){
      dir.create(  glue("{getwd()}/{station}/daily/{loop_month}"), recursive = T)
      write.table(data2, glue("{getwd()}/{station}/daily/{loop_month}/{day}_{loop_month}_2010_{hour}Z.txt"),
                  quote = F, sep = "\t",
                  row.names = F, col.names = T
      )
    }
    if(do_interpolation){
      data2 <- data2 %>%         
        interpol_orig(max_HGHT = 13000) %>% 
        filter(HGHT<12000)
    }
    
    if(check_wrong){
      data2 <- data2 %>% 
        slice(0)
    }
    if(check_thtv & !check_wrong){
      data2$THTV <- NA
    }
    if(check_drct & !check_wrong){
      data2$DRCT <- NA
      data2$DIR  <- NA
    }
    if(check_mixr & !check_wrong){
      
      data2$MIXR <- NA
    }
    
    if(check_height){
      data2 <- data2 %>% 
        slice(0)
    }
    outlist[[months_vec[loop_month]]][[glue("{day}_{hour}")]] <- data2
  }
 }
  return(outlist)
}

create_mean_vertical_profile <- function(sounding_data, timesteps, min_height, parameter = "THTA"){
  
  sel_data <- list()
  cntr <- 1
  
  
  for(j in names(sounding_data)){
    for(i in seq_along(sounding_data[[j]])){
      time_name <- names(sounding_data[[j]])[i]
      time_name <- unlist(strsplit(time_name, "_"))[2]
      
      if(as.numeric(time_name)%in%timesteps){
        sel_data[[cntr]] <- sounding_data[[j]][[i]]
        cntr <- cntr +1
      }
      
    }
  }
  
  
  pot_temp <- do.call(cbind, lapply(sel_data, FUN = function(x){
    if(min(x$HGHT,na.rm=T)==min_height){
      x[[parameter]]
    }
  })  
  )
  
  height <- sel_data[[1]]$HGHT
  
  mean_pot_temp <- apply(pot_temp, 1, mean)
  std_pot_temp  <- apply(pot_temp,1,sd)
  
  
  return(data.frame(height = height, 
                    theta = mean_pot_temp,
                    std = std_pot_temp)
  )
}

# DISCLAIMER
# There will be warnings. None are critical whatsoever.
#


# get_months_wyoming needs following input:

# year
# months (in numbers, can be 1:12)
# region; "europe" for Berlin
#         "mideast" for abu dhabi
# station: station number as per website
#         Berlin Lindenberg:  10393
#         Abu Dhabi:          41217
# write_monthly
#         writes text files with the measurement for each month
# write_daily
#         writes text files for each measurement
# It saves the monthly data in the workding directory + station_number
# It saves the "hourly" data in workking directory + station_number + daily + month
# leap
#         whether the year is a leap year or not. if TRUE and it is not a leap year, script will fail
# read_monthly_from_file
#         If you already downloaded and saved the monthly file, you can just read it from the disk

# Set the correct working directory if you want to save the data.
setwd("G:/15 InCityTakeOff/AbuDhabi")


abu_dhabi <-  get_months_wyoming(2010, 1:2, "mideast", 41217, write_monthly = FALSE, write_daily = FALSE)
berlin    <-  get_months_wyoming(2010, 1:3, "europe",  10393, write_monthly = FALSE, write_daily = FALSE,
                                 read_monthly_from_file = T)


# create_mean_vertical_profile needs following input
# sounding_data
#         the variable created by get_months_wyoming
# timesteps
#         at which timesteps should the evaluation be?
#         Berlin has 00Z, 06Z, 12Z, 18Z
#         Abu Dhabi has 00Z, 12Z
#         It can be one value (e.g. 0), or multiple ones (e.g. c(0,12))
# min_height
#         the minimal height of the measurements
#         Berlin is 115
#         Abu Dhabi is 27

abu_dhabi_10 <- create_mean_vertical_profile(abu_dhabi, 12, 27, parameter = "THTA")

# simple example plot
ggplot(abu_dhabi_10, aes(y = theta, x = height ))+
  geom_path() + 
  geom_ribbon(aes(ymin = theta-std, ymax = theta+std))+
  coord_flip()

# You can select every variable available
names(berlin$january$`1_0`)

berlin_10 <- create_mean_vertical_profile(berlin, c(6,12), 115, parameter = "UV2")

ggplot(berlin_10, aes(y = theta, x = height )) +
  geom_path() + 
  geom_ribbon(aes(ymin = theta-std, ymax = theta+std))+
  coord_flip()



# Write as CSV file:

write.csv2(abu_dhabi_10, "abu_dhabi_mean.csv", quote = F, row.names = F)


abu_dhabi_small <- abu_dhabi_10 %>% 
  filter(height<1500)

write.csv2(abu_dhabi_small, "abu_dhabi_mean_mean_1500.csv", quote = F, row.names = F)









