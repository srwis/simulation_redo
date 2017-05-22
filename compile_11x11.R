
library(jsonlite)

library(stringr)

library(dplyr)



compile <- function(cur.dir,csv){
  jsons <- list.files(path = cur.dir,
                      pattern = '*.json', recursive = TRUE)
  
  #create empty data.frame with variable names
  names = c("FILE", "BUSTED.LR", "BUSTED.SRV.LR", "BUSTED.omega3.MLE", "BUSTED.SRV.omega3.MLE", "BUSTED.omega3.prop",
            "BUSTED.SRV.omega3.prop", 'CV.SRV', 'BUSTED.P', 'BUSTED.SRV.P','BUSTED.AICc','BUSTED.SRV.AICc',
            'BUSTED.treelength' ,'BUSTED.SRV.treelength', 'Sites', 'Sequences', 
            'BUSTED.omega1.MLE','BUSTED.SRV.omega1.MLE', 'BUSTED.omega1.prop','BUSTED.SRV.omega1.prop',
            'BUSTED.omega2.MLE','BUSTED.SRV.omega2.MLE', 'BUSTED.omega2.prop','BUSTED.SRV.omega2.prop', 'SRV.alpha3.MLE',
            'SRV.alpha3.prop','SRV.alpha1.MLE','SRV.alpha1.prop','SRV.alpha2.MLE','SRV.alpha2.prop')
  
  classes = c("character", "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
              "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
              "numeric","numeric","numeric")
  df =read.table(text="", col.names = names, colClasses = classes)
  for (i in  seq(from=1, to=length(jsons), by=2)){
    filepath = paste(cur.dir,jsons[i], sep="")
    test = filepath %>% readLines() %>% gsub(x=.,pattern="nan",replacement ='"NA"') %>% fromJSON()
    
    
    
    FILE = jsons[i]
    Sites = length(test$profiles$unconstrained)
    tree_string = test$fits$`Unconstrained model`$`tree string`
    x= tree_string %>% str_replace_all("\\(",":") %>% str_replace_all("\\)",":") %>%     str_replace_all(",",":") %>% str_split(":")
    x= unlist(x)
    x =x[x !=""]
    br_len = matrix(x,ncol = 2,  byrow = T)
    colnames(br_len) = c("Branch", "Length")
    
    Sequences = sum(grepl("Node*", br_len[,1]) == FALSE)
    
    if (grepl("BUSTED-SRV",jsons[i])){
      filepath = paste(cur.dir,jsons[i], sep="")
      
      test = filepath %>% readLines() %>% gsub(x=.,pattern="nan",replacement ='"NA"') %>% fromJSON()      
      
      
      BUSTED.SRV.P = test$`test results`$p
      BUSTED.SRV.LR =test$`test results`$LR
      BUSTED.SRV.AICc = test$fits$`Unconstrained model`$`AIC-c`
      BUSTED.SRV.treelength = test$fits$`Unconstrained model`$`tree length`
      
      #OMEGA values for BUSTED.SRV
      BUSTED.SRV.omega3.MLE = test$fits$`Unconstrained model`$`rate distributions`$FG[3,1]
      BUSTED.SRV.omega3.prop = test$fits$`Unconstrained model`$`rate distributions`$FG[3,2]
      BUSTED.SRV.omega2.MLE = test$fits$`Unconstrained model`$`rate distributions`$FG[2,1]
      BUSTED.SRV.omega2.prop = test$fits$`Unconstrained model`$`rate distributions`$FG[2,2]
      BUSTED.SRV.omega1.MLE = test$fits$`Unconstrained model`$`rate distributions`$FG[1,1]
      BUSTED.SRV.omega1.prop = test$fits$`Unconstrained model`$`rate distributions`$FG[1,2]
      #ALPHA values for BUSTED.SRV
      SRV.alpha3.MLE = test$fits$`Unconstrained model`$`rate distributions`$SRV[3,1]
      SRV.alpha3.prop = test$fits$`Unconstrained model`$`rate distributions`$SRV[3,2]
      SRV.alpha2.MLE = test$fits$`Unconstrained model`$`rate distributions`$SRV[2,1]
      SRV.alpha2.prop = test$fits$`Unconstrained model`$`rate distributions`$SRV[2,2]
      SRV.alpha1.MLE = test$fits$`Unconstrained model`$`rate distributions`$SRV[1,1]
      SRV.alpha1.prop = test$fits$`Unconstrained model`$`rate distributions`$SRV[1,2]
      
      mom2 = SRV.alpha3.MLE^2*SRV.alpha3.prop+ SRV.alpha1.MLE^2*SRV.alpha1.prop+ SRV.alpha2.MLE^2*SRV.alpha2.prop
      mean= SRV.alpha3.MLE*SRV.alpha3.prop+ SRV.alpha1.MLE*SRV.alpha1.prop+ SRV.alpha2.MLE*SRV.alpha2.prop
      CV.SRV = sqrt(mom2-mean^2)/mean
      
    }
    if (grepl("BUSTED.json",jsons[i+1])){
      filepath = paste(cur.dir,jsons[i+1], sep="")
      
      test = filepath %>% readLines() %>% gsub(x=.,pattern="nan",replacement ='"NA"') %>% fromJSON()
      #print(filepath)    
      BUSTED.P = test$`test results`$p
      BUSTED.LR = test$`test results`$LR
      BUSTED.AICc = test$fits$`Unconstrained model`$`AIC-c`
      BUSTED.treelength = test$fits$`Unconstrained model`$`tree length`
      
      #OMEGA values for BUSTED
      BUSTED.omega3.MLE = test$fits$`Unconstrained model`$`rate distributions`$FG[3,1]
      BUSTED.omega3.prop = test$fits$`Unconstrained model`$`rate distributions`$FG[3,2]
      BUSTED.omega2.MLE = test$fits$`Unconstrained model`$`rate distributions`$FG[2,1]
      BUSTED.omega2.prop = test$fits$`Unconstrained model`$`rate distributions`$FG[2,2]
      BUSTED.omega1.MLE = test$fits$`Unconstrained model`$`rate distributions`$FG[1,1]
      BUSTED.omega1.prop = test$fits$`Unconstrained model`$`rate distributions`$FG[1,2]
      #ALPHA values for BUSTED
      #       BUSTED.alpha3.MLE = test$fits$`Unconstrained model`$`rate distributions`$SRV[3,1]
      #       BUSTED.alpha3.prop = test$fits$`Unconstrained model`$`rate distributions`$SRV[3,2]
      #       BUSTED.alpha2.MLE = test$fits$`Unconstrained model`$`rate distributions`$SRV[2,1]
      #       BUSTED.alpha2.prop = test$fits$`Unconstrained model`$`rate distributions`$SRV[2,2]
      #       BUSTED.alpha1.MLE = test$fits$`Unconstrained model`$`rate distributions`$SRV[1,1]
      #       BUSTED.alpha1.prop = test$fits$`Unconstrained model`$`rate distributions`$SRV[1,2]
      
    }
    #print(FILE)
    df[nrow(df)+1,] <- c(FILE, BUSTED.LR, BUSTED.SRV.LR, BUSTED.omega3.MLE, BUSTED.SRV.omega3.MLE, BUSTED.omega3.prop,
                         BUSTED.SRV.omega3.prop, CV.SRV, BUSTED.P, BUSTED.SRV.P,BUSTED.AICc,BUSTED.SRV.AICc,
                         BUSTED.treelength ,BUSTED.SRV.treelength, Sites, Sequences, 
                         BUSTED.omega1.MLE,BUSTED.SRV.omega1.MLE, BUSTED.omega1.prop,BUSTED.SRV.omega1.prop,
                         BUSTED.omega2.MLE,BUSTED.SRV.omega2.MLE, BUSTED.omega2.prop,BUSTED.SRV.omega2.prop, SRV.alpha3.MLE,
                         SRV.alpha3.prop,SRV.alpha1.MLE,SRV.alpha1.prop,SRV.alpha2.MLE,SRV.alpha2.prop)
    
  }
  df[,2:30]=as.numeric(unlist(df[,2:30]))
  write.csv(file = csv, x = df, row.names= F)
  # return(df)
}

simulation_inputs <- function(dir,csv){
  require("stringr")
  require("jsonlite")
  require("dplyr")
  list = list.files(path = dir, recursive = T, pattern ="^([^.]+)$")
  #set up the empty data frame
  #names = c("Sites" ,"Cat", "Omega 1 value", "Omega 2 value", "Omega 3 value", "Omega 1 prop", "Omega 2 prop", "Omega 3 prop",
            #"Alpha 1 value", "Alpha 2 value", "Alpha 3 value","Alpha 1 prop", "Alpha 2 prop", "Alpha 3 prop")
  names = c("Sites","File" , paste("Omega",seq(from=1, to = 11, by =1), "value" ), paste("Omega",seq(from=1, to = 11, by =1), "prop" ),
            paste("Alpha",seq(from=1, to = 11, by =1), "value" ), paste("Alpha",seq(from=1, to = 11, by =1), "prop" ),"CV.SRV")
  setup.tab = read.table(text = "",col.names = names)
  #loop thru each file to get info in correct format
  for(i in seq(from = 1, to= length(list))){
    x=readLines(paste(dir,list[i], sep = "/"))
        #making this a readable json
    
    x1 = x[2:(length(x)-1)]  %>% gsub(x=.,pattern="\\{",replacement ='\\[') %>% gsub(x=.,pattern ="\\}", replacement ="\\]")
    x1 = c(x[1],x1,x[length(x)])
    num_rates=as.numeric(str_extract(x1[4], "[0-9]."))
    
    end_1 = 6+(num_rates-2)
    start_2 = end_1+5
    end_2 = start_2+num_rates-2
    x1[c(6:end_1,start_2:end_2)]  = x1[c(6:end_1,start_2:end_2)] %>%  gsub(x=.,pattern ="\\]", replacement ="\\],")

    r= fromJSON(x1)
    omega_rates = r$`omega distribution`[,1]
    names(omega_rates)= paste("Omega",seq(from=1, to = r$`omega rate count`, by =1), "value" )
    omega_weights = r$`omega distribution`[,2]
    names(omega_weights) = paste("Omega",seq(from=1, to = r$`omega rate count`, by =1), "prop" )
    
    Alpha_rates = r$`alpha distribution`[,1]
    names(Alpha_rates)= paste("Alpha",seq(from=1, to = r$`alpha rate count`, by =1), "value" )
    Alpha_weights = r$`alpha distribution`[,2]
    names(Alpha_weights) = paste("Alpha",seq(from=1, to = r$`alpha rate count`, by =1), "prop" )
    
    
    
    
    
    mom2 = sum(Alpha_rates^2*Alpha_weights)
    
    mean = sum(Alpha_rates*Alpha_weights)
    
    CV.SRV = sqrt(mom2-mean^2)/mean
    setup.tab[nrow(setup.tab)+1,] = c(r$sites,list[i], omega_rates,omega_weights,Alpha_rates,Alpha_weights,CV.SRV)
    
  }
  #setup.tab =select(setup.tab, File,Sites,  Omega.1.value, Omega.1.prop, Omega.2.value, Omega.2.prop,Omega.3.value, Omega.3.prop,
                   # Alpha.1.value, Alpha.1.prop, Alpha.2.value, Alpha.2.prop,Alpha.3.value, Alpha.3.prop)
  #setup.tab[,3:14] = as.numeric(unlist(setup.tab[,3:14]))

  return(setup.tab)
  #write.csv(file = csv, x = setup.tab, row.names= F)
}

add_truth <- function(dat, truth){
  dat = dat %>% mutate(True.CV = 1337)
  dat = dat %>% mutate(True.omega3.value = 8008)
  dat = dat %>% mutate(True.omega1.value = 8008)
  dat = dat %>% mutate(True.omega2.value = 8008)
  dat = dat %>% mutate(True.alpha3.value = 8008)
  dat = dat %>% mutate(True.alpha2.value = 8008)
  dat = dat %>% mutate(True.alpha1.value = 8008)
  dat = dat %>% mutate(True.omega3.prop = 8008)
  dat = dat %>% mutate(True.omega1.prop = 8008)
  dat = dat %>% mutate(True.omega2.prop = 8008)
  dat = dat %>% mutate(True.alpha3.prop = 8008)
  dat = dat %>% mutate(True.alpha2.prop = 8008)
  dat = dat %>% mutate(True.alpha1.prop = 8008)
  for( i in seq(from = 1, to = nrow(truth), by = 1)){
    temp = which(str_detect(dat$FILE,truth$File[i]))
    dat$True.CV[temp] = truth$CV.SRV[i] 
    dat$True.omega3.value[temp] = truth$Omega.3.value[i]
    dat$True.omega2.value[temp] = truth$Omega.2.value[i]
    dat$True.omega1.value[temp] = truth$Omega.1.value[i]
    dat$True.alpha1.value[temp] = truth$Alpha.1.value[i]
    dat$True.alpha2.value[temp] = truth$Alpha.2.value[i]
    dat$True.alpha3.value[temp] = truth$Alpha.3.value[i]
    dat$True.omega3.prop[temp] = truth$Omega.3.prop[i]
    dat$True.omega2.prop[temp] = truth$Omega.2.prop[i]
    dat$True.omega1.prop[temp] = truth$Omega.1.prop[i]
    dat$True.alpha1.prop[temp] = truth$Alpha.1.prop[i]
    dat$True.alpha2.prop[temp] = truth$Alpha.2.prop[i]
    dat$True.alpha3.prop[temp] = truth$Alpha.3.prop[i]
  }
  return(dat)
}



process_dat <- function(dir, basename){
  temp = paste(dir,basename,"_Truth.csv", sep = "")
  simulation_inputs(dir,temp)
  truth = read.csv(temp, as.is = T)
  temp = paste(dir,basename,"_results.csv", sep = "")
  compile(dir,temp)
  dat = read.csv(temp, as.is = T)
  dat = add_truth(dat, truth)
  #dat = mutate(dat, Cat = str_extract(dat$FILE, "YesYes|YesNo|NoNo|NoYes"))
  dat$True.CV = round(dat$True.CV, 3)
  write.csv(file = paste(dir,basename,"_processed.csv", sep = ""),x = dat, row.names = F)
}



dir= "1000_Codons_11x11"
basename = "All_5_8_17"
all.dat = process_dat(dir, basename)
