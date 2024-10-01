#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

nombrito <- paste0("fabm_",as.character(args[1]))

if (require(stringr)==FALSE) install.packages("stringr")
library(stringr)

indir<-paste0(getwd(),"/")
outdir<-paste0(getwd(),"/")

#########################################################################
## COMPOSE & SAVE fabm.yaml with several size-classes per functional type
#########################################################################

## Define size-classes in an octave scale (2^n), for details: 
## "Dealing with size-spectra: some conceptual and mathematical problems" JM Blanco, F Echevarría, CM García - Scientia Marina 1994
clases<-c(-8:32)   # n
volumes<-2^clases  # biovolume in um^3
diameters<-((volumes*6)/pi)^(1/3) # equivalent spherical diameter (ESD)

## Define size range (bins) for each functional type:
clasesB1<-c(-5)
clasesP6<-c(-4:-3)
clasesP9<- c(-2:0) 
clasesP3<-  c(0:2)
clasesP7<-  c(2:7)  
clasesP8<-  c(2:7)
clasesP2<- c(4:13)
clasesP5<- c(6:15)
clasesP1<-c(11:20)
clasesP4<-c(11:20)
clasesZ6<-c(5:11)
clasesZ5<-c(9:18)
clasesZ4<-c(15:24)
clasesZ3<-c(22:31)

## Count how many types (all, phyto and zoo) fall in each bin:
all_groups<-c(clasesB1,
              clasesP1, clasesP2, clasesP3, clasesP4, clasesP5, clasesP6, clasesP7, clasesP8, clasesP9,
              clasesZ6, clasesZ5, clasesZ4, clasesZ3)
histo<-hist(all_groups, breaks=clases, plot=F)
clases_init<-histo$breaks[-1]
counts_init<-histo$counts

phyto_groups<-c(clasesP1, clasesP2, clasesP3, clasesP4, clasesP5, clasesP6, clasesP7, clasesP8, clasesP9)
histo<-hist(phyto_groups, breaks=clases, plot=F)
clases_init<-histo$breaks[-1]
counts_init_phyto<-histo$counts

zoo_groups<-c(clasesZ6, clasesZ5, clasesZ4, clasesZ3)
histo<-hist(zoo_groups, breaks=clases, plot=F)
clases_init<-histo$breaks[-1]
counts_init_zoo<-histo$counts

#cbind(counts_init,counts_init_phyto,counts_init_zoo)



## Initialization biomass (fabm) and reference state (gotm)
## if reference state is given (in gotm.yaml for relaxation) the value in fabm.yaml is not used
## Several options tested below
###############################################################################################

# Initialize each size-class to the same biomass c=2 mg m-3 (all phyto and zoo together)
# Note!!! B1 does not get the initialization value
c_init=2                                                    
# Initialize each type to the same biomass 
c_init_phyto=0.5
c_init_zoo=0.1
# Initialize each size-class (bin) to the same biomass, phyto and zoo separately
# realistic initialization values (similar to Serra-Pompei et al 2020 c=5 mg m-3)
c_init_phyto=1
c_init_zoo=0.5
# very small initialization values (in the order of magnitude Bruggeman & Kooijman 2007 c=0.0096 mg m-3 per type)
c_init_phyto=0.01
c_init_zoo=0.005
# intermediate initialization values (null net flow, on average, with migration in and out)
c_init_phyto=0.1
c_init_zoo=0.08

# Choose one
equal_initialization=FALSE
separate_phytozoo=TRUE

###############################################




## Read the fabm.yaml template that contains the 9PFTs and 4Zoo
## several text tags have been included to locate the different "lego bricks"
modulo_productor <- readLines(paste(indir,"fabm_multispectral_9PFTs_template.yaml",sep=""))
## Locate the commom part: lightspectral, bethicLayer, initialization for pelagic_base, PelagicCSYS and PelOxygen, B1, PelChem and CalciteDissolution
inicio<-grep(paste("start","base",sep="_"), modulo_productor)
final<-grep(paste("end","base",sep="_"), modulo_productor)



# START the fabm.yaml
sink(paste(outdir,nombrito,".yaml",sep=""))
# Write the common part
write.table(modulo_productor[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
            sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)



### Create BRICKS for PRODUCERS


####### P1 ########
old_name<-"P1"
old_long_name<-"diatoms"
clasesP<-clasesP1      # ESD = (((2^clasesP)*6)/pi)^(1/3)
clasesN<-str_replace(clasesP,pattern="-","m")   # the name of the instance cannot have "-", for names m5 = -5
# parameters to modify
parameters<-c("p_sum","p_srs","p_res","p_qlcPPY",   "p_qup",  "p_qun",   "p_qplc",  "p_qpcPPY", "p_qnlc",  "p_qncPPY", "p_quantum_yield")
old_value <-c(  "3.7", "0.06", "5.0",   "0.039",    "0.0046", "0.025",   "0.00057", "0.000786", "0.00687", "0.0126",    "0.504e-3")
digito<-c(3,3,2,3,  4,4, 6,6,6,6, 5)
# Extract module block (template) for Px
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)     
         # modulo_productor[inicio:final]

###################  
lower_limit<-2^clasesP      
upper_limit<-2^(clasesP+1)
low_diameter<-((lower_limit*6)/pi)^(1/3)

## Compute parameter values (in case allometry is taxon-dependent)

  # Maximum growth rate 
  #p_sum<- (10^0.54)*lower_limit^(-0.15) # (Tang 1995) d-1
  #p_sum<- (10^0.73)*lower_limit^(-0.17) # (Tang 1995) d-1 for PFG A (diatoms)
  p_sum<- 3.8*lower_limit^(-0.17)        # b from Tang 1995, a from Irwin 2006
  
  # Basal metabolism 
  # p_srs<- 10^(-0.063)*lower_limit^0.008 ?? (Ray 2001) d-1
  p_srs<- 0.063-0.008*log10(lower_limit)     # Shimoda 2016
  
  # Carbon content in cell, pgC cell-1 
  #carbon<-(10^(-0.29))*lower_limit^0.76    # Mullin 1996
  #carbon<-(10^(-0.665))*lower_limit^0.939  # Menden-Deuer and Lessard, 2000 (overall)
  #carbon<-(10^(-0.541))*lower_limit^0.811  # Menden-Deuer and Lessard, 2000 (diatoms)
  carbon<-(10^(-0.933))*lower_limit^0.881  # Menden-Deuer and Lessard, 2000 (diatoms>3000)
  
  # Minimum P:C quota (Grover 1989), fmol cell-1
  QPmin<-(10^(-1.04))*lower_limit^0.714 
  p_qplc<-(QPmin/carbon)*1e9*1e-12 # mmolP mgC-1
  # Maximal P:C quota (Grover 1989), fmol cell-1, divided by luxury storage to get optimum quota
  QPmax<-(10^(-0.29))*lower_limit^0.767 
  p_qpcPPY<-((QPmax/carbon)*1e9*1e-12)/2 # mmolP mgC-1  
  
  # Minimum N:C quota (Lichtman et al 2007), umol cell-1
  QNmin<-(1.36e-9)*lower_limit^0.77 
  p_qnlc<-(QNmin/carbon)*1e9*1e-3 # mmolN mgC-1
  # Maximal N:C quota (Montagnes & Franklin 2001), umol cell-1, divided by luxury storage to get optimum quota
  QNmax<-(4.64e-9)*lower_limit^0.81 
  p_qncPPY<-((QNmax/carbon)*1e9*1e-3)/2 # mmolN mgC-1    
  
  # Settling velocity, m d−1
  #p_res<-0.029*carbon^0.42       # Moloney and Field, 1991  
  #p_res<-0.036*lower_limit^0.43  # Durante et al, 2019  
  p_res<-0.024*lower_limit^0.37   # Durante et al, 2019 (cylinders)
  
      # P affinity  # Edwards 2012
      aP<-(10^(-8.1))*lower_limit^0.73              # L cell-1 d-1
      p_qup<-(aP/carbon)*1e-3*1e9/4                 # m3 mgC-1 d-1
      
      # N affinity # Edwards 2012
      aN<-(10^(-8.2))*lower_limit^0.75              # L cell-1 d-1
      p_qun<-(aN/carbon)*1e-3*1e9                   # m3 mgC-1 d-1
      
      ## Original values in BFM for quota
      #p_qplc<-rep(0.000432,length(lower_limit))
      #p_qpcPPY<-rep(0.000786,length(lower_limit))
      #p_qnlc<-rep(0.00687,length(lower_limit))
      #p_qncPPY<-rep(0.0126,length(lower_limit))   
      
      # The only allometry in use is optimal P:C quota, the rest is scaled as previously in BFM
      p_qplc<-p_qpcPPY*0.55
      p_qncPPY<-(p_qpcPPY*0.0126)/0.000786        
      p_qnlc<-p_qncPPY*0.55

  # Chl-a:C
  p_qlcPPY<-rep(0.02,length(lower_limit))
  
  # Photochemical efficiency
  p_quantum_yield<-rep(0.480e-3,length(lower_limit))
      
## For each instance (j): substitute old parameter value (i) by new one and write the brick in the fabm.yaml      
          for (j in c(1:length(clasesP))){    
            modulo_modificado<-modulo_productor
            diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
            nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesN[j],sep=""))
            modulo_modificado[inicio:final][3]<-nueva
            nueva<-str_replace(modulo_modificado[inicio:final][4], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
            modulo_modificado[inicio:final][4]<-nueva  
            
                      for (i in c(1:length(parameters))){
                        parametro<-get(parameters[i])
                      # Locate the parameter
                        filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
                        fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
                        cual<-which(!is.na(filas[,1])==T)
                      # Replace parameter value
                        nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
                        modulo_modificado[inicio:final][cual]<-nueva
                      } # end loop i parameters
            
            # Compute initialization values
            if (equal_initialization==T) {carbon<-c_init_phyto
            } else if (separate_phytozoo==T) { carbon<-c_init_phyto/counts_init_phyto[which(clases_init==clasesP[j])]
            } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
            nitro <-carbon*(0.1008/8)
            phosp <-carbon*(0.006288/8)
            silica <-carbon*(0.08/8)
            chloro<-carbon*(0.16/8)
            
            # Locate initialization chunk
            location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
            cual<-which(!is.na(location[,1])==T)
            # Replace initialization values
            nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "8.0", as.character(round(carbon,10)))
            modulo_modificado[inicio:final][cual+1]<-nueva
            nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.1008", as.character(round(nitro,10)))
            modulo_modificado[inicio:final][cual+2]<-nueva
            nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.006288", as.character(round(phosp,10)))
            modulo_modificado[inicio:final][cual+3]<-nueva
            nueva<-str_replace(modulo_modificado[inicio:final][cual+4], "0.08", as.character(round(silica,10)))
            modulo_modificado[inicio:final][cual+4]<-nueva
            nueva<-str_replace(modulo_modificado[inicio:final][cual+5], "0.16", as.character(round(chloro,10)))
            modulo_modificado[inicio:final][cual+5]<-nueva

            # Add to fabm.yaml
            write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
                        sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
              
          } # end loop j size_classes  
###################            

 
####### P4 ########
old_name<-"P4"
old_long_name<-"dinoflagellates"
clasesP<-clasesP4      # ESD = (((2^clasesP)*6)/pi)^(1/3)
clasesN<-str_replace(clasesP,pattern="-","m")
# parameters to modify
parameters<-c("p_sum","p_srs","p_res", "p_qlcPPY",    "p_qup",  "p_qun",   "p_qplc",    "p_qpcPPY", "p_qnlc",  "p_qncPPY",  "p_quantum_yield")
old_value <-c(  "1.5",  "0.1",  "2.5",    "0.015",    "0.0038", "0.025",   "0.0004288", "0.000786", "0.00687", "0.0126",    "0.617e-3")
digito<-c(3,3,2,3,  4,4, 6,6,6,6, 5)
# Extract module block for Px
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)
          # modulo_productor[inicio:final]
          
###################  
lower_limit<-2^clasesP    
upper_limit<-2^(clasesP+1)
  
  # Compute parameter values (in case allometry is taxon-dependent)
  
  # Maximum growth rate 
  #p_sum<- (10^0.54)*lower_limit^(-0.15)  # (Tang 1995) d-1
  p_sum<- 2.1*lower_limit^(-0.15)         # b from Tang 1995, a from Irwin 2006

  # Basal metabolism 
  #p_srs<- 10^(-0.063)*lower_limit^0.008 ?? (Ray 2001) d-1
  p_srs<- 0.063-0.008*log10(lower_limit) # Shimoda 2016
  
  # Carbon content per cell cell
  #carbon<-(10^(-0.29))*lower_limit^0.76    # Mullin 1996  pgC cell-1 
  #carbon<-(10^(-0.665))*lower_limit^0.939  # Menden-Deuer and Lessard, 2000 (overall)
  carbon<-(10^(-0.353))*lower_limit^0.864  # Menden-Deuer and Lessard, 2000 (dinoflagellates)

  # Minimum P:C quota (Grover 1989), fmol cell-1
  QPmin<-(10^(-1.04))*lower_limit^0.714 
  p_qplc<-(QPmin/carbon)*1e9*1e-12 # mmolP mgC-1
  # Maximal P:C quota (Grover 1989), fmol cell-1, divided by luxury storage to get optimum quota
  QPmax<-(10^(-0.29))*lower_limit^0.767 
  p_qpcPPY<-((QPmax/carbon)*1e9*1e-12)/2 # mmolP mgC-1  
  
  # Minimum N:C quota (Lichtman et al 2007), umol cell-1
  QNmin<-(1.36e-9)*lower_limit^0.77 
  p_qnlc<-(QNmin/carbon)*1e9*1e-3 # mmolN mgC-1
  # Maximal N:C quota (Montagnes & Franklin 2001), umol cell-1, divided by luxury storage to get optimum quota
  QNmax<-(4.64e-9)*lower_limit^0.81 
  p_qncPPY<-((QNmax/carbon)*1e9*1e-3)/2 # mmolN mgC-1  
  
  # Settling velocity, m d−1
  #p_res<-0.029*carbon^0.42        # Moloney and Field, 1991  
  p_res<-0.036*lower_limit^0.43    # Durante et al, 2019 
  
  # P affinity  # Edwards 2012
  aP<-(10^(-8.1))*lower_limit^0.73    # L cell-1 d-1
  p_qup<-(aP/carbon)*1e-3*1e9/3        # m3 mgC-1 d-1
  
  # N affinity # Edwards 2012
  aN<-(10^(-8.2))*lower_limit^0.75    # L cell-1 d-1
  p_qun<-(aN/carbon)*1e-3*1e9*(4/3)    # m3 mgC-1 d-1
  
  ## Original values in BFM
  #p_qplc<-rep(0.000432,length(lower_limit))
  #p_qpcPPY<-rep(0.000786,length(lower_limit))
  #p_qnlc<-rep(0.00687,length(lower_limit))
  #p_qncPPY<-rep(0.0126,length(lower_limit))   

  p_qplc<-p_qpcPPY*0.55
  p_qncPPY<-(p_qpcPPY*0.0126)/0.000786        
  p_qnlc<-p_qncPPY*0.55
  
  # Chl-a:C
  p_qlcPPY<-rep(0.02,length(lower_limit))
  
  # Photochemical efficiency
  p_quantum_yield<-rep(0.480e-3,length(lower_limit))
  
## For each instance (j): substitute old parameter value (i) by new one and write the brick in the fabm.yaml   
      for (j in c(1:length(clasesP))){    
        modulo_modificado<-modulo_productor
        diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
        nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesN[j],sep=""))
        modulo_modificado[inicio:final][3]<-nueva
        nueva<-str_replace(modulo_modificado[inicio:final][4], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
        modulo_modificado[inicio:final][4]<-nueva  
            
            for (i in c(1:length(parameters))){
              parametro<-get(parameters[i])
              # Locate the parameter
              filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
              fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
              cual<-which(!is.na(filas[,1])==T)
              # Replace parameter value
              nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
              modulo_modificado[inicio:final][cual]<-nueva
            } # end loop i parameters
            
        # Compute initialization values
        if (equal_initialization==T) {carbon<-c_init_phyto
        } else if (separate_phytozoo==T) { carbon<-c_init_phyto/counts_init_phyto[which(clases_init==clasesP[j])]
        } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
        nitro <-carbon*(0.1008/8)
        phosp <-carbon*(0.006288/8)
        chloro<-carbon*(0.16/8)
        # Locate initialization chunk
        location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
        cual<-which(!is.na(location[,1])==T)
        # Replace initialization values
        nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "8.0", as.character(round(carbon,10)))
        modulo_modificado[inicio:final][cual+1]<-nueva
        nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.1008", as.character(round(nitro,10)))
        modulo_modificado[inicio:final][cual+2]<-nueva
        nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.006288", as.character(round(phosp,10)))
        modulo_modificado[inicio:final][cual+3]<-nueva
        nueva<-str_replace(modulo_modificado[inicio:final][cual+4], "0.16", as.character(round(chloro,10)))
        modulo_modificado[inicio:final][cual+4]<-nueva

        # Add to fabm.yaml
        write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
                    sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
        
      } # end loop j size_classes  
###################            
  
  
  
####### P2 ########
old_name<-"P2"
old_long_name<-"prymnesiophyta"
clasesP<-clasesP2   # ESD = (((2^clasesP)*6)/pi)^(1/3)
clasesN<-str_replace(clasesP,pattern="-","m")
# parameters to modify
parameters<-c("p_sum","p_srs", "p_res","p_qlcPPY",   "p_qup",  "p_qun",   "p_qplc",   "p_qpcPPY", "p_qnlc",  "p_qncPPY",  "p_quantum_yield")
old_value <-c(  "2.6", "0.09",  "0.0",   "0.045",    "0.0033", "0.025",   "0.000352", "0.000556", "0.00687", "0.0126",    "0.563e-3")
digito<-c(3,3,2,3,  4,4, 6,6,6,6, 5)
# Extract module block for Px
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)
#modulo_productor[inicio:final]
  
###################  
lower_limit<-2^clasesP      
upper_limit<-2^(clasesP+1)
  
  # Compute parameter values (in case allometry is taxon-dependent)
  
  # Maximum growth rate 
  #p_sum<- (10^0.54)*lower_limit^(-0.15)  # (Tang 1995) d-1
  p_sum<- 1.2*lower_limit^(-0.10)         # Dut Cer 2020

  # Basal metabolism (Ray 2001) d-1
  # p_srs<- 10^(-0.063)*lower_limit^0.008 ??
  p_srs<- 0.063-0.008*log10(lower_limit) # Shimoda 2016
  
  # Carbon content in cell, pgC cell-1, Menden-Deuer and Lessard, 2000
  #carbon<-(10^(-0.29))*lower_limit^0.76    # Mullin 1996 pgC cell-1
  #carbon<-(10^(-0.665))*lower_limit^0.939  # Menden-Deuer and Lessard, 2000 (overall)
  carbon<-(10^(-0.642))*lower_limit^0.899  # Menden-Deuer and Lessard, 2000 (prymnesiophyceae)
  
  # Minimum P:C quota (Grover 1989), fmol cell-1
  QPmin<-(10^(-1.04))*lower_limit^0.714 
  p_qplc<-(QPmin/carbon)*1e9*1e-12 # mmolP mgC-1
  # Maximal P:C quota (Grover 1989), fmol cell-1, divided by luxury storage to get optimum quata
  QPmax<-(10^(-0.29))*lower_limit^0.767 
  p_qpcPPY<-((QPmax/carbon)*1e9*1e-12)/2 # mmolP mgC-1  
  
  # Minimum N:C quota (Lichtman et al 2007), umol cell-1
  QNmin<-(1.36e-9)*lower_limit^0.77 
  p_qnlc<-(QNmin/carbon)*1e9*1e-3 # mmolN mgC-1
  # Maximal N:C quota (Montagnes & Franklin 2001), umol cell-1, divided by luxury storage to get optimum quota
  QNmax<-(4.64e-9)*lower_limit^0.81 
  p_qncPPY<-((QNmax/carbon)*1e9*1e-3)/2 # mmolN mgC-1 
    
  # Settling velocity, m d−1
  #p_res<-0.029*carbon^0.42      # Moloney and Field, 1991  
  p_res<-0.036*lower_limit^0.43  # Durante et al, 2019 
  
  # P affinity  # Edwards 2012
  aP<-(10^(-8.1))*lower_limit^0.73    # L cell-1 d-1
  p_qup<-(aP/carbon)*1e-3*1e9/4       # m3 mgC-1 d-1
  
  # N affinity # Edwards 2012
  aN<-(10^(-8.2))*lower_limit^0.75    # L cell-1 d-1
  p_qun<-(aN/carbon)*1e-3*1e9         # m3 mgC-1 d-1
  
  ## Original values in BFM
  #p_qplc<-rep(0.000432,length(lower_limit))
  #p_qpcPPY<-rep(0.000786,length(lower_limit))
  #p_qnlc<-rep(0.00687,length(lower_limit))
  #p_qncPPY<-rep(0.0126,length(lower_limit))   
  
  p_qplc<-p_qpcPPY*0.55
  p_qncPPY<-(p_qpcPPY*0.0126)/0.000786        
  p_qnlc<-p_qncPPY*0.55  
  
  # Chl-a:C
  p_qlcPPY<-rep(0.02,length(lower_limit))
  # Photochemical efficiency
  p_quantum_yield<-rep(0.480e-3,length(lower_limit))

## For each instance (j): substitute old parameter value (i) by new one and write the brick in the fabm.yaml     
    for (j in c(1:length(clasesP))){    
      modulo_modificado<-modulo_productor
      diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
      nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesN[j],sep=""))
      modulo_modificado[inicio:final][3]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][4], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
      modulo_modificado[inicio:final][4]<-nueva  
      
          for (i in c(1:length(parameters))){
            parametro<-get(parameters[i])
            # Locate the parameter
            filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
            fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
            cual<-which(!is.na(filas[,1])==T)
            # Replace parameter value
            nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
            modulo_modificado[inicio:final][cual]<-nueva
          } # end loop i parameters
          
      # Compute initialization values
      if (equal_initialization==T) {carbon<-c_init_phyto
      } else if (separate_phytozoo==T) { carbon<-c_init_phyto/counts_init_phyto[which(clases_init==clasesP[j])]
      } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
      nitro <-carbon*(0.0252/2)
      phosp <-carbon*(0.001572/2)
      chloro<-carbon*(0.04/2)
      # Locate initialization chunk
      location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
      cual<-which(!is.na(location[,1])==T)
      # Replace initialization values
      nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "2.0", as.character(round(carbon,10)))
      modulo_modificado[inicio:final][cual+1]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.0252", as.character(round(nitro,10)))
      modulo_modificado[inicio:final][cual+2]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.001572", as.character(round(phosp,10)))
      modulo_modificado[inicio:final][cual+3]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+4], "0.04", as.character(round(chloro,10)))
      modulo_modificado[inicio:final][cual+4]<-nueva

      # Add to fabm.yaml
      write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
                  sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      
    } # end loop j size_classes  
###################            
  
  
####### P5 ########
old_name<-"P5"
old_long_name<-"coccolithophores"
clasesP<-clasesP5  # ESD = (((2^clasesP)*6)/pi)^(1/3)
clasesN<-str_replace(clasesP,pattern="-","m")
# parameters to modify
parameters<-c("p_sum","p_srs", "p_res","p_qlcPPY",   "p_qup",  "p_qun",   "p_qplc",   "p_qpcPPY",  "p_qnlc",  "p_qncPPY",  "p_quantum_yield")
old_value <-c(  "2.6", "0.09", "0.0", "0.045",       "0.0034", "0.025",   "0.000352", "0.000556",  "0.00687", "0.0126",    "0.447e-3 ")
digito<-c(3,3,2,3,  4,4, 6,6,6,6, 5)
# Extract module block for Px
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)
  #modulo_productor[inicio:final]
  
###################  
lower_limit<-2^clasesP      
upper_limit<-2^(clasesP+1)
  
  # Compute parameter values (in case allometry is taxon-dependent)
  
  # Maximum growth rate
  #p_sum<- (10^0.54)*lower_limit^(-0.15)  # Tang 1995, d-1
  #p_sum<- 2.1*lower_limit^(-0.15)        # b from Tang 1995, a from Irwin 2006
  p_sum<- 1.2*lower_limit^(-0.10)         # Dut Cer 2020

  # Basal metabolism 
  # p_srs<- 10^(-0.063)*lower_limit^0.008 ?? (Ray 2001) d-1
  p_srs<- 0.063-0.008*log10(lower_limit) # Shimoda 2016
  
  # Carbon content per cell
  #carbon<-(10^(-0.29))*lower_limit^0.76    # Mullin 1996  pgC cell-1 
  #carbon<-(10^(-0.665))*lower_limit^0.939  # Menden-Deuer and Lessard, 2000 (overall)
  carbon<-(10^(-0.642))*lower_limit^0.899  # Menden-Deuer and Lessard, 2000 (prymnesiophyceae) # check for coccos
  
  # Minimum P:C quota (Grover 1989), fmol cell-1
  QPmin<-(10^(-1.04))*lower_limit^0.714 
  p_qplc<-(QPmin/carbon)*1e9*1e-12 # mmolP mgC-1
  # Maximal P:C quota (Grover 1989), fmol cell-1, divided by luxury storage to get optimum quata
  QPmax<-(10^(-0.29))*lower_limit^0.767 
  p_qpcPPY<-((QPmax/carbon)*1e9*1e-12)/2 # mmolP mgC-1  
  
  # Minimum N:C quota (Lichtman et al 2007), umol cell-1
  QNmin<-(1.36e-9)*lower_limit^0.77 
  p_qnlc<-(QNmin/carbon)*1e9*1e-3 # mmolN mgC-1
  # Maximal N:C quota (Montagnes & Franklin 2001), umol cell-1, divided by luxury storage to get optimum quota
  QNmax<-(4.64e-9)*lower_limit^0.81 
  p_qncPPY<-((QNmax/carbon)*1e9*1e-3)/2 # mmolN mgC-1 
  
  # Settling velocity, m d−1
  #p_res<-0.029*carbon^0.42      # Moloney and Field, 1991  
  p_res<-0.036*lower_limit^0.43  # Durante et al, 2019 
  
  # P affinity  # Edwards 2012
  aP<-(10^(-8.1))*lower_limit^0.73    # L cell-1 d-1
  p_qup<-(aP/carbon)*1e-3*1e9/4       # m3 mgC-1 d-1
  
  # N affinity # Edwards 2012
  aN<-(10^(-8.2))*lower_limit^0.75    # L cell-1 d-1
  p_qun<-(aN/carbon)*1e-3*1e9         # m3 mgC-1 d-1
  
  ## Original values in BFM
  #p_qplc<-rep(0.000432,length(lower_limit))
  #p_qpcPPY<-rep(0.000786,length(lower_limit))
  #p_qnlc<-rep(0.00687,length(lower_limit))
  #p_qncPPY<-rep(0.0126,length(lower_limit))   
  
  p_qplc<-p_qpcPPY*0.55
  p_qncPPY<-(p_qpcPPY*0.0126)/0.000786        
  p_qnlc<-p_qncPPY*0.55  
  
  # Chl-a:C
  p_qlcPPY<-rep(0.02,length(lower_limit))
  # Photochemical efficiency
  p_quantum_yield<-rep(0.480e-3,length(lower_limit))
  
## For each instance (j): substitute old parameter value (i) by new one and write the brick in the fabm.yaml   
      for (j in c(1:length(clasesP))){    
        modulo_modificado<-modulo_productor
        diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
        nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesN[j],sep=""))
        modulo_modificado[inicio:final][3]<-nueva
        nueva<-str_replace(modulo_modificado[inicio:final][4], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
        modulo_modificado[inicio:final][4]<-nueva  
        
            for (i in c(1:length(parameters))){
              parametro<-get(parameters[i])
              # Locate the parameter
              filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
              fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
              cual<-which(!is.na(filas[,1])==T)
              # Replace parameter value
              nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
              modulo_modificado[inicio:final][cual]<-nueva
            } # end loop i parameters
            
        # Compute initialization values
        if (equal_initialization==T) {carbon<-c_init_phyto
        } else if (separate_phytozoo==T) { carbon<-c_init_phyto/counts_init_phyto[which(clases_init==clasesP[j])]
        } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
        nitro <-carbon*(0.0252/2)
        phosp <-carbon*(0.001572/2)
        chloro<-carbon*(0.04/2)
        # Locate initialization chunk
        location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
        cual<-which(!is.na(location[,1])==T)
        # Replace initialization values
        nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "2.0", as.character(round(carbon,10)))
        modulo_modificado[inicio:final][cual+1]<-nueva
        nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.0252", as.character(round(nitro,10)))
        modulo_modificado[inicio:final][cual+2]<-nueva
        nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.001572", as.character(round(phosp,10)))
        modulo_modificado[inicio:final][cual+3]<-nueva
        nueva<-str_replace(modulo_modificado[inicio:final][cual+4], "0.04", as.character(round(chloro,10)))
        modulo_modificado[inicio:final][cual+4]<-nueva
        
        # Add to fabm.yaml
        write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
                    sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
        
      } # end loop j size_classes  
###################            

  
####### P7 ########
old_name<-"P7"
old_long_name<-"chlorophyceae"
clasesP<-clasesP7  # ESD = (((2^clasesP)*6)/pi)^(1/3)
clasesN<-str_replace(clasesP,pattern="-","m")
# parameters to modify
parameters<-c("p_sum","p_srs", "p_res","p_qlcPPY",   "p_qup",  "p_qun",   "p_qplc",   "p_qpcPPY", "p_qnlc",  "p_qncPPY",  "p_quantum_yield")
old_value <-c(  "2.6", "0.09",   "0.0",  "0.010",    "0.0027", "0.025",   "0.000352", "0.000556", "0.00687", "0.0126",    "0.271e-3")
digito<-c(3,3,2,3,  4,4, 6,6,6,6, 5)
# Extract module block for Px
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)
  #modulo_productor[inicio:final]
  
###################  
lower_limit<-2^clasesP      
upper_limit<-2^(clasesP+1)
  
  # Compute parameter values (in case allometry is taxon-dependent)
  
  # Maximum growth rate
  #p_sum<- (10^0.54)*lower_limit^(-0.15)  # Tang 1995, d-1
  p_sum<- 1.1*lower_limit^(-0.08)         # Dut Cer 2020

  # Basal metabolism 
  # p_srs<- 10^(-0.063)*lower_limit^0.008 ?? (Ray 2001) d-1
  p_srs<- 0.063-0.008*log10(lower_limit) # Shimoda 2016
  
  # Carbon content per cell
  #carbon<-(10^(-0.29))*lower_limit^0.76    # Mullin 1996 pgC cell-1 
  #carbon<-(10^(-0.665))*lower_limit^0.939  # Menden-Deuer and Lessard, 2000 (overall)
  carbon<-(10^(-1.026))*lower_limit^1.088  # Menden-Deuer and Lessard, 2000 (chlorophyceae)
  
  # Minimum P:C quota (Grover 1989), fmol cell-1
  QPmin<-(10^(-1.04))*lower_limit^0.714 
  p_qplc<-(QPmin/carbon)*1e9*1e-12 # mmolP mgC-1
  # Maximal P:C quota (Grover 1989), fmol cell-1, divided by luxury storage to get optimum quota
  QPmax<-(10^(-0.29))*lower_limit^0.767 
  p_qpcPPY<-((QPmax/carbon)*1e9*1e-12)/2 # mmolP mgC-1  
  
  # Minimum N:C quota (Lichtman et al 2007), umol cell-1
  QNmin<-(1.36e-9)*lower_limit^0.77 
  p_qnlc<-(QNmin/carbon)*1e9*1e-3 # mmolN mgC-1
  # Maximal N:C quota (Montagnes & Franklin 2001), umol cell-1, divided by luxury storage to get optimum quota
  QNmax<-(4.64e-9)*lower_limit^0.81 
  p_qncPPY<-((QNmax/carbon)*1e9*1e-3)/2 # mmolN mgC-1   
  
  # Settling velocity, m d−1
  #p_res<-0.029*carbon^0.42      # Moloney and Field, 1991  
  p_res<-0.036*lower_limit^0.43  # Durante et al, 2019 
  
  # P affinity  # Edwards 2012
  aP<-(10^(-8.1))*lower_limit^0.73    # L cell-1 d-1
  p_qup<-(aP/carbon)*1e-3*1e9/4       # m3 mgC-1 d-1
  
  # N affinity # Edwards 2012
  aN<-(10^(-8.2))*lower_limit^0.75    # L cell-1 d-1
  p_qun<-(aN/carbon)*1e-3*1e9         # m3 mgC-1 d-1

  ## Original values in BFM
  #p_qplc<-rep(0.000432,length(lower_limit))
  #p_qpcPPY<-rep(0.000786,length(lower_limit))
  #p_qnlc<-rep(0.00687,length(lower_limit))
  #p_qncPPY<-rep(0.0126,length(lower_limit))   
  
  p_qplc<-p_qpcPPY*0.55
  p_qncPPY<-(p_qpcPPY*0.0126)/0.000786        
  p_qnlc<-p_qncPPY*0.55
  
  # Chl-a:C
  p_qlcPPY<-rep(0.02,length(lower_limit))  
  # Photochemical efficiency
  p_quantum_yield<-rep(0.480e-3,length(lower_limit))
  
## For each instance (j): substitute old parameter value (i) by new one and write the brick in the fabm.yaml   
    for (j in c(1:length(clasesP))){    
      modulo_modificado<-modulo_productor
      diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
      nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesN[j],sep=""))
      modulo_modificado[inicio:final][3]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][4], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
      modulo_modificado[inicio:final][4]<-nueva  
        
        for (i in c(1:length(parameters))){
          parametro<-get(parameters[i])
          # Locate the parameter
          filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
          fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
          cual<-which(!is.na(filas[,1])==T)
          # Replace parameter value
          nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
          modulo_modificado[inicio:final][cual]<-nueva
        } # end loop i parameters
        
      # Compute initialization values
      if (equal_initialization==T) {carbon<-c_init_phyto
      } else if (separate_phytozoo==T) { carbon<-c_init_phyto/counts_init_phyto[which(clases_init==clasesP[j])]
      } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
      nitro <-carbon*(0.0252/2)
      phosp <-carbon*(0.001572/2)
      chloro<-carbon*(0.04/2)
      # Locate initialization chunk
      location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
      cual<-which(!is.na(location[,1])==T)
      # Replace initialization values
      nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "2.0", as.character(round(carbon,10)))
      modulo_modificado[inicio:final][cual+1]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.0252", as.character(round(nitro,10)))
      modulo_modificado[inicio:final][cual+2]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.001572", as.character(round(phosp,10)))
      modulo_modificado[inicio:final][cual+3]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+4], "0.04", as.character(round(chloro,10)))
      modulo_modificado[inicio:final][cual+4]<-nueva
      
      # Add to fabm.yaml
      write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
                  sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      
    } # end loop j size_classes  
###################            

  
####### P8 ########
old_name<-"P8"
old_long_name<-"prasinophyceae"
clasesP<-clasesP8   # ESD = (((2^clasesP)*6)/pi)^(1/3)
clasesN<-str_replace(clasesP,pattern="-","m")
# parameters to modify
parameters<-c("p_sum","p_srs", "p_res","p_qlcPPY",   "p_qup",   "p_qun",   "p_qplc",   "p_qpcPPY", "p_qnlc",  "p_qncPPY",  "p_quantum_yield")
old_value <-c(  "2.6", "0.09",  "0.0",   "0.040",    "0.0038",  "0.025",   "0.000352", "0.000556", "0.00687", "0.0126",    "0.591e-3")
digito<-c(3,3,2,3,  4,4, 6,6,6,6, 5)
# Extract module block for Px
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)
  #modulo_productor[inicio:final]
  
###################  
lower_limit<-2^clasesP      
upper_limit<-2^(clasesP+1)
  
  # Compute parameter values (in case allometry is taxon-dependent)
  
  # Maximum growth rate
  #p_sum<- (10^0.54)*lower_limit^(-0.15)  # Tang 1995 d-1
  #p_sum<- 2.1*lower_limit^(-0.15)        # b from Tang 1995, a from Irwin 2006
  p_sum<- 1.1*lower_limit^(-0.08)         # Dut Cer 2020

  # Basal metabolism 
  # p_srs<- 10^(-0.063)*lower_limit^0.008 ?? (Ray 2001) d-1
  p_srs<- 0.063-0.008*log10(lower_limit) # Shimoda 2016
  
  # Carbon content per cell
  #carbon<-(10^(-0.29))*lower_limit^0.76    # Mullin 1996   pgC cell-1
  #carbon<-(10^(-0.665))*lower_limit^0.939  # Menden-Deuer and Lessard, 2000 (overall)
  carbon<-(10^(-0.545))*lower_limit^0.886  # Menden-Deuer and Lessard, 2000 (prasinophyceae)
  
  # Minimum P:C quota (Grover 1989), fmol cell-1
  QPmin<-(10^(-1.04))*lower_limit^0.714 
  p_qplc<-(QPmin/carbon)*1e9*1e-12 # mmolP mgC-1
  # Maximal P:C quota (Grover 1989), fmol cell-1, divided by luxury storage to get optimum quata
  QPmax<-(10^(-0.29))*lower_limit^0.767 
  p_qpcPPY<-((QPmax/carbon)*1e9*1e-12)/2 # mmolP mgC-1  
  
  # Minimum N:C quota (Lichtman et al 2007), umol cell-1
  QNmin<-(1.36e-9)*lower_limit^0.77 
  p_qnlc<-(QNmin/carbon)*1e9*1e-3 # mmolN mgC-1
  # Maximal N:C quota (Montagnes & Franklin 2001), umol cell-1, divided by luxury storage to get optimum quota
  QNmax<-(4.64e-9)*lower_limit^0.81 
  p_qncPPY<-((QNmax/carbon)*1e9*1e-3)/2 # mmolN mgC-1  
  
  # Settling velocity, m d−1
  #p_res<-0.029*carbon^0.42      # Moloney and Field, 1991  
  p_res<-0.036*lower_limit^0.43  # Durante et al, 2019 
  
  # P affinity  # Edwards 2012
  aP<-(10^(-8.1))*lower_limit^0.73    # L cell-1 d-1
  p_qup<-(aP/carbon)*1e-3*1e9/4       # m3 mgC-1 d-1
  
  # N affinity # Edwards 2012
  aN<-(10^(-8.2))*lower_limit^0.75    # L cell-1 d-1
  p_qun<-(aN/carbon)*1e-3*1e9         # m3 mgC-1 d-1
  
  ## Original values in BFM
  #p_qplc<-rep(0.000432,length(lower_limit))
  #p_qpcPPY<-rep(0.000786,length(lower_limit))
  #p_qnlc<-rep(0.00687,length(lower_limit))
  #p_qncPPY<-rep(0.0126,length(lower_limit))   
  
  p_qplc<-p_qpcPPY*0.55
  p_qncPPY<-(p_qpcPPY*0.0126)/0.000786        
  p_qnlc<-p_qncPPY*0.55
  
  # Chl-a:C
  p_qlcPPY<-rep(0.02,length(lower_limit))
  # Photochemical efficiency
  p_quantum_yield<-rep(0.480e-3,length(lower_limit))
  
## For each instance (j): substitute old parameter value (i) by new one and write the brick in the fabm.yaml   
    for (j in c(1:length(clasesP))){    
      modulo_modificado<-modulo_productor
      diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
      nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesN[j],sep=""))
      modulo_modificado[inicio:final][3]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][4], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
      modulo_modificado[inicio:final][4]<-nueva  
          
          for (i in c(1:length(parameters))){
            parametro<-get(parameters[i])
            # Locate the parameter
            filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
            fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
            cual<-which(!is.na(filas[,1])==T)
            # Replace parameter value
            nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
            modulo_modificado[inicio:final][cual]<-nueva
          } # end loop i parameters

      # Compute initialization values
      if (equal_initialization==T) {carbon<-c_init_phyto
      } else if (separate_phytozoo==T) { carbon<-c_init_phyto/counts_init_phyto[which(clases_init==clasesP[j])]
      } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
      nitro <-carbon*(0.0252/2)
      phosp <-carbon*(0.001572/2)
      chloro<-carbon*(0.04/2)
      # Locate initialization chunk
      location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
      cual<-which(!is.na(location[,1])==T)
      # Replace initialization values
      nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "2.0", as.character(round(carbon,10)))
      modulo_modificado[inicio:final][cual+1]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.0252", as.character(round(nitro,10)))
      modulo_modificado[inicio:final][cual+2]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.001572", as.character(round(phosp,10)))
      modulo_modificado[inicio:final][cual+3]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+4], "0.04", as.character(round(chloro,10)))
      modulo_modificado[inicio:final][cual+4]<-nueva
      
      # Add to fabm.yaml
      write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
                  sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      
    } # end loop j size_classes  
###################            
  


####### P3 ########
old_name<-"P3"
old_long_name<-"smalleukaryotes"
clasesP<-clasesP3   # ESD = (((2^clasesP)*6)/pi)^(1/3)
clasesN<-str_replace(clasesP,pattern="-","m")
# parameters to modify
parameters<-c("p_sum","p_srs", "p_res","p_qlcPPY",   "p_qup",  "p_qun",   "p_qplc",    "p_qpcPPY", "p_qnlc",  "p_qncPPY",  "p_quantum_yield")
old_value <-c(  "3.5", "0.1", "0.0", "0.010",       "0.0036",  "0.25",    "0.0004288", "0.00072",  "0.00687", "0.0126",    "1.043e-3")
digito<-c(3,3,2,3,  4,4, 6,6,6,6, 5)
# Extract module block for Px
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)
  #modulo_productor[inicio:final]
  
###################  
lower_limit<-2^clasesP      
upper_limit<-2^(clasesP+1)
  
  # Compute parameter values (in case allometry is taxon-dependent)
  
  # Maximum growth rate
  #p_sum<- (10^0.54)*lower_limit^(-0.15)  # Tang 1995, d-1
  #p_sum<- 2.1*lower_limit^(-0.15)        # b from Tang 1995, a from Irwin 2006
  p_sum<- 1.0*lower_limit^(-0.02)        # Dut Cer 2020

  # Basal metabolism (Ray 2001) d-1
  # p_srs<- 10^(-0.063)*lower_limit^0.008 ??
  p_srs<- 0.063-0.008*log10(lower_limit) # Shimoda 2016
  
  # Carbon content per cell
  #carbon<-(10^(-0.29))*lower_limit^0.76    # Mullin 1996 pgC cell-1 
  carbon<-(10^(-0.665))*lower_limit^0.939  # Menden-Deuer and Lessard, 2000 (overall)
  #carbon<-(10^(-1.694))*lower_limit^1.218  # Menden-Deuer and Lessard, 2000 (Chrysophytes)
  
  # Minimum P:C quota (Grover 1989), fmol cell-1
  QPmin<-(10^(-1.04))*lower_limit^0.714 
  p_qplc<-(QPmin/carbon)*1e9*1e-12 # mmolP mgC-1
  # Maximal P:C quota (Grover 1989), fmol cell-1, divided by luxury storage to get optimum quata
  QPmax<-(10^(-0.29))*lower_limit^0.767 
  p_qpcPPY<-((QPmax/carbon)*1e9*1e-12)/2 # mmolP mgC-1  
  
  # Minimum N:C quota (Lichtman et al 2007), umol cell-1
  QNmin<-(1.36e-9)*lower_limit^0.77 
  p_qnlc<-(QNmin/carbon)*1e9*1e-3 # mmolN mgC-1
  # Maximal N:C quota (Montagnes & Franklin 2001), umol cell-1, divided by luxury storage to get optimum quota
  QNmax<-(4.64e-9)*lower_limit^0.81 
  p_qncPPY<-((QNmax/carbon)*1e9*1e-3)/2 # mmolN mgC-1  
  
  # Settling velocity, m d−1
  #p_res<-0.029*carbon^0.42      # Moloney and Field, 1991  
  p_res<-0.036*lower_limit^0.43  # Durante et al, 2019 
  
  # P affinity  # Edwards 2012
  aP<-(10^(-8.1))*lower_limit^0.73      # L cell-1 d-1
  p_qup<-(aP/carbon)*1e-3*1e9/3         # m3 mgC-1 d-1
  
  # N affinity # Edwards 2012
  aN<-(10^(-8.2))*lower_limit^0.75      # L cell-1 d-1
  p_qun<-(aN/carbon)*1e-3*1e9*(4/3)     # m3 mgC-1 d-1
  
  ## Original values in BFM
  #p_qplc<-rep(0.000432,length(lower_limit))
  #p_qpcPPY<-rep(0.000786,length(lower_limit))
  #p_qnlc<-rep(0.00687,length(lower_limit))
  #p_qncPPY<-rep(0.0126,length(lower_limit))   
  
  p_qplc<-p_qpcPPY*0.55
  p_qncPPY<-(p_qpcPPY*0.0126)/0.000786        
  p_qnlc<-p_qncPPY*0.55
  
  # Chl-a:C
  p_qlcPPY<-rep(0.02,length(lower_limit))
  # Photochemical efficiency
  p_quantum_yield<-rep(0.480e-3,length(lower_limit))
  
## For each instance (j): substitute old parameter value (i) by new one and write the brick in the fabm.yaml   
    for (j in c(1:length(clasesP))){    
      modulo_modificado<-modulo_productor
      diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
      nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesN[j],sep=""))
      modulo_modificado[inicio:final][3]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][4], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
      modulo_modificado[inicio:final][4]<-nueva  
      
        for (i in c(1:length(parameters))){
          parametro<-get(parameters[i])
          # Locate the parameter
          filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
          fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
          cual<-which(!is.na(filas[,1])==T)
          # Replace parameter value
          nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
          modulo_modificado[inicio:final][cual]<-nueva
        } # end loop i parameters
        
      # Compute initialization values
      if (equal_initialization==T) {carbon<-c_init_phyto
      } else if (separate_phytozoo==T) { carbon<-c_init_phyto/counts_init_phyto[which(clases_init==clasesP[j])]
      } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
      nitro <-carbon*(0.0336/2.67)
      phosp <-carbon*(0.002096/2.67)
      chloro<-carbon*(0.053/2.67)
      # Locate initialization chunk
      location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
      cual<-which(!is.na(location[,1])==T)
      # Replace initialization values
      nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "2.67", as.character(round(carbon,10)))
      modulo_modificado[inicio:final][cual+1]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.0336", as.character(round(nitro,10)))
      modulo_modificado[inicio:final][cual+2]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.002096", as.character(round(phosp,10)))
      modulo_modificado[inicio:final][cual+3]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+4], "0.053", as.character(round(chloro,10)))
      modulo_modificado[inicio:final][cual+4]<-nueva

      # Add to fabm.yaml
      write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
                  sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      
    } # end loop j size_classes  
###################            

  
####### P9 ########
old_name<-"P9"
old_long_name<-"synechococcus"
clasesP<-clasesP9   # ESD = (((2^clasesP)*6)/pi)^(1/3)
clasesN<-str_replace(clasesP,pattern="-","m")
# parameters to modify
parameters<-c("p_sum","p_srs", "p_res","p_qlcPPY",   "p_qup",  "p_qun",   "p_qplc",  "p_qpcPPY", "p_qnlc",  "p_qncPPY", "p_quantum_yield")
old_value <-c(  "3.5", "0.1",  "0.0",   "0.014",    "0.0034",  "0.25", "0.0004288", "0.00072",  "0.00687", "0.0126",    "1.570e-3")
digito<-c(3,3,2,3,  4,4, 6,6,6,6, 5)
# Extract module block for Px
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)
  #modulo_productor[inicio:final]
  
###################  
lower_limit<-2^clasesP      
upper_limit<-2^(clasesP+1)
  
  # Compute parameter values (in case allometry is taxon-dependent)
  
  # Maximum growth rate
  #p_sum<- (10^0.54)*lower_limit^(-0.15)   # Tang 1995 d-1
  #p_sum<- 1.4*lower_limit^(-0.15)         # b from Tang 1995, a from Irwin 2006
  #p_sum<- 2.1*lower_limit^(-0.15)         # b from Tang 1995, a equal to other eukaryotes
  p_sum<- 1.0*lower_limit^(+0.06)          # Dut Cer 2020
      
  # Basal metabolism
  # p_srs<- 10^(-0.063)*lower_limit^0.008 ?? (Ray 2001) d-1
  p_srs<- 0.063-0.008*log10(lower_limit) # Shimoda 2016
  
  # Carbon content per cell
  #carbon<-(10^(-0.29))*lower_limit^0.76    # Mullin 1996 pgC cell-1
  #carbon<-(10^(-0.665))*lower_limit^0.939  # Menden-Deuer and Lessard, 2000 (overall)
  carbon<-(10^(-0.583))*lower_limit^0.860  # Menden-Deuer and Lessard, 2000 (protists<3000)
  
  # Minimum P:C quota (Grover 1989), fmol cell-1
  QPmin<-(10^(-1.04))*lower_limit^0.714 
  p_qplc<-(QPmin/carbon)*1e9*1e-12 # mmolP mgC-1
  # Maximal P:C quota (Grover 1989), fmol cell-1, divided by luxury storage to get optimum quata
  QPmax<-(10^(-0.29))*lower_limit^0.767 
  p_qpcPPY<-((QPmax/carbon)*1e9*1e-12)/2 # mmolP mgC-1  
  
  # Minimum N:C quota (Lichtman et al 2007), umol cell-1
  QNmin<-(1.36e-9)*lower_limit^0.77 
  p_qnlc<-(QNmin/carbon)*1e9*1e-3 # mmolN mgC-1
  # Maximal N:C quota (Montagnes & Franklin 2001), umol cell-1, divided by luxury storage to get optimum quota
  QNmax<-(4.64e-9)*lower_limit^0.81 
  p_qncPPY<-((QNmax/carbon)*1e9*1e-3)/2 # mmolN mgC-1  
  
  # Minimum N:C quota (Lichtman et al 2007), umol cell-1
  QNmin<-(1.36e-9)*lower_limit^0.77 
  p_qnlc<-(QNmin/carbon)*1e9*1e-3 # mmolN mgC-1
  # Maximal N:C quota (Montagnes & Franklin 2001), umol cell-1, divided by luxury storage to get optimum quota
  QNmax<-(4.64e-9)*lower_limit^0.81 
  p_qncPPY<-((QNmax/carbon)*1e9*1e-3)/2 # mmolN mgC-1  
  
  # Settling velocity, m d−1
  #p_res<-0.029*carbon^0.42      # Moloney and Field, 1991  
  p_res<-0.036*lower_limit^0.43  # Durante et al, 2019 
  
  # P affinity  # Edwards 2012
  aP<-(10^(-8.1))*lower_limit^0.73    # L cell-1 d-1
  p_qup<-(aP/carbon)*1e-3*1e9/2       # m3 mgC-1 d-1
  
  # N affinity # Edwards 2012
  aN<-(10^(-8.2))*lower_limit^0.75    # L cell-1 d-1
  p_qun<-(aN/carbon)*1e-3*1e9*(4/2)   # m3 mgC-1 d-1
  
  ## Original values in BFM
  #p_qplc<-rep(0.000432,length(lower_limit))
  #p_qpcPPY<-rep(0.000786,length(lower_limit))
  #p_qnlc<-rep(0.00687,length(lower_limit))
  #p_qncPPY<-rep(0.0126,length(lower_limit))   
  
  p_qplc<-p_qpcPPY*0.55
  p_qncPPY<-(p_qpcPPY*0.0126)/0.000786        
  p_qnlc<-p_qncPPY*0.55
  
  # Chl-a:C
  p_qlcPPY<-rep(0.02,length(lower_limit))
  # Photochemical efficiency
  p_quantum_yield<-rep(0.480e-3,length(lower_limit))
  
## For each instance (j): substitute old parameter value (i) by new one and write the brick in the fabm.yaml   
    for (j in c(1:length(clasesP))){    
      modulo_modificado<-modulo_productor
      diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
      nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesN[j],sep=""))
      modulo_modificado[inicio:final][3]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][4], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
      modulo_modificado[inicio:final][4]<-nueva  
      
        for (i in c(1:length(parameters))){
          parametro<-get(parameters[i])
          # Locate the parameter
          filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
          fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
          cual<-which(!is.na(filas[,1])==T)
          # Replace parameter value
          nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
          modulo_modificado[inicio:final][cual]<-nueva
        } # end loop i parameters
        
      # Compute initialization values
      if (equal_initialization==T) {carbon<-c_init_phyto
      } else if (separate_phytozoo==T) { carbon<-c_init_phyto/counts_init_phyto[which(clases_init==clasesP[j])]
      } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
      nitro <-carbon*(0.0336/2.67)
      phosp <-carbon*(0.002096/2.67)
      chloro<-carbon*(0.053/2.67)
      # Locate initialization chunk
      location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
      cual<-which(!is.na(location[,1])==T)
      # Replace initialization values
      nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "2.67", as.character(round(carbon,10)))
      modulo_modificado[inicio:final][cual+1]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.0336", as.character(round(nitro,10)))
      modulo_modificado[inicio:final][cual+2]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.002096", as.character(round(phosp,10)))
      modulo_modificado[inicio:final][cual+3]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+4], "0.053", as.character(round(chloro,10)))
      modulo_modificado[inicio:final][cual+4]<-nueva
      
      # Add to fabm.yaml
      write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
                  sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      
    } # end loop j size_classes  
###################            

  
####### P6 ########
old_name<-"P6"
old_long_name<-"prochlorococcus"
clasesP<-clasesP6   # ESD = (((2^clasesP)*6)/pi)^(1/3)
clasesN<-str_replace(clasesP,pattern="-","m")
# parameters to modify
parameters<-c("p_sum","p_srs", "p_res","p_qlcPPY",   "p_qup",  "p_qun",   "p_qplc",    "p_qpcPPY",  "p_qnlc",  "p_qncPPY",  "p_quantum_yield")
old_value <-c(  "3.5", "0.1", "0.0", "0.010",       "0.0027",   "0.25",   "0.0004288", "0.00072",  "0.00687",  "0.0126",    "1.489e-3")
digito<-c(3,3,2,3,  4,4, 6,6,6,6, 5)
# Extract module block for Px
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)
  #modulo_productor[inicio:final]
  
###################  
lower_limit<-2^clasesP      
upper_limit<-2^(clasesP+1)
  
  # Compute parameter values (in case allometry is taxon-dependent)
  
  # Maximum growth rate 
  #p_sum<- (10^0.54)*lower_limit^(-0.15)  # (Tang 1995) d-1
  #p_sum<- 1.0*lower_limit^(-0.15)        # b from Tang 1995, a from Irwin 2006
  #p_sum<- 2.1*lower_limit^(-0.15)        # b from Tang 1995, a equal to other eukaryotes
  p_sum<- 1.05*lower_limit^(+0.08)        # Dut Cer 2020
  
  # Basal metabolism (Ray 2001) d-1
  # p_srs<- 10^(-0.063)*lower_limit^0.008 ??
  p_srs<- 0.063-0.008*log10(lower_limit) # Shimoda 2016
  
  # Carbon content per cell
  #carbon<-(10^(-0.29))*lower_limit^0.76    # Mullin 1996 pgC cell-1 
  #carbon<-(10^(-0.665))*lower_limit^0.939  # Menden-Deuer and Lessard, 2000 (overall)
  carbon<-(10^(-0.583))*lower_limit^0.860  # Menden-Deuer and Lessard, 2000 (protists<3000)
  
  # Minimum P:C quota (Grover 1989), fmol cell-1
  QPmin<-(10^(-1.04))*lower_limit^0.714 
  p_qplc<-(QPmin/carbon)*1e9*1e-12 # mmolP mgC-1
  # Maximal P:C quota (Grover 1989), fmol cell-1, divided by luxury storage to get optimum quata
  QPmax<-(10^(-0.29))*lower_limit^0.767 
  p_qpcPPY<-((QPmax/carbon)*1e9*1e-12)/2 # mmolP mgC-1  
  
  # Minimum N:C quota (Lichtman et al 2007), umol cell-1
  QNmin<-(1.36e-9)*lower_limit^0.77 
  p_qnlc<-(QNmin/carbon)*1e9*1e-3 # mmolN mgC-1
  # Maximal N:C quota (Montagnes & Franklin 2001), umol cell-1, divided by luxury storage to get optimum quota
  QNmax<-(4.64e-9)*lower_limit^0.81 
  p_qncPPY<-((QNmax/carbon)*1e9*1e-3)/2 # mmolN mgC-1  
  
  # Settling velocity, m d−1
  #p_res<-0.029*carbon^0.42      # Moloney and Field, 1991  
  p_res<-0.036*lower_limit^0.43  # Durante et al, 2019 
  
  # P affinity  # Edwards 2012
  aP<-(10^(-8.1))*lower_limit^0.73    # L cell-1 d-1
  p_qup<-(aP/carbon)*1e-3*1e9/2       # m3 mgC-1 d-1
  
  # N affinity # Edwards 2012
  aN<-(10^(-8.2))*lower_limit^0.75    # L cell-1 d-1
  p_qun<-(aN/carbon)*1e-3*1e9*(4/2)   # m3 mgC-1 d-1
  
  ## Original values in BFM 
  #p_qplc<-rep(0.000432,length(lower_limit))
  #p_qpcPPY<-rep(0.000786,length(lower_limit))
  #p_qnlc<-rep(0.00687,length(lower_limit))
  #p_qncPPY<-rep(0.0126,length(lower_limit))   
  
  p_qplc<-p_qpcPPY*0.55
  p_qncPPY<-(p_qpcPPY*0.0126)/0.000786        
  p_qnlc<-p_qncPPY*0.55
  
  # Chl-a:C
  p_qlcPPY<-rep(0.02,length(lower_limit))
  # Photochemical efficiency
  p_quantum_yield<-rep(0.480e-3,length(lower_limit))
  
## For each instance (j): substitute old parameter value (i) by new one and write the brick in the fabm.yaml   
    for (j in c(1:length(clasesP))){    
      modulo_modificado<-modulo_productor
      diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
      nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesN[j],sep=""))
      modulo_modificado[inicio:final][3]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][4], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
      modulo_modificado[inicio:final][4]<-nueva  
      
          for (i in c(1:length(parameters))){
            parametro<-get(parameters[i])
            # Locate the parameter
            filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
            fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
            cual<-which(!is.na(filas[,1])==T)
            # Replace parameter value
            nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
            modulo_modificado[inicio:final][cual]<-nueva
          } # end loop i parameters

      # Compute initialization values
      if (equal_initialization==T) {carbon<-c_init_phyto
      } else if (separate_phytozoo==T) { carbon<-c_init_phyto/counts_init_phyto[which(clases_init==clasesP[j])]
      } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
      nitro <-carbon*(0.0336/2.67)
      phosp <-carbon*(0.002096/2.67)
      chloro<-carbon*(0.053/2.67)
      # Locate initialization chunk
      location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
      cual<-which(!is.na(location[,1])==T)
      # Replace initialization values
      nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "2.67", as.character(round(carbon,10)))
      modulo_modificado[inicio:final][cual+1]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.0336", as.character(round(nitro,10)))
      modulo_modificado[inicio:final][cual+2]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.002096", as.character(round(phosp,10)))
      modulo_modificado[inicio:final][cual+3]<-nueva
      nueva<-str_replace(modulo_modificado[inicio:final][cual+4], "0.053", as.character(round(chloro,10)))
      modulo_modificado[inicio:final][cual+4]<-nueva
      
      # Add to fabm.yaml
      write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
                  sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      
    } # end loop j size_classes  
###################            
  
  
  
  
### Create BRICKS for PREDATORS
  
####### Z6 ########
old_name<-"Z6"
old_long_name<-"Heterotrophic Nanoflagellates (HNAN)"
clasesP<-clasesZ6  # ESD = (((2^clasesP)*6)/pi)^(1/3)
# parameters to modify
parameters<-c("p_sum", "p_chuc")
old_value <-c("3.88",   "90")
digito<-c(2,0)

  
old_nprey=12
# Which are the potential preys? define here all prey choices non size-dependent, if it is not here, it is not eaten
preys<-c(clasesB1, clasesP2, clasesP3, clasesP5, clasesP6, clasesP7, clasesP8, clasesP9, clasesZ6)
names_preys<-c("B1", paste("P2",str_replace(clasesP2,pattern="-","m"), sep="_"),
               paste("P3",str_replace(clasesP3,pattern="-","m"), sep="_"),
               paste("P5",str_replace(clasesP5,pattern="-","m"), sep="_"),
               paste("P6",str_replace(clasesP6,pattern="-","m"), sep="_"),
               paste("P7",str_replace(clasesP7,pattern="-","m"), sep="_"),
               paste("P8",str_replace(clasesP8,pattern="-","m"), sep="_"),
               paste("P9",str_replace(clasesP9,pattern="-","m"), sep="_"),
               paste("Z6",str_replace(clasesZ6,pattern="-","m"), sep="_"))

################### 
# Extract module block for Zx
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)
   #modulo_productor[inicio:final]

lower_limit<-2^clasesP      
upper_limit<-2^(clasesP+1)
diameter<-((lower_limit*6)/pi)^(1/3)

## Compute parameter values (in case allometry is taxon-dependent)

# Maximum ingestion rate (Hansen 1997) d-1
p_sum<- 26*diameter^(-0.4)
# Half-saturation constant (Banas 2011) 3 uM N, not used finally
#p_chuc<-rep(3*(106/16)*12,length=length(diameter))
p_chuc<-rep(120,length=length(diameter))
#p_vum<-p_sum/p_chuc

## Find optimum prey size
diameter_preys<-(((2^preys)*6)/pi)^(1/3) # range(diameter_preys) 
#opt_prey_size<-0.65*diameter^0.56   # Hansen 1997, Banas 2011
opt_prey_size<-diameter/10           # Kiorboe 2008, Ward 2012
## Carbon conversion from Alcaraz et al 2003 (zooplankton<5mm WM) not used for now
#carbon<-0.0825+0.0780*lower_limit

## For each instance (j): substitute old parameter value (i) by new one,
##                        compute preferences and fill prey blocks,
##                        write the Z brick in the fabm.yaml

  for (j in c(1:length(clasesP))){    
    modulo_modificado<-modulo_productor
    diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
    # Replace module names
    nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesP[j],sep=""))
    modulo_modificado[inicio:final][3]<-nueva
    nueva<-str_replace(modulo_modificado[inicio:final][5], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
    modulo_modificado[inicio:final][5]<-nueva
    
      # Replace parameter values    
      for (i in c(1:length(parameters))){
        parametro<-get(parameters[i])
        # Locate the parameter
        filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
        fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
        cual<-which(!is.na(filas[,1])==T)
        # Replace parameter value
        nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
        modulo_modificado[inicio:final][cual]<-nueva
      } # end loop i parameters
      
    # Compute preferences, sd of the gaussian = 0.20
    prefs<-exp(-((log10(diameter_preys)-log10(opt_prey_size[j]))/0.20)^2)
    prefs[prefs<1e-2]<-0
    o<-order(2^preys)

    # Replace number of preys (nprey)
    filita<-str_locate(modulo_modificado[inicio:final], pattern="nprey")
    fila<-modulo_modificado[inicio:final][!is.na(filita[,1])]
    cual<-which(!is.na(filita[,1])==T)
    nueva<-str_replace(fila, as.character(old_nprey), as.character(length(names_preys[prefs!=0])))    
    modulo_modificado[inicio:final][cual]<-nueva
    
    ## For the next three blocks (preferences, calcifier and coupling) the template needs to have
    ## enough lines available to accommodate all preys, otherwise it does not complain but later fails at runtime.
    
    # Block of preferences
    new_lines<- paste("      suprey", c(1:length(names_preys[prefs!=0])),":   ",
                      round(prefs[prefs!=0],2),"     # [-]           Availability of ",
                      names_preys[prefs!=0]," to Z6_",clasesP[j], sep="")
    modulo_modificado[inicio:final][(cual+1):(cual+length(names_preys[prefs!=0]))]<-new_lines
    
    # Block of calcifier (1 if the prey is calcifier, 0 otherwise)
    boolean<-rep(0,length(names_preys[prefs!=0]))
    boolean[(substring(names_preys[prefs!=0],1,2)=="P5")]<-1
    new_lines<- paste("      isP2", c(1:length(names_preys[prefs!=0])),":   ",
                      boolean,"           # [-] identify P2 among the preys [1 for P2 otherwise 0]", sep="")
    filita<-str_locate(modulo_modificado[inicio:final], pattern="isP21:")
    fila<-modulo_modificado[inicio:final][!is.na(filita[,1])]
    cual<-which(!is.na(filita[,1])==T)    
    modulo_modificado[inicio:final][(cual):(cual-1+length(names_preys[prefs!=0]))]<-new_lines    
    
    # Block of coupling
    new_lines<- paste("      prey", c(1:length(names_preys[prefs!=0])),":   ",names_preys[prefs!=0], sep="")
    filita<-str_locate(modulo_modificado[inicio:final], pattern="      prey1:")
    fila<-modulo_modificado[inicio:final][!is.na(filita[,1])]
    cual<-which(!is.na(filita[,1])==T)    
    modulo_modificado[inicio:final][(cual):(cual-1+length(names_preys[prefs!=0]))]<-new_lines    

    # Compute initialization values
    if (equal_initialization==T) {carbon<-c_init_zoo
    } else if (separate_phytozoo==T) { carbon<-c_init_zoo/counts_init_zoo[which(clases_init==clasesP[j])]
    } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
    nitro <-carbon*(0.0167/1.0)
    phosp <-carbon*(0.00185/1.0)
    # Locate initialization chunk
    location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
    cual<-which(!is.na(location[,1])==T)
    # Replace initialization values
    nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "1.0", as.character(round(carbon,10)))
    modulo_modificado[inicio:final][cual+1]<-nueva
    nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.0167", as.character(round(nitro,10)))
    modulo_modificado[inicio:final][cual+2]<-nueva
    nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.00185", as.character(round(phosp,10)))
    modulo_modificado[inicio:final][cual+3]<-nueva

    # Add to fabm.yaml
    write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
                sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
    
  } # end loop j size_classes  
################### 



  
####### Z5 ########
old_name<-"Z5"
old_long_name<-"Microzooplankton"
clasesP<-clasesZ5  # ESD = (((2^clasesP)*6)/pi)^(1/3)
# parameters to modify
parameters<-c("p_sum","p_chuc")
old_value <-c( "2.71",   "20")
digito<-c(2,0)

old_nprey=12
# Which are the potential preys? define here all prey choices non size-dependent, if it is not here, it is not eaten
preys<-c(clasesB1,clasesP1,clasesP2,clasesP3,clasesP4,clasesP5,clasesP6,clasesP7,clasesP8,clasesP9,clasesZ5,clasesZ6)
names_preys<-c("B1", paste("P1",str_replace(clasesP1,pattern="-","m"), sep="_"),
               paste("P2",str_replace(clasesP2,pattern="-","m"), sep="_"),
               paste("P3",str_replace(clasesP3,pattern="-","m"), sep="_"),
               paste("P4",str_replace(clasesP4,pattern="-","m"), sep="_"),
               paste("P5",str_replace(clasesP5,pattern="-","m"), sep="_"),
               paste("P6",str_replace(clasesP6,pattern="-","m"), sep="_"),
               paste("P7",str_replace(clasesP7,pattern="-","m"), sep="_"),
               paste("P8",str_replace(clasesP8,pattern="-","m"), sep="_"),
               paste("P9",str_replace(clasesP9,pattern="-","m"), sep="_"),
               paste("Z5",str_replace(clasesZ5,pattern="-","m"), sep="_"),
               paste("Z6",str_replace(clasesZ6,pattern="-","m"), sep="_"))

################### 
# Extract module block for Zx
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)
#modulo_productor[inicio:final]

lower_limit<-2^clasesP      
upper_limit<-2^(clasesP+1)
diameter<-((lower_limit*6)/pi)^(1/3)

# Compute parameter values (in case allometry is taxon-dependent)
# Maximum ingestion rate (Hansen 1997) d-1
p_sum<- 26*diameter^(-0.4)
# Half-saturation constant (Banas 2011) 3 uM N, not used finally
#p_chuc<-rep(3*(106/16)*12,length=length(diameter))
p_chuc<-rep(30,length=length(diameter))
#p_vum<-p_sum/p_chuc

# Find optimum prey size
diameter_preys<-(((2^preys)*6)/pi)^(1/3) # range(diameter_preys) 
#opt_prey_size<-0.65*diameter^0.56   # Hansen 1997, Banas 2011
opt_prey_size<-diameter/10           # Kiorboe 2008, Ward 2012
# Carbon conversion from Alcaraz et al 2003 (zooplankton<5mm WM) not used for now
#carbon<-0.0825+0.0780*lower_limit

## For each instance (j): substitute old parameter value (i) by new one,
##                        compute preferences and fill prey blocks,
##                        write the Z brick in the fabm.yaml

for (j in c(1:length(clasesP))){    
  modulo_modificado<-modulo_productor
  diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
  # Replace module names
  nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesP[j],sep=""))
  modulo_modificado[inicio:final][3]<-nueva
  nueva<-str_replace(modulo_modificado[inicio:final][5], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
  modulo_modificado[inicio:final][5]<-nueva  
  
  # Replace parameter values    
  for (i in c(1:length(parameters))){
    parametro<-get(parameters[i])
    # Locate the parameter
    filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
    fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
    cual<-which(!is.na(filas[,1])==T)
    # Replace parameter value
    nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
    modulo_modificado[inicio:final][cual]<-nueva
  } # end loop i parameters
  
  # Compute preferences, sd of the gaussian = 0.20
  prefs<-exp(-((log10(diameter_preys)-log10(opt_prey_size[j]))/0.20)^2)
  prefs[prefs<1e-2]<-0
  o<-order(2^preys)

  # Replace number of preys
  filita<-str_locate(modulo_modificado[inicio:final], pattern="nprey")
  fila<-modulo_modificado[inicio:final][!is.na(filita[,1])]
  cual<-which(!is.na(filita[,1])==T)
  nueva<-str_replace(fila, as.character(old_nprey), as.character(length(names_preys[prefs!=0])))    
  modulo_modificado[inicio:final][cual]<-nueva
  
  ## For the next three blocks (preferences, calcifier and coupling) the template needs to have
  ## enough lines available to accommodate all preys, otherwise it does not complain but later fails at runtime.
  
  # Block of preferences
  new_lines<- paste("      suprey", c(1:length(names_preys[prefs!=0])),":   ",
                    round(prefs[prefs!=0],2),"     # [-]           Availability of ",
                    names_preys[prefs!=0]," to Z5_",clasesP[j], sep="")
  modulo_modificado[inicio:final][(cual+1):(cual+length(names_preys[prefs!=0]))]<-new_lines
  
  # Block of calcifier (1 if the prey is calcifier, 0 otherwise)
  boolean<-rep(0,length(names_preys[prefs!=0]))
  boolean[(substring(names_preys[prefs!=0],1,2)=="P5")]<-1
  new_lines<- paste("      isP2", c(1:length(names_preys[prefs!=0])),":   ",
                    boolean,"           # [-] identify P2 among the preys [1 for P2 otherwise 0]", sep="")
  filita<-str_locate(modulo_modificado[inicio:final], pattern="isP21:")
  fila<-modulo_modificado[inicio:final][!is.na(filita[,1])]
  cual<-which(!is.na(filita[,1])==T)    
  modulo_modificado[inicio:final][(cual):(cual-1+length(names_preys[prefs!=0]))]<-new_lines    
  
  # Block of coupling
  new_lines<- paste("      prey", c(1:length(names_preys[prefs!=0])),":   ",names_preys[prefs!=0], sep="")
  filita<-str_locate(modulo_modificado[inicio:final], pattern="      prey1:")
  fila<-modulo_modificado[inicio:final][!is.na(filita[,1])]
  cual<-which(!is.na(filita[,1])==T)    
  modulo_modificado[inicio:final][(cual):(cual-1+length(names_preys[prefs!=0]))]<-new_lines    
  
  # Compute initialization values
  if (equal_initialization==T) {carbon<-c_init_zoo
  } else if (separate_phytozoo==T) { carbon<-c_init_zoo/counts_init_zoo[which(clases_init==clasesP[j])]
  } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
  nitro <-carbon*(0.0167/1.0)
  phosp <-carbon*(0.00185/1.0)
  # Locate initialization chunk
  location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
  cual<-which(!is.na(location[,1])==T)
  # Replace initialization values
  nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "1.0", as.character(round(carbon,10)))
  modulo_modificado[inicio:final][cual+1]<-nueva
  nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.0167", as.character(round(nitro,10)))
  modulo_modificado[inicio:final][cual+2]<-nueva
  nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.00185", as.character(round(phosp,10)))
  modulo_modificado[inicio:final][cual+3]<-nueva
  
  # Add to fabm.yaml
  write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
              sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
  
} # end loop j size_classes  
################### 



####### Z4 ########
old_name<-"Z4"
old_long_name<-"Omnivorous Mesozooplankton"
clasesP<-clasesZ4  # ESD = (((2^clasesP)*6)/pi)^(1/3)
# parameters to modify
parameters<-c("p_sum") #,"p_vum")
old_value <-c( "2.0") #, "0.02")
digito<-c(2) #,4)

old_nprey=8
# Which are the potential preys? define here all prey choices non size-dependent, if it is not here, it is not eaten
preys<-c(clasesP1,clasesP2,clasesP4,clasesP5,clasesP7,clasesP8,clasesZ6,clasesZ5,clasesZ4)
names_preys<-c(paste("P1",str_replace(clasesP1,pattern="-","m"), sep="_"),
               paste("P2",str_replace(clasesP2,pattern="-","m"), sep="_"),
               paste("P4",str_replace(clasesP4,pattern="-","m"), sep="_"),
               paste("P5",str_replace(clasesP5,pattern="-","m"), sep="_"),
               paste("P7",str_replace(clasesP7,pattern="-","m"), sep="_"),
               paste("P8",str_replace(clasesP8,pattern="-","m"), sep="_"),
               paste("Z6",str_replace(clasesZ6,pattern="-","m"), sep="_"),
               paste("Z5",str_replace(clasesZ5,pattern="-","m"), sep="_"),
               paste("Z4",str_replace(clasesZ4,pattern="-","m"), sep="_"))

################### 
# Extract module block for Px
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)
#modulo_productor[inicio:final]

lower_limit<-2^clasesP      
upper_limit<-2^(clasesP+1)
diameter<-((lower_limit*6)/pi)^(1/3)

# Compute parameter values (in case allometry is taxon-dependent)
# Maximum ingestion rate (Hansen 1997) d-1
p_sum<- 26*diameter^(-0.4)
# Half-saturation constant (Banas 2011) 3 uM N
p_chuc<-rep(3*(106/16)*12,length=length(diameter)) # not used
p_vum<-p_sum/p_chuc                                # not used

# Find optimum prey size
diameter_preys<-(((2^preys)*6)/pi)^(1/3) # range(diameter_preys) 
#opt_prey_size<-0.65*diameter^0.56   # Hansen 1997, Banas 2011
opt_prey_size<-diameter/10           # Kiorboe 2008, Ward 2012
## Carbon conversion from Alcaraz et al 2003 (zooplankton<5mm WM) not used for now
#carbon<-0.0825+0.0780*lower_limit

## For each instance (j): substitute old parameter value (i) by new one,
##                        compute preferences and fill prey blocks,
##                        write the Z brick in the fabm.yaml

for (j in c(1:length(clasesP))){    
  modulo_modificado<-modulo_productor
  diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
  # Replace module names
  nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesP[j],sep=""))
  modulo_modificado[inicio:final][3]<-nueva
  nueva<-str_replace(modulo_modificado[inicio:final][5], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
  modulo_modificado[inicio:final][5]<-nueva
  
  # Replace parameter values    
  for (i in c(1:length(parameters))){
    parametro<-get(parameters[i])
    # Locate the parameter
    filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
    fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
    cual<-which(!is.na(filas[,1])==T)
    # Replace parameter value
    nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
    modulo_modificado[inicio:final][cual]<-nueva
  } # end loop i parameters
  
  # Compute preferences, sd of the gaussian = 0.2
  prefs<-exp(-((log10(diameter_preys)-log10(opt_prey_size[j]))/0.20)^2)
  prefs[prefs<1e-2]<-0
  o<-order(2^preys)

  # Replace number of preys
  filita<-str_locate(modulo_modificado[inicio:final], pattern="nprey")
  fila<-modulo_modificado[inicio:final][!is.na(filita[,1])]
  cual<-which(!is.na(filita[,1])==T)
  nueva<-str_replace(fila, as.character(old_nprey), as.character(length(names_preys[prefs!=0])))    
  modulo_modificado[inicio:final][cual]<-nueva
  
  ## For the next three blocks (preferences, calcifier and coupling) the template needs to have
  ## enough lines available to accommodate all preys, otherwise it does not complain but later fails at runtime.
  
  # Block of preferences
  new_lines<- paste("      suprey", c(1:length(names_preys[prefs!=0])),":   ",
                    round(prefs[prefs!=0],2),"     # [-]           Availability of ",
                    names_preys[prefs!=0]," to Z4_",clasesP[j], sep="")
  modulo_modificado[inicio:final][(cual+1):(cual+length(names_preys[prefs!=0]))]<-new_lines
  
  # Block of calcifier (1 if the prey is calcifier, 0 otherwise)
  boolean<-rep(0,length(names_preys[prefs!=0]))
  boolean[(substring(names_preys[prefs!=0],1,2)=="P5")]<-1
  new_lines<- paste("      isP2", c(1:length(names_preys[prefs!=0])),":   ",
                    boolean,"           # [-] identify P2 among the preys [1 for P2 otherwise 0]", sep="")
  filita<-str_locate(modulo_modificado[inicio:final], pattern="isP21:")
  fila<-modulo_modificado[inicio:final][!is.na(filita[,1])]
  cual<-which(!is.na(filita[,1])==T)    
  modulo_modificado[inicio:final][(cual):(cual-1+length(names_preys[prefs!=0]))]<-new_lines    
  
  # Block of coupling
  new_lines<- paste("      prey", c(1:length(names_preys[prefs!=0])),":   ",names_preys[prefs!=0], sep="")
  filita<-str_locate(modulo_modificado[inicio:final], pattern="      prey1:")
  fila<-modulo_modificado[inicio:final][!is.na(filita[,1])]
  cual<-which(!is.na(filita[,1])==T)    
  modulo_modificado[inicio:final][(cual):(cual-1+length(names_preys[prefs!=0]))]<-new_lines    
  
  # Compute initialization values
  if (equal_initialization==T) {carbon<-c_init_zoo
  } else if (separate_phytozoo==T) { carbon<-c_init_zoo/counts_init_zoo[which(clases_init==clasesP[j])]
  } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
  nitro <-carbon*(0.015/1.0)
  phosp <-carbon*(0.00167/1.0)
  # Locate initialization chunk
  location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
  cual<-which(!is.na(location[,1])==T)
  # Replace initialization values
  nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "1.0", as.character(round(carbon,10)))
  modulo_modificado[inicio:final][cual+1]<-nueva
  nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.015", as.character(round(nitro,10)))
  modulo_modificado[inicio:final][cual+2]<-nueva
  nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.00167", as.character(round(phosp,10)))
  modulo_modificado[inicio:final][cual+3]<-nueva
  
  # Add to fabm.yaml
  write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
              sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
  
} # end loop j size_classes  
################### 



####### Z3 ########
old_name<-"Z3"
old_long_name<-"Carnivorous Mesozooplankton"
clasesP<-clasesZ3  # ESD = (((2^clasesP)*6)/pi)^(1/3)
# parameters to modify
parameters<-c("p_sum") #,"p_vum")
old_value <-c("2.0") #,"0.008")
digito<-c(2) #,4)

old_nprey=3
# Which are the potential preys? define here all prey choices non size-dependent, if it is not here, it is not eaten
preys<-c(clasesP1,clasesP4,clasesZ5,clasesZ4,clasesZ3)
names_preys<-c(paste("P1",str_replace(clasesP1,pattern="-","m"), sep="_"),
               paste("P4",str_replace(clasesP4,pattern="-","m"), sep="_"),
               paste("Z5",str_replace(clasesZ5,pattern="-","m"), sep="_"),
               paste("Z4",str_replace(clasesZ4,pattern="-","m"), sep="_"),
               paste("Z3",str_replace(clasesZ3,pattern="-","m"), sep="_"))

################### 
# Extract module block for Px
inicio<-grep(paste("start",old_name,sep="_"), modulo_productor)
final<-grep(paste("end",old_name,sep="_"), modulo_productor)
#modulo_productor[inicio:final]

lower_limit<-2^clasesP      
upper_limit<-2^(clasesP+1)
diameter<-((lower_limit*6)/pi)^(1/3)

# Compute parameter values (in case allometry is taxon-dependent)
# Maximum ingestion rate (Hansen 1997) d-1
p_sum<- 26*diameter^(-0.4)
# Half-saturation constant (Banas 2011) 3 uM N
p_chuc<-rep(3*(106/16)*12,length=length(diameter)) # not used
p_vum<-p_sum/p_chuc                                # not used

# Find optimum prey size
diameter_preys<-(((2^preys)*6)/pi)^(1/3) # range(diameter_preys) 
#opt_prey_size<-0.65*diameter^0.56   # Hansen 1997, Banas 2011
opt_prey_size<-diameter/10           # Kiorboe 2008, Ward 2012
## Carbon conversion from Alcaraz et al 2003 (zooplankton<5mm WM) not used for now
#carbon<-0.0825+0.0780*lower_limit

## For each instance (j): substitute old parameter value (i) by new one,
##                        compute preferences and fill prey blocks,
##                        write the Z brick in the fabm.yaml

for (j in c(1:length(clasesP))){    
  modulo_modificado<-modulo_productor
  diameter_min<-((lower_limit[j]*6)/pi)^(1/3)
  # Replace module names
  nueva<-str_replace(modulo_modificado[inicio:final][3], old_name, paste(old_name,"_",clasesP[j],sep=""))
  modulo_modificado[inicio:final][3]<-nueva
  nueva<-str_replace(modulo_modificado[inicio:final][5], old_long_name, paste(old_long_name,"_",round(diameter_min,1),sep=""))
  modulo_modificado[inicio:final][5]<-nueva
  
  # Replace parameter values    
  for (i in c(1:length(parameters))){
    parametro<-get(parameters[i])
    # Locate the parameter
    filas<-str_locate(modulo_modificado[inicio:final], pattern=paste(parameters[i],":",sep=""))
    fila<-modulo_modificado[inicio:final][!is.na(filas[,1])]
    cual<-which(!is.na(filas[,1])==T)
    # Replace parameter value
    nueva<-str_replace(fila, old_value[i], as.character(round(parametro[j],digito[i])))
    modulo_modificado[inicio:final][cual]<-nueva
  } # end loop i parameters
  
  # Compute preferences, sd of the gaussian = 0.20
  prefs<-exp(-((log10(diameter_preys)-log10(opt_prey_size[j]))/0.20)^2)
  prefs[prefs<1e-2]<-0
  o<-order(2^preys)

  # Replace number of preys
  filita<-str_locate(modulo_modificado[inicio:final], pattern="nprey")
  fila<-modulo_modificado[inicio:final][!is.na(filita[,1])]
  cual<-which(!is.na(filita[,1])==T)
  nueva<-str_replace(fila, as.character(old_nprey), as.character(length(names_preys[prefs!=0])))    
  modulo_modificado[inicio:final][cual]<-nueva
  
  ## For the next three blocks (preferences, calcifier and coupling) the template needs to have
  ## enough lines available to accommodate all preys, otherwise it does not complain but later fails at runtime.
  
  # Block of preferences
  new_lines<- paste("      suprey", c(1:length(names_preys[prefs!=0])),":   ",
                    round(prefs[prefs!=0],2),"     # [-]           Availability of ",
                    names_preys[prefs!=0]," to Z3_",clasesP[j], sep="")
  modulo_modificado[inicio:final][(cual+1):(cual+length(names_preys[prefs!=0]))]<-new_lines
  
  # Block of calcifier (1 if the prey is calcifier, 0 otherwise)
  boolean<-rep(0,length(names_preys[prefs!=0]))
  boolean[(substring(names_preys[prefs!=0],1,2)=="P5")]<-1
  new_lines<- paste("      isP2", c(1:length(names_preys[prefs!=0])),":   ",
                    boolean,"           # [-] identify P2 among the preys [1 for P2 otherwise 0]", sep="")
  filita<-str_locate(modulo_modificado[inicio:final], pattern="isP21:")
  fila<-modulo_modificado[inicio:final][!is.na(filita[,1])]
  cual<-which(!is.na(filita[,1])==T)    
  modulo_modificado[inicio:final][(cual):(cual-1+length(names_preys[prefs!=0]))]<-new_lines    
  
  # Block of coupling
  new_lines<- paste("      prey", c(1:length(names_preys[prefs!=0])),":   ",names_preys[prefs!=0], sep="")
  filita<-str_locate(modulo_modificado[inicio:final], pattern="      prey1:")
  fila<-modulo_modificado[inicio:final][!is.na(filita[,1])]
  cual<-which(!is.na(filita[,1])==T)    
  modulo_modificado[inicio:final][(cual):(cual-1+length(names_preys[prefs!=0]))]<-new_lines    
  
  # Compute initialization values
  if (equal_initialization==T) {carbon<-c_init_zoo
  } else if (separate_phytozoo==T) { carbon<-c_init_zoo/counts_init_zoo[which(clases_init==clasesP[j])]
  } else {carbon<-c_init/counts_init[which(clases_init==clasesP[j])]}
  nitro <-carbon*(0.015/1.0)
  phosp <-carbon*(0.00167/1.0)
  # Locate initialization chunk
  location<-str_locate(modulo_modificado[inicio:final], pattern="initialization:")
  cual<-which(!is.na(location[,1])==T)
  # Replace initialization values
  nueva<-str_replace(modulo_modificado[inicio:final][cual+1], "1.0", as.character(round(carbon,4)))
  modulo_modificado[inicio:final][cual+1]<-nueva
  nueva<-str_replace(modulo_modificado[inicio:final][cual+2], "0.015", as.character(round(nitro,5)))
  modulo_modificado[inicio:final][cual+2]<-nueva
  nueva<-str_replace(modulo_modificado[inicio:final][cual+3], "0.00167", as.character(round(phosp,7)))
  modulo_modificado[inicio:final][cual+3]<-nueva

  # Add to fabm.yaml
  write.table(modulo_modificado[(inicio+1):(final-1)], file=paste(outdir,nombrito,".yaml",sep=""),
              sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
  
} # end loop j size_classes  
###################



sink()   



### Few things to do by hand to the newly created fabm.yaml (nombrito)
### - replace all occurrences of B1 by B1_m5
### - change initialization values of B1_m5 to 0.1 mgC m-3
### - replace all values of p_minfood to zero (this switches off the density-dependent predation in Z6 and Z5)
## Optional
## - PHYTO: all p_q10 to 2 and p_temp to 0 (no difference in temperature traits)
## - set compute_acdom, compute_anap, compute_bph and compute_bbc to TRUE



############################################################
## MIGRATION formally as relaxation towards reference values
############################################################

# Block to put in the gotm.yaml, in the section fabm:input:
relax_tau=157680000       # 5 years  relax_tau/60/60/24/365
redondear=10

sink(paste(outdir,str_replace(nombrito, "fabm", "gotm"),".yaml",sep=""))

# PHYTO
old_names<-c(rep("P1",length(clasesP1)),rep("P2",length(clasesP2)),rep("P3",length(clasesP3)),rep("P4",length(clasesP4)),rep("P5",length(clasesP5)),
             rep("P6",length(clasesP6)),rep("P7",length(clasesP7)),rep("P8",length(clasesP8)),rep("P9",length(clasesP9)))
clasesP<-c(clasesP1,clasesP2,clasesP3,clasesP4,clasesP5,clasesP6,clasesP7,clasesP8,clasesP9)# (((2^clasesP)*6)/pi)^(1/3)
clasesN<-str_replace(clasesP,pattern="-","m")
nombrecitos<-paste(old_names,"_",clasesN,sep="")
out_phyto<-nombrecitos

## Depending on the initialization choices:
## it splits the initial biomass for the bin into the types represented in the bin,
## and creates small blocks like this (where x is the elemental component):
##       name_instance/x:
##          constant_value: x.xxx   # [mg/m3]
##          relax_tau: 157680000    # [s]

      for (k in c(1:length(clasesN))){
        
        # Values
        if (equal_initialization==T) {carbon<-c_init_phyto
        } else if (separate_phytozoo==T) { carbon<-c_init_phyto/counts_init_phyto[which(clases_init==clasesP[k])]
        } else {carbon<-c_init/counts_init[which(clases_init==clasesP[k])]}
        carbon<-round(carbon,redondear)
        nitro <-round(carbon*(0.1008/8),redondear)
        phosp <-round(carbon*(0.006288/8),redondear)
        silica <-round(carbon*(0.08/8),redondear)
        chloro<-round(carbon*(0.16/8),redondear)        
        
        #CARBON
        datosCARBON  <- rbind(paste("      ",nombrecitos[k],"/c:", sep=""),
                             paste("         constant_value: ",carbon, "   # [mg/m3]", sep=""),
                             paste("         relax_tau: ",relax_tau, "   # [s]", sep=""))
        write.table(datosCARBON, file=paste(outdir,str_replace(nombrito, "fabm", "gotm"),".yaml",sep=""),
                    sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      
        #NITRO
        datosNITRO   <- rbind(paste("      ",nombrecitos[k],"/n:", sep=""),
                              paste("         constant_value: ",nitro, "   # [mmol/m3]", sep=""),
                              paste("         relax_tau: ",relax_tau, "   # [s]", sep=""))
        write.table(datosNITRO, file=paste(outdir,str_replace(nombrito, "fabm", "gotm"),".yaml",sep=""),
                    sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      
        #PHOSP
        datosPHOSP   <- rbind(paste("      ",nombrecitos[k],"/p:", sep=""),
                              paste("         constant_value: ",phosp, "   # [mmol/m3]", sep=""),
                              paste("         relax_tau: ",relax_tau, "   # [s]", sep=""))
        write.table(datosPHOSP, file=paste(outdir,str_replace(nombrito, "fabm", "gotm"),".yaml",sep=""),
                    sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
       
        #CHLORO
        datosCHLORO   <- rbind(paste("      ",nombrecitos[k],"/Chl:", sep=""),
                              paste("         constant_value: ",chloro, "   # [mg/m3]", sep=""),
                              paste("         relax_tau: ",relax_tau, "   # [s]", sep=""))
        write.table(datosCHLORO, file=paste(outdir,str_replace(nombrito, "fabm", "gotm"),".yaml",sep=""),
                    sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
        
        #SILICA
        if (k<=length(clasesP1)){
          datosSILICA  <- rbind(paste("      ",nombrecitos[k],"/s:", sep=""),
                                paste("         constant_value: ",silica, "   # [mmol/m3]", sep=""),
                                paste("         relax_tau: ",relax_tau, "   # [s]", sep=""))
          # Add to gotm.yaml
          write.table(datosSILICA, file=paste(outdir,str_replace(nombrito, "fabm", "gotm"),".yaml",sep=""),
                      sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
        }   
        } # end loop k phyto types  

        
# MICRO
old_names<-c(rep("Z6",length(clasesZ6)),rep("Z5",length(clasesZ5)))
clasesP<-c(clasesZ6,clasesZ5)
clasesN<-str_replace(clasesP,pattern="-","m")
nombrecitos<-paste(old_names,"_",clasesN,sep="")
out_micro<-nombrecitos

      for (k in c(1:length(clasesN))){
        
        # Values
        if (equal_initialization==T) {carbon<-c_init_zoo
        } else if (separate_phytozoo==T) { carbon<-c_init_zoo/counts_init_zoo[which(clases_init==clasesP[k])]
        } else {carbon<-c_init/counts_init[which(clases_init==clasesP[k])]}
        carbon<-round(carbon,redondear)
        nitro <-round(carbon*(0.0167/1.0),redondear)
        phosp <-round(carbon*(0.00185/1.0),redondear)        
        
        #CARBON
        datosCARBON  <- rbind(paste("      ",nombrecitos[k],"/c:", sep=""),
                              paste("         constant_value: ",carbon, "   # [mg/m3]", sep=""),
                              paste("         relax_tau: ",relax_tau, "   # [s]", sep=""))
        write.table(datosCARBON, file=paste(outdir,str_replace(nombrito, "fabm", "gotm"),".yaml",sep=""),
                    sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
        
        #NITRO
        datosNITRO   <- rbind(paste("      ",nombrecitos[k],"/n:", sep=""),
                              paste("         constant_value: ",nitro, "   # [mmol/m3]", sep=""),
                              paste("         relax_tau: ",relax_tau, "   # [s]", sep=""))
        write.table(datosNITRO, file=paste(outdir,str_replace(nombrito, "fabm", "gotm"),".yaml",sep=""),
                    sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
        
        #PHOSP
        datosPHOSP   <- rbind(paste("      ",nombrecitos[k],"/p:", sep=""),
                              paste("         constant_value: ",phosp, "   # [mmol/m3]", sep=""),
                              paste("         relax_tau: ",relax_tau, "   # [s]", sep=""))
        write.table(datosPHOSP, file=paste(outdir,str_replace(nombrito, "fabm", "gotm"),".yaml",sep=""),
                    sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
        } # end loop k micro types  

# MESO
old_names<-c(rep("Z4",length(clasesZ4)),rep("Z3",length(clasesZ3)))
clasesP<-c(clasesZ4,clasesZ3)
clasesN<-str_replace(clasesP,pattern="-","m")
nombrecitos<-paste(old_names,"_",clasesN,sep="")
out_meso<-nombrecitos

      for (k in c(1:length(clasesN))){

        # Values
        if (equal_initialization==T) {carbon<-c_init_zoo
        } else if (separate_phytozoo==T) { carbon<-c_init_zoo/counts_init_zoo[which(clases_init==clasesP[k])]
        } else {carbon<-c_init/counts_init[which(clases_init==clasesP[k])]}
        carbon<-round(carbon,redondear)
        nitro <-round(carbon*(0.015/1.0),redondear)
        phosp <-round(carbon*(0.00167/1.0),redondear)                
        
        #CARBON
        datosCARBON  <- rbind(paste("      ",nombrecitos[k],"/c:", sep=""),
                              paste("         constant_value: ",carbon, "   # [mg/m3]", sep=""),
                              paste("         relax_tau: ",relax_tau, "   # [s]", sep=""))
        write.table(datosCARBON, file=paste(outdir,str_replace(nombrito, "fabm", "gotm"),".yaml",sep=""),
                    sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
        
        #NITRO
        datosNITRO   <- rbind(paste("      ",nombrecitos[k],"/n:", sep=""),
                              paste("         constant_value: ",nitro, "   # [mmol/m3]", sep=""),
                              paste("         relax_tau: ",relax_tau, "   # [s]", sep=""))
        write.table(datosNITRO, file=paste(outdir,str_replace(nombrito, "fabm", "gotm"),".yaml",sep=""),
                    sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
        
        #PHOSP
        datosPHOSP   <- rbind(paste("      ",nombrecitos[k],"/p:", sep=""),
                              paste("         constant_value: ",phosp, "   # [mmol/m3]", sep=""),
                              paste("         relax_tau: ",relax_tau, "   # [s]", sep=""))
        write.table(datosPHOSP, file=paste(outdir,str_replace(nombrito, "fabm", "gotm"),".yaml",sep=""),
                    sep = "  ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      } # end loop k meso types  

  sink()



### BLOCK of state variables (one element) for phyto and zoo
### useful to get a limited amount of output
### in gotm.yaml, in  output: result: variables: