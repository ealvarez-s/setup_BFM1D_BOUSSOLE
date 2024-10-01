library(igraph)

## dev_1D_diverse
##### 60 PFTs
        if (Sys.info()['sysname']=="Windows") {OS<-"C:/"} else {OS<-paste("/Users/",Sys.info()['user'],"/", sep="")}
        o_dir<-"Datos/Res_C29_R2/"
        my_example<-"realdate_EXP-21_pmax_unimodal_aPN_quotaScaled_bps_photoNO_final_theta2_bX333_b_thetaM_long"
        archivo<-paste(my_example,"/fabm.yaml",sep="")
        filename <- paste(OS,o_dir,"RUNS_100PFTs", "/",archivo, sep="")
        texto<-readLines(filename)
        
        # Read numbers/names of consumers and preys, and build the adjacency matrix
        total_FTs<-c(grep("model: ogs/Phyto",texto),
                     grep("model: ogs/PelBac",texto),
                     grep("model: ogs/MicroZoo",texto),
                     grep("model: ogs/MesoZoo",texto))
        bloque<-c(texto[total_FTs-3],texto[total_FTs-2],texto[total_FTs-1])
        nombres_FTs<-bloque[grep("[A-Z][0-9]:|[A-Z][0-9]_",bloque)]
        nombres_FTs<-gsub(":", "", nombres_FTs)
        nombres_FTs<-gsub("[[:space:]]", "", nombres_FTs)
        
        inicio<-(str_locate(nombres_FTs, pattern="_")[,1])+1
        clasesG<-substring(nombres_FTs,inicio,inicio+2)
        clasesG<-as.numeric(str_replace(clasesG,"m","-"))        
        
        # Phyto large to small,B, Zoo small to large
        nombres_ordenados<-c(nombres_FTs[grep("P4",nombres_FTs)][order(clasesG[grep("P4",nombres_FTs)],decreasing=T)],
                             nombres_FTs[grep("P1",nombres_FTs)][order(clasesG[grep("P1",nombres_FTs)],decreasing=T)],
                             nombres_FTs[grep("P5",nombres_FTs)][order(clasesG[grep("P5",nombres_FTs)],decreasing=T)],
                             nombres_FTs[grep("P2",nombres_FTs)][order(clasesG[grep("P2",nombres_FTs)],decreasing=T)],
                             nombres_FTs[grep("P8",nombres_FTs)][order(clasesG[grep("P8",nombres_FTs)],decreasing=T)],
                             nombres_FTs[grep("P7",nombres_FTs)][order(clasesG[grep("P7",nombres_FTs)],decreasing=T)],
                             nombres_FTs[grep("P3",nombres_FTs)][order(clasesG[grep("P3",nombres_FTs)],decreasing=T)],
                             nombres_FTs[grep("P9",nombres_FTs)][order(clasesG[grep("P9",nombres_FTs)],decreasing=T)],
                             nombres_FTs[grep("P6",nombres_FTs)][order(clasesG[grep("P6",nombres_FTs)],decreasing=T)],
                             nombres_FTs[grep("B1",nombres_FTs)][order(clasesG[grep("B1",nombres_FTs)],decreasing=T)],
                             nombres_FTs[grep("Z6",nombres_FTs)][order(clasesG[grep("Z6",nombres_FTs)],decreasing=F)],
                             nombres_FTs[grep("Z5",nombres_FTs)][order(clasesG[grep("Z5",nombres_FTs)],decreasing=F)],
                             nombres_FTs[grep("Z4",nombres_FTs)][order(clasesG[grep("Z4",nombres_FTs)],decreasing=F)],
                             nombres_FTs[grep("Z3",nombres_FTs)][order(clasesG[grep("Z3",nombres_FTs)],decreasing=F)])
                              
        adj_matrix<-matrix(0,ncol=length(total_FTs),nrow=length(total_FTs))
        colnames(adj_matrix)<-nombres_ordenados
        rownames(adj_matrix)<-nombres_ordenados

        inicio<-(str_locate(nombres_ordenados, pattern="_")[,1])+1
        clasesG<-substring(nombres_ordenados,inicio,inicio+2)
        clasesG<-as.numeric(str_replace(clasesG,"m","-"))
        
        
        colorinchis<-c(rep("hotpink3",10),rep("lightpink2",10),
                       rep("darkgoldenrod4",10),rep("gold2",10),
                       rep("green4",6),rep("darkolivegreen2",6),
                       rep("darkslategray3",3),rep("deepskyblue4",3),rep("darkorchid3",2),
                       rep("black",1),
                       rep("grey20",7),
                       rep("grey40",10),
                       rep("grey60",10),
                       rep("grey80",10))

        # Fill the matrix one consumer type at a time
        instances<-grep("model: ogs/MicroZoo",texto)
        #length(instances)
        bloque<-c(texto[instances-3],texto[instances-2],texto[instances-1])
        cazadores<-bloque[grep("[A-Z][0-9]:|[A-Z][0-9]_",bloque)]
        cazadores<-gsub(":", "", cazadores)
        cazadores<-gsub("[[:space:]]", "", cazadores)
        
        # Microzoo
        for (i in c(1:length(instances))){
          cazador=cazadores[i]
          if (i!=length(instances)) {largo<-(instances[i+1]-instances[i])} else {largo=(grep("model: ogs/MesoZoo",texto)[1])-instances[i]}
          
          presas<-grep("prey",texto[instances[i]:(instances[i]+largo)])[-1]
          nombres<-presas[(1+(2*(length(presas)/3))):(length(presas))] 
          preferencias<-grep("suprey",texto[instances[i]:(instances[i]+largo)])
          
          donde<-str_locate(texto[instances[i]:(instances[i]+largo)][nombres], pattern=":")
          nombritos<-substring(texto[instances[i]:(instances[i]+largo)][nombres],donde[,1]+1, donde[,1]+9)
          nombritos<-gsub("[[:space:]]", "", nombritos)
          
          donde<-str_locate(texto[instances[i]:(instances[i]+largo)][preferencias], pattern=":")
          prefis<-substring(texto[instances[i]:(instances[i]+largo)][preferencias],donde[,1]+1, donde[,1]+9)
          prefis<-as.numeric(gsub("[[:space:]]", "", prefis))
          
          # Fill the matrix
          cuales<-match(nombritos,rownames(adj_matrix))
          adj_matrix[cuales,match(cazador,colnames(adj_matrix))]<-prefis }
        
        # Fill the matrix one consumer type at a time
        instances<-grep("model: ogs/MesoZoo",texto)
        #length(instances)
        bloque<-c(texto[instances-3],texto[instances-2],texto[instances-1])
        cazadores<-bloque[grep("[A-Z][0-9]:|[A-Z][0-9]_",bloque)]
        cazadores<-gsub(":", "", cazadores)
        cazadores<-gsub("[[:space:]]", "", cazadores)
        
        # Mesozoo
        for (i in c(1:length(instances))){
          cazador=cazadores[i]
          if (i!=length(instances)) {largo<-(instances[i+1]-instances[i])} else {largo=length(texto)-instances[i]}
          presas<-grep("prey",texto[instances[i]:(instances[i]+largo)])[-1]
          nombres<-presas[(1+(2*(length(presas)/3))):(length(presas))] 
          preferencias<-grep("suprey",texto[instances[i]:(instances[i]+largo)])
          
          donde<-str_locate(texto[instances[i]:(instances[i]+largo)][nombres], pattern=":")
          nombritos<-substring(texto[instances[i]:(instances[i]+largo)][nombres],donde[,1]+1, donde[,1]+10)
          nombritos<-gsub("[[:space:]]", "", nombritos)
          
          donde<-str_locate(texto[instances[i]:(instances[i]+largo)][preferencias], pattern=":")
          prefis<-substring(texto[instances[i]:(instances[i]+largo)][preferencias],donde[,1]+1, donde[,1]+9)
          prefis<-as.numeric(gsub("[[:space:]]", "", prefis))
          
          # Fill the matrix
          cuales<-match(nombritos,rownames(adj_matrix))
          adj_matrix[cuales,match(cazador,colnames(adj_matrix))]<-prefis }
        
        
## Plotting with igraph        
        #par(mfrow=c(1,1))
        #data<-adj_matrix
        #data[data!=0]<-1.0
        #network <- graph_from_adjacency_matrix(adj_matrix , mode='directed') #, weighted = TRUE)
        network <- graph_from_adjacency_matrix(adj_matrix , mode='undirected', weighted = TRUE)
        pesos<-as.vector(t(adj_matrix))
        pesos<-pesos[pesos!=0]
        
        
        plot.igraph(network, layout=layout.circle, main="circle",
                    vertex.size=log10(2^clasesG)+4, vertex.label.dist=0.0,
                    vertex.color=colorinchis, edge.arrow.size=0.5, vertex.label.cex=0.5,
                    vertex.label.col="white", edge.width=pesos)
        
        plot.igraph(network, layout=layout.fruchterman.reingold, main="fruchterman.reingold",
                    vertex.size=log10(2^clasesG)+4, vertex.label.dist=0.0,
                    vertex.color=colorinchis, edge.arrow.size=0.5, vertex.label.cex=0.1,
                    vertex.label.col="white", edge.width=pesos)

        plot.igraph(network, layout=layout.lgl, main="large graph layout",
                    vertex.size=log10(2^clasesG)+4, vertex.label.dist=0.0,
                    vertex.color=colorinchis, edge.arrow.size=0.5, vertex.label.cex=0.5,
                    vertex.label.col="white", edge.width=pesos)
        
########################## 
## Other options for plotting with igraph
## It opens interactive map    
      
        tkplot(network, canvas.width = 450, canvas.height = 450)
        
        tkplot(network, canvas.width = 600, canvas.height = 600,layout=layout.fruchterman.reingold, main="fruchterman.reingold",
               vertex.size=log10(2^clasesG)+4, vertex.label.dist=0.5,
               vertex.color=colorinchis, edge.arrow.size=0.5)          
        
        tkplot(network, canvas.width = 600, canvas.height = 600,
               layout=layout.lgl, main="large graph layout",
               vertex.size=log10(2^clasesG)+4, vertex.label.dist=0.0,
               vertex.color=colorinchis, edge.arrow.size=0.5, vertex.label.cex=0.5,
               vertex.label.col="white", edge.width=pesos)
        #handtune the placement of the vertices
        #query the coordinates by the tk_coords()
        #chordDiagram(adj_matrix, transparency = 0.5, directional = 0, scale=T)
###########################
        