####Figure S3D-E

library(ggplot2)
library(readxl)
axisSetting <- theme(  # remove grid line
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.border = element_blank(),
  axis.line.y = element_line(colour = "black"),
  axis.line.x = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size = 10, color = 'black', face = 'bold', family = "serif"),
  axis.title.y = element_blank(),
  axis.text.x = element_text(size = 8, color = 'black', face = 'bold', family = "serif"),
  axis.text.y = element_text(size = 8, color = 'black', family = "serif"),
  #legend.position ="none",
  legend.title=element_text(size = 8, color = 'black', face = 'bold', family = "serif"),
  plot.title = element_text(hjust = 0.5,color = 'black',face = 'bold', family = "serif")
)


# go_up&down
data<-read_excel("GO_2710up_down.xlsx",sheet = "sheet1")

ggplot(data,aes(reorder(Description, Ratio),Ratio,fill=ONTOLOGY))+
  geom_col()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black",size=8,face = 'bold', family = "serif"),
        axis.line.x = element_line(color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 8, color = 'black', face = 'bold', family = "serif"),
        legend.position = 'none')+
  coord_flip()+
  geom_segment(aes(y=0, yend=0,x=0,xend=40))+
  geom_text(data = data[which(data$Ratio>0),],aes(x=Description, y=-0.01, label=Description),
            hjust=1, size=3, family = "serif")+
  geom_text(data = data[which(data$Ratio<0),],aes(x=Description, y=0.01, label=Description),
            hjust=0, size=3,family = "serif")+
  geom_text(data = data[which(data$Ratio>0),],aes(label=pvalue),
            hjust=-0.1, size=3, color='darkred',family = "serif")+
  geom_text(data = data[which(data$Ratio<0),],aes(label=pvalue),
            hjust=1.1, size=3, color="darkred",family = "serif")+
  #scale_fill_manual(values = c("#5599ff",
  #                             "#ff8800"))+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  ylim(-0.1, 0.1)+
  labs(x='', y='Ratio')


#all
data<-read.csv("GO_2710all.csv")
ggplot(data=data, aes(x=Description,y=Ratio, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8)
Description_order=factor(as.integer(rownames(data)),labels=data$Description)

ggplot(data=data, aes(x=Description_order,y=Ratio, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) + 
  coord_flip() +  
  xlab("GO term") + 
  ylab("Ratio") + 
  theme_bw()+
  axisSetting 

