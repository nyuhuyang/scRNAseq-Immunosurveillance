library(readxl)
library(dplyr)
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
df <- read_excel("output/Cell_numbers.xlsx") %>% as.data.frame
rownames(df) = df[,1]
df = df[,-1]

fisher_p_value = c()
chisq_p_value = c()
ColSum <- colSums(df)

for(i in 1:nrow(df)){
        conting <- rbind(df[i,],ColSum-df[i,])
        FISH <- fisher.test(conting,conf.int = T)
        fisher_p_value[i] = FISH$p.value
        
        CHI = chisq.test(conting, correct = T)
        chisq_p_value[i] = CHI$p.value             
}

df$fisher_p_value = fisher_p_value
df$chisq_p_value = chisq_p_value
write.csv(df,paste0(path,"cell_numbers.csv"))
