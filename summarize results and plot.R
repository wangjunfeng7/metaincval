

Ref_bin <- c("Binary","Binary","Binary","Binary")
Ref_bin <- data.frame(Ref_bin,beta_x7_all,c(0,0,0,0))
names(Ref_bin) <- c("Variable","Beta", "AUC_ref")

Ref_con <- c("Continuous","Continuous","Continuous","Continuous")
Ref_con <- data.frame(Ref_con,beta_x7_all,c(0,0,0,0))
names(Ref_con) <- c("Variable","Beta", "AUC_ref")


for (j in 1:length(beta_x7_all)){
  Ref_bin$AUC_ref[j] <- performance[[j]]$auc_6_bin[1]
  Ref_con$AUC_ref[j] <- performance[[j]]$auc_6_con[1]
}














# continuous
con_est_0_all <- rbind(con_est_0_0_200_all,con_est_0_0_500_all,con_est_0_0_1000_all,
                       con_est_0_04_200_all,con_est_0_04_500_all,con_est_0_04_1000_all,
                       con_est_0_07_200_all,con_est_0_07_500_all,con_est_0_07_1000_all,
                       con_est_0_09_200_all,con_est_0_09_500_all,con_est_0_09_1000_all)



hline = aggregate(con_est_0_all$Ref,
                  by = list(con_est_0_all$Beta),
                  FUN = mean)

names(hline)[1] <- "Beta"


ggplot(con_est_0_all,aes(x= N.Study, y=Est, fill=Strategy)) +
  geom_boxplot()+
  geom_hline(data= hline, aes(yintercept=x),linetype="dashed", color = "red")+
  facet_grid(vars(Beta), vars(N.Sample),scales = "free")







con_est_0_all$SE <- (con_est_0_all$Est-con_est_0_all$Ref)^2


con_est_0_RMSE <- aggregate(con_est_0_all$SE,
                by = list(con_est_0_all$N.Sample,con_est_0_all$N.Study,con_est_0_all$Strategy,con_est_0_all$Beta),
                FUN = mean)


con_est_0_RMSE$RMSE <- con_est_0_RMSE$x^0.5
con_est_0_RMSE

names(con_est_0_RMSE) <- c("N.Sample","N.Study","Strategy","Beta","x","RMSE")

ggplot(con_est_0_RMSE, aes(x=N.Study, y=RMSE, group=Strategy)) + 
  geom_line(aes(color=Strategy),size=1,alpha=0.7)+
  geom_point(aes(color=Strategy),size=2,alpha=0.4)+
  facet_grid(vars(Beta), vars(N.Sample),scales = "free")





# binary


bin_est_0_all <- rbind(bin_est_0_0_200_all,bin_est_0_0_500_all,bin_est_0_0_1000_all,
                       bin_est_0_04_200_all,bin_est_0_04_500_all,bin_est_0_04_1000_all,
                       bin_est_0_07_200_all,bin_est_0_07_500_all,bin_est_0_07_1000_all,
                       bin_est_0_09_200_all,bin_est_0_09_500_all,bin_est_0_09_1000_all)



hline = aggregate(bin_est_0_all$Ref,
                  by = list(bin_est_0_all$Beta),
                  FUN = mean)

names(hline)[1] <- "Beta"


ggplot(bin_est_0_all,aes(x= N.Study, y=Est, fill=Strategy)) +
  geom_boxplot()+
  geom_hline(data= hline, aes(yintercept=x),linetype="dashed", color = "red")+
  facet_grid(vars(Beta), vars(N.Sample),scales = "free")







bin_est_0_all$SE <- (bin_est_0_all$Est-bin_est_0_all$Ref)^2


bin_est_0_RMSE <- aggregate(bin_est_0_all$SE,
                            by = list(bin_est_0_all$N.Sample,bin_est_0_all$N.Study,bin_est_0_all$Strategy,bin_est_0_all$Beta),
                            FUN = mean)


bin_est_0_RMSE$RMSE <- bin_est_0_RMSE$x^0.5
bin_est_0_RMSE

names(bin_est_0_RMSE) <- c("N.Sample","N.Study","Strategy","Beta","x","RMSE")

ggplot(bin_est_0_RMSE, aes(x=N.Study, y=RMSE, group=Strategy)) + 
  geom_line(aes(color=Strategy),size=1,alpha=0.7)+
  geom_point(aes(color=Strategy),size=2,alpha=0.4)+
  facet_grid(vars(Beta), vars(N.Sample),scales = "free")





  
reference <- read.csv("Delta.csv")
reference$Beta <- as.factor(reference$Beta)
reference$Correlation <- as.factor(reference$Correlation)


ggplot(reference, aes(x=Correlation, y=Reference, group=Beta)) + 
  geom_line(aes(color=Beta),size=1,alpha=0.7)+
  geom_point(aes(color=Beta),size=2,alpha=0.4)+
  facet_grid(cols=vars(Variable),scales = "free")+
  ylab("△AUC")




reference <- read.csv("References04.csv")
reference$Beta <- as.factor(reference$Beta)
reference$Correlation <- as.factor(reference$Correlation)


ggplot(reference, aes(x=Correlation, y=Reference, group=Beta)) + 
  geom_line(aes(color=Beta),size=1,alpha=0.7)+
  geom_point(aes(color=Beta),size=2,alpha=0.4)+
  facet_grid(cols=vars(Variable),scales = "free")+
  ylab("AUC (Reference Model)")








