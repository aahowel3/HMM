setwd("~/Documents/HMM")


data2state=read.csv("viterbi_nb_25_2state_validated.tsv",sep="",row.names = NULL)

data3state=read.csv("viterbi_nb_25_3state_validated_rerun_state1.tsv",sep="",row.names = NULL)

data3state=read.csv("viterbi_nb_25_3state_validated_state2.tsv",sep="",row.names = NULL)

data3state=read.csv("viterbi_nb_25_3state_validated_rerun_state3.tsv",sep="",row.names = NULL)

data3state=read.csv("viterbi_nb_25_3state_rerun_state1or3.tsv",sep="",row.names = NULL)



data2state$state_assignment="2 state"
data3state$state_assignment="3 state"


longform = rbind(data2state,data3state)

ggplot(longform, aes(value, fill=state_assignment)) + geom_histogram(position="dodge") + ggtitle("Detectable IESs") + theme(plot.title = element_text(hjust = 0.5))

