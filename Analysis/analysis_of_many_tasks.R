library(tidyverse)
library(Rfast)

raw_data <- read_delim("../Data/data.csv",
                       "&", escape_double = FALSE, trim_ws = TRUE)


###################################################################
#make a plot with a line for each number of tasks. Mode version with error bars
n_tasks = c(1,5,25,100,200,300,400,500,2000)
data1 <- raw_data %>% filter(number_of_tasks %in% n_tasks)  %>% 
  select(number_of_tasks, contains("mode_veridical?"))

temp <- data1 %>% select(contains("mode_veridical?")) %>% gather(variable, value) %>%
  separate(variable, into = c("x","y","z", "time"), sep="_")

#data is 0 or 1 so use prop.test to obtain confidence interval
GetLowerCI <- function(x,y){return(prop.test(x,y)$conf.int[1])}
GetTopCI <- function(x,y){return(prop.test(x,y)$conf.int[2])}

number_of_tasks <- rep(data1$number_of_tasks, dim(data1)[2]-1)
to_plot <- cbind(number_of_tasks, temp) %>% group_by(time, number_of_tasks) %>% 
  summarize(Samples=n(),Hits=sum(value, na.rm=T),Mean=mean(value, na.rm=T),Lower=GetLowerCI(Hits,Samples),Top=GetTopCI(Hits,Samples))
g <- ggplot(data = to_plot, aes(x = as.numeric(time), y = Mean, ymin = Lower, ymax = Top, color = as.factor(number_of_tasks))) + 
  geom_ribbon() + geom_line() + ylim(0,1) + guides(color = guide_legend(reverse = TRUE)) +
  ylab("proportion runs where mode strategy is veridical")
ggsave("plot.pdf", g)
