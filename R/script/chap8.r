## ----nikkeidata---------------------------------------------------------------
library(tidyverse)
dt <- read.csv("data/nikkei225.csv", stringsAsFactor = FALSE)
dt$date <- as.Date(dt$date, format = "%Y/%m/%d")
head(dt)


## ----TimeSeriesPlot, echo = FALSE, fig.cap = "Time Series Data of Nikkei 225 (Closed Price)", out.width="90%", fig.align = 'center'----
ggplot(dt, aes(x = date, y = close_price)) + 
  geom_line(color = "blue") +
  geom_vline(aes(xintercept = as.Date("2019/12/31")), linetype = 2) +
  geom_vline(aes(xintercept = as.Date("2020/01/16")), linetype = 2) +
  geom_vline(aes(xintercept = as.Date("2020/11/07")), linetype = 2) +
  annotate(
    geom = "rect",
    xmin = as.Date("2020/04/07"), xmax = as.Date("2020/05/25"), 
    ymin = -Inf, ymax = Inf,
    fill = "grey80", alpha = 0.4
  ) +
  annotate(
    geom = "text", 
    x = as.Date("2020/05/20"), y = 27500, 
    label = "First declaration of \n a state of emergency \n in Japan", 
    hjust = "center"
  ) +
  annotate(
    geom = "segment", 
    x = as.Date("2019/11/30"), y = 27500, 
    xend = as.Date("2019/12/31"), yend = 27500, 
    arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(
    geom = "text", 
    x = as.Date("2019/11/29"), y = 27500, 
    label = "First case of Covid-19 in China", 
    hjust = "right"
  ) +
  annotate(
    geom = "segment", 
    x = as.Date("2019/11/30"), y = 26500, 
    xend = as.Date("2020/01/16"), yend = 26500, 
    arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(
    geom = "text", 
    x = as.Date("2019/11/29"), y = 26500, 
    label = "First case of Covid-19 in Japan", 
    hjust = "right"
  ) +
  annotate(
    geom = "segment", 
    x = as.Date("2020/10/07"), y = 25500, 
    xend = as.Date("2020/11/07"), yend = 25500, 
    arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(
    geom = "text", 
    x = as.Date("2020/10/06"), y = 25500, 
    label = "U.S. Presidential Election", 
    hjust = "right"
  ) +
  scale_x_date(date_breaks = "2 months", date_labels = "%y/%m") +
  labs(x = "Date (Year/Month)", y = "Closed price") +
  theme_minimal() +
  theme(
    # setting: background
    plot.background = element_rect(
      #fill="#87CEEB50", 
      color = "transparent"
    ),
    
    # setting: plot
    panel.border = element_rect(color = "white", fill = NA),
    panel.background = element_rect(fill = "white"),   
    panel.grid = element_line(color = "grey80", size = 0.1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # setting: text
    plot.title = element_text(hjust=0.5,size=20),       
    plot.caption = element_text(size=8),       
    
    # setting: axis
    axis.text = element_text(color="black",size=8),    
    axis.title = element_text(size=10),
    axis.ticks.length = unit(0.25, "cm"),
    axis.ticks.x = element_line(size = 0.1),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    
    # setting: legend
    legend.title = element_text(size=8),               
    legend.text = element_text(size=8),                
    legend.key.size = unit(0.5,"cm"),
    legend.background = element_rect(color = "white"), 
    legend.position = "bottom"
  )


## ----AR1, eval = FALSE--------------------------------------------------------
## ar1 <- sarima(dt$close_price, p = 1, q = 0, d = 0, no.constant = TRUE)




## ----result_AR1---------------------------------------------------------------
sprintf(
  "The estimated beta is %1.5f (s.e. = %1.5f)", 
  ar1$fit$coef, sqrt(ar1$fit$var.coef)
)
sprintf("The estimated squared sigma is %1.2f", ar1$fit$sigma2)


## ----DFtest-------------------------------------------------------------------
library(tseries)
df1 <- adf.test(dt$close_price, k = 1)
sprintf("The DF stats is %1.4f (p-value = %1.4f)", df1$statistic, df1$p.value)


## ----ConstructLogReturn-------------------------------------------------------
dt <- dt %>% 
  arrange(date) %>% 
  mutate(lag_close_price = lag(close_price, n = 1)) %>% 
  mutate(log_return = log(close_price) - log(lag_close_price)) %>% 
  filter(!is.na(log_return))
head(dt[,c("date", "close_price", "lag_close_price", "log_return")])


## ----TimeSeriesPlot2, echo = FALSE, fig.cap = "Time Series Data of Nikkei 225 (Log Return)", out.width="90%", fig.align = 'center'----
ggplot(dt, aes(x = date, y = log_return)) + 
  geom_line(color = "blue") +
  geom_vline(aes(xintercept = as.Date("2019/12/31")), linetype = 2) +
  geom_vline(aes(xintercept = as.Date("2020/01/16")), linetype = 2) +
  geom_vline(aes(xintercept = as.Date("2020/11/07")), linetype = 2) +
  annotate(
    geom = "rect",
    xmin = as.Date("2020/04/07"), xmax = as.Date("2020/05/25"), 
    ymin = -Inf, ymax = Inf,
    fill = "grey80", alpha = 0.4
  ) +
  annotate(
    geom = "text", 
    x = as.Date("2020/05/20"), y = 0.08, 
    label = "First declaration of \n a state of emergency \n in Japan", 
    hjust = "center"
  ) +
  annotate(
    geom = "segment", 
    x = as.Date("2019/11/30"), y = 0.07, 
    xend = as.Date("2019/12/31"), yend = 0.07, 
    arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(
    geom = "text", 
    x = as.Date("2019/11/29"), y = 0.07, 
    label = "First case of Covid-19 in China", 
    hjust = "right"
  ) +
  annotate(
    geom = "segment", 
    x = as.Date("2019/11/30"), y = 0.06, 
    xend = as.Date("2020/01/16"), yend = 0.06, 
    arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(
    geom = "text", 
    x = as.Date("2019/11/29"), y = 0.06, 
    label = "First case of Covid-19 in Japan", 
    hjust = "right"
  ) +
  annotate(
    geom = "segment", 
    x = as.Date("2020/10/07"), y = 0.1, 
    xend = as.Date("2020/11/07"), yend = 0.1, 
    arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(
    geom = "text", 
    x = as.Date("2020/10/06"), y = 0.1, 
    label = "U.S. Presidential Election", 
    hjust = "right"
  ) +
  scale_x_date(date_breaks = "2 months", date_labels = "%y/%m") +
  labs(x = "Date (Year/Month)", y = "Log Return") +
  theme_minimal() +
  theme(
    # setting: background
    plot.background = element_rect(
      #fill="#87CEEB50", 
      color = "transparent"
    ),
    
    # setting: plot
    panel.border = element_rect(color = "white", fill = NA),
    panel.background = element_rect(fill = "white"),   
    panel.grid = element_line(color = "grey80", size = 0.1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # setting: text
    plot.title = element_text(hjust=0.5,size=20),       
    plot.caption = element_text(size=8),       
    
    # setting: axis
    axis.text = element_text(color="black",size=8),    
    axis.title = element_text(size=10),
    axis.ticks.length = unit(0.25, "cm"),
    axis.ticks.x = element_line(size = 0.1),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    
    # setting: legend
    legend.title = element_text(size=8),               
    legend.text = element_text(size=8),                
    legend.key.size = unit(0.5,"cm"),
    legend.background = element_rect(color = "white"), 
    legend.position = "bottom"
  )


## ----AR2, eval = FALSE--------------------------------------------------------
## ar2 <- sarima(dt$log_return, p = 1, q = 0, d = 0, no.constant = TRUE)




## ----result_AR2---------------------------------------------------------------
sprintf(
  "The estimated beta is %1.5f (s.e. = %1.5f)", 
  ar2$fit$coef, sqrt(ar2$fit$var.coef)
)
sprintf("The estimated squared sigma is %1.4f", ar2$fit$sigma2)


## ----DFtest2------------------------------------------------------------------
x <- dt[dt$date != as.Date("2019/01/04"), "log_return"]
df2 <- adf.test(x, k = 1)
sprintf("The DF stats is %1.4f (p-value = %1.4f)", df2$statistic, df2$p.value)

