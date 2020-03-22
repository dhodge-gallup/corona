# Modelling the Corona Epideic using a SIRD model
#
# author: Florian Nehonsky (florian@nehonsky.net)
# last change 2020-03-23
# published under GPLv3
#
# Models:
# see https://de.wikipedia.org/wiki/SIR-Modell
# see http://www.mathe.tu-freiberg.de/~wegert/Lehre/Seminar3/moehler.pdf
#
# Data:
# https://www.worldometers.info/coronavirus/coronavirus-cases/#recovered



library(tidyverse)
library(ggplot2)
library(scales)
library(lubridate)
library(deSolve)

setwd('C:/Users/Florian/Documents/corona')


calc_time <- 2 # calculation time frame in years, [0, calc_time]
obs_time <- 1 # observation time frame in years for plots, [0, obs_time]
use_solver <- 'mySolve' # 'mySolver', 'deSolve'; release measures only works with 'mySolver' but is slower

# ***************
# model specific parameters
N <- 8800000 # inhabitants
# 20.3.2020 Austria
# ~ 2400 recolgnized cases, assumption only 14% of total
today <- as.Date('2020-03-20')
detected_cases <- 2400
release_measures <- as.Date('2020-11-01')  # time measures (eg. lockdown) end
percentage_of_detection <- 0.14 # therefore, assumed darkfigure is (1 - percentage_of_detection)
nm_daily_growthrate <- 0.35     # growth rate without any measures (eg, lockdown)
wanted_daily_growthrate <- 0.12 # assumed growth rate with applied measures
                                # a growth rate of 0.30 .. I doubles within  2.7d
                                #                  0.25 .. I doubles within  3.1d
                                #                  0.20 .. I doubles within  3.8d
                                #                  0.15 .. I doubles within  5.0d
                                #                  0.10 .. I doubles within  7.3d
                                #                  0.07 .. I doubles within 10.3d 
obs_mortality_rate <- 0.0144    # observed mortality rate
obs_intensecare <- 0.02         # observed percentage of infected which needs intense care
                                # currently about 5%
                                # assumption: due to reports from china that mortality rate 
                                #             is about 1.44% instead of ~4-5% due to hidden cases (no symptoms)
                                #             we use a lower and more realistic figure of about 2%
obs_hospit <- 0.20              # observed percentage of infected which needs hospitalisation
# daily cure rate
# we calibrate on observed data, https://www.worldometers.info/coronavirus/coronavirus-cases/#recovered
# March 4th seems like a good point on the slope
#    why: slope is almost flat (epidemic was more or less only in china and china was getting it under control)
#         therefore recovery figures should be resilient
# total cases March 4th: 38 505
# recovered March 4th:   53 524
# recovered March 3th:   50 994
obs_curerate <- (53524 - 50994)/38505   # maybe the cure rate is better, probably the focus is currently on testing of newly infected people ?
av_intensecare_beds <- 2500             # number of available beds in intense care


param_string <- paste0('Date: ', as.Date(today),
                       ',\n N=', N, 
                       ', init_cases=', ceiling(detected_cases / percentage_of_detection) ,
                       ',\n gr=', wanted_daily_growthrate, 
                       ', cr=', round(obs_curerate,4),
                       ', dr=', obs_mortality_rate,
                       ',\n icp=', obs_intensecare,
                       ', rm= ', as.Date(release_measures))

filename <- paste0('./plots/',
                   as.character(today, format='%Y%m%d'),
                   '-g',round(wanted_daily_growthrate,4),
                   '-c',round(obs_curerate,4),
                   '-d',round(obs_mortality_rate,4),
                   '-rm',as.character(release_measures, format='%Y%m%d'))
# ***************
# initial values
# expected cases !
S0 <- N - (detected_cases / percentage_of_detection )  
I0 <- detected_cases / percentage_of_detection 
R0 <- 0
D0 <- 0

# assumption detected cases and expected case grow using the same dynamic
# neglecting cure and death for initial calibration
#
# a * I0_bar * S0_bar = I0_bar * nm_daily_growthrate
# a * S0_bar = nm_daily_growthrate
# a = nm_daily_growthrate / S0_bar
#
a0 <- (nm_daily_growthrate) / (N - detected_cases) # infection rate no measures
a1 <- (wanted_daily_growthrate) / (N - detected_cases) # infection rate wanted

# rel_m<- 30*release_measures # days from now lockdown released

# ***************
# ODE Parameters
a <- a1                     # infection rate, we are now in lockdown
b <- obs_curerate           # cure rate
c <- obs_mortality_rate * b # daily death rate

# # ***************
if (use_solver == 'mySolve') {
  rel_m <- as.numeric(release_measures - today)
  state <- tibble(i = 0,
                  S = S0,
                  I = I0,
                  R = R0,
                  D = D0)
  
  S_dot <- function(a, I, S) -a*I*S
  I_dot <- function(a, b, c, I, S) a*I*S - b*I - c*I
  R_dot <- function(b, I) b*I
  D_dot <- function(c, I) c*I
  
  
  # ***************
  # numerically solve ODE
  for (l in 0:(365*calc_time)) {
    if (l == rel_m) a <<- a0   # release lockdown when time has come
    
    temp <- state %>% slice(l + 1)
    S_l <- temp %>% select(S) %>% as.numeric()
    I_l <- temp %>% select(I) %>% as.numeric()
    R_l <- temp %>% select(R) %>% as.numeric()
    D_l <- temp %>% select(D) %>% as.numeric()
    new <- tibble(i = l + 1,
                  S = S_l + S_dot(a, I_l, S_l),
                  I = I_l + I_dot(a, b, c, I_l, S_l),
                  R = R_l + R_dot(b, I_l),
                  D = D_l + D_dot(c, I_l))
    state <- rbind(state, new)
  }
}

if (use_solver == 'deSolve') {
  parameters <- c(a = a1,
                  b = b,
                  c = c)
  initstate <- c(S = S0, 
                 I = I0, 
                 R = R0, 
                 D = D0)
  
  SIRD <- function(t, initstate, parameters) {
    with(as.list(c(initstate, parameters)),{
      # rate of change
      dS <- -a*S*I
      dI <- a*S*I - b*I - c*I
      dR <- b*I
      dD <- c*I
      # return the rate of change
      list(c(dS, dI, dR, dD))
    }) 
  }
  
  times <- seq(0, calc_time*365, by = 0.01)
  state  <- ode(y = state, times = times, func = SIRD, parms = parameters) %>% 
    as_tibble() %>% 
    select(i = time, everything())
}

# **************
# solution preparation and plots
state <- state %>% 
  select(Susceptibles = S,
         Infected = I,
         Recovered = R,
         Dead = D,
         everything()) %>% 
  mutate(Time = i + today) %>% 
  mutate(intensecare = Infected * obs_intensecare)

mortality_rate <- state %>% slice(365*obs_time) %>% select(Dead) %>% as.numeric() / state %>% slice(365*obs_time) %>% select(Recovered) %>% as.numeric()
max_infected <- state %>% summarise(max_infected = max(Infected)) %>% as.numeric()
total_died <- state %>% summarise(total_died = max(Dead)) %>% as.numeric()
percentage_immune <- state %>% summarise(total_recovered = max(Recovered)) %>% as.numeric() / N
total_recovered <- state %>% summarise(total_recovered = max(Recovered)) %>% as.numeric()

# transform to long format and plot
g <- state %>% 
  gather('Susceptibles', 'Infected', 'Recovered', 'Dead' , 'intensecare', key = type, value = Population) %>% 
  filter(type %in% c('Susceptibles', 'Infected', 'Recovered', 'Dead', 'intensecare')) %>% 
  ggplot(aes(x = Time, y = Population , group = type, color = type)) + 
  geom_line(lwd = 1) +
  labs(title = 'Dynamic of Corona Epidemic',
       subtitle = param_string)+ 
  scale_x_date(date_breaks = "months" , 
               date_labels = "%b-%y", 
               limits = c(today, today + years(obs_time))) + 
  scale_y_continuous(limits = c(0, N), breaks = seq(0, N, 500000)) + 
  geom_hline(yintercept = av_intensecare_beds, linetype="dashed", color = "red") + # draw line of available intense care beds
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    #size = 12, 
    angle = 45)#,
    #axis.text.y = element_text(face = "bold", color = "blue", 
    #                           size = 12, angle = 45)
  )
plot(g)
ggsave(filename = paste0(filename,'_total.png'),
       plot = g,
       device = png(),
       dpi = 'print')

# zoom into peak of infected
g1 <- g + coord_cartesian(xlim = c(today, today + years(obs_time)),
                          ylim = c(0, max_infected )) + 
  scale_y_continuous(limits = c(0, N), breaks = seq(0, N, 100000), labels = comma) +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    #size = 12, 
    angle = 45)#,
    #axis.text.y = element_text(face = "bold", color = "blue", 
    #                           size = 12, angle = 45)
  )
plot(g1)
ggsave(filename = paste0(filename,'_peak.png'),
       plot = g1,
       device = png(),
       dpi = 'print')

# zoom into peak of people needing intense care and 
g2 <- g1 + coord_cartesian(xlim = c(today, today + years(obs_time)),
                          ylim = c(0, max_infected * obs_intensecare )) + 
  scale_y_continuous(limits = c(0, N), breaks = seq(0, N, 2500), labels = comma) +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    #size = 12, 
    angle = 45)#,
    #axis.text.y = element_text(face = "bold", color = "blue", 
    #                           size = 12, angle = 45)
  )
plot(g2)
ggsave(filename = paste0(filename,'_IC.png'),
       plot = g2,
       device = png(),
       dpi = 'print')

# write summary
txt <- paste0('SIRD Model of Corona Epidemic
-----------------------------

Data:
-----
date                          ', today,'
detected cases                ', format(detected_cases, decimal.mark = '.', big.mark = ' ', nsmall = 0), '
rate of detection             ', round(percentage_of_detection*100, 4), '%
growth rate                   ', round(wanted_daily_growthrate*100, 4), '%
growth rate (no measures)     ', round(nm_daily_growthrate*100, 4), '%
observed mortality rate       ', round(obs_mortality_rate*100, 4), '%
observed hospitalisaton rate  ', round(obs_hospit*100, 4), '%
observed intense care rate    ', round(obs_intensecare*100, 4), '%
observed cure rate            ', round(obs_curerate*100, 4), '%
intense care beds             ', format(av_intensecare_beds, decimal.mark = '.', big.mark = ' ', nsmall = 0), '
release containment measures  ', release_measures ,'

ODE Parameters:
---------------
Population    ', format(N, decimal.mark = '.', big.mark = ' ', nsmall = 0), '
a             ', a1, '
b             ', b, '
c             ', c, '
S0            ', format(S0, decimal.mark = '.', big.mark = ' ', nsmall = 2), '
I0            ', format(I0, decimal.mark = '.', big.mark = ' ', nsmall = 2), '

Summary:
--------
Date of peak:                        ', state %>% filter(Infected == max_infected) %>% select(Time) %>% as.numeric %>% as.Date(origin='1970-01-01'), '
# of people which will be infected:  ', format(total_recovered + total_died, decimal.mark = '.', big.mark = ' ', nsmall = 0), ' 
# of people which will recover:      ', format(total_recovered, decimal.mark = '.', big.mark = ' ', nsmall = 0) ,'
# of people which will die:          ', format(total_died, decimal.mark = '.', big.mark = ' ', nsmall = 0) ,'
mortality rate:                      ', format(mortality_rate*100, decimal.mark = '.', big.mark = ' ', nsmall = 2), '%
maximum concurrent infected:         ', format(max_infected, decimal.mark = '.', big.mark = ' ', nsmall = 0), '
    thereof needing hospitalization: ', format(ceiling(max_infected * obs_hospit), decimal.mark = '.', big.mark = ' ', nsmall = 0), '
    thereof needing intense care:    ', format(ceiling(max_infected * obs_intensecare), decimal.mark = '.', big.mark = ' ', nsmall = 0))

fileConn <- file(paste0(filename,'.txt'))
writeLines(txt, fileConn)
close(fileConn)


