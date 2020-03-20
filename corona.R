library(tidyverse)
# using SIR Model
# see https://de.wikipedia.org/wiki/SIR-Modell
# see http://www.mathe.tu-freiberg.de/~wegert/Lehre/Seminar3/moehler.pdf

obs_time <- 1 # observation time fram in years

N <- 8000000 # inhabitants
nm_daily_growthrate <- 0.25 # growth rate without any measures (eg, lockdown)
wanted_daily_growthrate <- 0.12 # assumed growth rate with applied measures
release_measures <- 12 # time from now measures end


# 20.3.2020 Austria
# ~ 2400 recolgnized cases, assumption only 14% of total
detected_cases <- 2400
percentage_of_detection <- 0.14

# Startingpoint
# expected cases !
S0 <- N - (detected_cases / percentage_of_detection )  
I0 <- detected_cases / percentage_of_detection 
R0 <- 0
D0 <- 0

# assumption detected cases and expected case double within same time frame
# neglecting cure and death for initial calibration
#
# a * I0_bar * S0_bar = I0_bar * nm_daily_growthrate
# a * S0_bar = nm_daily_growthrate
# a = nm_daily_growthrate / S0_bar
#
a0 <- (nm_daily_growthrate) / (N - detected_cases) # infection rate no measures
a1 <- (wanted_daily_growthrate) / (N - detected_cases) # infection rate wanted

rel_m<- 30*release_measures # days from now lockdown released

a <- a1     # we are now in lockdown
b <- 1/12 #0.07   # cure rate - takes 7 to 17 days to cure
c <- 0.0012 # death rate


state <- tibble(i = 0,
                  S = S0,
                  I = I0,
                  R = R0,
                  D = D0)

S_dot <- function(a, I, S) -a*I*S
I_dot <- function(a, b, c, I, S) a*I*S - b*I - c*I
R_dot <- function(b, I) b*I
D_dot <- function(c, I) c*I


for (l in 0:(365*obs_time)) { #(365*2)
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

# to long format for ggplot2
state %>% gather('S', 'I', 'R', 'D' , key = type, value = pax) %>% 
  filter(type %in% c('S', 'I', 'R', 'D')) %>% 
  ggplot(aes(x = i, y = pax , group = type, color = type )) + geom_line() 

mortality_rate <- state %>% slice(365*obs_time) %>% select(D) %>% as.numeric() / state %>% slice(365*obs_time) %>% select(R) %>% as.numeric()
max_infected <- state %>% summarise(max_infected = max(I)) %>% as.numeric()
total_died <- state %>% summarise(total_died = max(D)) %>% as.numeric()
percentage_immune <- state %>% summarise(total_recovered = max(R)) %>% as.numeric() / N
total_recovered <- state %>% summarise(total_recovered = max(R)) %>% as.numeric()

# print summary
cat('Summary: \n
# of people which were infected: ', total_recovered + total_died, ' 
# of people which recovered: ', total_recovered ,'
# of people which died: ', total_died ,'
mortality rate: ', mortality_rate, '
maximum concurrent infected: ', max_infected, '
    thereof needing hospitalization: ', max_infected * 0.2, '
    thereof needing intense care: ', max_infected * 0.05
)