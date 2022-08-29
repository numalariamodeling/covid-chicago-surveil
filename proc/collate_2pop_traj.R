library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(argparser)

p <- arg_parser("Collate trajectories from realizations")
p <- add_argument(p, "--outname", help="Output filename")
p <- add_argument(p, "--rootid", type="numeric",
                  help="Root of serial numbers, e.g., 90000 correspond to ID 90001, 90002, ...")
p <- add_argument(p, "--n", type="numeric",
                  help="Number of realizations", default=500)


args <- parse_args(p)
ini <- args$rootid + 1
fin <- args$rootid + args$n

bigdf <- data.frame()
for (i in ini:fin) {
  df <- data.table::fread(paste0("/projects/b1139/covid-age-output/2pop-real/daily_output.", i))
  df$serial <- i
  bigdf <- bind_rows(bigdf, df)
}

bigdf$date <- bigdf$time + ymd("2020-02-13")
dates <- unique(bigdf$date) %>%
  sort

bigdf1 <- bigdf %>% group_by(serial, node) %>%
  mutate(new_sym = c(NA, diff(cumu_sym)), new_adm = c(NA, diff(cumu_adm)))

bigdf2 <- bigdf1 %>% group_by(serial, date) %>% 
  summarise(HOS = sum(HOS),
            hosp = sum(HOS) + sum(CRIT), 
            new_sym = sum(new_sym), 
            new_adm = sum(new_adm)) %>%
  group_by(serial)

bigdf2 <- bigdf2 %>%
  na.omit() %>%
  filter(date >= ymd("2020-03-01"))

# data.table::fwrite(bigdf2 %>% na.omit(), "out/simout-exp.csv")
data.table::fwrite(bigdf2 %>% na.omit(), paste0("out/", args$outname, ".csv"))

# bigdf2
