library(tidyverse)

files <- read_lines("data/read_counts/files.txt")
lines <- read_lines("data/read_counts/lines.txt") |> as.integer()

d <- tibble(file = files, lines = lines) |>
    mutate(id = str_match(file, "s_([:alnum:]+)_1")[, 2]) |>
    group_by(id) |>
    summarise(reads = sum(lines) / 4)

summary(d)

with(d, c(mean(reads), sd(reads)))
