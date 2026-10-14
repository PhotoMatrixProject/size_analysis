#Start from dataset dl from a pipeline, here i use area (most relevant from the perspective of competition between the images on a page)

## ============================================================
## 1. CREATE ISSUE-LEVEL DATA
## ============================================================

# A unique issue must be defined within its magazine and volume.
# Dates help distinguish identically numbered issues.
issue_ok <- !is.na(d_fil$magazineLabel) &
  !is.na(d_fil$issueLabel)

issue_key <- with(
  d_fil,
  paste(
    magazineLabel,
    ifelse(is.na(volumeLabel), "<no volume>", volumeLabel),
    issueLabel,
    startTime,
    endTime,
    sep = " | "
  )
)

summarise_groups <- function(data, key) {
  groups <- split(seq_len(nrow(data)), key)
  
  result <- lapply(names(groups), function(g) {
    ii <- groups[[g]]
    
    data.frame(
      group       = g,
      n_images    = length(ii),
      mean_area   = if (all(is.na(data$area[ii]))) NA_real_
      else mean(data$area[ii], na.rm = TRUE),
      mean_time   = mean(data$midTime[ii], na.rm = TRUE),
      midpoint    = mean(
        rowMeans(
          cbind(
            as.numeric(data$startTime[ii]),
            as.numeric(data$endTime[ii])
          ),
          na.rm = TRUE
        ),
        na.rm = TRUE
      )
    )
  })
  
  do.call(rbind, result)
}

issues <- summarise_groups(
  d_fil[issue_ok, ],
  issue_key[issue_ok]
)

issues$midpoint <- as.Date(
  issues$midpoint,
  origin = "1970-01-01"
)

issues$year <- as.integer(format(issues$midpoint, "%Y"))

summary(issues$n_images)
nrow(issues)

## ============================================================
## 2. ANNUAL AVERAGE NUMBER OF IMAGES PER ISSUE
## ============================================================

annual_mean <- aggregate(
  n_images ~ year,
  data = issues,
  FUN = mean
)

annual_median <- aggregate(
  n_images ~ year,
  data = issues,
  FUN = median
)

annual_n <- aggregate(
  n_images ~ year,
  data = issues,
  FUN = length
)

names(annual_mean)[2]   <- "mean_images"
names(annual_median)[2] <- "median_images"
names(annual_n)[2]      <- "n_issues"

annual <- Reduce(
  function(x, y) merge(x, y, by = "year"),
  list(annual_mean, annual_median, annual_n)
)

annual

plot(
  issues$year,
  issues$n_images,
  pch = 16,
  col = adjustcolor("grey40", alpha.f = 0.20),
  xlab = "Year",
  ylab = "Number of images per issue",
  main = "Number of images per issue over time"
)

points(
  annual$year,
  annual$mean_images,
  pch = 16,
  col = "navy"
)

# Flexible trend; quasipoisson is suitable for overdispersed counts
trend_spline <- glm(
  n_images ~ bs(year, df = 5),
  family = quasipoisson(link = "log"),
  data = issues
)

year_grid <- data.frame(
  year = seq(
    min(issues$year, na.rm = TRUE),
    max(issues$year, na.rm = TRUE),
    length.out = 500
  )
)

year_grid$predicted <- predict(
  trend_spline,
  newdata = year_grid,
  type = "response"
)

lines(
  year_grid$year,
  year_grid$predicted,
  col = "firebrick",
  lwd = 3
)

legend(
  "topleft",
  legend = c("Individual issues", "Annual mean", "Spline trend"),
  col = c("grey40", "navy", "firebrick"),
  pch = c(16, 16, NA),
  lty = c(NA, NA, 1),
  lwd = c(NA, NA, 3),
  bty = "n"
)

direction <- ifelse(
  diff(year_grid$predicted) >= 0,
  "up",
  "down"
)

runs <- rle(direction)
run_end <- cumsum(runs$lengths)
run_start <- c(1, head(run_end, -1) + 1)

trend_intervals <- data.frame(
  from = round(year_grid$year[run_start], 1),
  to = round(year_grid$year[pmin(run_end + 1, nrow(year_grid))], 1),
  direction = runs$values
)

trend_intervals


trend_linear <- glm(
  n_images ~ year,
  family = quasipoisson(link = "log"),
  data = issues
)

summary(trend_linear)

annual_percent_change <-
  100 * (exp(coef(trend_linear)["year"]) - 1)

annual_percent_change


## ============================================================
## 3. MAGAZINE-, VOLUME-, AND ISSUE-LEVEL SUMMARIES
## ============================================================

mag_ok <- !is.na(d_fil$magazineLabel)

magazine_key <- as.character(d_fil$magazineLabel)

magazines <- summarise_groups(
  d_fil[mag_ok, ],
  magazine_key[mag_ok]
)


vol_ok <- !is.na(d_fil$magazineLabel) &
  !is.na(d_fil$volumeLabel)

volume_key <- with(
  d_fil,
  paste(magazineLabel, volumeLabel, sep = " | ")
)

volumes <- summarise_groups(
  d_fil[vol_ok, ],
  volume_key[vol_ok]
)

# The issue-level table has already been created
head(issues)


correlation_results <- function(x) {
  keep <- complete.cases(x$n_images, x$mean_area)
  x <- x[keep, ]
  
  pearson <- cor.test(
    x$n_images,
    x$mean_area,
    method = "pearson"
  )
  
  spearman <- cor.test(
    x$n_images,
    x$mean_area,
    method = "spearman",
    exact = FALSE
  )
  
  data.frame(
    n_groups   = nrow(x),
    pearson_r  = unname(pearson$estimate),
    pearson_p  = pearson$p.value,
    spearman_r = unname(spearman$estimate),
    spearman_p = spearman$p.value
  )
}

correlations <- rbind(
  magazine = correlation_results(magazines),
  volume   = correlation_results(volumes),
  issue    = correlation_results(issues)
)

correlations

par(mfrow = c(1, 3))

plot(
  magazines$n_images,
  magazines$mean_area,
  log = "x",
  pch = 16,
  xlab = "Number of images (log scale)",
  ylab = "Average relative area",
  main = "Magazine level"
)

plot(
  volumes$n_images,
  volumes$mean_area,
  log = "x",
  pch = 16,
  xlab = "Number of images (log scale)",
  ylab = "Average relative area",
  main = "Volume level"
)

plot(
  issues$n_images,
  issues$mean_area,
  log = "x",
  pch = 16,
  xlab = "Number of images (log scale)",
  ylab = "Average relative area",
  main = "Issue level"
)

par(mfrow = c(1, 1))



## ============================================================
## SUMMARY OF MAGAZINE COVERAGE
## ============================================================

magazine_groups <- split(
  seq_len(nrow(d_fil)),
  d_fil$magazineLabel
)

magazine_summary <- lapply(names(magazine_groups), function(magazine) {
  
  ii <- magazine_groups[[magazine]]
  
  volume_values <- d_fil$volumeLabel[ii]
  volume_values <- volume_values[
    !is.na(volume_values) & volume_values != ""
  ]
  
  issue_ids <- paste(
    d_fil$volumeLabel[ii],
    d_fil$issueLabel[ii],
    d_fil$startTime[ii],
    d_fil$endTime[ii],
    sep = " | "
  )
  
  valid_issues <- !is.na(d_fil$issueLabel[ii]) &
    d_fil$issueLabel[ii] != ""
  
  issue_ids <- issue_ids[valid_issues]
  
  data.frame(
    magazine       = magazine,
    first_year     = min(
      as.integer(format(d_fil$startTime[ii], "%Y")),
      na.rm = TRUE
    ),
    last_year      = max(
      as.integer(format(d_fil$endTime[ii], "%Y")),
      na.rm = TRUE
    ),
    n_volumes      = length(unique(volume_values)),
    n_issues       = length(unique(issue_ids)),
    n_reproductions = length(ii)
  )
})

magazine_summary <- do.call(rbind, magazine_summary)

magazine_summary <- magazine_summary[
  order(
    magazine_summary$first_year,
    magazine_summary$magazine
  ),
]

row.names(magazine_summary) <- NULL

print(magazine_summary, row.names = FALSE)

write.table(
  magazine_summary,
  "magazine_summary.txt",
  row.names = FALSE,
  sep="\t",
  fileEncoding = "UTF-8"
)
