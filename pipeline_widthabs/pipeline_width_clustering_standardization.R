#setwd("/home/petr_turecek/size3")
library(splines)
library(rethinking)
library(stringi)

d <- read.csv("size_data_2026_10_03.csv")
nrow(d)
str(d)
range(d$width_mm, na.rm=TRUE)

d$startTime <- as.Date(d$startTime, format="%Y-%m-%dT%H:%M:%OSZ", tz="UTC")
d$endTime <- as.Date(d$endTime, format="%Y-%m-%dT%H:%M:%OSZ", tz="UTC")

origin_date <- as.Date("1900-01-01")
d$startTime_numeric <- as.integer(d$startTime - origin_date) + 1
d$endTime_numeric <- as.integer(d$endTime - origin_date) + 1
d$midTime <- d$startTime_numeric + (d$endTime_numeric - d$startTime_numeric)/2

d_fil <- subset(d, !is.na(startTime) & startTime >= as.Date("1900-01-01") & startTime <= as.Date("1950-12-31"))
range(d_fil$width_mm, na.rm=TRUE)
range(d_fil$midTime, na.rm=TRUE)

mtime <- mean(d_fil$midTime, na.rm=TRUE)
d_fil$midTime <- (d_fil$midTime - mtime)/365
range(d_fil$midTime, na.rm=TRUE)

K <- 5
knots <- seq(min(d_fil$midTime, na.rm=TRUE), max(d_fil$midTime, na.rm=TRUE), l=K)
B <- bs(d_fil$midTime, knots=knots[-c(1,K)], degree=3, intercept=TRUE)

d_fil$CountryLabel <- "Czechia"

dl <- data.frame(
  IDr=as.integer(as.factor(d_fil$item)),
  width=d_fil$width,
  height=d_fil$height,
  width_mm=d_fil$width_mm,
  page_width_mm=d_fil$page_width_mm,
  iss=as.integer(as.factor(d_fil$issueLabel)),
  vol=as.integer(as.factor(d_fil$volumeLabel)),
  mag=as.integer(as.factor(d_fil$magazineLabel)),
  country=as.integer(as.factor(d_fil$CountryLabel)),
  time=d_fil$midTime,
  B1=B[,1], B2=B[,2], B3=B[,3], B4=B[,4], B5=B[,5], B6=B[,6], B7=B[,7]
)

analysis_slug <- "width"
clustering_out_root <- "temporal_width_figures"
standardization_out_root <- "standardization_width"
value_col <- "width_mm"
value_label_plural <- "widths"
y_axis_label <- "Width in mm"
plot_measure_title <- "Absolute Width"
density_from <- 0
density_to <- 241
y_lim <- c(density_from,density_to)
y_ticks <- seq(0,240,by=40)
y_labs <- as.character(y_ticks)
default_blob_halfwidth <- c(2,9)
default_max_jump <- 10
to_display <- function(x) x
format_display_value <- function(x) paste0(round(x,1)," mm")


magazine_names <- levels(as.factor(d_fil$magazineLabel))
dl$magazine_name <- magazine_names[dl$mag]
dl$magazine_name_clean <- stri_trans_general(dl$magazine_name, "Latin-ASCII")

# X-axis labels in real years
xdates <- as.Date(c("1900-01-01","1910-01-01","1920-01-01","1930-01-01","1940-01-01","1950-01-01"))
xdati <- as.integer(xdates - origin_date) + 1
xdats <- (xdati - mtime)/365
labs <- c("1900","1910","1920","1930","1940","1950")

x_min <- min(dl$time, na.rm=TRUE)
x_max <- max(dl$time, na.rm=TRUE)

# ============================================================
# TEMPORAL DENSITY-MODE FIGURES + SUMMARIES
# ============================================================

fmt_tag <- function(x){
  x <- as.character(x)
  x <- gsub("\\.","p",x)
  x <- gsub("-","m",x)
  x
}

sanitize_name <- function(x){
  x <- stri_trans_general(x,"Latin-ASCII")
  x <- gsub("[^A-Za-z0-9_]+","_",x)
  x <- gsub("_+","_",x)
  x
}

choose_legend_pos <- function(d_known,xlim){
  xr <- diff(xlim)
  starts_late <- min(d_known$time,na.rm=TRUE) > xlim[1] + 0.15*xr
  spans_widely <- diff(range(d_known$time,na.rm=TRUE)) > 0.60*xr
  if(starts_late || spans_widely) "topleft" else "topright"
}

find_density_modes <- function(x, from=density_from, to=density_to, adjust=0.7, mode_frac=0.15, n=1024){
  x <- x[!is.na(x)]
  x <- x[x >= from & x <= to]
  if(length(x) < 5) return(data.frame(mode=numeric(0), height=numeric(0)))
  dens <- density(x, from=from, to=to, adjust=adjust, n=n)
  dy <- diff(dens$y)
  maxima <- which(diff(sign(dy)) == -2) + 1
  if(length(maxima) == 0) return(data.frame(mode=numeric(0), height=numeric(0)))
  mode_x <- dens$x[maxima]
  mode_y <- dens$y[maxima]
  keep <- mode_y > max(mode_y,na.rm=TRUE)*mode_frac
  data.frame(mode=mode_x[keep], height=mode_y[keep])
}

temporal_modes <- function(dl, magazine, window=3, step=0.5, adjust=0.7, mode_frac=0.10){
  d <- dl[dl$magazine_name == magazine & !is.na(dl[[value_col]]),]
  d <- d[d[[value_col]] >= density_from & d[[value_col]] <= density_to,]
  if(nrow(d) == 0) return(NULL)
  centers <- seq(min(d$time,na.rm=TRUE), max(d$time,na.rm=TRUE), by=step)
  out <- list()
  for(i in seq_along(centers)){
    lo <- centers[i] - window/2
    hi <- centers[i] + window/2
    dd <- d[d$time >= lo & d$time <= hi,]
    modes <- find_density_modes(dd[[value_col]], adjust=adjust, mode_frac=mode_frac)
    if(nrow(modes) > 0){
      modes$time <- centers[i]
      modes$n_window <- nrow(dd)
      out[[length(out)+1]] <- modes
    }
  }
  if(length(out) == 0) return(NULL)
  do.call(rbind,out)
}

track_modes <- function(tm, max_jump=0.5){
  tm <- tm[order(tm$time,tm$mode),]
  tm$track <- NA_integer_
  next_track <- 1L
  times <- sort(unique(tm$time))
  for(tt in times){
    ii <- which(tm$time == tt)
    if(tt == times[1]){
      tm$track[ii] <- seq_len(length(ii))
      next_track <- max(tm$track[ii]) + 1L
    } else {
      prev <- tm[tm$time < tt & !is.na(tm$track),]
      prev_time <- max(prev$time)
      prev <- prev[prev$time == prev_time,]
      used_prev <- rep(FALSE,nrow(prev))
      for(j in ii){
        dist <- abs(prev$mode - tm$mode[j])
        dist[used_prev] <- Inf
        if(length(dist) > 0 && min(dist) <= max_jump){
          k <- which.min(dist)
          tm$track[j] <- prev$track[k]
          used_prev[k] <- TRUE
        } else {
          tm$track[j] <- next_track
          next_track <- next_track + 1L
        }
      }
    }
  }
  tm
}

make_flat_capped_corridor <- function(x,y,w,ncap=50){
  n <- length(x)
  usr <- par("usr")
  pin <- par("pin")
  xupi <- diff(usr[1:2])/pin[1]
  yupi <- diff(usr[3:4])/pin[2]
  xr <- w*xupi/yupi
  thR <- seq(pi/2,-pi/2,length.out=ncap)
  thL <- seq(3*pi/2,pi/2,length.out=ncap)
  xR <- x[n] + xr[n]*cos(thR)
  yR <- y[n] + w[n]*sin(thR)
  xL <- x[1] + xr[1]*cos(thL)
  yL <- y[1] + w[1]*sin(thL)
  list(x=c(x,xR[-1],rev(x),xL[-1]), y=c(y+w,yR[-1],rev(y-w),yL[-1]))
}

make_param_tag <- function(window,step,mode_adjust,mode_frac,max_jump,min_track_length){
  paste0("w",fmt_tag(window),"_s",fmt_tag(step),"_adj",fmt_tag(mode_adjust),"_frac",fmt_tag(mode_frac),"_jump",fmt_tag(max_jump),"_min",fmt_tag(min_track_length))
}

get_temporal_mode_data <- function(dl, magazine, window=3, step=0.5, mode_adjust=0.7, mode_frac=0.10, max_jump=0.5, min_track_length=3){
  if(magazine == "ALL_MAGAZINES"){
    d <- dl
    dl_tmp <- dl
    dl_tmp$magazine_name <- "ALL_MAGAZINES"
    tm <- temporal_modes(dl_tmp, magazine="ALL_MAGAZINES", window=window, step=step, adjust=mode_adjust, mode_frac=mode_frac)
  } else {
    d <- dl[dl$magazine_name == magazine,]
    tm <- temporal_modes(dl, magazine=magazine, window=window, step=step, adjust=mode_adjust, mode_frac=mode_frac)
  }
  d_known <- d[!is.na(d[[value_col]]),]
  d_known <- d_known[d_known[[value_col]] >= density_from & d_known[[value_col]] <= density_to,]
  if(nrow(d_known) < 5) stop(paste0("Too few known ", value_label_plural, " for this magazine."))
  if(is.null(tm) || nrow(tm) == 0) stop("No temporal modes found.")
  tm <- track_modes(tm, max_jump=max_jump)
  track_tab <- table(tm$track)
  keep_tracks <- as.integer(names(track_tab)[track_tab >= min_track_length])
  tm <- tm[tm$track %in% keep_tracks,]
  if(nrow(tm) == 0) stop("No tracks survived min_track_length.")
  old_tracks <- sort(unique(tm$track))
  for(i in seq_along(old_tracks)) tm$track[tm$track == old_tracks[i]] <- i
  tm$rel_strength <- NA_real_
  for(tt in sort(unique(tm$time))){
    ii <- which(tm$time == tt)
    tm$rel_strength[ii] <- tm$height[ii]/sum(tm$height[ii], na.rm=TRUE)
  }
  tm$mode_display <- to_display(tm$mode)
  list(d=d, d_known=d_known, tm=tm)
}

make_temporal_mode_figure <- function(dl, magazine, outdir_root=clustering_out_root, window=3, step=0.5, mode_adjust=0.7, mode_frac=0.10, max_jump=0.5, min_track_length=3, raw_col="#0080801A", kde_fill="#00808040", blob_fill_alpha=0.28, blob_border_alpha=0.80, white_under_alpha=0.7, line_alpha=0.60, point_cex=0.95, blob_halfwidth=default_blob_halfwidth, xlim=NULL, show_legend=TRUE){
  param_tag <- make_param_tag(window,step,mode_adjust,mode_frac,max_jump,min_track_length)
  outdir <- paste0(outdir_root,"_",param_tag)
  dir.create(outdir, showWarnings=FALSE)
  obj <- get_temporal_mode_data(dl, magazine, window, step, mode_adjust, mode_frac, max_jump, min_track_length)
  d <- obj$d
  d_known <- obj$d_known
  tm <- obj$tm
  n_tracks <- length(unique(tm$track))
  base_cols <- hcl.colors(max(n_tracks,3), palette="Dark 3")[1:n_tracks]
  hmin <- min(tm$height)
  hmax <- max(tm$height)
  if(hmax - hmin < 1e-9){
    tm$halfwidth <- mean(blob_halfwidth)
  } else {
    tm$halfwidth <- blob_halfwidth[1] + (tm$height-hmin)/(hmax-hmin)*(blob_halfwidth[2]-blob_halfwidth[1])
  }
  if(is.null(xlim)) xlim <- range(dl$time, na.rm=TRUE)
  dens_all <- density(d_known[[value_col]], from=density_from, to=density_to, n=1024)
  safe_mag <- sanitize_name(magazine)
  file <- file.path(outdir, paste0(safe_mag,"_temporal_",analysis_slug,"_modes.png"))
  legend_pos <- choose_legend_pos(d_known,xlim)

  h <- 12
  png(file, height=h, width=h*1.6, units="cm", res=600, type="cairo")
  layout(matrix(c(1,2),nrow=1), widths=c(5,1))
  par(mar=c(3.5,4.2,2,1), mgp=c(2.0,0.65,0))
  plot(d$time, d[[value_col]], xlab="Year", ylab="", xaxt="n", yaxt="n", pch=16, cex=point_cex, col=raw_col, ylim=y_lim, xlim=xlim)
  axis(1, at=xdats, labels=labs, las=1)
  axis(2, at=y_ticks, labels=y_labs, las=2)
  mtext(y_axis_label, side=2, line=2.7)
  title(ifelse(magazine=="ALL_MAGAZINES","All magazines",magazine), adj=0)
  abline(h=y_ticks, lty=2, lwd=0.7)
  abline(v=xdats, lty=2, lwd=0.7)

  for(k in sort(unique(tm$track))){
    z <- tm[tm$track == k,]
    z <- z[order(z$time),]
    if(nrow(z) < 2) next
    xx <- seq(min(z$time), max(z$time), length.out=max(120,25*nrow(z)))
    yy <- approx(z$time,z$mode,xout=xx,rule=2)$y
    ww <- approx(z$time,z$halfwidth,xout=xx,rule=2)$y
    rb <- make_flat_capped_corridor(xx,yy,ww,ncap=60)
    polygon(rb$x, rb$y, col=adjustcolor("white",alpha.f=white_under_alpha), border="white")
    polygon(rb$x, rb$y, col=adjustcolor(base_cols[k],alpha.f=blob_fill_alpha), border=adjustcolor(base_cols[k],alpha.f=blob_border_alpha), lwd=1.2)
    lines(xx, yy, col=adjustcolor(base_cols[k],alpha.f=line_alpha), lwd=1.4)
  }

  if(show_legend && magazine!="ALL_MAGAZINES"){
    leg_centers <- tapply(tm$mode_display, tm$track, mean)
    legend(legend_pos, legend=paste0("C",seq_along(leg_centers),": ",format_display_value(as.numeric(leg_centers))), col=adjustcolor(base_cols,alpha.f=blob_border_alpha), lty=1, lwd=1.2, pch=21, pt.bg=adjustcolor(base_cols,alpha.f=blob_fill_alpha), pt.cex=1.1, bty="o", bg="white", box.col=NA, cex=0.75)
  }
  box()

  par(mar=c(3.5,1,2,0.5), mgp=c(2,0.65,0))
  plot(dens_all$y, dens_all$x, type="n", xlab="Density", ylab="", ylim=y_lim, yaxt="n", xaxs="i", xlim=c(0,max(dens_all$y)*1.10))
  axis(2, at=y_ticks, labels=FALSE)
  polygon(c(0,dens_all$y,0), c(dens_all$x[1],dens_all$x,dens_all$x[length(dens_all$x)]), col=kde_fill, border=NA)
  lines(dens_all$y, dens_all$x, lwd=1)
  abline(h=y_ticks, lty=2, lwd=0.7)
  for(k in sort(unique(tm$track))){
    z <- tm[tm$track == k,]
    idx <- sapply(z$mode,function(mm) which.min(abs(dens_all$x - mm)))
    lines(dens_all$y[idx], dens_all$x[idx], col=adjustcolor(base_cols[k],alpha.f=line_alpha), lwd=1.2)
    points(dens_all$y[idx], dens_all$x[idx], pch=21, bg=adjustcolor(base_cols[k],alpha.f=blob_fill_alpha), col=adjustcolor(base_cols[k],alpha.f=blob_border_alpha), cex=0.55 + 0.10*z$halfwidth, lwd=0.9)
  }
  dev.off()
  invisible(list(file=file, tm=tm, outdir=outdir))
}

make_raw_points_figure <- function(dl, magazine, outdir_root=clustering_out_root, window=3, step=0.5, mode_adjust=0.7, mode_frac=0.10, max_jump=0.5, min_track_length=3, raw_col="#0080801A", kde_fill="#00808040", point_cex=0.95, xlim=NULL){
  param_tag <- make_param_tag(window,step,mode_adjust,mode_frac,max_jump,min_track_length)
  outdir <- paste0(outdir_root,"_",param_tag)
  dir.create(outdir, showWarnings=FALSE)
  d <- dl[dl$magazine_name == magazine,]
  d_known <- d[!is.na(d[[value_col]]),]
  d_known <- d_known[d_known[[value_col]] >= density_from & d_known[[value_col]] <= density_to,]
  if(nrow(d_known) < 1) stop(paste0("No known ", value_label_plural, " for this magazine."))
  if(is.null(xlim)) xlim <- range(dl$time, na.rm=TRUE)
  dens_all <- density(d_known[[value_col]], from=density_from, to=density_to, n=1024)
  safe_mag <- sanitize_name(magazine)
  file <- file.path(outdir, paste0(safe_mag,"_temporal_",analysis_slug,"_points_only.png"))

  h <- 12
  png(file, height=h, width=h*1.6, units="cm", res=600, type="cairo")
  layout(matrix(c(1,2),nrow=1), widths=c(5,1))
  par(mar=c(3.5,4.2,2,1), mgp=c(2.0,0.65,0))
  plot(d$time, d[[value_col]], xlab="Year", ylab="", xaxt="n", yaxt="n", pch=16, cex=point_cex, col=raw_col, ylim=y_lim, xlim=xlim)
  axis(1, at=xdats, labels=labs, las=1)
  axis(2, at=y_ticks, labels=y_labs, las=2)
  mtext(y_axis_label, side=2, line=2.7)
  title(magazine, adj=0)
  abline(h=y_ticks, lty=2, lwd=0.7)
  abline(v=xdats, lty=2, lwd=0.7)
  box()

  par(mar=c(3.5,1,2,0.5), mgp=c(2,0.65,0))
  plot(dens_all$y, dens_all$x, type="n", xlab="Density", ylab="", ylim=y_lim, yaxt="n", xaxs="i", xlim=c(0,max(dens_all$y)*1.10))
  axis(2, at=y_ticks, labels=FALSE)
  polygon(c(0,dens_all$y,0), c(dens_all$x[1],dens_all$x,dens_all$x[length(dens_all$x)]), col=kde_fill, border=NA)
  lines(dens_all$y, dens_all$x, lwd=1)
  abline(h=y_ticks, lty=2, lwd=0.7)
  dev.off()
  invisible(list(file=file, outdir=outdir))
}

summarize_temporal_modes <- function(dl, magazine, window=3, step=0.5, mode_adjust=0.7, mode_frac=0.10, max_jump=0.5, min_track_length=3){
  obj <- get_temporal_mode_data(dl, magazine, window, step, mode_adjust, mode_frac, max_jump, min_track_length)
  d <- obj$d_known
  tm <- obj$tm
  tracks <- sort(unique(tm$track))
  track_out <- list()
  for(k in tracks){
    z <- tm[tm$track == k,]
    z <- z[order(z$time),]
    if(nrow(z) >= 2){
      mode_slope <- coef(lm(mode ~ time,data=z))[2]
      strength_slope <- coef(lm(rel_strength ~ time,data=z))[2]
    } else {
      mode_slope <- NA_real_
      strength_slope <- NA_real_
    }
    track_out[[length(track_out)+1]] <- data.frame(
      magazine=magazine,
      track=k,
      mean_mode=mean(z$mode,na.rm=TRUE),
      mean_mode_display=to_display(mean(z$mode,na.rm=TRUE)),
      first_mode=z$mode[1],
      last_mode=z$mode[nrow(z)],
      mode_change=z$mode[nrow(z)] - z$mode[1],
      first_mode_display=to_display(z$mode[1]),
      last_mode_display=to_display(z$mode[nrow(z)]),
      mode_change_display=to_display(z$mode[nrow(z)]) - to_display(z$mode[1]),
      mode_slope_per_year=mode_slope,
      mean_rel_strength=mean(z$rel_strength,na.rm=TRUE),
      first_rel_strength=z$rel_strength[1],
      last_rel_strength=z$rel_strength[nrow(z)],
      rel_strength_change=z$rel_strength[nrow(z)] - z$rel_strength[1],
      rel_strength_slope_per_year=strength_slope,
      first_time=min(z$time,na.rm=TRUE),
      last_time=max(z$time,na.rm=TRUE),
      n_windows=nrow(z)
    )
  }
  track_summary <- do.call(rbind,track_out)
  time_out <- list()
  for(tt in sort(unique(tm$time))){
    z <- tm[tm$time == tt,]
    z$rel_strength <- z$height/sum(z$height,na.rm=TRUE)
    weighted_mean_mode <- sum(z$mode*z$rel_strength,na.rm=TRUE)
    weighted_mean_mode_display <- to_display(weighted_mean_mode)
    large_threshold <- median(track_summary$mean_mode,na.rm=TRUE)
    large_share <- sum(z$rel_strength[z$mode >= large_threshold],na.rm=TRUE)
    time_out[[length(time_out)+1]] <- data.frame(
      magazine=magazine,
      time=tt,
      n_modes=nrow(z),
      weighted_mean_mode=weighted_mean_mode,
      weighted_mean_mode_display=weighted_mean_mode_display,
      large_threshold=large_threshold,
      large_threshold_display=to_display(large_threshold),
      large_share=large_share
    )
  }
  time_summary <- do.call(rbind,time_out)
  fit_wmean <- lm(weighted_mean_mode ~ time,data=time_summary)
  fit_wmean_display <- lm(weighted_mean_mode_display ~ time,data=time_summary)
  fit_lshare <- lm(large_share ~ time,data=time_summary)
  first_ts <- time_summary[which.min(time_summary$time),]
  last_ts <- time_summary[which.max(time_summary$time),]
  magazine_summary <- data.frame(
    magazine=magazine,
    n_observations=nrow(d),
    n_tracks=nrow(track_summary),
    first_time=min(time_summary$time,na.rm=TRUE),
    last_time=max(time_summary$time,na.rm=TRUE),
    first_weighted_mean_mode=first_ts$weighted_mean_mode,
    last_weighted_mean_mode=last_ts$weighted_mean_mode,
    weighted_mean_mode_change=last_ts$weighted_mean_mode - first_ts$weighted_mean_mode,
    weighted_mean_mode_slope_per_year=coef(fit_wmean)[2],
    first_weighted_mean_mode_display=first_ts$weighted_mean_mode_display,
    last_weighted_mean_mode_display=last_ts$weighted_mean_mode_display,
    weighted_mean_mode_change_display=last_ts$weighted_mean_mode_display - first_ts$weighted_mean_mode_display,
    weighted_mean_mode_slope_display_per_year=coef(fit_wmean_display)[2],
    first_large_share=first_ts$large_share,
    last_large_share=last_ts$large_share,
    large_share_change=last_ts$large_share - first_ts$large_share,
    large_share_slope_per_year=coef(fit_lshare)[2],
    mean_track_mode_change=mean(track_summary$mode_change,na.rm=TRUE),
    mean_abs_track_mode_change=mean(abs(track_summary$mode_change),na.rm=TRUE),
    mean_track_mode_change_display=mean(track_summary$mode_change_display,na.rm=TRUE),
    mean_abs_track_mode_change_display=mean(abs(track_summary$mode_change_display),na.rm=TRUE),
    mean_track_strength_change=mean(track_summary$rel_strength_change,na.rm=TRUE),
    mean_abs_track_strength_change=mean(abs(track_summary$rel_strength_change),na.rm=TRUE)
  )
  invisible(list(magazine_summary=magazine_summary, track_summary=track_summary, time_summary=time_summary, tm=tm))
}

run_temporal_mode_set <- function(dl, window=3, step=0.5, mode_adjust=0.7, mode_frac=0.10, max_jump=0.5, min_track_length=3, outdir_root=clustering_out_root, show_legend=FALSE){
  mags <- sort(unique(dl$magazine_name))
  mags <- mags[!is.na(mags)]
  mags <- c("ALL_MAGAZINES",mags)
  param_tag <- make_param_tag(window,step,mode_adjust,mode_frac,max_jump,min_track_length)
  outdir <- paste0(outdir_root,"_",param_tag)
  dir.create(outdir,showWarnings=FALSE)
  out_fig <- setNames(vector("list",length(mags)),mags)
  out_raw <- setNames(vector("list",length(mags)),mags)
  out_sum <- setNames(vector("list",length(mags)),mags)
  for(mag in mags){
    out_fig[[mag]] <- try(make_temporal_mode_figure(dl,magazine=mag,outdir_root=outdir_root,window=window,step=step,mode_adjust=mode_adjust,mode_frac=mode_frac,max_jump=max_jump,min_track_length=min_track_length,show_legend=show_legend),silent=TRUE)
    if(mag != "ALL_MAGAZINES"){
      out_raw[[mag]] <- try(make_raw_points_figure(dl,magazine=mag,outdir_root=outdir_root,window=window,step=step,mode_adjust=mode_adjust,mode_frac=mode_frac,max_jump=max_jump,min_track_length=min_track_length),silent=TRUE)
    }
    out_sum[[mag]] <- try(summarize_temporal_modes(dl,magazine=mag,window=window,step=step,mode_adjust=mode_adjust,mode_frac=mode_frac,max_jump=max_jump,min_track_length=min_track_length),silent=TRUE)
  }
  magazine_summary <- do.call(rbind,lapply(out_sum,function(x) if(!inherits(x,"try-error") && !is.null(x)) x$magazine_summary else NULL))
  track_summary <- do.call(rbind,lapply(out_sum,function(x) if(!inherits(x,"try-error") && !is.null(x)) x$track_summary else NULL))
  time_summary <- do.call(rbind,lapply(out_sum,function(x) if(!inherits(x,"try-error") && !is.null(x)) x$time_summary else NULL))
  if(!is.null(magazine_summary)){
    key_summary <- magazine_summary[,c("magazine","n_tracks","weighted_mean_mode_change_display","large_share_change","mean_track_mode_change_display","mean_abs_track_mode_change_display","mean_track_strength_change")]
    write.csv(magazine_summary,file.path(outdir,paste0("temporal_",analysis_slug,"_magazine_summary.csv")),row.names=FALSE)
    write.csv(key_summary,file.path(outdir,paste0("temporal_",analysis_slug,"_key_summary.csv")),row.names=FALSE)
  }
  if(!is.null(track_summary)) write.csv(track_summary,file.path(outdir,paste0("temporal_",analysis_slug,"_track_summary.csv")),row.names=FALSE)
  if(!is.null(time_summary)) write.csv(time_summary,file.path(outdir,paste0("temporal_",analysis_slug,"_time_summary.csv")),row.names=FALSE)
  failed_figures <- names(out_fig)[sapply(out_fig,inherits,"try-error")]
  failed_raw_figures <- names(out_raw)[sapply(out_raw,inherits,"try-error")]
  failed_summaries <- names(out_sum)[sapply(out_sum,inherits,"try-error")]
  if(length(failed_figures)>0) write.csv(data.frame(failed_figures=failed_figures),file.path(outdir,"failed_figures.csv"),row.names=FALSE)
  if(length(failed_raw_figures)>0) write.csv(data.frame(failed_raw_figures=failed_raw_figures),file.path(outdir,"failed_raw_figures.csv"),row.names=FALSE)
  if(length(failed_summaries)>0) write.csv(data.frame(failed_summaries=failed_summaries),file.path(outdir,"failed_summaries.csv"),row.names=FALSE)
  invisible(list(outdir=outdir,figures=out_fig,raw_figures=out_raw,summaries=out_sum,magazine_summary=magazine_summary,track_summary=track_summary,time_summary=time_summary))
}

# Final clustering parameters (3-year window)
candidate_params <- list(
  candidate_first = list(window=3, step=0.5, mode_adjust=0.7, mode_frac=0.10, max_jump=default_max_jump, min_track_length=3)
)

candidate_out <- list()
for(param_name in names(candidate_params)){
  pars <- candidate_params[[param_name]]
  candidate_out[[param_name]] <- run_temporal_mode_set(
    dl,
    window=pars$window,
    step=pars$step,
    mode_adjust=pars$mode_adjust,
    mode_frac=pars$mode_frac,
    max_jump=pars$max_jump,
    min_track_length=pars$min_track_length,
    show_legend=FALSE
  )
}

# ============================================================
# STANDARDIZATION ANALYSIS
# Yearly cross-magazine diversity, disruption, and within-magazine entropy
# ============================================================

std_n_bins <- 15
std_min_obs_year <- 5
std_min_mags_year <- 2
std_param_tag <- paste0("yearly_bins",std_n_bins)
std_outdir <- paste0(standardization_out_root,"_",std_param_tag)
dir.create(std_outdir,showWarnings=FALSE)

dl$std_value <- dl[[value_col]]
std_value_from <- density_from
std_value_to <- density_to

std_time_to_year <- function(x){
  as.integer(format(origin_date + round(x*365 + mtime - 1), "%Y"))
}
dl$year <- std_time_to_year(dl$time)
std_year_ticks <- seq(1900,1950,10)
std_breaks <- seq(std_value_from,std_value_to,length.out=std_n_bins+1)

std_js_divergence <- function(p,q){
  if(any(is.na(p)) || any(is.na(q))) return(NA_real_)
  if(sum(p) <= 0 || sum(q) <= 0) return(NA_real_)
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  kl_pm <- sum(ifelse(p > 0, p * log(p / m), 0))
  kl_qm <- sum(ifelse(q > 0, q * log(q / m), 0))
  0.5 * kl_pm + 0.5 * kl_qm
}

std_get_bin_probs <- function(x,breaks){
  x <- x[!is.na(x)]
  x <- x[x >= min(breaks) & x <= max(breaks)]
  if(length(x)==0) return(rep(NA_real_,length(breaks)-1))
  tab <- hist(x,breaks=breaks,plot=FALSE,include.lowest=TRUE,right=TRUE)$counts
  if(sum(tab)==0) return(rep(NA_real_,length(tab)))
  tab / sum(tab)
}

std_entropy_from_probs <- function(p, normalize=FALSE){
  p <- p[is.finite(p)]
  p <- p[p > 0]
  if(length(p) == 0) return(NA_real_)
  ent <- -sum(p * log(p))
  if(normalize) ent / 4 else ent
}

std_analysis_data <- function(dl){
  d <- dl[!is.na(dl$std_value) & !is.na(dl$magazine_name) & !is.na(dl$year),]
  d[d$std_value >= std_value_from & d$std_value <= std_value_to,]
}

std_within_magazine_entropy <- function(dl,n_bins=15,min_obs_year=5,breaks=std_breaks){
  d <- std_analysis_data(dl)
  mags <- sort(unique(d$magazine_name))
  years <- sort(unique(d$year))
  year_out <- list()
  for(mag in mags){
    for(yr in years){
      dd <- d[d$magazine_name == mag & d$year == yr,]
      if(nrow(dd) >= min_obs_year){
        probs <- std_get_bin_probs(dd$std_value,breaks=breaks)
        if(!all(is.na(probs))){
          year_out[[length(year_out)+1]] <- data.frame(
            magazine=mag,
            year=yr,
            n_obs=nrow(dd),
            entropy=std_entropy_from_probs(probs,FALSE),
            entropy_norm=std_entropy_from_probs(probs,TRUE)
          )
        }
      }
    }
  }
  yearly <- if(length(year_out)>0) do.call(rbind,year_out) else NULL
  if(is.null(yearly)) return(list(yearly=NULL,time_summary=NULL,magazine_summary=NULL,overall=NULL))
  time_summary <- do.call(rbind,lapply(split(yearly,yearly$year),function(z){
    data.frame(
      year=z$year[1],
      n_magazines=nrow(z),
      mean_entropy=mean(z$entropy,na.rm=TRUE),
      median_entropy=median(z$entropy,na.rm=TRUE),
      mean_entropy_norm=mean(z$entropy_norm,na.rm=TRUE),
      median_entropy_norm=median(z$entropy_norm,na.rm=TRUE),
      sd_entropy_norm=sd(z$entropy_norm,na.rm=TRUE),
      min_entropy_norm=min(z$entropy_norm,na.rm=TRUE),
      max_entropy_norm=max(z$entropy_norm,na.rm=TRUE)
    )
  }))
  time_summary <- time_summary[order(time_summary$year),]
  mag_summary <- do.call(rbind,lapply(split(yearly,yearly$magazine),function(z){
    z <- z[order(z$year),]
    if(nrow(z) >= 2){
      fit <- lm(entropy_norm ~ year,data=z,weights=n_obs)
      entropy_slope <- coef(fit)[2]
    } else {
      entropy_slope <- NA_real_
    }
    data.frame(
      magazine=z$magazine[1],
      n_years=nrow(z),
      n_observations=sum(z$n_obs,na.rm=TRUE),
      first_year=min(z$year,na.rm=TRUE),
      last_year=max(z$year,na.rm=TRUE),
      first_entropy_norm=z$entropy_norm[1],
      last_entropy_norm=z$entropy_norm[nrow(z)],
      entropy_norm_change=z$entropy_norm[nrow(z)] - z$entropy_norm[1],
      entropy_slope_per_year=entropy_slope,
      mean_entropy=mean(z$entropy,na.rm=TRUE),
      mean_entropy_norm=mean(z$entropy_norm,na.rm=TRUE)
    )
  }))
  mag_summary <- mag_summary[order(mag_summary$magazine),]
  overall <- data.frame(
    n_magazines=nrow(mag_summary),
    mean_entropy_slope_per_year=mean(mag_summary$entropy_slope_per_year,na.rm=TRUE),
    median_entropy_slope_per_year=median(mag_summary$entropy_slope_per_year,na.rm=TRUE),
    mean_entropy_norm=mean(mag_summary$mean_entropy_norm,na.rm=TRUE),
    median_entropy_norm=median(mag_summary$mean_entropy_norm,na.rm=TRUE)
  )
  list(yearly=yearly,time_summary=time_summary,magazine_summary=mag_summary,overall=overall)
}

std_entropy_out <- std_within_magazine_entropy(dl,n_bins=std_n_bins,min_obs_year=std_min_obs_year,breaks=std_breaks)
std_last_issue <- aggregate(time ~ magazine_name, data=dl[!is.na(dl$magazine_name) & !is.na(dl$time),], FUN=max, na.rm=TRUE)
names(std_last_issue) <- c("magazine","last_issue_time")
std_last_issue$last_issue_year <- std_time_to_year(std_last_issue$last_issue_time)
if(!is.null(std_entropy_out$magazine_summary)){
  std_entropy_out$magazine_summary <- merge(std_entropy_out$magazine_summary,std_last_issue,by="magazine",all.x=TRUE)
}
write.csv(std_entropy_out$yearly,file.path(std_outdir,"standardization_within_magazine_entropy_yearly.csv"),row.names=FALSE)
write.csv(std_entropy_out$time_summary,file.path(std_outdir,"standardization_within_magazine_entropy_time_summary.csv"),row.names=FALSE)
write.csv(std_entropy_out$magazine_summary,file.path(std_outdir,"standardization_within_magazine_entropy_magazine_summary.csv"),row.names=FALSE)
write.csv(std_entropy_out$overall,file.path(std_outdir,"standardization_within_magazine_entropy_overall.csv"),row.names=FALSE)

std_crossmag_js_yearly <- function(dl,n_bins=15,min_obs_year=5,min_mags_year=2,breaks=std_breaks){
  d <- std_analysis_data(dl)
  mags <- sort(unique(d$magazine_name))
  years <- sort(unique(d$year))
  comp_out <- list(); pair_out <- list(); time_out <- list()
  for(yr in years){
    probs_list <- list(); n_list <- c()
    for(mag in mags){
      dd <- d[d$magazine_name == mag & d$year == yr,]
      if(nrow(dd) >= min_obs_year){
        probs <- std_get_bin_probs(dd$std_value,breaks=breaks)
        if(!all(is.na(probs))){
          probs_list[[mag]] <- probs
          n_list[mag] <- nrow(dd)
          comp_out[[length(comp_out)+1]] <- data.frame(year=yr,magazine=mag,n_obs=nrow(dd),bin=seq_len(n_bins),bin_left=breaks[-length(breaks)],bin_right=breaks[-1],probability=probs)
        }
      }
    }
    if(length(probs_list) >= min_mags_year){
      mags_t <- names(probs_list)
      pairs <- combn(mags_t,2,simplify=FALSE)
      js_vals <- numeric(length(pairs))
      for(i in seq_along(pairs)){
        m1 <- pairs[[i]][1]; m2 <- pairs[[i]][2]
        js <- std_js_divergence(probs_list[[m1]],probs_list[[m2]])
        js_vals[i] <- js
        pair_out[[length(pair_out)+1]] <- data.frame(year=yr,magazine_1=m1,magazine_2=m2,n_1=n_list[m1],n_2=n_list[m2],js_divergence=js)
      }
      time_out[[length(time_out)+1]] <- data.frame(year=yr,n_magazines=length(probs_list),n_pairs=length(js_vals),mean_js=mean(js_vals,na.rm=TRUE),median_js=median(js_vals,na.rm=TRUE),sd_js=sd(js_vals,na.rm=TRUE),min_js=min(js_vals,na.rm=TRUE),max_js=max(js_vals,na.rm=TRUE))
    }
  }
  composition_long <- if(length(comp_out)>0) do.call(rbind,comp_out) else NULL
  pairwise_js <- if(length(pair_out)>0) do.call(rbind,pair_out) else NULL
  time_summary <- if(length(time_out)>0) do.call(rbind,time_out) else NULL
  overall <- NULL
  if(!is.null(time_summary) && nrow(time_summary) >= 2){
    fit_mean <- lm(mean_js ~ year,data=time_summary,weights=n_pairs)
    first <- time_summary[which.min(time_summary$year),]
    last <- time_summary[which.max(time_summary$year),]
    overall <- data.frame(first_year=first$year,last_year=last$year,first_mean_js=first$mean_js,last_mean_js=last$mean_js,mean_js_change=last$mean_js-first$mean_js,mean_js_slope_per_year=coef(fit_mean)[2])
  }
  list(composition_long=composition_long,pairwise_js=pairwise_js,time_summary=time_summary,overall=overall,breaks=breaks)
}

std_magazine_to_previous_magazines_pairwise_js <- function(dl,min_obs_year=5,breaks=std_breaks){
  d <- std_analysis_data(dl)
  years <- sort(unique(d$year))
  mags <- sort(unique(d$magazine_name))
  pair_out <- list()
  for(yr in years){
    prev_year_data <- d[d$year == yr - 1L,]
    if(nrow(prev_year_data) < min_obs_year) next
    prev_mags <- sort(unique(prev_year_data$magazine_name))
    for(mag_now in mags){
      dd_now <- d[d$magazine_name == mag_now & d$year == yr,]
      if(nrow(dd_now) < min_obs_year) next
      probs_now <- std_get_bin_probs(dd_now$std_value,breaks=breaks)
      if(all(is.na(probs_now))) next
      for(mag_prev in prev_mags){
        dd_prev <- prev_year_data[prev_year_data$magazine_name == mag_prev,]
        if(nrow(dd_prev) < min_obs_year) next
        probs_prev <- std_get_bin_probs(dd_prev$std_value,breaks=breaks)
        if(all(is.na(probs_prev))) next
        pair_out[[length(pair_out)+1]] <- data.frame(
          year=yr,
          magazine_current=mag_now,
          magazine_previous=mag_prev,
          n_current=nrow(dd_now),
          n_previous=nrow(dd_prev),
          js_divergence=std_js_divergence(probs_now,probs_prev)
        )
      }
    }
  }
  pairwise_prev_mag <- if(length(pair_out)>0) do.call(rbind,pair_out) else NULL
  if(is.null(pairwise_prev_mag)) return(list(pairwise_prev_mag=NULL,time_summary=NULL,overall=NULL))
  time_summary <- do.call(rbind,lapply(split(pairwise_prev_mag,pairwise_prev_mag$year),function(z){
    data.frame(
      year=z$year[1],
      n_pairs=nrow(z),
      n_magazines_current=length(unique(z$magazine_current)),
      n_magazines_previous=length(unique(z$magazine_previous)),
      mean_js=mean(z$js_divergence,na.rm=TRUE),
      median_js=median(z$js_divergence,na.rm=TRUE),
      sd_js=sd(z$js_divergence,na.rm=TRUE),
      min_js=min(z$js_divergence,na.rm=TRUE),
      max_js=max(z$js_divergence,na.rm=TRUE)
    )
  }))
  time_summary <- time_summary[order(time_summary$year),]
  overall <- NULL
  if(nrow(time_summary) >= 2){
    fit_mean <- lm(mean_js ~ year,data=time_summary,weights=n_pairs)
    first <- time_summary[which.min(time_summary$year),]
    last <- time_summary[which.max(time_summary$year),]
    overall <- data.frame(first_year=first$year,last_year=last$year,first_mean_js=first$mean_js,last_mean_js=last$mean_js,mean_js_change=last$mean_js-first$mean_js,mean_js_slope_per_year=coef(fit_mean)[2])
  }
  list(pairwise_prev_mag=pairwise_prev_mag,time_summary=time_summary,overall=overall)
}

std_cross_out <- std_crossmag_js_yearly(dl,n_bins=std_n_bins,min_obs_year=std_min_obs_year,min_mags_year=std_min_mags_year,breaks=std_breaks)
std_prev_pair_out <- std_magazine_to_previous_magazines_pairwise_js(dl,min_obs_year=std_min_obs_year,breaks=std_breaks)
write.csv(std_cross_out$composition_long,file.path(std_outdir,"standardization_crossmag_composition_bins.csv"),row.names=FALSE)
write.csv(std_cross_out$pairwise_js,file.path(std_outdir,"standardization_crossmag_pairwise_js.csv"),row.names=FALSE)
write.csv(std_cross_out$time_summary,file.path(std_outdir,"standardization_crossmag_yearly_js.csv"),row.names=FALSE)
write.csv(std_cross_out$overall,file.path(std_outdir,"standardization_crossmag_overall_js.csv"),row.names=FALSE)
write.csv(data.frame(bin=seq_len(length(std_cross_out$breaks)-1),bin_left=std_cross_out$breaks[-length(std_cross_out$breaks)],bin_right=std_cross_out$breaks[-1]),file.path(std_outdir,"standardization_crossmag_bin_breaks.csv"),row.names=FALSE)
write.csv(std_prev_pair_out$pairwise_prev_mag,file.path(std_outdir,"standardization_magazine_to_previous_year_magazines_pairwise_js.csv"),row.names=FALSE)
write.csv(std_prev_pair_out$time_summary,file.path(std_outdir,"standardization_previous_year_pairwise_tradition_disruption_js.csv"),row.names=FALSE)
write.csv(std_prev_pair_out$overall,file.path(std_outdir,"standardization_previous_year_pairwise_tradition_disruption_overall_js.csv"),row.names=FALSE)

std_pad_range <- function(x,pad=0.05,top_extra=0.10){
  x <- x[is.finite(x)]
  if(length(x)==0) return(c(0,1))
  r <- range(x,na.rm=TRUE)
  if(diff(r) == 0) rr <- r + c(-1,1)*max(0.01,abs(r[1])*pad) else rr <- r + c(-1,1)*diff(r)*pad
  rr[2] <- rr[2] + diff(rr)*top_extra
  rr
}

std_make_fit_ci <- function(df,y_col,weight_col=NULL){
  if(is.null(df) || nrow(df) < 2) return(NULL)
  tmp <- data.frame(year=df$year,y=df[[y_col]],w=if(is.null(weight_col)) rep(1,nrow(df)) else df[[weight_col]])
  fit <- lm(y ~ year,data=tmp,weights=w)
  x_new <- seq(min(tmp$year,na.rm=TRUE),max(tmp$year,na.rm=TRUE),length.out=200)
  pred <- predict(fit,newdata=data.frame(year=x_new),interval="confidence")
  list(fit=fit,x=x_new,pred=pred)
}

fit_prev_pair <- std_make_fit_ci(std_prev_pair_out$time_summary,"mean_js","n_pairs")
fit_entropy <- std_make_fit_ci(std_entropy_out$time_summary,"mean_entropy_norm","n_magazines")
fit_cross <- std_make_fit_ci(std_cross_out$time_summary,"mean_js","n_pairs")
std_y_values <- c(
  if(!is.null(std_cross_out$time_summary)) std_cross_out$time_summary$mean_js else NULL,
  if(!is.null(std_prev_pair_out$time_summary)) std_prev_pair_out$time_summary$mean_js else NULL,
  if(!is.null(std_entropy_out$time_summary)) std_entropy_out$time_summary$mean_entropy_norm else NULL,
  if(!is.null(fit_cross)) as.numeric(fit_cross$pred) else NULL,
  if(!is.null(fit_prev_pair)) as.numeric(fit_prev_pair$pred) else NULL,
  if(!is.null(fit_entropy)) as.numeric(fit_entropy$pred) else NULL
)
std_ylim <- std_pad_range(std_y_values)
col_cross <- "black"
col_prev <- "#D7301F"
col_entropy <- "#2166AC"

png(file.path(std_outdir,"standardization_crossmag_mean_js.png"),height=10,width=16,units="cm",res=600,type="cairo")
par(mar=c(3.5,3.5,2.5,1),mgp=c(2,0.65,0))
plot(NA,xlim=c(1900,1950),ylim=std_ylim,xlab="Year",ylab="Divergence / normalized entropy",
     main=paste0("Cross- and Within-Magazine Diversity per Year, ",plot_measure_title),
     xaxt="n",cex.main=0.92)
axis(1,at=std_year_ticks,labels=std_year_ticks,las=1)
abline(v=std_year_ticks,lty=2,lwd=0.7,col="grey80")

# Draw the previous-year pairwise divergence first, then entropy, then cross-magazine
# # diversity last, so the black line remains visually on top.
# if(!is.null(fit_prev_pair)){
#   polygon(c(fit_prev_pair$x,rev(fit_prev_pair$x)),c(fit_prev_pair$pred[,"lwr"],rev(fit_prev_pair$pred[,"upr"])),col=adjustcolor(col_prev,alpha.f=0.10),border=NA)
#   lines(fit_prev_pair$x,fit_prev_pair$pred[,"fit"],col=col_prev,lwd=2,lty=2)
# }
if(!is.null(fit_entropy)){
  polygon(c(fit_entropy$x,rev(fit_entropy$x)),c(fit_entropy$pred[,"lwr"],rev(fit_entropy$pred[,"upr"])),col=adjustcolor(col_entropy,alpha.f=0.10),border=NA)
  lines(fit_entropy$x,fit_entropy$pred[,"fit"],col=col_entropy,lwd=2,lty=2)
}
if(!is.null(fit_cross)){
  polygon(c(fit_cross$x,rev(fit_cross$x)),c(fit_cross$pred[,"lwr"],rev(fit_cross$pred[,"upr"])),col=adjustcolor(col_cross,alpha.f=0.10),border=NA)
  lines(fit_cross$x,fit_cross$pred[,"fit"],col=col_cross,lwd=2,lty=2)
}
# if(!is.null(std_prev_pair_out$time_summary) && nrow(std_prev_pair_out$time_summary) > 0) lines(std_prev_pair_out$time_summary$year,std_prev_pair_out$time_summary$mean_js,type="b",pch=17,lwd=1.2,col=col_prev)
if(!is.null(std_entropy_out$time_summary) && nrow(std_entropy_out$time_summary) > 0) lines(std_entropy_out$time_summary$year,std_entropy_out$time_summary$mean_entropy_norm,type="b",pch=15,lwd=1.2,col=col_entropy)
if(!is.null(std_cross_out$time_summary) && nrow(std_cross_out$time_summary) > 0) lines(std_cross_out$time_summary$year,std_cross_out$time_summary$mean_js,type="b",pch=16,lwd=1.2,col=col_cross)
legend("topright",legend=c("Cross-magazine diversity","Within-magazine entropy"),col=c(col_cross,col_entropy),lty=1,pch=c(16,15),lwd=1.4,bty="n",cex=0.80)
dev.off()

if(!is.null(std_entropy_out$magazine_summary) && nrow(std_entropy_out$magazine_summary) > 0){
  png(file.path(std_outdir,"standardization_within_magazine_entropy_slope.png"),height=12,width=16,units="cm",res=600,type="cairo")
  par(mar=c(11.5,5.1,4,1)+0.1, mgp=c(3.0,0.8,0))
  ord <- order(std_entropy_out$magazine_summary$entropy_slope_per_year)
  barplot(std_entropy_out$magazine_summary$entropy_slope_per_year[ord],
          names.arg=std_entropy_out$magazine_summary$magazine[ord],
          las=2,cex.names=0.55,
          ylab="Mean normalized entropy slope per year",
          main="Within-magazine entropy, ordered by entropy slope")
  abline(h=0,lty=2)
  dev.off()
}

# Print quick summaries
std_entropy_out$overall
std_cross_out$overall
std_prev_pair_out$overall


cor_data <- Reduce(function(x,y) merge(x,y,by="year",all=FALSE), list(
  std_entropy_out$time_summary[,c("year","mean_entropy_norm")],
  std_cross_out$time_summary[,c("year","mean_js")],
  std_prev_pair_out$time_summary[,c("year","mean_js")]
))

names(cor_data) <- c("year","entropy","cross_js","prev_year_js")

cor(cor_data$entropy, cor_data$cross_js, use="complete.obs")
cor(cor_data$cross_js, cor_data$prev_year_js, use="complete.obs")
