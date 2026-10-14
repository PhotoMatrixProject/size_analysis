library(magick)


export_image_paths <- function(root = ".", out_csv = "image_files_inventory.csv",
                               out_txt = "image_paths_for_R.txt",
                               extensions = c("png", "jpg", "jpeg", "tif", "tiff", "webp")) {
  
  pattern <- paste0("\\.(", paste(extensions, collapse = "|"), ")$")
  files <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  files <- sort(files)
  
  inventory <- data.frame(folder = dirname(files),
                          file = basename(files),
                          path = files,
                          r_path = paste0('"', files, '"'),
                          stringsAsFactors = FALSE)
  
  write.csv(inventory, out_csv, row.names = FALSE)
  writeLines(paste0("c(\n  ", paste(inventory$r_path, collapse = ",\n  "), "\n)"), out_txt)
  
  message("Found ", length(files), " image files.")
  message("Wrote: ", out_csv)
  message("Wrote: ", out_txt)
  
  invisible(inventory)
}


make_composite <- function(files, output, labels = LETTERS[seq_along(files)], ncol = NULL, nrow = NULL, byrow = TRUE,
                           label_size = 48, label_gravity = "northwest", label_location = "+35+25",
                           label_color = "black", label_weight = 700, scale_to = c("height", "width", "none"),
                           bg = "white") {
  
  scale_to <- match.arg(scale_to)
  stopifnot(length(files) == length(labels))
  
  n <- length(files)
  if (is.null(ncol) & is.null(nrow)) ncol <- n
  if (is.null(ncol)) ncol <- ceiling(n / nrow)
  if (is.null(nrow)) nrow <- ceiling(n / ncol)
  
  imgs <- Map(function(file, label) {
    img <- image_read(file)
    image_annotate(img, text = label, size = label_size, gravity = label_gravity,
                   location = label_location, color = label_color, weight = label_weight)
  }, files, labels)
  
  if (scale_to == "height") {
    target_h <- max(sapply(imgs, function(x) image_info(x)$height))
    imgs <- lapply(imgs, function(x) image_scale(x, paste0("x", target_h)))
  }
  
  if (scale_to == "width") {
    target_w <- max(sapply(imgs, function(x) image_info(x)$width))
    imgs <- lapply(imgs, function(x) image_scale(x, paste0(target_w)))
  }
  
  idx <- matrix(c(seq_len(n), rep(NA, nrow * ncol - n)), nrow = nrow, ncol = ncol, byrow = byrow)
  info <- image_info(imgs[[1]])
  blank <- image_blank(info$width, info$height, color = bg)
  
  row_imgs <- lapply(seq_len(nrow), function(i) {
    row <- lapply(idx[i, ], function(j) if (is.na(j)) blank else imgs[[j]])
    image_append(do.call(c, row), stack = FALSE)
  })
  
  max_w <- max(sapply(row_imgs, function(x) image_info(x)$width))
  row_imgs <- lapply(row_imgs, function(x) image_extent(x, geometry = paste0(max_w, "x", image_info(x)$height),
                                                        gravity = "west", color = bg))
  
  combined <- image_append(do.call(c, row_imgs), stack = TRUE)
  image_write(combined, output)
  invisible(combined)
}


imgs <- export_image_paths()

#Variances
make_composite(files = c("./pipeline_area/time_results/Var_Czech_logit_normal_step.png",
                         "./pipeline_area/time_results/Var_Czech_logit_normal_no_time.png"),
               labels = c("A", "B"),
               output = "Var_Czech_rel_width_logit_normal_spline_vs_no_time.png",label_size = 150, label_location = "+150+35")

#Flows
make_composite(files = c("./pipeline_area/time_results/CzechStep_logit_normal1.png",
                         "./pipeline_width_rel/time_results/CzechStep_rel_width_logit_normal1.png",
                         "./pipeline_widthabs/time_results/CzechStep_width_normal1.png"),
               labels = c("A", "B","C"),
               output = "Figure1_Flows.png",label_size = 150, label_location = "+150+35", ncol = 1)


#Divergences
make_composite(files = c("./pipeline_area/standardization_area_yearly_bins15/standardization_crossmag_mean_js.png",
                         "./pipeline_width_rel/standardization_rel_width_yearly_bins15/standardization_crossmag_mean_js.png",
                         "./pipeline_widthabs/standardization_width_yearly_bins15/standardization_crossmag_mean_js.png"),
               labels = c("A", "B","C"),
               output = "Figure3_Divergence.png",label_size = 150, label_location = "+150+35", ncol = 1)

#Particular_panels
make_composite(files = c("./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Topicuv_sbornik_literarni_a_umelecky_temporal_rel_width_modes.png",
                         "./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Drobne_umeni_temporal_rel_width_modes.png",
                         "./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Vytvarna_prace_temporal_rel_width_modes.png",
                         "./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Vytvarne_snahy_temporal_rel_width_modes.png"),
               labels = c("A", "B","C","D"),
               output = "Figure4_Panels1.png",label_size = 150, label_location = "+150+35", ncol = 2,byrow=T)


make_composite(files = c("./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Fotograficky_obzor_temporal_rel_width_modes.png",
                         "./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Ceskoslovenska_fotografie_temporal_rel_width_modes.png",
                         "./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Umeni_temporal_rel_width_modes.png",
                         "./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Vytvarna_vychova_temporal_rel_width_modes.png",
                         "./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Umeni_dneska_temporal_rel_width_modes.png",
                         "./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Architektura_temporal_rel_width_modes.png"),
               labels = c("A", "B","C","D","E","F"),
               output = "Figure5_Panels2.png",label_size = 150, label_location = "+150+35", ncol = 2,byrow=T)


# make_composite(files = c("./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Umeni_dneska_temporal_rel_width_modes.png",
#                          "./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Typografia_temporal_rel_width_modes.png",
#                "./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Architektura_temporal_rel_width_modes.png",
#                "./pipeline_width_rel/temporal_rel_width_figures_w3_s0p5_adj0p7_frac0p1_jump0p5_min3/Blok_temporal_rel_width_modes.png"),
#                labels = c("A", "B","C","D"),
#                output = "Figure6_Panels3.png",label_size = 150, label_location = "+150+35", ncol = 2,byrow=T)
