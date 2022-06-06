majiq <- readRDS("~/projects/AS/analysis/paper_gb/sim_1/MAJIQ/majiq2.rds")
majiq <- majiq |>
  mutate(
    as_type = dplyr::case_when(
      as.numeric(`Num. Junctions`) == 1 & ES == "False" & A5SS == "False" & A3SS == "False" & !is.na(`ir_coords`) ~ "RI",
      as.numeric(`Num. Junctions`) == 2 & as.numeric(`Num. Exons`) == 3 ~ "SE",
      as.numeric(`Num. Junctions`) == 4 & as.numeric(`Num. Exons`) == 4 ~ "MXE",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A3SS == "True" ~ "A3SS",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A5SS == "True" ~ "A5SS",
      T ~ "OTHER"
    )
  )
cols_coord <- c(
  "exons_coords", "junctions_coords", "ir_coords"
)
# majiq_e <- markEvent(
#   majiq, truths, cols_coord, split = T, split_by = "[;-]",
#   numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
#   error_distance = 1, as = T, gene = T, gene_symbol = T,
#   truths_col = c("X16", "X17"), mark_col = "idx_truths_as"
# )
saveRDS(majiq, "~/projects/AS/analysis/paper_gb/sim_1/MAJIQ/majiq2.rds")


majiq <- readRDS("~/projects/AS/analysis/paper_gb/sim_1/MAJIQ/majiq.rds")
majiq <- majiq |>
  mutate(
    as_type = dplyr::case_when(
      as.numeric(`Num. Junctions`) == 1 & ES == "False" & A5SS == "False" & A3SS == "False" & !is.na(`ir_coords`) ~ "RI",
      as.numeric(`Num. Junctions`) == 2 & as.numeric(`Num. Exons`) == 3 ~ "SE",
      as.numeric(`Num. Junctions`) == 4 & as.numeric(`Num. Exons`) == 4 ~ "MXE",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A3SS == "True" ~ "A3SS",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A5SS == "True" ~ "A5SS",
      T ~ "OTHER"
    )
  )
cols_coord <- c(
  "exons_coords", "junctions_coords", "ir_coords"
)
# majiq_e <- markEvent(
#   majiq, truths, cols_coord, split = T, split_by = "[;-]",
#   numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
#   error_distance = 1, as = T, gene = T, gene_symbol = T,
#   truths_col = c("X16", "X17"), mark_col = "idx_truths_as"
# )
saveRDS(majiq, "~/projects/AS/analysis/paper_gb/sim_1/MAJIQ/majiq.rds")




majiq <- readRDS("~/projects/AS/analysis/novel_junction/MAJIQ/majiq2.rds")
majiq <- majiq |>
  mutate(
    as_type = dplyr::case_when(
      as.numeric(`Num. Junctions`) == 1 & ES == "False" & A5SS == "False" & A3SS == "False" & !is.na(`ir_coords`) ~ "RI",
      as.numeric(`Num. Junctions`) == 2 & as.numeric(`Num. Exons`) == 3 ~ "SE",
      as.numeric(`Num. Junctions`) == 4 & as.numeric(`Num. Exons`) == 4 ~ "MXE",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A3SS == "True" ~ "A3SS",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A5SS == "True" ~ "A5SS",
      T ~ "OTHER"
    )
  )
# cols_coord <- c(
#   "exons_coords", "junctions_coords", "ir_coords"
# )
# majiq_e <- markEvent(
#   majiq, truths, cols_coord, split = T, split_by = "[;-]",
#   numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
#   error_distance = 1, as = T, gene = T, gene_symbol = T,
#   truths_col = c("X16", "X17"), mark_col = "idx_truths_as"
# )
saveRDS(majiq, "~/projects/AS/analysis/novel_junction/MAJIQ/majiq2.rds")


majiq <- readRDS("~/projects/AS/analysis/novel_junction/MAJIQ/majiq.rds")
majiq <- majiq |>
  mutate(
    as_type = dplyr::case_when(
      as.numeric(`Num. Junctions`) == 1 & ES == "False" & A5SS == "False" & A3SS == "False" & !is.na(`ir_coords`) ~ "RI",
      as.numeric(`Num. Junctions`) == 2 & as.numeric(`Num. Exons`) == 3 ~ "SE",
      as.numeric(`Num. Junctions`) == 4 & as.numeric(`Num. Exons`) == 4 ~ "MXE",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A3SS == "True" ~ "A3SS",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A5SS == "True" ~ "A5SS",
      T ~ "OTHER"
    )
  )
# cols_coord <- c(
#   "exons_coords", "junctions_coords", "ir_coords"
# )
# majiq_e <- markEvent(
#   majiq, truths, cols_coord, split = T, split_by = "[;-]",
#   numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
#   error_distance = 1, as = T, gene = T, gene_symbol = T,
#   truths_col = c("X16", "X17"), mark_col = "idx_truths_as"
# )
saveRDS(majiq, "~/projects/AS/analysis/novel_junction/MAJIQ/majiq.rds")




majiq <- readRDS("~/projects/AS/analysis/novel_splice_site/MAJIQ/majiq2.rds")
majiq <- majiq |>
  mutate(
    as_type = dplyr::case_when(
      as.numeric(`Num. Junctions`) == 1 & ES == "False" & A5SS == "False" & A3SS == "False" & !is.na(`ir_coords`) ~ "RI",
      as.numeric(`Num. Junctions`) == 2 & as.numeric(`Num. Exons`) == 3 ~ "SE",
      as.numeric(`Num. Junctions`) == 4 & as.numeric(`Num. Exons`) == 4 ~ "MXE",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A3SS == "True" ~ "A3SS",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A5SS == "True" ~ "A5SS",
      T ~ "OTHER"
    )
  )
# cols_coord <- c(
#   "exons_coords", "junctions_coords", "ir_coords"
# )
# majiq_e <- markEvent(
#   majiq, truths, cols_coord, split = T, split_by = "[;-]",
#   numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
#   error_distance = 1, as = T, gene = T, gene_symbol = T,
#   truths_col = c("X16", "X17"), mark_col = "idx_truths_as"
# )
saveRDS(majiq, "~/projects/AS/analysis/novel_splice_site/MAJIQ/majiq2.rds")


majiq <- readRDS("~/projects/AS/analysis/novel_splice_site/MAJIQ/majiq.rds")
majiq <- majiq |>
  mutate(
    as_type = dplyr::case_when(
      as.numeric(`Num. Junctions`) == 1 & ES == "False" & A5SS == "False" & A3SS == "False" & !is.na(`ir_coords`) ~ "RI",
      as.numeric(`Num. Junctions`) == 2 & as.numeric(`Num. Exons`) == 3 ~ "SE",
      as.numeric(`Num. Junctions`) == 4 & as.numeric(`Num. Exons`) == 4 ~ "MXE",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A3SS == "True" ~ "A3SS",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A5SS == "True" ~ "A5SS",
      T ~ "OTHER"
    )
  )
# cols_coord <- c(
#   "exons_coords", "junctions_coords", "ir_coords"
# )
# majiq_e <- markEvent(
#   majiq, truths, cols_coord, split = T, split_by = "[;-]",
#   numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
#   error_distance = 1, as = T, gene = T, gene_symbol = T,
#   truths_col = c("X16", "X17"), mark_col = "idx_truths_as"
# )
saveRDS(majiq, "~/projects/AS/analysis/novel_splice_site/MAJIQ/majiq.rds")




majiq <- readRDS("~/projects/AS/analysis/filter/MAJIQ/majiq2.rds")
majiq <- majiq |>
  mutate(
    as_type = dplyr::case_when(
      as.numeric(`Num. Junctions`) == 1 & ES == "False" & A5SS == "False" & A3SS == "False" & !is.na(`ir_coords`) ~ "RI",
      as.numeric(`Num. Junctions`) == 2 & as.numeric(`Num. Exons`) == 3 ~ "SE",
      as.numeric(`Num. Junctions`) == 4 & as.numeric(`Num. Exons`) == 4 ~ "MXE",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A3SS == "True" ~ "A3SS",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A5SS == "True" ~ "A5SS",
      T ~ "OTHER"
    )
  )
# cols_coord <- c(
#   "exons_coords", "junctions_coords", "ir_coords"
# )
# majiq_e <- markEvent(
#   majiq, truths, cols_coord, split = T, split_by = "[;-]",
#   numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
#   error_distance = 1, as = T, gene = T, gene_symbol = T,
#   truths_col = c("X16", "X17"), mark_col = "idx_truths_as"
# )
saveRDS(majiq, "~/projects/AS/analysis/filter/MAJIQ/majiq2.rds")


majiq <- readRDS("~/projects/AS/analysis/filter/MAJIQ/majiq.rds")
majiq <- majiq |>
  mutate(
    as_type = dplyr::case_when(
      as.numeric(`Num. Junctions`) == 1 & ES == "False" & A5SS == "False" & A3SS == "False" & !is.na(`ir_coords`) ~ "RI",
      as.numeric(`Num. Junctions`) == 2 & as.numeric(`Num. Exons`) == 3 ~ "SE",
      as.numeric(`Num. Junctions`) == 4 & as.numeric(`Num. Exons`) == 4 ~ "MXE",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A3SS == "True" ~ "A3SS",
      as.numeric(`Num. Junctions`) %in% c(2, 3) & A5SS == "True" ~ "A5SS",
      T ~ "OTHER"
    )
  )
# cols_coord <- c(
#   "exons_coords", "junctions_coords", "ir_coords"
# )
# majiq_e <- markEvent(
#   majiq, truths, cols_coord, split = T, split_by = "[;-]",
#   numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
#   error_distance = 1, as = T, gene = T, gene_symbol = T,
#   truths_col = c("X16", "X17"), mark_col = "idx_truths_as"
# )
saveRDS(majiq, "~/projects/AS/analysis/filter/MAJIQ/majiq.rds")








