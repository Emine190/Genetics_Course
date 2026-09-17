## ============================================================================
## Manual cTWAS-style locus fine-mapping — multi-gene + empirical prior pass
##
## Builds on a susie_rss()-based script. Two additions:
##
##  (A) run_ctwas_locus_region() now takes a VECTOR of candidate genes per
##      locus and fine-maps all of them jointly against the region's SNPs
##      in one susie_rss() call (one pseudo-SNP column per gene, all sharing
##      the same LD-derived correlation structure).
##
##  (B) A two-pass calibration: Pass 1 runs every locus with SuSiE's default
##      uniform prior_weights and records, on average, how much posterior
##      mass ends up on "Gene" rows vs "SNP" rows. Pass 2 re-runs every locus
##      using those empirical group-level probabilities as prior_weights.
##      This is NOT the same as cTWAS's own EM (which also estimates a
##      separate effect-size variance per group, and iterates within each
##      region rather than once across all regions) — but it's a reasonable,
##      transparent approximation of "let the data tell us whether genes or
##      SNPs deserve more prior credibility" without needing the ctwas
##      package's genome-wide machinery.
## ============================================================================

#install.packages(data.table)
#install.packages(susieR)
#install.packages(remotes)

#remotes::install_github("XingHua-Lab/ctwas", build_vignettes = FALSE)

library(data.table)
library(susieR)
library(ctwas)
library(Matrix)

setwd("~/Desktop/Feral_project")

## ----------------------------------------------------------------------
## Load shared inputs once
## ----------------------------------------------------------------------

expr_dt   <- fread("your_file_name_forexpression.txt")
pheno_dt  <- fread("your_pheno_file.txt")
geno_dt   <- fread("genotype_file.txt")
gemma_gwas <- fread("gemma_gwas_filehere.txt")
#e_qtl <- fread("eqtlfilehere.txt")

prep_gemma_zscores <- function(gemma_dt) {
  gwas_dt <- copy(as.data.table(gemma_dt))
  gwas_dt[, z := Beta / SE]
  gwas_dt[, .(id = SNP_ID, Trait = Trait, Beta = Beta, SE = SE, z = z)][!is.na(z)]
}

## ----------------------------------------------------------------------
## (A) Multi-gene joint locus fine-mapping
##     target_genes is now a character VECTOR, not a single string.
##     prior_weights, if supplied, must be a named vector: names should be
##     "SNP" and "Gene" (used as the group-level prior probability, applied
##     to every row of that type at this locus). If NULL, SuSiE's uniform
##     default is used (equivalent to Pass 1 / uncalibrated).
## ----------------------------------------------------------------------

run_ctwas_locus_region <- function(target_trait, target_genes,
                                    chr, start_pos, end_pos,
                                    gemma_dt, pheno_dt, geno_dt, expr_dt,
                                    n_gwas_samples = 500,
                                    group_prior_weights = NULL) {

  message(sprintf("--> Fine-mapping Locus [Chr %s: %s - %s] for %d gene(s): %s",
                   chr, start_pos, end_pos, length(target_genes),
                   paste(target_genes, collapse = ", ")))

  ## 1. GWAS Z-scores for this trait
  z_all <- prep_gemma_zscores(gemma_dt)
  z_sub <- z_all[Trait == target_trait]

  ## 2. SNPs in the physical window
  locus_geno_dt <- geno_dt[grepl(paste0("^0?", chr, ":"), SNP)]
  locus_geno_dt[, pos := as.numeric(tstrsplit(SNP, ":")[[2]])]
  locus_geno_dt <- locus_geno_dt[pos >= start_pos & pos <= end_pos]
  if (nrow(locus_geno_dt) < 2) {
    message("Not enough SNPs found in physical window.")
    return(NULL)
  }

  sample_cols_expr <- setdiff(names(expr_dt), "ID")
  sample_cols_geno <- setdiff(names(geno_dt), c("SNP", "pos"))
  samples_clean_expr <- gsub("[^0-9]", "", sample_cols_expr)
  samples_clean_geno <- gsub("[^0-9]", "", sample_cols_geno)
  common_samples <- intersect(samples_clean_expr, samples_clean_geno)
  idx_expr <- match(common_samples, samples_clean_expr)
  idx_geno <- match(common_samples, samples_clean_geno)

  X_mat <- t(as.matrix(locus_geno_dt[, ..sample_cols_geno]))[idx_geno, , drop = FALSE]
  colnames(X_mat) <- locus_geno_dt$SNP

  var_snps_pre <- apply(X_mat, 2, var, na.rm = TRUE)
  X_mat <- X_mat[, var_snps_pre > 0, drop = FALSE]
  if (ncol(X_mat) < 2) return(NULL)

  ## 3. Region-wide LD matrix (shared across every gene at this locus)
  ld_snps <- cor(X_mat, use = "pairwise.complete.obs")

  ## 4. Build ONE pseudo-SNP (gene) column per candidate gene
  gene_entries <- list()

  for (target_gene in target_genes) {
    gene_row <- expr_dt[ID %in% target_gene]
    if (nrow(gene_row) == 0) {
      message(sprintf("  [skip] %s not found in expression matrix.", target_gene))
      next
    }

    y_vec <- as.numeric(gene_row[1, ..sample_cols_expr])[idx_expr]
    valid_rows <- complete.cases(y_vec, X_mat)
    X_sub <- X_mat[valid_rows, , drop = FALSE]
    y_sub <- y_vec[valid_rows]

    var_snps <- apply(X_sub, 2, var, na.rm = TRUE)
    X_sub <- X_sub[, var_snps > 0, drop = FALSE]
    if (ncol(X_sub) < 2) {
      message(sprintf("  [skip] %s: no variable SNPs after filtering.", target_gene))
      next
    }

    weights <- apply(X_sub, 2, function(snp_col) coef(lm(y_sub ~ snp_col))[2])
    weights[is.na(weights)] <- 0

    locus_snps <- colnames(X_sub)
    z_snps_df <- z_sub[id %in% locus_snps]
    common_snps <- intersect(locus_snps, z_snps_df$id)
    if (length(common_snps) < 2) {
      message(sprintf("  [skip] %s: fewer than 2 SNPs overlap GWAS z-scores.", target_gene))
      next
    }

    ld_g <- ld_snps[common_snps, common_snps, drop = FALSE]
    w_vec <- weights[common_snps]
    w_denom <- sqrt(as.numeric(t(w_vec) %*% ld_g %*% w_vec))
    if (is.na(w_denom) || w_denom == 0) {
      message(sprintf("  [skip] %s: degenerate weight normalization.", target_gene))
      next
    }

    z_snps_vec <- z_snps_df[match(common_snps, id), z]
    names(z_snps_vec) <- common_snps
    gene_z <- sum(w_vec * z_snps_vec) / w_denom
    gene_ld_col <- as.vector(ld_g %*% w_vec) / w_denom

    gene_entries[[target_gene]] <- list(
      common_snps = common_snps, gene_z = gene_z, gene_ld_col = gene_ld_col,
      w_vec = w_vec, ld_g = ld_g
    )
  }

  if (length(gene_entries) == 0) {
    message("No candidate genes produced a usable weight model at this locus.")
    return(NULL)
  }

  ## 5. Assemble ONE joint z-vector / LD matrix spanning every SNP that
  ##    contributed to ANY gene's weight model, plus one row/col per gene.
  all_snps <- Reduce(union, lapply(gene_entries, function(e) e$common_snps))
  base_ld <- ld_snps[all_snps, all_snps, drop = FALSE]
  base_z  <- z_sub[match(all_snps, id), z]
  names(base_z) <- all_snps

  n_genes <- length(gene_entries)
  ld_joint <- rbind(
    cbind(base_ld, matrix(0, nrow = length(all_snps), ncol = n_genes)),
    cbind(matrix(0, nrow = n_genes, ncol = length(all_snps)),
          diag(n_genes))
  )
  gene_names <- names(gene_entries)
  colnames(ld_joint) <- rownames(ld_joint) <- c(all_snps, gene_names)

  for (g in gene_names) {
    e <- gene_entries[[g]]
    # gene-SNP correlation, placed only in the columns this gene's model used;
    # SNPs outside this gene's own model default to 0 correlation with it
    # (a simplifying approximation — see caveat below).
    col_vec <- setNames(rep(0, length(all_snps)), all_snps)
    col_vec[e$common_snps] <- e$gene_ld_col
    ld_joint[all_snps, g] <- col_vec
    ld_joint[g, all_snps] <- col_vec
  }

  z_joint <- c(base_z, setNames(sapply(gene_entries, function(e) e$gene_z), gene_names))
  type_vec <- setNames(c(rep("SNP", length(all_snps)), rep("Gene", n_genes)), names(z_joint))

  ## 6. Prior weights: uniform (Pass 1) or empirically calibrated (Pass 2)
  if (is.null(group_prior_weights)) {
    prior_weights <- NULL   # susie_rss default: uniform
  } else {
    prior_weights <- ifelse(type_vec == "Gene",
                             group_prior_weights["Gene"],
                             group_prior_weights["SNP"])
    names(prior_weights) <- names(z_joint)
  }

  ## 7. Fine-map
  susie_fit <- susieR::susie_rss(
    z = z_joint,
    R = ld_joint,
    n = n_gwas_samples,
    L = min(5, length(z_joint)),
    prior_weights = prior_weights,
    check_R = FALSE
  )

  res_dt <- data.table(
    Locus_Chr = chr,
    Locus_Start = start_pos,
    Locus_End = end_pos,
    Trait = target_trait,
    ID = names(z_joint),
    Type = type_vec[names(z_joint)],
    PIP = susie_fit$pip
  )

  res_dt[order(-PIP)]
}

## ----------------------------------------------------------------------
## Define your loci + candidate genes here.
## One row per locus; genes is a list-column so a locus can carry >1 gene.BERKAU:34:1632693-1642694, BER:2:3273757-3283757, KAU:1:166050552-166060552
## ----------------------------------------------------------------------

loci_table <- list(
  list(trait = "MEANXOTHER_D", chr = "1",  start = 166050552, end = 166060552,
       genes = c("ENSGALG00010003246"), n = 376)
  # add more loci here, e.g.:
  # list(trait = "OTHER_TRAIT", chr = "3", start = 1234000, end = 1244000,
  #      genes = c("ENSGALG..."), n = 66)
)

## ----------------------------------------------------------------------
## (B) Pass 1 — uncalibrated, uniform priors — run every locus once
## ----------------------------------------------------------------------

pass1_results <- rbindlist(lapply(loci_table, function(L) {
  run_ctwas_locus_region(
    target_trait = L$trait, target_genes = L$genes,
    chr = L$chr, start_pos = L$start, end_pos = L$end,
    gemma_dt = gemma_gwas, pheno_dt = pheno_dt, geno_dt = geno_dt, expr_dt = expr_dt,
    n_gwas_samples = L$n, group_prior_weights = NULL
  )
}), fill = TRUE)

## Empirical group-level prior: mean PIP achieved by each type across all
## loci in Pass 1. This is the stand-in for cTWAS's own estimated
## group_prior — a rough measure of "how often does this type of entry end
## up looking causal, on average, across the loci I actually have."
empirical_prior <- pass1_results[, .(mean_pip = mean(PIP)), by = Type]
cat("\nPass 1 empirical group priors:\n")
print(empirical_prior)

group_prior_weights <- setNames(empirical_prior$mean_pip, empirical_prior$Type)
# guard against a zero/NA prior collapsing a group entirely
group_prior_weights[is.na(group_prior_weights) | group_prior_weights <= 0] <- 0.01

## ----------------------------------------------------------------------
## Pass 2 — re-run every locus using the empirically calibrated priors
## ----------------------------------------------------------------------

pass2_results <- rbindlist(lapply(loci_table, function(L) {
  run_ctwas_locus_region(
    target_trait = L$trait, target_genes = L$genes,
    chr = L$chr, start_pos = L$start, end_pos = L$end,
    gemma_dt = gemma_gwas, pheno_dt = pheno_dt, geno_dt = geno_dt, expr_dt = expr_dt,
    n_gwas_samples = L$n, group_prior_weights = group_prior_weights
  )
}), fill = TRUE)

cat("\nPass 2 (calibrated) top hits:\n")
print(head(pass2_results[order(-PIP)], 15))

fwrite(pass1_results, "ctwas_manual_pass1_uncalibrated.txt", sep = "\t")
fwrite(pass2_results, "ctwas_manual_pass2_calibrated.txt", sep = "\t")

## ----------------------------------------------------------------------
## Caveats to keep in mind / mention in your methods section
## ----------------------------------------------------------------------
## 1. Gene-gene correlation within a locus is NOT modeled (off-diagonal
##    between two genes' pseudo-SNP rows is left at 0). If two candidate
##    genes at the same locus share overlapping weight SNPs, SuSiE won't
##    "know" they're correlated with each other — only with the SNPs.
##    For most 2-3-gene loci with mostly distinct eQTL SNPs this is a minor
##    approximation; flag it if your genes share most of their weight SNPs.
## 2. This calibration uses ONE shared prior per group (SNP vs Gene) across
##    ALL loci, estimated from the same small set of loci being tested —
##    real cTWAS estimates priors from many more regions than you're
##    fine-mapping, and iterates (niter1/niter2) rather than doing a single
##    pass. With only a handful of loci, treat Pass 2's calibration as a
##    rough sensitivity check against Pass 1, not a rigorously estimated
##    prior — the more loci you add to loci_table, the more defensible this
##    empirical step becomes.
## 3. No separate effect-size-variance term per group (group_prior_var in
##    real cTWAS) is estimated here — only inclusion probability.
