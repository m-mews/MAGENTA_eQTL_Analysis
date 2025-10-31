suppressMessages(library(data.table))
suppressMessages(library(magrittr))
suppressMessages(library(susieR))
suppressMessages(library(Cairo))
suppressMessages(library(tidyverse))

args = commandArgs(trailing = TRUE)
anc = args[1]
hgnc = args[2]
gene = args[3]
chr = args[4]
eqtl_n = args[5]
plink_ld_all = args[6]
base = args[7]

## Read eQTL file 

ss_name = paste0("/mnt/vstor/Data14/metabrain_lasso/magenta_susie/ieGenes/",anc,"/", hgnc, "_", gene, "_chr_", chr, "_", anc, "/MAGENTA_ieQTL_", gene, "_", hgnc, "_chr_", chr, "_", anc, ".txt.gz")

# Load in Summary statistics
ss_all <- fread(ss_name) %>%
        # ADD Z-score
        mutate(z = beta/se, z_ad = beta_ad/se_ad, z_ad_int = beta_ad_interaction/se_ad_interaction) %>%
	# Filter by Valid SNPs
        filter(!is.na(z_ad_int))



print('Total Variants used for SUSIE')
tvar_susie <- nrow(ss_all)
tvar_susie


# Derive Sample Size for SUSIE
eqtl_n_old <- ss_all$n[1]

print('Sample Size for SUSIE')
eqtl_n = as.numeric(eqtl_n)
eqtl_n

print('Sample Size OLD')
eqtl_n_old

# Load in PLINK LD ALL
ld_all <- read_tsv(plink_ld_all, col_names = FALSE)
ld_all_mat <- as.matrix(ld_all)

colnames(ld_all_mat) <- rownames(ld_all_mat)

print('Sanity check PLINK LD ALL should match total variants used for SUSIE')
tvar_plink_all <- nrow(ld_all)
tvar_plink_all
ifelse(tvar_susie == tvar_plink_all, TRUE, FALSE)

#################################################################################################
## Run SUSIE RSS - ALL VARIANTS
fitted_rss_eqtl_all <- susie_rss(ss_all$z_ad_int, ld_all_mat, n = eqtl_n, L = 10, estimate_residual_variance = TRUE, max_iter = 100)


print('SUSIE RSS Model Fit - ALL Variants')

###############################################################################################
# Write out summary statistics with PIP & convergence status included

ss_all$pip <- fitted_rss_eqtl_all$pip

###########################################################################################################

################################################################################################
## Save CS
cs_all <- if (!is.null(summary(fitted_rss_eqtl_all)$cs)) {
    summary(fitted_rss_eqtl_all)$cs

} else {
    data.frame(
        cs = NA,
        cs_log10bf = NA,
        cs_avg_r2 = NA,
        cs_min_r2 = NA,
        variable = NA
    )
}
cs_all$gene = gene
cs_all$hgnc = hgnc
cs_all$chr = chr
cs_all$anc = anc
cs_all$convergence = fitted_rss_eqtl_all$converged

### Save file

cs_all_file_name = paste0("/mnt/vstor/Data14/metabrain_lasso/magenta_susie/ieGenes/",anc , "/", hgnc, "_", gene, "_chr_", chr, "_", anc, "/MAGENTA_ieQTL_", gene, "_", hgnc, "_chr_", chr, "_", anc, "_susie_cs_corrected_n_est_resid_var.txt.gz")

fwrite(x = cs_all, file = cs_all_file_name, sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")

print('Credible Sets Saved All')
#################################################################################################
## Extract CS SNPs for OTTERs Rerun

if (!is.null(summary(fitted_rss_eqtl_all)$cs)) {
	print("Extracting Credible Sets")
	cs_snps <- unlist(strsplit(as.character(cs_all$variable), ","))
	ss_all$row <- rownames(ss_all)
	ss_all <- ss_all %>% mutate(CS = ifelse(row %in% cs_snps, 1, 0))
	cs_all_file_name = paste0("/mnt/vstor/Data14/metabrain_lasso/magenta_susie/ieGenes/",anc, "/", hgnc, "_", gene, "_chr_", chr, "_", anc, "/MAGENTA_ieQTL_", gene, "_", hgnc, "_chr_", chr, "_", anc, "_susie_pip_w_cs_corrected_n_est_resid_var.txt.gz")
	fwrite(x = ss_all, file = cs_all_file_name, sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")

	# Generate Credible Model Plot

	no_anno_title = paste0(anc , ": ", hgnc, " ADR-eGene SUSIE Fine-mapping")

	file_name_susie_pval <- paste0(
    	"/mnt/vstor/Data14/metabrain_lasso/magenta_susie/ieGenes/", anc, "/",
    	hgnc, "_", gene, "_chr_", chr, "_", anc, 
    	"/MAGENTA_ieQTL_", gene, "_", hgnc, "_chr_", chr, "_", anc, 
    	"_susie_rss_all_pval_ad_int_corrected_n_est_resid_var.png"
	)

	print("Plotting started PVAL AD INT")

	# Open CairoPNG device
	CairoPNG(file_name_susie_pval, width = 8, height = 6, units = "in", dpi = 150, bg = "white")

	# Create the ggplot object
	p <- ggplot(ss_all, aes(x = pos, y = -log10(pvalue_ad_interaction))) +
    	geom_jitter(aes(color = CS > 0), width = 0.1, height = 0) +  # Add jitter to prevent overlap
    	scale_color_manual(values = c("black", "magenta")) +  # Define color mapping
    	labs(
       		x = "Position", 
        	y = "-log10(pvalue_ad_interaction)",  # Fixed missing closing parenthesis in label
        	title = no_anno_title, 
        	subtitle = "All variants"
    	) +
    	theme_minimal() +  # Apply a minimal theme for aesthetics
    	theme(
        	# Bold axis labels
        	axis.title.x = element_text(face = "bold"),
        	axis.title.y = element_text(face = "bold"),
        	# Bold axis levels
        	axis.text.x = element_text(face = "bold"),
        	axis.text.y = element_text(face = "bold"),
        	# Bold title
        	plot.title = element_text(face = "bold"),
        	# Bold legend title
        	legend.title = element_text(face = "bold")
    	)

	# Print the ggplot object
	print(p)

	# Close the CairoPNG device
	dev.off()

	print("Plotting PVAL Done")

	# PIP Plot
	 file_name_susie_pip = paste0("/mnt/vstor/Data14/metabrain_lasso/magenta_susie/ieGenes/",anc, "/", hgnc, "_", gene, "_chr_", chr, "_", anc, "/MAGENTA_ieQTL_", gene, "_", hgnc, "_chr_", chr, "_", anc, "_susie_rss_all_pip_corrected_n_est_resid_var.png")

	CairoPNG(file_name_susie_pip, width = 8, height = 6, units = "in", dpi = 150, bg = "white")
	print("Plotting started PIP")

	# Create the ggplot object
	p <- ggplot(ss_all, aes(x = pos, y = pip)) +
        	geom_jitter(aes(color = CS > 0), width = 0.1, height = 0) +  # Add jitter to prevent overlap
        	scale_color_manual(values = c("black", "magenta")) +  # Define color mapping
        	labs(x = "Position", y = "PIP", title = no_anno_title, subtitle = "All variants") +  # Label the axes
        	theme_minimal() +  # Apply a minimal theme for aesthetics
        	theme(
            	# Bold axis labels
            	axis.title.x = element_text(face = "bold"),
            	axis.title.y = element_text(face = "bold"),
            	# Bold axis levels
            	axis.text.x = element_text(face = "bold"),
            	axis.text.y = element_text(face = "bold"),
            	# Bold title
            	plot.title = element_text(face = "bold"),
            	# Bold legend title
            	legend.title = element_text(face = "bold")
        	)

	# Explicitly print the ggplot object
	print(p)

	dev.off()
	print("Plotting done PIP")


} else {
	print("No Credible Sets")
	ss_all$CS <- 0
}


## Save Plot


susie_plot_all_name = paste0("/mnt/vstor/Data14/metabrain_lasso/magenta_susie/ieGenes/",anc, "/", hgnc, "_", gene, "_chr_", chr, "_", anc, "/MAGENTA_ieQTL_", gene, "_", hgnc, "_chr_", chr, "_", anc, "_susie_rss_all_corrected_n_est_resid_var.png")

# png(susie_plot_all_name, width = 8, height = 5, units = "in", res = 600)
CairoPNG(susie_plot_all_name, width = 8, height = 5, units = "in", dpi = 150)

susie_plot(fitted_rss_eqtl_all, y="PIP", add_legend = "topleft", xlab = "Variable", main = paste0(anc, ": ", hgnc, " (", gene, ")"), sub = paste0("Chr: ", chr), cex.main = 1.5, cex.sub = 1.25, cex.lab = 1.25, cex.axis = 1.5)
dev.off()

print('SUSIE RSS Plot All complete')

print('SUSIE RSS ALL Convergence Status')
fitted_rss_eqtl_all$converged

