# Slightly modify locuscompare to allow the user to specify the axis dimensions for all plots

require(ggplot2)
require(locuscomparer)
require(cowplot)

#if (FALSE) {
add_label = function(merged, snp){
	    merged$label = ifelse(merged$rsid %in% snp, merged$rsid, '')
    return(merged)
}

make_scatterplot_custom = function (merged, title1, title2, color, shape, size, xmax, ymax, legend = FALSE,
	      legend_position = c("bottomright", "topright", "topleft"))
{
	p = ggplot(merged, aes(logp1, logp2)) + geom_point(aes(fill = rsid,
	size = rsid, shape = rsid), alpha = 0.8) + geom_point(data = merged[merged$label !=
																	               "", ], aes(logp1, logp2, fill = rsid, size = rsid, shape = rsid)) +
        xlab(bquote(.(title1) ~ -log[10] * "(P)")) + ylab(bquote(.(title2) ~
								         -log[10] * "(P)")) + scale_fill_manual(values = color,
						              guide = "none") + scale_shape_manual(values = shape,
							              guide = "none") + scale_size_manual(values = size, guide = "none") +
        ggrepel::geom_text_repel(aes(label = label)) + theme_classic() + xlim(0, xmax) + ylim(0, ymax)
    if (legend == TRUE) {
	            legend_position = match.arg(legend_position)
            if (legend_position == "bottomright") {
		                legend_box = data.frame(x = 0.8, y = seq(0.4, 0.2,
									                 -0.05))
	            }
	            else if (legend_position == "topright") {
			                legend_box = data.frame(x = 0.8, y = seq(0.8, 0.6,
										                 -0.05))
		            }
		            else {
				                legend_box = data.frame(x = 0.2, y = seq(0.8, 0.6,
											                 -0.05))
			            }
			            p = ggdraw(p) + geom_rect(data = legend_box, aes(xmin = x,
										                 xmax = x + 0.05, ymin = y, ymax = y + 0.05), color = "black",
							                  fill = rev(c("blue4", "skyblue", "darkgreen", "orange",
										                       "red"))) + draw_label("0.8", x = legend_box$x[1] +
							                  0.05, y = legend_box$y[1], hjust = -0.3, size = 10) +
            draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2],
		                       hjust = -0.3, size = 10) + draw_label("0.4",
				                   x = legend_box$x[3] + 0.05, y = legend_box$y[3],
				                   hjust = -0.3, size = 10) + draw_label("0.2", x = legend_box$x[4] +
						               0.05, y = legend_box$y[4], hjust = -0.3, size = 10) +
            draw_label(parse(text = "r^2"), x = legend_box$x[1] +
		                       0.05, y = legend_box$y[1], vjust = -2, size = 10)
	        }
        return(p)
}

make_locuszoom_custom = function (metal, title, chr, color, shape, size, ymax, ylab_linebreak = FALSE)
{
	    p = ggplot(metal, aes(x = pos, logp)) + geom_point(aes(fill = rsid,
								           size = rsid, shape = rsid), alpha = 0.8) + geom_point(data = metal[metal$label !=
																	              "", ], aes(x = pos, logp, fill = rsid, size = rsid, shape = rsid)) +
	ylim(0,ymax) +
        scale_fill_manual(values = color, guide = "none") + scale_shape_manual(values = shape,
									               guide = "none") + scale_size_manual(values = size, guide = "none") +
        scale_x_continuous(labels = function(x) {
				               sprintf("%.1f", x/1e+06)
					               }) + ggrepel::geom_text_repel(aes(label = label)) + xlab(paste0("chr",
														               chr, " (Mb)")) + ylab(bquote(.(title) ~ -log[10] * "(P)")) +
        theme_classic() + theme(plot.margin = unit(c(0.5, 1,
						             0.5, 0.5), "lines"))
	    if (ylab_linebreak == TRUE) {
		            p = p + ylab(bquote(atop(.(title), -log[10] * "(P)")))
	        }
	    return(p)
}



make_combined_plot_custom = function (merged, title1, title2, ld, chr, xmax, ymax, snp = NULL, combine = TRUE,
			           legend = FALSE, legend_position = c("bottomright", "topright",
								              "topleft"), lz_ylab_linebreak = FALSE)
{
	    snp = get_lead_snp(merged, snp)
    color = assign_color(merged$rsid, snp, ld)
        shape = ifelse(merged$rsid == snp, 23, 21)
        names(shape) = merged$rsid
	    size = ifelse(merged$rsid == snp, 3, 2)
	    names(size) = merged$rsid
	        merged = add_label(merged, snp)
	        p1 = make_scatterplot_custom(merged, title1, title2, color, shape,
				              size, xmax, ymax, legend, legend_position)
		    metal1 = merged[, c("rsid", "logp1", "chr", "pos", "label")]
		    colnames(metal1)[which(colnames(metal1) == "logp1")] = "logp"
		        p2 = make_locuszoom_custom(metal1, title1, chr, color, shape, size, xmax,
					            lz_ylab_linebreak)
		        metal2 = merged[, c("rsid", "logp2", "chr", "pos", "label")]
			    colnames(metal2)[which(colnames(metal2) == "logp2")] = "logp"
			    p3 = make_locuszoom_custom(metal2, title2, chr, color, shape, size, ymax,
						        lz_ylab_linebreak)
			        if (combine) {
					        p2 = p2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
				        p4 = cowplot::plot_grid(p2, p3, align = "v", nrow = 2,
								            rel_heights = c(0.8, 1))
					        p5 = cowplot::plot_grid(p1, p4)
					        return(p5)
						    }
			        else {
					        return(list(locuscompare = p1, locuszoom1 = p2, locuszoom2 = p3))
				    }
}



locuscompare_custom = function (in_fn1, in_fn2, xmax, ymax, marker_col1 = "rsid", pval_col1 = "pval",
	      title1 = "eQTL", marker_col2 = "rsid", pval_col2 = "pval",
	          title2 = "GWAS", snp = NULL, population = "EUR", combine = TRUE,
	          legend = TRUE, legend_position = c("bottomright", "topright",
						             "topleft"), lz_ylab_linebreak = FALSE, genome = c("hg19",
							             "hg38"))
{
	    d1 = read_metal(in_fn1, marker_col1, pval_col1)
    d2 = read_metal(in_fn2, marker_col2, pval_col2)
        merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"),
		               all = FALSE)
        genome = match.arg(genome)
	    merged = get_position(merged, genome)
	    chr = unique(merged$chr)
	        if (length(chr) != 1)
			        stop("There must be one and only one chromosome.")
	        snp = get_lead_snp(merged, snp)
		    ld = retrieve_LD(chr, snp, population)
		    p = make_combined_plot_custom(merged, title1, title2, ld, chr, xmax, ymax, snp,
					           combine, legend, legend_position, lz_ylab_linebreak)
		        return(p)
}


