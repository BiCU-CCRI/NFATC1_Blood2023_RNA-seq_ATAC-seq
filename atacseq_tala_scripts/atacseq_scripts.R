#' rm_blacklisted_regions
#'
#' Removes blacklisted regions from peaks
#'
#' @param peak_gr peaks GRanges object
#' @param blacklisted_regions blacklisted regions GRanges object 
#'
#' @return GRanges object of peaks with blacklisted regions removed
#'
#' @examples rm_blacklisted_regions(peak_gr = peak_grlist[[sample_name]]
remove_blacklisted <- function(peak_gr=NULL, blacklisted_regions=blacklist){
  import::here(.from = IRanges, overlapsAny)  # IRanges::overlapsAny
  # removing regions overlapping with blacklisted
  # see https://biodatascience.github.io/compbio/bioc/ranges.html
  # If we just wanted to subset to the genes which overlap a given range, we can use overlapsAny:
  #g[overlapsAny(g, query[1])]
  #This is equivalent to the following:
  #g[g %over% query[1]] 
  
  blacklisted = sum(overlapsAny(peak_gr, blacklisted_regions))
  blacklisted_perc = (blacklisted/length(peak_gr))*100
  #not_blacklisted = sum(!overlapsAny(peak_gr, blacklisted_regions))
  message("Blacklisted regions regions: ", blacklisted, " (", round(blacklisted_perc,1), "%)")
  
  return(peak_gr[!overlapsAny(peak_gr, blacklisted_regions)])
}


#' add_flag
#'
#' Add flags to pheatmap
#' @description ADD more detailed description
#'
#' @param pheatmap refers to the original heatmap produced from the pheatmap() function
#' @param kept.labels should be a vector of labels you wish to show
#' @param repel.degree is a number in the range [0, 1], controlling how much the labels are spread out from one another
#' @param hiden.labels.character character to hide labels; e.g. if need to add spacing?!
#'
#' @return updated pheatmap
#'
#' @references original function from Z.Lin https://stackoverflow.com/questions/52599180/partial-row-labels-heatmap-r
#'
#' @examples add_flag(signifLog2Promoters_heatmap, kept.labels = highlight_genes, repel.degree = 0.2,hiden.labels.character = "")
add_flag <- function(pheatmap,
                     kept.labels,
                     repel.degree,
                     hiden.labels.character="") {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  # PR: hiden.labels.character - character to hide labels; e.g. if need to add spacing?!
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  # new.label$label <- ifelse(new.label$label %in% kept.labels, 
  #                           new.label$label, "")
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, hiden.labels.character)
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != hiden.labels.character)
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != hiden.labels.character],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != hiden.labels.character] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}