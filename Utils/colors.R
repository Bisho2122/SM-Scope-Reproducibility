library(colorspace)
library(ggsci)
# Color design assignment -------------------------------------------------

#Color paletes needed

# 1. Cell-types palette (for both coarse and fine-grained)
# 2. Metabotypes palette
# 3. General palette for :
    # 3.1 : Diverging colors for continuouos scales. The upper of which is for
            #gradient.
    # 3.2 : Deriving colors for one-time plots
            # 3.2.1 3 colors for 3 categories
            # 3.2.2 4 colors for 4 categories
            # 3.2.3 2 colors for 2 categories
            # 3.2.4 11 colors for 11 categories
            # 3.2.5 10 colors for 10 categories (This could be paired in brewer)

#Palettes:
# A = qualitative_hcl(n = 11, h = c(-311, 344), c = 90, l = 75)
# B = RColorBrewer::brewer.pal(11, "Set3") #Maybe increase saturation a bit
A = ggsci::pal_npg("nrc", alpha = 0.9)(10)
B = ggsci::pal_d3("category20", alpha = 0.9)(11)

# barplot(rep(1,length(A)), col = A, border = NA, main = "Cell types: fine 11")
# barplot(rep(1,length(B)), col = B, border = NA, main = "Metabotypes")


#Plot     #Palette
#1        A,C
#2        B
#3.1      2 from B
#3.2.1    3 from B
#3.2.2    4 from B
#3.2.3    2 from B
#3.2.4    B
#3.2.5    10 from B

my_cols = list(fine_cell_types = A,
               coarse_cell_types = A[1:7],
               metabotypes = B,
               de_misty_both = B[c(1:3)],
               corr_labels_histo = B[c(1:4)],
               highlight_colors = B[c(1,4)],
               regions_parents = B,
               spatialglue_clusters =B,
               transc_leiden_clusters = B)

scale_cols = list("div" = list("low" = B[1],
                               "mid" = "#F7F7F7",
                               "high" = B[4]),
                  "gradient" = list("low" = "lightgrey",
                                    "high" = B[4]))

# Function ---------------------------------------------------------------

assign_colors = function(labels, cols){
  if(length(labels) != length(cols)){
    stop("Number of labels must be the same as colors")
  }
  vec = cols
  names(vec) = labels
  return(vec)
}

