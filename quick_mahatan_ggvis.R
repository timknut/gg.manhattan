mtcars %>% 
  ggvis(~wt, ~mpg) %>% 
  layer_points(fill = ~factor(cyl))

install.packages("qqman")


require(qqman)
require(ggvis)


df <- gwasResults

df %>% 
  arrange(CHR, BP) %>% 
  mutate(bp_2 = cumsum(BP), chr_col = (CHR %% 2) == 0 ) %>%
  mutate(col2 = ifelse(chr_col, "black", "grey")) %>% 
  ggvis(~bp_2, ~-log10(P)) %>% 
  add_tooltip(function(df) df$bp_2)
  #layer_points(stroke := rep(c("red", "green")),11)
