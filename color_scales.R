color.gradient.discrete = function(color.low, color.high, n) {
  scales::seq_gradient_pal(low=color.low, high=color.high)(seq(0, 1, length.out = n)) 
  #ggplot2::scale_color_gradient ? doesn't work well :/ also not ggplot2::scale_color_gradientn
}
  
color.gradient.discrete.single = function(color, n, color.low="white", includeEnd=F) {
  if (includeEnd) return(color.gradient.discrete(color.low, color, n))
  else return(color.gradient.discrete(color.low, color, n+1)[-1])
}
  
color.gradient.discrete.divergent = function(color.low, color.high, n, color.middle = "white") {
  if (n %% 2 == 0) { #even number => don't include middle color
    c(color.gradient.discrete.single(color.low, n/2, color.middle, includeEnd=F) %>% rev(),
      color.gradient.discrete.single(color.high, n/2, color.middle, includeEnd=F)
    ) %>% return()
  } else { #odd number => round up but include middle color only once
    c(color.gradient.discrete.single(color.low, ceiling(n/2), color.middle, includeEnd=T) %>% rev(),
      color.gradient.discrete.single(color.high, ceiling(n/2), color.middle, includeEnd=T)
    ) %>% unique() %>% return()
  }
}
  
colors = color.gradient.discrete("blue", "red", n=threatLevels.n)
#scales::show_col(colors)