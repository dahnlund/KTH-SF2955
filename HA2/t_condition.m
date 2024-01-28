function cond = t_condition(breakpoints, T)
    
    cond = (sum(breakpoints < T) == length(breakpoints)) * (sum(breakpoints > 0) == length(breakpoints)) * sum(breakpoints == sort(breakpoints)) == length(breakpoints);

end