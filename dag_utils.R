library(dagitty)
library(ggdag)
library(cowplot)
library(data.table)
library(knitr)
library(boot)

plot_dag = function(dag, curves = T) {
  
  tdag = data.table(tidy_dagitty(dag)$data)
  tdag[!is.na(direction), curvature := 0]
  tdag[name == "L" & to == "S", curvature := ifelse(curves == T, -1, 0)]
  
  
  g = ggplot(tdag,
             aes(x = x, y = y,
                 xend = xend, yend = yend)) + 
    theme_dag() + 
    theme(panel.background = element_rect(fill = "white"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_dag_point(col = "white", show.legend = F) +
    geom_dag_text(col = "black") 
  
  if (any(tdag$curvature == -1,na.rm = T)) {
    g = g + geom_dag_edges_arc(curvature = tdag[!is.na(direction),
                                                 curvature])
  } else {
    g = g + geom_dag_edges()
  }
    
  
  for (k in grep("S$|S[0-9]$",tdag$name)) {
    g = g + geom_point(data = tdag[k],
                       pch = 1,
                       size = 10)
  }
  
  adjusted_vars = sapply(grep("adjusted",
                              strsplit(dag,"\n")[[1]],
                              value = T),
                         function(x)
                           strsplit(x," ")[[1]][1])
  
  for (v in adjusted_vars) {
    g = g + geom_point(data = tdag[name == v],
                       pch = 1, size = ifelse(nchar(v)<3,10,11))
  }
  
  g$coordinates$clip = "off"
  return(g)
}

plot_hists = function(x) {
  betaX = c()
  Estimates = c()
  for (k in 1:ncol(x)) {
    betaX = c(betaX,rep(colnames(x)[k], nrow(x)))
    Estimates = c(Estimates, x[,k])
  }
  
  if (min(x) >= 0) {
    xlim = c(0,quantile(x,.999))
  } else {
    xlim = c(quantile(x,1-.999),quantile(x,.999))
  }
  
  pdata = data.frame(betaX, Estimates)
  pdata$betaX = gsub("betax_","",pdata$betaX)
  
  p = ggplot(pdata, aes(x = Estimates,
                        group = betaX,
                        fill = betaX)) + 
    theme_classic() +
    geom_histogram(alpha = .5,
                   aes(y = ..density..),
                   position = 'identity',
                   bins = 30) +
    geom_vline(xintercept = 0, col = "green4", lwd = 1) + 
    xlab("Bias") + 
    theme(legend.position="top") + 
    xlim(xlim)
  p
  return(p)
}

plot_daghist = function(DAG,betas, curves = T, rel_heights = c(2,3)){
  plot_grid(
    plot_dag(DAG, curves = curves),
    plot_hists(betas),
    ncol = 1,
    rel_heights = rel_heights
  )
}

get_ipws = function(P,X, only_included = T, weights  = NULL) {
  if (class(X)[1] != "matrix")
    X = as.matrix(X)
  if (is.null(weights)) {
    selection_model = 
      glm(P ~ X, family = binomial)  
  } else {
    selection_model = 
      glm(P ~ X, family = binomial, weights = weights)
  }
  
  participation_probability = 
    predict(selection_model,type = "response")
  
  if(only_included == T)
    participation_probability = participation_probability[P]
  
  IPW = mean(P)/participation_probability
  return(IPW)
}


test_implications = function(DAG,dt) {
  library(rstanarm)
  icid = impliedConditionalIndependencies(DAG)
  tests = c()
  test_names = c()
  for (k in 1:length(icid)) {
    Yvar =  icid[[k]]$Y
    Xvar =  icid[[k]]$X
    Zvar =  icid[[k]]$Z
    Zvar = ifelse(length(Zvar) == 0,
                  "",
                  paste0(" + ", Zvar))
    f = as.formula(paste0(Yvar," ~ ", Xvar, Zvar))
    
    sf = stan_glm(f,
                  data = dt,
                  chains = 1)
    icid[[k]]$test = posterior_interval(sf, pars = Xvar)
    tests = rbind(tests,icid[[k]]$test)
    
    test_names = c(test_names,
                   paste0(Yvar,
                          " _||_ ",
                          Xvar,
                          gsub("\\+","| ",Zvar)))
  }
  rownames(tests) = test_names
  
  ys = 1:nrow(tests)
  cols = apply(tests,1,
               function(x) 
                 ifelse(mean(x)> .1,
                        "red",
                        "black"))
  par (mar=c(3,7,.5,.5), mgp=c(2,.7,0), tck=-.01)
  plot(0,
       type = "n",
       xlim = range(tests),
       ylim = range(ys),
       yaxt = "n", ylab = "",
       xlab = "coefficient")
  polygon(x = c(-.1,.1,.1,-.1,-.1),
          y = c(rep(par("usr")[3],2),
                rep(par("usr")[4],2),
                par("usr")[3]),
          col = adjustcolor("green3",alpha = .5),
          border = NA)
  abline(v = 0, lty = 2)
  segments(x0 = tests[,1],
           x1 = tests[,2],
           y0 = ys,
           lwd = 2,
           col = cols)
  axis(2,
       at = ys,
       labels = rownames(tests),
       las = 2)
}

sim_dag = function(DAG,
                   N = 1000,
                   show_sim_sequence = F,
                   sigma_S = 2) 
  {
  variables = unique(c(levels(edges(DAG)[,1]),
                       levels(edges(DAG)[,2])))
  
  dt = matrix(rnorm(N*length(variables)),
              ncol = length(variables))
  dt = apply(dt,2,scale)
  colnames(dt) = variables
  
  if ("S" %in% variables)
    dt[,"S"] = dt[,"S"]*sigma_S
  for (Sx in grep("S[0-9]",colnames(dt), value = T))
    dt[,Sx] = dt[,Sx]*sigma_S
  if ("U" %in% variables)
    dt[,"U"] = dt[,"U"]*1
  for (Ux in grep("U[0-9]",colnames(dt), value = T))
    dt[,Ux] = dt[,Ux]*1
  
  dag_egdes = edges(DAG)
  missing_edges = c()
  gen_model = ""
  for (k in 1:nrow(dag_egdes)) {
    if (!(dag_egdes[k,1] %in% dag_egdes[,2])) {
      cs = as.character(dag_egdes[k,1])
      ef = as.character(dag_egdes[k,2])
      dt[,ef] = scale(dt[,ef] + dt[,cs])
      gen_model = paste0(gen_model,
                         cs,
                         " -> ",
                         ef,
                         "\n")
    } else {
      missing_edges = c(missing_edges,k)
    }
  }
  
  missing_edgesb = c()
  for (k in missing_edges) { 
    if (!(dag_egdes[k,1] %in% dag_egdes[missing_edges,2])) {
      cs = as.character(dag_egdes[k,1])
      ef = as.character(dag_egdes[k,2])
      dt[,ef] = scale(dt[,ef] + dt[,cs])
      gen_model = paste0(gen_model,
                         cs,
                         " -> ",
                         ef,
                         "\n")
    } else {
      missing_edgesb = c(missing_edgesb,k)
    }
  }
  
  for (k in missing_edgesb) { 
      cs = as.character(dag_egdes[k,1])
      ef = as.character(dag_egdes[k,2])
      dt[,ef] = scale(dt[,ef] + dt[,cs])
      gen_model = paste0(gen_model,
                         cs,
                         " -> ",
                         ef,
                         "\n")
  }
  
  dt = data.table(dt)
  if ("S" %in% variables)
    dt[, P := S > 0] # same as inv.logit(S) > .5
  
  #enforce dependence of one selection variable on the other
  # i.e. if and only if S1 -> S2, S2 must be 0 if S1 is 0
  for (Sx in grep("S[0-9]",names(dt), value = T)) {
    P = sub("S","P",Sx)
    dt[[P]] = T
    dt[,(P) := get(Sx) > 0]
    SxP = dag_egdes[dag_egdes$w == Sx,"v"]
    if (length(grep("S",SxP))>0) {
      dt[get(grep("S",SxP,value = T)) < 0, (P) := F]
    }
  }
  
  if (show_sim_sequence == T)
    cat(gen_model)
  return(dt)
}