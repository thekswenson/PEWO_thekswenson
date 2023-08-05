#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument
if (length(args) < 2) {
  stop(paste0("The file results.csv must be supplied as 1st argument and the ",
              "workdir as 2nd argument (input file)."), call. = FALSE)
}

library(RColorBrewer)
library(grid)
library(ggplot2)
library(Cairo)

# Association software/parameters
# Need to be manually defined prior to analysis
# There should not be new placement software too
# often, so this is manageable.
#################################################

epa <- c("g")
epang_h1 <- c("g")
epang_h2 <- c("bigg")
epang_h3 <- c()
epang_h4 <- c()
pplacer <- c("ms", "sb", "mp")
rappas <- c("k", "o", "red", "ar")
rappas2 <- c("k", "o", "red", "ar", "filter", "sm")
apples <- c("meth", "crit")
appspam <- c("mode", "w", "pattern")

soft_params <- list(epa = epa, epang_h1 = epang_h1, epang_h2 = epang_h2,
                    epang_h3 = epang_h3, epang_h4 = epang_h4, pplacer = pplacer,
                    rappas = rappas, rappas2 = rappas2, apples = apples,
                    appspam = appspam)


# Associate read generator to paramaters
#################################################

partition <- list()
art <- list("art.co", "art.pl")

rgen_params <- list(PARTITION = partition, ART = art)

#helper functions
##################################################

#Compute all combinations of paramaters for the read generators that have data.
#Return a list of lists of pairs, where the first element of each pair is
#a column name, and the second is the value for that column.
#
#   EXAMPLE:
#     ((("rgen", "ART"), ("art.co", "ILLUMINA"), ("art.pl", "HS25")),
#      (("rgen", "ART"), ("art.co", "ILLUMINA"), ("art.pl", "HS25")),
#      (("rgen", "PARTITION")))
#
infer_read_generators <- function(data) {
  retlist <- list()
  rgenerators <- levels(data$rgen)
  for (rgen in rgenerators) {
    if (length(rgen_params[[rgen]])) {	#combine the second level parameters
      p2lists <- list()
      #collect the paramters and their values into a list of lists
      for (column in rgen_params[[rgen]]) {
        clist <- list()

        for (val in levels(data[, column])) {
          if (val != "") {
            clist <- append(clist, list(list(column, val)))
          }
        }
        p2lists <- append(p2lists, list(clist))	#yes, append() concatenates!
      }

      #now take the product of all the paramater lists
      combos <- expand.grid(p2lists)			#all combinations of parameters
      for (i in seq_len(nrow(combos))) {	#create the list of pairs for each row
        aslist <- list(list("rgen", rgen))
        for (pair in combos[i, ]) {
          aslist <- append(aslist, pair)
        }
        retlist[[length(retlist) + 1]] <- aslist
      }
    } else {
      retlist[[length(retlist) + 1]] <- list(list("rgen", rgen))
    }
  }
  return(retlist)
}

#Select the data that corresponds to the column/value pairs in cvpairs.
#These input columns/value pairs are specified the output of
#infer_read_generators().
#Return Dataset corresponding to cvpairs.
select_data_for_rgen <- function(data, cvpairs) {
  columns <- c()
  values <- c()
  for (pair in cvpairs) {			#widdle down the dataframe base on each pair
    columns <- c(columns, pair[[1]])
    prefix <- strsplit(pair[[1]], ".", fixed = TRUE)[1]

    values <- c(values, paste0(tail(prefix[[1]], n = 1), pair[[2]]))
    data <- data[data[, pair[[1]]] == pair[[2]], ]
  }
  return(new("Dataset", df = data, name = paste(values, collapse = "-"),
             formula = paste(columns, collapse = " + ")))
}


#classes
##################################################

setRefClass("Dataset",
            fields = list(
              df = "data.frame",       #the dataframe
              softname = "character",  #name of software
              name = "character",      #parameter combination string
              formula = "character"    #parameter combination formula
            ),
            methods = list(
              str = function() {
                return(paste(softname, name, sep = "_"))
              }
            ))


#main function
##################################################

analyze_data <- function(measure, measurestr) {
  #ND heatmaps per parameters
  alltables <- list()         #list of Datsets
  tablei <- 1
  for (softname in soft_analyzed) {
    #epang is treated separatly for each heuristic so it has a _ in the name
    softname_short <- strsplit(softname, "_")[[1]][1]

    print(paste0(measurestr, " heatmap for ", softname))
    #select data for current software
    if (softname_short != "epang") {
      current_soft_data <- data[data$software == softname_short, ]
    } else {
      heur <- substr(strsplit(softname, "_")[[1]][2], 2, 10)
      current_soft_data <- data[data$software == softname_short &
                                  data$h == as.numeric(heur), ]
    }
    #remove columns with only NA, meaning this parameter was not linked to the
    #current soft
    current_soft_data <- current_soft_data[, colSums(is.na(current_soft_data))
                                           != nrow(current_soft_data)]

    #build formulas dynamically for parameters
    formula_mean <- paste0(measure, " ~ pruning + r")
    formula_meanofmean <- paste0(measure, " ~ r")
    print(paste0("  soft_params: ", soft_params[softname][[1]]))
    if (length(soft_params[softname][[1]]) > 0) {  #is == 0 when no params
      for (j in seq_along(soft_params[softname][[1]])) {
        formula_mean <- paste0(formula_mean, " + ",
                               soft_params[softname][[1]][j])
        formula_meanofmean <- paste0(formula_meanofmean, " + ",
                                     soft_params[softname][[1]][j])
      }
    }
    print(paste0("  formula_mean: ", formula_mean))

    rgenerators <- infer_read_generators(current_soft_data)
    for (rgen in rgenerators) {
      #select data for the current read generator
      current <- select_data_for_rgen(current_soft_data, rgen)

      #aggregate as mean per pruning
      data_mean <- aggregate(as.formula(paste0(formula_mean, " + ",
                                               current$formula)),
                             current$df, mean)
      write.table(data_mean,
                  file = paste0(workdir, "/mean_per_pruning_", measurestr, "_",
                                current$str(), ".csv"),
                  quote = TRUE, sep = ";", dec = ".", row.names = FALSE,
                  col.names = TRUE)

      #aggregate as mean of means
      df_meanofmean <- aggregate(as.formula(formula_meanofmean), data_mean,
                                 mean)
      #order from best to wort parameters combination
      current$df <- df_meanofmean[order(df_meanofmean[[measure]]), ]
      current$df["software"] <- softname    #is this necessary?
      current$softname <- softname

      #register results

      alltables[[tablei]] <- current

      #ouputs results table per software
      print(paste0("CSV table for ", softname))
      write.table(current$df,
                  file = paste0(workdir, "/summary_table_", measurestr, "_",
                                current$str(), ".csv"),
                  quote = TRUE, sep = ";", dec = ".", row.names = FALSE,
                  col.names = TRUE)

      tablei <- tablei + 1
    }
  }

  #search for measure (e.g. e_ND, ND) min/max
  min_nd <- Inf
  max_nd <- 0
  for (table in alltables) {
    mi <- min(table$df[[measure]])
    ma <- max(table$df[[measure]])
    if (mi < min_nd) {
      min_nd <- mi
    }
    if (ma > max_nd) {
      max_nd <- ma
    }
  }

  #build all plots
  global_labeller <- labeller(
    .default = label_both
  )

  for (table in alltables) {
    softname_short <- strsplit(table$softname, "_")[[1]][1]
    params <- soft_params[table$softname][[1]]
    #if 0 parameter, build heatmap on fake x/y
    if (length(params) == 0) {
      table$df["none"] <- "none"
      params <- c(params, "none")
    }

    #if 1 parameter, build heatmap on fake y
    if (length(params) == 1) {
      table$df["none"] <- "none"
      params <- c(params, "none")
    }
    #if more than 2 parameters, build a facet_wrap combination
    wrap_string <- "~r"
    wrapcount <- 0
    if (length(params) > 2) {
      for (j in 3:length(params)) {
        wrap_string <- paste(wrap_string, params[j], sep = "+")
        wrapcount <- wrapcount + 1
      }
    }

    #build aes string from 2 first params + nd as fill
    #              ncol matches read length (param r)
    nrow <- length(unique(data[!is.na(data$r), ]$r))
    g <- ggplot(table$df, aes_string(x = sprintf("factor(%s)", params[1]),
                                     y = sprintf("factor(%s)", params[2])))
    g <- g + geom_tile(aes(fill = !!sym(measure))) #NOTE: measure is ND or e_ND
    g <- g + facet_wrap(as.formula(wrap_string), labeller = global_labeller,
                        nrow = nrow)
    g <- g + geom_text(aes(label = sprintf("%0.2f",
                                           round(!!sym(measure), digits = 2)),
                           size = 3))
    g <- g + scale_fill_distiller(limits = c(min_nd, max_nd),
                                  palette = "RdYlGn")
    g <- g + labs(title = paste("mean ", measurestr, ": ", table$str()),
                  x = paste("parameter: '", params[1], "'"),
                  y = paste("parameter: '", params[2], "'"))

    #(2 * parameter uniq value) * (combinations of 3rd to nth params) +
    # space for legend on the right
    columns <- 1
    if (length(params) > 2) {
      for (j in 3:length(params)) {
        columns <- columns * length(unique(table$df[[params[j]]]))
      }
    }
    svg_width <- 2 + (0.7 * length(unique(table$df[[params[1]]])) * columns)
    svg_height <- 1 + (0.7 * length(unique(table$df[[params[2]]])))
    CairoSVG(file = paste0(workdir, "/summary_plot_", measurestr, "_",
                           table$str(), ".svg"),
             width = svg_width, height = svg_height)
    print(g)
    dev.off()
  }
}

#load data
##################################################
data <- read.csv(file = args[1], sep = ";", header = TRUE)
workdir <- args[2]

#define list with software that were actually tested and remove them from
#soft_list and soft_param accordingly
#soft_list<-list("epa", "epang_h1", "epang_h2", "epang_h3", "epang_h4",
#                "pplacer", "rappas", "apples")
soft_analyzed <- levels(data$software)
#for epang, test which algorithms were tested
epang_idx <- match("epang", soft_analyzed)
if (!is.na(epang_idx)) {
  soft_analyzed <- soft_analyzed[-epang_idx]
  epang_algos <- unique(data[!is.na(data$h), ]$h)
  for (h in sort(epang_algos)) {
    soft_analyzed <- c(soft_analyzed, paste0("epang_h", h))
  }
}

analyze_data("nd", "ND")
analyze_data("e_nd", "eND")