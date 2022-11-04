#!/usr/bin/env Rscript

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 12) stop("Wrong number of input parameters: ''")

library(ape)
library(apTreeshape)
library(Brobdingnag)

#receive parameters
rho <- inputpar[1]        #0.0005 # sampling probability
d_t <- inputpar[2]        #0.01
time_scale <- inputpar[3] #8
b_rates <- inputpar[4]    #comma seperated
d_rates <- inputpar[5]    #comma seperated
mu <- inputpar[6]
timeFrom <- inputpar[7]   #0
timeTo <- inputpar[8]     #20
timeBy <- inputpar[9]     #0.1
non_negativity_cutoff <- inputpar[10]  #0
treeFile <- inputpar[11]  #mutation File
outFile <- inputpar[12]   #output file

#prepare parameters
rho = as.numeric(rho)
d_t = as.numeric(d_t)
time_scale = as.numeric(time_scale)
b_rate = as.numeric(strsplit(b_rates, ",")[[1]])
d_rate = as.numeric(strsplit(b_rates, ",")[[1]])
mu <- replicate(length(b_rate), as.numeric(mu))
#argument <- list(b_rate, d_rate, mu)
T_vector <- seq(from = as.numeric(timeFrom), to = as.numeric(timeTo), by = as.numeric(TimeBy))
non_negativity_cutoff = as.numeric(non_negativity_cutoff)


phy = mutationMatrix2Tree(treeFile)
outcome = ith.Fitness(phy, outFile, rho, d_t, time_scale, b_rate, d_rate, mu, T_vector, non_negativity_cutoff)
