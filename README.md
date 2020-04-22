# R-scripts
R scripts for paper "R package SamplingStrata: new developments and extension to Spatial Sampling"

This repository is intended to enable users to replicate what contained in the paper 
"R package SamplingStrata: new developments and extension to Spatial Sampling"
(authors: Marco Ballin and Giulio Barcaroli), published in arxiv.org.

Abstract

The R package SamplingStrata was developed in 2011 as an instrument to optimize the
design of stratified samples. The optimization is performed by considering the stratification
variables available in the sampling frame, and the precision constraints on target estimates of
the survey (Ballin & Barcaroli, 2014). The genetic algorithm at the basis of the optimization
step explores the universe of the possible alternative stratifications determining for each of
them the best allocation, that is the one of minumum total size that allows to satisfy the
precision constraints: the final optimal solution is the one that ensures the global minimum
sample size. One fundamental requirement to make this approach feasible is the possibility to
estimate the variability of target variables in generated strata; in general, as target variable
values are not available in the frame, but only proxy ones, anticipated variance is calculated
by modelling the relations between target and proxy variables. In case of spatial sampling, it
is important to consider not only the total model variance, but also the co-variance derived
by the spatial auto-correlation. The last release of SamplingStrata enables to consider both
components of variance, thus allowing to harness spatial auto-correlation in order to obtain
more efficient samples.
