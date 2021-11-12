############################################################
## Readme.txt
############################################################

Reference: Bayesian personalized treatment selection strategies that integrate predictive with 
           prognostic determinants
Junsheng Ma, Fracesco C. Stingo and Brian P. Hobbs

This readme.txt file provides details of the R scrips of "HCBPP.R" and "Lasso.R", which were 
written to implement the treatment selection methods proposed in the above referred paper. 
The key functions for model estimation and summary measure calculations are listed in the R script of
"Functions.R", and the example data are included as an R object of "LGGdata.rda".

1. The LGG data as decribed in Section 4.1
   "outcom"    = an ordinal outcome variable with three-level of responses (0,1,2). Large values
                 indicate better results.
   "trt"       = the treatment variable. Values of 0 and 1 represent for non-targeted
                 and targeted treatment, respectively.
   "gene.norm" = the standarized protein expression data.

2. Methods
   (a) "HCBPP.R" is for the proposed method of BPP with hierachical clustering method.
   (b) "Lasso.R" is for the penalized regression approach of Lasso. 

3. Key functions

       mymultt: posterior density of the prognostic features used in the proposed Bayesian 
                predictive modeling approaches.
                input: measures of the prognostic features.
               output: density

   con.cluster: for the proposed Bayesian predictive modeling approaches
                this function calculates the predictive utility for both treatments
                    input: similarity measures
                   output: the recommended treatment

            LR: for the methods of Lasso and ridge regression; Similarly, this function
                calculates the predictive utility for both treatments
                    input: alpha, the indicator of Lasso and Ridge regression
                   output: the recommended treatment

        PreUtBPP: this function calculates the summary measures as described in Section 4.2.
                    input: predicted treatment utilities
                   output: the summary measure

4. Required R library packages
     library(NMF);
     library(ConsensusClusterPlus);
     library(glmnetcr);

5. Run the example R script of "ExampleRcode.R".
     Output:  Empirical summary measure and CPO counts.

   