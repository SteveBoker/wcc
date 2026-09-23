# ------- Version 0.4.2 -------

2026-09-23  Steven Boker  <smb3u@virginia.edu>

    Bug fixes: 
    1. Added warning to wccPlot when requested startwindow or endwindow is out of range.

    Manual edits and additions:
    1. wccTestSurrogateDyads
    
    Additions:
    1. Added function wccTestSurrogateDyads which is just a convenience wrapper for
       wccFindDyadParam.  This will be used when people have found a set of parameters
       and now want to test their data against surrogates.
    
# ------- Version 0.4.0 -------

2026-07-12  Steven Boker  <smb3u@virginia.edu>

    Bug fixes:
    1. Compiled and R versions of wccCalc now agree on last elapsed time calculation when wInc != 1
    
    Manual edits and additions:
    1. wccCalcBatch
    1. wccVectorFieldCalc (beta release)
    1. wccVectorFieldPlot (beta release)
    
    Additions:
     from Aaron Peikert
    1. Variety of speed-ups 
    2. Addition of two new methods for calculating a WCC grid: cumulative sum in R and in C.  
       Note that both of these methods are faster than previous versions, but fail when there are any
       missing values in either time series.
    3. Reuse of cumulative wcc calc results for pairings of surrogate series.
     from Steve Boker
    1. wccVectorFieldCalc()  This function calculates the direction of an isobar on a wcc grid.
       It is a beta function that is intended to help identify the elapsed time of the start and 
       end of intervals of synchronization.
    2. wccVectorFieldPlot()  This plots the results of wccVectorFieldCalc() as a slope field.
    
# ------- Version 0.3.2 -------

2026-04-15  Steven Boker  <smb3u@virginia.edu>

    Bug fixes:
    1. The compiled version of windcross returned NA for any window with any missing values in either time series.  
        It now calculates correlation on pairwise complete observations.
        
    2. Manual edits and additions
        a. Added defaults into argument descriptions.
        b. Note that timeseries 1 and 2 can be swapped.
        c. wcc-package referenced wcc instead of wccCalc.
        d. Improved lspan description in wccPeakPick
        e. noted "Max" and "Min" for maxMean and maxVar in wccAggregate.
        f. Better definition of what is considered to be a missing window.
        g. Better definition of the returned dataframe from wccFindDyadParam.
        h. better documentation for the arguments to wccHeatMap.
        i. Better recommendations for the use of wccPlot.
        
    Additions:
    1. wccVectorField
        a. Add the wccVectorField function.
        b. Calculate the forward time derivative for each point in a wcc matrix.
    2. wccPlotVectorField
        a. Add the wccPlotVectorField function
        b. Plot a heatmap of the vector field.
        c. Plot a traditional vector field

    4. Add single parameter distribution testing.
    
    5. Vignettes
    
# ------- Version 0.3.1 ------- 

2026-04-09  Steven Boker    <smb3u@virginia.edu>

    1. wccPlot.R:
        a. call wccCalc and wccPeakPick from within the function rather than loading data from disk.
        b. Change the arguments list to match the arguments to wccAggregate.
        c. Remove the pdffile argument and always plot to the current device.
        d. Change default to windcross=TRUE
        d. Add argument error checking.
        e. Add a legend.
    
    2. wccAggregate.R:
        a. Remove the option to save the result of wccCalc and wccPeakPick to disk.
        b. Remove argument savefile
        c. Add argument embedD=9
        d. Change defaults to samplespersecond=1 and windcross=TRUE
        e. Add an aggregation of the 2nd derivative of the peakpicking lag timeseries.
        f. Add "d2lagMean", "d2lagVar" to wccAggNames

    3. wccFindDyadParam.R:
        a. Change defaults to samplespersecond=1 and windcross=TRUE
        b. Add argument embedD=9
        c. Change call to wccAggregate to include embedD
        d. Add "d2lagMeanKS", "d2lagVaKSr", "d2lagMeanQdiff", "d2lagVarQdiff" to wccTestNames
        e. Calculate values for "d2lagMeanKS", "d2lagVaKSr", "d2lagMeanQdiff", and "d2lagVarQdiff" 
    
    4. wccSurrogateDyads.R:
        a. Change defaults to windcross=TRUE
        c. Add argument embedD=9
        d. Change call to wccAggregate to include embedD
        
    5. wccCalc.R:
        a. Change defaults to windcross=TRUE
    
# ------- Version 0.3.0 ------- 

2026-04-03  Steven Boker    <smb3u@virginia.edu>

    1. windcross.c
        a. Rewrote to directly connect to R.
    2. wccCalc.R
        a. Added call to windcross.

