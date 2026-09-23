#
#   Copyright 2001-2026 by the individuals mentioned in the source code history
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# ---------------------------------------------------------------------
# Program: wccFindDyadParam.R
#  Author: Steven Boker
#    Date: Mon Feb 2 13:30:58 EST 2026
#
# This function receives a set of parameters to optimize and two arrays of
#   timeseries.  It repeatedly calls wccSurrogateDyads() and wccAggregate() 
#   for each parameter choice and returns the parameters for which real results are
#   most different than surrogate distributions for the results defined by opt1, opt2, and opt3.
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- Mon Feb 2 13:31:03 EST 2026
#      Created wccFindDyadParam.R
#
# ---------------------------------------------------------------------


wccTestSurrogateDyads <- function(inArray1=NA, inArray2=NA, wMax=50, tMax=50, wInc=1, tInc=1,
                             Lsize=8, pspan=.25, type="Max", nSurrogates=NA, samplespersecond=1, method=c("c", "cumr", "cumc", "r"), embedD=9, ...) {
    # Deprecation: allow old windcross argument
    dots <- list(...)
    testFrame <- wccFindDyadParam(inArray1=inArray1, inArray2=inArray2, wMaxvector=c(wMax), tMaxvector=c(tMax), wIncvector=c(wInc), tIncvector=c(tInc),
                                 Lsizevector=c(Lsize), pspanvector=c(pspan), type=type, nSurrogates=nSurrogates, samplespersecond=samplespersecond, 
                                 method=method, embedD=embedD, ...) 
    return(testFrame)
}