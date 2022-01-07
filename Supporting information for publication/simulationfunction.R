##' This function is adapted from OrgMassSpecR and computes theoretical HDX spectra
##' given the incorporation and charge. Note the method is stochastic to there
##' will be some slight difference between different runs.
##' 
##' @title Compute theoretical HDX Spectra 
##' @param sequence The sequence for which we are computing the theoretical spectra
##' @param incorp The deuterium incorporation level at which to compute the spectra
##' @param charge The charge state of the desired spectra.
##' @param custom A list. This allows custom additions to the spectra for example
##' modifications
##' @param fast Should optimised C++ code be used.
##' @return An object of class spectra with the desired theoretical spectra
##' @md
##' 
##' @rdname hdx-distributions
isotopicDistributionHDX <- function(sequence,
                                    incorp = 0,
                                    charge = 1, 
                                    custom = list(code = NULL, elements = NULL),
                                    fast = TRUE) {
    
    if(length(custom$elements != 0)) {
        custom_elements <- c(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0)
        custom_elements[names(custom$elements)] <- custom$elements
    }
    
    if(charge < 0 | charge > 8) stop("charge must be between 1 and 8")
    
    # create vector for sequences and holder for element annotations
    seq_vector <- strsplit(sequence, split = "")[[1]]
    x <- c(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0)
    
    # Update x with amino acid information
    for(i in seq.int(length(seq_vector))) {
        if(seq_vector[i] == "A") x <- x + c(C = 3, H = 5, N = 1, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "R") x <- x + c(C = 6, H =12, N = 4, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "N") x <- x + c(C = 4, H = 6, N = 2, O = 2, S = 0, P = 0)
        if(seq_vector[i] == "D") x <- x + c(C = 4, H = 5, N = 1, O = 3, S = 0, P = 0)
        if(seq_vector[i] == "C") x <- x + c(C = 3, H = 5, N = 1, O = 1, S = 1, P = 0) 
        if(seq_vector[i] == "E") x <- x + c(C = 5, H = 7, N = 1, O = 3, S = 0, P = 0)
        if(seq_vector[i] == "Q") x <- x + c(C = 5, H = 8, N = 2, O = 2, S = 0, P = 0)
        if(seq_vector[i] == "G") x <- x + c(C = 2, H = 3, N = 1, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "H") x <- x + c(C = 6, H = 7, N = 3, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "I") x <- x + c(C = 6, H =11, N = 1, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "L") x <- x + c(C = 6, H =11, N = 1, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "K") x <- x + c(C = 6, H =12, N = 2, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "M") x <- x + c(C = 5, H = 9, N = 1, O = 1, S = 1, P = 0)
        if(seq_vector[i] == "F") x <- x + c(C = 9, H = 9, N = 1, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "P") x <- x + c(C = 5, H = 7, N = 1, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "S") x <- x + c(C = 3, H = 5, N = 1, O = 2, S = 0, P = 0)
        if(seq_vector[i] == "T") x <- x + c(C = 4, H = 7, N = 1, O = 2, S = 0, P = 0)
        if(seq_vector[i] == "W") x <- x + c(C =11, H =10, N = 2, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "Y") x <- x + c(C = 9, H = 9, N = 1, O = 2, S = 0, P = 0)
        if(seq_vector[i] == "V") x <- x + c(C = 5, H = 9, N = 1, O = 1, S = 0, P = 0)
        
        if(length(custom$elements != 0))
            if(seq_vector[i] == custom$code) x <- x + custom_elements    
    }
    
    ## add N-terminal H and C-terminal OH
    elements <- x + c(C = 0, H = 2, N = 0, O = 1, S = 0, P = 0) 
    
    # get the number of exchangeable amides
    num_exch_sites <- exchangeableAmides(sequence)
    
    ## run simulation
    sim <- replicate(1000, expr = simulationHDX(elements = elements,
                                                    incorp = incorp,
                                                    num_exch_sites = num_exch_sites,
                                                    charge = charge))

    ## bin ions
    b <- seq(from = min(sim) - (1 / (2 * charge)),
             to = max(sim) + 1, 
             by = 1 / charge)
    bins <- cut(sim, breaks = b)
    mz <- round(tapply(sim, bins, mean), digits = 2)
    mz <- as.numeric(mz[!is.na(mz)])
    
    intensity <- as.vector(table(bins))
    intensity <- as.numeric(intensity[intensity != 0])
    
    # generate spectra object
    spec <- DataFrame(
        msLevel = c(1L),
        charge = charge,
        sequence = sequence)
    
    spec$mz <- list(mz)
    spec$intensity <-  list(intensity)
    
    # construct spectra object
    sps <- Spectra(spec)
    
    return(sps)
    
}

##' Computes the number of exchangeable amides based on the sequnece
##' @title Compute exchangeable amides.
##' @param sequence The sequence of the peptide
##' @return Returns a numeric indicating the number of exchangeable amides
##' @md
##' 
##' @rdname hdx-distributions
exchangeableAmides <- function(sequence) {
    
    n <- length(sequence)
    x <- vector(mode = "numeric", length = n)
    
    for(i in 1:n) {
        seq_vector <- strsplit(as.character(sequence[i]), split = "")[[1]]
        x[i] <- length(na.omit(sub("P", NA, seq_vector))) - 2
    }
    
    return(x)	
}

##' This function performs simulation for deuterium incorporation
##' @title Simulate possible hdx m/z's
##' @param elements The elements/residues used in the simulation
##' @param incorp The desired deuterium incoperation level
##' @param num_exch_sites The number of exchangeable amides
##' @param charge The charge of the Spectra
##' @return A vector of mz indicating the isotopic expansion of a hdx spectra
##' @md
##' 
##' @rdname hdx-distributions
simulationHDX <- function(elements,
                          incorp,
                          num_exch_sites,
                          charge) {
    
    mz <- vector(mode = "numeric")
    
    ## mass of carbons
    mc <- sum(sample(c(12.0000000, 13.0033548378), 
                     size = elements["C"], 
                     replace = TRUE, 
                     prob = c(0.9893, 0.0107)))   
    
    ## mass of unexchanged (non-backbone) hydrogens
    mh <- sum(sample(c(1.0078250321, 2.0141017780), 
                     size = elements["H"] - num_exch_sites,
                     replace = TRUE, 
                     prob = c(0.999885, 0.000115)))
    
    ## mass of exchanged (backbone) hydrogens
    md <- sum(sample(c(1.0078250321, 2.0141017780),
                     size = num_exch_sites,
                     replace = TRUE,
                     prob = c(1 - incorp, incorp)))
    
    ## mass of nitrogens
    mn <- sum(sample(c(14.0030740052, 15.0001088984), 
                     size = elements["N"], 
                     replace = TRUE, 
                     prob = c(0.99632, 0.00368)))
    
    ## mass of oxygens
    mo <- sum(sample(c(15.9949146221, 16.99913150, 17.9991604), 
                     size = elements["O"], 
                     replace = TRUE, 
                     prob=c(0.99757, 0.00038, 0.00205)))
    
    ## mass of sulfers
    ms <- sum(sample(c(31.97207069, 32.97145850, 33.96786683, 35.96708088), 
                     size=elements["S"], 
                     replace = TRUE, 
                     prob = c(0.9493, 0.0076, 0.0429, 0.0002)))
    
    ## mass of charge
    mch <- sum(sample(c(1.0072764522, 2.0135531981), 
                      size = charge, 
                      replace = TRUE, 
                      prob=c(0.999885, 0.000115)))
    
    ## m/z of molecule
    if(charge != 0) mz <- sum(mc, mh, md, mn, mo, ms, mch) / charge 
    else mz <- mz <- sum(mc, mh, md, mn, mo, ms, mch)		
    return(mz)	
}