################################################################################
#' @title NEON Soil Inorganic N Concentrations and Net N Transformation Rates

#' @author Samantha Weintraub \email{sweintraub@battelleecology.org}

#' @description Calculate soil extractable inorganic nitrogen concentrations and
#' net N transformation rates for NEON L1 data. It is recommended to use the
#' neonUtilities package to download data prior to running this function.

#' @param kclInt A data frame containing soil masses and kcl volumes used in kcl
#'   extractions. Data product table name is ntr_internalLab
#' @param kclIntBlank A data frame containing information needed to link kcl
#'   extraction samples to procedural blanks. Data product table name is
#'   ntr_internalLabBlanks
#' @param kclExt A data frame containing inorganic N concentrations measured in
#'   kcl extractions and blanks. Data product table name is ntr_externalLab
#' @param soilMoist A data frame containing soil moisture values. Data product
#'   table name is sls_soilMoisture
#' @param dropConditions An optional list of sampleCondition or dataQF values for which to exclude
#'   ammonium and nitrate concentration measurements
#' @param dropFlagged An option to exclude ammonium and nitrate concentration measurements for
#'   samples with external lab quality flags. Defaults to F
#' @param keepAll An option to keep all variables and blank info used in calculations. 
#'   Defaults to F, meaning only sample information (not blanks), critical input variables, 
#'   and calculated outputs will be included in the output data frame. 
#' @return A data frame of soil inorganic N concentrations in micrograms per gram
#'   in t-initial and t-final soil samples, as well as net N transformation rates
#'   for t-final cores. Rows for blank samples are not included when keepAll = F (default).

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords soil, nitrogen, mineralization, nitrification

#' @examples
#' \dontrun{
#' Load data to R using loadByProduct (neonUtilties)
#' 
#' NTR <- loadByProduct(site = "GUAN", dpID = "DP1.10080.001", 
#' package = "basic", check.size = F)
#' 
#' out <- def.calc.ntrans(kclInt = NTR$ntr_internalLab, kclIntBlank = NTR$ntr_internalLabBlanks, 
#' kclExt = NTR$ntr_externalLab, soilMoist = NTR$sls_soilMoisture)
#' 
#' # If data is downloaded to a computer
#' 
#' df1 <- "path/to/data/ntr_internalLab"
#' df2 <- "path/to/data/ntr_internalLabBlanks"
#' df3 <- "path/to/data/ntr_ntr_externalLab"
#' df4 <- "path/to/data/sls_soilMoisture"
#'
#' out <- def.calc.ntrans(kclInt = df1, kclIntBlank = df2, kclExt =df3, soilMoist = df4,
#' dropConditions = c("deprecatedMethod", "other"), dropFlagged = T)
#' }

#' @seealso Currently none

#' @export

# changelog and author contributions / copyrights
#   Samantha Weintraub (2017-11-22)
#     original creation
#   Samantha Weintraub (2018-04-20)
#     minor updates to allow for multiple dropConditions
#   Samantha Weintraub (2019-03-29)
#     bug fixes
################################################################################

# Function
def.calc.ntrans <- function(kclInt,
                            kclIntBlank,
                            kclExt,
                            soilMoist,
                            dropConditions,
                            dropFlagged = FALSE,
                            keepAll = FALSE
){
  
  # check for missing datasets
  null.check = sapply(list(kclInt, kclIntBlank, kclExt, soilMoist), is.null)
  nullDSs = c('kclInt', 'kclIntBlank', 'kclExt', 'soilMoist')[null.check]
  if (length(nullDSs) > 0) {
    stop(paste0(paste(nullDSs, collapse = ', '), ' dataset(s) missing.'))
  }
  
  # join the internal and external lab data
  suppressWarnings(suppressMessages(combinedDF <- left_join(kclExt, kclInt, 
                                                            by = c("sampleID",
                                                                   "kclSampleID", 
                                                                   "domainID", 
                                                                   "siteID", 
                                                                   "plotID", 
                                                                   "namedLocation", 
                                                                   "collectDate"))))
  
  # set N data to NA based on sample condition or dataQF values (optional)
  if(!missing(dropConditions)) {
    # conditions and NEON quality flags in external lab data
    combinedDF$kclAmmoniumNConc[combinedDF$sampleCondition.x %in% dropConditions] <- NA
    combinedDF$kclNitrateNitriteNConc[combinedDF$sampleCondition.x %in% dropConditions] <- NA
    combinedDF$kclAmmoniumNConc[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.x)] <- NA
    combinedDF$kclNitrateNitriteNConc[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.x)] <- NA
    # conditions and NEON quality flags in internal lab data
    combinedDF$kclAmmoniumNConc[combinedDF$sampleCondition.y %in% dropConditions] <- NA
    combinedDF$kclNitrateNitriteNConc[combinedDF$sampleCondition.y %in% dropConditions] <- NA
    combinedDF$kclAmmoniumNConc[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.y)] <- NA
    combinedDF$kclNitrateNitriteNConc[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.y)] <- NA
    
    # how many values got set to NA with extneral lab condition filtering
    if (any(combinedDF$sampleCondition.x %in% dropConditions)) {
      num1 <-
        length(combinedDF$sampleID[combinedDF$sampleCondition.x %in% dropConditions])
      warning1 <-
        paste(
          'warning:',
          num1,
          'records had concentration values set to NA due to external lab sample conditions',
          sep = " "
        )
      print(warning1)
    } 
    
    # how many values got set to NA with internal lab condition filtering
    if (any(combinedDF$sampleCondition.y %in% dropConditions)) {
      num1a <-
        length(combinedDF$sampleID[combinedDF$sampleCondition.y %in% dropConditions])
      warning1a <-
        paste(
          'warning:',
          num1a,
          'records had concentration values set to NA due to NEON lab sample conditions',
          sep = " "
        )
      print(warning1a)
    } 
    
    # how many values got set to NA with external lab dataQF filtering
    if (any(grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.x))) {
      num2 <-
        length(combinedDF$sampleID[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.x)])
      warning2 <-
        paste(
          'warning:',
          num2,
          'records had concentration values set to NA due to data quality issues',
          sep = " "
        )
      print(warning2)
    } 
    
    # how many values got set to NA with internal lab dataQF filtering
    if (any(grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.y))) {
      num2a <-
        length(combinedDF$sampleID[grepl(paste(dropConditions, collapse = "|"), combinedDF$dataQF.y)])
      warning2a <-
        paste(
          'warning:',
          num2a,
          'records had concentration values set to NA due to data quality issues',
          sep = " "
        )
      print(warning2a)
    } else {
      combinedDF
    }
  }
  
  # set concentration data to NA based on ammonium or nitrate quality flags (optional)
  if(dropFlagged) {
    combinedDF$kclAmmoniumNConc[combinedDF$ammoniumNQF %in% c("1", "2")] <- NA
    combinedDF$kclNitrateNitriteNConc[combinedDF$nitrateNitriteNQF %in% c("1", "2")] <- NA
    
    # compile list of how many ammonium values got set to NA with filtering
    if (any(combinedDF$ammoniumNQF %in% c("1", "2"))) {
      num3 <-
        length(combinedDF$kclAmmoniumNConc[combinedDF$ammoniumNQF %in% c("1", "2")])
      warning3 <-
        paste(
          'warning:',
          num3,
          'records had ammonium concentrations set to NA due to the QF value',
          sep = " "
        )
      print(warning3)
    } 
    
    # compile list of how many nitrate + nitrite values got set to NA with filtering
    if (any(combinedDF$nitrateNitriteNQF %in% c("1", "2"))) {
      num4 <-
        length(combinedDF$kclNitrateNitriteNConc[combinedDF$nitrateNitriteNQF %in% c("1", "2")])
      warning4 <-
        paste(
          'warning:',
          num4,
          'records had nitrate + nitrite concentrations set to NA due to the QF value',
          sep = " "
        )
      print(warning4)
    } 
  } 
  
  # add blank info & values
  combinedDF <-  suppressWarnings(suppressMessages(combinedDF %>%
                                                     mutate(kclReferenceID = toupper(kclReferenceID),
                                                            incubationPairID = ifelse(is.na(sampleID), NA, substr(sampleID, 1, nchar(sampleID) - 9))) %>%
                                                     left_join(select(kclIntBlank, kclReferenceID, kclBlank1ID, kclBlank2ID, kclBlank3ID), 
                                                               by = "kclReferenceID") %>%
                                                     mutate(blank1NH4 = as.numeric(kclAmmoniumNConc[match(kclBlank1ID, kclSampleID)]),
                                                            blank2NH4 = as.numeric(kclAmmoniumNConc[match(kclBlank2ID, kclSampleID)]),
                                                            blank3NH4 = as.numeric(kclAmmoniumNConc[match(kclBlank3ID, kclSampleID)]),
                                                            blank1NO3 = as.numeric(kclNitrateNitriteNConc[match(kclBlank1ID, kclSampleID)]),
                                                            blank2NO3 = as.numeric(kclNitrateNitriteNConc[match(kclBlank2ID, kclSampleID)]),
                                                            blank3NO3 = as.numeric(kclNitrateNitriteNConc[match(kclBlank3ID, kclSampleID)]),
                                                            blankNH4mean = rowMeans(data.frame(blank1NH4, blank2NH4, blank3NH4), na.rm = TRUE), 
                                                            blankNO3mean = rowMeans(data.frame(blank1NO3, blank2NO3, blank3NO3), na.rm = TRUE),
                                                            kclAmmoniumNBlankCor = ifelse(as.numeric(kclAmmoniumNConc) - blankNH4mean < 0, 0, as.numeric(kclAmmoniumNConc) - blankNH4mean),
                                                            kclNitrateNitriteNBlankCor = ifelse(as.numeric(kclNitrateNitriteNConc) - blankNO3mean < 0, 0, as.numeric(kclNitrateNitriteNConc) - blankNO3mean)) %>%
                                                     left_join(select(soilMoist, sampleID, soilMoisture, dryMassFraction), by = "sampleID") %>%
                                                     mutate(soilDryMass = soilFreshMass * dryMassFraction, 
                                                            soilAmmoniumNugPerGram = kclAmmoniumNBlankCor * (kclVolume / 1000) / soilDryMass * 1000,
                                                            soilNitrateNitriteNugPerGram = kclNitrateNitriteNBlankCor * (kclVolume / 1000) / soilDryMass * 1000,
                                                            soilInorganicNugPerGram = soilAmmoniumNugPerGram + soilNitrateNitriteNugPerGram)))
  
  # count how many samples are missing moisture values
  samples <- combinedDF[!grepl("BREF", combinedDF$kclSampleID),]
  if (any(is.na(samples$soilMoisture))) {
    num5 <-
      sum(is.na(samples$soilMoisture))
    warning5 <-
      paste(
        'warning:',
        num5,
        'records were missing soil moisture values',
        sep = " "
      )
    print(warning5)
  }
  
  # create wide (cast) version of the df in order to calculate net rates with incubationPairID and nTransBoutType
  combinedDFforCast <- combinedDF %>%
    filter(!sampleID == "") 
  
  cast1 <- data.table::dcast(data.table::setDT(combinedDFforCast), incubationPairID ~ nTransBoutType,
                             value.var = c("incubationLength","soilInorganicNugPerGram", "soilNitrateNitriteNugPerGram"),
                             fun = mean, na.rm = T)
  
  # calculate net rates
  cast1 <- cast1 %>%
    mutate(netInorganicNugPerGram = soilInorganicNugPerGram_tFinal - soilInorganicNugPerGram_tInitial,
           netNitrateNitriteNugPerGram =  soilNitrateNitriteNugPerGram_tFinal - soilNitrateNitriteNugPerGram_tInitial,
           netNminugPerGramPerDay = netInorganicNugPerGram / incubationLength_tFinal,
           netNitugPerGramPerDay = netNitrateNitriteNugPerGram / incubationLength_tFinal)
  
  # attach net rates onto combined df
  combinedDF <- suppressWarnings(suppressMessages(combinedDF %>%
                                                    left_join(select(cast1, incubationPairID, netNminugPerGramPerDay, netNitugPerGramPerDay), by = "incubationPairID") %>%
                                                    mutate(netNminugPerGramPerDay = ifelse(nTransBoutType == "tInitial", NA, netNminugPerGramPerDay),
                                                           netNitugPerGramPerDay = ifelse(nTransBoutType == "tInitial", NA, netNitugPerGramPerDay))))
  
  # set NaN to NA
  combinedDF[is.na(combinedDF)] <- NA
  
  # determine whether to keep all variable or just a subset
  if (keepAll) {
    combinedDF
  } else {
    combinedDF <-
      subset(
        combinedDF,
        select = c(
          "plotID",
          "collectDate",
          "nTransBoutType",
          "sampleID",
          "incubationPairID",
          "incubationLength",
          "soilFreshMass",
          "dryMassFraction",
          "soilDryMass",
          "kclVolume",
          "kclAmmoniumNBlankCor",
          "kclNitrateNitriteNBlankCor",
          "soilAmmoniumNugPerGram",
          "soilNitrateNitriteNugPerGram",
          "netNminugPerGramPerDay",
          "netNitugPerGramPerDay"
        )
      )
    combinedDF <- combinedDF[!is.na(combinedDF$nTransBoutType), ]# Drop records that have no plotID (blanks)
    combinedDF <- combinedDF[order(combinedDF$sampleID),]
  }
  
  # Round all numeric variables to 3 digits
  combinedDF %>% mutate_if(is.numeric, round, digits=3)
}