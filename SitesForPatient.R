get_args <- function() {
    suppressMessages(library(argparse))
    
    p <- ArgumentParser(formatter_class='argparse.RawTextHelpFormatter')
    p$add_argument("-f", "--freeze", type="character", nargs=1,
        default="hg18", help="hg18, etc")
    p$add_argument("-p", "--patient", type="character", nargs=1,
        required=TRUE, help="subject name, e.g. pFR03")
    p$add_argument("-t", "--time", type="character", nargs=1,
        default="all", help="time, e.g. m12, all for all timepoints")
    p$add_argument("-c", "--cell", type="character", nargs=1,
        default="all", help="cell type, e.g. Tcell, all for all cell types")
    p$add_argument("-g", "--group", type="character", nargs=1,
        #default="intsites_miseq.read", 
        #default="intsitesdb.read", 
        default="dryga-test-db", 
        help="mysql groupname")
    
    p$parse_args(commandArgs(trailingOnly=TRUE))
}

args <- get_args()
print(t(as.data.frame(args)), quote=FALSE)

## libs
libs <- c("stats", "methods", "dplyr", "RMySQL", "GenomicRanges")
null <- suppressMessages(sapply(libs, library, character.only=TRUE))

## options
options(stringsAsFactors=F)
options(dplyr.width = Inf)
#' increase output width to console width
wideScreen <- function(howWide=as.numeric(strsplit(system('stty size', intern=T), ' ')[[1]])[2]) {
    options(width=as.integer(howWide))
}
if( interactive() ) wideScreen()

## check if file exist
## permission ~/.my.cnf should be 600
## stopifnot(file.exists("~/.my.cnf"))
## stopifnot(file.info("~/.my.cnf")$mode == as.octmode("600"))

## initialize connection to database
## ~/.my.cnf must be present
null <- sapply(dbListConnections(MySQL()), dbDisconnect)
#dbConn <- dbConnect(MySQL(), default.file="~/.my.cnf", group=args$group) 
dbConn <- dbConnect(MySQL(), default.file="~/.my.cnf", group="intsitesdb.read")
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)

## get sampleInfo from specimen_management.gtsp
sql <- sprintf("SELECT * FROM specimen_management.gtsp WHERE patient = '%s'",
               args$patient)
sampleInfo <- suppressWarnings( dbGetQuery(dbConn, sql) )
stopifnot( nrow(sampleInfo)>=1 )
colnames(sampleInfo) <- tolower(colnames(sampleInfo))
stopifnot(length(unique(sampleInfo$trial))==1)

patientInfo <- (dplyr::select(sampleInfo,
                              gtsp=specimenaccnum,
                              trial,
                              patient,
                              timepoint,
                              celltype,
                              vcn,
                              sampleprepmethod,
                              seqmethod) %>%
                dplyr::filter(timepoint %in% args$time | "all" %in% args$time) %>%
                dplyr::filter(celltype %in% args$cell | "all" %in% args$cell) )

similarGTSP <- sprintf("'^(%s)'",
                       paste(paste0(patientInfo$gtsp,"-"), collapse="|"))

## get singlehit sites, gtsp, posid, estAbund
message("\nGet unique sites")
sql <- sprintf("SELECT DISTINCT *
FROM samples JOIN sites
ON samples.sampleID = sites.sampleID
JOIN pcrbreakpoints
ON pcrbreakpoints.siteID = sites.siteID 
WHERE samples.sampleName REGEXP %s AND
samples.refGenome = '%s' ;", similarGTSP, args$freeze)
message(sql)

null <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), default.file="~/.my.cnf", group=args$group) 

sites.uniq <- suppressWarnings( dbGetQuery(dbConn, sql) ) 
sites.uniq <- sites.uniq[, !duplicated(colnames(sites.uniq))]

sites.uniq <- dplyr::mutate(sites.uniq,
                            gtsp=sub("-\\d+$", "", sampleName),
                            replicate=sub(".*-", "", sampleName),
                            start=ifelse(strand=="+", position, breakpoint),
                            end=ifelse(strand=="+", breakpoint, position),
                            enzyme="SHEAR",
                            SeqType="Illumina")

sites.uniq <- dplyr::left_join(sites.uniq, patientInfo, by="gtsp")

needed <- c("gtsp", "sampleName", "replicate",
            "siteID", "timepoint", "patient", "timepointDay", "celltype",
            "VCN", "Trial", "enzyme", "SeqType",
            "chr", "strand", "start", "end", "count",
            "refGenome")

sites.uniq <- sites.uniq[, tolower(colnames(sites.uniq)) %in% tolower(needed)]

sites.uniq <- dplyr::select(sites.uniq,
                            chr,
                            start,
                            end,
                            strand,
                            count,
                            sampleName=gtsp,
                            replicate,
                            siteID,
                            Trial=trial,
                            patient,
                            timepoint,
                            celltype,
                            VCN=vcn,
                            enzyme,
                            SeqType)

sites.uniq.gr <- makeGRangesFromDataFrame(sites.uniq,
                                          seqnames.field="chr",
                                          start.field="start",
                                          end.field="end",
                                          strand.field="strand",
                                          keep.extra.columns=TRUE)

fileName <- paste(args$patient,args$time,args$cell,"RData", sep=".")
save(sites.uniq.gr, file=fileName)

message("\n", length(sites.uniq.gr), " records written to ", fileName)

#### exit ####
q()

