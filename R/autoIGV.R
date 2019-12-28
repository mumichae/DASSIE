# authors: Litovchenko, Mertes
# 
# Interacting with IGV over R (AutoIGV)
#
require("data.table")
AutoIGV_DEBUG <- FALSE

#'
#' toogles the debuging mode for AutoIGV on and off
#' 
set_AutoIGV_DEBUG <- function(debugging = FALSE){
	stopifnot(debugging)
	AutoIGV_DEBUG <<- debugging
}

#' Function to produce IGV session files and screenshorts automatically
#' @param coords data frame with 6 columns: chr, start, end, maximum, 
#'        coverage.only and suffix.to.file. Maximum is the maximum of the
#'        coverage tracks to be displayd. If you want to use auto maximum, 
#'        put -1 to maximum column. Column coverage.only indicates,
#'        whatever or not only coverage tracks should be displayed: FALSE, if you want
#'        reads to be displayed too and TRUE otherwise. Column 
#'        suffix.to.file contains suffixes to be appended to files names
#' @param indents vector of the nuumber of nucleoteds to be dysplayed from the
#'        left and rigth of region of interest. Default value is c(0, 0)
#' @param bams.paths vector, containing paths to bam files. The order 
#'        matters: tracks will appear in igv according to the order of 
#'        this vector
#' @param refs.paths vector, containing paths to reference files. The order
#'        also matters
#' @param coverage.only boolean, indicating, whatever or not only coverage 
#'        tracks will be diplayed
#' @param folder.name name of the folder to put results in
#' @return folder with one IGV session and 1 png per line of coords data 
#'         frame
#' @author Maria Litovchenko
#' @export
#' @example
#' source("autoIGV.R")
#' bams.dir <- "/s/project/mitoMultiOmics/bam_sorted"
#' refs.dir <- "/s/project/mitoMultiOmics/tmp_litovchenko/igv_tracks"
#' coords <- data.frame(chr = c("chr1", "chr2"), start = c(12000, 12000),
#'                      end = c(14000, 14000), maximum = c(-1, -1), 
#'                      coverage.only = c(TRUE, TRUE), 
#'                      suffix.to.file = c("A", "B"), stringsAsFactors = F)
#' bams.paths <- c(file.path(bams.dir, "49470_exome_merged_srt.bam"), 
#'                 file.path(bams.dir, "62524_exome_merged_srt.bam"), 
#'                 file.path(bams.dir, "71724_exome_merged_srt.bam"))
#' refs.paths <- c(file.path(refs.dir, "DESeq2_probes_AgilentV5.bed"),
#'                 file.path(refs.dir, "ED_exons.bed"), 
#'                 file.path(refs.dir, "49470_DESeq2_HomoDel.bed"))
#' folder.name <- "testAutoIGV"
#' max.matrix <- data.table(chr = c("chr1", "chr2"), start = c(12000, 12000),
#'                      end = c(14000, 14000), 
#'                      x49470_exome_merged_srt.bam = c(100, 40),
#'                      x62524_exome_merged_srt.bam =  c(150, 50),
#'                      x71724_exome_merged_srt.bam = c(200, 60))
#' setnames(max.matrix, c("x49470_exome_merged_srt.bam", 
#'                        "x62524_exome_merged_srt.bam",
#'                        "x71724_exome_merged_srt.bam"),
#'                       c("49470_exome_merged_srt.bam", 
#'                         "62524_exome_merged_srt.bam",
#'                         "71724_exome_merged_srt.bam"))   
#' AutoIGV(coords, indents = c(100, 100), bams.paths, refs.paths, 
#'         max.matrix, folder.name)
AutoIGV <- function(coords, indents = c(0, 0), bams.paths, refs.paths = data.table(paths = ""[0], name = ""[0]),
                    max.matrix = data.table(), folder.name = tempdir()) {
    # TODO make sure that column names exists
    # directory, where all session files and images will be stored
    # add /, so all files, that will be created, will be saved in a directory
    dir.path <- paste0(folder.name, "/IGV_session/")
    if(!file.exists(dir.path)){
        dir.create(dir.path, TRUE, TRUE)
    }
    
    # convert bams/refs.paths to table if needed
    if(is.vector(bams.paths)){
        bams.paths <- data.table(paths = bams.paths)
        bams.paths[, names := basename(bams.paths[, paths])]
    }
    
    if(is.vector(refs.paths)){
        refs.paths <- data.table(paths = refs.paths)
        refs.paths[, names := basename(refs.paths[, paths])]
    }
    
    # print session for each entry
    apply(coords, 1, ProduceIGVSession, indents, bams.paths, refs.paths, 
          max.matrix, dir.path)
    print("Produced session files")
    
    # produce batch script
    batch.path <- file.path(dir.path, "IGV_batch.txt")
    session.paths <- paste(dir.path, coords$suffix.to.file, ".xml", sep = "")
    # create batch file
    sapply(session.paths, PrintIGVbatch, bams.paths, batch.path)
    print("Produced batch script")
  
    print("Started to take pictures. Please wait, it can take some time")
    cat("exit", file =  batch.path, append = TRUE)
  
    # run the batch file
    system2("igv", c("-b", batch.path))
    print("Finished!")
}

#' ProduceIGVSession
#' Produces session file for 1 entry of coords data frame
#' @param region.params row of coords data frame 
#' @param indents vector of the number of nucleotids to be displayed from
#'        the left and rigth of region of interest. Default value is 
#'        c(0, 0)
#' @param bams.paths vector, containing paths to bam files. The order 
#'        matters: tracks will appear in igv according to the order
#'        of this vector
#' @param refs.paths vector, containing paths to reference files. 
#'        The order also matters.
#' @param dir.path name of the folder to put results in
#' @author Maria Litovchenko
#' @return folder with one IGV session and 1 png per line of coords data frame
#' @export
#' @example
ProduceIGVSession <- function(region.params, indents = c(0, 0), bams.paths, 
                              refs.paths = data.table(paths = ""[0], name = ""[0]), 
							  max.matrix = data.table(), dir.path = tempdir()) {
    file.path <- paste(dir.path, region.params['suffix.to.file'], ".xml", 
                       sep = "")
    PrintIGVSessionHeader(region.params, indents, file.path)
    PrintIGVSessionResourses(c(bams.paths[, paths], refs.paths[, paths]), 
                             file.path)
    # in case of coverage only I need to estimate the best heigth of the 
    # tracks. I want to produce images with overall heigth <= 1200 pixels.
    # So the optimal heigth of each covarage track in case of 
    # coverage.only = T will be (1200 - 200)/number of bams. I am doing
    # -200, because we need to fit reference in the picture
    bams.length <- nrow(bams.paths)
    refs.length <- nrow(refs.paths)
    coverage.only <- ifelse(region.params[['coverage.only']] == FALSE,  
                            -1, floor((1200 - 45 * refs.length) / bams.length))
    
    # print track informations
    if (nrow(max.matrix) == 0) {
        bams.paths[, maximum := region.params[['maximum']]]
    } else {
        bams.maxs <- max.matrix[chr == region.params[['chr']] &
                    start == region.params[['start']] &
                    end == region.params[['end']]][, -3:-1, with = FALSE]
        bams.maxs <- unlist(as.data.frame(bams.maxs))
        bams.paths[, maximum := bams.maxs]
    }
    apply(bams.paths, 1, PrintIGVTrack, coverage.only, file.path)
    PrintAnnotationHeader(file.path, 45 * refs.length)
    sapply(refs.paths[, paths], PrintAnnotationTrack, file.path)
    SetHeigthAndClose(bams.paths, file.path)
	
	return(file.path)
}

#' PrintIGVSessionHeader
#' Prints session header with coordinates and genome to the session file
#' @param coords vector with 5 items: chr, start, end, coverage.only,
#'        maximum and suffix.to.file
#' @param indents vector of the number of nucleoteds to be dysplayed from the
#'        left and rigth of region of interest. Default value is c(0, 0)
#' @param file.path path to the session file to print in
#' @author Maria Litovchenko
#' @export
PrintIGVSessionHeader <- function(coords, indents = c(0, 0), file.path) {
    # make coordinates in a form, understandable by IGV
    coordinates <- paste(sep = "",
                        coords[['chr']], ":",
                        as.numeric(coords[['start']]) - as.numeric(indents[1]),
                        "-",
                        as.numeric(coords[['end']]) + as.numeric(indents[2]))
    # coordinates <- paste(sep = "",
    #                      coords[['gene_name']])
    # usual header for IGV session
    header <- paste(sep = "",
                   '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n',
                   '  <Session genome="hg19" locus="', coordinates, 
                   '" version="5">\n')
    # print to file
    cat(header, file = file.path, append = FALSE)
}

#' PrintIGVSessionResourses
#' Prints session resourses with paths to files
#' @param paths.to.files vector, containing paths to bam files. 
#'        The order matters: tracks will appear in igv according to the 
#'        order of this vector
#' @param file.path path to the session file to print in
#' @author Maria Litovchenko
#' @export
PrintIGVSessionResourses <- function(paths.to.files, file.path) {  
    # prefix, which signifies resourse string
    resourse.preffix <- '    <Resource path="'
    resourse.suffix  <- '" />'
    all.paths <- paste(resourse.preffix, paths.to.files, 
                       resourse.suffix, sep = "")    
    # combine all resources
    all.paths <- paste(all.paths, collapse = "\n", sep = "\n")
    cat('  <Resources>\n',  file = file.path, append = TRUE)
    # write it to the session file
    cat(all.paths, file = file.path, append = TRUE)
    cat('\n',  file = file.path, append = TRUE)
    cat('  </Resources>\n', file = file.path, append = TRUE)
}

#' PrintIGVTrack
#' Prints tracks of bams to the session file
#' @param bam.paths path to the bam file
#' @param file.path path to the session file to print in
#' @param coverage.only parameter, which equals to -1, if reads are also
#'        displayed and the heigth of the track otherrwise
#' @author Maria Litovchenko
#' @export
PrintIGVTrack <- function(bams.paths, coverage.only, file.path) {
    # collect needed information
    if(AutoIGV_DEBUG) message(bams.paths)
	
    sample.name <- bams.paths[['names']]
    path.to.bam <- bams.paths[['paths']]
    maxim <- bams.paths[['maximum']]
    
    max.y <- -1
    autoScale <- "true"
    if (maxim != -1) {
        max.y <- maxim
        autoScale <- "false"
     }
    track.heigth <- ifelse(coverage.only == -1, 2550, coverage.only)
    
    # create tracks
    track.xml <- paste(sep="",
                       '  <Panel height="TRACKHEIGTH" name="PANELID" width="2541">\n',
                       '    <Track altColor="0,0,178" autoScale="', autoScale, '" color="175,175,175"\n',
                       '        displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10"',
                       '        id="PATHTOSAMPLE_coverage" name="SAMPLENAME Coverage" \n',
                       '        showReference="false" snpThreshold="0.2" \n',
                       '        sortable="true" visible="true"', 
                       ifelse(coverage.only == -1, "", paste(' height="TRACKHEIGTH"', sep ='')),' >\n',
                       '      <DataRange baseline="0.0" drawBaseline="true" \n',
                       '          flipAxis="false" maximum="', max.y, '" minimum="0.0" type="LINEAR"/> \n',
                       '    </Track>\n')
    # if also reads are requested
    if (coverage.only == -1) {
        track.xml <- paste(track.xml,  
                           '    <Track altColor="0,0,178" autoScale="false" color="0,0,178" \n',
                           '        displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" \n',
                           '        id="PATHTOSAMPLE" name="SAMPLENAME" showSpliceJunctions="false" \n',
                           '        sortable="true" visible="true">',
                           '      <RenderOptions colorByTag="" colorOption="UNEXPECTED_PAIR" flagUnmappedPairs="false" \n',
                           '          groupByTag="" maxInsertSize="1000" minInsertSize="50" shadeBasesOption="QUALITY" \n',
                           '          shadeCenters="true" showAllBases="false" sortByTag="" />\n',
                           '    </Track>\n')
    }
    track.xml <- paste(track.xml, '  </Panel>\n', sep = "")
  
    # print the lines to the session file
    track.xml <- gsub("TRACKHEIGTH", track.heigth, track.xml)
    track.xml <- gsub("PATHTOSAMPLE", path.to.bam, track.xml)
    track.xml <- gsub("PANELID", paste("Panel1429518", sample(100000:999999, 1), sep = ""), track.xml)
    track.xml <- gsub("SAMPLENAME", sample.name, track.xml)
  
    cat(track.xml, file = file.path, append = TRUE)
}

#' PrintAnnotationHeader
#' Prints header of the annotation part to the session file.
#' @param file.path path to the session file to print in
#' @param heigth heigth of the track
#' @author Maria Litovchenko
#' @export
PrintAnnotationHeader <- function(file.path, heigth) {
  header <- paste('\    <Panel height="', heigth, '"', sep = "")
  header <- paste(paste(header, 'name="FeaturePanel" width="2541">'),
                  paste('        <Track altColor="0,0,178"', 
                        'autoScale="false"', 'color="0,0,178"', 
                        'displayMode="COLLAPSED"', 
                        'featureVisibilityWindow="-1"', 
                        'fontSize="10"',
                        'id="Reference sequence"', 
                        'name="Reference sequence"',
                        'sortable="false"', 'visible="true"/>'),
                  paste('        <Track altColor="0,0,178"', 
                        'autoScale="false"',
                        'clazz="org.broad.igv.track.FeatureTrack"', 
                        'color="0,0,178"',
                        'colorScale="ContinuousColorScale;0.0;308.0;255,255,255;0,0,178"',
                        'displayMode="






D"', 
                        'featureVisibilityWindow="-1"',
                        'fontSize="10"', 'height="35"', 'id="hg19_genes"',
                        'name="RefSeq Genes"', 'renderer="BASIC_FEATURE"',
                        'sortable="false"', 'visible="true"', 
                        'windowFunction="count">'),
                  paste('             <DataRange baseline="0.0"',
                        'drawBaseline="true"', 'flipAxis="false"',
                        'maximum="308.0"', 'minimum="0.0"',
                        'type="LINEAR"/>'), 
                  '        </Track>\n', sep = "\n")
  cat(header, file = file.path, append = TRUE)
}

#' PrintAnnotationTrack
#' Prints annotation tracks to the session file
#' @param path.to.ann path to the annotation file
#' @param file.path path to the session file to print in
#' @author Maria Litovchenko
#' @export
PrintAnnotationTrack <- function(path.to.ann, file.path) {
  ann.track <- paste(paste('        <Track altColor="0,0,178"',
                           'autoScale="false"', 
                           'clazz="org.broad.igv.track.FeatureTrack"',
                           'color="0,0,178"',
                           'colorScale="ContinuousColorScale;0.0;1110.0;255,255,255;0,0,178"',
                           'displayMode="COLLAPSED"', 
                           'featureVisibilityWindow="-1"',
                           'fontSize="10"', 'height="10"', 
                           'id="ANNPATH"', 'name="ANNNAME"', 
                           'renderer="BASIC_FEATURE"',
                           'sortable="false"', 'visible="true"', 
                           'windowFunction="count">'),
                     paste('            <DataRange baseline="0.0"', 
                           'drawBaseline="true"',
                           'flipAxis="false"', 'maximum="1110.0"',
                           'minimum="0.0"', 'type="LINEAR"/>'),
                           '        </Track>\n', sep = "\n")
  ann.name <- basename(path.to.ann)
  ann.track <- gsub("ANNPATH", path.to.ann, ann.track)
  ann.track <- gsub("ANNNAME", ann.name, ann.track)
  cat(ann.track, file = file.path, append = TRUE)
}

#' SetHeigthAndClose
#' Sets up equal heiths to the IGV tracks
#' @param bams.paths paths to the bams files
#' @param file.path path to the session file to print in
#' @author Maria Litovchenko
#' @export
SetHeigthAndClose <- function(bams.paths, file.path) {
  div.factors <- seq(0, 0.994, length.out = nrow(bams.paths) + 2)
  div.factors[1] <- 0.005
  div.factors <- paste(div.factors, collapse = ",")
  closing <- paste('    </Panel>',
                   paste('    <PanelLayout ', 'dividerFractions="',
                   div.factors, '"/>', sep = ""),
                   '</Session>', sep = "\n")
  cat(closing, file = file.path, append = TRUE)
}

# PrintIGVbatch
#' Prints IGV batch script file for further png profuction
#' @param bams.paths paths to the bams files
#' @param batch.path path to the batch file to print in
#' @return IGV batch template file
#' @examples 
#' @seealso 
#' @keywords benchmark
#' @author Maria Litovchenko
#' @export
PrintIGVbatch <- function(session.path, bams.paths, batch.path) {
  command <- paste('load ', session.path, "\n", 
                   'maxPanelHeight 100\n', sep = "")
  cat(command, file = batch.path)
  
  #preffix <- "collapse"
  #bams.collapse <- paste(preffix, bams.paths[,names])
  #bams.collapse <- paste(paste(bams.collapse, collapse = "\n"), "\n")
  #cat(bams.collapse, file = batch.path, append = TRUE)
  
  snapshot.name <- strsplit(basename(session.path), "\\.") [[1]][1]
  snapshot.name <- paste(snapshot.name, ".png", sep = "")
  snapshot.path <- paste(dirname(session.path), "/", snapshot.name, sep = "")
  command <- paste('snapshot', snapshot.path, "\n")
  cat(command, file = batch.path, append = TRUE)
}

