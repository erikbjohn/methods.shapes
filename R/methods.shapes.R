#' \code{methods.shapes} package
#'
#' methods.shapes
#'
#' See the README on
#'
#' @docType package
#' @name methods.shapes
#' @importFrom dplyr %>%
#' @importFrom data.table :=
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' @title clean.shape
#'
#' @description Reprojects shape file to proj.env and turns factor data columns to strings.
#' @param x spatialPolygonsDataFrame
#' @param proj.env.name projection environment Example: denver
#' @keywords transform, project, spatialPolygonsDataFrame
#' @export
#' @import rgdal
#'      sp
clean.shape <- function(x, proj.env.name=NULL){
    proj.env <- projection.value(proj.name = proj.env.name)
    pc <- spTransform(x, CRS(proj.env))
    #pc <- clgeo_Clean(pc) # Useful for holes
    #pc <- gBuffer(pc, byid=TRUE, width=0)
    # For raster union
    if ('id' %in% names(pc@data)){
        pc@data$id <- as.integer(pc@data$id)
    }
    ids.factor <- which(sapply(pc@data, class) == 'factor')
    if (length(ids.factor) > 0){
        for (id.factor in ids.factor){
            pc@data[, id.factor] <- as.character(pc@data[,id.factor])
        }
    }
    return(pc)
}
#' @title explode.shape.union
#'
#' @description When the spatialPolygon has only a length of 1, this will bring everything up so shapes are accessible.
#' @param b spatialPolygonsDataFrame
#' @keywords transform, project, spatialPolygonsDataFrame
#' @export
#' @import rgdal
#'     sp
#'     methods
explode.shape.union <- function(b){
    polys <- b@polygons[[1]]@Polygons
    pl <- vector("list", length(polys))
    for (i in 1:length(polys)) { pl[i] <- Polygons(list(polys[[i]]), i) }
    b.spolys <- SpatialPolygons(pl)
    row.ids <- sapply(slot(b.spolys, "polygons"), function(i) slot(i, "ID"))
    b.exploded <- SpatialPolygonsDataFrame(b.spolys, data.frame(FID=as.numeric(row.ids)))
    return(b.exploded)
}
#' @title export.ogr
#'
#' @description Simple function to export as an ogr file.
#' @param shp spatialPolygonsDataFrame
#' @param shp.name output .shp file root name
#' @keywords export.ogr
#' @export
#' @import rgdal
#'     sp
export.ogr <- function(shp, shp.name){
    if(!dir.exists('EsriShapes')) dir.create('EsriShapes')
    writeOGR(shp, ".", paste0("EsriShapes/", shp.name), driver="ESRI Shapefile")
}
#' @title line.orthogonal
#'
#' @description Finds orthogonal line for the streetview direction aiming.
#' @param midpoint spatialPointsDataframe
#' @param shp.line spatialLinesDataframe
#' @param proj.env.name Projection envrionment name Example: denver
#' @keywords streetview, orthogonal line
#' @export
#' @import rgdal
#'     sp
#'     geosphere
#' @importFrom data.table as.data.table
#'     data.table rbindlist
line.orthogonal <- function(midpoint, shp.line, proj.env.name=NULL){
    proj.env <- projection.value(proj.name = proj.env.name)
    proj.wgs84 <- projection.value(proj.name = 'wgs84')
    midpoint.wgs <- spTransform(midpoint, CRS(proj.wgs84))
    line.wgs <- spTransform(shp.line,  CRS(proj.wgs84))
    line.ends <- slot(line.wgs@lines[[1]]@Lines[[1]], 'coords')
    line.bearing <- bearing(line.ends[1,], line.ends[2,])
    m.ray <- round(as.matrix(rbindlist(lapply(c(line.bearing-90, line.bearing+90),
                                              function(x) as.data.table(destPoint(midpoint.wgs, x, d=140))))),8)
    line.ray <- Lines(list(Line(cbind(m.ray[,1], m.ray[,2]))), ID='a')
    sp.line.ray <- SpatialLines(list(line.ray), proj4string = CRS(proj.wgs84))
    df.line.ray <- data.frame(ID=row.names(sp.line.ray))
    row.names(df.line.ray) <- row.names(sp.line.ray)
    sp.line.ray <- SpatialLinesDataFrame(sp.line.ray, data=df.line.ray)
    # Export to qgis and check
    #line.wgs@data$ID <- row.names(line.wgs@data)
    #row.names(line.wgs) <- line.wgs@data$ID
    l <- list()
    l$line <- spTransform(sp.line.ray, CRS(proj.env))
    return(l$line)
}
#' @title shp.2.lines
#'
#' @description Used by streetview to break apart a shapefile to its component outline.
#' @param shp spatialPolygonsDataFrame
#' @param proj.env.name Projection envrionment name Example: denver
#' @keywords streetview, shape.2.lines, shape to line
#' @export
#' @import rgdal
#'     sp
shp.2.lines <- function(shp, proj.env.name=NULL){
    proj.env <- projection.value(proj.name = proj.env.name)
    sl <- as(shp, 'SpatialLines')
    coord.1 <- slot(sl@lines[[1]]@Lines[[1]], 'coords')
    coord.1.length <- dim(coord.1)[1]
    coord.1 <- rbind(coord.1, coord.1[1,])
    Ls <- list()
    for (i in 1:coord.1.length){
        L1 <- Line(coord.1[i:(i+1),1:2])
        Ls[[i]] <- Lines(list(L1), ID = as.character(i))
    }
    SL1 = SpatialLines(Ls)
    sp::proj4string(SL1) <- proj.env
    SL1 <- spTransform(SL1, CRS(proj.env))
    df <- data.frame(ID = row.names(SL1))
    SL1 <- SpatialLinesDataFrame(SL1, df)
    SL1@data$LineLength <- sapply(SL1@lines, function(x) LineLength(Line(slot(slot(x, 'Lines')[[1]], 'coords'))))
    SL1 <- SL1[SL1@data$LineLength>0,]
    return(SL1)
}
#' @title shapes.check.identical
#'
#' @description Check to see if there are any identical shapes and consolidate.
#' @param DT SpatialPolygonsDataframe
#' @keywords consolidate shapes
#' @export
#' @import rgdal
#'     sp
shapes.check.identical <- function(DT){
    DT@data$id <-  sapply(slot(DT, "polygons"), function(x) slot(x, "ID"))
    poly.slot <-  sapply(slot(DT, "polygons"), function(x) slot(x, "Polygons"))
    coords <- sapply(poly.slot, function(x) slot(x, "coords"), simplify=FALSE)
    i.coord <- 1
    n.coord <- length(coords)
    coord.unique.locs <- 1
    while (i.coord < n.coord){
        coord.base <- coords[i.coord]
        ident <- TRUE
        while(ident == TRUE & i.coord < n.coord){
            i.coord <- i.coord +1
            ident <- identical(coord.base, coords[i.coord])
        }
        if(ident==TRUE){
            coord.unique.locs
        } else {
            coord.unique.locs <- c(coord.unique.locs, i.coord)
        }
    }
    DT.ids <-  DT@data$id[coord.unique.locs]
    DT <- DT[DT@data$id %in% DT.ids,]
    return(DT)
}
#' @title shapes.coords2points
#'
#' @description Convert data table with lats and longs to spatial points.
#' @param DT data.table$long, data.table$lat
#' @param proj.env.name Projection envrionment name (Example: denver)
#' @keywords consolidate shapes
#' @export
#' @import rgdal
#'     sp
#'     pkg.data.paths
#' @importFrom dplyr select
shapes.coords2points <- function(DT, proj.env.name = NULL){
    # Input file:
    ## DT: data.table with lat, long
    # Output file:
    ## spatialPointsDataFrame
    proj.env <- projection.value(proj.name = proj.env.name)
    proj.wgs84 <- projection.value(proj.name = 'wgs84')
    lat <- NULL; long <- NULL
    coords <- cbind(Longitude = as.numeric(as.character(DT$long)),
                    Latitude = as.numeric(as.character(DT$lat)))
    Location.pts <- SpatialPointsDataFrame(coords, dplyr::select(DT, -lat, -long),
                                           proj4string = CRS(proj.wgs84))
    return(spTransform(Location.pts, CRS(proj.env)))
}
#' @title shapes.extent
#'
#' @description Returns the extent of any shapefile in WGS84 format.
#' @param shp spatialPolygonsDataFrame
#' @keywords consolidate shapes
#' @export
#' @import rgdal
#'     sp
#' @importFrom raster extent
shapes.extent <- function(shp){
    proj.wgs84 <- projection.value(proj.name='wgs84')
    shp <- spTransform(shp, CRS(proj.wgs84))
    shapes.extent <- extent(shp)
    return(shapes.extent)
}
#' @title shapes.points.2.buffer
#' @description Create a 1000 foot buffer around spatial points/parcel intersection this is currently only used for drug treatment locations.
#' @param sp spatialPointsDataframe for lat/long
#' @param shp spatialPolygonsDataFrame the intersecting polygon
#' @keywords points 2 buffer
#' @export
#' @import rgdal
#'     sp
#'     rgeos
#' @importFrom stats na.omit
shapes.points.2.Buffer <- function(sp, shp){
    BufferB <- gBuffer(sp, byid = TRUE, width = 305)
    sp.shp <- over(sp, shp[,"schednum"])
    sp.shp <- na.omit(unique(sp.shp))
    if (nrow(sp.shp)>0){
        BufferA <- gBuffer(shp[shp@data$schednum %in% sp.shp$schednum,],
                           byid = TRUE, width=305)
        #writeOGR(BufferA, '../Maps/', "BufferA", driver="ESRI Shapefile")
        # B. If not intersected, create 1000 meter buffer around point
        BufferFinal <- gUnion(BufferA, BufferB)
    } else {
        BufferFinal <- BufferB
    }
    # C. Union A and B
    # Store map
    return(BufferFinal)
}
#' @title projection.value
#'
#' @description Returns projection details from l.pkg
#' @param path.root root path Default: ~/Dropbox/pkg.data
#' @param proj.name projection name Example: 'wgs84'
#' @param list.only echos list of all possible projections
#' @keywords projection name
#' @export
#' @import pkg.data.paths
#'     knitr
#' @importFrom data.table as.data.table
#'     data.table rbindlist
#'
#' @importFrom utils str
projection.value <- function(path.root = NULL, proj.name = NULL, list.only= FALSE){
    file.name <- NULL; sys.path <- NULL
    pkg.paths <- pkg.data.paths::paths(path.root, str.pkg.name = 'methods.shapes')
    l.pkg.path <- pkg.paths[file.name=='l.pkg.rdata', sys.path]
    load(l.pkg.path)
    # Check if wgs84 already exists
    check.wgs84 <- length(l.pkg[names(l.pkg)=='wgs84']) > 0
    if (check.wgs84 == 0){
        l.proj <- '+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
        names(l.proj) <- 'wgs84'
        l.pkg <- c(l.pkg, l.proj)
        save(l.pkg, file=l.pkg.path)
    }
    if (list.only == TRUE){ # Just return list of possible projections
        dt.pkg <- rbindlist(mapply(function(x, y) data.table(proj.name = x, projection.value = y), names(l.pkg), l.pkg, SIMPLIFY = FALSE))
        knitr::kable(dt.pkg)
    } else {
        if (is.null(proj.name)){
            proj.name <- 'denver'
            cat('No projection name set. Using', proj.name)
        }
        l.proj <- l.pkg[names(l.pkg)==proj.name]
        while (length(l.proj)==0){
            cat(paste0('Projection ', proj.name, ' not found. Available projections are: \n'))
            cat(paste0(str(l.pkg)))
            proj.continue <- readline(prompt = cat(paste0('Assign name ',  proj.name, ' to a projection from the list (Example: wgs84) or type "1" for new projection')))
            if (proj.continue != 1){ # Assign additional name to existing projection
                l.proj <- l.pkg[names(l.pkg)==proj.continue]
                if (length(l.proj)>0){
                    names(l.proj) <- proj.name
                    l.pkg <- c(l.pkg, l.proj)
                }
            } else { # New projection
                l.proj <- readline(prompt = cat('Enter projection value. (Example: +init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0)'))
                names(l.proj) <- proj.name
                l.pkg <- c(l.pkg, l.proj)
            }
            save(l.pkg, file=l.pkg.path)
        }
        return(l.proj[[1]])
    }
}

# shapes.points.2.shape <- function(shapes.points){
#     shapes.points2shape <- gBuffer(shapes.points, byid = TRUE, width=5)
#     return(shapes.points2shape)
# }
# funShp.intersect <- function(){
#     print('Using full parcel outline intersect from QGIS')
#     if (!exists('shp.intersect')){
#         if (file.exists(shp.intersect.location)){
#             load(shp.intersect.location)
#         } else {
#             shp.intersect <- readOGR(dsn = 'CleanData/EsriShapes/',
#                                      layer = 'shapes.parcels.outlines', verbose = FALSE)
#             shp.intersect@data$parcel.id <- shp.intersect@data$parcel_id
#             shp.intersect@data$outline.id <- shp.intersect@data$outline_id
#             proj.orig <- proj4string(shp.intersect) # Original projections
#             proj.new <- proj.env
#             shp.intersect <- shp.intersect %>%
#                 spTransform(CRS(proj.new))
#             save(shp.intersect, file=shp.intersect.location)
#         }
#     }
#     return(shp.intersect)
# }
# funShapes.DT <- function(DT, DT.name, debug.check){
#     # Purpose: Assign data to actual shapes (costar.rent -> shapes.costar.rent)
#     # Assumes: Have already assigned to parcel/outline/point the funLocation.assign.rawData
#     # In:
#     ## DT - data.table such as costar.rent with locations already found. (costar.rent)
#     ## DT.name - data.table name as string (costar.rent)
#     # Out:
#     ## shapes.DT - shape file with all location.ids assigned to shapes
#     if(missing(debug.check)) debug.check <- FALSE
#     shapes.DT.location <- paste0('shapes.', DT.name, '.location')
#     shapes.DT.name <- paste0('shapes.', DT.name)
#     if (!exists(shapes.DT.name)){
#         if (file.exists(get(shapes.DT.location))){
#             shapes.DT <- readRDS(get(shapes.DT.location))
#         } else {
#             # parcels
#             DT.parcel.ids <- unique(DT[location.source=='parcels', location.id])
#             shapes.parcels <- funShapes.parcels()
#             shapes.DT.parcels <- shapes.parcels[shapes.parcels@data$parcel.id %in% DT.parcel.ids, ]
#             shapes.DT.parcels <- gBuffer(shapes.DT.parcels, byid=TRUE, width=0)
#             shapes.DT.parcels@data$parcel.area <- round(gArea(shapes.DT.parcels, byid=TRUE),2)
#             if (debug.check==TRUE){ # Examine the parcel areas (in sq feet)
#                 hist(shapes.DT.parcels@data$area)
#                 # The big ones should be erased and we should try to geocode to building outline (focus on mmed?)
#                 # Also, should do this in the funLocation.assign function (check right away after string match)
#                 big.DT.parcels <- shapes.DT.parcels[shapes.DT.parcels@data$area > 50000,]
#                 writeOGR(big.DT.parcels, '.', 'CleanData/EsriShapes/mmed/parcels.big', driver='ESRI Shapefile')
#                 writeOGR(shapes.parcels[shapes.parcels@data$parcel.id == '261123000001',], '.',
#                          'CleanData/EsriShapes/mmed/parcels.really.big', driver='ESRI Shapefile')
#             }
#
#             # outlines
#             DT.outlines.ids <- unique(DT[location.source=='outlines', location.id])
#             shapes.outlines <- funShapes.outlines()
#             shapes.DT.outlines <- shapes.outlines[shapes.outlines@data$outline.id %in% DT.outlines.ids, ]
#             shapes.DT.outlines@data$outline.area <- round(gArea(shapes.DT.outlines, byid=TRUE),2)
#
#             # points
#             DT.points.ids <- unique(DT[location.source=='points', location.id])
#             points.address <- funPoints.address()
#             points.DT.address <- points.address[location.id %in% DT.points.ids, .(location.id, long, lat)]
#             coords <- cbind(Longitude = as.numeric(as.character(points.DT.address$long)),
#                             Latitude = as.numeric(as.character(points.DT.address$lat)))
#             points.DT  <- SpatialPointsDataFrame(coords, dplyr::select(points.DT.address, -long, -lat),
#                                                  proj4string = CRS("+init=epsg:4326"))
#             points.DT <-  spTransform(points.DT, CRS(proj.env))
#             shapes.DT.points <- gBuffer(points.DT, byid=TRUE, width=1)
#             shapes.DT.points@data$point.id <- shapes.DT.points@data$location.id
#             shapes.DT.points@data$location.id <- NULL
#
#             if(debug.check==TRUE){
#                 plot(shapes.DT.parcels[shapes.DT.parcels@data$parcel.id=='162871601',], col='green')
#                 plot(shapes.DT.outlines[shapes.DT.outlines@data$outline.id=='352088',], col='grey', add=TRUE)
#             }
#
#             # If parcels exist, then assign to shapes.DT
#             if (length(shapes.DT.parcels)>0){
#                 shapes.DT <- shapes.DT.parcels
#
#                 # Dissolve parcels to underlying building outlines if exist (use full outline shapefile)
#                 print('Intersecting selected parcels with full building outlines')
#                 shapes.DT.intersect <- intersect(shapes.DT.parcels, shapes.outlines)
#                 shapes.DT.intersect.ids <- sapply(slot(shapes.DT.intersect, "polygons"), function(x) slot(x, "ID"))
#
#                 if (debug.check==TRUE){
#                     # Not sure I really want to dissolve right here.
#                     writeOGR(shapes.DT.intersect, '.', 'CleanData/EsriShapes/mmed/shapes.mmed.intersect', driver='ESRI Shapefile')
#                     shapes.DT.intersect.dissolve <- gUnaryUnion(shapes.DT.intersect, id = shapes.DT.intersect@data$parcel.id)
#                     DT.dissolve.parcel.ids <- sapply(slot(shapes.DT.intersect.dissolve, "polygons"), function(x) slot(x, "ID"))
#                     newdata <- data.frame(data.id = 1:length(DT.dissolve.parcel.ids), parcel.id = DT.dissolve.parcel.ids)
#                     row.names(newdata) <- row.names(shapes.DT.intersect.dissolve)
#                     shapes.DT.intersect.dissolve <- SpatialPolygonsDataFrame(shapes.DT.intersect.dissolve, data=newdata)
#                     polygons.n <- sapply(shapes.DT.intersect.dissolve@polygons, function(p) length(p@Polygons))
#                     polygons.fix <- which(polygons.n > 1)
#                     shapes.DT.intersect.dissolve[polygons.fix,]
#                     # Find those with multiple polygons
#                     writeOGR(shapes.DT.intersect.dissolve, '.', 'CleanData/EsriShapes/mmed/outlines.parcels.intersect.dissolve', driver='ESRI Shapefile')
#
#                     # Areas defined in sqfeet
#                     # http://www.canorml.org/news/A_SUMMARY_OF_THE_MEDICAL_MARIJUANA_REGULATION_AND_SAFETY_ACT
#                     # CULTIVATION SIZE LIMITATIONS  The maximum allowable size is 1 acre (43,560 sq ft) outdoors (Type 3)
#                     # or 22,000 sq ft indoors (Type 3A and 3B licenses).
#                     # The DFA is directed to limit the number of Type 3, 3A and 3B licenses.  (AB 243, 19332(g)).
#
#                     DT.dissolve.area <- sapply(slot(shapes.DT.intersect.dissolve, "polygons"), function(x) slot(x, "area"))
#                     # Relatively large - good
#                     plot(shapes.DT.parcels[shapes.DT.parcels@data$parcel.id=='160512290',], col='red')
#                     plot(shapes.DT.intersect[shapes.DT.intersect@data$parcel.id=='160512290',], col='yellow', add=TRUE)
#
#                     # Biggest one
#                     biggest.location.id <- DT.dissolve.parcel.ids[which.max(DT.dissolve.area)]
#                     biggest.address <- mmed[location.id==DT.dissolve.parcel.ids[228], completeAddress][1]
#                     plot(shapes.DT.parcels[shapes.DT.parcels@data$parcel.id==biggest.location.id,], col='red')
#                     plot(shapes.DT.intersect.dissolve[which.max(DT.dissolve.area),], col='blue', add=TRUE)
#                     x$area_sqkm <- area(x) / 1000000
#                 }
#             }
#
#             # If outlines 1) no shapes.DT exists ) yes shapes.DT
#             if (length(shapes.DT.outlines)>0){ # yes outlines
#                 if (!exists('shapes.DT')) { # No parcels
#                     shapes.DT <- shapes.DT.outlines
#                 } else { # yes shapes.DT
#                     shapes.DT <- union(shapes.DT.outlines, shapes.DT)
#                 }
#             }
#
#             # If points 1) no shapes.DT exists ) yes shapes.DT
#             if (length(shapes.DT.points)>0){ # yes outlines
#                 if (!exists('shapes.DT')) { # No parcels
#                     shapes.DT <- shapes.DT.points
#                 } else { # yes shapes.DT
#                     shapes.DT <- union(shapes.DT.points, shapes.DT)
#                 }
#             }
#
#             # Save individual and union to shapefiles
#             ogr.directory <- paste0('CleanData/EsriShapes/', DT.name)
#             if (!dir.exists(ogr.directory)) dir.create(ogr.directory)
#             ogr.location.union <- paste0(ogr.directory,'/shapes.', DT.name)
#             ogr.location.parcels <- paste0(ogr.directory,'/shapes.', DT.name, '.parcels.new')
#             ogr.location.outlines <- paste0(ogr.directory,'/shapes.', DT.name, '.outlines.new')
#             ogr.location.points <- paste0(ogr.directory,'/shapes.', DT.name, '.points.new')
#
#             writeOGR(shapes.DT, '.', ogr.location.union, driver='ESRI Shapefile')
#             writeOGR(shapes.DT.parcels, '.', ogr.location.parcels, driver='ESRI Shapefile')
#             writeOGR(shapes.DT.outlines, '.', ogr.location.outlines, driver='ESRI Shapefile')
#             writeOGR(shapes.DT.points, '.', ogr.location.points, driver='ESRI Shapefile')
#
#             saveRDS(shapes.DT, file=get(shapes.DT.location))
#         }
#     }
#     return(shapes.DT)
# # }
# shapes.get.kml.attributes <- function(kmlfile, ignoreAltitude = FALSE){
#     # Added by me
#     re <- 'Type :.+?<\\/p>'
#     kml <- paste(readLines(kmlfile, encoding = "UTF-8"), collapse = " ")
#     mtchs <- gregexpr(re, kml)[[1]]
#     attr(mtchs, "match.length")
#     schType <- vector()
#     for (i in 1:(length(mtchs))) {
#         temp <-substr(kml, mtchs[i], (mtchs[i] + attr(mtchs, "match.length")[i]))
#         schType <- append(schType,temp)
#     }
#
#
#     re <- 'District.+?<br\\/>'
#     mtchs <- gregexpr(re, kml)[[1]]
#     attr(mtchs, "match.length")
#     district <- vector()
#     for (i in 1:(length(mtchs))) {
#         temp <-substr(kml, mtchs[i], (mtchs[i] + attr(mtchs, "match.length")[i]))
#         district <- append(district,temp)
#     }
#
#     attributes <- cbind(schType, district)
# }
# shapes.union.all.2 <- function(x){
#     seq.start <- seq(1, length(x), 2)
#     seq.end <- c(seq.start[2:length(seq.start)], length(x))
#     seq.start <- seq.start+1
#     seq.start[1] <- 1
#     seqs <- data.table(seq.start, seq.end)
#     seqs <- seqs[seq.start <= seq.end & !is.na(seq.end) & !(seq.start==1 & seq.end==1)]
#     shape.osm<- list()
#     for (iSeq in 1:nrow(seqs)){
#         print(iSeq/nrow(seqs))
#         seq.range <- seqs[iSeq]
#         row.start <- seq.range$seq.start
#         row.end <- seq.range$seq.end
#         # print(iSeq)
#         # print(row.start)
#         shape.osm[[iSeq]] <- x[[row.start]]
#         seq.range <- (seq.range$seq.start+1):seq.range$seq.end
#         seq.range <- sapply(seq.range, function(x) min(x, row.end))
#         for (iSeq.sub in seq.range){
#             shape.osm[[iSeq]] <- rbind(shape.osm[[iSeq]], x[[iSeq.sub]], makeUniqueIDs=TRUE)
#         }
#     }
#     return(shape.osm)
# }
# funShapes.union.all.iter <- function(shapes.osm){
#     while (length(shapes.osm) > 1){
#         print(paste(length(shapes.osm), 'to be joined'))
#         shapes.osm <- funShapes.union.all.2(shapes.osm, nGroup = 2)
#     }
#     saveRDS(shapes.osm, file='RawData/osm/shape.osm.rds')
#
#     return(shapes.osm)
# }
