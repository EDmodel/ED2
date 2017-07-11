#==========================================================================================#
#==========================================================================================#
#     This function generates the projection table for the UTM zone.                       #
#------------------------------------------------------------------------------------------#
utm.projcs <<- function(utm.zone,south){
   #----- Find the values that change depending on the UTM zone. --------------------------#
   merid    = sprintf("%6.1f",-177 + (utm.zone - 1) * 6)
   utm.epsg = 32600 + 100 * south + utm.zone
   #---------------------------------------------------------------------------------------#



   #----- Populate the strings. -----------------------------------------------------------#
   labzone            = paste("\"WGS 84 / UTM Zone ",utm.zone,ifelse(south,"S","N"),"\""
                             , sep=""
                             )#end paste
   spheroid           = paste( "SPHEROID[\"WGS 84\", 6378137.0, 298.257223563, "
                             , "AUTHORITY[\"EPSG\",\"7030\"]]"
                             , sep=""
                             )#end paste
   datum              = paste( "DATUM[\"World Geodetic System 1984\", "
                             , spheroid
                             , ", AUTHORITY[\"EPSG\",\"6326\"]]"
                             , sep = ""
                             )#end paste
   primem             = "PRIMEM[\"Greenwich\", 0.0, AUTHORITY[\"EPSG\",\"8901\"]]"
   degree             = "UNIT[\"degree\", 0.017453292519943295]"
   axis.east          = "AXIS[\"Geodetic longitude\", EAST]"
   axis.north         = "AXIS[\"Geodetic latitude\", NORTH]"
   geogcs             = paste( "GEOGCS[\"WGS 84\", ",datum,", ",primem,", ",degree,", "
                             , axis.east,", ",axis.north,", AUTHORITY[\"EPSG\",\"4326\"]]"
                             , sep = ""
                             )#end paste
   projection         = "PROJECTION[\"Transverse_Mercator\", AUTHORITY[\"EPSG\",\"9807\"]]"
   central.meridian   = paste("PARAMETER[\"central_meridian\", ",trim(merid),"]",sep="")
   latitude.of.origin = "PARAMETER[\"latitude_of_origin\", 0.0]"
   scale.factor       = "PARAMETER[\"scale_factor\", 0.9996]"
   false.easting      = "PARAMETER[\"false_easting\", 500000.0]"
   false.northing     = "PARAMETER[\"false_northing\", 10000000.0]"
   unit               = "UNIT[\"m\", 1.0]"
   axis.easting       = "AXIS[\"Easting\", EAST]"
   axis.northing      = "AXIS[\"Northing\", NORTH]"
   utm.epsg           = paste("AUTHORITY[\"EPSG\",\"",utm.epsg,"\"]",sep="")
   #---------------------------------------------------------------------------------------#



   #------ Build the string ready to write to the projection file. ------------------------#
   projcs             = paste( "PROJCS[",labzone,", ",geogcs,", ",projection,", "
                             , central.meridian,", ",latitude.of.origin,", "
                             , scale.factor,", ",false.easting,", ",false.northing,", "
                             , unit,", ",axis.easting,", ",axis.northing,", ",utm.epsg,"]"
                             , sep = ""
                             )#end paste
   return(projcs)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
