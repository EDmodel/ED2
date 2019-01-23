#==========================================================================================#
#==========================================================================================#
#     This function generates the projection table for the UTM zone.                       #
#------------------------------------------------------------------------------------------#
utm.projcs <<- function(utm.zone,south){
   #----- Find the values that change depending on the UTM zone. --------------------------#
   merid          = sprintf("%6.1f",-177 + (utm.zone - 1) * 6)
   nors           = ifelse(south,"S","N")
   utm.epsg       = 32600 + 100 * south + utm.zone
   false.northing = sprintf("%.1f",10000000 * south)
   #---------------------------------------------------------------------------------------#



   #----- Populate the strings. -----------------------------------------------------------#
   labzone            = paste0("\"WGS_1984_UTM_Zone_",utm.zone,nors,"\"")
      gcslab          = "GCS_WGS_1984"
      spheroid        = "SPHEROID[\"WGS_1984\",6378137,298.257223563]"
      datum           = paste0("DATUM[\"D_WGS_1984\",",spheroid,"]")
      primem          = "PRIMEM[\"Greenwich\",0]"
      unit            = "UNIT[\"Degree\",0.017453292519943295]"
   geogcs             = paste0("GEOGCS["
                              ,paste(gcslab,spheroid,datum,primem,unit,sep=",")
                              ,"]"
                              )#end paste0
   projection         = "PROJECTION[\"Transverse_Mercator\"]"
   latitude.of.origin = "PARAMETER[\"latitude_of_origin\",0]"
   central.meridian   = paste0("PARAMETER[\"central_meridian\", ",trim(merid),"]")
   scale.factor       = "PARAMETER[\"scale_factor\", 0.9996]"
   false.easting      = "PARAMETER[\"false_easting\", 500000.0]"
   false.northing     = paste0("PARAMETER[\"false_northing\",",false.northing,"]")
   unit               = "UNIT[\"m\", 1.0]"
   #---------------------------------------------------------------------------------------#



   #------ Build the string ready to write to the projection file. ------------------------#
   projcs = paste0( "PROJCS["
                  , paste( labzone
                         , geogcs
                         , projection
                         , latitude.of.origin
                         , central.meridian
                         , scale.factor
                         , false.easting
                         , false.northing
                         , unit
                         , sep = ","
                         )#end paste
                  , "]"
                  )#end paste0
   return(projcs)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
