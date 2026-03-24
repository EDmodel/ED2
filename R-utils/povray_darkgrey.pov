//========================================================================================//
//      polygon.pov                                                                       //
//      This file is a template for plotting ED-2.2 using povray.  The cohorts must be    //
//      appended to the end of the script, using the following standard:                  //
//      tree(dbh,pft,x,y)                                                                 //
//      PFT numbers are the same as ED.                                                   //
//----------------------------------------------------------------------------------------//



//----- Define some common settings. -----------------------------------------------------//
#version 3.6;
global_settings{assumed_gamma 1.0}
#default{ finish{ ambient 0.1 diffuse 0.9 }}
#include "colors.inc"
#include "textures.inc"
//----------------------------------------------------------------------------------------//



//----------------------------------------------------------------------------------------//
//     This macro define the plant based on DBH and PFT.                                  //
//----------------------------------------------------------------------------------------//
#macro plant(dbh,ipft,xtree,ztree)
   //---- Define the allometric traits (based on tropical trees. -------------------------//
   #local mdbh    = min(dbh,96.2578);
   #local height  = 61.7 * (1. - exp (- 0.0352 * pow(mdbh,0.694) ) );
   #local cwidth  = sqrt( 1.12573 * pow(mdbh,1.05212) / pi);
   #local clength = 0.3106775 * pow(height,1.098);
   #local c23     = height - 1/3 * clength;
   #local crown   = height - 0.5 * clength;
   #local c13     = height - 2/3 * clength;
   #local bole    = height - clength;
   #local rbh     = 0.01 * mdbh;
   //-------------------------------------------------------------------------------------//



   //-------------------------------------------------------------------------------------//
   //      Define the colour of the PFTs.                                                 //
   //-------------------------------------------------------------------------------------//
   #switch (ipft)
      #case (1)
         //------ C4 Grass. --------------------------------------------------------------//
         #local pftcol   = pigment { color rgb <1.00,0.84,0.00> };
         #local pftshape = 1;
         #break
      #case (2)
         //------ Early tropical. --------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.50,1.00,0.00> };
         #local pftshape = 1;
         #break
      #case (3)
         //------ Mid tropical. ----------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.27,0.55,0.00> };
         #local pftshape = 1;
         #break
      #case (4)
         //------ Late tropical. --------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.00,0.31,0.00> };
         #local pftshape = 1;
         #break
      #case (5)
         //------ Temperate C3 grass. ----------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.67,0.51,1.00> };
         #local pftshape = 1;
         #break
      #case (6)
         //------ Northern pine. ---------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.00,0.75,1.00> };
         #local pftshape = 2;
         #break
      #case (7)
         //------ Southern pine. ---------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.28,0.82,0.80> };
         #local pftshape = 2;
         #break
      #case (8)
         //------ Late conifer. ----------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.15,0.25,0.55> };
         #local pftshape = 2;
         #break
      #case (9)
         //------ Early hardwood. --------------------------------------------------------//
         #local pftcol   = pigment { color rgb <1.00,0.55,0.00> };
         #local pftshape = 1;
         #break
      #case (10)
         //------ Mid hardwood. ----------------------------------------------------------//
         #local pftcol   = pigment { color rgb <1.00,0.27,0.00> };
         #local pftshape = 1;
         #break
      #case (11)
         //------ Late hardwood. ---------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.70,0.13,0.13> };
         #local pftshape = 1;
         #break
      #case (12)
         //------ C3 Crop. ---------------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.33,0.10,0.55> };
         #local pftshape = 1;
         #break
      #case (13)
         //------ C3 pasture. ------------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.75,0.24,1.00> };
         #local pftshape = 1;
         #break
      #case (14)
         //------ C4 Crop. ---------------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.72,0.53,0.04> };
         #local pftshape = 1;
         #break
      #case (15)
         //------ C4 pasture. ------------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.94,0.90,0.55> };
         #local pftshape = 1;
         #break
      #case (16)
         //------ C3 grass. --------------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.80,0.75,0.44> };
         #local pftshape = 1;
         #break
      #case (17)
         //------ Araucaria. -------------------------------------------------------------//
         #local pftcol   = pigment { color rgb <0.31,0.58,0.80> };
         #local pftshape = 2;
         #break
   #end // switch (ipft)
   //-------------------------------------------------------------------------------------//



   //-------------------------------------------------------------------------------------//
   //     Decide what to do based on whether this is a flowering plant or a conifer.      //
   //-------------------------------------------------------------------------------------//
   #switch(pftshape)
      #case(1)
         //----- Flowering plant, use an ellipsoid. --------------------------------------//
         merge{
            //----- The stem. ------------------------------------------------------------//
            cylinder{ <0,0,0>, <0,1,0>, 1.0
                      texture{ pigment{ color rgb <0.64,0.32,0.16>}
                               finish { phong 1 }
                             } // end of texture
                      scale <rbh,height,rbh>
                      translate<xtree,0,ztree>
                    } // end of sphere
            //----------------------------------------------------------------------------//



            //----- First crown. ---------------------------------------------------------//
            sphere{ <0,0,0>, 1.0
                    texture{ pigment{ pftcol  }
                             finish { phong 1 }
                           } // end texture
                    scale <0.25*cwidth,2*clength/3,0.25*cwidth>
                    translate<xtree,c13,ztree>
                  } // end sphere
            //----- Second crown. --------------------------------------------------------//
            sphere{ <0,0,0>, 1.0
                    texture{ pigment{ pftcol  }
                             finish { phong 1 }
                           } // end texture
                    scale <0.5*cwidth,2*clength/3,0.5*cwidth>
                    translate<xtree,crown,ztree>
                  } // end sphere
            //----- Third crown. ---------------------------------------------------------//
            sphere{ <0,0,0>, 1.0
                    texture{ pigment{ pftcol  }
                             finish { phong 1 }
                           } // end texture
                    scale <cwidth,2*clength/3,cwidth>
                    translate<xtree,c23,ztree>
                  } // end sphere
            //----------------------------------------------------------------------------//
         } // end merge
         //-------------------------------------------------------------------------------//
         #break
      #case(2)
         //----- Conifer, use a cone... --------------------------------------------------//
         //-------------------------------------------------------------------------------//
         //----- Flowering plant, use an ellipsoid. --------------------------------------//
         merge{
            //----- The stem. ------------------------------------------------------------//
            cylinder{ <0,0,0>, <0,1,0>, 1.0
                      texture{ pigment{ color rgb <0.64,0.32,0.16>}
                               finish { phong 1 }
                             } // end of texture
                      scale <rbh,crown,rbh>
                      translate<xtree,0,ztree>
                    } // end of sphere
            //----------------------------------------------------------------------------//



            //----- First crown. ---------------------------------------------------------//
            cone{ <0,0,0>, 1.0, <0,clength/3,0>, 1/3
                    texture{ pigment{ pftcol  }
                             finish { phong 1 }
                           } // end texture
                    scale <cwidth,1,cwidth>
                    translate<xtree,bole,ztree>
                  } // end cone
            //----- Second crown. --------------------------------------------------------//
            cone{ <0,0,0>, 0.8, <0,clength/3,0>, 1/3
                    texture{ pigment{ pftcol  }
                             finish { phong 1 }
                           } // end texture
                    scale <cwidth,1,cwidth>
                    translate<xtree,c13,ztree>
                  } // end cone
            //----- First crown. ---------------------------------------------------------//
            cone{ <0,0,0>, 0.6, <0,clength/3,0>, 0.0
                    texture{ pigment{ pftcol  }
                             finish { phong 1 }
                           } // end texture
                    scale <cwidth,1,cwidth>
                    translate<xtree,c23,ztree>
                  } // end cone
            //----------------------------------------------------------------------------//
         } // end merge
         //-------------------------------------------------------------------------------//
         #break
   #end // switch(pftshape)
   //-------------------------------------------------------------------------------------//
#end // macro plant
//----------------------------------------------------------------------------------------//


//----- Set the camera -------------------------------------------------------------------//
camera{ location  <-450.0 ,350.0 , -450.0>
        look_at   <   0.0 ,  0.0 ,    0.0>
        right x*image_width/image_height
        angle 45 }
//----------------------------------------------------------------------------------------//


//----- Fiat lux... ----------------------------------------------------------------------//
light_source{<1500,3000,-2500> color rgb<0.66,0.66,0.66>}
//----------------------------------------------------------------------------------------//


//----- Ground ---------------------------------------------------------------------------//
box{ <-205,-1,-205>,<205,0,205>
     texture{ pigment{ color rgb<0.96,0.76,0.60>}
              normal { bumps 0.75 scale 0.25 }
              finish { phong 0.1 }
            } // end texture
   } // end box
//----------------------------------------------------------------------------------------//
