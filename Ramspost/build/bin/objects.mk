#Makefile objects.mk

# Define main source.

MAIN = $(RPOST_DRIVER)/rpost_main.f90
MAINOBJ = rpost_main.o


# Define objects.

OBJECTS =                           \
            abort_run.o             \
            an_header.o             \
            brams_data.o            \
            comp_lib.o              \
            charutils.o             \
            dateutils.o             \
            dted.o                  \
            eenviron.o              \
            error_mess.o            \
            interp_lib.o            \
            leaf_coms.o             \
            micro_coms.o            \
            misc_coms.o             \
            numutils.o              \
            parlib.o                \
            polarst.o               \
            rams_read_header.o      \
            rfvec.o                 \
            rcio.o                  \
            rconstants.o            \
            rgrad.o                 \
            rnamel.o                \
            rnumr.o                 \
            rout_coms.o             \
            rpost_coms.o            \
            rpost_dims.o            \
            rpost_filelist.o        \
            rpost_misc.o            \
            rsys.o                  \
            rutil.o                 \
            scratch_coms.o          \
            soil_coms.o             \
            somevars.o              \
            therm_lib.o             \
            tmpname.o               \
            utils_c.o               \
            utils_f.o               \
            variables.o             \
            vformat_brams3.3.o
