#Makefile objects.mk

# Define main source.

MAIN = $(RPOST_DRIVER)/rpost_main.f90
MAINOBJ = rpost_main.o


# Define objects.

OBJECTS =                           \
            abort_run.o             \
            an_header.o             \
            charutils.o             \
            dateutils.o             \
            dted.o                  \
            eenviron.o              \
            error_mess.o            \
            interp_lib.o            \
            parlib.o                \
            polarst.o               \
            rams_read_header.o      \
            rfvec.o                 \
            rcio.o                  \
            rconstants.o            \
            rgrad.o                 \
            rnamel.o                \
            rnumr.o                 \
            rpost_misc.o            \
            rutil.o                 \
            somevars.o              \
            therm_lib.o             \
            tmpname.o               \
            utils_c.o               \
            utils_f.o               \
            variables.o             \
            vformat_brams3.3.o
