# ==============================================================================
# 1. PETSc Directory & Variables
# ==============================================================================
PETSC_DIR = /opt/homebrew/Cellar/petsc/3.24.2
EXTRA_LIB_PATH = -L/opt/homebrew/opt/openblas/lib -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
APP = fea_phd

# ==============================================================================
# 2. DEFAULT TARGET
# ==============================================================================
all: debug ${APP}

# ==============================================================================
# 3. Include PETSc Rules
# ==============================================================================
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# ==============================================================================
# 4. Source Files (EXPLICIT ORDER)
# ==============================================================================
SOURCES = mod_globals.F90 \
          mod_mesh.F90 \
          mod_io.F90 \
          mod_stiffness.F90 \
          mod_solver.F90 \
          mod_post.F90 \
          main.F90
		  
OBJECTS = $(SOURCES:.F90=.o)

# ==============================================================================
# 5. Debug Rules
# ==============================================================================
debug:
	@echo "----------------------------------------"
	@echo "   DEBUG MODE: Compiling in Order"
	@echo "----------------------------------------"

# ==============================================================================
# 6. Linking Rule
# ==============================================================================
${APP}: ${OBJECTS}
	@echo "--> Linking Application..."
	-${FLINKER} -o ${APP} ${OBJECTS} ${PETSC_LIB} ${EXTRA_LIB_PATH}
	@echo "--> SUCCESS! Executable created: ./${APP}"
	@echo "----------------------------------------"
	${RM} *.o *.mod

clean_app:
	@echo "--> Cleaning..."
	${RM} ${APP} *.o *.mod