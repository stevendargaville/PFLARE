# Just import the pflare definitions
import pflare_defs
# And register all the types with PETSc
pflare_defs.py_PCRegister_PFLARE()
# Define CF splitting variables
CF_PMISR_DDC = 0
CF_PMIS=1
CF_PMIS_DIST2=2
CF_AGG=3
CF_PMIS_AGG=4
