# ----------------------------
# PHYS1 - Physics 1 solver for the TI-84 Plus CE
#   make                      -> bin/PHYS1.8xp  (solver + Ask AI over USB)
#   make NAME=PHYSL AI=0      -> bin/PHYSL.8xp  (solver only, no USB)
# Run "make clean" between the two.
# ----------------------------

NAME ?= PHYS1
AI ?= 1
DESCRIPTION = "Physics 1 solver"
COMPRESSED = NO

CFLAGS = -Wall -Wextra -Oz -DUSE_CLAUDE=$(AI)
CXXFLAGS = -Wall -Wextra -Oz

# ----------------------------

include $(shell cedev-config --makefile)
