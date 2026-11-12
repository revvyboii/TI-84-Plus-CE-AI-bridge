# ----------------------------
# PHYS - Physics 1 solver for the TI-84 Plus CE
#   make NAME=PHYS1 PART=1            motion, forces, energy
#   make NAME=PHYS2 PART=2            momentum, collisions, changing forces
#   make NAME=PHYS3 PART=3            rotation, gravity, springs, fluids
#   add AI=0 for a build with no USB / Ask AI code (e.g. NAME=PHL1 PART=1 AI=0)
# Run "make clean" between builds.
# ----------------------------

NAME ?= PHYS1
PART ?= 1
AI ?= 1
DESCRIPTION = "Physics 1 solver"
COMPRESSED = YES

CFLAGS = -Wall -Wextra -Wno-unused-function -Wno-unused-const-variable -Oz -DUSE_CLAUDE=$(AI) -DPART=$(PART)
CXXFLAGS = -Wall -Wextra -Oz

# ----------------------------

include $(shell cedev-config --makefile)
