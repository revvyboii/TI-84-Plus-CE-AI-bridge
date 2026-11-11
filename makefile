# ----------------------------
# CLAUDE - chat with Claude from a TI-84 Plus CE
# Build with the CE C toolchain:  make
# Output: bin/CLAUDE.8xp
# ----------------------------

NAME = CLAUDE
DESCRIPTION = "Chat with Claude over USB"
COMPRESSED = NO

CFLAGS = -Wall -Wextra -Oz
CXXFLAGS = -Wall -Wextra -Oz

# ----------------------------

include $(shell cedev-config --makefile)
