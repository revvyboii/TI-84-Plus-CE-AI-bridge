# CHEM: Chemistry 1 solver for the TI-84 Plus CE

Offline (no USB, no AI). Same screens and keys as the PHYS solvers.

## Put it on the calculator

1. Click **CHEM.8xp** in this folder, then the download button (Download raw file).
2. In WebTILP: **Choose File** → CHEM.8xp → **Send Selected Files**.
3. Run it from **PRGM → CHEM**.

It needs the CE C libraries (clibs), the same ones PHYS already uses, and about 63 KB of free RAM.

## What it solves (23 problem types)

| Category | Topics |
| --- | --- |
| Measurement | Density, Percent error, Temperature units |
| Atoms + moles | Grams/moles/particles, Percent by mass, Average atomic mass, Empirical formula |
| Solutions | Molarity, Dilution |
| Gases | Ideal gas PV = nRT, Combined gas law, Gas at STP, Partial pressures |
| Reactions | Stoichiometry, Limiting reactant, Percent yield |
| Heat | q = m c dT, Calorimetry, Heat of reaction |
| Light + electrons | Light and photons, Hydrogen atom (Bohr) |
| Acids + bases | pH and pOH, Titration |

## Keys

- **ZOOM** on the home screen: *Help me pick* asks a few questions, opens the right topic and walks you box by box.
- **GRAPH** in a topic (or **WINDOW** on the home screen) opens the periodic table:
  ALPHA + a letter jumps (ALPHA C = carbon, again = Cl, Ca, ...), ENTER adds one atom,
  digits set how many, DEL removes, TRACE shows % composition, GRAPH puts the molar mass in the box.
- **LEFT** makes a value negative, **ENTER** on a solved value shows the steps, **ZOOM** in a topic shows tips.

## Building

`make -C chem` with the CE C toolchain.
PC tests: `gcc -DHOST_TEST -o chem_test chem/src/main.c -lm && ./chem_test`
