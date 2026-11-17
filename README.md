# TI-84-Plus-CE-AI-bridge

Programs for the TI-84 Plus CE. To put one on the calculator: open its `.8xp` file below,
click the download button, then send it with WebTILP (**Choose File** → **Send Selected Files**).

| File | What it is |
| --- | --- |
| `phys1/PHYS1.8xp` | Physics 1 solver: motion, forces, energy (+ Ask AI) |
| `phys1/PHYS2.8xp` | Physics 1 solver: momentum, collisions, center of mass, Help me pick (+ Ask AI) |
| `phys1/PHYS3.8xp` | Physics 1 solver: rotation, gravity, springs, fluids (+ Ask AI) |
| `phys1/PHL1.8xp` – `PHL3.8xp` | The same three solvers, offline (no USB) |
| `chem/CHEM.8xp` | Chemistry 1 solver with a periodic table, offline ([details](chem/README.md)) |

All of them need the CE C libraries (clibs). Every push also builds fresh copies of all of these:
the **Actions** tab → latest run → *send-these-to-your-calculator*.

Source: `calc/` (CLAUDE chat program), `phys1/src/main.c`, `chem/src/main.c`.
