# ../example

# <FILE> = ac.cc
# loesen der Allen Cahn Gleichung in 3D mit Neumann-Randdaten auf uniformen
# parallel verteilter Triangulierung mittels parallelem PCG Verfahren

# <FILE> = ac_adaptive.cc
# loesen der Allen Cahn Gleichung in 3D mit Neumann-Randdaten auf adaptiv verfeinerter und
# parallel verteilter Triangulierung mittels parallelem PCG Verfahren

# <FILE> = poisson.cc
# loesen der Poisson-Gleichung in 2D auf uniformen parallel verteilter Triangulierung mittels 
# parallelem PCG Verfahren


# compilieren mit Optimierungs-Flags
mpicxx <FILE> -O3 -funroll-loops -DNDEBUG -std=c++11 -lparmetis -lmetis -o <OUTPUT>

# ausfuehren mit 4 Prozessen
mpirun -np 4 <OUTPUT>

# compilieren im Debugg-Modus
mpicxx <FILE> -ggdb -std=c++11 -DDEBUG -lparmetis -lmetis -o <OUTPUT>

# ausfuehren im Debug-Modus mit 4 Prozessen
mpirun -np 4 xterm -e gdb <OUTPUT> -ex run
