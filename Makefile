hydro_solver: hydro_solver.c
	gcc -o hydro_solver hydro_solver.c -I.

output.log: hydro_solver
	./hydro_solver > output.log

output.pdf: output.log
	python plot_output.py

# Basically an alias
run: output.log
	python plot_output.py

