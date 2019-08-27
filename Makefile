hydro_solver: hydro_solver.c
	gcc -o hydro_solver hydro_solver.c -I.

hydro_solver_second_order: hydro_solver_second_order.c
	gcc -o hydro_solver_second_order hydro_solver_second_order.c -I.

output.log: hydro_solver
	./hydro_solver > output.log

output_second_order.log: hydro_solver_second_order
	./hydro_solver_second_order > output_second_order.log

output.pdf: output.log
	python plot_output.py output.log

output_second_order.pdf: output_second_order.log
	python plot_output.py output_second_order.log; mv output.pdf output_second_order.pdf

# Basically an alias
first: output.log
	python plot_output.py output.log

second: output_second_order.log
	python plot_output.py output_second_order.log; mv output.pdf output_second_order.pdf


