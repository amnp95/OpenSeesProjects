clear
conda deactivate

python3 modelCreator.py


mpirun -np 8 OpenSeesMP model.tcl

