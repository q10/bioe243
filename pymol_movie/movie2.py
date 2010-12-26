from glob import glob
from pymol import cmd

file_list = glob("mov*.pdb")

for file in file_list:
   cmd.color("white")
   cmd.load(file,"mov")
   cmd.hide("everything")
   cmd.show("nb_spheres")
   cmd.color("white")
   cmd.recolor

cmd.mset("1 -%d"%len(file_list))
