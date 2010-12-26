from glob import glob
from pymol import cmd
import re

fn_re = re.compile(r'mov(\d+)\.pdb')
file_list = glob("mov*.pdb")
file_list.sort(key=lambda f:int(fn_re.match(f).groups()[0]))

for file in file_list:
   cmd.color("white")
   cmd.load(file,"mov")
   cmd.hide("everything")
   cmd.show("spheres")
   cmd.show("lines")
   cmd.color("white")
   cmd.recolor

cmd.mset("1 -%d"%len(file_list))
#cmd.mpng("mov")
