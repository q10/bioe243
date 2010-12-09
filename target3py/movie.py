import string
import random
from pymol.cgo import *
from pymol import cmd
from target3mod import *

cluster = SCCluster()
box = 10.0

axes = [
       LINEWIDTH, 2.0,
       BEGIN, LINES,
       COLOR, 0.2, 1.0, 0.2,

       VERTEX, 0.0, 0.0, 0.0,
       VERTEX, box, 0.0, 0.0,
      
       #VERTEX, box.0, 0.0, 0.0,
       #VERTEX, box, box, 0.0,

       #VERTEX, box, box, 0.0,
       #VERTEX, 0.0, 100.0, 0.0,

       #VERTEX, box, box, 0.0,
       #VERTEX, box, box, box,

       #VERTEX, 0.0, 100.0, 0.0,
       #VERTEX, 0.0, 100.0, 100.0,

       #VERTEX, box, 0.0, 0.0,
       #VERTEX, box, 0.0, box,

       #VERTEX, 0.0, 0.0, box,
       #VERTEX, 0.0, box, box,

       #VERTEX, 0.0, box, box,
       #VERTEX, box, box, box,

       #VERTEX, box, box, box,
       #VERTEX, box, 0.0, box,

       #VERTEX, box, 0.0, box,
       #VERTEX, 0.0, 0.0, box,

       COLOR, 1.0, 0.2, 0.2,
       VERTEX, 0.0, 0.0, 0.0,
       VERTEX, 0.0, box, 0.0,

       COLOR, 0.2, 0.2, 1.0,
       VERTEX, 0.0, 0.0, 0.0,
       VERTEX, 00, 0.0, box,
       END
       ]



c=0
for a in xrange(0,100):
   cluster.MD_step()
   balls = []
   for i in xrange(0,cluster.N):
       balls.extend([COLOR,
                     1.0,
                     1.0,
                     1.0,
                     SPHERE,
                     cluster.particles[i].x,
                     cluster.particles[i].y,
                     cluster.particles[i].z,
                     0.5,])
   obj = axes + balls
   cmd.load_cgo(obj,'cgo01',c)
   c = c + 1

pdb_list = [
"HETATM 1 X UNK 1 %8.3f%8.3f%8.3f 1.00 10.00\n"%(3.2,0,0), "HETATM 2 Y UNK 2 %8.3f%8.3f%8.3f 1.00 10.00\n"%(0,3.2,0), "HETATM 3 Z UNK 3 %8.3f%8.3f%8.3f 1.00 10.00\n"%(0,0,3.2),
           ]
cmd.read_pdbstr(string.join(pdb_list,''),'lab2')
cmd.hide('(lab2)')
cmd.label('lab2','name')
cmd.color('white','lab2')
#cmd.zoom('cgo01')
cmd.clip('far',-5)
