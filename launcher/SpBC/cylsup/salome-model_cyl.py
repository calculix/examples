# -*- coding: utf-8 -*-

###
### This file is parametric model.
### Let's define variables

R1=2 #Outer radius
R2=1.25 #Hole
L=24 #Distance between the holes
t=1 #thickness

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy



###
### GEOM component
###

#Let's correct numbers using variables

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Cylinder_1 = geompy.MakeCylinderRH(R1, t)

Translation_1 = geompy.MakeTranslation(Cylinder_1, 0, L, 0)
Box_1 = geompy.MakeBoxDXDYDZ(2*R1, L, t)
Translation_2 = geompy.MakeTranslation(Box_1, -R1, 0, 0)
Fuse_1 = geompy.MakeFuseList([Cylinder_1, Translation_1, Translation_2], True, True)
Cylinder_2 = geompy.MakeCylinderRH(R2, t)
Cut_1 = geompy.MakeCutList(Fuse_1, [Cylinder_2], True)
geompy.TranslateDXDYDZ(Cylinder_2, 0, L, 0)
Cut_2 = geompy.MakeCutList(Cut_1, [Cylinder_2], True)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudy( Translation_2, 'Translation_2' )
geompy.addToStudy( Fuse_1, 'Fuse_1' )
geompy.addToStudy( Cylinder_2, 'Cylinder_2' )
geompy.addToStudy( Cut_1, 'Cut_1' )
geompy.addToStudy( Cut_2, 'Cut_2' )


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)
