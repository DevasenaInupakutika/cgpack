#$Id: debug.py 9 2014-12-01 09:55:21Z mexas $

from paraview.simple import *

reader=ImageReader(FilePrefix= "z9end.raw" )
reader.DataByteOrder=1
reader.DataExtent=[1,80,1,80,1,640]
reader.DataScalarType=5

view = GetRenderView()
Show()
