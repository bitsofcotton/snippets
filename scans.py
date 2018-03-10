#! /usr/bin/env python
# coding: UTF-8

import os.path
from gimpfu import *

def roundup_scan(timg, tdrawable, resolution = 72, ratio = 4, hilight = 207, radius = 8, strength = 1, thresho = 16, savecopy = TRUE):
	pdb.gimp_levels(tdrawable, 0, 0, hilight, 1, 0, 255)
	drawable = pdb.gimp_image_get_active_layer(timg)
	pdb.gimp_image_scale(timg,
						 float(pdb.gimp_image_width(timg))  / ratio,
						 float(pdb.gimp_image_height(timg)) / ratio)
	pdb.plug_in_unsharp_mask(timg, drawable, radius, strength, thresho)
	pdb.gimp_image_set_resolution(timg, resolution, resolution)
	
	if savecopy:
		root, ext = os.path.splitext(timg.filename)
		pdb.file_jpeg_save(timg, tdrawable,
						   root+".jpeg",
						   root+".jpeg",
						   0.8, 0, 0, 0, "", 0, 0, 0, 0)

register(
    "python-fu-roundup-scan",
    "Round up after scanning.",
    "Round up after scanning.",
    "kazunobu",
    "Public Domain",
    "2011",
   	"<Image>/Image/After scan",
    "GRAY*",
    [
		(PF_INT, "resolution", "resulution?", 72),
		(PF_FLOAT, "ratio", "1/n, n?", 4),
		(PF_INT, "hilight", "high light?", 207),
		(PF_INT, "radius", "unsharp radius?", 8),
		(PF_INT, "strength", "  strength? (0-5)", 1),
		(PF_INT, "thresho", "  threshold? (0-255)", 16),
		(PF_BOOL, "savecopy", "Save into JPEG?", TRUE),
	],
    [],
    roundup_scan)

main()
