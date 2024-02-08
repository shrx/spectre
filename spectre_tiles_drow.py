## draw Polygons Svg by drawsvg #####
from spectre import buildSpectreTiles,get_color_array, SPECTRE_POINTS, Mystic_SPECTRE_POINTS, Edge_a,Edge_b, N_ITERATIONS, print_trot_inv_prof, trot_inv
from spectre import get_transformation_min, get_transformation_max
from time import time
import drawsvg

start = time()
spectreTiles = buildSpectreTiles(N_ITERATIONS,Edge_a,Edge_b)
time1 = time()-start
print(f"supertiling loop took {round(time1, 4)} seconds")
transformation_min = get_transformation_min()
transformation_max = get_transformation_max()
print(f"transformation range is {transformation_min} to {transformation_max}")

start = time()
def flattenPts(lst): # drowsvg
    return [item for sublist in lst for item in sublist] # drowsvg

SPECTRE_SHAPE = drawsvg.Lines(*flattenPts([p for p in SPECTRE_POINTS]), stroke="black", stroke_width=0.5,close=True) # drowsvg
Mystic_SPECTRE_SHAPE = drawsvg.Lines(*flattenPts([p for p in Mystic_SPECTRE_POINTS]), stroke="black",   stroke_width=0.5, close=True) # drowsvg

viewWidth = transformation_max - transformation_min
svgContens = drawsvg.Drawing(viewWidth,
                    viewWidth,
                    ) # @TODO: ajust to polygons X-Y min and max. 
svgContens.view_box = (transformation_min,transformation_min,viewWidth,viewWidth)
SvgContens_drowSvg_transform_scaleY = 0 # spectreTiles["Delta"].transformations[0][0,0]
num_tiles = 0 # drowswvg
def drawPolygon2Svg(T, label): #drowsvg
    """
    T: transformation matrix
    label: label of shape type
    """
    global num_tiles,svgContens,SvgContens_drowSvg_transform_scaleY
    num_tiles += 1
    color_array = get_color_array(T, label) # drowsvg
    degAngle, _scaleY = trot_inv(T) 
    # @TODO: It depends on the initial data being one of the following
    #         [[-1. 0. 0.][ 0. 1. 0.]] [[1. 0. 0.] [0. 1. 0.]] first 
    # @TODO: However, it is not good to depend on the order of data . repare trot_inv(T) .
    if SvgContens_drowSvg_transform_scaleY == 0:
        SvgContens_drowSvg_transform_scaleY = 1 if T[0,0] > 0 else (-1 if T[0,0] < 0 else 0)
    
    fill = f"rgb({int(round(color_array[0]* 255, 0))}, {int(round(color_array[1]* 255,0))}, {int(round(color_array[2]* 255,0))})"
    stroke_f = "gray" # tile stroke color
    stroke_w = 0.1 if (fill[0] != 0) | (fill[1] != 0) | (fill[2] != 0) else 0 # tile stroke width
    shape = SPECTRE_SHAPE if label != "Gamma2" else Mystic_SPECTRE_SHAPE  # geometric points used.
    # print(f"transform-matrix,{T[0,0]},{T[1,0]},{T[0,1]},{T[1,1]},{T[0,2]},{T[1,2]}")

    svgContens.append(drawsvg.Use(
        shape,
        0, 0,
        transform=f"translate({T[0,2]},{T[1,2]}) rotate({degAngle}) scale(1,{SvgContens_drowSvg_transform_scaleY})",
        # transform=f"matrix({T[0,0]} {T[1,0]} {T[0,1]} {T[1,1]} {T[0,2]} {T[1,2]})",
        fill=fill,
        stroke=stroke_f,
        stroke_width=stroke_w))

spectreTiles["Delta"].drawPolygon(drawPolygon2Svg) # updates num_tiles
saveFileName = f"spectre_tile{Edge_a:.1f}-{Edge_b:.1f}_{N_ITERATIONS}-{num_tiles}useRef.svg"
svgContens.save_svg(saveFileName)
time4 = time()-start
print(f"drowsvg: SVG drawing took {round(time4, 4)} seconds, generated {num_tiles} tiles")
print_trot_inv_prof()
print("drowsvg: drawPolygon save to " + saveFileName)
print(f"drowsvg: total processing time {round(time1+time4, 4)} seconds, {round(1000000*(time1+time4)/num_tiles, 4)} μs/tile")
