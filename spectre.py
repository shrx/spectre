#!/usr/bin/python3
import numpy as np
from time import time
import matplotlib.pyplot as plt

import drawsvg

## configlation
#* increase this number for larger tilings.
N_ITERATIONS = 5
#* shape Edge_ration tile(Edge_a, Edge_b)
Edge_a = 20.0 / (np.sqrt(3) + 1.0)
Edge_b = 20.0 - Edge_a

# COLOR_MAP_orig
COLOR_MAP = {
	'Gamma': np.array((255, 255, 255),'f')/255.,
	'Gamma1': np.array((255, 255, 255),'f')/255.,
	'Gamma2': np.array((255, 255, 255),'f')/255.,
	'Delta': np.array((220, 220, 220),'f')/255.,
	'Theta': np.array((255, 191, 191),'f')/255.,
	'Lambda': np.array((255, 160, 122),'f')/255.,
	'Xi': np.array((255, 242, 0),'f')/255.,
	'Pi': np.array((135, 206, 250),'f')/255.,
	'Sigma': np.array((245, 245, 220),'f')/255.,
	'Phi': np.array((0, 255, 0),'f')/255.,
	'Psi': np.array((0, 255, 255),'f')/255.
}

# COLOR_MAP_mystics 
COLOR_MAP = {
	'Gamma': np.array((196, 201, 169),'f')/255.,
	'Gamma1': np.array((196, 201, 169),'f')/255.,
	'Gamma2': np.array((156, 160, 116),'f')/255.,
	'Delta': np.array((247, 252, 248),'f')/255.,
	'Theta': np.array((247, 252, 248),'f')/255.,
	'Lambda': np.array((247, 252, 248),'f')/255.,
	'Xi': np.array((247, 252, 248),'f')/255.,
	'Pi': np.array((247, 252, 248),'f')/255.,
	'Sigma': np.array((247, 252, 248),'f')/255.,
	'Phi': np.array((247, 252, 248),'f')/255.,
	'Psi': np.array((247, 252, 248),'f')/255.
}

# COLOR_MAP_pride
COLOR_MAP = {
    "Gamma":  np.array((255, 255, 255),'f')/255.,
    "Gamma1": np.array(( 97,  57,  21),'f')/255.,
    "Gamma2": np.array(( 64,  64,  64),'f')/255.,
    "Delta":  np.array((  2, 129,  33),'f')/255.,
    "Theta":  np.array((  0,  76, 255),'f')/255.,
    "Lambda": np.array((118,   0, 136),'f')/255.,
    "Xi":     np.array((229,   0,   0),'f')/255.,
    "Pi":     np.array((255, 175, 199),'f')/255.,
    "Sigma":  np.array((115, 215, 238),'f')/255.,
    "Phi":    np.array((255, 141,   0),'f')/255.,
    "Psi":    np.array((255, 238,   0),'f')/255.
}
## end of configilation.

TILE_NAMES = ["Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Phi", "Psi"]

PI = np.pi

IDENTITY = np.array([[1,0,0],
                     [0,1,0]], 'float32')

def get_spectre_points(edge_a, edge_b):
    a = edge_a
    a_sqrt3_d2 = a * np.sqrt(3)/2 # a*sin(60 deg)
    a_d2 = a * 0.5  # a* cos(60 deg)

    b = edge_b
    b_sqrt3_d2 = b * np.sqrt(3) / 2 # b*sin(60 deg)
    b_d2 = b * 0.5 # b* cos(60 deg)

    spectre_points = np.array([
		(0                        , 0                            ), #// 1: - b
		(a                        , 0                            ), #// 2: + a
		(a +     a_d2             , 0 - a_sqrt3_d2               ), #// 3: + ~a
		(a +     a_d2 + b_sqrt3_d2, 0 - a_sqrt3_d2 +         b_d2), #// 4: + ~b
		(a +     a_d2 + b_sqrt3_d2, 0 - a_sqrt3_d2 +     b + b_d2), #// 5: + b
		(a + a + a_d2 + b_sqrt3_d2, 0 - a_sqrt3_d2 +     b + b_d2), #// 6: + a
		(a + a + a +    b_sqrt3_d2,                      b + b_d2), #// 7: + ~a
		(a + a + a                ,                  b + b       ), #// 8: - ~b 
		(a + a + a    - b_sqrt3_d2,                  b + b - b_d2), #// 9: - ~b
		(a + a + a_d2 - b_sqrt3_d2,     a_sqrt3_d2 + b + b - b_d2), #// 10: +~a
		(a +     a_d2 - b_sqrt3_d2,     a_sqrt3_d2 + b + b - b_d2), #// 11: -a
		(        a_d2 - b_sqrt3_d2,     a_sqrt3_d2 + b + b - b_d2), #// 12: -a
		(0            - b_sqrt3_d2,                  b + b - b_d2), #// 13: -~a
		(0                        ,                      b       )  #// 14: +~b
    ], 'float32')
    # print(spectre_points)
    return spectre_points
   
SPECTRE_POINTS = get_spectre_points(Edge_a, Edge_b) # tile(Edge_a, Edge_b)
Mystic_SPECTRE_POINTS = get_spectre_points(Edge_b, Edge_a) # tile(Edge_b, Edge_a)
SPECTRE_QUAD = SPECTRE_POINTS[[3,5,7,11],:]

def flattenPts(lst): # drowsvg
    return [item for sublist in lst for item in sublist] # drowsvg

SPECTRE_SHAPE = drawsvg.Lines(*flattenPts([p for p in SPECTRE_POINTS]), stroke="black", stroke_width=0.5,close=True) # drowsvg
Mystic_SPECTRE_SHAPE = drawsvg.Lines(*flattenPts([p for p in Mystic_SPECTRE_POINTS]), stroke="black",   stroke_width=0.5, close=True) # drowsvg
num_tiles = 0 # drowswvg
def drawPolygon(drawing, T, label, f, s, w): #drowsvg
    """
    drawing: drawing to draw on
    T: transformation matrix
    label: label of shape type
    f: tile fill color
    s: tile stroke color
    w: tile stroke width
    """
    global num_tiles
    num_tiles += 1

    fill = f"rgb({int(round(f[0]* 255, 0))}, {int(round(f[1]* 255,0))}, {int(round(f[2]* 255,0))})"
    stroke_f = s
    stroke_w = w if s else 0
    shape = Mystic_SPECTRE_SHAPE if label == "Gamma2" else SPECTRE_SHAPE # geometric points used.

    drawing.append(drawsvg.Use(
        shape,
        0, 0,
        transform=f"matrix({T[0,0]} {T[1,0]} {T[0,1]} {T[1,1]} {T[0,2]} {T[1,2]})",
        fill=fill,
        stroke=stroke_f,
        stroke_width=stroke_w))

# Matrix * point
def mul(A, B):
    AB = A.copy()
    AB[:,:2] = A[:,:2].dot(B[:,:2]) 
    AB[:,2] += A[:,:2].dot(B[:,2])
    return AB

class Tile:
    def __init__(self, label):
        """
        _: NO list of Tile coordinate points
        label: Tile type used for shapes coloring
        """
        self.label = label
        self.quad = SPECTRE_QUAD.copy()

    def draw(self, polygons, tile_transformation=IDENTITY.copy()):
        vertices = (SPECTRE_POINTS if self.label != "Gamma2" else Mystic_SPECTRE_POINTS).dot(tile_transformation[:,:2].T) + tile_transformation[:,2]
        polygons.append((vertices, self.label))

    def drawPolygon(self, drawing, tile_transformation=IDENTITY):
        return drawPolygon(drawing, tile_transformation, self.label, COLOR_MAP[self.label], "black", 0.1)

class MetaTile:
    def __init__(self, tiles=[], transformations=[], quad=SPECTRE_QUAD.copy()):
        """
        tiles: list of Tiles(No points)
        transformations: list of transformation matrices
        quad: MetaTile quad points
        """
        self.tiles = tiles
        self.transformations = transformations
        self.quad = quad

    def draw(self, polygons, transformation=IDENTITY.copy()):
        """
        recursively expand MetaTiles down to Tiles and draw those
        """
        for tile, trsf in zip(self.tiles, self.transformations):
           tile.draw(polygons, mul(transformation, trsf))

    def drawPolygon(self, drawing, metatile_transformation=IDENTITY.copy()):
        """
        recursively expand MetaTiles down to Tiles and draw those
        """
        # TODO: parallelize?
        for tile, trsf in zip(self.tiles, self.transformations):
           tile.drawPolygon(drawing,  mul(metatile_transformation, trsf))
                            
def buildSpectreBase():
    ttrans = np.array([[1,0,SPECTRE_POINTS[8,0]],
                       [0,1,SPECTRE_POINTS[8,1]]])
    trot = np.array([[np.cos(PI/6),-np.sin(PI/6),0.],
                     [np.sin(PI/6), np.cos(PI/6),0.]],'float32')
    trsf = mul(ttrans, trot)
    tiles = {label: (Tile(label) ) for label in TILE_NAMES if label != "Gamma"}
    # special rule for Mystic == Gamma == Gamma1 + Gamma2
    tiles["Gamma"] = MetaTile(tiles=[Tile("Gamma1"),Tile("Gamma2")],
                                     transformations=[IDENTITY.copy(),trsf],
                                     quad=SPECTRE_QUAD.copy())
    return tiles

def buildSupertiles(input_tiles):
    """
    iteratively build on current system of tiles
    input_tiles = current system of tiles, initially built with buildSpectreBase()
    """
    # First, use any of the nine-unit tiles in "tiles" to obtain a
    # list of transformation matrices for placing tiles within supertiles.
    quad = input_tiles["Delta"].quad

    transformations = [IDENTITY.copy()]
    total_angle = 0
    trot = IDENTITY.copy()
    transformed_quad = quad
    for _angle, _from, _to in ((   PI/3, 3, 1),
                               (     0., 2, 0),
                               (   PI/3, 3, 1),
                               (   PI/3, 3, 1),
                               (     0., 2, 0),
                               (   PI/3, 3, 1),
                               (-2*PI/3, 3, 3)):
        if _angle != 0:
            total_angle += _angle
            trot = np.array([[1, 0,0],[0,1,0]])*np.cos(total_angle) \
                  +np.array([[0,-1,0],[1,0,0]])*np.sin(total_angle)
            transformed_quad = quad.dot(trot[:,:2].T) # + trot[:,2]
        last_trsf = transformations[-1]
        ttrans = IDENTITY.copy()
        ttrans[:,2] = last_trsf[:,:2].dot(quad[_from,:]) + last_trsf[:,2] \
                     -transformed_quad[_to,:]
        transformations.append(mul(ttrans, trot))

    R = np.array([[-1,0,0],[ 0,1,0]], 'float32')
    transformations = [ mul(R, trsf) for trsf in transformations ]

    # Now build the actual supertiles, labeling appropriately.
    super_quad = quad[[2,1,2,1],:]
    for i,itrsf in enumerate([6,5,3,0]):
        trsf = transformations[itrsf]
        super_quad[i,:] = trsf[:,:2].dot(super_quad[i,:]) + trsf[:,2]

    tiles = {}
    for label, substitutions in (("Gamma",  ("Pi",  "Delta", None,  "Theta", "Sigma", "Xi",  "Phi",    "Gamma")),
                                 ("Delta",  ("Xi",  "Delta", "Xi",  "Phi",   "Sigma", "Pi",  "Phi",    "Gamma")),
                                 ("Theta",  ("Psi", "Delta", "Pi",  "Phi",   "Sigma", "Pi",  "Phi",    "Gamma")),
                                 ("Lambda", ("Psi", "Delta", "Xi",  "Phi",   "Sigma", "Pi",  "Phi",    "Gamma")),
                                 ("Xi",     ("Psi", "Delta", "Pi",  "Phi",   "Sigma", "Psi", "Phi",    "Gamma")),
                                 ("Pi",     ("Psi", "Delta", "Xi",  "Phi",   "Sigma", "Psi", "Phi",    "Gamma")),
                                 ("Sigma",  ("Xi",  "Delta", "Xi",  "Phi",   "Sigma", "Pi",  "Lambda", "Gamma")),
                                 ("Phi",    ("Psi", "Delta", "Psi", "Phi",   "Sigma", "Pi",  "Phi",    "Gamma")),
                                 ("Psi",    ("Psi", "Delta", "Psi", "Phi",   "Sigma", "Psi", "Phi",    "Gamma"))):
        tiles[label] = MetaTile(tiles=[input_tiles[subst] for subst in substitutions if subst],
                     transformations=[trsf for subst, trsf in zip(substitutions, transformations) if subst],
                     quad=super_quad)
    return tiles

start = time()
tiles = buildSpectreBase()
for _ in range(N_ITERATIONS):
    tiles = buildSupertiles(tiles)
time1 = time()-start
print(f"supertiling loop took {round(time1, 4)} seconds")

start = time()
polygons = []
tiles["Delta"].draw(polygons)
time2 = time()-start
print(f"matplotlib.pyplot: tile recursion loop took {round(time2, 4)} seconds, generated {len(polygons)} tiles")

start = time()
plt.figure(figsize=(8, 8))
plt.axis('equal')
for pts,label in polygons:
    # plt.text((pts[1,0] + pts[7,0])/2, (pts[1,1] + pts[7,1])/2, label, fontsize=8, color='gray')
    plt.fill(pts[:,0],pts[:,1],facecolor=COLOR_MAP[label])
    plt.plot(pts[:,0],pts[:,1],color='gray',linewidth=0.2)

saveFileName = f"spectre_tile{Edge_a:.1f}-{Edge_b:.1f}_{N_ITERATIONS}-{len(polygons)}pts.svg"
print("matplotlib.pyplot: file save to " + saveFileName)
plt.savefig(saveFileName)
time3 = time()-start
print(f"matplotlib.pyplot SVG drawing took {round(time3, 4)} seconds")
print(f"matplotlib.pyplot total processing time {round(time1+time2+time3, 4)} seconds, {round(1000000*(time1+time2+time3)/len(polygons), 4)} μs/tile")

start = time()
d = drawsvg.Drawing(8000, 8000, origin="center") # @TODO: ajust to polygons X-Y min and max. 
tiles["Delta"].drawPolygon(d) # updates num_tiles
saveFileName = f"spectre_tile{Edge_a:.1f}-{Edge_b:.1f}_{N_ITERATIONS}-{num_tiles}useRef.svg"
d.save_svg(saveFileName)
time4 = time()-start
print(f"drowsvg: SVG drawing took {round(time4, 4)} seconds, generated {num_tiles} tiles")
print("drowsvg: drawPolygon save to " + saveFileName)
print(f"drowsvg: total processing time {round(time1+time4, 4)} seconds, {round(1000000*(time1+time4)/len(polygons), 4)} μs/tile")

plt.show()
