#!/usr/bin/python3
import numpy as np
## configlation
#* increase this number for larger tilings.
N_ITERATIONS = 3
#* shape Edge_ration tile(Edge_a, Edge_b)
Edge_a = 20.0 / (np.sqrt(3) + 1.0)
Edge_b = 20.0 - Edge_a
## end of configilation.

TILE_NAMES = ["Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Phi", "Psi"]

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

IDENTITY = np.array([[1,0,0],[0,1,0]], 'float32') # == trot(0)

# Rotation matrix for Affine transform
trot_memo = {
     0:  np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0]]),
     30: np.array([[np.sqrt(3)/2, -0.5, 0.0], [0.5, np.sqrt(3)/2, 0.0]]),
     60: np.array([[0.5, -np.sqrt(3)/2, 0.0], [np.sqrt(3)/2, 0.5, 0.0]]),
     120: np.array([[-0.5, -np.sqrt(3)/2, 0.0], [np.sqrt(3)/2, -0.5, 0.0]]),
     180: np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
     240: np.array([[-0.5, np.sqrt(3)/2, 0.0], [-np.sqrt(3)/2, -0.5, 0.0]]),
}
def trot(degAngle):
    """
    degAngle: integer degree angle 
    """
    global trot_memo
    if degAngle not in trot_memo:
        ang = np.deg2rad(degAngle)
        c = np.cos(ang)
        s = np.sin(ang)
        trot_memo[degAngle] = np.array([[c, -s, 0],[s, c, 0]])
        print(f"trot_memo[{degAngle}]={trot_memo[degAngle]}")
    return trot_memo[degAngle].copy()

def trot_inv(T):
    """
    T: rotation matrix for Affine transform
    """
    degAngle1 = int(np.round(np.rad2deg(np.arctan2(T[1, 0], T[0, 0]))))
    if degAngle1 == -180:
        degAngle1 = 180
    degAngle2 = int(np.round(np.rad2deg(np.arctan2(-T[0, 1], T[1, 1]))))
    if (degAngle1 == degAngle2): # self validate angle
        scaleY = 1
    elif (degAngle1 == (-degAngle2)):
        scaleY = 1
    elif (degAngle1 == (180 - degAngle2)) or (degAngle2 == (180 - degAngle1)):
        scaleY = -1
    elif (degAngle1 == (degAngle2 - 180)) or (degAngle2 == (degAngle1 - 180)):
        scaleY = -1
    else:
        scaleY = -1
        print(f"ValueError at trot_inv: degAngle1={degAngle1}, degAngle2={degAngle2} T={T}")
        # raise ValueError("trot_inv: degAngle1.abs != degAngle2.abs")
        
    return (degAngle1, scaleY)

# Matrix * point
def transPt(trsf, quad):
    trPt = (trsf[:,:2].dot(quad) + trsf[:,2])
    # print(f"at transPt={trPt}")
    return trPt

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
        self.quad = SPECTRE_QUAD

    def forEachTile(self, doProc, tile_transformation=IDENTITY):
        # print(f"at Tile.drawPolygon {self.label} angle={trot_inv(tile_transformation)} tile_transformation={tile_transformation}")
        return doProc(tile_transformation, self.label)

class MetaTile:
    def __init__(self, tiles=[], transformations=[], quad=SPECTRE_QUAD):
        """
        tiles: list of Tiles(No points)
        transformations: list of transformation matrices
        quad: MetaTile quad points
        """
        self.tiles = tiles
        self.transformations = transformations
        self.quad = quad

    def forEachTile(self, doProc, transformation=IDENTITY):
        """
        recursively expand MetaTiles down to Tiles and draw those
        """
        # TODO: parallelize?
        for tile, trsf in zip(self.tiles, self.transformations):
           tile.forEachTile(doProc, (mul(transformation, trsf)))
                            
def buildSpectreBase():
    tiles = {label: (Tile(label) ) for label in TILE_NAMES if label != "Gamma"}
    # special rule for Mystic == Gamma == Gamma1 + Gamma2
    tiles["Gamma"] = MetaTile(tiles=[Tile("Gamma1"),
                                     Tile("Gamma2")
                              ],
                              transformations=[
                                         IDENTITY.copy(),
                                         mul(np.array([
                                             [1,0,SPECTRE_POINTS[8,0]],
                                             [0,1,SPECTRE_POINTS[8,1]]
                                         ]), trot(30))
                              ],
                              quad=SPECTRE_QUAD.copy())
    # print(f"at buildSpectreBase: tiles[Gamma]={tiles['Gamma'].transformations}")
    return tiles

def get_transformation_range():
    global transformation_min_X,transformation_min_Y,transformation_max_X,transformation_max_Y
    return (transformation_min_X,transformation_min_Y,transformation_max_X,transformation_max_Y)
    
def buildSupertiles(input_tiles):
    """
    iteratively build on current system of tiles
    input_tiles = current system of tiles, initially built with buildSpectreBase()
    """
    # First, use any of the nine-unit tiles in "tiles" to obtain a
    # list of transformation matrices for placing tiles within supertiles.
    quad = input_tiles["Delta"].quad

    total_angle = 0
    rotation =  trot(total_angle) # IDENTITY.copy() #
    transformations =  [rotation.copy()] # [IDENTITY.copy()]
    transformed_quad = quad
    for _angle, _from, _to in ((  60, 3, 1),
                               (   0, 2, 0),
                               (  60, 3, 1),
                               (  60, 3, 1),
                               (   0, 2, 0),
                               (  60, 3, 1),
                               (-120, 3, 3)):
        if _angle != 0:
            total_angle += _angle
            rotation = trot(total_angle)
            transformed_quad = np.array([transPt(rotation, quad1) for quad1 in quad]) ### quad.dot(rotation[:,:2].T) # + trot[:,2]
        ttrans = IDENTITY.copy()
        ttrans[:,2] = transPt(transformations[-1], quad[_from]) - transformed_quad[_to,:]
        transformations.append(mul(ttrans, rotation))

    R = np.array([[-1.0, 0.0, 0.0],[0.0, 1.0, 0.0]]) # @TODO: Not trot(180).  Instead of rotating 180 degrees, get a mirror image.
    transformations = [(mul(R, trsf)) for trsf in transformations ] # @TODO Note that mul(trsf, R) is not commutible
    # @TODO: TOBE auto update svg transform.translate scaleY. failed by (SvgContens_drowSvg_transform_scaleY=spectreTiles["Delta"].transformations[0][0,0]) 

    # print(f"transformations={[transformations[i] for i in [6,5,3,0]]}")
    # Now build the actual supertiles, labeling appropriately.
    super_quad =  np.array([
        transPt(transformations[6], quad[2]),
        transPt(transformations[5], quad[1]),
        transPt(transformations[3], quad[2]),
        transPt(transformations[0], quad[1]) 
    ])
    # print(f"super_quad={super_quad}")

    tiles = {label: MetaTile(tiles=[input_tiles[subst] for subst in substitutions if subst],
                     transformations=[trsf for subst, trsf in zip(substitutions, transformations) if subst],
                     quad=super_quad
                     ) for label, substitutions in (
                         ("Gamma",  ("Pi",  "Delta", None,  "Theta", "Sigma", "Xi",  "Phi",    "Gamma")),
                         ("Delta",  ("Xi",  "Delta", "Xi",  "Phi",   "Sigma", "Pi",  "Phi",    "Gamma")),
                         ("Theta",  ("Psi", "Delta", "Pi",  "Phi",   "Sigma", "Pi",  "Phi",    "Gamma")),
                         ("Lambda", ("Psi", "Delta", "Xi",  "Phi",   "Sigma", "Pi",  "Phi",    "Gamma")),
                         ("Xi",     ("Psi", "Delta", "Pi",  "Phi",   "Sigma", "Psi", "Phi",    "Gamma")),
                         ("Pi",     ("Psi", "Delta", "Xi",  "Phi",   "Sigma", "Psi", "Phi",    "Gamma")),
                         ("Sigma",  ("Xi",  "Delta", "Xi",  "Phi",   "Sigma", "Pi",  "Lambda", "Gamma")),
                         ("Phi",    ("Psi", "Delta", "Psi", "Phi",   "Sigma", "Pi",  "Phi",    "Gamma")),
                         ("Psi",    ("Psi", "Delta", "Psi", "Phi",   "Sigma", "Psi", "Phi",    "Gamma"))
                      )}
    return tiles

transformation_min_X = np.inf
transformation_min_Y = np.inf
transformation_max_X = -np.inf
transformation_max_Y = -np.inf
def update_transformation_range(T, _label): # drowsvg
    """
    T: transformation matrix
    label: unused label string
    """
    global transformation_min_X, transformation_min_Y, transformation_max_X, transformation_max_Y
    transformation_min_X = min(transformation_min_X, T[0,2]) # drowsvg
    transformation_min_Y = min(transformation_min_Y, T[1,2]) # drowsvg
    transformation_max_X = max(transformation_max_X, T[0,2]) # drowsvg
    transformation_max_Y = max(transformation_max_Y, T[1,2]) # drowsvg
    return

#### main process ####
def buildSpectreTiles(n_ITERATIONS,edge_a,edge_b):
    global SPECTRE_POINTS, Mystic_SPECTRE_POINTS, SPECTRE_QUAD

    SPECTRE_POINTS = get_spectre_points(edge_a, edge_b) # tile(Edge_a, Edge_b)
    Mystic_SPECTRE_POINTS = get_spectre_points(edge_b, edge_a) # tile(Edge_b, Edge_a)
    SPECTRE_QUAD = SPECTRE_POINTS[[3,5,7,11],:]
    tiles = buildSpectreBase()
    for _ in range(n_ITERATIONS):
        tiles = buildSupertiles(tiles)

    tiles["Delta"].forEachTile(update_transformation_range) # scan all Tile

    global transformation_min_X, transformation_min_Y, transformation_max_X, transformation_max_Y
    transformation_min_X = int(np.floor(transformation_min_X - Edge_a * 3 - Edge_b * 3))
    transformation_min_Y = int(np.floor(transformation_min_Y - Edge_a * 3 - Edge_b * 3))
    transformation_max_X = int(np.ceil(transformation_max_X + Edge_a * 3 + Edge_b * 3))
    transformation_max_Y = int(np.ceil(transformation_max_Y + Edge_a * 3 + Edge_b * 3))

    return tiles


### drawing parameter data
# Color map from Figure 5.3
COLOR_MAP = {
	'Gamma': np.array((203, 157, 126),'f')/255.,
	'Gamma1': np.array((203, 157, 126),'f')/255.,
	'Gamma2': np.array((203, 157, 126),'f')/255.,
	'Delta': np.array((163, 150, 133),'f')/255.,
	'Theta': np.array((208, 215, 150),'f')/255.,
	'Lambda': np.array((184, 205, 178),'f')/255.,
	'Xi': np.array((211, 177, 144),'f')/255.,
	'Pi': np.array((218, 197, 161),'f')/255.,
	'Sigma': np.array((191, 146, 126),'f')/255.,
	'Phi': np.array((228, 213, 167),'f')/255.,
	'Psi': np.array((224, 223, 156),'f')/255.
}

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

trot_inv_prof = {
    # -180: 0, # to be 0, becaluse angle -180=>180
    -150: 0, # Gamma2
    -120: 0,
    -90: 0, # Gamma2
    -60: 0,
    -30: 0,  # Gamma2
    0: 0,
    30: 0,  # Gamma2
    60: 0,
    90: 0, # Gamma2
    120: 0,
    150: 0,  # Gamma2
    180: 0,
    360: 0 # Gamma2 total
}
def print_trot_inv_prof():
    global trot_inv_prof
    print("transformation rotation profile(angle: count)={")
    for angle, count in (sorted(trot_inv_prof.items())):
        print(f"\t{angle}: {count},")
    print("}")
    return trot_inv_prof

def get_color_array(tile_transformation, label):
    global trot_inv_prof
    angle, _scale = trot_inv(tile_transformation)
    trot_inv_prof[angle] += 1
    if (label == 'Gamma2'):
        trot_inv_prof[360] += 1
        return np.array([0.25,0.25,0.25])
    else :
        rgb = {
                # -180: (  0,   0, 1.0), # sangle -180 == 180
                -120: (0.9, 0.8,   0),
                -60:  (0.9, 0.4, 0.4),
                0:    (1.0,   0,   0),
                60:   (0.4, 0.4, 0.9),
                120:  (  0, 0.8, 0.9),
                180:  (  0,   0, 1.0)
        }[angle]
        if rgb:
            return np.array(rgb, 'f')
        else:
            print(f"Inalid color {rgb} {label}, {tile_transformation}")
    return COLOR_MAP[label]

