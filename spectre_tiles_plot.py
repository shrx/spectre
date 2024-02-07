# draw Polygons Svg by matplotlib #####
from spectre import buildSpectreTiles,get_color_array, SPECTRE_POINTS, Mystic_SPECTRE_POINTS, Edge_a,Edge_b, N_ITERATIONS, print_trot_inv_prof
from time import time
import matplotlib.pyplot as plt

start = time()
spectreTiles = buildSpectreTiles(N_ITERATIONS,Edge_a,Edge_b)
time1 = time()-start

print(f"supertiling loop took {round(time1, 4)} seconds")

start = time()
plt.figure(figsize=(8, 8))
plt.axis('equal')

num_tiles = 0
def plotVertices(tile_transformation, label):
    """
    T: transformation matrix
    label: label of shape type
    """
    global num_tiles
    num_tiles += 1
    vertices = (SPECTRE_POINTS if label != "Gamma2" else Mystic_SPECTRE_POINTS).dot(tile_transformation[:,:2].T) + tile_transformation[:,2]
    color_array = get_color_array(tile_transformation, label)
    # plt.text((vertices[1,0] + vertices[7,0])/2, (vertices[1,1] + vertices[7,1])/2, label, fontsize=8, color='gray')
    plt.fill(vertices[:,0],vertices[:,1],facecolor=color_array)
    plt.plot(vertices[:,0],vertices[:,1],color='gray',linewidth=0.2)

spectreTiles["Delta"].drawPolygon(plotVertices)
time2 = time()-start
print(f"matplotlib.pyplot: tile recursion loop took {round(time2, 4)} seconds, generated {num_tiles} tiles")
print_trot_inv_prof()

start = time()
saveFileName = f"spectre_tile{Edge_a:.1f}-{Edge_b:.1f}_{N_ITERATIONS}-{num_tiles}pts.svg"
print("matplotlib.pyplot: file save to " + saveFileName)
plt.savefig(saveFileName)
time3 = time()-start
print(f"matplotlib.pyplot SVG drawing took {round(time3, 4)} seconds")
print(f"matplotlib.pyplot total processing time {round(time1+time2+time3, 4)} seconds, {round(1000000*(time1+time2+time3)/num_tiles, 4)} μs/tile")

plt.show()
