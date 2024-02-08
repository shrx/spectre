Python script for generating tilings of the weakly chiral aperiodic monotile Tile(1,1) "Spectre".

Code ported from JavaScript from the web app [1] provided [2] by the authors of the original research paper [3].

[1]: https://cs.uwaterloo.ca/~csk/spectre/app.html

[2]: https://cs.uwaterloo.ca/~csk/spectre/

[3]: https://arxiv.org/abs/2305.17743

![Rendered tiling.](./spectre.svg)


* USAGE

    * When drawing with drowsvg the command is : 
       ```python spectre_tiles_drow.py```
    * When drawing with mathplot.plot, the command is : 
       ```python spectre_tiles_plot.py```
    * when customization;
        To ensure that the same pattern is visible no matter which command you use to draw the spectre tile,
        the customization related to the drawing is embedded in the ```spectre.py```

* CHANGES

    * Made it possible to compare the drawing speed between the path drawing process of all polygons by mathplotlib and the two polygon reference processes via transform by drowsvg.
    * Made it possible to draw spectre tile(edge_a, edge_b) at any ratio.
    * split mathplot.plot and drowsvg
    * In order to reduce the size of the SVG file, 
      the Transform of DrawSVG replaced the matrix with 6 floating-point numbers 
      with a translate with 2 floating-point numbers and a rotate and scale expansion with 3 integers. 


![Rendered tiling ratio sqrt(3)  tile(7.3, 12.7)](./spectre_tile7.3-12.7_3-559useRef.svg)
