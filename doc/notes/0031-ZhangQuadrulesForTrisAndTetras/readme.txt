
quadrature rules from the PHG code, http://lsec.cc.ac.cn/phg/index_en.htm;

 @article{zhang2009set,
 title={A set of symmetric quadrature rules on triangles and tetrahedra},
 author={Zhang, Linbo and Cui, Tao and Liu, Hui and others},
 journal={J. Comput. Math},
 volume={27},
 number={1},
 pages={89--96},
 year={2009}
 }

Code wurde mit __float128 kompiliert.
Die Tetra-Regeln scheinen OK zu sein; die Dreieck-Regeln haben irgendeinen Bug, wie es scheint.

Von mir modifiziert: quad_test.c, schreibt nun BoSSS-Quadregeln in Textfile.

