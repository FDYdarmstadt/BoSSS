using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Platform.LinAlg;

namespace VoronoiTests.Grid
{
    static class GridShapes
    {
        public static Vector[] LShape()
        {
            double a = 1;
            Vector[] LShapedPolygon = new[] {
                    new Vector(-a,a),
                    new Vector(a,a),
                    new Vector(a,-a),
                    new Vector(0,-a),
                    new Vector(0,0),
                    new Vector(-a,0)
                };
            return LShapedPolygon;
        }

        public static Vector[] Rectangle(double width, double height)
        {
            Vector[] polygonBoundary = new Vector[]
            {
                new Vector(-width / 2, height / 2),
                new Vector(width / 2, height / 2),
                new Vector(width / 2, -height / 2),
                new Vector(-width / 2, -height / 2)
            };
            return polygonBoundary;
        }
    }
}
