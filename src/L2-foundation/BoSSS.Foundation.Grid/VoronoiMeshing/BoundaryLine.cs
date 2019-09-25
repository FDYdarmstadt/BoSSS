using System;
using System.Collections.Generic;
using BoSSS.Platform.LinAlg;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class BoundaryLine : Line
    {
        public double Length()
        {
            double length = (End.Position - Start.Position).Abs();
            return length;
        }

        public static BoundaryLine[] ToLines( Vector[] polygon)
        {
            BoundaryLine[] lines = new BoundaryLine[polygon.Length];
            BoundaryLine line;
            for (int i = 0; i < polygon.Length - 1; ++i)
            {
                line = new BoundaryLine
                {
                    Start = new Vertex()
                    {
                        Position = polygon[i]
                    },
                    End = new Vertex()
                    {
                        Position = polygon[i + 1]
                    },
                };
                lines[i] = line;
            }
            line = new BoundaryLine { Start = new Vertex(), End = new Vertex() };
            line.Start.Position = polygon[polygon.Length - 1];
            line.End.Position = polygon[0];
            lines[lines.Length - 1] = line;

            return lines;
        }

        public static BoundaryLine GetReverse(BoundaryLine line)
        {
            BoundaryLine reverse = new BoundaryLine
            {
                Start = line.End,
                End = line.Start
            };
            return reverse;
        }
    }
}
