using BoSSS.Platform.LinAlg;
using ilPSP;
using System;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Recomposer
{
    class BoundaryTransformation : Transformation
    {
        public BoundaryTransformation(BoundaryLine source, BoundaryLine target)
            : base(2)
        {
            if(source.Length() != target.Length())
            {
                throw new Exception("There is something rotten in the state of the boundary lengths");
            }

            this.Matrix = GetRotationMatrixFrom(
                source.Start.Position - source.End.Position, 
                target.Start.Position - target.End.Position);
            this.AffineTransformation = 
                (target.Start.Position + target.End.Position - source.Start.Position - source.End.Position) / 2;
        }

        static MultidimensionalArray GetRotationMatrixFrom(Vector source, Vector target)
        {
            source.Normalize();
            target.Normalize();
            return GetRotationMatrixFromNormal(source, target); 
        }

        //https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        static MultidimensionalArray GetRotationMatrixFromNormal(Vector source, Vector target)
        {
            double crossProduct = source.CrossProduct2D(target);
            double dotProduct = source * target;

            MultidimensionalArray rotation = MultidimensionalArray.Create(2, 2);
            rotation[0, 0] = dotProduct;
            rotation[0, 1] = -crossProduct;
            rotation[1, 0] = dotProduct;
            rotation[1, 1] = crossProduct;
            return rotation;
        }
    }
}
