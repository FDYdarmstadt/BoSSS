using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid
{
    class GridComparer<T> : IEqualityComparer<IGrid>
        where T : IGrid
    {
        readonly Func<T, T, bool> CheckEquality;

        public GridComparer(Func<T, T, bool> CheckEquality)
        {
            this.CheckEquality = CheckEquality;
        }

        public bool Equals(IGrid x, IGrid y)
        {
            bool isEqual;
            if (x is T X && y is T Y)
            {
                isEqual = CheckEquality(X, Y);
            }
            else
            {
                isEqual = false;
            }
            return isEqual;
        }

        public int GetHashCode(IGrid obj)
        {
            throw new NotImplementedException();
        }
    }
}
