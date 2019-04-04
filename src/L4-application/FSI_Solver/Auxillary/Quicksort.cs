using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FSI_Solver
{
    class Quicksort // copied from Wikipedia
    {
        internal static void Main(int Leftelement, int RightElement, ref int[] Data)
        {
            if (Leftelement < RightElement)
            {
                int division = Divide(Leftelement, RightElement, ref Data);
                Main(Leftelement, division - 1, ref Data);
                Main(division + 1, RightElement, ref Data);
            }
        }

        private static int Divide(int Leftelement, int RightElement, ref int[] Data)
        {
            int i = Leftelement;
            int j = RightElement - 1;
            int pivot = Data[RightElement];

            do
            {
                while (Data[i] <= pivot && i < RightElement)
                    i += 1;

                while (Data[j] >= pivot && j > Leftelement)
                    j -= 1;

                if (i < j)
                {
                    int z = Data[i];
                    Data[i] = Data[j];
                    Data[j] = z;
                }

            } while (i < j);

            if (Data[i] > pivot)
            {
                int z = Data[i];
                Data[i] = Data[RightElement];
                Data[RightElement] = z;
            }
            return i; 
        }
    }
}
