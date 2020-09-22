using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using static BoSSS.Foundation.XDG.LevelSetTracker;
using System.Threading;
using BoSSS.Platform;
using System.ComponentModel;
using System.Collections;
using BoSSS.Foundation.Quadrature;

namespace BoSSS.Foundation.XDG.Quadrature
{
    struct CombinedID
    {
        public int LevSet0;
        public JumpTypes Jmp0;
        public int LevSet1;
        public JumpTypes Jmp1;

        public bool Equals(CombinedID otherID)
        {
            if ((LevSet0 == otherID.LevSet0)
                && (Jmp0 == otherID.Jmp0)
                && (LevSet1 == otherID.LevSet1)
                && (Jmp1 == otherID.Jmp1))
            {
                return true;
            }
            else if ((LevSet0 == otherID.LevSet1)
                && (Jmp0 == otherID.Jmp1)
                && (LevSet1 == otherID.LevSet0)
                && (Jmp1 == otherID.Jmp0))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }

    class CombinedLevelSet
    {
        CombinedID iD;

        public LevelSet LevelSet;

        public CombinedLevelSet(LevelSet levelSet, CombinedID iD)
        {
            this.LevelSet = levelSet;
            this.iD = iD;
        }

        public bool Equals(CombinedID otherID)
        {
            return iD.Equals(otherID);
        }
    }

    class CombinedLevelSetData
    {
        int degree;

        GridData gridData;

        LevelSetData[] originalLevelSetDatas;

        CombinedLevelSet[] combinedLevelSets;
        
        LevelSetData[] combinedLevelSetDatas;

        /// <summary>
        /// For now only 2 combined LevelSets
        /// </summary>
        /// <param name="datas"></param>
        /// <param name="degree"></param>
        public CombinedLevelSetData(LevelSetData[] datas, int degree)
        {
            this.degree = degree;
            this.originalLevelSetDatas = datas;
            gridData = datas[0].GridDat;
            
            CombineOriginalLevelSets();
            ExtractLevelSetData();
        }

        void CombineOriginalLevelSets()
        {
            //--A | -+ B
            //___________
            //+-D | ++C
            combinedLevelSets = new CombinedLevelSet[4];
            combinedLevelSets[0] = Combine(0, JumpTypes.OneMinusHeaviside, 1, JumpTypes.OneMinusHeaviside, "A");
            combinedLevelSets[1] = Combine(0, JumpTypes.OneMinusHeaviside, 1, JumpTypes.Heaviside, "B");
            combinedLevelSets[2] = Combine(0, JumpTypes.Heaviside, 1, JumpTypes.Heaviside, "C");
            combinedLevelSets[3] = Combine(0, JumpTypes.Heaviside, 1, JumpTypes.OneMinusHeaviside, "D");
        }

        CombinedLevelSet Combine(int levSetIndice0, JumpTypes jmp0, int levSetIndice1, JumpTypes jmp1, string name) 
        {
            LevelSet levSet0 = originalLevelSetDatas[levSetIndice0].LevelSet.As<LevelSet>();
            LevelSet levSet1 = originalLevelSetDatas[levSetIndice1].LevelSet.As<LevelSet>();

            LevelSet combination = CombineLevelSet(levSet0, jmp0, levSet1, jmp1, name);
            CombinedID id = new CombinedID
            {
                LevSet0 = levSetIndice0,
                Jmp0 = jmp0,
                LevSet1 = levSetIndice1,
                Jmp1 = jmp1
            };
            return new CombinedLevelSet(combination, id);
        }

        LevelSet CombineLevelSet(LevelSet levSet0, JumpTypes jmp0, LevelSet levSet1, JumpTypes jmp1, string name)
        {
            double sign0 = ToDouble(jmp0);
            double sign1 = ToDouble(jmp1);
            double Min(Vector position, double[] levSetValues, int cellIndex)
            {
               return Math.Min(sign0 * levSetValues[0], sign1 * levSetValues[1]);
            };

            LevelSet combination = new LevelSet(new Basis(gridData, degree), name);
            combination.Clear();
            combination.ProjectFunction(1, Min, new CellQuadratureScheme(), levSet0, levSet1);
            return combination;
        }

        double ToDouble(JumpTypes type)
        {
            if (type == JumpTypes.Heaviside)
            {
                return 1.0;
            }
            else if (type == JumpTypes.OneMinusHeaviside)
            {
                return -1.0;
            }
            else
            {
                throw new NotSupportedException();
            }
        }

        void ExtractLevelSetData()
        {
            //Hack di Hack
            string[,,,] speciesTable = new string[2,2,2,2];
            for (int i = 0; i < 16; ++i)
            {
                BitArray j = new BitArray(new int[] { i });
                int[] bits = j.Cast<bool>().Select(bit => bit ? 1 : 0).ToArray();
                speciesTable[bits[0], bits[1], bits[2], bits[3]] = "wayne";
            }
            LevelSetTracker tracker = new LevelSetTracker(gridData, 
                XQuadFactoryHelper.MomentFittingVariants.Saye, 
                0, 
                speciesTable, 
                combinedLevelSets[0].LevelSet,
                combinedLevelSets[1].LevelSet,
                combinedLevelSets[2].LevelSet,
                combinedLevelSets[3].LevelSet);
            tracker.UpdateTracker();
            tracker.PushStacks();
            combinedLevelSetDatas = tracker.DataHistories.Select(hist => hist[0]).ToArray();
        }

        public LevelSetData[] GetCombinedLevelSetDatas()
        {
            return combinedLevelSetDatas;
        }

        public int GetCombinedLevelSetIndice(int levSet0, JumpTypes jmp0, int levSet1, JumpTypes jmp1)
        {
            CombinedID id = new CombinedID
            {
                LevSet0 = levSet0,
                Jmp0 = jmp0,
                LevSet1 = levSet1,
                Jmp1 = jmp1
            };
            for(int i = 0; i < combinedLevelSets.Count(); ++i)
            {
                if (combinedLevelSets[i].Equals(id))
                {
                    return i;
                }
            }
            throw new Exception("combination of LevelSets was not found");
        }
    }
}
