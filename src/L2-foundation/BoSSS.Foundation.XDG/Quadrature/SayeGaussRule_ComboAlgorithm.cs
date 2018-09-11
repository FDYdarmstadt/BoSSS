using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using ilPSP;

namespace BoSSS.Foundation.XDG.Quadrature
{
    public abstract class SayeComboIntegrand<S, T>
        : SayeIntegrand<S, T>
        where S : IPsi
        where T : SayeArgument<S>
        
    {
        /// <summary>
        /// Returns surface & volume quadrature nodes, respectively.
        /// </summary>
        /// <param name="Cell"></param>
        /// <param name="arg"></param>
        /// <returns></returns>
        public QuadRule[] ComboEvaluate(int Cell, T arg)
        {
            //Setup algorithm
            //----------------------------------------------------------------------
            arg.Surface = false;
            cell = Cell;
            TreeNode<T> recursionTree = new TreeNode<T>();
            TreeNode<T> fullSpace = recursionTree.AddChild(arg);

            //Build Integrand
            //----------------------------------------------------------------------
            SayeRecursion(fullSpace);

            //Evaluate Integrand
            //----------------------------------------------------------------------
            //Fill nodesAndWeights
            recursionTree.UnrollFunc(ComboIntegrandEvaluation);
            //ConvertToQuadRule
            QuadRule[] toRuleThemAll = GetComboQuadRule();
            return toRuleThemAll;
        }

        private QuadRule[] GetComboQuadRule()
        {
            QuadRule[] rulez = new QuadRule[2];
            rulez[0] = VolumeRule.GetQuadRule();
            rulez[1] = SurfRule.GetQuadRule();
            return rulez;
        }

        ISayeQuadRule SurfRule;
        ISayeQuadRule VolumeRule; 

        private void ComboIntegrandEvaluation(TreeNode<T> node)
        {
            if (node.level == 1)
            {
                //Calculate Surface and Volume Rule
                T nodeArg = node.Value;
                SayeQuadRule newRule;
                switch (nodeArg.Status)
                {
                    case SayeArgument<S>.Mode.Standart:
                        SetStandartComboNodes(node);
                        break;
                    case SayeArgument<S>.Mode.GaussQuadrature:
                        //Surface must be 0;
                        throw new NotImplementedException();
                        break;
                    case SayeArgument<S>.Mode.LowOrderQuadrature:
                        throw new NotImplementedException();
                        break;
                    case SayeArgument<S>.Mode.DomainIsEmpty:
                        break;
                    default:
                        throw new NotSupportedException();
                }
            }
            else
            {
                IntegrandEvaluation(node);
            }
        }

        private void SetStandartComboNodes(TreeNode<T> node)
        {
            T arg = node.Value;

            foreach (TreeNode<T> childNode in node.Children)
            {
                T childArg = childNode.Value;
                foreach (Tuple<MultidimensionalArray, double>
                    integrationNode in childArg.NodesAndWeights.IntegrationNodes)
                {
                    ComboIntegrand(integrationNode.Item1, integrationNode.Item2, arg);
                }
            }
        }

        private void ComboIntegrand(MultidimensionalArray X, double X_weight, T arg)
        {
            //Calculate the set of roots and union with the interval endpoints ( line 1)
            int heightDirection = arg.HeightDirection;

            //Sort R into ascending order such that ... (line 2) : Implemented with a sortedList
            ISortedBoundedList<double> roots = new SayeSortedList(arg.n * 2);
            double[] bounds = GetBoundaries(arg, heightDirection);
            roots.SetBounds(bounds[0], bounds[1]);
            RestrictToBound(X, bounds[0], arg.HeightDirection);
            for (int i = 0; i < arg.n; ++i)
            {
                S psi_i = arg.PsiAndS[i].Item1;
                double[] rootsInPsi = FindRoots(psi_i, X, heightDirection, bounds, this.cell);
                roots.Add(rootsInPsi);
            }

            //For j = 1 to l - 1 do (line 4)
            bool xIsUnchanged = true;
            for (int j = 0; j < roots.Count() - 1; ++j)
            {
                //Define L and x_c(line 5)
                double L = roots[j + 1] - roots[j];
                NodeSet x_c = NodeOnRay(X, heightDirection, roots[j] - roots[0] + L / 2.0);

                //If s_i * psi_i >= 0 for all i (line 6)
                bool updateIntegrand = true;
                foreach (Tuple<S, int> psiAndS in arg.PsiAndS)
                {
                    S psi_i = psiAndS.Item1;
                    int s_i = psiAndS.Item2;
                    updateIntegrand &= s_i * EvaluateAt(psi_i, x_c, cell) >= 0;
                }
                //Update I = ...(line 7)
                if (updateIntegrand)
                {
                    SayeQuadRule newRule = BuildQuadRule(x_c, X_weight, heightDirection, L);
                    bool deriveFromExistingNode = xIsUnchanged;
                    arg.NodesAndWeights.AddRule(newRule, deriveFromExistingNode);
                    xIsUnchanged = false;
                }
            }
            if (xIsUnchanged)
            {
                arg.NodesAndWeights.RemoveActiveNode();
            }
        }
    }
}
