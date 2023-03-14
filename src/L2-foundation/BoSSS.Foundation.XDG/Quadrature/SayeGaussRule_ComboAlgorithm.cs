﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using ilPSP;

namespace BoSSS.Foundation.XDG.Quadrature
{
    abstract class SayeComboIntegrand<S, T>
        : SayeIntegrand<S, T>
        where T : SayeArgument<S>    
    {
        //The surface quadrature rule is saved here.
        ISayeQuadRule SurfRule;

        /// <summary>
        /// Returns surface and volume quadrature nodes, respectively.
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
            SayeRecursion(fullSpace,0);

            //Evaluate Integrand
            //----------------------------------------------------------------------
            //Fill nodesAndWeights

            SurfRule = GetEmptyQuadratureRule();
            recursionTree.UnrollFunc(ComboIntegrandEvaluation);

            //Convert data to QuadRule. The volume rule is saved in SayeArgument as in the single algorithm. 
            //The surface rule is saved SurfRule.
            ISayeQuadRule ruleZ = GetEmptyQuadratureRule();
            foreach (TreeNode<T> childNode in recursionTree.Children) {
                foreach (SayeQuadRule rule in childNode.Value.NodesAndWeights.GetRules()) {
                    ruleZ.AddRule(rule);
                }
            }

            QuadRule[] rulezOfKrom = new QuadRule[2];
            rulezOfKrom[0] = ruleZ.GetQuadRule(); //Volume rule
            rulezOfKrom[1] = SurfRule.GetQuadRule();    //Surface rule
            
            //Handle empty QuadRule
            if (rulezOfKrom[0].NoOfNodes == 0)
            {
                Debug.WriteLine("Evaluation of cell {0} returned empty volume quadrule", cell);
                rulezOfKrom[0] = CreateZeroQuadrule();
            }
            if(rulezOfKrom[1].NoOfNodes == 0)
            {
                Debug.WriteLine("Evaluation of cell {0} returned empty surface quadrule", cell);
                rulezOfKrom[1] = CreateZeroQuadrule();
            }


            return rulezOfKrom;
        }
        
        private void ComboIntegrandEvaluation(TreeNode<T> node)
        {
            if (node.level == 1)
            {
                //Calculate Surface and Volume Rule
                T nodeArg = node.Value;
                switch (nodeArg.Status)
                {
                    case SayeArgument<S>.Mode.Standard:
                        SetStandartComboNodes(node);
                        break;
                    case SayeArgument<S>.Mode.GaussQuadrature:
                        //No need for surface quadNodes
                        SayeQuadRule newRule = SetGaussQuadratureNodes(nodeArg);
                        nodeArg.NodesAndWeights.AddRule(newRule);
                        break;
                    case SayeArgument<S>.Mode.LowOrderQuadrature:
                        SetLowOrderNodes(node);
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

        private void SetLowOrderNodes(TreeNode<T> node) {
            T arg = node.Value;
            NodeSet c = arg.GetCellCenter();
            SayeQuadRule surf = BuildSurfaceQuadRule(c, 2, arg.HeightDirection, cell);
            SurfRule.AddRule(surf);
            SayeQuadRule volume =  SetLowOrderQuadratureNodes(arg);
            arg.NodesAndWeights.AddRule(volume);
        }

        //Evaluate volume and surface integrand
        private void ComboIntegrand(MultidimensionalArray X, double X_weight, T arg)
        {
            //Calculate the set of roots and union with the interval endpoints ( line 1)
            int heightDirection = arg.HeightDirection;

            //Sort R into ascending order such that ... (line 2) : Implemented with a sortedList
            ISortedBoundedList<double> roots = new SayeSortedList(2 + arg.n, 1.0e-14);
            double[] bounds = GetBoundaries(arg, heightDirection);
            roots.SetBounds(bounds[0], bounds[1]);
            RestrictToBound(X, bounds[0], arg.HeightDirection);
            for (int i = 0; i < arg.n; ++i)
            {
                S psi_i = arg.PsiAndS[i].Item1;
                double[] rootsInPsi = FindRoots(psi_i, X, heightDirection, bounds, this.cell);
                roots.Add(rootsInPsi);
            }

            //Volume evaluation
            //-----------------------------------------------------------------------------------
            //For j = 1 to l - 1 do (line 4)
            for (int j = 0; j < roots.Count() - 1; ++j)
            {
                //Define L and x_c(line 5)
                double L = roots[j + 1] - roots[j];
                NodeSet x_c = NodeOnRay(X, heightDirection, roots[j] - bounds[0] + L / 2.0);

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
                    //Update 
                    SayeQuadRule newRule = BuildQuadRule(x_c, X_weight, heightDirection, L);
                    arg.NodesAndWeights.AddRule(newRule);
                }
            }
            //Surface evaluation
            //-----------------------------------------------------------------------------------
            for (int i = 1; i < roots.Count() - 1; ++i) 
            {
                X[0, heightDirection] = roots[i];
                SayeQuadRule surfaceQuadNode = BuildSurfaceQuadRule(X, X_weight, heightDirection, this.cell);
                SurfRule.AddRule(surfaceQuadNode);
            }
        }
    }
}