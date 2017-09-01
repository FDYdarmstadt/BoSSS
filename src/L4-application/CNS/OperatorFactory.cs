using System;
using BoSSS.Foundation;
using CNS.BoundaryCondition;

namespace CNS {

    abstract class OperatorFactory {

        protected CNSControl control;

        protected BoundaryConditionMap boundaryMap;

        protected string[] argumentOrdering;

        public OperatorFactory(CNSControl control, BoundaryConditionMap boundaryMap) {
            this.control = control;
            this.boundaryMap = boundaryMap;
            argumentOrdering = new string[control.Config.VariableOrdering.Count];
            control.Config.VariableOrdering.Keys.CopyTo(argumentOrdering, 0);
        }

        public SpatialDifferentialOperator GetConvectiveOperator() {
            SpatialDifferentialOperator op =
                new SpatialDifferentialOperator(argumentOrdering, argumentOrdering);
            AddConvectiveComponents(op);
            op.Commit();

            return op;
        }

        public SpatialDifferentialOperator GetDiffusiveOperator() {
            SpatialDifferentialOperator op =
                new SpatialDifferentialOperator(argumentOrdering, argumentOrdering);
            AddDiffusiveComponents(op);
            op.Commit();

            return op;
        }

        public SpatialDifferentialOperator GetCombinedOperator() {
            SpatialDifferentialOperator op =
                new SpatialDifferentialOperator(argumentOrdering, argumentOrdering);
            AddConvectiveComponents(op);
            AddDiffusiveComponents(op);
            op.Commit();

            return op;
        }

        abstract protected void AddConvectiveComponents(SpatialDifferentialOperator op);

        abstract protected void AddDiffusiveComponents(SpatialDifferentialOperator op);
    }
}
