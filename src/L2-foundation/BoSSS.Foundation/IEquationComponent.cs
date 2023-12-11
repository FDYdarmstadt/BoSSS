/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System.Collections.Generic;
using BoSSS.Platform;
using System;
using BoSSS.Foundation.Grid;
using ilPSP;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;

namespace BoSSS.Foundation {

    /// <summary>
    /// Properties that all kinds of objects that form equations have in common.
    /// </summary>
    public interface IEquationComponent {
        /// <summary>
        /// Defines the order in which numerical values are provided as Arguments for the
        /// flux functions (domain variables).
        /// E.g., if this property returns the array {"u","v","w"},
        /// the values of u,v,w  are put into into the <c>Uin</c>,<c>Uout</c>-arguments 
        /// of the flux function when they are called.
        /// </summary>
        IList<string> ArgumentOrdering { get; }

        /// <summary>
        /// Defines the order in which numerical values of the 
        /// parameter variables are sorted for this flux;
        /// </summary>
        IList<string> ParameterOrdering { get; }
    }


    /// <summary>
    /// imposes constraints (e.g. DG polynomial degree grater than 0) on the DG polynomial degree of a discretization
    /// </summary>
    public interface IDGdegreeConstraint : IEquationComponent {

        /// <summary>
        /// Indicate wheter a specfic combination of DG degrees is a valid combination (e.g. with respect to numerical stability, a degree 0 is invalid for an SIP discretization)
        /// </summary>
        /// <param name="DomainDegreesPerVariable">
        /// - index correlates with <see cref="IEquationComponent.ArgumentOrdering"/>
        /// - content: DG polynomial degree for respective variable
        /// </param>
        /// <param name="CodomainDegree">
        ///  DG polynomial degree for codomain respective variable
        /// </param>        
        /// <returns>
        /// - true: combination of DG degrees is OK;
        /// - false if not
        /// </returns>
        bool IsValidDomainDegreeCombination(int[] DomainDegreesPerVariable, int CodomainDegree);
    }

    /// <summary>
    /// Interface for equation components that require parameters which are based on the current state 
    /// </summary>
    public interface IParameterHandling : IEquationComponent {

        /// <summary>
        /// Update of parameter fields used by this operator;
        /// (alternatively, <see cref="IDifferentialOperator.ParameterUpdates"/> can be used.) 
        /// </summary>
        /// <param name="Arguments">
        /// input, the current state of the approximate solution at which the operator should be evaluated or linearized.
        /// sequence correlates with <see cref="IEquationComponent.ArgumentOrdering"/>.
        /// </param>
        /// <param name="Parameters">
        /// output, sequence correlates with <see cref="IEquationComponent.ParameterOrdering"/>.
        /// </param>
        void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters);

        /// <summary>
        /// Factory for the allocation of storage for storing the parameters for this component
        /// (alternatively, <see cref="IDifferentialOperator.ParameterFactories"/> can be used.)
        /// </summary>
        /// <param name="Arguments">
        /// input, the currently used DG fields to store approximate solution.
        /// (required if e.g. the DG degree of the parameters somehow should depend on the DG degree of the arguments/domain variables).
        /// sequence correlates with <see cref="IEquationComponent.ArgumentOrdering"/>.
        /// </param>
        /// <returns></returns>
        DGField[] MyParameterAlloc(DGField[] Arguments);
    }

    /// <summary>
    /// Additional hints for checking the implementation
    /// </summary>
    public interface IEquationComponentChecking : IEquationComponent { 
        
        /// <summary>
        /// Only for performance measurements of the vectorized implementations (e.g. <see cref="INonlinVolumeForm_V"/>, <see cref="IVolumeForm_UxV"/>, etc.)
        /// against the base implementation (<see cref="IVolumeForm"/>, <see cref="IEdgeForm"/>):
        /// If true, the base implementation will be used even if a vectorized version is provided.
        /// </summary>
        bool IgnoreVectorizedImplementation { get; }
    }


    /// <summary>
    /// Interface for equation components which require e.g. grid and/or problem-dependent coefficients,
    /// e.g. cell length scales;
    /// </summary>
    public interface IEquationComponentCoefficient : IEquationComponent {

        /// <summary>
        /// Passes various coefficients to the equation component.
        /// </summary>
        /// <param name="cs">
        /// </param>
        /// <param name="DomainDGdeg">
        /// DG polynomial order of trial/domain variables/arguments; ordering corresponds with <see cref="IEquationComponent.ArgumentOrdering"/>;
        /// </param>
        /// <param name="TestDGdeg">
        /// DG polynomial degree of test/codomain variable
        /// </param>
        void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg);
    }
    
    /// <summary>
    /// Set of various custom resp. predefined coefficients, 
    /// </summary>
    /// <remarks>
    /// By encapsulation of the arguments of <see cref="IEquationComponentCoefficient.CoefficientUpdate"/> in a separate
    /// class it is easy tho change/add variables without updating each and every interface implementation.
    /// </remarks>
    public class CoefficientSet {

        /// <summary>
        /// Reference to grid
        /// </summary>
        public IGridData GrdDat;

        /// <summary>
        /// length scales for cells (e.g. for computing penalty parameters or local CFL numbers),
        /// typically (in the DG case) the values from <see cref="Grid.Classic.GridData.CellData.CellLengthScale"/>
        /// </summary>
        public MultidimensionalArray CellLengthScales;

        /// <summary>
        /// length scales for edges (e.g. for computing penalty parameters)
        /// typically (in the DG case) the values from <see cref="Grid.Classic.GridData.EdgeData.h_min_Edge"/>
        /// </summary>
        public MultidimensionalArray EdgeLengthScales;

        /// <summary>
        /// scalar on the homotopy path a nonlinear solver
        /// </summary>
        public double HomotopyValue;

        /// <summary>
        /// Species Subgrid
        /// </summary>
        public System.Collections.BitArray SpeciesSubGrdMask;

        /// <summary>
        /// collection of user-defined objects
        /// </summary>
        public Dictionary<string, object> UserDefinedValues = new Dictionary<string, object>();


    }


    /// <summary>
    /// Interface for components that provide their derivative,
    /// required for Newton solvers
    /// (see <see cref="DifferentialOperator.GetJacobiOperator"/>).
    /// </summary>
    public interface ISupportsJacobianComponent {

        /// <summary>
        /// A collection of components which in sum for the derivative of this component;
        /// For (bi-) linear components, this is usually the flux itself.
        /// </summary>
        IEquationComponent[] GetJacobianComponents(int SpatialDimension);
    }

    /// <summary>
    /// defines a nonlinear source term.
    /// ```math 
    /// s(\vec{U},v) = \int_{\Omega} f(\vec{U}) v \ \mathrm{dV}
    /// ```
    /// </summary>
    public interface INonlinearSource : IEquationComponent {
        
        /// <summary>
        /// the point-wise source term <em>f</em>;
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x">
        /// global coordinated of quadrature nodes for each cell
        /// 1st index: cell index
        /// 2nd index: quadrature node index
        /// 3rd index: spatial dimension
        /// </param>
        /// <param name="U"></param>
        /// <param name="IndexOffset"></param>
        /// <param name="Lenght"></param>
        /// <param name="Output">
        /// 
        /// </param>
        /// <param name="FirstCellInd">
        /// index of the first of the <paramref name="Lenght"/> cells to process
        /// </param>
        void Source(double time, 
                    MultidimensionalArray x, 
                    MultidimensionalArray[] U,
                    int IndexOffset, int FirstCellInd, int Lenght,
                    MultidimensionalArray Output);
    }

    /// <summary>
    /// specifies a nonlinear flux function; Equations 
    /// containing that kind of flux can only be treated explicit.
    /// </summary>
    public interface INonlinearFlux : IEquationComponent {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="time">actual time</param>
        /// <param name="x">
        /// 1st index: local edge index, with some offset;
        /// 2nd index: node index;
        /// 3rd index: spatial dimension index
        /// </param>
        /// <param name="normal">
        /// 1st index: edge index;
        /// 2nd index: quadrature node´s index
        /// 3rd index: spatial dimension index;
        /// </param>
        /// <param name="Uin">
        /// 1st index: edge index, with some offset;
        /// 2nd index: node index;
        /// </param>
        /// <param name="IndexOffset">
        /// an index offset into the first dimension of 
        /// <paramref name="Output"/>, <paramref name="x"/> and each <paramref name="Uin"/>-entry;
        /// </param>
        /// <param name="Lenght">
        /// number of edges to process
        /// </param>
        /// <param name="Output">
        /// accumulator (!) for output values;
        /// 1st index: local edge index, with some offset;
        /// 2nd index: node index;
        /// </param>
        /// <param name="EdgeTags"></param>
        /// <param name="EdgeTagsOffset"></param>
        /// <param name="normalFliped"></param>
        /// <param name="jEdge">index of first edge</param>
        void BorderEdgeFlux(double time, int jEdge,
                            MultidimensionalArray x, 
                            MultidimensionalArray normal, bool normalFliped,
                            byte[] EdgeTags, int EdgeTagsOffset,
                            MultidimensionalArray[] Uin,
                            int IndexOffset, int Lenght,
                            MultidimensionalArray Output);



        /// <summary>
        /// 
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="normal">
        /// edge normal's;
        /// 1st index: local edge index, with some offset;
        /// 2nd index: node index;
        /// 3rd index: spatial dimension index
        /// </param>
        /// <param name="Uin"></param>
        /// <param name="Uout"></param>
        /// <param name="Offset"></param>
        /// <param name="Lenght"></param>
        /// <param name="Output"></param>
        /// <param name="jEdge">index of first edge</param>
        void InnerEdgeFlux(double time, int jEdge,
                           MultidimensionalArray x,
                           MultidimensionalArray normal, 
                           MultidimensionalArray[] Uin,
                           MultidimensionalArray[] Uout,
                           int Offset, int Lenght,
                           MultidimensionalArray Output);


        /// <summary>
        /// 
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x">
        /// nodes position in global space;
        /// 1st index: cell index with some offset;
        /// 2nd index: node index within cell;
        /// 3rd index: spatial dimension;
        /// </param>
        /// <param name="U">
        /// 1st index: Field index as defined by 
        /// 2nd index: cell index with some offset;
        /// 3rd index: node index within cell;
        /// </param>
        /// <param name="Offset">
        /// 1st cell index into <paramref name="x"/> that must be evaluated.
        /// Evaluation has to be done only for the cells from index <paramref name="Offset"/> to 
        /// <paramref name="Offset"/>+<paramref name="Length"/>-1;
        /// </param>
        /// <param name="Length">
        /// Number of cells to evaluate; NOT EQUAL to the 1st length of <paramref name="x"/> or <paramref name="Output"/>!!!
        /// </param>
        /// <param name="Output">
        /// On exit, the flux vector field value for the corresponding nodes in <paramref name="x"/> should be 
        /// <em>accumulated</em> here;<br/>
        /// 1st index: cell index with some offset;
        /// 2nd index: node index;
        /// 3rd index: spatial dimension;
        /// All indices correspond with <paramref name="x"/>;
        /// </param>
        void Flux(double time, 
                  MultidimensionalArray x,
                  MultidimensionalArray[] U,
                  int Offset, int Length,
                  MultidimensionalArray Output);
    }


    /// <summary>
    /// specifies a nonlinear flux function; Equations 
    /// containing that kind of flux can only be treated explicit.
    /// It differs from <see cref="INonlinearFlux"/>
    /// </summary>
    public interface INonlinearFluxEx : IEquationComponent {


        /// <summary>
        /// Numerical flux at edges on the boundary of the physical domain
        /// </summary>
        /// <param name="time">actual time</param>
        /// <param name="x">
        /// 1st index: local edge index, with some offset;
        /// 2nd index: node index;
        /// 3rd index: spatial dimension index
        /// </param>
        /// <param name="normal">
        /// 1st index: local edge index, with some offset;
        /// 2nd index: node index;
        /// 3rd index: spatial dimension index
        /// </param>
        /// <param name="Uin">
        /// 1st index: edge index, with some offset;
        /// 2nd index: node index;
        /// </param>
        /// <param name="UinMean"></param>
        /// <param name="IndexOffset">
        /// an index offset into the first dimension of 
        /// <paramref name="Output"/>, <paramref name="x"/> and each <paramref name="Uin"/>-entry;
        /// </param>
        /// <param name="Lenght">
        /// number of edges to process
        /// </param>
        /// <param name="Output">
        /// accumulator (!) for output values;
        /// 1st index: local edge index, with some offset;
        /// 2nd index: node index;
        /// </param>
        /// <param name="jEdge">
        /// local edge index of first edge to process
        /// </param>
        /// <param name="EdgeTags"></param>
        /// <param name="EdgeTagsOffset"></param>
        /// <param name="normalFliped"></param>
        void BorderEdgeFlux(double time, int jEdge,
                            MultidimensionalArray x,
                            MultidimensionalArray normal, bool normalFliped, 
                            byte[] EdgeTags, int EdgeTagsOffset,
                            MultidimensionalArray[] Uin,
                            MultidimensionalArray[] UinMean,
                            int IndexOffset, int Lenght,
                            MultidimensionalArray Output);


        /// <summary>
        /// 
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="normal"></param>
        /// <param name="Uin"></param>
        /// <param name="Uout"></param>
        /// <param name="UinMean"></param>
        /// <param name="UoutMean"></param>
        /// <param name="Offset"></param>
        /// <param name="Lenght"></param>
        /// <param name="Output"></param>
        /// <param name="jEdge">
        /// local edge index of first edge to process
        /// </param>
        void InnerEdgeFlux(double time, int jEdge,
                           MultidimensionalArray x,
                           MultidimensionalArray normal, 
                           MultidimensionalArray[] Uin,
                           MultidimensionalArray[] Uout,
                           MultidimensionalArray[] UinMean,
                           MultidimensionalArray[] UoutMean,
                           int Offset, int Lenght,
                           MultidimensionalArray Output);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="U"></param>
        /// <param name="Offset"></param>
        /// <param name="Length"></param>
        /// <param name="Output"></param>
        /// <param name="jCell">
        /// local cell index of the first (vectorized evaluation !) cell
        /// </param>
        void Flux(double time,
                  MultidimensionalArray x,
                  MultidimensionalArray[] U,
                  int Offset, int Length,
                  MultidimensionalArray Output,
                  int jCell);

    }

    /// <summary>
    /// 
    /// </summary>
    public struct InParams {
        /// <summary>
        /// spatial position in space, at which flux should be evaluated
        /// </summary>
        public double[] X;

        /// <summary>
        /// normal vector, pointing from IN- to OUT-cell;<br/>
        /// index: spatial direction.
        /// </summary>
        public double[] normal;

        /// <summary>
        /// minimum length of any line within edge simplex
        /// </summary>
        public double h_min_edge;

        /// <summary>
        /// maximum length of any line within edge simplex
        /// </summary>
        public double h_max_edge;

        /// <summary>
        /// minimum length of any line within in IN - cell;
        /// </summary>
        public double h_min_in;

        /// <summary>
        /// maximum length of any line within in IN - cell;
        /// </summary>
        public double h_max_in;

        /// <summary>
        /// minimum length of any line within in OUT - cell;
        /// (invalid for border edges);
        /// </summary>
        public double h_min_out;

        /// <summary>
        /// maximum length of any line within in OUT - cell;
        /// (invalid for border edges);
        /// </summary>
        public double h_max_out;

        /// <summary>
        /// volume (to be more exact: D - dimensional measure) of the IN - cell
        /// </summary>
        public double VolumeIn;

        /// <summary>
        /// volume (to be more exact: D - dimensional measure) of the OUT - cell
        /// </summary>
        public double VolumeOut;


        /// <summary>
        /// area (to be more exact: (D-1) - dimensional measure) of boundary the IN - cell
        /// </summary>
        public double SurfaceAreaIn;

        /// <summary>
        /// area (to be more exact: (D-1) - dimensional measure) of boundary the OUT - cell
        /// </summary>
        public double SurfaceAreaOut;
        
        /// <summary>
        /// area (to be more exact: (D-1) - dimensional measure) of the edge
        /// </summary>
        public double EdgeArea;

        /// <summary>
        /// local cell index for IN - cell;
        /// </summary>
        public int jCellIn;

        /// <summary>
        /// local cell index for OUT - cell;
        /// (invalid for border edges);
        /// </summary>
        public int jCellOut;

        /// <summary>
        /// local edge index
        /// </summary>
        public int jEdge;

        /// <summary>
        /// edge tag
        /// </summary>
        public byte EdgeTag;

        /// <summary>
        /// values of parameters in IN - cell; values are sorted according to <see cref="IEquationComponent.ParameterOrdering"/>;
        /// </summary>
        public double[] ParameterValuesIn;

        /// <summary>
        /// values of parameters in OUT - cell; values are sorted according to <see cref="IEquationComponent.ParameterOrdering"/>;
        /// </summary>
        public double[] ParameterValuesOut;
                

        /// <summary>
        /// true if integrating over an edge that is on the boundary of the subgrid
        /// </summary>
        public bool SubGridBoundary;
    }

    /// <summary>
    /// These Flags Control, whether Certain Terms are evaluated during quadrature of the Forms
    /// Multiple Flags can be Defined by using "bitwise and" i.e. "|".
    /// E.g.
    /// TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV
    /// </summary>
    [Flags]
    public enum TermActivationFlags {
        /// <summary>
        /// Don't evaluate Anything
        /// </summary>
        None = 0,

        /// <summary>
        /// Trial- x TestFunction
        /// </summary>
        UxV = 0x1,

        /// <summary>
        /// Gradient of TrialFunction x TestFunction
        /// </summary>
        GradUxV = 0x2,

        /// <summary>
        /// TrialFunction Gradient of TestFunction
        /// </summary>
        UxGradV = 0x4,

        /// <summary>
        /// Gradient of TrialFunction x Gradient of TestFunction
        /// </summary>
        GradUxGradV = 0x8,

        /// <summary>
        /// TestFunction Only
        /// This is required e.g. for SourceTerms and BoundaryValues
        /// </summary>
        V = 0x10,

        /// <summary>
        /// Gradient of TestFunction Only
        /// This is required e.g. for SourceTerms and BoundaryValues
        /// </summary>
        GradV = 0x20,

        /// <summary>
        /// All flags are activated
        /// </summary>
        AllOn = 0xFF
    }

    /// <summary>
    /// parameter value structure for interior edges.
    /// </summary>
    public struct CommonParams {

        /// <summary>
        /// normal vector
        /// </summary>
        public Vector Normal;

        /// <summary>
        /// Quadrature node in global coordinates
        /// </summary>
        public Vector X;

        /// <summary>
        /// parameter values on IN-cell
        /// </summary>
        public double[] Parameters_IN;

        /// <summary>
        /// parameter values on OUT-cell
        /// </summary>
        public double[] Parameters_OUT;

        /// <summary>
        /// edge index (local on current MPI rank), i.e. 0th index into <see cref="ILogicalEdgeData.CellIndices"/>;
        /// - the IN-cell is: <see cref="ILogicalEdgeData.CellIndices"/>[iEdge, 0] == <see cref="jCellIn"/>
        /// - the OUT-cell is: <see cref="ILogicalEdgeData.CellIndices"/>[iEdge, 1] == <see cref="jCellOut"/>
        /// - access via <see cref="GridDat"/>, <see cref="IGridData.iLogicalEdges"/>, <see cref="ILogicalEdgeData.CellIndices"/>
        /// </summary>
        public int iEdge;

        /// <summary>
        /// Index of IN-cell (local geometrical index).
        /// </summary>
        public int jCellIn;

        /// <summary>
        /// Index of OUT-cell (local geometrical index).
        /// </summary>
        public int jCellOut;

        /// <summary>
        /// reference to grid data structure.
        /// </summary>
        public IGridData GridDat;

        /// <summary>
        /// spatial dimension
        /// </summary>
        public int D {
            get {
                Debug.Assert(X.Dim == Normal.Dim);
                Debug.Assert(X.Dim == GridDat.SpatialDimension);
                return X.Dim;
            }
        }

        /// <summary>
        /// Physical time.
        /// </summary>
        public double time;

        /// <summary>
        /// edge tag for the respective edge (see <see cref="BoSSS.Foundation.Grid.IGeometricalEdgeData.EdgeTags"/>)
        /// For interior edges, this is 
        /// - equal to 0 on ordinary edges
        /// - greater or equal to <see cref="Grid.Classic.GridCommons.FIRST_PERIODIC_BC_TAG"/> on periodic edges
        /// </summary>
        public byte EdgeTag;
    }

    /// <summary>
    /// parameter value structure for boundary edges.
    /// </summary>
    public struct CommonParamsBnd {

        /// <summary>
        /// normal vector
        /// </summary>
        public Vector Normal;

        /// <summary>
        /// Quadrature node in global coordinates
        /// </summary>
        public Vector X;

        /// <summary>
        /// parameter values on IN-cell
        /// </summary>
        public double[] Parameters_IN;

        /// <summary>
        /// edge tag for the respective edge (see <see cref="BoSSS.Foundation.Grid.IGeometricalEdgeData.EdgeTags"/>).
        /// </summary>
        public byte EdgeTag;

        /// <summary>
        /// edge index (local on current MPI rank)
        /// </summary>
        public int iEdge;

        /// <summary>
        /// Index of IN-cell (local geometrical index).
        /// </summary>
        public int jCellIn {
            get {
                return GridDat.iGeomEdges.CellIndices[iEdge, 0];
            }
        }

        /// <summary>
        /// Physical time.
        /// </summary>
        public double time;

        /// <summary>
        /// reference to grid data structure.
        /// </summary>
        public IGridData GridDat;

        /// <summary>
        /// spatial dimension
        /// </summary>
        public int D {
            get {
                Debug.Assert(X.Dim == Normal.Dim);
                Debug.Assert(X.Dim == GridDat.SpatialDimension);
                return X.Dim;
            }
        }
    }


    /// <summary>
    /// Defines a general _edge term_, i.e. a form
    /// ```math 
    /// a(\vec{U},v) = \int_{\partial K} f(\vec{U}) g(v) n  \mathrm{dS} .
    /// ```
    /// </summary>
    public interface IInnerEdgeForm : IEquationComponent {

        /// <summary>
        /// The form which is integrated over interior edges
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="_uIN">value of trial function on IN-cell</param>
        /// <param name="_Grad_uIN">gradient of trial function on IN-cell</param>
        /// <param name="_vIN">value of test function on IN-cell</param>
        /// <param name="_Grad_vIN">gradient of test function on IN-cell</param>
        /// <param name="_uOUT">value of trial function on OUT-cell</param>
        /// <param name="_Grad_uOUT">gradient of trial function on OUT-cell</param>
        /// <param name="_vOUT">value of test function on OUT-cell</param>
        /// <param name="_Grad_vOUT">gradient of test function on OUT-cell</param>
        /// <returns></returns>
        double InnerEdgeForm(ref CommonParams inp,
            double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT,
            double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT);
    }

    /// <summary>
    /// Defines a general _edge term_, i.e. a form
    /// ```math 
    /// a(\vec{U},v) = \int_{\partial K} f(\vec{U}) g(v) n  \mathrm{dS} .
    /// ```
    /// </summary>
    public interface IBoundaryEdgeForm : IEquationComponent {
        /// <summary>
        /// The form which is integrated over boundary edges.
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="_uA">value of trial function on IN-cell</param>
        /// <param name="_Grad_uA">gradient of trial function on IN-cell</param>
        /// <param name="_vA">value of test function on IN-cell</param>
        /// <param name="_Grad_vA">gradient of test function on IN-cell</param>
        /// <returns></returns>
        double BoundaryEdgeForm(ref CommonParamsBnd inp,
            double[] _uA, double[,] _Grad_uA,
            double _vA, double[] _Grad_vA);
    }

    /// <summary>
    /// Defines a general complete _edge term_ consisting of boundary and internal component, i.e. a form
    /// ```math 
    /// a(\vec{U},v) = \int_{\partial K} f(\vec{U}) g(v) n  \mathrm{dS} .
    /// ```
    /// </summary>
    public interface IEdgeForm : IInnerEdgeForm, IBoundaryEdgeForm {

        /// <summary>
        /// Activation Flags For Boundary Edges
        /// <see cref="TermActivationFlags"/>
        /// </summary>
        TermActivationFlags BoundaryEdgeTerms { get; }

        /// <summary>
        /// Activation Flags for Inner Edges 
        /// <see cref="TermActivationFlags"/>
        /// </summary>
        TermActivationFlags InnerEdgeTerms { get; }
    }



    /// <summary>
    /// parameter value structure for interior edges.
    /// </summary>
    public struct CommonParamsVol {
        /// <summary>
        /// cell index.
        /// </summary>
        public int jCell;

        /// <summary>
        /// reference to grid data structure.
        /// </summary>
        public IGridData GridDat;

        /// <summary>
        /// Physical time.
        /// </summary>
        public double time;

        /// <summary>
        /// value of parameter fields.
        /// </summary>
        public double[] Parameters;

        /// <summary>
        /// node in global coordinates
        /// </summary>
        public Vector Xglobal;

        /// <summary>
        /// Spatial dimension.
        /// </summary>
        public int D {
            get {
                Debug.Assert(Xglobal.Dim == GridDat.SpatialDimension);
                return Xglobal.Dim;
            }
        }
    }

    /// <summary>
    /// defines a volume term.
    /// ```math 
    /// a(\vec{U},v) = \int_{\Omega} f(\vec{U}) g(v)   \mathrm{dX}
    /// ```
    /// </summary>
    public interface IVolumeForm : IEquationComponent {

        /// <summary>
        /// <see cref="TermActivationFlags"/>
        /// </summary>
        TermActivationFlags VolTerms { get; }


        /// <summary>
        /// bilinear form on volume integrals.
        /// </summary>
        /// <param name="cpv"></param>
        /// <param name="GradU">gradient of trial function</param>
        /// <param name="GradV">gradient of test function</param>
        /// <param name="U">value of trial function</param>
        /// <param name="V">value of test function</param>
        /// <returns></returns>
        double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV);
    }


    /// <summary>
    /// parameters on which volume forms 
    /// (e.g. <see cref="IVolumeForm_UxV"/>, <see cref="IVolumeForm_GradUxV"/>, <see cref="IVolumeForm_UxGradV"/>, <see cref="IVolumeForm_GradUxGradV"/>)
    /// may depend.
    /// </summary>
    public struct VolumFormParams {
        
        /// <summary>
        /// first cell
        /// </summary>
        public int j0;

        /// <summary>
        /// number of cells
        /// </summary>
        public int Len;

        /// <summary>
        /// reference to the current grid 
        /// </summary>
        public IGridData GridDat;

        /// <summary>
        /// Physical time.
        /// </summary>
        public double time;

        /// <summary>
        /// Values of parameter fields at quadrature nodes<br/>
        /// array index: parameter variable, as specified by the parameter mapping (see <see cref="IEquationComponent.ParameterOrdering"/>); <br/>
        /// for each multidimensional array: 
        ///  - 1st index: cell 
        ///  - 2nd index: node 
        /// </summary>
        public MultidimensionalArray[] ParameterVars;

        /// <summary>
        /// quadrature nodes in global coordinates
        ///  - 1st index: cell 
        ///  - 2nd index: quad node 
        ///  - 3rd index: spatial direction 
        /// </summary>
        public MultidimensionalArray Xglobal;
    }

    /// <summary>
    /// a bi-linear form of the type
    /// ```math 
    ///    a(U,v) = \sum_{l} \int_{\Gamma_{\mathrm{int}}} 
    ///              u_l  {f}_{l}(\vec{x}) v
    ///           \ \mathrm{dV}
    /// ```.
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable) and 
    /// \f$ U = (u_0, \ldots, u_{L-1} ) \f$ denotes the trial functions (correspond to domain variable, defined by the 
    /// argument ordering <see cref="IEquationComponent.ArgumentOrdering"/>).
    /// </summary>
    public interface IVolumeForm_UxV : IVolumeForm {

        /// <summary>
        /// the values of \f$ {f}_{l}(\vec{x}) \f$
        /// </summary>
        /// <param name="prm">parameters on which the form may depend</param>
        /// <param name="UxV">
        /// Output: the values of \f$ f_l(\vec{x}) \f$ 
        /// - 1st index: cell index 
        /// - 2nd index: quadrature node 
        /// - 3rd index: correlates with argument ordering, i.e. index \f$ l \f$, of trial function; see <see cref="IEquationComponent.ArgumentOrdering"/> 
        /// </param>
        void Form(ref VolumFormParams prm, MultidimensionalArray UxV);

    }

    /// <summary>
    /// a bi-linear form of the type
    /// ```math 
    ///    a(U,v) = \sum_{l} \int_{\Gamma_{\mathrm{int}}} 
    ///               v \vec{f}_{l}(\vec{x}) \cdot \nabla u_l
    ///           \ \mathrm{dV}
    /// ```.
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable) and 
    /// \f$ U = (u_0, \ldots, u_{L-1} ) \f$ denotes the trial functions (correspond to domain variable, defined by the 
    /// argument ordering <see cref="IEquationComponent.ArgumentOrdering"/>).
    /// </summary>
    public interface IVolumeForm_GradUxV : IVolumeForm {

        /// <summary>
        /// the values of \f$ \vec{f}_{l}(\vec{x}) \f$
        /// </summary>
        /// <param name="prm">parameters on which the form may depend</param>
        /// <param name="GradUxV">
        /// Output: the values of \f$ \vec{f}_{l}(\vec{x}) \f$ 
        ///  - 1st index: cell index 
        ///  - 2nd index: quadrature node 
        ///  - 3rd index: correlates with argument ordering, i.e. index \f$ l \f$, of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> 
        ///  - 4th index: component index of \f$ \vec{f}_{l}(\vec{x}) \f$, i.e. correlates with 
        ///            component index of trial function gradient \f$ \nabla u_l \f$
        /// </param>
        void Form(ref VolumFormParams prm, MultidimensionalArray GradUxV);

    }

    /// <summary>
    /// a bi-linear form of the type
    /// ```math 
    ///    a(U,v) = \sum_{l} \int_{\Gamma_{\mathrm{int}}} 
    ///               \nabla v \cdot \vec{f}_{l}(\vec{x}) u_l
    ///           \ \mathrm{dV}
    /// ```.
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable) and 
    /// \f$ U = (u_0, \ldots, u_{L-1} ) \f$ denotes the trial functions (correspond to domain variable, defined by the 
    /// argument ordering <see cref="IEquationComponent.ArgumentOrdering"/>).
    /// </summary>
    public interface IVolumeForm_UxGradV : IVolumeForm {

        /// <summary>
        /// the values of \f$ \vec{f}_{l}(\vec{x}) \f$
        /// </summary>
        /// <param name="prm">parameters on which the form may depend</param>
        /// <param name="UxGradV">
        ///  - Output: the values of \f$ \vec{f}_{l}(\vec{x}) \f$ 
        ///  - 1st index: cell index 
        ///  - 2nd index: quadrature node 
        ///  - 3rd index: correlates with argument ordering, i.e. index \f$ l \f$, of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> 
        ///  - 4th index: component index of \f$ \vec{f}_{l}(\vec{x}) \f$, i.e. correlates with 
        ///            component index of test function gradient \f$ \nabla v \f$
        /// </param>
        void Form(ref VolumFormParams prm, MultidimensionalArray UxGradV);

    }

    /// <summary>
    /// a bi-linear form of the type
    /// ```math 
    ///    a(U,v) = \sum_{l} \int_{\Gamma_{\mathrm{int}}} 
    ///              \nabla v^T \cdot \vec{\vec{f}}_{l}(\vec{x}) \cdot \nabla u_l
    ///           \ \mathrm{dV}
    /// ```.
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable) and 
    /// \f$ U = (u_0, \ldots, u_{L-1} ) \f$ denotes the trial functions (correspond to domain variable, defined by the 
    /// argument ordering <see cref="IEquationComponent.ArgumentOrdering"/>).
    /// </summary>
    public interface IVolumeForm_GradUxGradV : IVolumeForm {

        /// <summary>
        /// the values of \f$ \vec{\vec{f}}_{l}(\vec{x}) \f$
        /// </summary>
        /// <param name="prm">parameters on which the form may depend</param>
        /// <param name="GradUxGradV">
        /// Output: the values of \f$ \vec{\vec{f}}_{l}(\vec{x}) \f$ <br/>
        /// 1st index: cell index <br/>
        /// 2nd index: quadrature node <br/>
        /// 3rd index: correlates with argument ordering, i.e. index \f$ l \f$, of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> <br/>
        /// 4th index: row index of \f$ \vec{\vec{f}}_{l}(\vec{x}) \f$, i.e. correlates with 
        ///            component index of test function gradient \f$ \nabla v \f$
        /// 5th index: column index of \f$ \vec{\vec{f}}_{l}(\vec{x}) \f$, i.e. correlates with 
        ///            component index of trial function gradient \f$ \nabla u_l \f$
        /// </param>
        void Form(ref VolumFormParams prm, MultidimensionalArray GradUxGradV);

    }

    /// <summary>
    /// a linear form of the type
    /// ```math 
    ///    a(v) = \int_{\Gamma_{\mathrm{int}}} 
    ///              v \ f(\vec{x})
    ///           \ \mathrm{dV}
    /// ```
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable).
    /// </summary>
    public interface IVolumeSource_V : IVolumeForm {

        /// <summary>
        /// the values of \f$ f(\vec{x}) \f$
        /// </summary>
        /// <param name="prm">parameters on which the form may depend</param>
        /// <param name="V">
        /// Output: the values of \f$ f(\vec{x}) \f$ <br/>
        /// 1st index: cell index <br/>
        /// 2nd index: quadrature node <br/>
        /// </param>
        void Form(ref VolumFormParams prm, MultidimensionalArray V);
    }


    /// <summary>
    /// a linear form of the type
    /// ```math 
    ///    a(v) = \int_{\Gamma_{\mathrm{int}}} 
    ///              \nabla v^T \cdot \vec{f}(\vec{x})
    ///           \ \mathrm{dV}
    /// ```
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable).
    /// </summary>
    public interface IVolumeSource_GradV : IVolumeForm {

        /// <summary>
        /// the values of \f$ \vec{f}(\vec{x}) \f$
        /// </summary>
        /// <param name="prm">parameters on which the form may depend</param>
        /// <param name="GradV">
        /// Output: the values of \f$ \vec{f}(\vec{x}) \f$ <br/>
        /// 1st index: cell index <br/>
        /// 2nd index: quadrature node <br/>
        /// 3rd index: component index/spatial direction
        ///            of \f$ \vec{f}(\vec{x}) \f$, i.e. correlates with 
        ///            component index of test function gradient \f$ \nabla v \f$
        /// </param>
        void Form(ref VolumFormParams prm, MultidimensionalArray GradV);
    }


    /// <summary>
    /// a non-linear form of the type
    /// ```math 
    ///    a(U,v) = \int_{\Gamma_{\mathrm{int}}} 
    ///              u_l  {f}(\vec{x},U) v
    ///             \ \mathrm{dV}
    /// ```
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable) and 
    /// \f$ U = (u_0, \ldots, u_{L-1} ) \f$ denotes the trial functions (correspond to domain variable, defined by the 
    /// argument ordering <see cref="IEquationComponent.ArgumentOrdering"/>).
    /// </summary>
    public interface INonlinVolumeForm_V : IVolumeForm {

        /// <summary>
        /// the function \f$ {f}(\vec{x},U,\nabla U) \f$
        /// </summary>
        /// <param name="prm">parameters on which the form may depend</param>
        /// <param name="f">
        /// Output: the values of \f$ {f}(\vec{x},U,\nabla U) \f$<br/>
        /// 1st index: cell index <br/>
        /// 2nd index: quadrature node <br/>
        /// </param>
        /// <param name="U">
        /// Input: the values of \f$ U \f$<br/>
        /// 1st index: correlates with argument ordering, i.e. index \f$ l \f$, of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> <br/>
        /// 2nd index: cell index <br/>
        /// 3rd index: quadrature node <br/>
        /// </param>
        /// <param name="GradU">
        /// Input: the values of \f$ \nabla U \f$ <br/>
        /// 1st index: correlates with argument ordering, i.e. index \f$ l \f$ , of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> <br/>
        /// 2nd index: cell index <br/>
        /// 3rd index: quadrature node <br/>
        /// 4th index: spatial direction of derivative
        /// </param>
        void Form(ref VolumFormParams prm, MultidimensionalArray[] U, MultidimensionalArray[] GradU, MultidimensionalArray f);
    }

    /// <summary>
    /// a non-linear form of the type
    /// ```math 
    ///    a(U,v) = \int_{\Gamma_{\mathrm{int}}} 
    ///              \nabla v^T \cdot \vec{f}(\vec{x},U,\nabla U)
    ///             \ \mathrm{dV}
    /// ```
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable) and 
    /// \f$ U = (u_0, \ldots, u_{L-1} ) \f$  denotes the trial functions (correspond to domain variable, defined by the 
    /// argument ordering <see cref="IEquationComponent.ArgumentOrdering"/>).
    /// </summary>
    public interface INonlinVolumeForm_GradV : IVolumeForm {

        /// <summary>
        /// the values of \f$ \vec{f}(\vec{x},U,\nabla U) \f$ 
        /// </summary>
        /// <param name="prm">parameters on which the form may depend</param>
        /// <param name="U">
        /// Input: the values of \f$ U \f$ <br/>
        /// 1st index: correlates with argument ordering, i.e. index \f$ l \f$ , of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> <br/>
        /// 2nd index: cell index <br/>
        /// 3rd index: quadrature node <br/>
        /// </param>
        /// <param name="GradU">
        /// Input: the values of \f$ \nabla U \f$ <br/>
        /// 1st index: correlates with argument ordering, i.e. index \f$ l \f$ , of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> <br/>
        /// 2nd index: cell index <br/>
        /// 3rd index: quadrature node <br/>
        /// 4th index: spatial direction of derivative
        /// </param>
        /// <param name="f">
        /// Output: the values of \f$ \vec{f}(\vec{x},U,\nabla U) \f$  <br/>
        /// 1st index: cell index <br/>
        /// 2nd index: quadrature node <br/>
        /// 3rd index: vector-component of \f$ \vec{f} \f$  <br/>
        /// </param>
        void Form(ref VolumFormParams prm, MultidimensionalArray[] U, MultidimensionalArray[] GradU, MultidimensionalArray f);
    }



    /// <summary>
    /// parameters on which edge forms (e.g. <see cref="IEdgeForm_UxV"/>, <see cref="IEdgeform_GradUxV"/>, <see cref="IEdgeform_UxGradV"/>, <see cref="IEdgeSource_V"/>, <see cref="IEdgeSource_GradV"/>)
    /// may depend.
    /// </summary>
    public struct EdgeFormParams {
        
        /// <summary>
        /// Quadrature nodes in global coordinates. 
        /// - 1st index: edge
        /// - 2nd index: quadrature node
        /// - 3rd index: spatial direction
        /// </summary>
        public MultidimensionalArray Nodes;

        /// <summary>
        /// Edge normals at quadrature nodes. 
        /// - 1st index: edge
        /// - 2nd index: quadrature node
        /// - 3rd index: spatial direction
        /// </summary>
        public MultidimensionalArray Normals;

        /// <summary>
        /// Values of parameter fields at quadrature nodes, for the IN-cell:
        /// array index: parameter variable, as specified by the parameter mapping (see <see cref="IEquationComponent.ParameterOrdering"/>); 
        /// for each multidimensional array: 
        /// - 1st index: cell
        /// - 2nd index: node 
        /// </summary>
        public MultidimensionalArray[] ParameterVars_IN;

        /// <summary>
        /// Values of parameter fields at quadrature nodes, for the OUT-cell:
        /// array index: parameter variable, as specified by the parameter mapping (see <see cref="IEquationComponent.ParameterOrdering"/>); 
        /// for each multidimensional array:
        /// - 1st index: cell
        /// - 2nd index: node 
        /// </summary>
        public MultidimensionalArray[] ParameterVars_OUT;

       
        /// <summary>
        /// first edge
        /// </summary>
        public int e0;

        /// <summary>
        /// number of edges
        /// </summary>
        public int Len;

        /// <summary>
        /// Physical time.
        /// </summary>
        public double time;

        /// <summary>
        /// reference to the current grid 
        /// </summary>
        public IGridData GridDat;
    }

    /// <summary>
    /// Inner edge integrand of <see cref="IEdgeForm_UxV"/>
    /// </summary>
    public interface IInnerEdgeform_UxV : IInnerEdgeForm {
        
        /// <summary>
        /// the values of \f$ {f}_{i \ j \ l}(\vec{x})\f$   on interior edges on \f$ \Gamma_{\mathrm{int}} \f$ .
        /// </summary>
        /// <param name="efp"></param>
        /// <param name="UxV">
        /// output: the values of \f$ {f}_{i \ j \ l}(\vec{x}) \f$:
        /// - 1st index: edge index
        /// - 2nd index: quadrature node
        /// - 3rd index: in and out - coefficients with respect to trial function, i.e. index \f$ i \f$ 
        ///              ('U': index 0 corresponds to IN-cell, index 1 corresponds to OUT-cell). 
        /// - 4th index: in and out - coefficients with respect to test function, i.e. index \f$ j \f$ 
        ///              ('V': index 0 corresponds to IN-cell, index 1 corresponds to OUT-cell). 
        /// - 5th index: correlates with argument ordering, i.e. index \f$ l \f$ , of trial function; 
        ///              see <see cref="IEquationComponent.ArgumentOrdering"/>
        /// </param>
        void InternalEdge_UxV(ref EdgeFormParams efp, MultidimensionalArray UxV);
    }

    /// <summary>
    /// Boundary edge integrand of <see cref="IEdgeForm_UxV"/>
    /// </summary>
    public interface IBoundaryEdgeform_UxV : IBoundaryEdgeForm {
        
        /// <summary>
        /// the values of \f$ {f}_{0 \ 0 \ l}(\vec{x}) \f$   on boundary edges on \f$ \partial \Omega \f$ .
        /// </summary>
        /// <param name="efp"></param>
        /// <param name="UxV">
        /// output: the values of \f$ {f}_{0 \ 0 \ l}(\vec{x}) \f$ :<br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// 3rd index: correlates with argument ordering (of trial function; see <see cref="IEquationComponent.ArgumentOrdering"/> <br/>
        /// </param>
        void BoundaryEdge_UxV(ref EdgeFormParams efp, MultidimensionalArray UxV);
    }


    /// <summary>
    /// a bi-linear form of the type
    /// ```math 
    ///    a(U,v) = \sum_{l} \oint_{\Gamma_{\mathrm{int}}} 
    ///              u^\mathrm{in}_l  {f}_{0 \ 0 \ l}(\vec{x}) v^\mathrm{in} 
    ///            + u^\mathrm{in}_l  {f}_{0 \ 1 \ l}(\vec{x}) v^\mathrm{out} 
    ///            + u^\mathrm{out}_l {f}_{1 \ 0 \ l}(\vec{x}) v^\mathrm{in} 
    ///            + u^\mathrm{out}_l {f}_{1 \ 1 \ l}(\vec{x}) v^\mathrm{out} 
    ///           \ \mathrm{dS}
    ///         +
    ///          \sum_{l} \oint_{\partial \Omega} 
    ///              u^\mathrm{in}_\gamma {f}_{0 \ 0 \ l}(\vec{x}) v^\mathrm{in} 
    ///           \ \mathrm{dS}          
    /// ```
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable) and 
    /// \f$ U = (u_0, \ldots, u_{L-1} ) \f$  denotes the trial functions (correspond to domain variable, defined by the 
    /// argument ordering <see cref="IEquationComponent.ArgumentOrdering"/>).
    /// </summary>
    public interface IEdgeForm_UxV : IInnerEdgeform_UxV, IBoundaryEdgeform_UxV, IEdgeForm {


       
    }


    /// <summary>
    /// a bi-linear form of the type
    /// ```math 
    ///    a(U,v) = \sum_{l} \oint_{\Gamma_{\mathrm{int}}} 
    ///              \nabla u^\mathrm{in}_l  \cdot \vec{f}_{0 \ 0 \ l}(\vec{x}) v^\mathrm{in} 
    ///            + \nabla u^\mathrm{in}_l  \cdot \vec{f}_{0 \ 1 \ l}(\vec{x}) v^\mathrm{out} 
    ///            + \nabla u^\mathrm{out}_l \cdot \vec{f}_{1 \ 0 \ l}(\vec{x}) v^\mathrm{in} 
    ///            + \nabla u^\mathrm{out}_l \cdot \vec{f}_{1 \ 1 \ l}(\vec{x}) v^\mathrm{out} 
    ///           \ \mathrm{dS}
    ///         +
    ///          \sum_{l} \oint_{\partial \Omega} 
    ///              \nabla u^\mathrm{in}_\gamma \cdot \vec{f}_{0 \ 0 \ l}(\vec{x}) v^\mathrm{in} 
    ///           \ \mathrm{dS}          
    /// ```
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable) and 
    /// \f$ U = (u_0, \ldots, u_{L-1} ) \f$  denotes the trial functions (correspond to domain variable, defined by the 
    /// argument ordering <see cref="IEquationComponent.ArgumentOrdering"/>).
    /// </summary>
    public interface IEdgeform_GradUxV : IEdgeForm, IInnerEdgeform_GradUxV, IBoundaryEdgeform_GradUxV {

    }

    /// <summary>
    /// Inner part of <see cref="IEdgeform_GradUxV"/>
    /// </summary>
    public interface IInnerEdgeform_GradUxV : IInnerEdgeForm {
        /// <summary>
        /// the values of \f$ \vec{f}_{i \ j \ l}(\vec{x}) \f$   on interior edges on \f$ \Gamma_{\mathrm{int}} \f$ .
        /// </summary>
        /// <param name="efp"></param>
        /// <param name="GradUxV">
        /// output: the values of \f$ \vec{f}_{i \ j \ l}(\vec{x}) \f$ :<br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// 3rd index: in and out - coefficients with respect to trial function, i.e. index \f$ i \f$ 
        ///            ('U': index 0 corresponds to IN-cell, index 1 corresponds to OUT-cell). <br/>
        /// 4th index: in and out - coefficients with respect to test function, i.e. index \f$ j \f$ 
        ///            ('V': index 0 corresponds to IN-cell, index 1 corresponds to OUT-cell). <br/>
        /// 5th index: correlates with argument ordering, i.e. index \f$ l \f$ , of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> <br/>
        /// 6th index: spatial direction of trial ('U') function gradient  <br/>
        /// </param>
        void InternalEdge_GradUxV(ref EdgeFormParams efp, MultidimensionalArray GradUxV);
    }

    /// <summary>
    /// Boundary part of <see cref="IEdgeform_GradUxV"/>
    /// </summary>
    public interface IBoundaryEdgeform_GradUxV : IBoundaryEdgeForm { 
        
        /// <summary>
        /// the values of \f$ \vec{f}_{0 \ 0 \ l}(\vec{x}) \f$   on boundary edges on \f$ \partial \Omega\f$ .
        /// </summary>
        /// <param name="efp"></param>
        /// <param name="GradUxV">
        /// output: the values of \f$ \vec{f}_{0 \ 0 \ l}(\vec{x}) \f$ :<br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// 3rd index: correlates with argument ordering, i.e. index \f$ l \f$ , of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> <br/>
        /// 4th index: spatial direction of trial function gradient <br/>
        /// </param>
        void BoundaryEdge_GradUxV(ref EdgeFormParams efp, MultidimensionalArray GradUxV);
    }

    /// <summary>
    /// a bi-linear form of the type
    /// ```math 
    ///    a(U,v) = \sum_{l} \oint_{\Gamma_{\mathrm{int}}} 
    ///              u^\mathrm{in}_l  \vec{f}_{0 \ 0 \ l}(\vec{x}) \cdot \nabla v^\mathrm{in} 
    ///            + u^\mathrm{in}_l  \vec{f}_{0 \ 1 \ l}(\vec{x}) \cdot \nabla v^\mathrm{out} 
    ///            + u^\mathrm{out}_l \vec{f}_{1 \ 0 \ l}(\vec{x}) \cdot \nabla v^\mathrm{in} 
    ///            + u^\mathrm{out}_l \vec{f}_{1 \ 1 \ l}(\vec{x}) \cdot \nabla v^\mathrm{out} 
    ///           \ \mathrm{dS}
    ///         +
    ///          \sum_{l} \oint_{\partial \Omega} 
    ///              u^\mathrm{in}_\gamma \vec{f}_{0 \ 0 \ l}(\vec{x}) \cdot \nabla v^\mathrm{in} 
    ///           \ \mathrm{dS}          
    /// ```
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable) and 
    /// \f$ U = (u_0, \ldots, u_{L-1} ) \f$  denotes the trial functions (correspond to domain variable, defined by the 
    /// argument ordering <see cref="IEquationComponent.ArgumentOrdering"/>).
    /// </summary>
    public interface IEdgeform_UxGradV : IEdgeForm, IInnerEdgeform_UxGradV, IBoundaryEdgeform_UxGradV {

    }

    /// <summary>
    /// <see cref="IInnerEdgeform_UxGradV"/>
    /// </summary>
    public interface IInnerEdgeform_UxGradV : IInnerEdgeForm { 
        /// <summary>
        /// the values of \f$ \vec{f}_{i \ j \ l}(\vec{x}) \f$   on interior edges on \f$ \Gamma_{\mathrm{int}} \f$ .
        /// </summary>
        /// <param name="efp"></param>
        /// <param name="UxGradV">
        /// Output: the value of \f$ \vec{f}_{i \ j \ l}(\vec{x}) \f$  <br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// 3rd index: in and out - coefficients with respect to trial function, i.e. index \f$ i \f$ 
        ///            ('U': index 0 corresponds to IN-cell, index 1 corresponds to OUT-cell). <br/>
        /// 4th index: in and out - coefficients with respect to test function, i.e. index \f$ j \f$ 
        ///            ('V': index 0 corresponds to IN-cell, index 1 corresponds to OUT-cell). <br/>
        /// 5th index: correlates with argument ordering, i.e. index \f$ l \f$ , of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> <br/>
        /// 6th index: spatial direction of test function ('V') gradient  <br/>
        /// </param>
        void InternalEdge_UxGradV(ref EdgeFormParams efp, MultidimensionalArray UxGradV);
    }

    /// <summary>
    /// <see cref="IEdgeform_UxGradV"/>
    /// </summary>
    public interface IBoundaryEdgeform_UxGradV : IBoundaryEdgeForm { 
        /// <summary>
        /// the values of \f$ \vec{f}_{0 \ 0 \ l}(\vec{x}) \f$   on boundary edges on \f$ \partial \Omega\f$ .
        /// </summary>
        /// <param name="efp"></param>
        /// <param name="UxGradV">
        /// Output: values of \f$ \vec{f}_{0 \ 0 \ l}(\vec{x}) \f$ .<br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// 3rd index: correlates with argument ordering, i.e. index \f$ l \f$ , of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> <br/>
        /// 4th index: spatial direction of test function ('V') gradient <br/>
        /// </param>
        void BoundaryEdge_UxGradV(ref EdgeFormParams efp, MultidimensionalArray UxGradV);
    }

    /// <summary>
    /// a bi-linear form of the type
    /// ```math 
    ///    a(U,v) = \sum_{l} \oint_{\Gamma_{\mathrm{int}}} 
    ///              (u^\mathrm{in}_l)^T   \cdot \vec{\vec{f}}_{0 \ 0 \ l}(\vec{x}) \cdot \nabla v^\mathrm{in} 
    ///            + (u^\mathrm{in}_l)^T   \cdot \vec{\vec{f}}_{0 \ 1 \ l}(\vec{x}) \cdot \nabla v^\mathrm{out} 
    ///            + (u^\mathrm{out}_l)^T  \cdot \vec{\vec{f}}_{1 \ 0 \ l}(\vec{x}) \cdot \nabla v^\mathrm{in} 
    ///            + (u^\mathrm{out}_l)^T  \cdot \vec{\vec{f}}_{1 \ 1 \ l}(\vec{x}) \cdot \nabla v^\mathrm{out} 
    ///           \ \mathrm{dS}
    ///         +
    ///          \sum_{l} \oint_{\partial \Omega} 
    ///              u^\mathrm{in}_\gamma \vec{f}_{0 \ 0 \ l}(\vec{x}) \cdot \nabla v^\mathrm{in} 
    ///           \ \mathrm{dS}          
    /// ```
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable) and 
    /// \f$ U = (u_0, \ldots, u_{L-1} ) \f$  denotes the trial functions (correspond to domain variable, defined by the 
    /// argument ordering <see cref="IEquationComponent.ArgumentOrdering"/>).
    /// </summary>
    public interface IEdgeform_GradUxGradV : IEdgeForm, IInnerEdgeform_GradUxGradV, IBoundaryEdgeform_GradUxGradV {

    }

    /// <summary>
    /// Inner part of <see cref="IEdgeform_GradUxGradV"/>
    /// </summary>
    public interface IInnerEdgeform_GradUxGradV : IInnerEdgeForm {
        /// <summary>
        /// the values of \f$ \vec{f}_{i \ j \ l}(\vec{x}) \f$   on interior edges on \f$ \Gamma_{\mathrm{int}} \f$ .
        /// </summary>
        /// <param name="efp"></param>
        /// <param name="GradUxGradV">
        /// Output: the value of \f$ \vec{f}_{i \ j \ l}(\vec{x}) \f$  <br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// 3rd index: in and out - coefficients with respect to trial function, i.e. index \f$ i \f$ 
        ///            ('U': index 0 corresponds to IN-cell, index 1 corresponds to OUT-cell). <br/>
        /// 4th index: in and out - coefficients with respect to test function, i.e. index \f$ j \f$ 
        ///            ('V': index 0 corresponds to IN-cell, index 1 corresponds to OUT-cell). <br/>
        /// 5th index: correlates with argument ordering, i.e. index \f$ l\f$ , of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> <br/>
        /// 6th index: spatial direction of test function ('U') gradient  <br/>
        /// 7th index: spatial direction of test function ('V') gradient  <br/>
        /// </param>
        void InternalEdge_GradUxGradV(ref EdgeFormParams efp, MultidimensionalArray GradUxGradV);
    }

    /// <summary>
    /// Boundary part of <see cref="IEdgeform_GradUxGradV"/>
    /// </summary>
    public interface IBoundaryEdgeform_GradUxGradV : IBoundaryEdgeForm {
        /// <summary>
        /// the values of \f$ \vec{f}_{0 \ 0 \ l}(\vec{x})\f$   on boundary edges on \f$ \partial \Omega\f$ .
        /// </summary>
        /// <param name="efp"></param>
        /// <param name="GradUxGradV">
        /// Output: values of \f$ \vec{f}_{0 \ 0 \ l}(\vec{x})\f$ .<br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// 3rd index: correlates with argument ordering, i.e. index \f$ l\f$ , of trial function; 
        ///            see <see cref="IEquationComponent.ArgumentOrdering"/> <br/>
        /// 4th index: spatial direction of test function ('U') gradient <br/>
        /// 5th index: spatial direction of test function ('V') gradient <br/>
        /// </param>
        void BoundaryEdge_GradUxGradV(ref EdgeFormParams efp, MultidimensionalArray GradUxGradV);
    }

    /// <summary>
    /// a linear form of the type
    /// ```math 
    ///    a(v) = 
    ///         \oint_{\Gamma_{\mathrm{int}}} f_{0}(\vec{x}) v^\mathrm{in} + f_{1}(\vec{x}) v^\mathrm{out}  \ \mathrm{dS}
    ///         + \oint_{\partial \Omega} f_{0}(\vec{x}) v^\mathrm{in}  \ \mathrm{dS}
    /// ```
    /// where <em>v</em> denotes the test function.
    /// </summary>
    public interface IEdgeSource_V : IEdgeForm, IInnerEdgeSource_V, IBoundaryEdgeSource_V {
    }

    /// <summary>
    /// <see cref="IEdgeSource_V"/>
    /// </summary>
    public interface IInnerEdgeSource_V : IInnerEdgeForm {
        /// <summary>
        /// the point-wise source term \f$ f_{i}(\vec{x})\f$  on interior edges \f$ \Gamma_{\mathrm{int}}\f$ .
        /// </summary>
        /// <param name="efp">parameters on which  \f$ f\f$  may depend on.</param>
        /// <param name="V">
        /// Output: the value of \f$ f(\vec{x})\f$ , which will be multiplied with the test function gradient. <br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// 3rd index: in and out - coefficients with respect to test function, i.e. index \f$ i\f$ 
        ///            ('V': index 0 corresponds to IN-cell, index 1 corresponds to OUT-cell). <br/>
        /// </param>
        void InternalEdge_V(ref EdgeFormParams efp, MultidimensionalArray V);
    }

    /// <summary>
    /// <see cref="IEdgeSource_V"/>
    /// </summary>
    public interface IBoundaryEdgeSource_V : IBoundaryEdgeForm {
        /// <summary>
        /// the point-wise source term \f$ f_{0}(\vec{x}) \f$  on boundary edges \f$ \partial \Omega\f$ .
        /// </summary>
        /// <param name="efp">parameters on which  \f$ f\f$  may depend on.</param>
        /// <param name="V">
        /// Output: the value of \f$ f(\vec{x})\f$ , which will be multiplied with the test function gradient. <br/>
        /// The value of the test function.<br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// </param>
        void BoundaryEdge_V(ref EdgeFormParams efp, MultidimensionalArray V);
    }

    /// <summary>
    /// a linear form of the type
    /// ```math 
    ///    a(v) = \oint_{\Gamma_{\mathrm{int}}} 
    ///                     \vec{f}_{0}(\vec{x}) \cdot \nabla v^\mathrm{in} 
    ///                   + \vec{f}_{1}(\vec{x}) \cdot \nabla v^\mathrm{out} \ \mathrm{dS}
    ///           + \oint_{\partial \Omega} \vec{f}_{0}(\vec{x}) \cdot \nabla v^\mathrm{in} \ \mathrm{dS}
    /// ```
    /// where <em>v</em> denotes the test function.
    /// </summary>
    public interface IEdgeSource_GradV : IEdgeForm, IInnerEdgeSource_GradV, IBoundaryEdgeSource_GradV {
    }

    /// <summary>
    /// <see cref="IEdgeSource_GradV"/>
    /// </summary>
    public interface IInnerEdgeSource_GradV : IInnerEdgeForm {
        /// <summary>
        /// the point-wise source term \f$ \vec{f}_{i}(\vec{x})\f$  on interior edges \f$ \Gamma_{\mathrm{int}}\f$ .
        /// </summary>
        /// <param name="efp">parameters on which  \f$ \vec{f}\f$  may depend on.</param>
        /// <param name="GradV">
        /// Output: the value of \f$ \vec{f}_{i}(\vec{x})\f$ , which will be multiplied with the test function gradient.
        /// - 1st index: edge index
        /// - 2nd index: quadrature node
        /// - 3rd index: in and out - coefficients with respect to test function, i.e. index \f$ i\f$ 
        ///              ('V': index 0 corresponds to IN-cell, index 1 corresponds to OUT-cell).
        /// - 4th index: spatial direction of test function ('V') gradient <br/>
        /// </param>
        void InternalEdge_GradV(ref EdgeFormParams efp, MultidimensionalArray GradV);
    }

    /// <summary>
    /// <see cref="IEdgeSource_GradV"/>
    /// </summary>
    public interface IBoundaryEdgeSource_GradV : IBoundaryEdgeForm {
        /// <summary>
        /// the point-wise source term \f$ \vec{f}_{0}(\vec{x})\f$  on boundary edges \f$ \partial \Omega\f$ .
        /// </summary>
        /// <param name="efp">parameters on which  \f$ \vec{f}\f$  may depend on.</param> 
        /// <param name="GradV">
        /// Output: the value of \f$ \vec{f}_{0}(\vec{x})\f$ , which will be multiplied with the test function gradient. <br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// 3rd index: spatial direction of test function ('V') gradient <br/>
        /// </param>
        void BoundaryEdge_GradV(ref EdgeFormParams efp, MultidimensionalArray GradV);
    }

    /// <summary>
    /// a form of the type
    /// ```math 
    ///    a(U,v) = \oint_{\Gamma_{\mathrm{int}}} 
    ///              {f}^\mathrm{in} (\vec{x}, U^\mathrm{in}, U^\mathrm{out}, \nabla U^\mathrm{in}, \nabla U^\mathrm{out}) v^\mathrm{in} 
    ///            + {f}^\mathrm{out}(\vec{x}, U^\mathrm{in}, U^\mathrm{out}, \nabla U^\mathrm{in}, \nabla U^\mathrm{out}) v^\mathrm{out} 
    ///           \ \mathrm{dS}
    ///         +
    ///          \oint_{\partial \Omega} 
    ///              {f}_{0}(\vec{x}, U^\mathrm{in}, \nabla U^\mathrm{in}) v^\mathrm{in} 
    ///           \ \mathrm{dS}          
    /// ```
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable) and 
    /// $` U = (u_0, \ldots, u_{L-1} ) `$  denotes the trial functions (correspond to domain variable, defined by the 
    /// argument ordering <see cref="IEquationComponent.ArgumentOrdering"/>).
    /// </summary>
    public interface INonlinEdgeForm_V : IEdgeForm, INonlinInnerEdgeForm_V, INonlinBoundaryEdgeForm_V {
    }

    /// <summary>
    /// <see cref="INonlinEdgeForm_V"/>
    /// </summary>
    public interface INonlinInnerEdgeForm_V : IInnerEdgeForm {
        /// <summary>
        /// the values of \f$ {f}^{*}(\ldots) \f$ on interior edges on \f$ \Gamma_{\mathrm{int}}\f$ .
        /// </summary>
        /// <param name="efp"></param>
        /// <param name="fin">
        /// - output: the values of \f$ {f}^\mathrm{in}(\ldots)\f$ :
        /// - 1st index: edge index
        /// - 2nd index: quadrature node
        /// </param>
        /// <param name="fot">
        /// output: the values of \f$ {f}^\mathrm{out}(\ldots)\f$
        /// </param>
        /// <param name="Uin">
        /// - input: the values of \f$ U^\mathrm{in}\f$ 
        /// - 1st index: component index of \f$ U\f$
        /// - 2nd index: edge index
        /// - 3rd index: quadrature node
        /// </param>
        /// <param name="Uout">
        /// analog to <paramref name="Uin"/>, for 'out'-values.<br/>
        /// </param>
        /// <param name="GradUin">
        /// input: the values of \f$ \nabla U^\mathrm{in}\f$ <br/>
        /// 1st index: component index of \f$ U\f$ <br/>
        /// 2nd index: edge index<br/>
        /// 3rd index: quadrature node<br/>
        /// 4th index: spatial direction of derivative<br/>
        /// </param>
        /// <param name="GradUout">
        /// analog to <paramref name="GradUin"/>, for 'out'-values.
        /// </param>
        void NonlinInternalEdge_V(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout,
            MultidimensionalArray fin, MultidimensionalArray fot);

    }

    /// <summary>
    /// <see cref="INonlinEdgeForm_V"/>
    /// </summary>
    public interface INonlinBoundaryEdgeForm_V : IBoundaryEdgeForm {

        /// <summary>
        /// the values of \f$ {f}^\mathrm{in}(\ldots)\f$   on boundary edges on \f$ \partial \Omega\f$ .
        /// </summary>
        /// <param name="efp"></param>
        /// <param name="fin">
        /// output: the values of \f$ {f}^\mathrm{in}(\ldots)\f$ :<br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// </param>
        /// <param name="Uin">
        /// input: the values of \f$ U^\mathrm{in}\f$ <br/>
        /// 1st index: component index of \f$ U\f$ <br/>
        /// 2nd index: edge index<br/>
        /// 3rd index: quadrature node<br/>
        /// </param>
        /// <param name="GradUin">
        /// input: the values of \f$ \nabla U^\mathrm{in}\f$ <br/>
        /// 1st index: component index of \f$ U\f$ <br/>
        /// 2nd index: edge index<br/>
        /// 3rd index: quadrature node<br/>
        /// 4th index: spatial direction of derivative<br/>
        /// </param>
        void NonlinBoundaryEdge_V(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin,
            MultidimensionalArray fin);
    }


    /// <summary>
    /// a form of the type
    /// ```math 
    ///    a(U,v) = \oint_{\Gamma_{\mathrm{int}}} 
    ///              \vec{f}^\mathrm{in}  (\vec{x}, U^\mathrm{in}, U^\mathrm{out}, \nabla U^\mathrm{in}, \nabla U^\mathrm{out}) \cdot \nabla v^\mathrm{in} 
    ///            + \vec{f}^\mathrm{out} (\vec{x}, U^\mathrm{in}, U^\mathrm{out}, \nabla U^\mathrm{in}, \nabla U^\mathrm{out}) \cdot \nabla v^\mathrm{out} 
    ///           \ \mathrm{dS}
    ///         +
    ///          \oint_{\partial \Omega} 
    ///              \vec{f}^\mathrm{in}  (\vec{x}, U^\mathrm{in}, \nabla U^\mathrm{in}) \cdot \nabla v^\mathrm{in} 
    ///           \ \mathrm{dS}          
    /// ```
    /// where <em>v</em> denotes the test function (corresponds to co-domain variable) and 
    /// \f$ U = (u_0, \ldots, u_{L-1} )\f$  denotes the trial functions (correspond to domain variable, defined by the 
    /// argument ordering <see cref="IEquationComponent.ArgumentOrdering"/>).
    /// </summary>
    public interface INonlinEdgeForm_GradV : IEdgeForm, INonlinInnerEdgeForm_GradV, INonlinBoundaryEdgeForm_GradV {
    }

    /// <summary>
    /// <see cref="INonlinEdgeForm_GradV"/>
    /// </summary>
    public interface INonlinInnerEdgeForm_GradV : IInnerEdgeForm {
        /// <summary>
        /// the values of \f$ \vec{f}^{*}(\ldots)\f$   on interior edges on \f$ \Gamma_{\mathrm{int}}\f$ .
        /// </summary>
        /// <param name="efp"></param>
        /// <param name="fIN">
        /// output: the values of \f$ {f}^\mathrm{in}(\ldots)\f$ :<br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// 3rd index: component index of  \f$ \vec{f}^\mathrm{in}(\ldots)\f$ . <br/>
        /// </param>
        /// <param name="fOT">
        /// like <paramref name="fIN"/>, for OUT-cells.
        /// </param>
        /// <param name="Uin">
        /// input: the values of \f$ U^\mathrm{in}\f$ <br/>
        /// 1st index: component index of \f$ U\f$ <br/>
        /// 2nd index: edge index<br/>
        /// 3rd index: quadrature node<br/>
        /// </param>
        /// <param name="Uout">
        /// analog to <paramref name="Uin"/>, for 'out'-values.<br/>
        /// </param>
        /// <param name="GradUin">
        /// input: the values of \f$ \nabla U^\mathrm{in}\f$ <br/>
        /// 1st index: component index of \f$ U\f$ <br/>
        /// 2nd index: edge index<br/>
        /// 3rd index: quadrature node<br/>
        /// 4th index: spatial direction of derivative<br/>
        /// </param>
        /// <param name="GradUout">
        /// analog to <paramref name="GradUin"/>, for 'out'-values.
        /// </param>
        void NonlinInternalEdge_GradV(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout,
            MultidimensionalArray fIN, MultidimensionalArray fOT);
    }


    /// <summary>
    /// <see cref="INonlinEdgeForm_GradV"/>
    /// </summary>
    public interface INonlinBoundaryEdgeForm_GradV : IBoundaryEdgeForm {
        /// <summary>
        /// the values of \f$ \vec{f}_{0}(\ldots)\f$   on boundary edges on \f$ \partial \Omega\f$ .
        /// </summary>
        /// <param name="efp"></param>
        /// <param name="f">
        /// output: the values of \f$ {f}_{0}(\ldots)\f$ :<br/>
        /// 1st index: edge index<br/>
        /// 2nd index: quadrature node<br/>
        /// 3rd index: component index of  \f$ \vec{f}_{i}(\ldots)\f$ . <br/>
        /// </param>
        /// <param name="Uin">
        /// input: the values of \f$ U^\mathrm{in}\f$ <br/>
        /// 1st index: component index of \f$ U\f$ <br/>
        /// 2nd index: edge index<br/>
        /// 3rd index: quadrature node<br/>
        /// </param>
        /// <param name="GradUin">
        /// input: the values of \f$ \nabla U^\mathrm{in}\f$ 
        /// - 1st index: component index of \f$ U\f$ 
        /// - 2nd index: edge index
        /// - 3rd index: quadrature node
        /// - 4th index: spatial direction of derivative
        /// </param>
        void NonlinBoundaryEdge_GradV(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin,
            MultidimensionalArray f);
    }
}
