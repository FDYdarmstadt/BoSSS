using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools {

    
    /// <summary>
    /// Driver options that control level-set evolution
    /// </summary>
    public enum LevelSetEvolution {

        /// <summary>
        /// turn level-set evolution off
        /// </summary>
        None,

        /// <summary>
        /// Prescribed level-set, the initial value is set every time.
        /// </summary>
        Prescribed,

        /// <summary>
        /// Prescribed level-set (wave-like) projected from imported amplitude values.
        /// </summary>
        PrescribedLSwave,

        /// <summary>
        /// The level-set is constructed from a graph which describes the interface,
        /// e.g. height or radius.
        /// The graph is described by a Fourier series;
        /// Can only be used for special cases, e.g. 2D fluid layers or bubbles/droplets
        /// </summary>
        /// <remarks>
        /// For details:
        /// <see cref="BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.FourierEvolver"/>.
        /// </remarks>
        Fourier,

        /// <summary>
        /// Extension velocity based (fast) marching algorithm.
        ///
        /// In contrast to other methods, this should be reliable for general interface topologies
        /// and it is also known to work for contact line setups 
        /// and should also work for non-material interfaces with a velocity jump.
        /// </summary>
        /// <remarks>
        /// For details:
        /// <see cref="BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.FastMarchingEvolver"/>.
        /// </remarks>
        FastMarching,

        /// <summary>
        /// Level set evolution where the level set field is treated by a scalar convection.
        ///
        /// This is known to be not very precise, especially for high viscosity ratios and/or
        /// non-material interfaces.
        /// This is because, due to the mismatch of approximation spaces -- i.e. the 
        /// velocity is XDG, while the level-set-field is DG -- a large error is induced
        /// when the XDG velocity is projected onto the DG space.
        /// </summary>
        ScalarConvection,

        /// <summary>
        /// The Level-Set is moved by an extension velocity idea using a PDE (cf. Utz et. al 2017, Utz and Kummer 2017)
        /// ReInitialization is done by fast marching, no ReInit on Cut Cells!
        /// -> TODO
        /// </summary>
        ExtensionVelocity,

        /// <summary>
        /// Phasefield level set by Cahn-Hilliard equation
        /// </summary>
        Phasefield,

        /// <summary>
        /// Level Set is moved by advection particles.
        /// This evolution features three subtypes <see cref="BoSSS.Application.SemiLagrangianLevelSetTestSuite.LagrangianMode"/>
        /// </summary>
        /// <remarks>
        /// 
        /// </remarks>
        SemiLagrangianLevelSet,

        /// <summary>
        /// Spline Level Set
        /// </summary>
        SplineLS,
                
        /// <summary>
        /// An extension velocity computed from a Stokes equation; 
        /// In contrast to other methods, this should be reliable for general interface topologies
        /// and it should work for contact line setups as well as non-material interfaces with a velocity jump.
        /// </summary>
        /// <remarks>
        /// For details:
        /// <see cref="BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.StokesExtensionEvolver"/>.
        /// </remarks>
        StokesExtension,

        /// <summary>
        /// An extension velocity computed from a Laplace equation; 
        /// Same as <see cref="StokesExtension"/> but lacks artificial pressure
        /// </summary>
        /// <remarks>
        /// For details:
        /// <see cref="BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.StokesExtensionEvolver"/>.
        /// </remarks>
        LaplaceExtension,

        /// <summary>
        /// A level set formulation for the surface of rigid objects. The Movement depends only on the position and orientation of the object.
        /// </summary>
        RigidObject,

        /// <summary>
        /// Parameterized Level Set
        /// </summary>
        /// <remarks>
        /// For details:
        /// <see cref="BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.ParameterizedInterfaceEvolver"/>.
        /// </remarks>
        ParameterizedLevelSet,
    }
}
