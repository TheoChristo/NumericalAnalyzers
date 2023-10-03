using System;
using System.Diagnostics;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.NumericalAnalyzers.Logging;
using System.Collections;
using MGroup.LinearAlgebra.Iterative;
using System.Collections.Generic;
using System.Linq;

namespace MGroup.NumericalAnalyzers.Dynamic
{
	public class PseudoTransientAnalyzer2 : INonLinearParentAnalyzer, IStepwiseAnalyzer
	{
		private const string TIME = TransientLiterals.TIME;
		private const string CURRENTTIMESTEP = "Current timestep";
		private const string CURRENTSOLUTION = "Current solution";
		private const string PREVIOUSSOLUTION = "Previous solution";
		private const string FIRSTORDERSOLUTION = "First order derivative of solution";
		private const string SECONDORDERSOLUTION = "Second order derivative of solution";

		/// <summary>
		/// This class implements a Pseudo-Transient Analyzer
		/// Authors: George Stavroulakis, Theofilos Christodoulou.
		/// </summary>
		private readonly double timeStep;

		/// <summary>
		/// This class makes the appropriate arrangements for the solution of linear dynamic equations
		/// according to implicit Newmark method
		/// Authors: George Stavroulakis, Theofilos Christodoulou.
		/// </summary>
		private readonly double totalTime;
		//private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ITransientAnalysisProvider provider;
		private IGlobalVector rhs;
		private int currentStep;
		private DateTime start, end;
		private GenericAnalyzerState currentState;
		private IList<IterativeStatistics> analysisStatistics;
		private IGlobalVector[] solutions;
		private IGlobalVector zeroOrderDerivativeSolutionOfPreviousStep;
		private bool newMarkVelocity;
		//Newmark velocity
		private readonly double beta;
		private readonly double gamma;
		private readonly double a0;
		private readonly double a1;
		private readonly double a2;
		private readonly double a3;
		private readonly double a4;
		private readonly double a5;
		private readonly double a6;
		private readonly double a7;

		//extra Generalized alpha parameters
		private readonly double am, af;
		private readonly double a0N;
		private readonly double a2N;
		private readonly double a3N;
		private readonly double a6N;
		private readonly double a7N;







		/// <summary>
		/// Creates an instance that uses a specific problem type and an appropriate child analyzer for the construction of the system of equations arising from the actual physical problem.
		/// </summary>
		/// <param name="model">Instance of the model to be solved.</param>
		/// <param name="provider">Instance of the problem type to be solver.</param>
		/// <param name="childAnalyzer">Instance of the child analyzer that will handle the solution of the system of equations.</param>
		/// <param name="timeStep">Instance of the time step of the method that will be initialized.</param>
		/// <param name="totalTime">Instance of the total time of the method that will be initialized.</param>
		private PseudoTransientAnalyzer2(IAlgebraicModel algebraicModel, ITransientAnalysisProvider provider,
			IChildAnalyzer childAnalyzer, double timeStep, double totalTime, int currentStep, double alpha, double delta, double am, double af, bool newMarkVelocity)
		{
			this.algebraicModel = algebraicModel;
			this.provider = provider;
			this.ChildAnalyzer = childAnalyzer;
			this.timeStep = timeStep;
			this.totalTime = totalTime;
			this.currentStep = currentStep;
			this.ChildAnalyzer.ParentAnalyzer = this;
			this.analysisStatistics = Enumerable.Range(0, Steps).Select(x => new IterativeStatistics() { AlgorithmName = "Pseudo-transient analyzer" }).ToArray();
			this.newMarkVelocity = newMarkVelocity;
			if (newMarkVelocity)
			{
				this.beta = alpha;
				this.gamma = delta;
				a0 = 1 / (alpha * timeStep * timeStep);
				a1 = delta / (alpha * timeStep);
				a2 = 1 / (alpha * timeStep);
				a3 = (1 / (2 * alpha)) - 1;
				a4 = (delta / alpha) - 1;
				a5 = timeStep * 0.5 * ((delta / alpha) - 2);
				a6 = timeStep * (1 - delta);
				a7 = delta * timeStep;
			}
			else
			{
			    this.beta = alpha;
				this.gamma = delta;
				this.am = am;
				this.af = af;
				a0N = 1 / (alpha * timeStep * timeStep);
				a0 = (1 - am) / (alpha * timeStep * timeStep);
				a1 = delta * (1 - af) / (alpha * timeStep);
				a2N = 1 / (alpha * timeStep);
				a2 = (1 - am) / (alpha * timeStep);
				a3N = (1 / (2 * alpha)) - 1;
				a3 = ((1 - am) / (2 * alpha)) - 1;
				a4 = ((delta - delta * af) / alpha) - 1;
				a5 = (1 - af) * timeStep * 0.5 * ((delta / alpha) - 2);
				a6N = timeStep * (1 - delta);
				a7N = delta * timeStep;
			}


		}

		public IAnalysisWorkflowLog[] Logs => null;

		public IGlobalVector CurrentAnalysisResult { get => ChildAnalyzer.CurrentAnalysisResult.Copy(); }

		public ImplicitIntegrationAnalyzerLog ResultStorage { get; set; }

		public IChildAnalyzer ChildAnalyzer { get; }

		public int CurrentStep { get => currentStep; }

		public int Steps { get => (int)(totalTime / timeStep); }

		GenericAnalyzerState IAnalyzer.CurrentState
		{
			get => currentState;
			set
			{
				currentState = value;
				currentStep = (int)currentState.StateValues[CURRENTTIMESTEP];

				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = false;
				currentState.StateVectors[PREVIOUSSOLUTION].CheckForCompatibility = false;
				currentState.StateVectors[FIRSTORDERSOLUTION].CheckForCompatibility = false;
				currentState.StateVectors[SECONDORDERSOLUTION].CheckForCompatibility = false;

				solutions[0].CopyFrom(currentState.StateVectors[CURRENTSOLUTION]);
				zeroOrderDerivativeSolutionOfPreviousStep.CopyFrom(currentState.StateVectors[PREVIOUSSOLUTION]);
				solutions[(int)DifferentiationOrder.First].CopyFrom(currentState.StateVectors[FIRSTORDERSOLUTION]);
				solutions[(int)DifferentiationOrder.Second].CopyFrom(currentState.StateVectors[SECONDORDERSOLUTION]);

				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = true;
				currentState.StateVectors[PREVIOUSSOLUTION].CheckForCompatibility = true;
				currentState.StateVectors[FIRSTORDERSOLUTION].CheckForCompatibility = true;
				currentState.StateVectors[SECONDORDERSOLUTION].CheckForCompatibility = true;
			}
		}

		public IList<IterativeStatistics> AnalysisStatistics => analysisStatistics;

		/// <summary>
		/// Makes the proper solver-specific initializations before the solution of the linear system of equations. This method MUST be called before the actual solution of the aforementioned system
		/// </summary>
		public void BuildMatrices()
		{
			algebraicModel.LinearSystem.Matrix = provider.GetMatrix(DifferentiationOrder.Zero);
		}

		/// <summary>
		/// Calculates inertia forces and damping forces.
		/// </summary>
		public IGlobalVector GetOtherRhsComponents(IGlobalVector currentSolution) => algebraicModel.CreateZeroVector();

		/// <summary>
		/// Initializes the models, the solvers, child analyzers, builds the matrices, assigns loads and initializes right-hand-side vectors.
		/// </summary>
		public void Initialize(bool isFirstAnalysis = true)
		{
			if (isFirstAnalysis)
			{
				// Connect data structures of model is called by the algebraic model
				algebraicModel.OrderDofs();
			}

			BuildMatrices();
			InitializeInternalVectors();
			InitializeRhs();

			if (ChildAnalyzer == null)
			{
				throw new InvalidOperationException("Pseudo-Transient analyzer must contain an embedded analyzer.");
			}

			ChildAnalyzer.Initialize(isFirstAnalysis);
		}

		private void SolveCurrentTimestep()
		{
			Debug.WriteLine("Pseudo-Transient Analyzer step: {0}", currentStep);

			IGlobalVector rhsVector = provider.GetRhs(currentStep * timeStep);
			ChildAnalyzer.CurrentAnalysisLinearSystemRhs.CopyFrom(rhsVector);

			InitializeRhs();
			CalculateRhsImplicit();

			start = DateTime.Now;
			ChildAnalyzer.Initialize(false);
			ChildAnalyzer.Solve();
			analysisStatistics[currentStep] = ChildAnalyzer.AnalysisStatistics;
			end = DateTime.Now;
			Debug.WriteLine("Pseudo-Transient Analyzer elapsed time: {0}", end - start);
		}

		/// <summary>
		/// Solves the linear system of equations by calling the corresponding method of the specific solver attached during construction of the current instance
		/// </summary>
		public void Solve()
		{
			for (int i = 0; i < Steps; ++i)
			{
				SolveCurrentTimestep();
				AdvanceStep();
			}
		}

		/// <summary>
		/// Calculates the right-hand-side of the implicit dyanmic method. This will be used for the solution of the linear dynamic system.
		/// </summary>
		private void CalculateRhsImplicit()
		{
			ChildAnalyzer.CurrentAnalysisLinearSystemRhs.CopyFrom(rhs);
		}

		private void InitializeInternalVectors()
		{
			rhs = algebraicModel.CreateZeroVector();
			solutions = new IGlobalVector[3];
			for (int i = 0; i < 3; i++)
			{
				solutions[i] = algebraicModel.CreateZeroVector();
			}

		}

		private void InitializeRhs()
		{
			rhs.CopyFrom(ChildAnalyzer.CurrentAnalysisLinearSystemRhs);
		}

		private void UpdateResultStorages(DateTime start, DateTime end)
		{
			if (ResultStorage != null)
			{
				foreach (var l in ChildAnalyzer.Logs)
				{
					ResultStorage.StoreResults(start, end, l);
				}
			}
		}

		GenericAnalyzerState CreateState()
		{
			currentState = new GenericAnalyzerState(this,
				new (string, IGlobalVector)[]
				{
					(CURRENTSOLUTION, solutions[0]),
					(PREVIOUSSOLUTION, zeroOrderDerivativeSolutionOfPreviousStep),
					(FIRSTORDERSOLUTION, solutions[(int)DifferentiationOrder.First]),
					(SECONDORDERSOLUTION, solutions[(int)DifferentiationOrder.Second]),
				},
				new[]
				{
					(TIME, (double)currentStep * (double) timeStep),
					(CURRENTTIMESTEP, (double)currentStep),
				});

			return currentState;
		}

		IHaveState ICreateState.CreateState() => CreateState();
		GenericAnalyzerState IAnalyzer.CreateState() => CreateState();

		/// <summary>
		/// Solves the linear system of equations of the current timestep
		/// </summary>
		void IStepwiseAnalyzer.Solve()
		{
			SolveCurrentTimestep();
		}

		public void AdvanceStep()
		{
			Debug.WriteLine("Advancing step");
			UpdateHigherOrderDerivatives();
		    UpdateResultStorages(start, end);
			currentStep++;
		}

		private void UpdateHigherOrderDerivatives()
		{
			if (newMarkVelocity)
			{
				zeroOrderDerivativeSolutionOfPreviousStep.CopyFrom(solutions[0]);
				solutions[0].CopyFrom(ChildAnalyzer.CurrentAnalysisResult);

				var secondOrderDerivativeOfSolutionOfPreviousStep = solutions[(int)DifferentiationOrder.Second].Copy();

				solutions[(int)DifferentiationOrder.Second] = solutions[0].Subtract(zeroOrderDerivativeSolutionOfPreviousStep);
				solutions[(int)DifferentiationOrder.Second].LinearCombinationIntoThis(a0, solutions[(int)DifferentiationOrder.First], -a2);
				solutions[(int)DifferentiationOrder.Second].AxpyIntoThis(secondOrderDerivativeOfSolutionOfPreviousStep, -a3);

				solutions[(int)DifferentiationOrder.First].AxpyIntoThis(secondOrderDerivativeOfSolutionOfPreviousStep, a6);
				solutions[(int)DifferentiationOrder.First].AxpyIntoThis(solutions[(int)DifferentiationOrder.Second], a7);
			}
			else
			{
				zeroOrderDerivativeSolutionOfPreviousStep.CopyFrom(solutions[0]);
				solutions[0].CopyFrom(ChildAnalyzer.CurrentAnalysisResult);

				var secondOrderDerivativeOfSolutionOfPreviousStep = solutions[(int)DifferentiationOrder.Second].Copy();

				solutions[(int)DifferentiationOrder.Second] = solutions[0].Subtract(zeroOrderDerivativeSolutionOfPreviousStep);
				solutions[(int)DifferentiationOrder.Second].LinearCombinationIntoThis(a0N, solutions[(int)DifferentiationOrder.First], -a2N);
				solutions[(int)DifferentiationOrder.Second].AxpyIntoThis(secondOrderDerivativeOfSolutionOfPreviousStep, -a3N);

				solutions[(int)DifferentiationOrder.First].AxpyIntoThis(secondOrderDerivativeOfSolutionOfPreviousStep, a6N);
				solutions[(int)DifferentiationOrder.First].AxpyIntoThis(solutions[(int)DifferentiationOrder.Second], a7N);
			}

		}

		public class Builder
		{
			private readonly double timeStep;
			private readonly double totalTime;
			private readonly IChildAnalyzer childAnalyzer;
			private readonly IAlgebraicModel algebraicModel;
			private readonly ITransientAnalysisProvider provider;
			private int currentStep = 0;
			//Newmark velocity
			private double beta = 0.25;
			private double gamma = 0.5;
			//extra Generalised alpha parameters
			private double am = 0;
			private double af = 0;
			private bool newMarkVelocity;


			public Builder(IAlgebraicModel algebraicModel, ITransientAnalysisProvider provider,
				IChildAnalyzer childAnalyzer, double timeStep, double totalTime, bool NewmarkVelocity, int currentStep = 0)
			{
				this.algebraicModel = algebraicModel;
				this.provider = provider;
				this.childAnalyzer = childAnalyzer;
				this.currentStep = currentStep;

				this.timeStep = timeStep;
				this.totalTime = totalTime;
				this.newMarkVelocity = NewmarkVelocity;
			}

			public void SetNewmarkParameters(double beta, double gamma, bool allowConditionallyStable = false)
			{
				if (!allowConditionallyStable)
				{
					if (gamma < 0.5)
					{
						throw new ArgumentException(
						"Newmark delta has to be bigger than 0.5 to ensure unconditional stability.");
					}

					if (beta < 0.25)
					{
						throw new ArgumentException(
						"Newmark alpha has to be bigger than 0.25 to ensure unconditional stability.");
					}
				}
				if (gamma < 0.5)
				{
					throw new ArgumentException("Newmark delta has to be bigger than 0.5.");
				}

				double aLimit = 0.25 * Math.Pow(0.5 + gamma, 2);
				if (beta < aLimit)
				{
					throw new ArgumentException($"Newmark alpha has to be bigger than {aLimit}.");
				}

				this.gamma = gamma;
				this.beta = beta;
			}

			public void SetNewmarkParametersForCentralDifferences()
			{
				gamma = 0.5;
				beta = 0.0;
			}

			public void SetNewmarkParametersForConstantAcceleration()
			{
				gamma = 0.5;
				beta = 0.25;
			}

			public void SetNewmarkParametersForLinearAcceleration()
			{
				gamma = 0.5;
				beta = 1.0 / 6.0;
			}

			//generalized alpha extra paarameters 
			public void SetSpectralRadius(double spectralRadius)
			{
				this.am = (2 * spectralRadius - 1) / (spectralRadius + 1);
				this.af = spectralRadius / (spectralRadius + 1);
				this.beta = 0.25 * Math.Pow(1 - am + af, 2);
				this.gamma = 0.5 - am + af;
			}


			public PseudoTransientAnalyzer2 Build()
			{
				
					return new PseudoTransientAnalyzer2(algebraicModel, provider, childAnalyzer, timeStep, totalTime, currentStep, beta, gamma, am, af, newMarkVelocity);
				
			}
		}
	}
}
