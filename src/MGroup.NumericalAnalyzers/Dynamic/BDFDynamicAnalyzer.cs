using System;
using System.Collections.Generic;
using System.Diagnostics;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.DataStructures;

namespace MGroup.NumericalAnalyzers.Dynamic
{
	/// <summary>
	/// This class makes the appropriate arrangements for the solution of linear dynamic equations
	/// according to the Backward Differentiation Formula
	/// Authors: Orestis Papas, Theofilos Christodoulou
	/// </summary>
	public class BDFDynamicAnalyzer : INonLinearParentAnalyzer, IStepwiseAnalyzer
	{
		private const string CURRENTTIMESTEP = "Current timestep";
		private const string CURRENTSOLUTION = "Current solution";
		private const string SOLUTION_N_1 = "Previous solution (n-1)";
		private const string SOLUTION_N_2 = "Previous solution (n-2)";
		private const string SOLUTION_N_3 = "Previous solution (n-3)";
		private const string SOLUTION_N_4 = "Previous solution (n-4)";
		private const string SOLUTION_N_5 = "Previous solution (n-5)";
		private const string FIRSTORDERSOLUTION = "First order derivative of solution";
		
		private readonly double timeStep;

		private readonly double totalTime;

		private int currentTimeStep = 0;

		private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ISolver solver;
		private readonly ITransientAnalysisProvider provider;
		private IGlobalVector rhs;
		private IGlobalVector solution;
		private IGlobalVector[] solutionOfPreviousStep;
		private IGlobalVector firstOrderDerivativeOfSolution;
		private IGlobalVector firstOrderDerivativeOfSolutionForRhs;
		private IGlobalVector firstOrderDerivativeComponentOfRhs;
		private DateTime start, end;
		private int currentStep;
		private GenericAnalyzerState currentState;
		private int lastUpdatedStep = -1;

		/// <summary>
		/// Creates an instance that uses a specific problem type and an appropriate child analyzer for the construction of the system of equations arising from the actual physical problem
		/// </summary>
		/// <param name="model">Instance of the model to be solved</param>
		/// <param name="solver">Instance of the solver that will handle the solution of the system of equations</param>
		/// <param name="provider">Instance of the problem type to be solver</param>
		/// <param name="childAnalyzer">Instance of the child analyzer that will handle the solution of the system of equations</param>
		/// <param name="timeStep">Instance of the time step of the method that will be initialized</param>
		/// <param name="totalTime">Instance of the total time of the method that will be initialized</param>
		/// <param name="bdfOrder">Order of the scheme [1,5]</param>

		private BDFDynamicAnalyzer(IModel model, IAlgebraicModel algebraicModel, ISolver solver, ITransientAnalysisProvider provider,
			IChildAnalyzer childAnalyzer, double timeStep, double totalTime, int bdfOrder, int currentTimeStep)
		{
			this.model = model;
			this.algebraicModel = algebraicModel;
			this.solver = solver;
			this.provider = provider;
			this.ChildAnalyzer = childAnalyzer;
			this.currentTimeStep = currentTimeStep;
			this.timeStep = timeStep;
			this.totalTime = totalTime;
			this.ChildAnalyzer.ParentAnalyzer = this;
			this.BDFOrder = bdfOrder;

			if (BDFOrder < 1 || BDFOrder > 5)
			{
				throw new ArgumentException("Wrong BDF order. Must be in [1,5]");
			}
		}

		public int BDFOrder { get; }

		public IAnalysisWorkflowLog[] Logs => null;

		public ImplicitIntegrationAnalyzerLog ResultStorage { get; set; }

		public IChildAnalyzer ChildAnalyzer { get; }

		public int CurrentStep { get => currentTimeStep; }

		public int Steps { get => (int)(totalTime / timeStep); }

		GenericAnalyzerState IAnalyzer.CurrentState
		{
			get => currentState;
			set
			{
				currentState = value;
				currentStep = (int)currentState.StateValues[CURRENTTIMESTEP];
				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = false;
				if (solutionOfPreviousStep.Length > 0) currentState.StateVectors[SOLUTION_N_1].CheckForCompatibility = false;
				if (solutionOfPreviousStep.Length > 1) currentState.StateVectors[SOLUTION_N_2].CheckForCompatibility = false;
				if (solutionOfPreviousStep.Length > 2) currentState.StateVectors[SOLUTION_N_3].CheckForCompatibility = false;
				if (solutionOfPreviousStep.Length > 3) currentState.StateVectors[SOLUTION_N_4].CheckForCompatibility = false;
				if (solutionOfPreviousStep.Length > 4) currentState.StateVectors[SOLUTION_N_5].CheckForCompatibility = false;
				currentState.StateVectors[FIRSTORDERSOLUTION].CheckForCompatibility = false;
				
				solution.CopyFrom(currentState.StateVectors[CURRENTSOLUTION]);

				solutionOfPreviousStep[0].CopyFrom(currentState.StateVectors[SOLUTION_N_1]);
				if (solutionOfPreviousStep.Length > 1)
				{
					solutionOfPreviousStep[1].CopyFrom(currentState.StateVectors[SOLUTION_N_2]);
				}
				if (solutionOfPreviousStep.Length > 2)
				{
					solutionOfPreviousStep[2].CopyFrom(currentState.StateVectors[SOLUTION_N_3]);
				}
				if (solutionOfPreviousStep.Length > 3)
				{
					solutionOfPreviousStep[3].CopyFrom(currentState.StateVectors[SOLUTION_N_4]);
				}
				if (solutionOfPreviousStep.Length > 4)
				{
					solutionOfPreviousStep[4].CopyFrom(currentState.StateVectors[SOLUTION_N_5]);
				}

				firstOrderDerivativeOfSolution.CopyFrom(currentState.StateVectors[FIRSTORDERSOLUTION]);
				
				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = true;
				if (solutionOfPreviousStep.Length > 0) currentState.StateVectors[SOLUTION_N_1].CheckForCompatibility = true;
				if (solutionOfPreviousStep.Length > 1) currentState.StateVectors[SOLUTION_N_2].CheckForCompatibility = true;
				if (solutionOfPreviousStep.Length > 2) currentState.StateVectors[SOLUTION_N_3].CheckForCompatibility = true;
				if (solutionOfPreviousStep.Length > 3) currentState.StateVectors[SOLUTION_N_4].CheckForCompatibility = true;
				if (solutionOfPreviousStep.Length > 4) currentState.StateVectors[SOLUTION_N_5].CheckForCompatibility = true;
				currentState.StateVectors[FIRSTORDERSOLUTION].CheckForCompatibility = true;
			}
		}

		GenericAnalyzerState CreateState()
		{
			currentState = new GenericAnalyzerState(this,
				new[]
				{
					(CURRENTSOLUTION, solution),
					(SOLUTION_N_1, solutionOfPreviousStep[0]),
					(SOLUTION_N_2, solutionOfPreviousStep.Length > 1 ? solutionOfPreviousStep[1] : null),
					(SOLUTION_N_3, solutionOfPreviousStep.Length > 2 ? solutionOfPreviousStep[2] : null),
					(SOLUTION_N_4, solutionOfPreviousStep.Length > 3 ? solutionOfPreviousStep[3] : null),
					(SOLUTION_N_5, solutionOfPreviousStep.Length > 4 ? solutionOfPreviousStep[4] : null),
					(FIRSTORDERSOLUTION, firstOrderDerivativeOfSolution),
				},
				new[]
				{
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
		/// <summary>
		/// Makes the proper solver-specific initializations before the solution of the linear system of equations. This method MUST be called before the actual solution of the aforementioned system
		/// </summary>
		public void BuildMatrices()
		{
			var coeffs = GetTransientAnalysisCoefficients();
			provider.LinearCombinationOfMatricesIntoEffectiveMatrixNoOverwrite(coeffs);
		}

		private TransientAnalysisCoefficients GetTransientAnalysisCoefficients()
		{
			int bdfOrderInternal = Math.Min(currentTimeStep + 1, BDFOrder);
			double timeStepNumerator;
			switch (bdfOrderInternal)
			{
				case 1:
					timeStepNumerator = 1;
					break;
				case 2:
					timeStepNumerator = 3.0 / 2.0;
					break;
				case 3:
					timeStepNumerator = 11.0 / 6.0;
					break;
				case 4:
					timeStepNumerator = 25.0 / 12.0;
					break;
				case 5:
					timeStepNumerator = 137.0 / 60.0;
					break;
				default:
					throw new ArgumentException("Wrong BDF Order");

			}

			return new TransientAnalysisCoefficients
			{
				SecondOrderDerivativeCoefficient = 0d,
				FirstOrderDerivativeCoefficient = timeStepNumerator / timeStep,
				ZeroOrderDerivativeCoefficient = 1,
			};
		}


		/// <summary>
		/// Calculates inertia forces and damping forces.
		/// </summary>
		public IGlobalVector GetOtherRhsComponents(IGlobalVector currentSolution)
		{
			return currentSolution;
		}

		/// <summary>
		/// Initializes the models, the solvers, child analyzers, builds the matrices, assigns loads and initializes right-hand-side vectors.
		/// </summary>
		public void Initialize(bool isFirstAnalysis = true)
		{//TODO: get initial conditions!
			if (isFirstAnalysis)
			{
				//provider.GetProblemDofTypes();
				model.ConnectDataStructures();
				algebraicModel.OrderDofs();
			}

			BuildMatrices();

			provider.AssignRhs();

			InitializeInternalVectors();

			InitializeRhs();

			if (ChildAnalyzer == null)
			{
				throw new InvalidOperationException("Newmark analyzer must contain an embedded analyzer.");
			}

			ChildAnalyzer.Initialize(isFirstAnalysis);
		}

		private void SolveCurrentTimestep()
		{
			Debug.WriteLine("BDF step: {0}", currentTimeStep);

			IGlobalVector rhsVector = provider.GetRhs(currentTimeStep * timeStep);
			solver.LinearSystem.RhsVector = rhsVector;

			if (currentTimeStep + 1 <= BDFOrder)
			{
				BuildMatrices();
			}

			InitializeRhs();
			CalculateRhsImplicit(currentTimeStep * timeStep);

			start = DateTime.Now;
			ChildAnalyzer.Initialize(false);
			ChildAnalyzer.Solve();
			end = DateTime.Now;
			Debug.WriteLine("BDF elapsed time: {0}", end-start);
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

		public void AdvanceStep()
		{
			Debug.WriteLine("Advancing step");

			UpdateSolutionDerivatives();
			UpdateResultStorages(start, end);

			currentTimeStep++;
		}

		/// <summary>
		/// Calculates the right-hand-side of the implicit dyanmic method. This will be used for the solution of the linear dynamic system.
		/// </summary>
		private void CalculateRhsImplicit(double time)
		{
			double a2 = 1 / timeStep;
			int bdfOrderInternal = Math.Min(currentTimeStep + 1, BDFOrder);

			double[] rhsFactors = new double[bdfOrderInternal];
			switch (bdfOrderInternal)
			{
				case 1:
					rhsFactors[0] = 1; // T_n
					break;
				case 2:
					rhsFactors[0] = 4.0 / 2.0; // T_n == n+1
					rhsFactors[1] = -1.0 / 2.0; // T_n-1 == n
					break;
				case 3:
					rhsFactors[0] = 18.0 / 6.0; // T_n == n+2
					rhsFactors[1] = -9.0 / 6.0; // T_n-1 == n+1
					rhsFactors[2] = 2.0 / 6.0; // T_n-2 == n
					break;
				case 4:
					rhsFactors[0] = 48.0 / 12.0; // T_n == n+3
					rhsFactors[1] = -36.0 / 12.0; // T_n-1 == n+2
					rhsFactors[2] = 16.0 / 12.0; // T_n-2 == n+1
					rhsFactors[3] = -3.0 / 12.0; // T_n-3 == n
					break;
				case 5:
					rhsFactors[0] = 300.0 / 60.0; // T_n == n+4
					rhsFactors[1] = -300.0 / 60.0; // T_n-1 == n+3
					rhsFactors[2] = 200.0 / 60.0; // T_n-2 == n+2
					rhsFactors[3] = -75.0 / 60.0; // T_n-4 == n+1
					rhsFactors[4] = 12.0 / 60.0; // T_n-5 == n
					break;
				default: throw new ArgumentException("Wrong BDF Order");
			}

			//var solutionTerm = solver.LinearSystem.Solution.Scale(rhsFactors[0]);
			var solutionTerm = solution.Scale(rhsFactors[0]);
			for (int bdfTerm = 1; bdfTerm < bdfOrderInternal; bdfTerm++)
			{
				solutionTerm.AddIntoThis(solutionOfPreviousStep[bdfTerm - 1].Scale(rhsFactors[bdfTerm]));
			}

			firstOrderDerivativeComponentOfRhs = provider.FirstOrderDerivativeMatrixVectorProduct(solutionTerm);
			firstOrderDerivativeComponentOfRhs.AddIntoThis(provider.FirstOrderDerivativeMatrixVectorProduct(provider.GetFirstOrderDerivativeVectorFromBoundaryConditions(time)));
			firstOrderDerivativeComponentOfRhs.ScaleIntoThis(a2);

			IGlobalVector rhsResult = firstOrderDerivativeComponentOfRhs;
			rhsResult.AddIntoThis(rhs);

			solver.LinearSystem.RhsVector = rhsResult;
		}

		private void InitializeInternalVectors()
		{
			firstOrderDerivativeOfSolutionForRhs = algebraicModel.CreateZeroVector();
			firstOrderDerivativeComponentOfRhs = algebraicModel.CreateZeroVector();
			firstOrderDerivativeOfSolution = algebraicModel.CreateZeroVector();
			rhs = algebraicModel.CreateZeroVector();

			if (solver.LinearSystem.Solution != null)
			{
				solution = solver.LinearSystem.Solution.Copy();
			}
			else
			{
				solution = algebraicModel.CreateZeroVector();
			}

			solutionOfPreviousStep = new IGlobalVector[BDFOrder];
			for (int i = 0; i < BDFOrder; i++)
			{
				solutionOfPreviousStep[i] = algebraicModel.CreateZeroVector();
			}

			solutionOfPreviousStep[0].CopyFrom(solution);
		}

		private void InitializeRhs()
		{
			TransientAnalysisCoefficients coeffs = GetTransientAnalysisCoefficients();
			provider.ProcessRhs(coeffs, solver.LinearSystem.RhsVector);
			rhs.CopyFrom(solver.LinearSystem.RhsVector);
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

		private void UpdateSolutionDerivatives()
		{
			int bdfOrderInternal = Math.Min(currentTimeStep + 1, BDFOrder);

			double[] firstOrderDerivativeFactors = new double[bdfOrderInternal+1];
			switch (bdfOrderInternal)
			{
				case 1:
					firstOrderDerivativeFactors[0] = 1; // T_n+1
					firstOrderDerivativeFactors[1] = -1; // T_n
					break;
				case 2:
					firstOrderDerivativeFactors[0] = 3.0 / 2.0; // T_n+2
					firstOrderDerivativeFactors[1] = -4.0 / 2.0; // T_n+1
					firstOrderDerivativeFactors[2] = 1.0 / 2.0; // T_n
					break;
				case 3:
					firstOrderDerivativeFactors[0] = 11.0 / 6.0; // T_n+3
					firstOrderDerivativeFactors[1] = -18.0 / 6.0; // T_n+2
					firstOrderDerivativeFactors[2] = 9.0 / 6.0; // T_n+1
					firstOrderDerivativeFactors[3] = -2.0 / 6.0; // T_n
					break;
				case 4:
					firstOrderDerivativeFactors[0] = 25.0 / 12.0; // T_n+4
					firstOrderDerivativeFactors[1] = -48.0 / 12.0; // T_n+3
					firstOrderDerivativeFactors[2] = 36.0 / 12.0; // T_n+2
					firstOrderDerivativeFactors[3] = -16.0 / 12.0; // T_n+1
					firstOrderDerivativeFactors[4] = 3.0 / 12.0; // T_n
					break;
				case 5:
					firstOrderDerivativeFactors[0] = 137.0 / 60.0; // T_n+5
					firstOrderDerivativeFactors[1] = -300.0 / 60.0; // T_n+4
					firstOrderDerivativeFactors[2] = 300.0 / 60.0; // T_n+3
					firstOrderDerivativeFactors[3] = -200.0 / 60.0; // T_n+2
					firstOrderDerivativeFactors[4] = 75.0 / 60.0; // T_n+1
					firstOrderDerivativeFactors[5] = -12.0 / 60.0; // T_n
					break;
				default: throw new ArgumentException("Wrong BDF Order");
			}

			firstOrderDerivativeOfSolution.Clear();
			for (int i = 1; i < bdfOrderInternal + 1; i++)
			{
				firstOrderDerivativeOfSolution.AddIntoThis(solutionOfPreviousStep[i - 1].Scale(firstOrderDerivativeFactors[i]));
			}

			firstOrderDerivativeOfSolution.AddIntoThis(solver.LinearSystem.Solution.Scale(firstOrderDerivativeFactors[0]));

			for (int jj = Math.Min(currentTimeStep, BDFOrder - 1); jj >= 1; jj--)
			{
				solutionOfPreviousStep[jj] = solutionOfPreviousStep[jj - 1].Copy();
			}

			solutionOfPreviousStep[0].CopyFrom(solution);
			solution.CopyFrom(solver.LinearSystem.Solution);

		}

		public class Builder
		{
			private readonly double timeStep;
			private readonly double totalTime;
			private readonly IChildAnalyzer childAnalyzer;
			private readonly IModel model;
			private readonly IAlgebraicModel algebraicModel;
			private readonly ISolver solver;
			private readonly ITransientAnalysisProvider provider;
			private readonly int bdfOrder, currentTimeStep;

			public Builder(IModel model, IAlgebraicModel algebraicModel, ISolver solver, ITransientAnalysisProvider provider,
				IChildAnalyzer childAnalyzer, double timeStep, double totalTime, int bdfOrder, int currentTimeStep = 0)
			{
				this.model = model;
				this.algebraicModel = algebraicModel;
				this.solver = solver;
				this.provider = provider;
				this.childAnalyzer = childAnalyzer;

				this.currentTimeStep = currentTimeStep;
				this.timeStep = timeStep;
				this.totalTime = totalTime;
				this.bdfOrder = bdfOrder;
			}

			public BDFDynamicAnalyzer Build()
				=> new BDFDynamicAnalyzer(model, algebraicModel, solver, provider, childAnalyzer, timeStep, totalTime, bdfOrder, currentTimeStep);
		}
	}
}
