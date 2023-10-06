using System;
using System.Diagnostics;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.NonLinear;
using System.Collections.Generic;

namespace MGroup.NumericalAnalyzers.Discretization.NonLinear
{
	public class LoadControlAnalyzer2Tolerance : NonLinearAnalyzerBase
	{
		private bool stopIfNotConverged = false;

		//public List<IGlobalVector> SolverRHSs { get; set; } = new List<IGlobalVector>();
		//public List<IGlobalVector> SolverSolutions { get; set; } = new List<IGlobalVector>();

		//private int minimumIters = 10;
		//private int maximumIters = 100;
		private IGlobalVector uPlusduPrevious;

		/// <summary>
		/// This class solves the linearized geoemtrically nonlinear system of equations according to Newton-Raphson's load control incremental-iterative method.
		/// </summary>
		/// <param name="model">Instance of the model that will be solved</param>
		/// <param name="solver">Instance of the solver that will solve the linear system of equations</param>
		/// <param name="provider">Instance of the problem type to be solved</param>
		/// <param name="subdomainUpdaters">Instance that updates constraints, right-hand-side vector, updates and resets state</param>
		/// <param name="numIncrements">Number of total load increments</param>
		/// <param name="maxIterationsPerIncrement">Number of maximum iterations within a load increment</param>
		/// <param name="numIterationsForMatrixRebuild">Number of iterations for the rebuild of the siffness matrix within a load increment</param>
		/// <param name="residualTolerance">Tolerance for the convergence criterion of the residual forces</param>
		private LoadControlAnalyzer2Tolerance(IAlgebraicModel algebraicModel, ISolver solver, INonLinearProvider provider,
			int numIncrements, int maxIterationsPerIncrement, int numIterationsForMatrixRebuild, double residualTolerance, bool stopIfNotConverged)
			: base(algebraicModel, solver, provider, numIncrements, maxIterationsPerIncrement,
				numIterationsForMatrixRebuild, residualTolerance)
		{
			this.stopIfNotConverged = stopIfNotConverged;
			analysisStatistics.AlgorithmName = "Load control analyzer";
		}


		/// <summary>
		/// Solves the nonlinear equations and calculates the displacements vector.
		/// </summary>
		public override void Solve()
		{
			bool notConverged = false;

			analysisStatistics.NumIterationsRequired = 0;
			InitializeLogs();

			DateTime start = DateTime.Now;
			UpdateInternalVectors();
			for (int increment = 0; increment < numIncrements; increment++)
			{
				double errorNorm = 0;
				ClearIncrementalSolutionVector();
				UpdateRhs(increment);
				UpdateMatrices();

				double firstError = 0;
				int iteration = 0;
				for (iteration = 0; iteration < maxIterationsPerIncrement; iteration++)
				{
					analysisStatistics.NumIterationsRequired++;

					if (iteration == maxIterationsPerIncrement - 1)
					{
						notConverged = true;
						analysisStatistics.ResidualNormRatioEstimation = errorNorm;
						if (stopIfNotConverged)
						{
							return;
						}
						else
						{
							break;
						}
					}

					if (double.IsNaN(errorNorm))
					{
						notConverged = true;
						analysisStatistics.ResidualNormRatioEstimation = errorNorm;
						if (stopIfNotConverged)
						{
							return;
						}
						else
						{
							break;
						}
					}

					solver.Solve();

					////outuput
					//SolverRHSs.Add(solver.LinearSystem.RhsVector.Copy());
					//SolverSolutions.Add(solver.LinearSystem.Solution.Copy());

					IGlobalVector internalRhsVector = CalculateInternalRhs(increment, iteration);
					double residualNormCurrent = UpdateResidualForcesAndNorm(increment, iteration, internalRhsVector);
					errorNorm = globalRhsNormInitial != 0 ? residualNormCurrent / globalRhsNormInitial : 0;

					if (iteration == 0)
					{
						firstError = errorNorm;
					}

					if (IncrementalDisplacementsLog != null)
					{
						IncrementalDisplacementsLog.StoreDisplacements(uPlusdu);
					}

					//if (iteration < minimumIters) { errorNorm = 1; }
					if (iteration == 0)
					{
						errorNorm = 1;
						uPlusduPrevious= uPlusdu.Copy();
					}
					else
					{
						var correction = uPlusdu.Subtract(uPlusduPrevious);
						errorNorm = correction.Norm2()/uPlusdu.Norm2();
						uPlusduPrevious=uPlusdu.Copy();
					}

					if (errorNorm < residualTolerance)
					{
						if (analysisStatistics.ResidualNormRatioEstimation < errorNorm)
						{
							analysisStatistics.ResidualNormRatioEstimation = errorNorm;
						}

						if (IncrementalLog != null)
						{
							IncrementalLog.LogTotalDataForIncrement(increment, iteration, errorNorm, uPlusdu, internalRhsVector);
						}

						break;
					}

					if ((iteration + 1) % numIterationsForMatrixRebuild == 0)
					{
						provider.Reset();
						parentAnalyzer.BuildMatrices();
					}
				}

				if (TotalDisplacementsPerIterationLog != null)
				{
					TotalDisplacementsPerIterationLog.StoreDisplacements(uPlusdu);
				}

				Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", iteration, firstError, errorNorm);
				SaveMaterialStateAndUpdateSolution();
			}

			analysisStatistics.HasConverged = !notConverged;
			DateTime end = DateTime.Now;
			StoreLogResults(start, end);
		}

		private void UpdateMatrices()
		{
			provider.Reset();
			parentAnalyzer.BuildMatrices();
		}

		protected double UpdateResidualForcesAndNorm(int currentIncrement, int iteration, IGlobalVector internalRhs)
		{
			solver.LinearSystem.RhsVector.CopyFrom(lastRhs);
			solver.LinearSystem.RhsVector.SubtractIntoThis(internalRhs);
			return provider.CalculateRhsNorm(solver.LinearSystem.RhsVector);
		}

		private new void UpdateRhs(int step)
		{
			lastRhs.AddIntoThis(rhsIncrement);
			

			solver.LinearSystem.RhsVector.CopyFrom(lastRhs);

			// edw upologizetai to state kai lamavametai upopsin kai i dinamiki suneisfora logw allagis tou (u_(t-1) ara kai tou u_dot(t) kai tou u_ddot(t))
			//alla afta ginontai meta to update tou lastRhs pou idi exei ginei kai den epireazetai apo afta
			IGlobalVector internalRhs = provider.CalculateResponseIntegralVector(u);
			provider.ProcessInternalRhs(u, internalRhs);
			if (parentAnalyzer != null)
			{
				IGlobalVector otherRhsComponents = parentAnalyzer.GetOtherRhsComponents(u);
				internalRhs.AddIntoThis(otherRhsComponents);
			}

			solver.LinearSystem.RhsVector.SubtractIntoThis(internalRhs);
		}

		public class Builder : NonLinearAnalyzerBuilderBase
		{
			private bool stopIfNotConverged = true;

			public Builder(IAlgebraicModel algebraicModel, ISolver solver, INonLinearProvider provider, int numIncrements, bool stopIfNotConverged = true)
				: base(algebraicModel, solver, provider, numIncrements)
			{
				MaxIterationsPerIncrement = 1000;
				NumIterationsForMatrixRebuild = 1;
				ResidualTolerance = 1E-3;
				this.stopIfNotConverged = stopIfNotConverged;
			}

			public LoadControlAnalyzer2Tolerance Build() => new LoadControlAnalyzer2Tolerance(algebraicModel, solver, provider,
				numIncrements, maxIterationsPerIncrement, numIterationsForMatrixRebuild, residualTolerance, stopIfNotConverged);
		}
	}
}
