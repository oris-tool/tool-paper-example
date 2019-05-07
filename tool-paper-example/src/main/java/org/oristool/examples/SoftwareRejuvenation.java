package org.oristool.examples;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.oristool.analyzer.log.AnalysisLogger;
import org.oristool.math.OmegaBigDecimal;
import org.oristool.math.domain.DBMZone;
import org.oristool.math.expression.Expolynomial;
import org.oristool.math.expression.Variable;
import org.oristool.math.function.GEN;
import org.oristool.math.function.PartitionedGEN;
import org.oristool.models.pn.PostUpdater;
import org.oristool.models.pn.Priority;
import org.oristool.models.stpn.MarkingExpr;
import org.oristool.models.stpn.RewardRate;
import org.oristool.models.stpn.TransientSolution;
import org.oristool.models.stpn.TransientSolutionViewer;
import org.oristool.models.stpn.trans.RegTransient;
import org.oristool.models.stpn.trees.DeterministicEnablingState;
import org.oristool.models.stpn.trees.StochasticTransitionFeature;
import org.oristool.petrinet.Marking;
import org.oristool.petrinet.MarkingCondition;
import org.oristool.petrinet.PetriNet;
import org.oristool.petrinet.Place;
import org.oristool.petrinet.Transition;
import org.oristool.simulator.Sequencer;
import org.oristool.simulator.TimeSeriesRewardResult;
import org.oristool.simulator.rewards.ContinuousRewardTime;
import org.oristool.simulator.rewards.RewardEvaluator;
import org.oristool.simulator.stpn.STPNSimulatorComponentsFactory;
import org.oristool.simulator.stpn.TransientMarkingConditionProbability;

public class SoftwareRejuvenation {	
	
    public static void main(String[] args) {
    	
    	// Parameters
    	BigDecimal timeLimit = new BigDecimal("1344"); 
    	BigDecimal timeStep = new BigDecimal("0.05");
    	String unavailability = "Down>0||Detected>0||Rej>0";
        String unreliability = "Down>0";
        String stop = "Down>0";
        int simulationRunsNumber = 100;
    	
    	// Model specification
    	PetriNet pn = new PetriNet();
        Marking m = new Marking();        
        build(pn, m);
        
        // Model evaluation
        analyze(pn, m, timeLimit, timeStep, unavailability);
        analyzeWithStopCondition(pn, m, timeLimit, timeStep, unreliability, stop);

        // Model simulation
        simulate(pn,m, timeLimit, timeStep, unavailability, simulationRunsNumber);
    }
    
    public static void build(PetriNet net, Marking marking) {

		// Places
		Place Detected = net.addPlace("Detected");
		Place Down = net.addPlace("Down");
		Place Rej = net.addPlace("Rej");
		Place Up = net.addPlace("Up");
		Place Wait = net.addPlace("Wait");
		
		// Transitions
		Transition clock = net.addTransition("clock");
		Transition detect = net.addTransition("detect");
		Transition fail = net.addTransition("fail");
		Transition rejuvenate = net.addTransition("rejuvenate");
		Transition repair = net.addTransition("repair");

		// Preconditions and postconditions
		net.addPostcondition(repair, Up);
		net.addPrecondition(Down, detect);
		net.addPostcondition(clock, Rej);
		net.addPrecondition(Rej, rejuvenate);
		net.addPrecondition(Detected, repair);
		net.addPostcondition(detect, Detected);
		net.addPostcondition(rejuvenate, Wait);
		net.addPrecondition(Up, fail);
		net.addPostcondition(fail, Down);
		net.addPrecondition(Wait, clock);

		// clock is a DET transition with firing time 168 and with update function "Up=0;Down=0;Detected=0"
		clock.addFeature(new PostUpdater("Up=0;Down=0;Detected=0", net));
		clock.addFeature(StochasticTransitionFeature.newDeterministicInstance(new BigDecimal("168"), MarkingExpr.from("1", net)));
		clock.addFeature(new Priority(0)); // Not relevant for the model
		
		// detect is a GEN transition with uniform distribution over [0,4] and with update function "Wait=0;"
		detect.addFeature(new PostUpdater("Wait=0", net));
		detect.addFeature(StochasticTransitionFeature.newUniformInstance(new BigDecimal("0"), new BigDecimal("4")));

		// The failure time is observed to be: 
		// lower than 72 h (3 days) with probability 0.001
		// lower than 144 h (6 days) with probability 0.006
		// lower than 216 h (9 days) with probability 0:016
		// larger than 216 h (9 days) with probability 0.984 and mean value equal to 672 h (28 days)

		// fail is a GEN transition with piece-wise distribution represented over 4 intervals
		List<GEN> fail_gens = new ArrayList<>();

		// 1st interval: uniform distribution over [0,72]
		DBMZone fail_d_0 = new DBMZone(Variable.X);
		Expolynomial fail_e_0 = Expolynomial.fromString("0.0000139"); // constant
		fail_d_0.setCoefficient(Variable.X, Variable.TSTAR, new OmegaBigDecimal("72")); // t - tstar <= 72
		fail_d_0.setCoefficient(Variable.TSTAR, Variable.X, new OmegaBigDecimal("0"));  // t - tstart >= 0 (i.e., tstar - t <= 0)
		GEN fail_gen_0 = new GEN(fail_d_0, fail_e_0);
		fail_gens.add(fail_gen_0);

		// 2nd interval: uniform distribution over [72,144]
		DBMZone fail_d_1 = new DBMZone(Variable.X);
		Expolynomial fail_e_1 = Expolynomial.fromString("0.0000694"); // constant
		fail_d_1.setCoefficient(Variable.X, Variable.TSTAR, new OmegaBigDecimal("144")); // t - tstar <= 144
		fail_d_1.setCoefficient(Variable.TSTAR, Variable.X, new OmegaBigDecimal("-72")); // t - tstar >= 72
		GEN fail_gen_1 = new GEN(fail_d_1, fail_e_1);
		fail_gens.add(fail_gen_1);

		// 3rd interval: uniform distribution over [144,216]
		DBMZone fail_d_2 = new DBMZone(Variable.X);
		Expolynomial fail_e_2 = Expolynomial.fromString("0.000139"); // constant
		fail_d_2.setCoefficient(Variable.X, Variable.TSTAR, new OmegaBigDecimal("216"));  // t - tstar <= 216
		fail_d_2.setCoefficient(Variable.TSTAR, Variable.X, new OmegaBigDecimal("-144")); // t - tstar >= 144
		GEN fail_gen_2 = new GEN(fail_d_2, fail_e_2);
		fail_gens.add(fail_gen_2);

		// 4th interval: shifted exponential distribution over [216,infty)
		DBMZone fail_d_3 = new DBMZone(Variable.X);
		Expolynomial fail_e_3 = Expolynomial.fromString("0.00347 * Exp[-0.00219 x]"); // exponential
		fail_d_3.setCoefficient(Variable.X, Variable.TSTAR, OmegaBigDecimal.POSITIVE_INFINITY); // t - tstar < infty
		fail_d_3.setCoefficient(Variable.TSTAR, Variable.X, new OmegaBigDecimal("-216"));       // t - tstar >= 216
		GEN fail_gen_3 = new GEN(fail_d_3, fail_e_3);
		fail_gens.add(fail_gen_3);

		PartitionedGEN fail_pFunction = new PartitionedGEN(fail_gens);
		StochasticTransitionFeature fail_feature = StochasticTransitionFeature.of(fail_pFunction);
		fail.addFeature(fail_feature);

		// rejuvenate is a GEN transition with uniform distribution over [0,2] and with update function "Up=1"
		rejuvenate.addFeature(new PostUpdater("Up=1", net));
		rejuvenate.addFeature(StochasticTransitionFeature.newUniformInstance(new BigDecimal("0"), new BigDecimal("2")));
		
		// repair is a GEN transition with uniform distribution over [4,24] and with update function "Up=1"
		repair.addFeature(new PostUpdater("Wait=1", net));
		repair.addFeature(StochasticTransitionFeature.newUniformInstance(new BigDecimal("4"), new BigDecimal("24")));
		
		// Initial marking
		marking.setTokens(Detected, 0);
		marking.setTokens(Down, 0);
		marking.setTokens(Rej, 0);
		marking.setTokens(Up, 1);
		marking.setTokens(Wait, 1);
	}
	
    public static void analyze(PetriNet net, Marking marking, BigDecimal timeLimit, BigDecimal timeStep, String rewardString) {
    	
        // Configure the regenerative engine for transient analysis (error per epoch equal to 0)
        RegTransient analysis = RegTransient.builder()
        		.greedyPolicy(timeLimit, BigDecimal.ZERO)
                .timeStep(timeStep)
                .build();

        // Evaluate transient probabilities of markings
        TransientSolution<DeterministicEnablingState, Marking> solution = analysis.compute(net, marking);
        
        // Evaluate transient unavailability
        TransientSolution<DeterministicEnablingState, RewardRate> transientUnavailability = 
        		TransientSolution.computeRewards(false, solution, rewardString);

        // Plot results
        new TransientSolutionViewer(transientUnavailability);

        // Evaluate cumulative unavailability
        TransientSolution<DeterministicEnablingState, RewardRate> cumulativeUnavailability = 
        		TransientSolution.computeRewards(true, solution, rewardString);

        // Plot results
        new TransientSolutionViewer(cumulativeUnavailability);
    }

    public static void analyzeWithStopCondition(PetriNet net, Marking marking, BigDecimal timeLimit, BigDecimal timeStep, String rewardString, String stopString) {

        // Configure the regenerative engine for transient analysis with stop condition (error per epoch equal to 0)
        RegTransient analysisWithStopCondition = RegTransient.builder()
        		.greedyPolicy(timeLimit, BigDecimal.ZERO)
                .timeStep(timeStep)
                .stopOn(MarkingCondition.fromString(stopString))
                .build();

        // Evaluate transient probabilities of markings
        TransientSolution<DeterministicEnablingState, Marking> solution = analysisWithStopCondition.compute(net, marking);

        // Evaluate transient unreliability
        TransientSolution<DeterministicEnablingState, RewardRate> transientUnreliability = 
        		TransientSolution.computeRewards(false, solution, stopString);
        
        // Plot results
        new TransientSolutionViewer(transientUnreliability);
    }
    
    public static void simulate(PetriNet net, Marking initialMarking, BigDecimal timeLimit, BigDecimal timeStep, String rewardString, int runsNumber) {

    	Sequencer s = new Sequencer(net, initialMarking, new STPNSimulatorComponentsFactory(), new AnalysisLogger() {
			@Override
			public void log(String message) { }
			@Override
			public void debug(String string) { }
		});
        
    	// Derive the number of time points
    	int samplesNumber = (timeLimit.divide(timeStep)).intValue() + 1;
    	
    	// Create a reward (which is a sequencer observer)
    	TransientMarkingConditionProbability reward = new TransientMarkingConditionProbability(
        		s, new ContinuousRewardTime(timeStep), samplesNumber, MarkingCondition.fromString(rewardString));
        
    	// Create a reward evaluator (which is a reward observer)
    	RewardEvaluator rewardEvaluator = new RewardEvaluator(reward, runsNumber);
    	
    	// Run simulation
        s.simulate();
        
        // Get simulation results
        TimeSeriesRewardResult result = (TimeSeriesRewardResult)reward.evaluate();
        
        // Plot results
        new TransientSolutionViewer(getTransientSolutionFromSimulatorResult(result, rewardString, initialMarking, timeLimit, timeStep));
    }

    public static TransientSolution<Marking, RewardRate> getTransientSolutionFromSimulatorResult(
    		TimeSeriesRewardResult result, String rewardString, Marking initialMarking, BigDecimal timeLimit, BigDecimal timeStep) {

        RewardRate rewardRate = RewardRate.fromString(rewardString);
    	List<Marking> regenerations = new ArrayList<>(Arrays.asList(initialMarking));
        List<RewardRate> columnStates = new ArrayList<>();
        columnStates.add(rewardRate);
        
        TransientSolution<Marking, RewardRate> solution = new TransientSolution<>(timeLimit, timeStep, regenerations, columnStates, initialMarking);
        
        List<Marking> mrkTmp = new ArrayList<>(result.getMarkings());
        TransientSolution<Marking, Marking> tmpSolution = new TransientSolution<>(timeLimit, timeStep, regenerations, mrkTmp, initialMarking);
        
        for (int t = 0; t < tmpSolution.getSolution().length; t++) {
            for (int i = 0; i < mrkTmp.size(); i++) {
                tmpSolution.getSolution()[t][0][i] = result.getTimeSeries(mrkTmp.get(i))[t].doubleValue();
            }
        }
        
        // Evaluate the reward
        TransientSolution<Marking, RewardRate> rewardTmpResult = TransientSolution.computeRewards(false, tmpSolution, rewardRate);
        for (int t = 0; t < solution.getSolution().length; t++) {
            solution.getSolution()[t][0][columnStates.indexOf(rewardRate)] 
                    = rewardTmpResult.getSolution()[t][0][0];
        }
        
        return solution;
    }
    
}
