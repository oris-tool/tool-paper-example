package org.oristool.examples;

import java.math.BigDecimal;
import java.util.ArrayList;
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
import org.oristool.simulator.stpn.STPNSimulatorComponentsFactory;
import org.oristool.simulator.stpn.TransientMarkingConditionProbability;

public class SoftwareRejuvenation {	
	
    public static void main(String[] args) {
        
    	// Model specification
    	PetriNet pn = new PetriNet();
        Marking m = new Marking();        
        build(pn, m);
        
        // Model evaluation
        analyze(pn, m);
        analyzeWithStopCondition(pn, m);

        // Model simulation
        simulate(pn,m);
    }
    
    public static void simulate(PetriNet net, Marking marking) {

    	Sequencer s = new Sequencer(net, marking, new STPNSimulatorComponentsFactory(), new AnalysisLogger() {
			@Override
			public void log(String message) { }
			@Override
			public void debug(String string) { }
		});
        
        TransientMarkingConditionProbability reward = new TransientMarkingConditionProbability(
        		s, new ContinuousRewardTime(new BigDecimal("0.01")), 100, MarkingCondition.fromString("Down>0 || Detected>0 || Rej>0"));
        
        s.addObserver(reward);
        s.simulate();
        
        TimeSeriesRewardResult result = (TimeSeriesRewardResult)reward.evaluate();    	
    }

    public static void analyzeWithStopCondition(PetriNet net, Marking marking) {

        // Configuration of the regenerative engine for transient analysis with stop condition: time limit=1344, error 0 (per epoch), step size=0.005
        RegTransient analysisWithStopCondition = RegTransient.builder()
        		.greedyPolicy(new BigDecimal("1344"), BigDecimal.ZERO)
                .timeStep(new BigDecimal("0.05"))
                .stopOn(MarkingCondition.fromString("Down>0"))
                .build();

        // Transient analysis of the model: evaluation of transient probabilities of markings
        TransientSolution<DeterministicEnablingState, Marking> solutionWithStopCondition = analysisWithStopCondition.compute(net, marking);

        // Reward property: transient unreliability
        TransientSolution<DeterministicEnablingState, RewardRate> transientUnreliability = 
        		TransientSolution.computeRewards(false, solutionWithStopCondition, "Down>0");
        new TransientSolutionViewer(transientUnreliability);
    }

    public static void analyze(PetriNet net, Marking marking) {
    	
        // Configuration of the regenerative engine for transient analysis: time limit=1344, error 0 (per epoch), step size=0.005
        RegTransient analysis = RegTransient.builder()
        		.greedyPolicy(new BigDecimal("1344"), BigDecimal.ZERO)
                .timeStep(new BigDecimal("0.05"))
                .build();

        // Transient analysis of the model: evaluation of transient probabilities of markings
        TransientSolution<DeterministicEnablingState, Marking> solution = analysis.compute(net, marking);
        
        // Reward property: transient unavailability
        TransientSolution<DeterministicEnablingState, RewardRate> transientUnavailability = 
        		TransientSolution.computeRewards(false, solution, "Down>0 || Detected>0 || Rej>0");
        new TransientSolutionViewer(transientUnavailability);

        // Reward property: cumulative unavailability
        TransientSolution<DeterministicEnablingState, RewardRate> cumulativeUnavailability = 
        		TransientSolution.computeRewards(true, solution, "Down>0 || Detected>0 || Rej>0");
        new TransientSolutionViewer(cumulativeUnavailability);
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

		clock.addFeature(new PostUpdater("Up=0;Down=0;Detected=0", net));
		clock.addFeature(StochasticTransitionFeature.newDeterministicInstance(new BigDecimal("168"), MarkingExpr.from("1", net)));
		clock.addFeature(new Priority(0));
		
		detect.addFeature(new PostUpdater("Wait=0;", net));
		detect.addFeature(StochasticTransitionFeature.newUniformInstance(new BigDecimal("0"), new BigDecimal("4")));
		
		List<GEN> fail_gens = new ArrayList<>();

		DBMZone fail_d_0 = new DBMZone(new Variable("x"));
		Expolynomial fail_e_0 = Expolynomial.fromString("0.0000139");
		fail_d_0.setCoefficient(new Variable("x"), new Variable("t*"), new OmegaBigDecimal("72"));
		fail_d_0.setCoefficient(new Variable("t*"), new Variable("x"), new OmegaBigDecimal("0"));
		GEN fail_gen_0 = new GEN(fail_d_0, fail_e_0);
		fail_gens.add(fail_gen_0);

		DBMZone fail_d_1 = new DBMZone(new Variable("x"));
		Expolynomial fail_e_1 = Expolynomial.fromString("0.0000694");
		fail_d_1.setCoefficient(new Variable("x"), new Variable("t*"), new OmegaBigDecimal("144"));
		fail_d_1.setCoefficient(new Variable("t*"), new Variable("x"), new OmegaBigDecimal("-72"));
		GEN fail_gen_1 = new GEN(fail_d_1, fail_e_1);
		fail_gens.add(fail_gen_1);

		DBMZone fail_d_2 = new DBMZone(new Variable("x"));
		Expolynomial fail_e_2 = Expolynomial.fromString("0.000139");
		fail_d_2.setCoefficient(new Variable("x"), new Variable("t*"), new OmegaBigDecimal("216"));
		fail_d_2.setCoefficient(new Variable("t*"), new Variable("x"), new OmegaBigDecimal("-144"));
		GEN fail_gen_2 = new GEN(fail_d_2, fail_e_2);
		fail_gens.add(fail_gen_2);

		DBMZone fail_d_3 = new DBMZone(new Variable("x"));
		Expolynomial fail_e_3 = Expolynomial.fromString("0.00347 * Exp[-0.00219 x]");
		fail_d_3.setCoefficient(new Variable("x"), new Variable("t*"), OmegaBigDecimal.POSITIVE_INFINITY);
		fail_d_3.setCoefficient(new Variable("t*"), new Variable("x"), new OmegaBigDecimal("-216"));
		GEN fail_gen_3 = new GEN(fail_d_3, fail_e_3);
		fail_gens.add(fail_gen_3);

		PartitionedGEN fail_pFunction = new PartitionedGEN(fail_gens);
		StochasticTransitionFeature fail_feature = StochasticTransitionFeature.of(fail_pFunction);
		fail.addFeature(fail_feature);

		rejuvenate.addFeature(new PostUpdater("Up=1", net));
		rejuvenate.addFeature(StochasticTransitionFeature.newUniformInstance(new BigDecimal("0"), new BigDecimal("2")));
		
		repair.addFeature(new PostUpdater("Wait=1", net));
		repair.addFeature(StochasticTransitionFeature.newUniformInstance(new BigDecimal("4"), new BigDecimal("24")));
		
		// Initial marking
		marking.setTokens(Detected, 0);
		marking.setTokens(Down, 0);
		marking.setTokens(Rej, 0);
		marking.setTokens(Up, 1);
		marking.setTokens(Wait, 1);
	}
}
