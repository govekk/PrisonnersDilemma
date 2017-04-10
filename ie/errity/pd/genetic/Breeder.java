package ie.errity.pd.genetic;

import ie.errity.pd.Prisoner;
import ie.errity.pd.Rules;

import javax.swing.*;
import java.awt.*;
import java.util.BitSet;
import java.util.Random;


/**
 * Provides a means of evolving a population of 
 *{@link  ie.errity.pd.Prisoner Prisoners} 
 *via a genetic algorithm
 * @author	Andrew Errity 99086921
 * @author      Sherri Goings (modified 4/4/2017)
 */
public class Breeder extends JPanel
{
    private Prisoner curPopulation[];
    private double mutateP, crossP; //mutate, cross probabilities
    private int selection, selParam, seed; // selection method, any associated selection parameter, random number seed
    private int popSize;
    private Random rand;
	
    /**
     *Create a new Genetic Breeder
     */
    public Breeder()
    {	
	//defaults
	mutateP = .001;
	crossP = .95;
	selection = 0;
	selParam = 1;
	seed = -1;
	rand = new Random(); //uses current time as seed
    }
	
    /**
     *Create a new Genetic Breeder using {@link  ie.errity.pd.Rules Rules} given
     *@param r1 parameters for the breeder
     */
    public Breeder(Rules r1)
    {	
	//init based on rules	
	mutateP = r1.getMutateP();
	crossP = r1.getCrossP();
	selection = r1.getSelection();
	selParam = r1.getSelectionParameter();
	seed = r1.getSeed();
	if (seed==-1)
	    rand = new Random();
	else
	    rand = new Random(seed); //uses seed set in rules for repeatability
    }

	
    /**
     *Breeds the next generation (panmictic mating) of an array of 
     *{@link  ie.errity.pd.Prisoner Prisoners} 
     *@param c	initial population (raw fitness of population must be calcualted previously)
     *@return the next generation
     */
    public Prisoner[] Breed(Prisoner[] c) {
		curPopulation = c;	//population to breed
		popSize = curPopulation.length;

		Prisoner Selected[] = new Prisoner[popSize]; // parent pop after selection


		// Select parents for next gen
		// ***ADD CODE (make separate functions) to perform each of the types of selection***
		// selection is an int that determines the method to use
		// selParam is an int that can be used as a parameter for a selection method if required

		// One silly selection method which uses rand and the selParam is given below
		// Be sure to use "rand" class variable here as your random number generator (i.e. don't create a new one!)

		// selection method 0, use parameter as a threshold percentage.  If a random number is less than param,
		// select the best individual from the population.  Otherwise select a random individual.
		// In this method a param of 0 gives random selection, and a param of 100 gives best-wins-all selection.
		if (selection == 0) {
			// find index of most fit individual
			double maxFit = 0;
			int indexBest = 0;
			for (int i=0; i<popSize; i++) {
				if (curPopulation[i].getScore() > maxFit) {
					maxFit = curPopulation[i].getScore();
					indexBest = i;
				}
			}

			// select according to description above for this method
			for (int i=0; i<popSize; i++) {
			int selIndex = 0;
			if (rand.nextInt(100) < selParam)  // rand less than param, select best individual
				{
				selIndex = indexBest;
				}
			else  // otherwise select random individual
				{
				selIndex = rand.nextInt(popSize);
				}
			Selected[i] = (Prisoner)curPopulation[selIndex].clone();
			}
		} else if (selection == 1) {
			// if elitism, choose selParam highest fitness individuals to not mutate
			Prisoner[] sortedPrisoners = sortPrisoners(curPopulation);
			if (selParam > 0) {
				for (int i = 0; i < selParam; i++) {
					Selected[i] = new Prisoner(sortedPrisoners[i].getStrat());
				}
				for (int i = selParam; i < popSize; i++) {
					// Currently sets all non-elites to "ALLD"
					Selected[i] = new Prisoner("ALLD");
				}
			}
			else {
				for (int i = 0; i < popSize; i++) {
					Selected[i] = new Prisoner("ALLD");
				}
			}
			repaint();
			return Selected;
		} else {  // any other selection method fill pop with always cooperate
			for (int i=0; i<popSize; i++)
			Selected[i] = new Prisoner("ALLC");
		}

		return Variation(Selected);
	}

	/**
	 * Scales the fitness scores of a population using sigma scaling
	 * @param curPopulation: a population of prisoners
	 * @return the same population of prisoners with their scores scaled
	 */
	private Prisoner[] SigmaScaling(Prisoner[] curPopulation) {
    	// get average fitness
		float avgScore = 0;
		for (int i = 0; i < popSize; i++) {
			avgScore += curPopulation[i].getScore();
		}
		avgScore = avgScore / popSize;

		// get std dev
		float sumSquares = 0;
		for (int i = 0; i < popSize; i++) {
			sumSquares += (curPopulation[i].getScore() - avgScore) * (curPopulation[i].getScore() - avgScore);
		}
		float stdDev = sumSquares/popSize;

		// set scores to scaled fitness
		for (int i = 0; i < popSize; i++) {
			int scaledFit = Math.round(curPopulation[i].getScore()/(2*stdDev));
			curPopulation[i].setScore(scaledFit);
		}
		return curPopulation;
	}

	/**
	 * Selects next generation using fitness proportional selection (roulette wheel)
	 * @param curPopulation: the current population from which to select parents
	 * @param toSelect: the number of prisoners to select for the next generation
	 * @return the next generation
	 */
	private Prisoner[] FitPropSelect(Prisoner[] curPopulation, int toSelect){
		return curPopulation;
	}

	/**
	 * Performs mutation and crossover on a given population of prisoners
	 * @param Selected: the population to be varied
	 * @return varied population of same size
	 */
	private Prisoner[] Variation(Prisoner[] Selected) {
		//Crossover & Mutate each pair of selected parents
		BitSet Offspring[] = new BitSet[2];  // temporarily holds 2 children during crossover/mutation
		for (int d=0; d<popSize; d+=2) {
			// in case of odd population, just mutate and replace last individual
			if (d+1 >= popSize) {
				Offspring[0] = Genetic.mutate(Selected[d].getStrat(), mutateP, rand);
				Selected[d] = new Prisoner(Offspring[0]);
			}
			else {

				if(rand.nextDouble() <= crossP) //Cross Over
					Offspring = Genetic.crossover(Selected[d].getStrat(),Selected[d+1].getStrat(), rand);
				else //clones
				{
					Offspring[0] = (BitSet)Selected[d].getStrat().clone();
					Offspring[1] = (BitSet)Selected[d+1].getStrat().clone();
				}

				//Mutation
				Offspring[0] = Genetic.mutate(Offspring[0],mutateP, rand);
				Offspring[1] = Genetic.mutate(Offspring[1],mutateP, rand);

				//Replacement - we are done with parents d & d+1, so just replace with children without
				// creating an entire new array
				Selected[d] = new Prisoner(Offspring[0]);
				Selected[d+1] = new Prisoner(Offspring[1]);
			}
		}
		// pass on children pop to be parents of next gen
		curPopulation = Selected;
		repaint();	//update display (if any)
		return curPopulation; //return the bred population
	}

    /**
     *Responsible for updating the graphical representation
     */
    public void paintComponent(Graphics g) 
    {
        super.paintComponent(g); //paint background

	//Get limits of viewable area
      	Insets insets = getInsets();
      	int x0 = insets.left;
      	int y0 = insets.top;
	int currentWidth = getWidth() - insets.left - insets.right;
	int currentHeight = getHeight() - insets.top - insets.bottom;
	    
	//Display a series of rectangles, representing the players
	for(int i = 0; i < popSize; i++)
	    {
	    	g.setColor(curPopulation[i].getColor());	
	    	g.fillRect((x0*2)+((currentWidth/popSize)*(i)),(currentHeight/4)+y0,(currentWidth/popSize),currentHeight/2);
	    }
    }
	
    /**
     *Set the {@link  ie.errity.pd.Rules Rules}
     *@param r new {@link  ie.errity.pd.Rules Rules}
     */
    public void setRules(Rules r) {
        mutateP = r.getMutateP();
        crossP = r.getCrossP();
        selection = r.getSelection();
        selParam = r.getSelectionParameter();
        seed = r.getSeed();
        if (seed > -1)
            rand.setSeed(seed);
    }

    /**
     * Sorts an array of prisoners by fitness using mergesort
     * @param prisoners (an array of prisoners)
     * @return array of prisoners sorted in descending order of fitness
     */
    private Prisoner[] sortPrisoners(Prisoner[] prisoners){
        int numPrisoners = prisoners.length;
        if (numPrisoners <= 1) return prisoners;
        else {
            int leftSize;
            int rightSize = numPrisoners/2;
            if (numPrisoners % 2 == 1) {
                leftSize = rightSize + 1;
            } else {
                leftSize = rightSize;
            }
            Prisoner[] leftPrisoners = new Prisoner[leftSize];
            Prisoner[] rightPrisoners = new Prisoner[rightSize];
            for (int i = 0; i < leftSize; i++) {
                leftPrisoners[i] = prisoners[i];
            }
            for (int i = 0; i < rightSize; i++) {
                rightPrisoners[i] = prisoners[i+leftSize];
            }
            leftPrisoners = sortPrisoners(leftPrisoners);
            rightPrisoners = sortPrisoners(rightPrisoners);
            int j = 0;
            Prisoner[] sortedPrisoners = new Prisoner[popSize];
            for (int i = 0; i < leftSize; i++) {
                while ((j < rightSize) && leftPrisoners[i].getScore() < rightPrisoners[j].getScore()) {
                    sortedPrisoners[i+j] = rightPrisoners[j];
                    j++;
                }
                sortedPrisoners[i+j] = leftPrisoners[i];
            }
            for (int k = j; k < rightSize; k++) {
            	sortedPrisoners[leftSize + k] = rightPrisoners[k];
			}
            return sortedPrisoners;
        }
    }

    /**
     *Reset the breeder
     */
    public void clear()
    {
	popSize = 0;
	repaint();
    }
}
	
