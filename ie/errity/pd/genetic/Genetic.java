package ie.errity.pd.genetic;

import ie.errity.pd.*;

import java.util.BitSet;
import java.util.Random;

/**
 *Provides Genetic operations
 *@author Andrew Errity 99086921
 * @author      Sherri Goings (modified 4/4/2017) to take rand as parameter so can repeat identical runs
 */
public class Genetic
{	

    /**
     *Mate two parents using random, one point crossover
     *@param parenta first parent (<code>BitSet</code> representation)
     *@param parentb second parent (<code>BitSet</code> representation)
     *@return an array containing the two children (<code>BitSet</code> 
     *representation)
     */
    public static BitSet[] crossover(BitSet parenta,BitSet parentb, Random rand)
    {
	BitSet child1 = new BitSet(71);
	BitSet child2 = new BitSet(71);
		
	//One point splicing
	int slicePoint = rand.nextInt(71); //rnd num between 0-70
	BitSet a = (BitSet)parenta.clone();
	a.clear(slicePoint,71);
	BitSet b = (BitSet)parenta.clone();
	b.clear(0,slicePoint);
	BitSet c = (BitSet)parentb.clone();
	c.clear(slicePoint,71);
	BitSet d = (BitSet)parentb.clone();
	d.clear(0,slicePoint);
		
	//Combine start of p1 with end of p2
	child1.or(a);
	child1.or(d);
	//Combine start of p2 with end of p1		
	child2.or(c);
	child2.or(b);
		
	//Return the children
	BitSet[] offspring = {child1, child2};
	return offspring;
    }
	
    /**
     *Mutate (Flip a bit in the bitset) with probability mProb
     *@param original the entity to be mutated
     *@param mProb the probability of a bit being mutated
     *@return the (possibly) mutated entity
     */
    public static BitSet mutate(BitSet original, double mProb, Random rand)
    {
	for(int m = 0; m < 71; m++)
	    {
		//Small possibility a bit copied from parent to child is mutated
		double test = rand.nextDouble();
	
		if(test <= mProb) {
		    original.flip(m);
		}
	    }
	//Return the (possibly) strategy
	return original;
    }
	
	
}
