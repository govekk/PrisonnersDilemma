package ie.errity.pd;

import java.util.BitSet;
import java.util.Random;

/**
 *Iterated Prisoner's Dilemma Tournament
 *<BR>In such a tournament every player plays the IPD against every other player
 *and themselves.
 *@author Andrew Errity 99086921
 * @author      Sherri Goings (modified 4/4/2017)
 */
public class Tournament
{
    private Prisoner Players[];
    private int numPlayers;
    private int Results[];
    private Rules rules;
    boolean done;
	
    /**
     *Create new Iterated Prisoner's Dilemma Tournament
     *@param plyrs an array containing the competitors
     *@param r the rules to govern the Tournament
     */
    public Tournament(Prisoner plyrs[], Rules r)
    {
	Players = plyrs;
	rules = r;
	numPlayers = Players.length;
	Results = new int[numPlayers];
	for(int i = 0; i < numPlayers; i++)
	    Results[i] = 0;
	done = false;
    }
	
    /**
     *Play the tournament
     *@return the Tournament results, an array of payoffs corresponding to the 
     *player's scores
     */
    public void Play()
    {
	int i;
	Game g;
	Prisoner opponent;
		
	// each individual plays 1 or more games and adds the score from each to their total fitness
	for(i = 0; i < numPlayers; i++)
	    {
		// get the user setting for evaluation method and play the appropriate opponents, current options are
		// 0 - play everyone including yourself
		// 1 - play once against always defect
		// 2 - play once against always cooperate
		// 3 - play once against tit for tat
		// 4 - play once against tit for two tat
		// 5 - play one game against each of the 4 fixed strategies
		// 6 - play 5 games where for each game the opponent is chosen randomly from the 4 fixed strategies
		int evalMethod = rules.getEvaluationMethod();

		//Every player plays themselves and every other player	
		if (evalMethod == 0) {
		    
		    // set opponent to cur player i, then play against every other player in pop, updating their scores
		    // use clone so cur opponent's score is not updated, only each of the other players in the pop
		    opponent = (Prisoner)Players[i].clone();
		    for(int j = 0; j < numPlayers; j++) {
			g = new Game(opponent,Players[j],rules);
			g.Play();
		    }
		}

		// play 1 game of each player against a fixed player that always defects 
		else if (evalMethod == 1) {
		    opponent = new Prisoner("ALLD");
		    g = new Game(opponent,Players[i],rules);
		    g.Play();
		}

		// play 1 game of each player against a fixed player that always cooperates
		else if (evalMethod == 2) {
		    opponent = new Prisoner("ALLC");
		    g = new Game(opponent,Players[i],rules);
		    g.Play();
		}

		// play 1 game of each player against a fixed player that uses tit-for-tat
		else if (evalMethod == 3) {
		    opponent = new Prisoner("TFT");
		    g = new Game(opponent,Players[i],rules);
		    g.Play();
		}

		// play 1 game of each player against a fixed player that uses tit-for-two-tat
		else if (evalMethod == 4) {
		    opponent = new Prisoner("TFTT");
		    g = new Game(opponent,Players[i],rules);
		    g.Play();
		}
		
		// play 1 game of each player against each of the 4 fixed player strategies
		else if (evalMethod == 5) {
		    opponent = new Prisoner("TFTT");
		    g = new Game(opponent,Players[i],rules);
		    g.Play();
		    opponent = new Prisoner("TFT");
		    g = new Game(opponent,Players[i],rules);
		    g.Play();
		    opponent = new Prisoner("ALLD");
		    g = new Game(opponent,Players[i],rules);
		    g.Play();
		    opponent = new Prisoner("ALLC");
		    g = new Game(opponent,Players[i],rules);
		    g.Play();
		}

		// play 5 games, each against a random strategy from the 4 fixed strategies
		else if (evalMethod == 6) {
		    Random rand = new Random(rules.getSeed());
		    
		    for (int j=0; j<5; j++) {
			int strategy = rand.nextInt(4);

			if (strategy == 0)
			    opponent = new Prisoner("TFTT");
			else if (strategy == 1)
			    opponent = new Prisoner("TFT");
			else if (strategy == 2)
			    opponent = new Prisoner("ALLD");
			else
			    opponent = new Prisoner("ALLC");

			g = new Game(opponent,Players[i],rules);
			g.Play();
		    }
		}

		// silly workaround to fact that when scores are displayed it is assumed that pop_size games were played
		// so the total score is divided by pop_size, in this case only 1 game was played so multiply that 1
		// score by pop_size, as opposed to playing pop_size games when the outcome of each will be identical
		if (evalMethod >= 1 && evalMethod <= 4)
		    Players[i].setScore(Players[i].getScore()*numPlayers);
		else if (evalMethod == 5) // played 4 games
		    Players[i].setScore(Players[i].getScore()*(int)(numPlayers/4.0));
		else if (evalMethod == 6) // played 5 games
		    Players[i].setScore(Players[i].getScore()*(int)(numPlayers/5.0));

		
	    }
	done = true;
		
    }
	
    /**
     *Returns the Tournament results
     *@return the Tournament results
     */
    public int[] getResults(){return Results;}
	
	
    /**
     *Returns index of lowest scoring individual
     *@return index of lowest scoring individual
     */
    public int minResult() 
    {
	int index = -1;
	if(done)
	    {
		int min = Results[0];
		index = 0;
		for(int i = 1; i < numPlayers; i++)
		    if(Results[i] < min)
			{
			    min = Results[i];
			    index = i;
			}
	    }
	return index;
    }
	
    /**
     *Returns index of highest scoring individual
     *@return index of highest scoring individual
     */
    public int maxResult()
    {
	int index = -1;
	if(done)
	    {
		int max = Results[0];
		index = 0;
		for(int i = 1; i < numPlayers; i++)
		    if(Results[i] > max)
			{
			    max = Results[i];
			    index = i;
			}
	    }
	return index;
    }
	
    /**
     *Returns the average score
     *@return the average score
     */
    public double meanResult()
    {
	double mean = 0;
	if(done)
	    {
		int total = 0;
		for(int i = 0; i < numPlayers; i++)
		    total += Results[i];
		mean = total/numPlayers;
	    }	
	return mean;
    }
	
	
}
