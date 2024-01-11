import java.io.*;
import java.util.*;

public class Main {
    // meta variable
    static String fileName = "Capital_Cities10.txt";

    // CHANGE CITY AND PRIZEGOAL
    static String begin = "";
    static String end = "";
    static double budget = 0; // budget in miles
    static int n = 48;
    private static double remainingBudget;

    // static variables to be tweaked by user
    static final int TRIALS = 10;
    static final int NUM_AGENTS = 5;
    static final double W = 100.0; // constant value to update the reward table
    static double alpha = 0.125; // .125 learning rate
    static double gamma = 0.35; // .35 discount factor
    static final double delta = 1; // power for Q value
    static final double beta = 2; // power for distance
    static double q0 = 0.8; // coefficient for exploration and exploitation

    // **************************** greedy epsilon
    static double explor_rate = 1.0;  // default 1.0
    static double max_explor = 1.0;  // default 1.0
    static double min_explore = 0.01;  // default 0.01
    static double decay_rate = 0.01;  // default 0.01
    // *************************************

    // flags for graph (do not touch)
    static final int UNVISITED = 0;
    static final int VISITED = 1;
    static final int LAST_VISIT = 2;

    // pre-initialization parameters (do not touch)
    static LinkedList<CityNode> arrCities;
    static ArrayList<String> nameList;
    static ArrayList<String> allNameList;
    static Graph sGraph;
    static double[][] Q;
    static double[][] R;
    static int statesCt;
    static int total_prize = 0;
    static double total_wt = 0;
    static ArrayList<Integer> route;
    static ArrayList<Integer> routeMax;// keep track of best case
    static int prizeMax = Integer.MIN_VALUE;// keep track of best case
    static int episodeMax = 0;// keep track of episode of best case
    // Create a HashMap with Integer keys and String array values
    static HashMap<Integer, ArrayList<Integer>> listMax = new HashMap<>();

    // index of cities
    static int endCity;
    static int beginCity;

    // variables for extra credits
    static long randomSeed = 12345; // random seed for inaccessible edge, can change to any long
    static double missingProb = 0.0; // probability of an inaccessible edge

    public static void main(String[] args) throws IOException {
        // variables for time tracking
        long startTime, endTime;
        double totalTime;

        askForUserInputs();

        System.out.println("\n========== BC-PC-TSP MARL algorithm ==========");
        // initialize list
        initList();
        // initialize Graph
        initGraph();
        // initialize Tables
        initStatics();


        startTime = System.nanoTime();

        printTable("Q_Table_Before",Q);
        printTable("R_Table_Before",R);

        // MARL algorithm
        learnQ();
        traverseQ();

        printTable("Q_Table_After",Q);
        printTable("R_Table_After",R);
        endTime = System.nanoTime();


        totalTime = (double) (endTime - startTime) / 1000000;

//        System.out.println("Algorithm took " + totalTime + "ms to process.");

        System.out.println();
        System.out.println("Value of P^m: " + (prizeMax - arrCities.get(0).pop));
        System.out.println("R^m Route: " + makeRouteString(true));
        System.out.println("Episode P^m and R^m: " + episodeMax);

        // Display Q_table
        for (int i = 0; i < Q.length; i++) {
            for (int j = 0; j < Q[i].length; j++) {
                System.out.print(Q[i][j] + " ");
            }
            System.out.println(); // Move to the next line after printing each row
        }
    }

    public static void printTable(String text, double[][] array){
        // CSV file path
        String csvFilePath = text + "_output.csv";

        // Write the 2D array to the CSV file
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(csvFilePath))) {
            for (int i = 0; i < array.length; i++) {
                for (int j = 0; j < array[i].length; j++) {
                    writer.write(String.valueOf(array[i][j]));

                    // Add a comma if it's not the last element in the row
                    if (j < array[i].length - 1) {
                        writer.write(",");
                    }
                }
                // Move to the next line after each row
                writer.newLine();
            }

//            System.out.println("2D array has been written to " + csvFilePath);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private static String routeString(ArrayList<Integer> routeArray) {
        String ret = "";

        for (int i = 0; i <routeArray.size();i++){
            ret+=arrCities.get(routeArray.get(i)+1).name;
            if(i!=routeArray.size()-1){
                ret+="("+arrCities.get(routeArray.get(i)+1).originalIndex+")";
            }else{
                ret+="("+arrCities.get(routeArray.get(i)+1).originalIndex+")";
            }

            if(i!=routeArray.size()-1){
                ret+=", ";
            }

        }
        return ret;
    }
    private static String makeRouteString(boolean value) {
        String ret = "";
        for (int i = 0; i <routeMax.size();i++){
            ret+=arrCities.get(routeMax.get(i)).name;
            if(i!=routeMax.size()-1){
                ret+="("+arrCities.get(routeMax.get(i)).originalIndex+")";
            }else{
                ret+="("+arrCities.get(routeMax.get(i)).originalIndex+")";
            }
            
            if(i!=routeMax.size()-1){
                ret+=", ";
            }
            
        }
        return ret;
    }



    static void askForUserInputs() {
        Scanner scanner = new Scanner(System.in);
//        System.out.print("Enter the start city: ");
        begin = "CarsonCity,NV";
//        begin = scanner.nextLine();
//        System.out.print("Enter the end city: ");
        end = "CarsonCity,NV";
//        end = scanner.nextLine();
//        System.out.print("Enter the budget in miles: ");
        budget = 6000;
//        budget = scanner.nextInt();
        System.out.println("Starting and Ending in " + begin + " with Range: "+budget);
    }

    /*
     * Documentation for initList()
     * (1) attempts to read in file (throws error if file not found)
     * (2) converts text file info into cityNode object and places into array for
     * use in creating graph
     * (3) optional scanner functionality for user inputted begin and end points
     * (4) restructure arrayList to make BEGIN city the first node, and END city the
     * last node
     */
    static void initList() {
        // (1)
        File towns = new File("src/" + fileName);
        arrCities = new LinkedList<>();
        nameList = new ArrayList<>();
        allNameList = new ArrayList<>();

        try {
            Scanner scan = new Scanner(towns);
            // (2)
            int originIndex=0;
            while (scan.hasNextLine()) {
                String name = scan.next();
                double lat = scan.nextDouble();
                double lon = scan.nextDouble();
                int pop = scan.nextInt();

                arrCities.add(new CityNode(name, lat, lon, pop));
                nameList.add(name.toLowerCase());
                allNameList.add(name.toLowerCase());

                arrCities.get(arrCities.size()-1).originalIndex=originIndex;
                originIndex++;
            }
            scan.close();
        }

        catch (FileNotFoundException p) {
            System.out.println("FILE NOT FOUND.");
            System.exit(0);
        }

        // (3)
        String startCity, endCity;
        startCity = begin.toLowerCase();
        endCity = end.toLowerCase();

        // (4)
        if (nameList.contains(startCity) && nameList.contains(endCity)) {
            if (startCity.equalsIgnoreCase(endCity)) {
                int sIndex = nameList.indexOf(startCity);
                CityNode t1 = new CityNode(arrCities.get(sIndex));
                t1.originalIndex = arrCities.get(sIndex).originalIndex;
                arrCities.remove(sIndex);
                nameList.remove(sIndex);

                arrCities.add(0, t1);
                arrCities.add(t1);
            } else {
                int sIndex = nameList.indexOf(startCity);
                CityNode t1 = new CityNode(arrCities.get(sIndex));
                t1.originalIndex = arrCities.get(sIndex).originalIndex;
                arrCities.remove(sIndex);
                nameList.remove(sIndex);

                int fIndex = nameList.indexOf(endCity);
                CityNode t2 = new CityNode(arrCities.get(fIndex));
                t2.originalIndex = arrCities.get(fIndex).originalIndex;
                arrCities.remove(fIndex);

                arrCities.add(0, t1);
                arrCities.add(t2);
            }
        }

        else {
            System.out.println("Cannot find city.");
            System.exit(0);
        }
    }

    /*
     * Documentation for initGraph()
     * This function simply initializes the graph with requisite flags,
     * prize and weight values for each node,
     * and a name for each node (debug purposes only)
     */
    static void initGraph() {
        // we can have different keys for difference cases
        Random rand = new Random(randomSeed);
        sGraph = new Graph();
        int size = arrCities.size();
        sGraph.Init(size);

        for (int i = 0; i < size; i++) {
            sGraph.setMark(i, UNVISITED);
            sGraph.setName(i, arrCities.get(i).name);
            sGraph.setPrize(i, arrCities.get(i).pop);
            for (int j = 0; j < i + 1; j++) {
                if (i == j) {
                    sGraph.setEdge(i, j, 0);
                } else {
                    // randomly mark some path as inaccessible
                    if (rand.nextDouble() > missingProb) {
                        sGraph.setEdge(i, j, CityNode.getDistance(arrCities.get(i), arrCities.get(j)));
                        sGraph.setEdge(j, i, CityNode.getDistance(arrCities.get(j), arrCities.get(i)));
                    } else {
                        sGraph.setEdge(i, j, Double.MAX_VALUE);
                        sGraph.setEdge(j, i, Double.MAX_VALUE);
                    }
                }
            }
        }
        sGraph.setMark(sGraph.getLastNode(), LAST_VISIT);
        sGraph.constructShortestPath();
    }

    /*
     * Documentation for initStatics()
     * This function simply initializes the 2d arrays used for Q-Learning
     * Reward of edge (i, j), initially -w(i,j)/pv
     */
    static void initStatics() {
        statesCt = sGraph.n()-1;
        Q = new double[statesCt][statesCt];
        R = new double[statesCt][statesCt];

        for (int i = 0; i < statesCt; i++)
            for (int j = 0; j < statesCt; j++) {
                if(i==j){
                    R[i][j] = -999999;
                    Q[i][j] = -9999999;
                } else {
                    R[i][j] = sGraph.weight(i, j) / sGraph.getPrize(j) * -1;
                Q[i][j] = (sGraph.getPrize(i) + sGraph.getPrize(j)) / sGraph.weight(i, j);
//                    Q[i][j] =  sGraph.getPrize(j) / sGraph.weight(i, j);
//                Q[i][j] = 0;
                }

            }
    }

    /*
     * Documentation for learnQ()
     * learnQ() is the primary function for the program
     *
     * (1) All the m agents are initially located at the starting node s with zero
     * collected prizes.
     * [line 236-242]
     *
     * (2) Then each independently follows the action rule to move to the next node
     * to collect prizes and
     * collaboratively updates the Q-value.
     * [line 252-262]
     *
     * (3) When an agent can no longer find a feasible unvisited node to move to due
     * to its insufficient budget,
     * it terminates and goes to t.
     * [line 246]
     *
     * (4) Otherwise, it moves to the next node and collects the prize and continues
     * the prize-collecting process.
     * [line 248-265]
     *
     * (5) Then, it finds among the m routes the one with the maximum collected
     * prizes and updates the
     * reward value and Q-value of the edges that belong to this route.
     * [line 269-281]
     */
    static void learnQ() {
        for (int i = 0; i < TRIALS; i++) {
            // for every episode, change learning rate and discount factor and epsilon?
            //setHypers(i);
            System.out.println("************** Episode: " + i + "***************");
            Agent[] aList = new Agent[NUM_AGENTS];
            for (int j = 0; j < NUM_AGENTS; j++) {
                Graph newGraph = new Graph(sGraph);
                Agent a = new Agent(statesCt, budget, newGraph);
                aList[j] = a;
            }

            while (!allQuotaMet(aList)) {
                for (int j = 0; j < NUM_AGENTS; j++) {
                    Agent aj = aList[j];

                    if (!aj.isDone) {

                        int nextState = getNextStateFromCurState(aj, aj.curState, 1 - q0 * (TRIALS - i) / TRIALS);
//                        System.out.println(nextState);
                        if (nextState == aj.getLastNode()) {
                            aj.isDone = true;
                        }
//                        double maxQ = maxQ(aj, nextState);
//                        double maxQ = findFesibleMaxQ(nextState);
                        // System.out.println("current state    "+ aj.curState);
                        // updating Q-Table
//                        Q[aj.curState][nextState] = (1 - alpha) * Q[aj.curState][nextState] + alpha * gamma * maxQ;
                        // updating current agent
                        aj.indexPath.add(nextState);
//                        System.out.println("current state: " + aj.curState + "   Next State:  "+ nextState);
                        aj.total_wt += aj.shortestPath(aj.curState, nextState);
                        aj.total_prize += aj.getTotalPrize(nextState);
                        aj.setAgentMark(nextState, VISITED);
                        aj.curState = nextState;
                    }
                }
            }

            int mostFitIndex = findHighestPrize(aList);

            Agent jStar = aList[mostFitIndex];
            ArrayList<Integer> path = jStar.indexPath;
            jStar.resetAgentMarks();

            int lastDestination = path.getLast();
            path.add(0,lastDestination);
//            System.out.println(lastDestination);
            for (int v = 0; v < path.size()-1; v++) {
//                double q = Q[path.get(v)][path.get(v + 1)];
//                double maxQ = maxQ(jStar, path.get(v + 1));
                double maxQ = findMaxQ(path.get(v+1)-1);
//                System.out.println("-----------updated ---------------");
                //R-Table updated
//                R[path.get(v)][path.get(v + 1)] += (2*W / jStar.total_prize);
                // NO R-Table
                //make new reward 200-prize of city / total prizes of route

                // distance between cities = sGraph.shortestPath(v, v+1)
                // jstar total prizes = jStar.total_prize
                // ************** jStar prize subtract last reward ******************
                int lastIndex = jStar.indexPath.getLast();
                double lastPrize = sGraph.getPrize(lastIndex);
//                System.out.println(" jStar last Prize: "+lastPrize);
                double jStarPrizeTotal = jStar.total_prize - lastPrize;
//                System.out.println(" jStar new Prizes: "+jStarPrizeTotal);

                // ******************* reward Table
//                R[(path.get(v))-1][(path.get(v + 1))-1] += (W / jStarPrizeTotal);
//                double reward = R[(path.get(v))-1][(path.get(v + 1))-1];

                // ******************* no reward Table
                //this worked
//                double reward = 2*jStar.total_prize/sGraph.shortestPath(v, v+1);
                //trying this
                double reward = 10-(sGraph.getPrize(path.get(v+1)-1)/jStarPrizeTotal);

                //Q-Table updated
                double beforeQValue = Q[(path.get(v))-1][(path.get(v + 1))-1];
//                double firstPart = (1 - alpha) * Q[path.get(v)-1][path.get(v + 1)-1];
//                double secondPart = alpha * (reward + gamma * maxQ);
                Q[(path.get(v))-1][(path.get(v + 1))-1] = (1 - alpha) * Q[path.get(v)-1][path.get(v + 1)-1] + alpha *
                        (reward + gamma * maxQ);
//                System.out.println(firstPart);
//                System.out.println(secondPart);
//                System.out.println("Before: " + beforeQValue + "    After: " + Q[(path.get(v))-1][(path.get(v + 1))-1]);
//                System.out.println(path);
//                Q[path.get(v)][path.get(v + 1)] = (1 - alpha) * q
//                        + alpha * (R[path.get(v)][path.get(v + 1)] + gamma * maxQ);
            }
            //new logic:
            if(aList[mostFitIndex].total_prize > prizeMax){
                System.out.println("Better Path Found ----------------------");
                prizeMax = aList[mostFitIndex].total_prize;
                routeMax = path;
                episodeMax = i;
            }
//            //# update explore rate
//                    explor_rate = min_explore + (max_explor - min_explore) * Math.exp(-decay_rate * i);
////            System.out.println(explor_rate);
//            // update date Q-table with highest route again.
            //This is working... but try to change the prizemax to the correct number
//            double newPrizeMax = prizeMax - sGraph.getPrize(routeMax.getLast());
//
//            for (int v = 0; v < routeMax.size()-1; v++) {
//                double maxQ = findMaxQ(routeMax.get(v+1)-1);
////                R[(routeMax.get(v))-1][(routeMax.get(v + 1))-1] += (W / newPrizeMax);
////                double reward = R[(routeMax.get(v))-1][(routeMax.get(v + 1))-1];
//
//                double reward = 10-(sGraph.getPrize(routeMax.get(v+1)-1)/newPrizeMax);
//
//                //Q-Table updated
//                double beforeQValue = Q[(routeMax.get(v))-1][(routeMax.get(v + 1))-1];
//                Q[(routeMax.get(v))-1][(routeMax.get(v + 1))-1] = (1 - alpha) * Q[routeMax.get(v)-1][routeMax.get(v + 1)-1] + alpha *
//                        (reward + gamma * maxQ);
////                System.out.println("Before: " + beforeQValue + "    After: " + Q[(routeMax.get(v))-1][(routeMax.get(v + 1))-1]);
//            }

        }

    }

    /*
     * Documentation for allQuotaMet()
     * This function checks if every agent has no extra budget to travel other than
     * goes to the destination
     * if every agent in aList does, return true, and break out of while loop. else,
     * return false, while loop continues
     */
    static boolean allQuotaMet(Agent[] aList) {
        int trueCounter = 0;

        for (int i = 0; i < aList.length; i++)
            if (aList[i].isDone)
                trueCounter++;

        if (trueCounter == aList.length)
            return true;
        else
            return false;
    }

    /*
     * Documentation for findHighestPrize
     * This function find the agent in the list with the highest prize collected
     * it then returns its index, and becomes agent jStar
     */
    static int findHighestPrize(Agent[] aList) {
        double runningHigh = Double.NEGATIVE_INFINITY;
        int index = 0;
        for (int j = 0; j < aList.length; j++) {
            System.out.println("total prize for: " + j + " " +aList[j].total_prize);
            if (aList[j].total_prize > runningHigh) {
                runningHigh = aList[j].total_prize;
                index = j;
            }
        }
        return index;
    }

    /*
     * Documentation for reset()
     * resets pertinent variables and graph flags for use in final traversal (see
     * traverse())
     */
    static void reset() {
        total_wt = 0;
        total_prize = 0;
        remainingBudget = 0.0;
        route = new ArrayList<>();
        for (int i = 0; i < statesCt; i++)
            sGraph.setMark(i, UNVISITED);
        sGraph.setMark(sGraph.getLastNode(), LAST_VISIT);
    }

    /*
     * Documentation for getNextStateFromCurState()
     * This function returns the next state for a specific agent according to action
     * rule in BC-PC-TSP
     * The q0 passed to this function will change according to different trails,
     * from small to large
     * At the very end, we won't do exploration at all
     */
    static int getNextStateFromCurState(Agent aj, int s, double q0) {
        ArrayList<Integer> feasible = new ArrayList<>();
        for (int i = 1; i < aj.getLastNode(); i++) {
            if (aj.getMark(i) == UNVISITED
                    && aj.shortestPath(s, i) + aj.shortestPath(i, aj.getLastNode()) < aj.budget - aj.total_wt) {
                feasible.add(i);
            }
        }

        if (feasible.size() == 0) {
            return aj.getLastNode();
        } else {
            Random rand = new Random();
            if (rand.nextDouble() > q0) {
                // Exploration
//                System.out.println("Exploration "+ q0);
                double[] prob = new double[feasible.size()];
                double total = 0;
                for (int i = 0; i < feasible.size(); i++) {
                    int u = feasible.get(i);
                    prob[i] = Math.pow(Q[s][u], delta) * aj.getPrize(u) / Math.pow(aj.weight(s, u), beta);
                    total += prob[i];
                }
                // uniform the distribution
                for (int i = 0; i < feasible.size(); i++) {
                    prob[i] /= total;
                }

                double target = rand.nextDouble();
                int idx = -1;
                while (target > 0) {
                    idx++;
                    target -= prob[idx];
                }
                return feasible.get(idx);
            } else {
                // Exploitation
//                System.out.println("NOT!!!!!!!!!!!!");
//                System.out.println(" _______________Exploitation");
//                int maxIdx = -1;
//                double curMax = Double.NEGATIVE_INFINITY;
//                for (int i = 0; i < feasible.size(); i++) {
//                    int u = feasible.get(i);
//                    double val = Math.pow(Q[s][u], delta) * aj.getPrize(u) / Math.pow(aj.weight(s, u), beta);
//
//                    if (val > curMax) {
//                        maxIdx = i;
//                        curMax = val;
//                    }
//                }
//                ################## MINE !!!!!!!!!!!!!!!!
                double maxNumber = -99999.9;
                int maxIndex = -1;

                // Find the index with the highest number element in each row
                for (int col = 0; col < Q[s].length; col++) {
                    if (Q[s][col] >= maxNumber && feasible.contains(col)) {
                        maxIndex = col;  // Update the index if a higher value is found
                        maxNumber = Q[s][maxIndex];
                    }
                }
//                System.out.println(maxIndex);
                return maxIndex;
            }
        }
    }


    /*
     * Documentation for getFeasibleSet
     * Given the current node s wherein the traveling salesman is located and his
     * current available budget B,
     * the set of s’s neighbor nodes that the salesman can travel to while still
     * having enough budgets to go to
     * destination node t is called node s's feasible set.
     * This function will be used across all algorithms
     */
    private static ArrayList<Integer> getFeasibleSet(int s) {
        ArrayList<Integer> feasible = new ArrayList<>();
        for (int i = 1; i < sGraph.getLastNode(); i++) {
            if (sGraph.getMark(i) == UNVISITED &&
                    sGraph.shortestPath(s, i) + sGraph.shortestPath(i, sGraph.getLastNode()) < budget - total_wt) {
                feasible.add(i);
            }
        }
        return feasible;
    }

    /*
     * Documentation for traverse() and DFSGreed_P()
     * This function handles the final traversal using the Q-Table as reference
     * DFSGreed_P finds the index of the largest Q-Value in the table for the
     * current node
     * then visits it, and repeats until the prizeGoal is met.
     * Includes printouts for debugging
     */
    private static void traverseQ() {
        reset(); // makes all static variables reset from previous algorithms...
        qTableTotal();
        remainingBudget = budget - total_wt; // updates remaining budget
    }
    static double findMaxQ(int columnIndex) {
        //find max Q from Q-Table
        double max = 0;
        for (int i = 1; i < Q.length; i++) {
            if (Q[i][columnIndex] > max) {
                max = Q[i][columnIndex];
            }
        }
        return max;
    }
    static void qTableTotal() {
        int beginIndex = 0;
        double qTableBudget = budget;
        double qTableReward = 0;
        ArrayList<Integer> qTableRoute = new ArrayList<>();
        ArrayList<Integer> feasible = new ArrayList<>();
        for(int i = 0;i<10;i++){
            feasible.add(i);
        }

        for (int i = 0; i < allNameList.size(); i++) {
            String name = allNameList.get(i);
            if (name.equals(begin.toLowerCase())) {
                beginIndex = i;
            }
        }
        System.out.println("------------- Q-Table ---------------");
//        System.out.println(qTableBudget);
        qTableRoute.add(beginIndex);
        int cityIndex = beginIndex;
        feasible.remove(cityIndex);
        //find max Q
        while(!feasible.isEmpty()){

            int maxIndex = -1;
            double maxNumber = 0.0;
            // Find the index with the highest number element in each row
            for (int col = 0; col < Q[cityIndex].length; col++) {
                double home = sGraph.shortestPath(cityIndex+1, col+1) +
                        sGraph.shortestPath(beginIndex+1, col+1);
                // System.out.println("Home " + home);
                if (Q[cityIndex][col] > maxNumber && feasible.contains(col)
                        && qTableBudget >= home) {
                    maxIndex = col;  // Update the index if a higher value is found
                    maxNumber = Q[cityIndex][maxIndex];
                }
            }
//            System.out.println("done!!! " + feasible + "    " + maxIndex + "    " + qTableBudget);
            if(maxIndex == -1){
                break;
            }
//            System.out.println(maxIndex);

            //subtract from budget
//            System.out.printf(
//                    "Going from %-18s to %-18s was %-5.2fmiles collecting $%-5d with a ratio of $%.7f/miles\n",
//                    sGraph.getName(cityIndex+1), sGraph.getName(maxIndex+1),
//                    sGraph.shortestPath(cityIndex+1, maxIndex+1),
//                    sGraph.getPrize(maxIndex+1),
//                    sGraph.getPrize(maxIndex+1) / sGraph.shortestPath(cityIndex+1, maxIndex+1));


            qTableBudget -= sGraph.shortestPath(cityIndex+1, maxIndex+1);
            qTableReward += sGraph.getPrize(maxIndex+1);
            cityIndex = maxIndex;
            for(int i = 0;i < feasible.size();i++) {
                if(feasible.get(i) == maxIndex){
                    feasible.remove(i);
                    qTableRoute.add(maxIndex);
                    break;
                }
//                double home = sGraph.shortestPath(maxIndex+1, i+1) + sGraph.shortestPath(beginIndex+1, i+1)
//                if()

            }

        }
        //get last city before home
        int cityBeforeHome = qTableRoute.getLast();
        //add home back into route
        qTableRoute.add(qTableRoute.get(0));
        //print out last city before home to home
        System.out.printf(
                "Going from %-18s to %-18s was %-5.2fmiles.\n",
                sGraph.getName(cityBeforeHome+1), sGraph.getName(beginCity),
                sGraph.shortestPath(cityBeforeHome+1, beginCity));

        qTableBudget -= sGraph.shortestPath(cityBeforeHome+1, beginCity);

//        System.out.println("Total distance: " + (budget-qTableBudget));
        System.out.println("Collected Prize: " + qTableReward);
        System.out.println("Route: " + routeString(qTableRoute));
        System.out.println(qTableRoute);
//        System.out.println("Remaining Budget: " + qTableBudget);

    }



    }