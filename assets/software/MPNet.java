import java.io.*;
import java.util.*;
import ilog.concert.*;
import ilog.cplex.*;

public class MPNet {

    // uses CPLEX
    public static boolean DEBUG = false;
    public static boolean SILENT = false;
    static int GAP = 9999;

    public static void main(String[] args) {

        Long seed = new Long(487641078);
        Random generator = new Random(seed);
        int maxstates = 0;

        if (args.length < 2 || args.length > 9) {
            System.out.println("----------- MPNet -----------");
            System.out.println("Software for computing the (hardwired and softwired) maximum parsimony scores of a phylogenetic network");
            System.out.println("Uses CPLEX");
            System.out.println("----------- USAGE -----------");
            System.out.println("java MPNet network.tree sequences.fasta [options]");
            System.out.println("network.tree should contain at least one network in e-newick format");
            System.out.println("sequences.fasta should contain, on each line, a taxon name followed by a space and a character state, or a sequence of character states");
            System.out.println("---------- OPTIONS ----------");
            System.out.println("--nolabels\t hides taxon labels");
            System.out.println("--nostates\t hides character states");
            System.out.println("--softwired\t only compute the softwired parsimony score, not the hardwired one");
            System.out.println("--hardwired\t only compute the hardwired parsimony score, not the softwired one");
            //System.out.println("--approx\t compute an approximation (faster)");
            System.out.println("--rand k\t use character states randomly chosen from 1 to k");
            System.out.println("--silent k\t do not show intermediate results");
            System.out.println("-Djava.library.path=[path to cplex.jar]\t to tell java the location of cplex.jar");
            return;
        }

        boolean nolabels = false;
        boolean nostates = false;
        boolean onlysoftwired = false;
        boolean onlyhardwired = false;
        boolean rand = false;
        int num_states = -1;
        boolean relax = false;
        for (int i = 1; i < args.length; i++) {
            if (args[i].equals("--nolabels")) {
                nolabels = true;
            } else if (args[i].equals("--nostates")) {
                nostates = true;
            } else if (args[i].equals("--softwired")) {
                onlysoftwired = true;
            } else if (args[i].equals("--hardwired")) {
                onlyhardwired = true;
            } else if (args[i].equals("--silent")) {
                SILENT = true;
            } else if (args[i].equals("--approx")) {
                relax = true;
            } else if (args[i].equals("--rand")) {
                rand = true;
                try {
                    num_states = Integer.parseInt(args[i + 1]);
                } catch (NumberFormatException nfe) {
                    System.out.println("** Integer number of states required");
                    return;
                }
                i++;
            } else if (i > 1) {
                System.out.println("Unknown option: " + args[i]);
                return;
            }
        }

        String netwerkFile = args[0];
        if(!SILENT) {
        System.out.println("\\ ** Reading e-newick from " + netwerkFile + "...");
        }
        
        File file = new File(netwerkFile);
        BufferedReader reader = null;
        String newick = null;
        Vector<String> newicks = new Vector();
        try {
            reader = new BufferedReader(new FileReader(file));
            while ((newick = reader.readLine()) != null) {
                newicks.add(newick);
            }
            reader.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

        Vector<Integer> num_taxa = new Vector();
        Vector<Integer> num_retic = new Vector();
        Vector<Integer> hardwired_score = new Vector();
        Vector<Integer> softwired_score = new Vector();
        Vector<Double> hardwired_time = new Vector();
        Vector<Double> softwired_time = new Vector();

        // loop through input networks
        for (String n : newicks) {
            if (newicks.indexOf(n) > 0 & !SILENT) {
                System.out.println("***** Results so far *****");
                System.out.println("***** Random seed: " + seed);
                System.out.println("***** Numbers of taxa: " + num_taxa.toString());
                System.out.println("***** Numbers of reticulations: " + num_retic.toString());
                System.out.println("***** Hardwired parsimony scores: " + hardwired_score.toString());
                System.out.println("***** Softwired parsimony scores: " + softwired_score.toString());
                System.out.println("***** Computation time hardwired parsimony score: " + hardwired_time.toString());
                System.out.println("***** Computation time softwired parsimony score: " + softwired_time.toString());
            }

            System.out.println("\\ ** Processing network " + (newicks.indexOf(n) + 1) + " out of " + newicks.size());

//        System.out.println("\\ ** Read the following:");
//        System.out.println("\\ " + newick);

            Long hardwiredTime = new Long(0);
            Long softwiredTime = new Long(0);

            Network.TAXON_LABELS = new Vector();
            Network.STATE_LABELS = new Vector();
            if (!SILENT) {
                System.out.println("\\ ** Parsing...");
            }
            Network N = Network.newick2netwerk(n);

            Vector<Vector<Character>> allStates = new Vector();
            Vector<String> taxa = new Vector();

            if (!rand) {
                // read character states from file
                String characterFile = args[1];
                if (!SILENT) {
                    System.out.println("\\ ** Reading character data from " + characterFile + "...");
                }

                file = new File(characterFile);
                String record = null;
                try {
                    reader = new BufferedReader(new FileReader(file));
                    Vector<Character> states = new Vector();
                    while ((record = reader.readLine()) != null) {
                        if (record.length() == 0 || record.startsWith("//")) {
                            continue; // ignore comments  and empty lines
                        }
                        if (record.startsWith(">")) {
                            if (!states.isEmpty()) {
                                allStates.add(states);
                                states = new Vector();
                            }
                            String[] data = record.split(" ");
                            taxa.add(data[0].substring(1));
                            if (data.length > 1) {
                                for (char a : data[1].toCharArray()) {
                                    states.add(a);
                                }
                            }
                        } else {
                            String[] data = record.split(" ");
                            if (data.length > 1) {
                                if (!states.isEmpty()) {
                                    allStates.add(states);
                                    states = new Vector();
                                }
                                taxa.add(data[0]);
                                for (char a : data[1].toCharArray()) {
                                    states.add(a);
                                }
                            } else {
                                for (char a : data[0].toCharArray()) {
                                    states.add(a);
                                }
                            }
                        }
                    }
                    if (!states.isEmpty()) {
                        allStates.add(states);
                    }
                    reader.close();
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    return;
                } catch (IOException e) {
                    e.printStackTrace();
                    return;
                }
            } else {
                // assign random states
                if (!SILENT) {
                    System.out.println("\\ ** Assigning random character states");
                }
                for (int i = 0; i < Network.TAXON_LABELS.size(); i++) {
                    Vector<Character> states = new Vector();
                    int c = generator.nextInt(num_states) + 1;
                    char cc = (char) (c + 65);
                    states.add(cc);
                    allStates.add(states);
                    taxa.add(Network.TAXON_LABELS.elementAt(i));
                }
            }

            if (taxa.size() != Network.TAXON_LABELS.size()) {
                System.out.println("Error: character and network do not have the same number of taxa.");
                System.out.println("Network taxa: " + Network.TAXON_LABELS.toString());
                System.out.println("Character taxa: " + taxa.toString());
                return;
            } else if (!taxa.containsAll(Network.TAXON_LABELS)) {
                // System.out.println(taxa.toString());
                // System.out.println(Network.TAXON_LABELS);
                System.out.println("Error: character and network do not have identical taxon sets.");
                System.out.println("Network taxa: " + Network.TAXON_LABELS.toString());
                System.out.println("Character taxa: " + taxa.toString());
                return;
            }

            int ss = 0; // softwired parsimony score
            int hs = 0; // hardwired parsimony score
            String uOutFile = netwerkFile + ".hardwiredPS.dot";
            String uPDFFile = uOutFile + ".pdf";
            String rOutFile = netwerkFile + ".softwiredPS.dot";
            String rPDFFile = rOutFile + ".pdf";
            String eol = System.getProperty("line.separator");
         
            // construct vector with all vertices
            // this also numbers the vertices
            N.cleanNetwork();
            int[] num = new int[1];
            num[0] = Network.TAXON_LABELS.size() + 1;
            Vector<Network> vertices = N.getVertices(num);

            if (!MPNet.SILENT) {
                System.out.println("\\ ** Network has " + Network.TAXON_LABELS.size() + " taxa.");
            }
            if (!MPNet.SILENT) {
                System.out.println("\\ ** Network has " + vertices.size() + " vertices.");
            }

            // construct vector with all edges
            Vector<Vector<Network>> edges = N.getEdges();
            if (!MPNet.SILENT) {
                System.out.println("\\ ** Network has " + edges.size() + " edges.");
            }

            // construct vector with all reticulations
            // this also numbers the vertices
            N.cleanNetwork();
            num[0] = Network.TAXON_LABELS.size() + 1;
            Vector<Network> reticulations = N.getReticulations(num);
            if (!MPNet.SILENT) {
                System.out.println("\\ ** Network has " + reticulations.size() + " reticulations.");
            }

            for (int state_index = 0; state_index < allStates.elementAt(0).size(); state_index++) {

                int hcs = 0; // hardwired parsimony score of this character
                int scs = 0; // softwired parsimony score of this character

                if (!SILENT && allStates.elementAt(0).size() > 1) {
                    System.out.println("\\ ** Processing character " + (state_index + 1) + " out of " + allStates.elementAt(0).size() + "...");
                }

                // clear existing character states
                N.clearStates();
                //N.cleanNetwork();
                Network.STATE_LABELS = new Vector();

                // add character data to network
                int k = N.setCharacterStates(taxa, allStates, state_index);
                if (k > maxstates) {
                    maxstates = k;
                }

                if (Network.STATE_LABELS.isEmpty()) {
                    System.out.println("\\ ** No states.");
                    continue;
                }

                if (Network.STATE_LABELS.size() == 1) {
                    System.out.println("\\ ** Only one state.");
                    continue;
                }
                
                if (!MPNet.SILENT) {
                    System.out.println("\\ ** Character has " + k + " states.");
                }

                // compute parsimony scores
                
                // first softwired
                if (!onlyhardwired) {

                    if (!SILENT) {
                        System.out.println("** ----- Solving softwired ILP with CPLEX ------");
                    }
                    Long startingTime = System.currentTimeMillis();
                    
                    // compute softwired parsimony score
                    computePS(true, relax, N, vertices, edges, reticulations);
                    if(relax) {
                        // in case some states are still -1 we assing them a value
                        N.roundStates(true);
                    }
                    scs = N.getScore(true);
                    N.resetSeen();
                    N.finaliseStates(state_index, true);
                    N.resetSeen();
                    
                    //update time
                    softwiredTime += System.currentTimeMillis() - startingTime;
                    if (!SILENT) {
                        System.out.println("** ----- Finished solving softwired ILP with CPLEX -----");
                    }
                }

                // now hardwired
                int rel_hcs = 0;
                if (!onlysoftwired) {

                    Long startingTime = System.currentTimeMillis();
                    if (!SILENT) {
                        System.out.println("** ----- Solving hardwired ILP with CPLEX ------");
                    }
                    
                    // compute hardwired parsimony score
                    N.clearInternalStates();
                    computePS(false, relax, N, vertices, edges, reticulations);
                    N.roundStates(false);
                    N.resetSeen();
                    hcs = N.getScore(false);
                    N.resetSeen();
                    N.finaliseStates(state_index, false);
                    N.resetSeen();
                    
                    //update time
                    hardwiredTime += System.currentTimeMillis() - startingTime;
                    if (!SILENT) {
                        System.out.println("** ----- Finished solving hardwired ILP with CPLEX -----");
                    }
                }

                ss += scs;
                hs += hcs;
            }

            // compute edge support
            N.computeRetEdgeSupport();
            
            // output
            Vector<String> networkStrings;
            if (!onlyhardwired) {

                if (newicks.size() == 1) {
                    // write output network to file
                    networkStrings = N.toDot(nolabels, nostates, true);
                    try {
                        BufferedWriter out = new BufferedWriter(new FileWriter(rOutFile));
                        for (String line : networkStrings) {
                            out.write(line + eol);
                        }
                        out.close();
                        if (!SILENT) {
                            System.out.println("** Network for softwired parsimony score in DOT format has been written to: " + rOutFile);
                        }

                    } catch (IOException e) {
                        System.out.println("** Could not write output network to file.");
                    }


                    try {
                        String line;
                        Process p = Runtime.getRuntime().exec("dot -Tpdf " + rOutFile + " -O");
                        BufferedReader bre = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                        while ((line = bre.readLine()) != null) {
                            if (!SILENT) {
                                System.out.println(line);
                            }
                        }
                        bre.close();
                        p.waitFor();
                        System.out.println("** Network for softwired parsimony score in PDF format has been written to: " + rPDFFile);

                    } catch (Exception err) {
                        System.out.println("** Could not convert network to PDF format.");
                    }
                }
            }

            if (!onlysoftwired) {

                if (newicks.size() == 1) {

                    networkStrings = N.toDot(nolabels, nostates, false);
                    try {
                        BufferedWriter out = new BufferedWriter(new FileWriter(uOutFile));
                        for (String line : networkStrings) {
                            out.write(line + eol);
                        }
                        out.close();
                    } catch (IOException e) {
                        return;
                    }

                    if (!SILENT) {
                        System.out.println("** Network for hardwired parsimony score in DOT format has been written to: " + uOutFile);
                    }

                    try {
                        String line;
                        Process p = Runtime.getRuntime().exec("dot -Tpdf " + uOutFile + " -O ");
                        BufferedReader bre = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                        while ((line = bre.readLine()) != null) {
                            System.out.println(line);
                        }
                        bre.close();
                        p.waitFor();
                        if (!SILENT) {
                            System.out.println("** Network for hardwired parsimony score in PDF format has been written to: " + uPDFFile);
                        }

                    } catch (Exception err) {
                        System.out.println("** Could not convert network to PDF format.");
                    }
                }
            }

            N.cleanNetwork();
            if (!SILENT) {
                System.out.println("***** Number of taxa: " + Network.TAXON_LABELS.size());
            }
            if (!SILENT) {
                System.out.println("***** Number of reticulations: " + reticulations.size());
                System.out.println("***** Number of characters: " + allStates.elementAt(0).size());
                System.out.println("***** Number of character states: " + maxstates);
            }

            if (!onlysoftwired & !SILENT) {
                System.out.println("***** Hardwired Parsimony Score: " + hs);
            }
            if (!onlyhardwired & !SILENT) {
                System.out.println("***** Softwired Parsimony Score: " + ss);
            }

            if (!SILENT) {
                System.out.println("***** Computation time for Hardwired Parsimony Score " + hardwiredTime / 1000 + " seconds.");
                System.out.println("***** Computation time for Softwired Parsimony Score: " + softwiredTime / 1000 + " seconds.");
            }
            num_taxa.add(Network.TAXON_LABELS.size());
            num_retic.add(reticulations.size());
            hardwired_score.add(hs);
            softwired_score.add(ss);
            hardwired_time.add(hardwiredTime / 1000.0);
            softwired_time.add(softwiredTime / 1000.0);

        }

        if (newicks.size() > 1) {
            System.out.println("***** Finished all networks");
            System.out.println("***** Random seed: " + seed);
            System.out.println("***** Numbers of taxa: " + num_taxa.toString());
            System.out.println("***** Numbers of reticulations: " + num_retic.toString());
            System.out.println("***** Hardwired parsimony scores: " + hardwired_score.toString());
            System.out.println("***** Softwired parsimony scores: " + softwired_score.toString());
            System.out.println("***** Computation time hardwired parsimony score: " + hardwired_time.toString());
            System.out.println("***** Computation time softwired parsimony score: " + softwired_time.toString());
        }

        // this can be use to compute average computation times per run
        if (newicks.size() == 60) {
            for (int run = 0; run < 6; run++) {
                double hwsum = 0;
                double swsum = 0;
                for (int i = run * 10; i < run * 10 + 10; i++) {
                    hwsum += hardwired_time.elementAt(i);
                    swsum += softwired_time.elementAt(i);
                }
                System.out.println("***** Avg computation time hardwired run " + (run + 1) + ": " + (hwsum / 10.0));
                System.out.println("***** Avg computation time softwired run " + (run + 1) + ": " + (swsum / 10.0));
            }
        }
    }
    
    public static double computePS(boolean softwired, boolean relax, Network N, Vector<Network> vertices, Vector<Vector<Network>> edges, Vector<Network> reticulations) {
        double ps = -1; // the parsimony score
        
        // Vector of Strings with ILP formulation
        Vector<String> ILPStrings = new Vector();
        String eol = System.getProperty("line.separator");
        
        // generate ILP for CPLEX
        int k = Network.STATE_LABELS.size();
        ILPStrings.add("Minimize");
        // objective function
        String obj = "";
        boolean gotone = false;
        for (int e = 0; e < edges.size(); e++) {
            Network v = edges.elementAt(e).elementAt(1);
            if (v.state == GAP) {
                continue;
            }
            if (gotone) {
                obj += " + ";
            }
            gotone = true;
            obj += "c_" + e;
        }
        ILPStrings.add(obj);
        ILPStrings.add("Subject To");
        // edge constraints
        for (int e = 0; e < edges.size(); e++) {
            // we only create variables for internal vertices
            Network u = edges.elementAt(e).elementAt(0);
            Network v = edges.elementAt(e).elementAt(1);
            boolean retedge = (v.parents.size() > 1);
            if (retedge & softwired) {
                //reticulation edge
                if (v.isLeaf && v.state != GAP) {
                    int s = v.state;
                    String xus = "x_" + u.number + "," + s;
                    String ye = "y_" + e;
                    // first constraint
                    ILPStrings.add("c_" + e + " + " + xus + " - " + ye + " >= 0");
                    // second constraint
                    ILPStrings.add("c_" + e + " - " + xus + " - " + ye + " >= -2");
                } else if (v.isLeaf) {
                    // v is a gap
                    // no constraint
                } else {
                    for (int s = 0; s < k; s++) {
                        String xus = "x_" + u.number + "," + s;
                        String xvs = "x_" + v.number + "," + s;
                        String ye = "y_" + e;
                        // first constraint
                        ILPStrings.add("c_" + e + " + " + xus + " - " + xvs + " - " + ye + " >= -1");
                        // second constraint
                        ILPStrings.add("c_" + e + " - " + xus + " + " + xvs + " - " + ye + " >= -1");
                    }
                }
            } else {
                if (v.isLeaf && v.state != GAP) {
                    int s = v.state;
                    String xus = "x_" + u.number + "," + s;
                    // first constraint
                    ILPStrings.add("c_" + e + " + " + xus + " >= 1");
                    // second constraint
                    ILPStrings.add("c_" + e + " - " + xus + " >= -1");
                } else if (v.isLeaf) {
                    // v is a gap
                    // no constraint
                } else {
                    for (int s = 0; s < k; s++) {
                        String xus = "x_" + u.number + "," + s;
                        String xvs = "x_" + v.number + "," + s;
                        // first constraint
                        ILPStrings.add("c_" + e + " + " + xus + " - " + xvs + " >= 0");
                        // second constraint
                        ILPStrings.add("c_" + e + " - " + xus + " + " + xvs + " >= 0");
                    }
                }
            }
        }
        if (softwired) {
            // reticulation constraints
            for (Network ret : reticulations) {
                boolean first = true;
                String con = "";
                for (int e = 0; e < edges.size(); e++) {
                    Network v = edges.elementAt(e).elementAt(1);
                    if (v != ret) {
                        continue;
                    }
                    if (v.state == GAP) {
                        continue;
                    }
                    if (first) {
                        con += "y_" + e;
                    } else {
                        con += " + y_" + e;
                    }
                    first = false;
                }
                con += " = 1";
                ILPStrings.add(con);
            }
        }
        // vertex constraints
        for (Network v : vertices) {
            if (v.isLeaf) {
                continue;
            }
            String con = "";
            for (int s = 0; s < k; s++) {
                if (s > 0) {
                    con += " + ";
                }
                con += "x_" + v.number + "," + s;
            }
            con += " = 1";
            ILPStrings.add(con);
        }
        // bounds
        if (relax) {
            ILPStrings.add("Bounds");
            for (Network vertex : vertices) {
                if (vertex.isLeaf) {
                    continue;
                }
                for (int s = 0; s < k; s++) {
                    ILPStrings.add("0 <= x_" + vertex.number + "," + s + " <= 1");
                }
            }
            for (int e = 0; e < edges.size(); e++) {
                Network v = edges.elementAt(e).elementAt(1);
                if (v.state != GAP) {
                    ILPStrings.add("0 <= c_" + e + " <= 1");
                }
            }
            ILPStrings.add("Binary");
            for (int e = 0; e < edges.size(); e++) {
                Network v = edges.elementAt(e).elementAt(1);
                boolean retedge = (v.parents.size() > 1);
                if (retedge & softwired & v.state != GAP) {
                    ILPStrings.add("y_" + e);
                }
            }
        } else {
            ILPStrings.add("Binary");
            for (Network vertex : vertices) {
                if (vertex.isLeaf) {
                    continue;
                }
                for (int s = 0; s < k; s++) {
                    ILPStrings.add("x_" + vertex.number + "," + s);
                }
            }
            for (int e = 0; e < edges.size(); e++) {
                Network v = edges.elementAt(e).elementAt(1);
                if (v.state != GAP) {
                    ILPStrings.add("c_" + e);
                }
                boolean retedge = (v.parents.size() > 1);
                if (retedge & softwired & v.state != GAP) {
                    ILPStrings.add("y_" + e);
                }
            }
        }
        ILPStrings.add("End");

        // write ILP to file
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter("ILP.tmp"));
            for (String line : ILPStrings) {
                out.write(line + eol);
            }
            out.close();
        } catch (IOException e) {
            return -1;
        }
        
        // ----- run CPLEX -----
        try {

            IloCplex cplex = new IloCplex();

            //! fileName is the name of the file where your ILP is
            cplex.importModel("ILP.tmp");

            //! uncomment this to suppress visual output from cplex
            cplex.setOut(null);

            //! this is the solving bit
            if (cplex.solve()) {
                IloNumVar[] var = parse(cplex);

                //! read the optimal solution
                double x[] = cplex.getValues(var);
                if (relax) {
                    // give each vertex a state with maximum x_{v,s} value
                    for (int loop = 0; loop < x.length; loop++) {
                        String varname = var[loop].getName();
                        String[] splitVarName = varname.split("_");
                        if (splitVarName[0].equals("x")) {
                            // this is an x-variable (indicating character state)
                            String[] doubleSplit = splitVarName[1].split(",");
                            int vertexNum = Integer.parseInt(doubleSplit[0]);
                            int stateValue = Integer.parseInt(doubleSplit[1]);
                            Network v = N.getVertex(vertexNum);
                            // give the corresponding vertex the right state
                            if(x[loop] > v.stateDouble) {
                                v.state = stateValue;
                                v.stateDouble = x[loop];
                            }
                        } else {
                            // this is a c-variable: skip
                        }
                    }
                } else {
                    for (int loop = 0; loop < x.length; loop++) {
                        //System.out.println(var[loop].getName() + " = " + x[loop]);
                        String varname = var[loop].getName();
                        int intvalue = (int) Math.round(x[loop]);
                        if (intvalue == 1) {
                            String[] splitVarName = varname.split("_");
                            if (splitVarName[0].equals("y")) {
                                // this is a y-variable (indicating the switching)
                                            /*
                                 int edgenum = Integer.parseInt(splitVarName[1]);
                                 Vector<Network> edge = softwiredILP.edges.elementAt(edgenum);
                                 Network u = edge.elementAt(0);
                                 Network v = edge.elementAt(1);
                                 int index = v.parents.indexOf(u);
                                 // record that this edge was switched on for this character
                                 Vector<Boolean> used = new Vector();
                                 used.setSize(v.parents.size());
                                 used.set(index, true);
                                 v.retEdgeUsed.set(state_index,used);
                                 */
                            } else if (splitVarName[0].equals("x")) {
                                // this is an x-variable (indicating character state)
                                String[] doubleSplit = splitVarName[1].split(",");
                                int vertexNum = Integer.parseInt(doubleSplit[0]);
                                int stateValue = Integer.parseInt(doubleSplit[1]);
                                // give the corresponding vertex the right state
                                N.setState(vertexNum, stateValue);
                            } else {
                                // this is a c-variable: skip
                            }
                        }
                    }
                }
            }

            //! this gets the objective function value, rounded to an int
            ps = cplex.getObjValue();

            //! this deallocates the CPLEX resources
            cplex.end();

        } catch (IloException e) {
            System.out.println("Something went wrong with CPLEX.");
            System.exit(0);
        }
        
        return ps;
    }

    private static IloNumVar[] parse(IloCplex cplex) throws IloException {
        HashSet<IloNumVar> vars = new HashSet<IloNumVar>();
        Iterator it = cplex.iterator();
        IloLinearNumExpr expr;
        IloLinearNumExprIterator it2;
        while (it.hasNext()) {
            IloAddable thing = (IloAddable) it.next();
            if (thing instanceof IloRange) {
                expr = (IloLinearNumExpr) ((IloRange) thing).getExpr();
                it2 = expr.linearIterator();
                while (it2.hasNext()) {
                    vars.add(it2.nextNumVar());
                }
            } else if (thing instanceof IloObjective) {
                expr = (IloLinearNumExpr) ((IloObjective) thing).getExpr();
                it2 = expr.linearIterator();
                while (it2.hasNext()) {
                    vars.add(it2.nextNumVar());
                }
            } else if (thing instanceof IloSOS1) {
                vars.addAll(Arrays.asList(((IloSOS1) thing).getNumVars()));
            } else if (thing instanceof IloSOS2) {
                vars.addAll(Arrays.asList(((IloSOS2) thing).getNumVars()));
            } else if (thing instanceof IloLPMatrix) {
                vars.addAll(Arrays.asList(((IloLPMatrix) thing).getNumVars()));
            }
        }
        IloNumVar[] varray = vars.toArray(new IloNumVar[1]);
        return varray;
    }
}

class Network {

    static int MAX_RET = 0; // for printing purposes
    static int GAP = 9999;
    static Vector<String> TAXON_LABELS = new Vector();
    static Vector<Character> STATE_LABELS = new Vector();
    Vector<Network> children;
    Vector<Network> parents;
    int state;
    double stateDouble;
    Vector<Character> softwiredStates;
    Vector<Character> hardwiredStates;
    Vector<Double> retEdgeSupport;
    boolean isLeaf;
    String label;
    Vector TreeVertices;
    int aafComp;
    boolean isRoot;
    int retNum;
    boolean seen;
    public int number;

    public Network() {
        isLeaf = false;
        label = null;
        parents = new Vector();
        children = new Vector();
        TreeVertices = new Vector();
        aafComp = -1;
        retNum = -1;
        seen = false;
        isRoot = false;
        number = 0;
        state = -1;
        stateDouble = -1;
        hardwiredStates = new Vector();
        softwiredStates = new Vector();
        retEdgeSupport = new Vector();
    }

    public static Network newick2netwerk(String newick) {
        if (newick.endsWith(";")) {
            int lastclosepar = newick.lastIndexOf(")");
            newick = newick.substring(0, lastclosepar + 1);
        } else {
            return null;
        }
        Network N = newick2netwerk(newick, new Vector());
        N.isRoot = true;
        N.cleanNetwork();

        // suppress indegree-1 outdegree-1
        N.suppress();

        return N;
    }
    
    public boolean setState(int vertex_num, int state) {
        if(this.number == vertex_num) {
            this.state = state;
            return true;
        }
        for (Network child : children) {
            boolean cb = child.setState(vertex_num, state);
            if(cb) {
                return true;
            }
        }
        return false;
    }
    
    public Network getVertex(int vertex_num) {
        if(this.number == vertex_num) {
            return this;
        }
        for (Network child : children) {
            Network v = child.getVertex(vertex_num);
            if(v != null) {
                return v;
            }
        }
        return null;
    }

    public void suppress() {
        for (Network child : children) {
            child.suppress();
            if (child.children.size() == 1 && child.parents.size() == 1) {
                // indegree-1 outdegree-1
                // suppress
                Network grandchild = child.children.elementAt(0);
                children.setElementAt(grandchild, children.indexOf(child));
                grandchild.parents.setElementAt(this, grandchild.parents.indexOf(child));
            }
        }
    }

    public static Network newick2netwerk(String newick, Vector<Network> reticulations) {
        int lastclosepar = newick.lastIndexOf(")");
        int lasthash = newick.lastIndexOf("#");
        int lastcolon = newick.lastIndexOf(":");

        // get rid of weights
        if (lastcolon > lastclosepar & lastcolon > lasthash) {
            return newick2netwerk(newick.substring(0, lastcolon), reticulations);
        }

        Network N = new Network();

        if (newick.startsWith("(")) {
            if (lastclosepar < newick.length() - 1 && newick.charAt(lastclosepar + 1) == '#') {
                // a new reticulation
                reticulations.add(N);
                N.retNum = new Integer(newick.substring(lastclosepar + 3, newick.length()));
                Network child = newick2netwerk(newick.substring(0, lastclosepar + 1), reticulations);
                N.children.add(child);
                child.parents.add(N);
                return N;
            } else {
                // split vertex
                int openpar = 0;
                int closepar = 0;
                int start = 1;
                Vector<String> childrenNewick = new Vector();
                for (int i = 0; i < newick.length(); i++) {
                    if (newick.charAt(i) == '(') {
                        openpar++;
                    }
                    if (newick.charAt(i) == ')') {
                        closepar++;
                    }
                    if ((openpar == closepar + 1) && (newick.charAt(i) == ',')) {
                        childrenNewick.add(newick.substring(start, i));
                        start = i + 1;
                    }
                    if (i == newick.length() - 1) {
                        childrenNewick.add(newick.substring(start, i));
                    }
                }

                for (String childNewick : childrenNewick) {
                    Network child = newick2netwerk(childNewick, reticulations);
                    N.children.add(child);
                    child.parents.add(N);
                }
                return N;
            }

        } else {
            if (newick.startsWith("#H")) {
                // a reticulation
                N.retNum = Integer.parseInt(newick.substring(2, newick.length()));
                for (Network reticulation : reticulations) {
                    if (reticulation.retNum == N.retNum) {
                        // an existing reticulation
                        N.children.add(reticulation);
                        reticulation.parents.add(N);
                        return N;
                    }
                }
            } else {
                // a leaf
                if (newick.contains("#")) {
                    // a reticulation leaf
                    int hash = newick.indexOf("#");
                    N.retNum = Integer.parseInt(newick.substring(hash + 2, newick.length()));

                    // check if we've already seen this reticulation
                    for (Network reticulation : reticulations) {
                        if (reticulation.retNum == N.retNum) {
                            // an existing reticulation
                            N.children.add(reticulation);
                            reticulation.parents.add(N);
                            return N;
                        }
                    }

                    // apparently this is a new reticulation
                    reticulations.add(N);
                    Network child = new Network();
                    child.isLeaf = true;
                    String lab = newick.substring(0, hash);
                    child.label = lab;
                    N.children.add(child);
                    child.parents.add(N);
                    TAXON_LABELS.add(lab);
                } else {
                    // a normal leaf
                    N.isLeaf = true;
                    N.label = newick;
                    TAXON_LABELS.add(newick);
                }
            }
        }
        return N;
    }

    public String toString() {
        String output;
        // returns eNewick string of the network
        if (isLeaf) {
            return label + ":1.0";
        }

        String childString1 = ((Network) children.elementAt(0)).toString();

        if (parents.size() > 1) {
            // reticulation
            if (seen) {
                output = "#H" + retNum;
            } else {
                MAX_RET++;
                retNum = MAX_RET;
                seen = true;
                if (childString1.startsWith("(")) {
                    output = childString1 + "#H" + retNum;
                } else {
                    output = "(" + childString1 + ")" + "#H" + retNum;
                }
            }
            return output;
        }

        output = "(" + childString1;

        for (int i = 1; i < children.size(); i++) {
            output += "," + ((Network) children.elementAt(i)).toString();
        }
        output += "):1.0";
        if (parents.isEmpty()) {
            output += ";";
        }
        return output;
    }
    
    public void computeRetEdgeSupport() {
        retEdgeSupport = new Vector();
        for(Network parent : parents) {
            retEdgeSupport.add(0.0);
        }
        if (parents.size() > 1) {
            for (int i = 0; i < softwiredStates.size(); i++) {
                char c = softwiredStates.elementAt(i);
                int numParentsWithC = 0;
                for(Network parent : parents) {
                    char d = parent.softwiredStates.elementAt(i);
                    if (c == d) {
                        numParentsWithC++;
                    }
                }
                //double score = 1.0 / (numParentsWithC * softwiredStates.size());
                for (int p = 0; p < parents.size(); p++) {
                    Network parent = parents.elementAt(p);
                    char d = parent.softwiredStates.elementAt(i);
                    if (c == d & numParentsWithC == 1) {
                    //if (c == d) {
                        retEdgeSupport.setElementAt(retEdgeSupport.elementAt(p) + 1, p);
                    }
                }
            }
            // normalise
            int total = 0;
            for (int p = 0; p < parents.size(); p++) {
                total += retEdgeSupport.elementAt(p);
            }
            if (total != 0) {
                for (int p = 0; p < parents.size(); p++) {
                    retEdgeSupport.setElementAt(retEdgeSupport.elementAt(p) / total, p);
                }
            }
        }
        // recurse
        for(Network child : children) {
            child.computeRetEdgeSupport();
        }
    }

    public int setCharacterStates(Vector<String> taxa, Vector<Vector<Character>> states, int state_index) {
        // returns the number of states
        if (isLeaf) {
            for (String taxon : taxa) {
                if (label.equals(taxon)) {
                    char s = states.elementAt(taxa.indexOf(taxon)).elementAt(state_index);
                    if (s != '-' & s != '?') {
                        if (!STATE_LABELS.contains(s)) {
                            STATE_LABELS.add(s);
                        }
                        this.state = STATE_LABELS.indexOf(s);
                    } else {
                        this.state = GAP;
                    }
                }
            }
        } else {
            for (Network child : children) {
                child.setCharacterStates(taxa, states, state_index);
            }
        }
        return STATE_LABELS.size();
    }

    public void resetSeen() {
        seen = false;
        if (!isLeaf) {
            for (Network child : children) {
                child.resetSeen();
            }
        }
    }

    public void cleanNetwork() {
        MAX_RET = 0;
        seen = false;
        retNum = -1;
        number = 0;
        if (!isLeaf) {
            for (Network child : children) {
                child.cleanNetwork();
            }
        }
    }

    public int getScore(boolean softwired) {
        // returns the parsimony score
        if (this.state == GAP) {
            return 0;
        }

        int s_h = 0; // hardwired score
        int s_s = 0; // softwired score
        boolean one_is_same = false;
        for (Network parent : parents) {
            if (parent.state != this.state) {
                s_h++;
            } else {
                one_is_same = true;
            }
        }
        if (!one_is_same && parents.size() > 0) {
            s_s++;
        }
        this.seen = true;
        for (Network child : children) {
            if (!child.seen) {
                int childscore = child.getScore(softwired);
                s_s += childscore;
                s_h += childscore;
            }
        }
        if (softwired) {
            return s_s;
        } else {
            return s_h;
        }
    }

    public Vector<String> toDot(boolean nolabels, boolean nostates, boolean softwired) {
        // returns vector with the network in dot format
        Vector<String> out = new Vector();
        out.add("strict digraph G {");
        int[] num = new int[1];
        num[0] = 1000;
        this.cleanNetwork();
        out.addAll(this.nodes2dot(num, nolabels, nostates, softwired));
        out.addAll(this.arcs2dot(softwired));
        out.add("}");
        this.cleanNetwork();
        return out;
    }

    public Vector<String> nodes2dot(int num[], boolean nolabels, boolean nostates, boolean softwired) {
        Vector<String> out = new Vector();
        if (number != 0) {
            return out; // already visited
        }
        if (isLeaf) {
            // this is a leaf
            number = TAXON_LABELS.indexOf(label) + 1;
            //System.out.println(number + " [shape=circle, width=0.3, label=\"" + label + " (" + STATE_LABELS.elementAt(state) + ")\"" + "];");
            //System.out.println(number + " [shape=circle, width=0.3, label=\"" + STATE_LABELS.elementAt(state) + "\"" + "];");
            String slabel = "";
            if (softwired) {
                for (char cstate : softwiredStates) {
                    if (cstate == GAP) {
                        slabel += '-';
                    } else {
                        slabel += cstate;
                    }
                }
            } else {
                for (char cstate : hardwiredStates) {
                    if (cstate == GAP) {
                        slabel += '-';
                    } else {
                        slabel += cstate;
                    }
                }
            }
            if (nolabels) {
                if (!nostates) {
                    out.add(number + " [shape=box, width=0.2, label=\"" + slabel + "\"" + "];");
                } else {
                    out.add(number + " [shape=point];");
                }
            } else {
                if (!nostates) {
                    out.add(number + " [shape=box, width=0.2, label=\"" + label + "\\n" + slabel + "\"" + "];");
                } else {
                    out.add(number + " [shape=none, label=\"" + label + "\"];");
                }
            }
        } else {
            number = num[0];
            //System.out.println(number + " [shape=point];");
            //System.out.println(number + " [shape=circle, width=0.3, label=\"" + STATE_LABELS.elementAt(state) + "\"" + "];");
            String slabel = "";
            if (softwired) {
                for (char cstate : softwiredStates) {
                    slabel += cstate;
                }
            } else {
                for (char cstate : hardwiredStates) {
                    slabel += cstate;
                }
            }
            if (!nostates) {
                out.add(number + " [shape=box, width=0.2, label=\"" + slabel + "\"" + "];");
            } else {
                out.add(number + " [shape=point];");
            }
            for (Network child : children) {
                num[0]++;
                out.addAll(child.nodes2dot(num, nolabels, nostates, softwired));
            }
        }
        return out;
    }

    public Vector<String> arcs2dot(boolean softwired) {
        // returns the arcs in dot format
        Vector<String> out = new Vector();
        for (Network child : children) {
            int intlabel = -1;
            double doublelabel = -1;
            /*
             if(child.parents.size()>1 & softwired & MPNet.USE_CPLEX) {
             // find the fraction of characters using this edge
             int index = child.parents.indexOf(this);
             doublelabel = 0.0;
             for(int s = 0; s < numchar; s++) {
             Boolean b = child.retEdgeUsed.elementAt(s).elementAt(index);
             if(b!=null && b) {
             doublelabel++;
             }
             }
             doublelabel = doublelabel / numchar;
             * */
            if (child.parents.size() > 1 & softwired) {
                // find the fraction of characters that could be inherited over this edge
                /*
                doublelabel = 0;
                for (int i = 0; i < softwiredStates.size(); i++) {
                    if ((softwiredStates.elementAt(i) == child.softwiredStates.elementAt(i)) | (child.softwiredStates.elementAt(i) == '-')) {
                        doublelabel++;
                    }
                }
                doublelabel = doublelabel / numchar;
                */
                doublelabel = child.retEdgeSupport.elementAt(child.parents.indexOf(this));
                doublelabel = Math.round(doublelabel*10.0) / 10.0;
            } else {
                // find the number of changes on this edge
                intlabel = countChanges(this, child, softwired);
            }
            if (doublelabel != -1) {
                out.add(number + " -> " + child.number + "[color=blue,label=" + doublelabel + "]");
            } else if (softwiredStates.size() == 1 | hardwiredStates.size() == 1) {
                if (intlabel > 0) {
                    out.add(number + " -> " + child.number + "[color=red]");
                } else {
                    out.add(number + " -> " + child.number + "[color=black]");
                }
            } else {
                //out.add(number + " -> " + child.number + "[colorscheme=dark28,color=" + (changes+1) + "]");
                if (intlabel > 0) {
                    out.add(number + " -> " + child.number + "[color=red,label=" + intlabel + "]");
                } else {
                    out.add(number + " -> " + child.number + "[color=black]");
                }
            }
        }
        seen = true;
        for (Network child : children) {
            if (!child.seen) {
                out.addAll(child.arcs2dot(softwired));
            }
        }
        return out;
    }

    public int countChanges(Network vertex, Network child, boolean softwired) {
        int changes = 0;
        if (softwired) {
            for (int i = 0; i < child.softwiredStates.size(); i++) {
                char childstate = child.softwiredStates.elementAt(i);
                if (childstate == '-') {
                    continue;
                }
                boolean change = true;
                for (Network parent : child.parents) {
                    if (parent.softwiredStates.elementAt(i).equals(childstate)) {
                        change = false;
                    }
                }
                if (change) {
                    changes++;
                }
            }
        } else {
            for (int i = 0; i < vertex.hardwiredStates.size(); i++) {
                char mystate = vertex.hardwiredStates.elementAt(i);
                char childstate = child.hardwiredStates.elementAt(i);
                if (childstate == '-') {
                    continue;
                }
                if (childstate != mystate) {
                    changes++;
                }
            }
        }
        return changes;
    }

    public Vector<Network> getVertices(int num[]) {
        Vector out = new Vector();
        if (number != 0) {
            return out; // already visited
        }
        if (isLeaf) {
            // this is a leaf
            number = TAXON_LABELS.indexOf(label) + 1;
            out.add(this);
        } else {
            number = num[0];
            out.add(this);
            for (Network child : children) {
                num[0]++;
                out.addAll(child.getVertices(num));
            }
        }
        return out;
    }

    public Vector<Network> getReticulations(int num[]) {
        Vector out = new Vector();
        if (number != 0) {
            return out; // already visited
        }
        if (parents.size() > 1) {
            out.add(this);
        }
        if (isLeaf) {
            number = TAXON_LABELS.indexOf(label) + 1;
        } else {
            number = num[0];
            for (Network child : children) {
                num[0]++;
                out.addAll(child.getReticulations(num));
            }
        }
        return out;
    }

    public Vector<Vector<Network>> getEdges() {
        Vector out = new Vector();
        for (Network child : children) {
            Vector edge = new Vector();
            edge.add(this);
            edge.add(child);
            out.add(edge);
        }
        seen = true;
        for (Network child : children) {
            if (!child.seen) {
                out.addAll(child.getEdges());
            }
        }
        return out;
    }

    public void roundStates(boolean softwired) {
        // rounding procedure for states
        if (state == -1) {
            if (softwired) {
                // for the softwired parsimony score, we make the state equal to the state of at least one parent
                // we choose the parent whose state is equal to the state of a maximum number of children
                if (parents.size() > 0) {
                    Network bestParent = null;
                    int bestnum = 0;
                    for (Network parent : parents) {
                        int num = 0;
                        for (Network child : children) {
                            if (parent.state == child.state) {
                                num++;
                            }
                        }
                        if (num >= bestnum) {
                            bestnum = num;
                            bestParent = parent;
                        }
                    }
                    state = bestParent.state;
                } else {
                    // the root we give state 0 (if it's not yet integral)
                    state = 0;
                }
            } else {
                // for the hardwired parsimony score, we round to the nearest integer
                state = (int) Math.rint(stateDouble);

                // might want to do something more clever here

            }
        }
        for (Network child : children) {
            child.roundStates(softwired);
        }
    }

    public void clearStates() {
        state = -1;
        stateDouble = -1;
        //if (!STATE_LABELS.isEmpty()) {
        //    STATE_LABELS = new Vector();
        //}
        for (Network child : children) {
            child.clearStates();
        }
    }
    
    public void clearInternalStates() {
        if(!isLeaf) {
            state = -1;
            stateDouble = -1;
        }
        //if (!STATE_LABELS.isEmpty()) {
        //    STATE_LABELS = new Vector();
        //}
        for (Network child : children) {
            child.clearInternalStates();
        }
    }

    public void finaliseStates(int index, boolean softwired) {
        seen = true;
        char s;
        if (state == GAP) {
            s = '-';
        } else {
            s = STATE_LABELS.elementAt(state);
        }
        if (softwired) {
            if (index >= softwiredStates.size()) {
                softwiredStates.add(s);
            } else {
                softwiredStates.setElementAt(s, index);
            }
        } else {
            if (index >= hardwiredStates.size()) {
                hardwiredStates.add(s);
            } else {
                hardwiredStates.setElementAt(s, index);
            }
        }
        for (Network child : children) {
            if (!child.seen) {
                child.finaliseStates(index, softwired);
            }
        }
    }

    public Network cloneLeaf() {
        Network leaf = new Network();
        leaf.label = label;
        leaf.isLeaf = true;
        return leaf;
    }

    public Network getTree(int index) {
        // returns tree obtained by using the reticulation edge specified by index for each reticulation
        Network tree = new Network();
        if (isLeaf) {
            return (Network) this.cloneLeaf();
        }
        for (Network child : children) {
            if (child.parents.size() == 1 || child.parents.indexOf(this) == index) {
                Network treeChild = child.getTree(index);
                tree.children.add(treeChild);
                treeChild.parents.add(tree);
            }
        }
        return tree;
    }
}