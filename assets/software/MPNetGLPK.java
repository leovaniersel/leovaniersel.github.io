import java.io.*;
import java.util.*;

public class MPNetGLPK {

    public static boolean DEBUG = false;
    public static boolean SILENT = false;

    public static void main(String[] args) {

        Long seed = new Long(487641078);
        Random generator = new Random(seed);
        int maxstates = 0;

        if (args.length < 2 || args.length > 9) {
            System.out.println("----------- MPNetGLPKGLPK -------");
            System.out.println("Software for computing the (hardwired and softwired) maximum parsimony scores of a phylogenetic network");
            System.out.println("This version uses GLPK to solve generated ILPs");
            System.out.println("----------- USAGE -----------");
            System.out.println("java MPNetGLPK network.tree sequences.fasta [options]");
            System.out.println("network.tree should contain at least one network in e-newick format");
            System.out.println("sequences.fasta should contain, on each line, a taxon name followed by a space and a character state, or a sequence of character states");
            System.out.println("---------- OPTIONS ----------");
            System.out.println("--nolabels\t hides taxon labels");
            System.out.println("--nostates\t hides character states");
            System.out.println("--softwired\t only compute the softwired parsimony score, not the hardwired one");
            System.out.println("--hardwired\t only compute the hardwired parsimony score, not the softwired one");
            System.out.println("--rand k\t use character states randomly chosen from 1 to k");
            System.out.println("--silent k\t do not show intermediate results");
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
        System.out.println("\\ ** Reading e-newick from " + netwerkFile + "...");

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

            Netwerk.TAXON_LABELS = new Vector();
            Netwerk.STATE_LABELS = new Vector();
            if (!SILENT) {
                System.out.println("\\ ** Parsing...");
            }
            Netwerk N = Netwerk.newick2netwerk(n);

//        System.out.println("\\ ** Parsed the following network in e-newick format:");
//        System.out.println("\\ " + N.toString());
//        N.cleanNetwork();

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
                //System.out.println("\\ ** Read the following character data:");

            } else {
                // assign random states
                if (!SILENT) {
                    System.out.println("\\ ** Assigning random character states");
                }
                for (int i = 0; i < Netwerk.TAXON_LABELS.size(); i++) {
                    Vector<Character> states = new Vector();
                    int c = generator.nextInt(num_states) + 1;
                    char cc = (char) (c + 65);
                    states.add(cc);
                    allStates.add(states);
                    taxa.add(Netwerk.TAXON_LABELS.elementAt(i));
                }
                //System.out.println("\\ ** Assigned the following character data:");

            }

//            for (int i = 0; i < allStates.size(); i++) {
//                System.out.println("\\ Taxon " + taxa.elementAt(i) + " has states " + allStates.elementAt(i).toString());
//            }

            if (taxa.size() != Netwerk.TAXON_LABELS.size()) {
                System.out.println("Error: character and network do not have the same number of taxa.");
                System.out.println("Network taxa: " + Netwerk.TAXON_LABELS.toString());
                System.out.println("Character taxa: " + taxa.toString());
                return;
            } else if (!taxa.containsAll(Netwerk.TAXON_LABELS)) {
                // System.out.println(taxa.toString());
                // System.out.println(Netwerk.TAXON_LABELS);
                System.out.println("Error: character and network do not have identical taxon sets.");
                System.out.println("Network taxa: " + Netwerk.TAXON_LABELS.toString());
                System.out.println("Character taxa: " + taxa.toString());
                return;
            }

            int ss = 0; // softwired parsimony score
            int hs = 0; // hardwired parsimony score
            String hsolfile = "hsol.tmp";
            String uOutFile = netwerkFile + ".hardwiredPS.dot";
            String uPDFFile = uOutFile + ".pdf";
            String ssolfile = "ssol.tmp";
            String rOutFile = netwerkFile + ".softwiredPS.dot";
            String rPDFFile = rOutFile + ".pdf";
            String eol = System.getProperty("line.separator");

            for (int state_index = 0; state_index < allStates.elementAt(0).size(); state_index++) {

                int hcs = 0; // hardwired parsimony score of this character
                int scs = 0; // softwired parsimony score of this character

                if (!SILENT && allStates.elementAt(0).size() > 1) {
                    System.out.println("\\ ** Processing character " + (state_index + 1) + " out of " + allStates.elementAt(0).size() + "...");
                }

                // clear existing character states
                N.clearStates();
                N.cleanNetwork();
                Netwerk.STATE_LABELS = new Vector();

                // add character data to network
                int k = N.setCharacterStates(taxa, allStates, state_index);
                if (k > maxstates) {
                    maxstates = k;
                }

                if (Netwerk.STATE_LABELS.isEmpty()) {
                    System.out.println("\\ ** No states.");
                    continue;
                }

                if (Netwerk.STATE_LABELS.size() == 1) {
                    System.out.println("\\ ** Only one state.");
                    continue;
                }

                // construct softwired ILP
                ILP softwiredILP = N.toILP(true, relax);
                // write to file "sILP.tmp"
                try {
                    BufferedWriter out = new BufferedWriter(new FileWriter("sILP.tmp"));
                    for (String line : softwiredILP.strings) {
                        out.write(line + eol);
                    }
                    out.close();
                } catch (IOException e) {
                    return;
                }

                // construct hardwired ILP
                ILP hardwiredILP = N.toILP(false, relax);
                // write to file "hILP.tmp"
                try {
                    BufferedWriter out = new BufferedWriter(new FileWriter("hILP.tmp"));
                    for (String line : hardwiredILP.strings) {
                        out.write(line + eol);
                    }
                    out.close();
                } catch (IOException e) {
                }

                // now we want to solve the ILPs with GLPK
                if (!onlyhardwired) {
                    if (!SILENT) {
                        System.out.println("** ----- Solving softwired ILP with GLPK ------");
                    }
                    Long startingTime = System.currentTimeMillis();
                    try {
                        String line;

                        boolean windows = System.getProperty("os.name").startsWith("Windows");
                        Process p;
                        if (windows) {
                            p = Runtime.getRuntime().exec("glpsol.exe --cpxlp -w ssol.tmp sILP.tmp");
                        } else {
                            p = Runtime.getRuntime().exec("./glpsol --lp --pcost -w ssol.tmp sILP.tmp");
                        }
                        BufferedReader bri = new BufferedReader(new InputStreamReader(p.getInputStream()));
                        BufferedReader bre = new BufferedReader(new InputStreamReader(p.getErrorStream()));

                        while ((line = bri.readLine()) != null) {
                            if (!SILENT) {
                                System.out.println(line);
                            }
                        }
                        while ((line = bre.readLine()) != null) {
                            if (!SILENT) {
                                System.out.println(line);
                            }
                        }
                        bri.close();
                        bre.close();
                        p.waitFor();
                    } catch (Exception err) {
                        err.printStackTrace();
                        return;
                    }
                    if(!SILENT) {
                        System.out.println("** ----- Finished solving softwired ILP with GLPK -----");
                    }
                    
                    softwiredTime += System.currentTimeMillis() - startingTime;
                }

                // now hardwired
                if (!onlysoftwired) {
                    Long startingTime = System.currentTimeMillis();
                    if (!SILENT) {
                        System.out.println("** ----- Solving hardwired ILP with GLPK ------");
                    }
                    try {
                        String line;
                        boolean windows = System.getProperty("os.name").startsWith("Windows");
                        Process p;
                        if (windows) {
                            p = Runtime.getRuntime().exec("glpsol.exe --cpxlp -w hsol.tmp hILP.tmp");
                        } else {
                            p = Runtime.getRuntime().exec("./glpsol --lp --pcost -w hsol.tmp hILP.tmp");
                        }

                        BufferedReader bri = new BufferedReader(new InputStreamReader(p.getInputStream()));
                        BufferedReader bre = new BufferedReader(new InputStreamReader(p.getErrorStream()));

                        while ((line = bri.readLine()) != null) {
                            if (!SILENT) {
                                System.out.println(line);
                            }
                        }
                        while ((line = bre.readLine()) != null) {
                            if (!SILENT) {
                                System.out.println(line);
                            }
                        }
                        bri.close();
                        bre.close();
                        p.waitFor();
                    } catch (Exception err) {
                        err.printStackTrace();
                        return;
                    }
                    hardwiredTime += System.currentTimeMillis() - startingTime;

                    if (!SILENT) {
                        System.out.println("** ----- Finished solving hardwired ILP with GLPK -----");
                    }
                }

                // after running the ILP solver with GLPK, we want to read its output

                // first the softwired ILP
                if (!SILENT) {
                    System.out.println("** Trying to read solution of softwired ILP from " + ssolfile + "...");
                }
                file = new File(ssolfile);
                String record = null;
                Vector<String> lines = new Vector();
                if (!onlyhardwired) {
                    Long startingTime = System.currentTimeMillis();
                    try {
                        reader = new BufferedReader(new FileReader(file));
                        while ((record = reader.readLine()) != null) {
                            if (record.length() == 0 || record.startsWith("//")) {
                                continue; // ignore comments  and empty lines
                            }
                            String[] data = record.split(" ");
                            if (data.length > 1) {
                                lines.add(data[1]);
                            } else {
                                lines.add(data[0]);
                            }
                        }
                        reader.close();
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                        return;
                    } catch (IOException e) {
                        e.printStackTrace();
                        return;
                    }
                    for (int v = 0; v < softwiredILP.vertices.size(); v++) {
                        Netwerk vertex = softwiredILP.vertices.elementAt(v);
                        int index = lines.size() - softwiredILP.vertices.size() + v;
                        Double doublevalue = Double.parseDouble(lines.elementAt(index));
                        try {
                            int intvalue = Integer.parseInt(lines.elementAt(index));
                            softwiredILP.setValue(v, intvalue);
                            if (intvalue == 1 & vertex != null) {
                                vertex.state = softwiredILP.states.elementAt(v);
                            }
                        } catch (NumberFormatException nfe) {
                            System.out.println("** Warning: nonintegral statevariable: " + doublevalue);
                            vertex.state = -2;
                        }
                    }

                    if (DEBUG) {
                        for (int i = 0; i < softwiredILP.variables.size(); i++) {
                            String var = softwiredILP.variables.elementAt(i);
                            if (var == null) {
                                continue;
                            }
                            int val = softwiredILP.values.elementAt(i);
                            System.out.println("Variable " + var + " has value " + val);
                        }
                    }

                    N.resetSeen();
                    scs = N.getScore(true);
                    N.resetSeen();
                    N.finaliseStates(state_index, true);
                    N.resetSeen();
                    //N.clearStates();
                    //N.resetSeen();

                    softwiredTime += System.currentTimeMillis() - startingTime;
                }

                // now the hardwired ILP
                if (!onlysoftwired) {
                    Long startingTime = System.currentTimeMillis();
                    if (!SILENT) {
                        System.out.println("** Trying to read solution of hardwired ILP from " + hsolfile + "...");
                    }
                    file = new File(hsolfile);
                    record = null;
                    lines = new Vector();
                    try {
                        reader = new BufferedReader(new FileReader(file));
                        while ((record = reader.readLine()) != null) {
                            if (record.length() == 0 || record.startsWith("//")) {
                                continue; // ignore comments  and empty lines
                            }
                            // where the solution is seems to depend on whether we solve an LP or ILP
                            // for an ILP only the solution is given
                            // for an LP three values are given, and the solution seems to be the second one
                            String[] data = record.split(" ");
                            if (data.length > 1) {
                                lines.add(data[1]);
                            } else {
                                lines.add(data[0]);
                            }
                        }
                        reader.close();
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                        return;
                    } catch (IOException e) {
                        e.printStackTrace();
                        return;
                    }
                    for (int v = 0; v < hardwiredILP.vertices.size(); v++) {
                        Netwerk vertex = hardwiredILP.vertices.elementAt(v);
                        if (vertex != null) {
                            int index = lines.size() - hardwiredILP.vertices.size() + v;
                            Double doublevalue = Double.parseDouble(lines.elementAt(index));
                            try {
                                int intvalue = Integer.parseInt(lines.elementAt(index));
                                if (intvalue == 1) {
                                    vertex.state = hardwiredILP.states.elementAt(v);
                                }
                            } catch (NumberFormatException nfe) {
                                System.out.println("** Warning: nonintegral statevariable: " + doublevalue);
                                vertex.state = -1;
                            }
                        }
                    }

                    N.resetSeen();
                    hcs = N.getScore(false);
                    N.resetSeen();
                    N.finaliseStates(state_index, false);
                    N.resetSeen();
                    //N.clearStates();
                    //N.resetSeen();

                    hardwiredTime += System.currentTimeMillis() - startingTime;
                }

                ss += scs;
                hs += hcs;
            }
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
            int[] num = new int[1];
            if (!SILENT) {
                System.out.println("***** Number of taxa: " + Netwerk.TAXON_LABELS.size());
            }
            num[0] = Netwerk.TAXON_LABELS.size()
                    + 1;
            Vector<Netwerk> reticulations = N.getReticulations(num);
            if (!SILENT) {
                System.out.println("***** Number of reticulations: " + reticulations.size());
                System.out.println("***** Number of characters: " + allStates.elementAt(0).size());
                System.out.println("***** Number of character states: " + maxstates);
            }
            if (!onlysoftwired & !SILENT) {
                System.out.println("***** Hardwired Parsimony Score: " + hs);
            }
            if (!onlyhardwired & !SILENT) {
                System.out.println("***** Softwired Parsimony Score is at most: " + ss);
            }
            if (!SILENT) {
                System.out.println("***** Computation time for Hardwired Parsimony Score " + hardwiredTime / 1000 + " seconds.");
                System.out.println("***** Computation time for Softwired Parsimony Score: " + softwiredTime / 1000 + " seconds.");
            }

            num_taxa.add(Netwerk.TAXON_LABELS.size());
            num_retic.add(reticulations.size());
            hardwired_score.add(hs);

            softwired_score.add(ss);

            hardwired_time.add(hardwiredTime
                    / 1000.0);
            softwired_time.add(softwiredTime
                    / 1000.0);
        }

        if (newicks.size()
                > 1) {
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

        if (false) {
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
}

class ILP {

    public final Vector<String> strings;
    public final Vector<Netwerk> vertices;
    public final Vector<Integer> states;
    public final Vector<String> variables;
    public final Vector<Vector<Netwerk>> edges;
    public Vector<Integer> values;

    public ILP(Vector<String> ILPStrings, Vector<Netwerk> ILPvertices, Vector<Vector<Netwerk>> ILPedges, Vector<Integer> ILPstates, Vector<String> ILPvariables) {
        this.strings = ILPStrings;
        this.vertices = ILPvertices;
        this.states = ILPstates;
        this.variables = ILPvariables;
        this.values = new Vector();
        this.edges = ILPedges;
    }

    public void setState(int vertex_num, int state) {
        for (Netwerk vertex : vertices) {
            if (vertex == null) {
                continue;
            }
            if (vertex.number == vertex_num) {
                vertex.state = state;
                return;
            }
        }
        System.out.println("Warning: vertex with number " + vertex_num + " not found!");
    }

    public void setValue(int index, int value) {
        if (index >= values.size()) {
            values.setSize(index + 1);
        }
        this.values.set(index, value);
    }
}

class Netwerk {

    static int MAX_RET = 0; // for printing purposes
    static int GAP = 9999;
    static Vector<String> TAXON_LABELS = new Vector();
    static Vector<Character> STATE_LABELS = new Vector();
    Vector<Netwerk> children;
    Vector<Netwerk> parents;
    int state;
    double stateDouble;
    Vector<Character> softwiredStates;
    Vector<Character> hardwiredStates;
    Vector<Vector<Boolean>> retEdgeUsed;
    boolean isLeaf;
    String label;
    Vector TreeVertices;
    int aafComp;
    boolean isRoot;
    int retNum;
    boolean seen;
    public int number;

    public Netwerk() {
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
        retEdgeUsed = new Vector();
    }

    public static Netwerk newick2netwerk(String newick) {
        if (newick.endsWith(";")) {
            int lastclosepar = newick.lastIndexOf(")");
            newick = newick.substring(0, lastclosepar + 1);
        } else {
            return null;
        }
        Netwerk N = newick2netwerk(newick, new Vector());
        N.isRoot = true;
        N.cleanNetwork();

        // suppress indegree-1 outdegree-1
        N.suppress();

        return N;
    }

    public void suppress() {
        for (Netwerk child : children) {
            child.suppress();
            if (child.children.size() == 1 && child.parents.size() == 1) {
                // indegree-1 outdegree-1
                // suppress
                Netwerk grandchild = child.children.elementAt(0);
                children.setElementAt(grandchild, children.indexOf(child));
                grandchild.parents.setElementAt(this, grandchild.parents.indexOf(child));
            }
        }
    }

    public static Netwerk newick2netwerk(String newick, Vector<Netwerk> reticulations) {
        int lastclosepar = newick.lastIndexOf(")");
        int lasthash = newick.lastIndexOf("#");
        int lastcolon = newick.lastIndexOf(":");

        // get rid of weights
        if (lastcolon > lastclosepar & lastcolon > lasthash) {
            return newick2netwerk(newick.substring(0, lastcolon), reticulations);
        }

        Netwerk N = new Netwerk();

        if (newick.startsWith("(")) {
            if (lastclosepar < newick.length() - 1 && newick.charAt(lastclosepar + 1) == '#') {
                // a new reticulation
                reticulations.add(N);
                N.retNum = new Integer(newick.substring(lastclosepar + 3, newick.length()));
                Netwerk child = newick2netwerk(newick.substring(0, lastclosepar + 1), reticulations);
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
                    Netwerk child = newick2netwerk(childNewick, reticulations);
                    N.children.add(child);
                    child.parents.add(N);
                }
                return N;
            }

        } else {
            if (newick.startsWith("#H")) {
                // a reticulation
                N.retNum = Integer.parseInt(newick.substring(2, newick.length()));
                for (Netwerk reticulation : reticulations) {
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
                    for (Netwerk reticulation : reticulations) {
                        if (reticulation.retNum == N.retNum) {
                            // an existing reticulation
                            N.children.add(reticulation);
                            reticulation.parents.add(N);
                            return N;
                        }
                    }

                    // apparently this is a new reticulation
                    reticulations.add(N);
                    Netwerk child = new Netwerk();
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

        String childString1 = ((Netwerk) children.elementAt(0)).toString();

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
            output += "," + ((Netwerk) children.elementAt(i)).toString();
        }
        output += "):1.0";
        if (parents.isEmpty()) {
            output += ";";
        }
        return output;
    }

    public int setCharacterStates(Vector<String> taxa, Vector<Vector<Character>> states, int state_index) {
        // returns the number of states
        this.retEdgeUsed.add(new Vector());
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
            for (Netwerk child : children) {
                child.setCharacterStates(taxa, states, state_index);
            }
        }
        return STATE_LABELS.size();
    }

    public void resetSeen() {
        seen = false;
        if (!isLeaf) {
            for (Netwerk child : children) {
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
            for (Netwerk child : children) {
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
        for (Netwerk parent : parents) {
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
        for (Netwerk child : children) {
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
            for (Netwerk child : children) {
                num[0]++;
                out.addAll(child.nodes2dot(num, nolabels, nostates, softwired));
            }
        }
        return out;
    }

    public Vector<String> arcs2dot(boolean softwired) {
        // returns the arcs in dot format
        int numchar = softwiredStates.size();
        Vector<String> out = new Vector();
        for (Netwerk child : children) {
            int intlabel = -1;
            double doublelabel = -1;
            /*
             if(child.parents.size()>1 & softwired & MPNetGLPK.USE_CPLEX) {
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
                doublelabel = 0;
                for (int i = 0; i < softwiredStates.size(); i++) {
                    if ((softwiredStates.elementAt(i) == child.softwiredStates.elementAt(i)) | (child.softwiredStates.elementAt(i) == '-')) {
                        doublelabel++;
                    }
                }
                doublelabel = doublelabel / numchar;
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
        for (Netwerk child : children) {
            if (!child.seen) {
                out.addAll(child.arcs2dot(softwired));
            }
        }
        return out;
    }

    public int countChanges(Netwerk vertex, Netwerk child, boolean softwired) {
        int changes = 0;
        if (softwired) {
            for (int i = 0; i < child.softwiredStates.size(); i++) {
                char childstate = child.softwiredStates.elementAt(i);
                if (childstate == '-') {
                    continue;
                }
                boolean change = true;
                for (Netwerk parent : child.parents) {
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

    public Vector<Netwerk> getVertices(int num[]) {
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
            for (Netwerk child : children) {
                num[0]++;
                out.addAll(child.getVertices(num));
            }
        }
        return out;
    }

    public Vector<Netwerk> getReticulations(int num[]) {
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
            for (Netwerk child : children) {
                num[0]++;
                out.addAll(child.getReticulations(num));
            }
        }
        return out;
    }

    public Vector<Vector<Netwerk>> getEdges() {
        Vector out = new Vector();
        for (Netwerk child : children) {
            Vector edge = new Vector();
            edge.add(this);
            edge.add(child);
            out.add(edge);
        }
        seen = true;
        for (Netwerk child : children) {
            if (!child.seen) {
                out.addAll(child.getEdges());
            }
        }
        return out;
    }

    public ILP toILP(boolean softwired, boolean relax) {
        // outputs Vector of Strings with ILP formulation
        Vector<String> ILPStrings = new Vector();
        // and a Vector with the vertices of the variables in the ILP formulation
        Vector<Netwerk> ILPVertices = new Vector();
        // and a Vector with the states of the variables in the ILP variables
        Vector<Integer> ILPStates = new Vector();
        // a vector with the variable names
        Vector<String> ILPVariables = new Vector();

        // construct vector with all vertices
        // this also numbers the vertices
        cleanNetwork();
        int[] num = new int[1];
        num[0] = TAXON_LABELS.size() + 1;
        Vector<Netwerk> vertices = getVertices(num);

        if (softwired & !MPNetGLPK.SILENT) {
            System.out.println("\\ ** Network has " + TAXON_LABELS.size() + " taxa.");
        }
        if (softwired & !MPNetGLPK.SILENT) {
            System.out.println("\\ ** Network has " + vertices.size() + " vertices.");
        }

        // construct vector with all edges
        Vector<Vector<Netwerk>> edges = getEdges();
        if (softwired & !MPNetGLPK.SILENT) {
            System.out.println("\\ ** Network has " + edges.size() + " edges.");
        }

        // construct vector with all reticulations
        // this also numbers the vertices
        cleanNetwork();
        num[0] = TAXON_LABELS.size() + 1;
        Vector<Netwerk> reticulations = getReticulations(num);
        if (softwired & !MPNetGLPK.SILENT) {
            System.out.println("\\ ** Network has " + reticulations.size() + " reticulations.");
        }

        // generate ILP for CPLEX
        int k = Netwerk.STATE_LABELS.size();
        if (softwired & !MPNetGLPK.SILENT) {
            System.out.println("\\ ** Character has " + k + " states.");
        }

        ILPStrings.add("Minimize");
        // objective function
        String obj = "";
        boolean gotone = false;
        for (int e = 0; e < edges.size(); e++) {
            Netwerk v = edges.elementAt(e).elementAt(1);
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
            Netwerk u = edges.elementAt(e).elementAt(0);
            Netwerk v = edges.elementAt(e).elementAt(1);
            boolean retedge = (v.parents.size() > 1);
            if (retedge & softwired) {
                //reticulation edge
                if (v.isLeaf && v.state != GAP) {
                    int s = v.state;
                    String xus = "x_" + u.number + "," + s;
                    if (!ILPVariables.contains(xus)) {
                        ILPVertices.add(u);
                        ILPStates.add(s);
                        ILPVariables.add(xus);
                    }
                    String ye = "y_" + e;
                    if (!ILPVariables.contains(ye)) {
                        ILPVariables.add(ye);
                        ILPVertices.add(null);
                        ILPStates.add(null);
                    }
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
                        if (!ILPVariables.contains(xus)) {
                            ILPVertices.add(u);
                            ILPStates.add(s);
                            ILPVariables.add(xus);
                        }
                        String xvs = "x_" + v.number + "," + s;
                        if (!ILPVariables.contains(xvs)) {
                            ILPVertices.add(v);
                            ILPStates.add(s);
                            ILPVariables.add(xvs);
                        }

                        String ye = "y_" + e;
                        if (!ILPVariables.contains(ye)) {
                            ILPVariables.add(ye);
                            ILPVertices.add(null);
                            ILPStates.add(null);
                        }

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
                    if (!ILPVariables.contains(xus)) {
                        ILPVertices.add(u);
                        ILPStates.add(s);
                        ILPVariables.add(xus);
                    }
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
                        if (!ILPVariables.contains(xus)) {
                            ILPVertices.add(u);
                            ILPStates.add(s);
                            ILPVariables.add(xus);
                        }
                        String xvs = "x_" + v.number + "," + s;
                        if (!ILPVariables.contains(xvs)) {
                            ILPVertices.add(v);
                            ILPStates.add(s);
                            ILPVariables.add(xvs);
                        }
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
            for (Netwerk ret : reticulations) {
                boolean first = true;
                String con = "";
                for (int e = 0; e < edges.size(); e++) {
                    Netwerk v = edges.elementAt(e).elementAt(1);
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
        for (Netwerk v : vertices) {
            if (v.isLeaf) {
                continue;
            }
            String con = "";
            for (int s = 0; s < k; s++) {
                String xvs = "x_" + v.number + "," + s;
                if (!ILPVariables.contains(xvs)) {
                    ILPVertices.add(v);
                    ILPStates.add(s);
                    ILPVariables.add(xvs);
                }

                if (s > 0) {
                    con += " + ";
                }
                con += "x_" + v.number + "," + s;
            }
            con += " = 1";
            ILPStrings.add(con);
        }
        // integrality constraints
        ILPStrings.add("Binary");
        for (Netwerk vertex : vertices) {
            if (vertex.isLeaf) {
                continue;
            }
            for (int s = 0; s < k; s++) {
                // should not be necessary
                String xvs = "x_" + vertex.number + "," + s;
                if (!ILPVariables.contains(xvs)) {
                    ILPVertices.add(vertex);
                    ILPStates.add(s);
                    ILPVariables.add(xvs);
                }

                ILPStrings.add("x_" + vertex.number + "," + s);
            }
        }
        for (int e = 0; e < edges.size(); e++) {
            Netwerk v = edges.elementAt(e).elementAt(1);
            if (v.state != GAP) {
                ILPStrings.add("c_" + e);
            }
            boolean retedge = (v.parents.size() > 1);
            if (retedge & softwired & v.state != GAP) {
                ILPStrings.add("y_" + e);
            }
        }
        ILPStrings.add("End");

        ILP ilp = new ILP(ILPStrings, ILPVertices, edges, ILPStates, ILPVariables);
        return ilp;
    }

    public void roundStates(boolean softwired) {
        // rounding procedure for states
        if (state == -1) {
            if (softwired) {
                // for the softwired parsimony score, we make the state equal to the state of at least one parent
                // we choose the parent whose state is equal to the state of a maximum number of children
                if (parents.size() > 0) {
                    Netwerk bestParent = null;
                    int bestnum = 0;;
                    for (Netwerk parent : parents) {
                        int num = 0;
                        for (Netwerk child : children) {
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
        for (Netwerk child : children) {
            child.roundStates(softwired);
        }
    }

    public void clearStates() {
        state = -1;
        stateDouble = -1;
        //if (!STATE_LABELS.isEmpty()) {
        //    STATE_LABELS = new Vector();
        //}
        for (Netwerk child : children) {
            child.clearStates();
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
        for (Netwerk child : children) {
            if (!child.seen) {
                child.finaliseStates(index, softwired);
            }
        }
    }

    public Netwerk cloneLeaf() {
        Netwerk leaf = new Netwerk();
        leaf.label = label;
        leaf.isLeaf = true;
        return leaf;
    }

    public Netwerk getTree(int index) {
        // returns tree obtained by using the reticulation edge specified by index for each reticulation
        Netwerk tree = new Netwerk();
        if (isLeaf) {
            return (Netwerk) this.cloneLeaf();
        }
        for (Netwerk child : children) {
            if (child.parents.size() == 1 || child.parents.indexOf(this) == index) {
                Netwerk treeChild = child.getTree(index);
                tree.children.add(treeChild);
                treeChild.parents.add(tree);
            }
        }
        return tree;
    }
}