
import java.io.*;
import java.util.*;
import ilog.concert.*;
import ilog.cplex.*;

// by Leo van Iersel (2013)
public class NonbinaryCycleKiller {

    static boolean DEBUG = false;
    static String EOL = System.getProperty("line.separator");
    static List<String> TAXA = new LinkedList();
    static String MODE = "maf";
    static boolean ONECORE = false;
    static boolean REDUCTION = true;
    static boolean MAX = false;

    public static void dfvs(String[] args) {

        String file = args[0];
        BufferedReader reader;
        String newick = null;
        List<List<Integer>> outlists = new LinkedList();
        try {
            reader = new BufferedReader(new FileReader(file));
            while ((newick = reader.readLine()) != null && newick.length() > 0) {
                String[] split = newick.split(" ");
                List<Integer> outlist = new LinkedList();
                for (int i = 0; i < split.length; i++) {
                    outlist.add(Integer.parseInt(split[i]));
                }
                outlists.add(outlist);
            }
            reader.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error reading file: " + file);
            return;
        } catch (IOException e) {
        }
        int[] weights = new int[outlists.size()];
        for (int i = 0; i < outlists.size(); i++) {
            weights[i] = 1;;
        }
        List<List<Integer>> inlists = new LinkedList();
        for (int i = 0; i < outlists.size(); i++) {
            inlists.add(new LinkedList());
        }
        for (int i = 0; i < outlists.size(); i++) {
            for (int x : outlists.get(i)) {
                inlists.get(x).add(new Integer(i));
            }
        }

        List<Integer> dfvs = solveDFVS(outlists, inlists, weights, new LinkedList());
        System.out.println("Directed feedback vertex set:");
        System.out.println(dfvs.toString());
    }

    public static void main(String[] args) {
        System.out.println("");
        System.out.println("----------------------------------------------------------------------------");
        System.out.println("");
        System.out.println("Nonbinary Cycle Killer version 1.3");
        System.out.println("");
        System.out.println("-------------- USAGE --------------------------------------------------------");
        System.out.println("");
        System.out.println("java -Djava.library.path=libpath -classpath jarpath;. NonbinaryCycleKiller test.tree");
        System.out.println("");
        System.out.println("test.tree \t text file containing two trees in newick format, one per line");
        System.out.println("libpath\t\t path to CPLEX shared library on your computer");
        System.out.println("jarpath\t\t path to cplex.jar on your computer");
        System.out.println("");
        System.out.println("-------------- OPTIONS ------------------------------------------------------");
        System.out.println("");
        System.out.println("--maf\t\t (default) use MAF FPT algorithm to find an initial agreement forest");
        System.out.println("--mafapp\t use MAF approximation algorithm to find an initial agreement forest");
        System.out.println("--rspr\t\t use rSPR FPT algorithm to find an initial agreement forest");
        System.out.println("--rsprapp\t use rSPR approximation algorithm to find an initial agreement forest");
        System.out.println("--given af.txt \t use given initial agreement forest in file af.txt");
        System.out.println("--nored \t do not use cluster reduction in MAF");
        System.out.println("--onecore \t let CPLEX use only one core");
        System.out.println("--dfvs\t\t only find directed feedback vertex set of a given digraph");
        System.out.println("--max\t\t merge components at the end to ensure maximality");
        System.out.println("");
        System.out.println("-----------------------------------------------------------------------------");
        System.out.println("");
        System.out.println("By Leo van Iersel, Steven Kelk, Nela Lekic, Celine Scornavacca and Leen Stougie (2013)");
        System.out.println("http://homepages.cwi.nl/~iersel/cyclekiller/");
        System.out.println("");
        System.out.println("-----------------------------------------------------------------------------");
        System.out.println("");

        Long startingTime = System.currentTimeMillis();

        if (args.length == 0) {
            System.out.println("Error: at least one argument required");
            return;
        }

        String afFile = null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--dfvs")) {
                dfvs(args);
                return;
            } else if (args[i].equals("--maf")) {
                MODE = "maf";
            } else if (args[i].equals("--mafapp")) {
                MODE = "mafapp";
            } else if (args[i].equals("--rspr")) {
                MODE = "rspr";
            } else if (args[i].equals("--rsprapp")) {
                MODE = "rsprapp";
            } else if (args[i].equals("--onecore")) {
                ONECORE = true;
            } else if (args[i].equals("--nored")) {
                REDUCTION = false;
            } else if (args[i].equals("--max")) {
                MAX = true;
            } else if (args[i].equals("--given")) {
                MODE = "given";
                if(args.length >= i) {
                    afFile = args[i+1];
                    i++;
                } else {
                    System.out.println("Error: no agreement forest file specified.");
                    return;
                }
            } else if(i > 0) {
                System.out.println("Unknown option: " + args[i]);
                return;
            }
        }

        List<int[]> afComponents = new LinkedList();
        List<Baum> trees = new LinkedList();
        
        String treeFile = args[0];
        int d = treeFile.lastIndexOf('.');
        String outFile = treeFile.substring(0, d) + "-network.tree";
        String aafFile = treeFile.substring(0, d) + "-aaf.txt";

        // read in the input trees
        File file = new File(treeFile);
        BufferedReader reader;
        String newick = null;
        try {
            reader = new BufferedReader(new FileReader(file));
            while ((newick = reader.readLine()) != null && newick.length() > 0) {
                Baum tree = Baum.newickToBaum(newick);
                trees.add(tree);
            }
            reader.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error reading file: " + treeFile);
            return;
        } catch (IOException e) {
        }

        if (trees.size() < 2) {
            System.out.println("Error: not enough input trees!");
            return;
        }
        
        if (trees.size() > 2) {
            System.out.println("Error: too many input trees!");
            return;
        }
        
        if( !trees.get(0).isBinary() && trees.get(1).isBinary()) {
            System.out.println("First tree is not binary.");
            System.out.println("Swapping first and second tree.");
            Baum tree = trees.get(1);
            trees.set(1,trees.get(0));
            trees.set(0,tree);
        }

        if( (MODE.equals("rspr") | MODE.equals("rsprapp")) &&  !trees.get(0).isBinary()) {
            System.out.println("Error: rSPR can only be used when at least one of the input trees is binary!");
            return;
        }

        List<Baum> leaves = trees.get(0).getLeaves();
        for (Baum leaf : leaves) {
            if (!TAXA.contains(leaf.label)) {
                TAXA.add(leaf.label);
            } else {
                System.out.println("Error: first tree contains duplicate taxa!");
                return;
            }
            // label the leaf by an int
            leaf.number = TAXA.indexOf(leaf.label);
        }

        leaves = trees.get(1).getLeaves();
        List<String> tree1taxa = new LinkedList();
        for (Baum leaf : leaves) {
            if (!TAXA.contains(leaf.label)) {
                System.out.println("Error: trees do not have same taxa!");
                return;
            }
            if (!tree1taxa.contains(leaf.label)) {
                tree1taxa.add(leaf.label);
            } else {
                System.out.println("Error: second tree contains duplicate taxa!");
                return;
            }
            // label the leaf by an int
            leaf.number = TAXA.indexOf(leaf.label);
        }

        List<Integer> cluster1 = trees.get(0).getCluster();
        List<Integer> cluster2 = trees.get(1).getCluster();
        if (cluster1.size() != cluster2.size()) {
            System.out.println("Error: trees do not have the same number of taxa!");
            return;
        }

        // print the trees
        System.out.println("Relabelled trees:");
        for (Baum tree : trees) {
            int t = trees.indexOf(tree);
            System.out.println("Tree " + (t + 1) + ": " + tree.toString());
        }

        if (MODE.equals("rspr") | MODE.equals("rsprapp")) {
            
            if (MODE.equals("rspr")) {
                System.out.println("");
                System.out.println("---------------------------------------------------------------");
                System.out.println("Using rSPR FPT algorithm to find an initial agreement forest...");
                System.out.println("---------------------------------------------------------------");
                System.out.println("");
            } else {

                System.out.println("");
                System.out.println("-------------------------------------------------------------------------");
                System.out.println("Using rSPR approximation algorithm to find an initial agreement forest...");
                System.out.println("-------------------------------------------------------------------------");
                System.out.println("");
            }
            
            String f2String = "";
            try {
                String line;
                OutputStream stdin = null;
                InputStream stderr = null;
                InputStream stdout = null;

                String cmd = null;
                boolean windows = System.getProperty("os.name").startsWith("Windows");
                if (windows) {
                    cmd = "rspr";
                } else {
                    cmd = "./rspr";
                }
                if (MODE.equals("rsprapp")) {
                    cmd += " -approx";
                }
                // launch the command and grab stdin/stdout and stderr
                Process process = Runtime.getRuntime().exec(cmd);
                stdin = process.getOutputStream();
                stderr = process.getErrorStream();
                stdout = process.getInputStream();

                for(Baum T : trees) {
                    T.relabel("p","p_RELABELLED");
                    String n = T.toStringLabels() + EOL;
                    //System.out.println(n);
                    T.relabel("p_RELABELLED","p");
                    stdin.write(n.getBytes());
                    stdin.flush();
                }
                
                stdin.close();

                // clean up if any output in stdout
                BufferedReader brCleanUp = new BufferedReader(new InputStreamReader(stdout));
                while ((line = brCleanUp.readLine()) != null) {
                    System.out.println(line);
                    if(MODE.equals("rspr") & line.startsWith("F2: ")) {
                        f2String = line;
                    } else if(MODE.equals("rsprapp") & line.startsWith("approx F2: ")) {
                        f2String = line.substring(8);
                    }
                }
                brCleanUp.close();

                // clean up if any output in stderr
                brCleanUp = new BufferedReader(new InputStreamReader(stderr));
                while ((line = brCleanUp.readLine()) != null) {
                    System.out.println ("[Error:] " + line);
                }
                brCleanUp.close();
            } catch (Exception err) {
                err.printStackTrace();
            }
            
            System.out.println("");
            System.out.println("--------------");
            System.out.println("rSPR finished!");
            System.out.println("--------------");
            System.out.println("");
            
            String[] split = f2String.split(" ");
            List<String> rsprTaxa = new LinkedList();
            for (int i = 1; i < split.length; i++) {
                if(split[i].equals("p")) {
                    // this is the dummy root in a separate component
                    System.out.println("Agreement forest computed by rSPR contained a single component p. Ignoring this component.");
                    System.out.println("");
                    continue;
                }
                if(split[i].contains("<")) {
                    split[i] = split[i].replace('<','(').replace('>',')');
                    //System.out.println(split[i]);
                }
                Baum T = Baum.newickToBaum(split[i]);
                List<Baum> L = T.getLeaves();
                for(Baum l : L) {
                    if(l.label.equals("p")) {
                        L.remove(l);
                        break;
                    }
                }
                T.relabel("p_RELABELLED","p");
                int[] comp = new int[L.size()];
                for (int j = 0; j < L.size(); j++) {
                    String s = L.get(j).label;
                    if(!TAXA.contains(s)) {
                        System.out.println("Error: rSPR output contains incorrect label!");
                        System.out.println("Try using --maf instead");
                        return;
                    }
                    if(rsprTaxa.contains(s)) {
                        System.out.println("Error: rSPR output contains duplicate taxa!");
                        System.out.println("Try using --maf instead");
                        return;
                    } else {
                        rsprTaxa.add(s);
                    }
                    comp[j] = TAXA.indexOf(s);
                }
                afComponents.add(Baum.sortCluster(comp));
            }
            if(rsprTaxa.size() < TAXA.size()) {
                System.out.println("Error: rSPR output does not contain all taxa!");
                System.out.println("Agreement forest found by rSPR contains " + rsprTaxa.size() + " taxa while input trees have " + TAXA.size() + " taxa.");
                System.out.println("Try using --maf instead");
                System.out.println("");
                return;
            }
            
        } else if (MODE.equals("maf") | MODE.equals("mafapp")) {
            // use the MAF program to find an agreement forest
            if (MODE.equals("maf")) {
                System.out.println("");
                System.out.println("--------------------------------------------------------------");
                System.out.println("Using MAF FPT algorithm to find an initial agreement forest...");
                System.out.println("--------------------------------------------------------------");
                System.out.println("");
            } else {
                System.out.println("");
                System.out.println("------------------------------------------------------------------------");
                System.out.println("Using MAF approximation algorithm to find an initial agreement forest...");
                System.out.println("------------------------------------------------------------------------");
                System.out.println("");
            }
            
            String[] mafargs;
            if (MODE.equals("mafapp")) {
                if(REDUCTION) {
                    mafargs = new String[2];
                    mafargs[0] = treeFile;
                    mafargs[1] = "--nofpt";
                } else {
                    mafargs = new String[3];
                    mafargs[0] = treeFile;
                    mafargs[1] = "--nofpt";
                    mafargs[2] = "--nored";
                }
                afFile = "app_maf.txt";
            } else {
                if(REDUCTION) {
                    mafargs = new String[1];
                    mafargs[0] = treeFile;
                } else {
                    mafargs = new String[2];
                    mafargs[0] = treeFile;
                    mafargs[1] = "--nored";
                }
                afFile = "maf.txt";
            }
            MAF.main(mafargs);
            
            System.out.println("");
            System.out.println("-----------------------");
            System.out.println("MAF algorithm finished!");
            System.out.println("-----------------------");
            System.out.println("");
            
        }
        
        if (MODE.equals("maf") | MODE.equals("mafapp") | MODE.equals("given")) {
            
            if(MODE.equals("given")) {
                System.out.println("");
                System.out.println("-------------------------------------------------");
                System.out.println("Reading initial agreement forest from: " + afFile);
                System.out.println("-------------------------------------------------");
                System.out.println("");
            }
            
            // read in the agreement forest
            File file2 = new File(afFile);
            afComponents = new LinkedList();
            BufferedReader reader2 = null;
            try {
                reader2 = new BufferedReader(new FileReader(file2));
                String text = null;
                while ((text = reader2.readLine()) != null) {
                    if (!text.equals("")) {
                        String[] split = text.split(" ");
                        int[] comp = new int[split.length];
                        for (int i = 0; i < split.length; i++) {
                            comp[i] = TAXA.indexOf(split[i]);
                        }
                        afComponents.add(Baum.sortCluster(comp));
                    }
                }

                if (reader2 != null) {
                    reader2.close();
                }
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        // print AF components
        System.out.println("Agreement forest:");
        for (int c = 0; c < afComponents.size(); c++) {
            System.out.print("Component " + (c + 1) + ": ");
            int[] component = afComponents.get(c);
            for (int i = 0; i < component.length; i++) {
                System.out.print(component[i] + " ");
            }
            System.out.print(EOL);
        }
        
        // get the components of the agreement forest in tree-form, and simultaneously check if the agreement forest is correct
        System.out.println("Checkin if the agreement forest is correct...");
        List<Baum> afBaumComponents = getBaumComponents(afComponents,trees);
        if(afBaumComponents == null) {
            System.out.println("Error: incorrect agreement forest!");
        } else {
            System.out.println("Agreement forest is correct.");
        }
        
        // convert the agreement forest to an input graph for DFVS
        DFVSInstance instance = afToGraph(afBaumComponents,trees,true);
        int[][] adjacencyMatrix = instance.adjacencyMatrix;
        int[] weights = instance.weights;
        List<Baum> intAFVertices = instance.intAFVertices;
        int n1 = intAFVertices.size();
        List<Baum> nonRootAFVertices = instance.nonRootAFVertices;
       
        // solve DFVS
        System.out.println("Solving DFVS...");
        List<Integer> DFVS = solveDFVS(adjacencyMatrix, weights);

        if (DFVS != null) {
            System.out.println("Optimal DFVS found!");
        } else {
            System.out.println("Error: no DFVS found!");
            return;
        }

        if(DEBUG) {
            System.out.println("DFVS:" + DFVS.toString());
        } else {
            System.out.println("DFVS size: " + DFVS.size());
        }

        // now we need to remove the correseponding vertices and edges from the agreement forest

        // first we find the corresponding vertices and edges of the agreement forest
        List<Baum> DFVSVertices = new LinkedList();
        List<Baum> DFVSEdges = new LinkedList();
        for (int i : DFVS) {
            if (i < n1) {
                DFVSVertices.add(intAFVertices.get(i));
            } else {
                DFVSEdges.add(nonRootAFVertices.get(i - n1));
            }
        }

        // remove the DFVS vertices from the AF
        for (Baum v : DFVSVertices) {
            while (!v.children.isEmpty()) {
                Baum child = v.children.get(0);
                v.children.remove(0);
                child.parent = null;
                if (!DFVSVertices.contains(child)) {
                    afBaumComponents.add(child);
                }
            }
            afBaumComponents.remove(v);
        }

        // remove the DFVS edges from the AF
        for (Baum v : DFVSEdges) {
            Baum u = v.parent;
            u.children.remove(v);
            v.parent = null;
            if (!DFVSVertices.contains(v)) {
                afBaumComponents.add(v);
            }
        }
        
        // clean up the AAF
        for (Baum tree : afBaumComponents) {
            int i = afBaumComponents.indexOf(tree);
            tree.removeOutdegreeZero();
            tree = tree.removeFakeRoots();
            afBaumComponents.set(i, tree);
        }

        // we found this acyclic agreement forest (in tree form):
        System.out.println("Resulting acyclic agreement forest (in tree form):");
        for (Baum tree : afBaumComponents) {
            System.out.println(tree.toString());
        }

        // get the partition form of the AAF
        List<int[]> aafComponents = new LinkedList();
        for (Baum tree : afBaumComponents) {
            List<Integer> cluster = tree.getCluster();
            int[] intcluster = new int[cluster.size()];
            for (int i = 0; i < cluster.size(); i++) {
                intcluster[i] = cluster.get(i);
            }
            aafComponents.add(intcluster);
        }
        
        Long time = System.currentTimeMillis() - startingTime;
        System.out.println("Computation time so far: " + time / 1000 + " seconds.");
        
        if(MAX) {
            // try merging components until a maximal AAF is obtained
            System.out.println("Trying to merge components...");
            boolean success = true;
            boolean didone = false;
            while(success) {
                for(int i = 0; i < aafComponents.size() - 1; i++) {
                    for(int j = i+1; j < aafComponents.size(); j++) {
                        // try merging component i and j
                        int[] compi = aafComponents.get(i);
                        int[] compj = aafComponents.get(j);
                        int[] merged = new int[compi.length + compj.length];
                        System.arraycopy(compi, 0, merged, 0, compi.length);
                        System.arraycopy(compj, 0, merged, compi.length, compj.length);
                        List<int[]> newComponents = new LinkedList();
                        for(int k = 0; k < i; k++) {
                            newComponents.add(aafComponents.get(k));
                        }
                        newComponents.add(merged);
                        for(int k = i+1; k < aafComponents.size(); k++) {
                            if(k!=j) {
                                newComponents.add(aafComponents.get(k));
                            }
                        }
                        // clone the trees
                        List<Baum> treeClones = new LinkedList();
                        for(Baum tree : trees) {
                            treeClones.add(tree.clone());
                        }
                        // get the components of the agreement forest in tree-form, and simultaneously check if the agreement forest is correct
                        List<Baum> newBaumComponents = getBaumComponents(newComponents, treeClones);
                        if (newBaumComponents == null) {
                            continue;
                        }

                        // convert the agreement forest to an input graph for DFVS
                        DFVSInstance newInstance = afToGraph(newBaumComponents, treeClones, false);
                        int[][] newmatrix = newInstance.adjacencyMatrix;

                        if(isAcyclic(newmatrix)) {
                            System.out.println("Merging component " + (i+1) + " and " + (j+1) + "...");
                            aafComponents = newComponents;
                            didone = true;
                            j--;
                        }
                    }
                }
                success = false;
            }
            if(didone) {
                System.out.println("No further components can be merged.");
            } else {
                System.out.println("No components can be merged.");
            }
        }
  
        System.out.println("Acyclic agreement forest (as a partition):");
        for (int c = 0; c < aafComponents.size(); c++) {
            System.out.print("Component " + (c + 1) + ": ");
            int[] component = aafComponents.get(c);
            for (int i = 0; i < component.length; i++) {
                System.out.print(component[i] + " ");
            }
            System.out.print(EOL);
        }
        System.out.println("Acyclic agreement forest (as a partition, with the original taxon tabels):");
        for (int c = 0; c < aafComponents.size(); c++) {
            System.out.print("Component " + (c + 1) + ": ");
            int[] component = aafComponents.get(c);
            for (int i = 0; i < component.length; i++) {
                System.out.print(TAXA.get(component[i]) + " ");
            }
            System.out.print(EOL);
        }

        // write AAF to file
        // use the original taxa names
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(aafFile));
            for (int[] component : aafComponents) {
                for (int i = 0; i < component.length; i++) {
                    if (i > 0) {
                        out.write(" ");
                    }
                    out.write(TAXA.get(component[i]));
                }
                out.write(EOL);
            }
            out.close();
        } catch (IOException e) {
            System.out.println("Error writing to file: " + aafFile);
            return;
        }

        time = System.currentTimeMillis() - startingTime;
        System.out.println("Computation time so far: " + time / 1000 + " seconds.");
        
        System.out.println("Constructing network...");

        Network N = aaf2network(afBaumComponents, trees);
        System.out.println("Network found:");
        String networkString = N.toString();
        System.out.println(networkString);

        System.out.println("Checking if network displays the trees...");
        if (N.displays(trees)) {
            System.out.println("Network displays both input trees!");
        } else {
            System.out.println("Error: network does NOT display the input trees!");
        }

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
            out.write(networkString);
            out.close();
        } catch (IOException e) {
            System.out.println("Error writing to file: " + outFile);
            return;
        }
        System.out.println("----------------------------------------------------------------- ");
        System.out.println("Acyclic agreement forest saved to " + aafFile);
        System.out.println("Number of components: " + aafComponents.size());
        System.out.println("Network saved to " + outFile);
        System.out.println("Reticulation number: " + (aafComponents.size() - 1));
        time = System.currentTimeMillis() - startingTime;
        System.out.println("Total computation time: " + time / 1000 + " seconds.");
    }
    
    static class DFVSInstance {
        List<Baum> intAFVertices;
        List<Baum> nonRootAFVertices;
        int[][] adjacencyMatrix;
        int[] weights;
        
        public DFVSInstance(List<Baum> i, List<Baum> n, int[][] a, int[] w) {
            intAFVertices = i;
            nonRootAFVertices = n;
            adjacencyMatrix = a;
            weights = w;
        }
    }
    
    private static DFVSInstance afToGraph(List<Baum> afBaumComponents, List<Baum> trees, boolean talk) {
        
        // we need to map each vertex of each component to a vertex of each tree
        // and we need to map each edge of each component to a vertex of each tree (the head of the first edge it corresponds to)
       
        // get the internal vertices of the AF
        List<Baum> intAFVertices = new LinkedList();
        for (Baum tree : afBaumComponents) {
            intAFVertices.addAll(tree.getInternalVertices());
        }
        int n1 = intAFVertices.size();
        if(talk) {
            System.out.println("The agreement forest has " + n1 + " internal vertices.");
        }
        // map each internal vertex of the AF to a vertex of each tree
        for (Baum v : intAFVertices) {
            for (Baum tree : trees) {
                v.vertexMapsTo.add(tree.getLCA(v.getCluster()));
            }
        }

        // we identify each edge with its head, which is a non-root vertex
        List<Baum> nonRootAFVertices = new LinkedList();
        for (Baum tree : afBaumComponents) {
            nonRootAFVertices.addAll(tree.getNonRootVertices());
        }
        int n2 = nonRootAFVertices.size();
        if(talk) {
            System.out.println("The agreement forest has " + n2 + " edges.");
        }
        int n = n1 + n2;
        if(talk) {
            System.out.println("The input graph for DFVS gets " + n + " vertices.");
        }

        // map each non-root vertex (i.e. edge) of the AF to a vertex of each tree
        for (Baum v : nonRootAFVertices) {
            Baum u = v.parent;
            for (Baum tree : trees) {
                Baum uMapsTo = tree.getLCA(u.getCluster());
                List<Integer> vcluster = v.getCluster();
                for (Baum child : uMapsTo.children) {
                    if (child.hasDecendantFrom(vcluster)) {
                        v.edgeMapsTo.add(child);
                    }
                }
            }
        }

        // now we construct the input graph for DFVS
        int[][] adjacencyMatrix = new int[n][n];
        // there is an arc from v to e if v is the tail of e
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                if (intAFVertices.get(i) == nonRootAFVertices.get(j).parent) {
                    adjacencyMatrix[i][n1 + j] = 1;
                    adjacencyMatrix[n1 + j][i] = -1;
                }
            }
        }

        // there is an arc from e to v if v is reachable from the head of e in either tree
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                for (int k = 0; k < trees.size(); k++) {
                    Baum v = intAFVertices.get(i).vertexMapsTo.get(k);
                    Baum e = nonRootAFVertices.get(j).edgeMapsTo.get(k);
                    if (e.hasDecendant(v)) {
                        adjacencyMatrix[i][n1 + j] = -1;
                        adjacencyMatrix[n1 + j][i] = 1;
                    }
                }
            }
        }

        // now we need to define the weights of the vertices of the graph
        int[] weights = new int[n];
        for (int i = 0; i < n1; i++) {
            weights[i] = intAFVertices.get(i).children.size() - 1;
        }
        for (int j = 0; j < n2; j++) {
            weights[n1 + j] = 1;
        }
        
        return new DFVSInstance(intAFVertices,nonRootAFVertices,adjacencyMatrix,weights);
    }
    
    private static List<Baum> getBaumComponents(List<int[]> afComponents, List<Baum> trees) {
       
        List<Baum> afBaumComponents = new LinkedList();

        // label the edges of the trees by the AF components
        List<Integer> allNumbers = new LinkedList();
        for (int[] component : afComponents) {
            for(Integer i : component) {
                if(allNumbers.contains(i)) {
                    // System.out.println("Error: incorrect agreement forest!");
                    return null;
                } else {
                    allNumbers.add(i);
                }
            }
        }
        if(allNumbers.size() != TAXA.size()) {
            // System.out.println("Error: incorrect agreement forest!");
            return null;
        }
        for (int[] component : afComponents) {
            for (Baum tree : trees) {
                boolean suc = tree.addComponent(component, afComponents.indexOf(component));
                if(!suc) {
                    // System.out.println("Error: incorrect agreement forest!");
                    return null;
                }
            }
        }
        // System.out.println("Agreement forest is correct!");
        
        if (DEBUG) {
            for (Baum tree : trees) {
                String tFile = "T" + (trees.indexOf(tree) + 1) + ".dot";
                tree.toFile(tFile, afComponents.size());
            }
        }

        // construct a tree for each component of the AF
        for (int k = 0; k < afComponents.size(); k++) {
            List<List<Integer>> clusters = new LinkedList();
            if (afComponents.get(k).length == 1) {
                List<Integer> cluster = new LinkedList();
                cluster.add(afComponents.get(k)[0]);
                clusters.add(cluster);
            } else {
                for (Baum tree : trees) {
                    List<List<Integer>> newClusters = tree.getClustersOfComponent(k);
                    for (List<Integer> cluster : newClusters) {
                        List<Integer> sortedCluster = Baum.sortCluster(cluster);
                        if (!clusters.contains(sortedCluster)) {
                            clusters.add(sortedCluster);
                        }
                    }
                }
            }
            afBaumComponents.add(Baum.clustersToBaum(clusters));
        }
        
        // here are the components that we found
//        System.out.println("The trees for each component of the agreement forest:");
//        for (Baum tree : afBaumComponents) {
//            System.out.println(tree.toString());
//        }
        // give the leaves correct labels
        for (Baum tree : afBaumComponents) {
            tree.setLabels();
        }
        
        return afBaumComponents;
    }

    public static Network aaf2network(List<Baum> aafComponents, List<Baum> trees) {
        List<Network> components = new LinkedList();
        for (Baum baum : aafComponents) {
            Network component = baum.toNetwork();
            // map each vertex of the AAF to a vertex of each tree
            for (Network v : component.getVertices()) {
                for (Baum tree : trees) {
                    v.treeVertices.add(tree.getLCA(v.getCluster()));
                }
            }
            components.add(component);
        }

        // construct the inheritance graph
        int n = components.size();
        int[][] IG = new int[n][n];
        for (int i = 0; i < n; i++) {
            Network component = components.get(i);
            List<Integer> cluster = component.getCluster();
            Baum v1 = component.treeVertices.get(0);
            Baum v2 = component.treeVertices.get(1);
            for (int j = 0; j < n; j++) {
                if (j == i) {
                    continue;
                }
                Network component2 = components.get(j);
                Baum w1 = component2.treeVertices.get(0);
                Baum w2 = component2.treeVertices.get(1);
                for (Baum c1 : v1.children) {
                    if (!c1.hasDecendantFrom(cluster)) {
                        continue;
                    }
                    if (c1.hasDecendant(w1)) {
                        IG[i][j] = 1;
                        IG[j][i] = -1;
                    }
                }
                for (Baum c2 : v2.children) {
                    if (!c2.hasDecendantFrom(cluster)) {
                        continue;
                    }
                    if (c2.hasDecendant(w2)) {
                        IG[i][j] = 1;
                        IG[j][i] = -1;
                    }
                }
            }
        }
        if (DEBUG) {
            // print the inheritance graph
            System.out.println("Adjacency matrix of the inheritance graph:");
            for (int i = 0; i < n; i++) {
                System.out.print("[");
                for (int j = 0; j < n; j++) {
                    System.out.print(IG[i][j]);
                    if (j < n - 1) {
                        System.out.print(" ");
                    }
                }
                System.out.print("]" + EOL);
            }
        }
        // find a topological ordering of the vertices of the inheritance graph
        List<Integer> topOrd = getTopOrd(IG);
        if (DEBUG) {
            System.out.println("Topological ordering: " + topOrd.toString());
        }
        // contruct an initial network
        Network N = components.get(topOrd.get(0));
        // add a fake root
        Network root = new Network();
        root.children.add(N);
        N.parents.add(root);
        N = root;
        List<Integer> networkTaxa;
        for (int i = 1; i < n; i++) {
            networkTaxa = N.getCluster();
            Network component = components.get(topOrd.get(i));
            // hang component in network N
            // the corresponding vertex in the first tree
            Baum v1 = component.treeVertices.get(0);
            // search upward for a vertex whose cluster contains network taxa
            v1 = v1.searchUpwards(networkTaxa);
            // and in the second tree
            Baum v2 = component.treeVertices.get(1);
            v2 = v2.searchUpwards(networkTaxa);
            // we need to find v1 and v2 in network N
            // first v1
            List<Integer> cluster1 = v1.getCluster();
            cluster1.retainAll(networkTaxa);
            Network u1 = N.children.get(0).findVertexWithCluster(cluster1, 0);
            // now v2
            List<Integer> cluster2 = v2.getCluster();
            cluster2.retainAll(networkTaxa);
            Network u2 = N.children.get(0).findVertexWithCluster(cluster2, 1);
            // subdivide the edge entering u1 by w1
            Network w1 = new Network();
            Network p1 = u1.parents.get(0);
            u1.parents.set(0, w1);
            p1.children.set(p1.children.indexOf(u1), w1);
            w1.parents.add(p1);
            w1.children.add(u1);
            // subdivide the edge entering u2 by w2
            Network w2 = new Network();
            Network p2 = u2.parents.get(0);
            u2.parents.set(0, w2);
            p2.children.set(p2.children.indexOf(u2), w2);
            w2.parents.add(p2);
            w2.children.add(u2);
            // hang the component below a reticulation below w1 and w2
            Network r = new Network();
            r.parents.add(w1);
            r.parents.add(w2);
            w1.children.add(r);
            w2.children.add(r);
            r.children.add(component);
            component.parents.add(r);
        }
        // check if we can safely remove the fake root
        if (N.children.size() == 1) {
            boolean remove = true;
            Network child = N.children.get(0);
            for (Network grandchild : child.children) {
                if (grandchild.parents.size() > 1) {
                    remove = false;
                }
            }
            if (remove) {
                // remove the fake root
                N = N.children.get(0);
                N.parents.clear();
            }
        }
        return N;
    }

    public static List<Integer> getTopOrd(int[][] A) {
        List<Integer> topOrd = new LinkedList();
        int n = A.length;
        boolean[] removed = new boolean[n];
        while (topOrd.size() < n) {
            for (int i = 0; i < n; i++) {
                if (removed[i]) {
                    continue;
                }
                boolean source = true;
                for (int j = 0; j < n; j++) {
                    if (removed[j]) {
                        continue;
                    }
                    if (A[i][j] == -1) {
                        source = false;
                        break;
                    }
                }
                if (source) {
                    topOrd.add(i);
                    removed[i] = true;
                    break;
                }
            }
        }
        return topOrd;
    }

    public static String clusterToString(List<Integer> cluster) {
        String out = "[";
        for (int x : cluster) {
            out += TAXA.get(x) + ", ";
        }
        out = out.substring(0, out.length() - 2) + "]";
        return out;
    }

    public static List<Integer> solveDFVS(int[][] adjacencyMatrix, int[] weights) {
        int n = weights.length;
        // make outlists and inlists
        List<List<Integer>> outlists = new LinkedList();
        List<List<Integer>> inlists = new LinkedList();
        for (int i = 0; i < n; i++) {
            List<Integer> outlist = new LinkedList();
            List<Integer> inlist = new LinkedList();
            for (int j = 0; j < n; j++) {
                if (adjacencyMatrix[i][j] == 1) {
                    outlist.add(j);
                }
                if (adjacencyMatrix[i][j] == -1) {
                    inlist.add(j);
                }
            }
            outlists.add(outlist);
            inlists.add(inlist);
        }

        if (DEBUG) {
            System.out.println("Searching for directed cycles...");
        }
        List<List<Integer>> dirCycles = getDirectedCycles(outlists, inlists, new LinkedList());
        if (dirCycles.isEmpty()) {
            if (DEBUG) {
                System.out.println("No directed cycles!");
            }
            return new LinkedList();
        } else {
            if (DEBUG) {
                for (List<Integer> dirCycle : dirCycles) {
                    System.out.println("Directed cycle: " + dirCycle.toString());
                }
            }
        }

        return solveDFVS(outlists, inlists, weights, dirCycles);
    }
    
    public static boolean isAcyclic(int[][] adjacencyMatrix) {
        int n = adjacencyMatrix.length;
        // make outlists and inlists
        List<List<Integer>> outlists = new LinkedList();
        List<List<Integer>> inlists = new LinkedList();
        for (int i = 0; i < n; i++) {
            List<Integer> outlist = new LinkedList();
            List<Integer> inlist = new LinkedList();
            for (int j = 0; j < n; j++) {
                if (adjacencyMatrix[i][j] == 1) {
                    outlist.add(j);
                }
                if (adjacencyMatrix[i][j] == -1) {
                    inlist.add(j);
                }
            }
            outlists.add(outlist);
            inlists.add(inlist);
        }

        List<List<Integer>> dirCycles = getDirectedCycles(outlists, inlists, new LinkedList());
        if (dirCycles.isEmpty()) {
            return true;
        } else {
            return false;
        }
    }

    public static List<Integer> solveDFVS(List<List<Integer>> outlists, List<List<Integer>> inlists, int[] weights, List<List<Integer>> directedCycles) {
        int n = weights.length;

        while (true) {

            List<Integer> DFVS = new LinkedList();
            // make the ILP formulation
            List<String> ILP = new LinkedList();
            // first the objective function
            ILP.add("Minimize");
            StringBuffer obj = new StringBuffer();
            for (int i = 0; i < n; i++) {
                obj.append(" + ").append(weights[i]).append(" x_").append(i);
            }
            ILP.add(obj.substring(3));

            // now the constraints
            ILP.add("Subject To");
            for (List<Integer> directedCycle : directedCycles) {
                obj = new StringBuffer();
                for (int v : directedCycle) {
                    obj.append(" + x_").append(v);
                }
                obj.append(" >= 1");
                ILP.add(obj.substring(3));
            }


//            ILP.add("Bounds");
//            for (int i = 0; i < n; i++) {
//                ILP.add("0 <= x_" + i + " <= 1");
//            }

            // integrality constraints
            ILP.add("Binary");
            for (int i = 0; i < n; i++) {
                ILP.add(" x_" + i);
            }

            ILP.add("End");

            // write ILP formulation to file
            try {
                BufferedWriter out = new BufferedWriter(new FileWriter("ilp.tmp"));
                for (String s : ILP) {
                    out.write(s + EOL);
                }
                out.close();
            } catch (IOException e) {
                return null;
            }

            if (DEBUG) {
                System.out.println("Solving ILP with CPLEX...");
            }

            // solve the ILP by CPLEX
            try {

                IloCplex cplex = new IloCplex();

                // set the maximum number of threads to use
                if(ONECORE) {
                    cplex.setParam(IloCplex.IntParam.Threads, 1);
                }

                // set the time limit
                // cplex.setParam(IloCplex.DoubleParam.TiLim, 120);

                //! filename is the name of the file where your ILP is
                cplex.importModel("ilp.tmp");

                // read initial solution
                // cplex.readMIPStart("mipstart.tmp");

                //! uncomment this to suppress visual output from cplex
                cplex.setOut(null);

                //! this is the solving bit
                cplex.solve();

                IloNumVar[] var = parse(cplex);

                //! read the optimal solution
                double[] x = cplex.getValues(var);

                for (int loop = 0; loop < x.length; loop++) {
                    String varname = var[loop].getName();
                    String[] splitVarName = varname.split("_");
                    if (splitVarName[0].equals("x")) {

                        if (x[loop] < 0.99) {
                            // this variable is set to 0. Skip!
                            continue;
                        }

                        int v = Integer.parseInt(splitVarName[1]);
                        DFVS.add(v);
                    }
                }

                //! this gets the objective function value, rounded to an int
                int opt = (int) Math.round(cplex.getObjValue());

                // System.out.println("Optimal solution to ILP has value: " + opt);

                // save solution for future reference
                // cplex.writeMIPStart("mipstart.tmp");

                //! this deallocates the CPLEX resources
                cplex.end();

            } catch (IloException e) {
                System.out.println("Something went wrong with CPLEX.");
                System.out.print(e.getMessage());
                System.exit(0);
            }

            if (DEBUG) {
                System.out.println("Candidate for DFVS: " + DFVS.toString());
            }

            if (DEBUG) {
                System.out.println("Searching for directed cycles...");
            }
            List<List<Integer>> dirCycles = getDirectedCycles(outlists, inlists, DFVS);
            if (dirCycles.isEmpty()) {
                System.out.println("No directed cycles left!");
                return DFVS;
            }

            if (DEBUG) {
                for (List<Integer> dirCycle : dirCycles) {
                    System.out.println("Directed cycle: " + dirCycle.toString());
                }
            }
            directedCycles.addAll(dirCycles);
        }
    }

    public static List<List<Integer>> getDirectedCycles(List<List<Integer>> out, List<List<Integer>> in, List<Integer> U) {
        List<List<Integer>> directedCycles = new LinkedList();
        int n = out.size();
        // first construct the set of vertices
        List<Integer> V = new LinkedList();
        for (int i = 0; i < n; i++) {
            V.add(i);
        }
        V.removeAll(U);
        // now copy the outlists and inlists
        List<List<Integer>> outlists = new LinkedList();
        List<List<Integer>> inlists = new LinkedList();
        for (int i = 0; i < n; i++) {
            List<Integer> outlist = new LinkedList();
            List<Integer> inlist = new LinkedList();
            outlist.addAll(out.get(i));
            inlist.addAll(in.get(i));
            outlist.removeAll(U);
            inlist.removeAll(U);
            outlists.add(outlist);
            inlists.add(inlist);
        }
        // repeatedly remove all indegree-0 and outdegree-0 vertices
        boolean success;
        while (V.size() > 1) {
            success = false;
            for (int i = 0; i < V.size(); i++) {
                Integer v = V.get(i);
                if (outlists.get(v).isEmpty() | inlists.get(v).isEmpty()) {
                    //System.out.println("Removing vertex " + v);
                    for (List<Integer> outlist : outlists) {
                        outlist.remove(v);
                    }
                    for (List<Integer> inlist : inlists) {
                        inlist.remove(v);
                    }
                    V.remove(v);
                    i--;
                    success = true;
                }
            }
            if (!success & V.size() > 1) {
                // there is at least one directed cycle
                // moreover, each vertex is in at least one directed cycle
                List<Integer> directedCycle = new LinkedList();
                Integer v = V.get(0);
                Integer w = -1;
                // start searching for a directed cycle at v
                Stack<Integer> stack = new Stack();
                int[] prev = new int[n];
                stack.push(v);
                while (!stack.isEmpty()) {
                    w = stack.get(0);
                    stack.remove(0);
                    for (Integer y : outlists.get(w)) {
                        if (prev[y] == 0) {
                            prev[y] = w;
                            stack.push(y);
                        }
                        if (y.equals(v)) {
                            directedCycle.add(w);
                            directedCycle.add(v);
                            stack.clear();
                            break;
                        }
                    }
                }
                // backtracking
                w = directedCycle.get(0);
                while (prev[w] != v.intValue()) {
                    w = prev[w];
                    directedCycle.add(0, w);
                }
                directedCycles.add(directedCycle);
                //System.out.println("Directed cycle found: " + directedCycle.toString());

                // remove all vertices of the directed cycle
                V.removeAll(directedCycle);
                for (List<Integer> outlist : outlists) {
                    outlist.removeAll(directedCycle);
                }
                for (List<Integer> inlist : inlists) {
                    inlist.removeAll(directedCycle);
                }
            }
        }

        return directedCycles;
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

class Baum {

    List<Baum> children;
    Baum parent;
    String label;
    int number;
    int printnum;
    List<Integer> components;
    List<Baum> vertexMapsTo;
    List<Baum> edgeMapsTo;
    
    @Override
    public Baum clone() {
        Baum clone = new Baum();
        clone.number = number;
        clone.printnum = printnum;
        clone.label = label;
        for(Baum child : children) {
            Baum childClone = child.clone();
            clone.children.add(childClone);
            childClone.parent = clone;
        }
        return clone;
    }

    public Baum() {
        children = new LinkedList();
        parent = null;
        label = null;
        number = -1;
        printnum = -1;
        components = new LinkedList();
        vertexMapsTo = new LinkedList();
        edgeMapsTo = new LinkedList();
    }

    public Baum(String l) {
        children = new LinkedList();
        parent = null;
        label = l;
        number = -1;
        printnum = -1;
        components = null;
        vertexMapsTo = new LinkedList();
        edgeMapsTo = new LinkedList();
    }

    public Baum(int k) {
        children = new LinkedList();
        parent = null;
        label = null;
        number = k;
        printnum = -1;
        components = null;
        vertexMapsTo = new LinkedList();
        edgeMapsTo = new LinkedList();
    }

    public boolean isBinary() {
        if (children.size() == 1 | children.size() > 2) {
            return false;
        }
        for (Baum child : children) {
            if (!child.isBinary()) {
                return false;
            }
        }
        return true;
    }
    
    public Baum searchUpwards(List<Integer> taxa) {
        List<Integer> cluster = this.getCluster();
        for (Integer x : cluster) {
            if (taxa.contains(x)) {
                return this;
            }
        }
        return parent.searchUpwards(taxa);
    }

    public List<List<Integer>> getClusters() {
        List<List<Integer>> clusters = new LinkedList();
        clusters.add(this.getCluster());
        for (Baum child : children) {
            clusters.addAll(child.getClusters());
        }
        return clusters;
    }

    public void setLabels() {
        if (children.isEmpty()) {
            label = NonbinaryCycleKiller.TAXA.get(number);
            return;
        }
        for (Baum child : children) {
            child.setLabels();
        }
    }

    public Network toNetwork() {
        if (children.isEmpty()) {
            return new Network(label, number);
        }
        Network N = new Network();
        for (Baum child : children) {
            Network networkChild = child.toNetwork();
            N.children.add(networkChild);
            networkChild.parents.add(N);
        }
        return N;
    }

    public List<String> toDot(int s) {
        List<String> out = new LinkedList();
        out.add("strict digraph G {");
        int[] num = new int[1];
        num[0] = 1000;
        out.addAll(this.nodes2dot(num, s));
        out.addAll(this.arcs2dot(s));
        out.add("}");
        clearPrintNum();
        return out;
    }
    
    public void clearPrintNum(){
        printnum = -1;
        for(Baum child : children) {
            child.clearPrintNum();
        }
    }
    

    public List<String> nodes2dot(int num[], int s) {
        List<String> out = new LinkedList();
        if (children.isEmpty()) {
            // this is a leaf
            printnum = num[0];
            out.add(printnum + " [shape=none, label=\"" + label + "\"];");
        } else {
            printnum = num[0];
            out.add(printnum + " [shape=point, width=0, color=black];");
            for (Baum child : children) {
                num[0]++;
                out.addAll(child.nodes2dot(num, s));
            }
        }
        return out;
    }

    public List<String> arcs2dot(int s) {
        // s is the number of components
        // if there are more than 12 components we do not have enough colours to colour them
        List<String> out = new LinkedList();
        for (Baum child : children) {
            int c = children.indexOf(child);
            if (components.get(c) == -1) {
                out.add(printnum + " -> " + child.printnum + "[color=grey]");
            } else {
                if (s <= 12) {
                    out.add(printnum + " -> " + child.printnum + "[colorscheme=paired12, color=" + (components.get(c)+1) + ", label=\"" + (components.get(c)+1) + "\"]");
                } else {
                    out.add(printnum + " -> " + child.printnum + "[label=\"" + (components.get(c)+1) + "\"]");
                }
            }
        }
        for (Baum child : children) {
            out.addAll(child.arcs2dot(s));
        }
        return out;
    }
    
    public static Baum clustersToBaum(List<List<Integer>> clusters) {

        if (clusters.size() == 1) {
            List<Integer> cluster = clusters.get(0);
            if (cluster.size() == 1) {
                return new Baum(cluster.get(0).intValue());
            }
        }

        Baum tree = new Baum();

        // find the maximal clusters
        List<List<Integer>> maxClusters = new LinkedList();
        for (List<Integer> cluster : clusters) {
            boolean max = true;
            for (List<Integer> cluster2 : clusters) {
                if (!cluster.equals(cluster2) && isContainedIn(cluster, cluster2)) {
                    max = false;
                }
            }
            if (max) {
                maxClusters.add(cluster);
            }
        }

        // the maximal clusters become children of the root
        for (List<Integer> cluster : maxClusters) {
            Baum child = clustersToBaum(restrictTo(clusters, cluster));
            tree.children.add(child);
            child.parent = tree;
        }

        return tree;
    }

    public static boolean isContainedIn(List<Integer> cluster1, List<Integer> cluster2) {
        for (int taxon : cluster1) {
            if (!cluster2.contains(taxon)) {
                return false;
            }
        }
        return true;
    }

    public static List<List<Integer>> restrictTo(List<List<Integer>> clusters, List<Integer> cluster) {
        List<List<Integer>> resClusters = new LinkedList();
        for (List<Integer> cluster2 : clusters) {
            if (isContainedIn(cluster2, cluster)) {
                if (!cluster2.equals(cluster) | cluster2.size() == 1) {
                    resClusters.add(cluster2);
                }
            }
        }
        return resClusters;
    }

    public List<Baum> getLeaves() {
        List<Baum> leaves = new LinkedList();
        if (children.isEmpty()) {
            leaves.add(this);
            return leaves;
        }
        for (Baum child : children) {
            leaves.addAll(child.getLeaves());
        }
        return leaves;
    }

    public List<List<Integer>> getClustersOfComponent(int k) {
        List<List<Integer>> clusters = new LinkedList();
        if (children.isEmpty()) {
            return clusters;
        }
        for (Baum child : children) {
            int i = children.indexOf(child);
            if (components.get(i).intValue() == k) {
                List<Integer> cluster = child.getComponentCluster(k);
                if (!clusters.contains(cluster)) {
                    clusters.add(cluster);
                }
            }
            for (List<Integer> cluster : child.getClustersOfComponent(k)) {
                if (!clusters.contains(cluster)) {
                    clusters.add(cluster);
                }
            }
        }
        return clusters;
    }

    public List<Integer> getComponentCluster(int k) {
        List<Integer> cluster = new LinkedList();
        if (children.isEmpty()) {
            cluster.add(number);
            return cluster;
        }
        for (Baum child : children) {
            int i = children.indexOf(child);
            if (components.get(i).intValue() == k) {
                cluster.addAll(child.getComponentCluster(k));
            }
        }
        return sortCluster(cluster);
    }

    public List<Baum> getInternalVertices() {
        List<Baum> vertices = new LinkedList();
        if (children.isEmpty()) {
            return vertices;
        }
        vertices.add(this);
        for (Baum child : children) {
            vertices.addAll(child.getInternalVertices());
        }
        return vertices;
    }

    public List<Baum> getNonRootVertices() {
        List<Baum> vertices = new LinkedList();
        if (parent != null) {
            vertices.add(this);
        }
        for (Baum child : children) {
            vertices.addAll(child.getNonRootVertices());
        }
        return vertices;
    }

    public boolean isLCA(List<Integer> taxa) {
        List<Integer> cluster = getCluster();
        if (!isContainedIn(taxa, cluster)) {
            return false;
        }
        for (Baum child : children) {
            List<Integer> childCluster = child.getCluster();
            if (isContainedIn(taxa, childCluster)) {
                return false;
            }
        }
        return true;
    }

    public List<Integer> getCluster() {
        List<Integer> cluster = new LinkedList();
        if (children.isEmpty()) {
            cluster.add(number);
            return cluster;
        }
        for (Baum child : children) {
            cluster.addAll(child.getCluster());
        }
        return sortCluster(cluster);
    }

    public Baum getLCA(List<Integer> taxa) {
        if (isLCA(taxa)) {
            return this;
        }
        for (Baum child : children) {
            Baum rec = child.getLCA(taxa);
            if (rec != null) {
                return rec;
            }
        }
        return null;
    }

    public boolean hasDecendantFrom(List<Integer> taxa) {
        List<Integer> cluster = getCluster();
        for (int x : cluster) {
            if (taxa.contains(x)) {
                return true;
            }
        }
        return false;
    }

    public boolean hasDecendant(Baum v) {
        if (this.equals(v)) {
            return true;
        }
        for (Baum child : children) {
            if (child.hasDecendant(v)) {
                return true;
            }
        }
        return false;
    }
    
    public void relabel(String label1, String label2) {
        if(label != null && label.equals(label1)) {
            label = label2;
        }
        for(Baum child : children) {
            child.relabel(label1,label2);
        }
    }
    
    public void removeOutdegreeZero() {
        for(int c = 0; c < children.size(); c++) {
            List<Integer> cluster = children.get(c).getCluster();
            boolean phantom = true;
            for(Integer x : cluster) {
                if(!x.equals(-1)) {
                    phantom = false;
                }
            }
            if(phantom){
                children.remove(c);
                c--;
            }
        }
        for(Baum child : children) {
            child.removeOutdegreeZero();
        }
    }
    
    public Baum removeFakeRoots() {
        if(children.size() == 1) {
            children.get(0).parent = null;
            return children.get(0).removeFakeRoots();
        } else {
            return this;
        }
    }
    
    public void toFile(String filename, int s) {
        List<String> dot = this.toDot(s);
        // write to file
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(filename));
            for (String str : dot) {
                out.write(str + "\n");
            }
            out.close();
        } catch (IOException e) {
            return;
        }
        // transform dot to PDF
        try {
            String line;
            Process p = Runtime.getRuntime().exec("dot -Tpdf " + filename + " -O ");
            BufferedReader bre = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            while ((line = bre.readLine()) != null) {
                System.out.println(line);
            }
            bre.close();
            p.waitFor();
            //System.out.println("Converted to PDF in " + filename + ".pdf");
        } catch (Exception err) {
        }
    }

    public boolean addComponent(int[] component, int k) {
        if (children.isEmpty()) {
            return true;
        }
        while (components.size() < children.size()) {
            components.add(-1);
        }
        for (Baum child : children) {
            // this outgoing edge belongs to the component if at least one, but not all, of the taxa in the component can be reached from this edge
            int below = 0;
            for (int taxon : component) {
                if (child.isAncestorOf(taxon)) {
                    below++;
                }
            }
            if (below > 0 && below < component.length) {
                if(components.get(children.indexOf(child)) != -1) {
                    return false;
                }
                components.set(children.indexOf(child), k);
            }
            // recurse
            boolean childsuc = child.addComponent(component, k);
            if(!childsuc) {
                return false;
            }
        }
        return true;
    }

    public boolean isAncestorOf(int taxon) {
        if (children.isEmpty() && number == taxon) {
            return true;
        }
        for (Baum child : children) {
            if (child.isAncestorOf(taxon)) {
                return true;
            }
        }
        return false;
    }

    public static Baum newickToBaum(String newick) {
        Baum tree;
        int lastclosepar = newick.lastIndexOf(")");
        int lastcolon = newick.lastIndexOf(":");

        // get rid of semicolon
        if (newick.endsWith(";")) {
            newick = newick.substring(0, newick.length() - 1);
        }

        // get the edge length
        Double length = 0.0;
        if (lastcolon != -1 && (lastcolon > lastclosepar | lastclosepar == -1)) {
            length = Double.parseDouble(newick.substring(lastcolon + 1, newick.length()));
            newick = newick.substring(0, lastcolon);
        }

        if (newick.startsWith("(")) {
            // internal vertex
            tree = new Baum();
            // remove everything after the last closing parenthesis (bootstrap value)
            newick = newick.substring(0, lastclosepar + 1);
            int openpar = 0;
            int closepar = 0;
            int start = 1;
            List<String> childrenNewick = new LinkedList();
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
                Baum child = newickToBaum(childNewick);
                tree.children.add(child);
                child.parent = tree;
            }
        } else {
            // a leaf
            tree = new Baum(newick);
        }
        return tree;
    }

    public boolean equals(Baum v) {
        if (this == v) {
            return true;
        }
        if (this.label == null | v.label == null) {
            return false;
        }

        if (this.label.equals(v.label)) {
            return true;
        }
        return false;
    }

    public static List<Integer> sortCluster(List<Integer> cluster) {
        int a = 0;
        for (int i = 0; i < cluster.size() - 1; i++) {
            for (int j = i + 1; j < cluster.size(); j++) {
                if (cluster.get(i) > cluster.get(j)) {
                    // swap
                    a = cluster.get(i);
                    cluster.set(i, cluster.get(j));
                    cluster.set(j, a);
                }
            }
        }
        return cluster;
    }

    public static int[] sortCluster(int[] cluster) {
        int a = 0;
        for (int i = 0; i < cluster.length - 1; i++) {
            for (int j = i + 1; j < cluster.length; j++) {
                if (cluster[i] > cluster[j]) {
                    // swap
                    a = cluster[i];
                    cluster[i] = cluster[j];
                    cluster[j] = a;
                }
            }
        }
        return cluster;
    }

    public String toStringLabels() {
        if (children.isEmpty()) {
            return "" + label;
        }
        String newick = "(";
        for (Baum child : children) {
            newick += child.toStringLabels();
            newick += ",";
        }
        // remove the last comma
        newick = newick.substring(0, newick.length() - 1);
        newick += ")";
        if (parent == null) {
            // add a semicolon
            newick += ";";
        }
        return newick;
    }
    
    @Override
    public String toString() {
        if (children.isEmpty()) {
            return "" + number;
        }
        String newick = "(";
        for (Baum child : children) {
            newick += child.toString();
            newick += ",";
        }
        // remove the last comma
        newick = newick.substring(0, newick.length() - 1);
        newick += ")";
        if (parent == null) {
            // add a semicolon
            newick += ";";
        }
        return newick;
    }
}

class Network {

    List<Network> children;
    List<Network> parents;
    String label;
    static int RET = 0;
    int ret;
    int number;
    boolean seen;
    List<Baum> treeVertices;

    public Network() {
        children = new LinkedList();
        parents = new LinkedList();
        treeVertices = new LinkedList();
        label = "";
        ret = -1;
        number = -1;
        seen = false;
    }

    public Network findVertexWithCluster(List<Integer> cluster, int index) {
        List<Integer> myCluster = this.getCluster(index);
        if (cluster.containsAll(myCluster) && myCluster.containsAll(cluster) && parents.size() == 1) {
            return this;
        }
        for (Network child : children) {
            Network v = child.findVertexWithCluster(cluster, index);
            if (v != null) {
                return v;
            }
        }
        return null;
    }

    public Network(String l, int i) {
        children = new LinkedList();
        parents = new LinkedList();
        treeVertices = new LinkedList();
        label = l;
        number = i;
        ret = -1;
        seen = false;
    }

    public List<List<Integer>> getClusters(int s) {
        List<List<Integer>> clusters = new LinkedList();
        clusters.add(this.getCluster(s));
        for (Network child : children) {
            if (child.parents.size() > 1 & child.parents.indexOf(this) != s) {
                continue;
            }
            List<List<Integer>> childClusters = child.getClusters(s);
            for (List<Integer> cluster : childClusters) {
                if (!clusters.contains(cluster)) {
                    clusters.add(cluster);
                }
            }
        }
        return clusters;
    }

    public List<Integer> getCluster(int s) {
        List<Integer> cluster = new LinkedList();
        if (children.isEmpty()) {
            cluster.add(this.number);
        }
        for (Network child : children) {
            if (child.parents.size() > 1 & child.parents.indexOf(this) != s) {
                continue;
            }
            cluster.addAll(child.getCluster(s));
        }
        cluster = sortCluster(cluster);
        return cluster;
    }

    public static List<Integer> sortCluster(List<Integer> cluster) {
        int a = 0;
        for (int i = 0; i < cluster.size() - 1; i++) {
            for (int j = i + 1; j < cluster.size(); j++) {
                if (cluster.get(i) > cluster.get(j)) {
                    // swap
                    a = cluster.get(i);
                    cluster.set(i, cluster.get(j));
                    cluster.set(j, a);
                }
            }
        }
        return cluster;
    }

    public boolean displays(List<Baum> trees) {
        List<List<Integer>> tree1Clusters = trees.get(0).getClusters();
        //System.out.println("Tree 1: " + tree1Clusters.toString());
        List<List<Integer>> networkTree1Clusters = this.getClusters(0);
        //System.out.println("Network tree 1: " + networkTree1Clusters.toString());
        List<List<Integer>> tree2Clusters = trees.get(1).getClusters();
        //System.out.println("Tree 2: " + tree2Clusters.toString());
        List<List<Integer>> networkTree2Clusters = this.getClusters(1);
        //System.out.println("Network tree 2: " + networkTree1Clusters.toString());
        // tree1Clusters and networkTree1Clusters should be compatible
        for (List<Integer> cluster1 : tree1Clusters) {
            for (List<Integer> cluster2 : networkTree1Clusters) {
                if (!compatible(cluster1, cluster2)) {
                    System.out.println("Incompatible clusters:");
                    System.out.println("In network: " + NonbinaryCycleKiller.clusterToString(cluster2));
                    System.out.println("In tree 1: " + NonbinaryCycleKiller.clusterToString(cluster1));
                    return false;
                }
            }
        }
        // tree2Clusters and networkTree2Clusters should be compatible
        for (List<Integer> cluster1 : tree2Clusters) {
            for (List<Integer> cluster2 : networkTree2Clusters) {
                if (!compatible(cluster1, cluster2)) {
                    System.out.println("Incompatible clusters:");
                    System.out.println("In network: " + NonbinaryCycleKiller.clusterToString(cluster2));
                    System.out.println("In tree 2: " + NonbinaryCycleKiller.clusterToString(cluster1));
                    return false;
                }
            }
        }
        return true;
    }

    public static boolean compatible(List<Integer> cluster1, List<Integer> cluster2) {
        if (disjoint(cluster1, cluster2)) {
            return true;
        }
        if (cluster2.containsAll(cluster1)) {
            return true;
        }
        if (cluster1.containsAll(cluster2)) {
            return true;
        }
        return false;
    }

    public static boolean disjoint(List<Integer> cluster1, List<Integer> cluster2) {
        for (Integer x : cluster1) {
            for (Integer y : cluster2) {
                if (x.equals(y)) {
                    return false;
                }
            }
        }
        return true;
    }

    public String toString() {
        String output = toStringRec();
        unseen();
        return output;
    }

    public List<Integer> getCluster() {
        List<Integer> cluster = getClusterRec();
        unseen();
        return cluster;
    }

    public List<Integer> getClusterRec() {
        List<Integer> cluster = new LinkedList();
        if (seen) {
            return cluster;
        }
        seen = true;
        if (children.isEmpty()) {
            cluster.add(this.number);
        }
        for (Network child : children) {
            cluster.addAll(child.getClusterRec());
        }
        return cluster;
    }

    public List<Network> getVertices() {
        List<Network> vertices = getVerticesRec();
        unseen();
        return vertices;
    }

    public List<Network> getVerticesRec() {
        List<Network> vertices = new LinkedList();
        if (seen) {
            return vertices;
        }
        seen = true;
        vertices.add(this);
        for (Network child : children) {
            vertices.addAll(child.getVerticesRec());
        }
        return vertices;
    }

    public void unseen() {
        seen = false;
        for (Network child : children) {
            child.unseen();
        }
    }

    public String toStringRec() {
        String output;
        // returns eNewick string of the network
        if (children.isEmpty()) {
            return label;
        }

        String childString1 = children.get(0).toStringRec();

        if (parents.size() > 1) {
            // reticulation
            if (seen) {
                output = "#H" + ret;
            } else {
                RET++;
                ret = RET;
                seen = true;
                if (childString1.startsWith("(")) {
                    output = childString1 + "#H" + ret;
                } else {
                    output = "(" + childString1 + ")" + "#H" + ret;
                }
            }
            return output;
        }

        output = "(" + childString1;

        for (int i = 1; i < children.size(); i++) {
            output += "," + children.get(i).toStringRec();
        }
        output += ")";
        if (parents.isEmpty()) {
            output += ";";
        }
        return output;
    }
}

class MAF {

    public static int UB = 0;
    public static boolean DEBUG = false;
    public static boolean DEBUGG = false;
    public static boolean SILENT = true;
    public static boolean APP_PRUNING = false;
    public static boolean CLUSTER_REDUCTION = true;

    public static void main(String[] args) {
        Tree T1;
        Tree T2;

        Long startingTime = System.currentTimeMillis();

        boolean fpt = true;
        boolean fine = true;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--nofpt")) {
                fpt = false;
            } else if (args[i].equals("--info")) {
                DEBUG = true;
            } else if (args[i].equals("--nored")) {
                CLUSTER_REDUCTION = false;
            } else if (args[i].equals("--silent")) {
                SILENT = true;
            } else if (args[i].equals("--prune")) {
                APP_PRUNING = true;
            } else if (args[i].equals("--it")) {
                if (args.length > i + 1) {
                    try {
                        UB = Integer.parseInt(args[i + 1]);
                    } catch (NumberFormatException nfe) {
                        fine = false;
                    }
                    i++;
                } else {
                    fine = false;
                }
            } else if (args[i].equals("--help")) {
                fine = false;
            } else if (i > 0) {
                System.out.println("Unknown option: " + args[i]);
                fine = false;
            }
        }

        if (args.length == 0 | !fine) {
            System.out.println("-------------------------USAGE----------------------------------------");
            System.out.println("java MAF trees.txt [options]");
            System.out.println("trees.txt should contain two trees in newick format on two lines");
            System.out.println("------------------------OPTIONS---------------------------------------");
            System.out.println("--nofpt\t apply only the approximation algorithm");
            System.out.println("--info\t show all steps of the algorithm");
            System.out.println("--it k\t search for forest of size 1,...,k before branch and bound");
            System.out.println("--prune\t use the approximation algorithm for pruning");
            System.out.println("--nored\t no cluster reductions");
            System.out.println("----------------------------------------------------------------------");
            return;
        }

        // parse the trees
        String treeFile = args[0];
        File file = new File(treeFile);
        BufferedReader reader;
        String newick1 = "";
        String newick2 = "";
        try {
            reader = new BufferedReader(new FileReader(file));
            newick1 = reader.readLine();
            newick2 = reader.readLine();
            reader.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error: file not found: " + treeFile);
            return;
        } catch (IOException e) {
        }
        T1 = Tree.newickToTree(newick1);
        T1 = T1.clean();
        T2 = Tree.newickToTree(newick2);
        T2 = T2.clean();
        System.out.println("First input tree:");
        System.out.println(T1.toString());
        //T1.toFile("T1.dot", 0);
        System.out.println("Second input tree:");
        System.out.println(T2.toString());
        //T2.toFile("T2.dot", 0);

        if (newick1.length() == 0) {
            System.out.println("Error: first input tree is empty");
            return;
        }
        if (newick2.length() == 0) {
            System.out.println("Error: second input tree is empty");
            return;
        }

        // check if both trees have the same leaves
        List<Tree> X1 = T1.getCluster();
        List<Tree> X2 = T2.getCluster();
        int n = X1.size();
        Tree T1res;
        Tree T2res;
        if (!allDifferent(X1)) {
            System.out.println("Error: first tree has duplicate leaves.");
            return;
        }
        if (!allDifferent(X2)) {
            System.out.println("Error: second tree has duplicate leaves.");
            return;
        }
        if (clustersEqual(X1, X2)) {
            System.out.println("Trees have " + n + " leaves.");
            T1res = T1;
            T2res = T2;
        } else {
            System.out.println("Error: trees do not have identical leaf sets!");

            System.out.println("Restricting to common leaves...");
            T1res = T1.clone();
            T2res = T2.clone();
            for (Tree taxon : X1) {
                if (!T2.contains(taxon)) {
                    T1res.deleteLeaf(taxon);
                }
            }
            for (Tree taxon : X2) {
                if (!T1.contains(taxon)) {
                    T2res.deleteLeaf(taxon);
                }
            }
            // clean up
            T1res = T1res.clean();
            T2res = T2res.clean();

            // check if trees now have same number of taxa
            List<Tree> X1res = T1res.getCluster();
            List<Tree> X2res = T2res.getCluster();
            n = X1res.size();
            if (n == X2res.size()) {
                System.out.println("Restricted trees have " + n + " leaves");
            } else {
                System.out.println("Error: trees still have unequal numbers of taxa");
                return;
            }
        }

        if (n == 0) {
            System.out.println("Error: no leaves");
            return;
        } else if (n == 1) {
            System.out.println("Error: only one leaf");
            return;
        }

        // run the approximation algorithm
        System.out.println("Running the approximation algorithm...");
        Tree T2app = approxMAF(T1res, T2res);
        List<List<Tree>> app_forest = T2app.getForest();
        System.out.println("Approximation algorithm finished");

        // the size of the agreement forest
        int app = app_forest.size();

        // check if the solution is correct
        System.out.println("Checking if the constructed forest is an agreement forest of the original trees...");
        boolean correct = checkForest(app_forest, T1, T2); // this also labels T1 and T2 by the component numbers

        if (!correct) {
            System.out.println("Error: constructed forest is NOT an agreement forest!");
        } else {
            System.out.println("Constructed forest is an agreement forest!");
        }
        Long approxTime = System.currentTimeMillis() - startingTime;
        //System.out.println("Computation took " + approxTime/1000 + " seconds.");
        // write agreement forest to file
        //System.out.println("Graphical representation of agreement forest saved as app_forest.dot");
        //T2app.toFile("app_forest.dot", 0);
        // write T1 with component numbers to file
        //System.out.println("T1 edge-coloured by agreement forest saved as app_T1forest.dot");
        //T1.toFile("app_T1forest.dot", app);
        // write T2 with component numbers to file
        //System.out.println("T2 edge-coloured by agreement forest saved as app_T2forest.dot");
        //T2.toFile("app_T2forest.dot", app);

        // print agreement forest to file
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter("app_maf.txt"));
            for (List<Tree> component : app_forest) {
                for (Tree x : component) {
                    out.write(x.label + " ");
                }
                out.write("\n");
            }
            out.close();
        } catch (IOException e) {
            return;
        }
        System.out.println("Agreement forest saved as app_maf.txt");
        System.out.println("Approximate agreement forest has " + app + " components");

        if (fpt) {
            System.out.println("---------------------------------------------------------");
            // remove the component labels from the trees
            T1.clearComponents();
            T2.clearComponents();

            // run the FPT algorithm
            System.out.println("Running the FPT algorithm...");
            // using upper bound given by the approximation algorithm
            // lower bound is the approximation divided by 4
            int lb = (int) Math.ceil(app * 1.0 / 4);
            Tree T2fpt;
            if (lb != app) {
                System.out.println("Optimal agreement forest has at least " + lb + " and at most " + app + " components.");
                T2fpt = fptMAF(T1res, T2res, lb, app, false);
            } else {
                System.out.println("Optimal agreement forest has " + app + " components.");
                T2fpt = T2res;
            }
            System.out.println("FPT algorithm finished.");
            List<List<Tree>> opt_forest;
            if (T2fpt != null) {
                opt_forest = T2fpt.getForest();
            } else {
                opt_forest = app_forest;
                T2fpt = T2app;
                System.out.println("Approximate agreement forest turned out to be optimal.");
            }
            int opt = opt_forest.size();

            // check if the solution is correct
            System.out.println("Checking if the constructed forest is an agreement forest of the original trees");
            correct = checkForest(opt_forest, T1, T2); // this also labels T1 and T2 by the component numbers

            if (!correct) {
                System.out.println("Error: constructed forest is NOT an agreement forest!");
            } else {
                System.out.println("Constructed forest is an agreement forest!");
            }

            System.out.println("(The approximate agreement forest had an approximation factor of " + app * 1.0 / opt + ")");

            // write agreement forest to file
            //System.out.println("Graphical representation of agreement forest saved as forest.dot");
            //T2fpt.toFile("forest.dot", 0);
            // write T1 with component numbers to file
            //System.out.println("T1 edge-coloured by agreement forest saved as T1forest.dot");
            //T1.toFile("T1forest.dot", opt);
            // write T2 with component numbers to file
            //System.out.println("T2 edge-coloured by agreement forest saved as T2forest.dot");
            //T2.toFile("T2forest.dot", opt);

            // print agreement forest to file
            try {
                BufferedWriter out = new BufferedWriter(new FileWriter("maf.txt"));
                for (List<Tree> component : opt_forest) {
                    for (Tree x : component) {
                        out.write(x.label + " ");
                    }
                    out.write("\n");
                }
                out.close();
            } catch (IOException e) {
                return;
            }
            System.out.println("Agreement forest saved as maf.txt");
            Long fptTime = System.currentTimeMillis() - startingTime;
            System.out.println("Total computation time: " + fptTime / 1000 + " seconds.");
            System.out.println("Maximum agreement forest has " + opt + " components.");
        }
    }

    public static boolean clustersEqual(List<Tree> X1, List<Tree> X2) {
        List<Tree> common = new LinkedList();
        for (Tree x1 : X1) {
            if (clusterContainsTaxon(X2, x1)) {
                common.add(x1);
            }
        }
        return (X1.size() == X2.size() & X1.size() == common.size());
    }

    public static boolean allDifferent(List<Tree> X) {
        for (Tree x : X) {
            for (Tree y : X) {
                if (X.indexOf(x) == X.indexOf(y)) {
                    continue;
                }
                if (x.equals(y)) {
                    return false;
                }
            }
        }
        return true;
    }

    public static boolean clusterContainsTaxon(List<Tree> cluster, Tree taxon) {
        for (Tree x : cluster) {
            if (x.equals(taxon)) {
                return true;
            }
        }
        return false;
    }

    private static void printSolution(Tree sol, Tree Tree1, Tree Tree2) {

        sol = sol.decollapse();
        sol = sol.clean();

        List<List<Tree>> forest = sol.getForest();
        int val = forest.size();

        boolean correct = checkForest(forest, Tree1, Tree2);

        if (correct) {
            //sol.toFile("forest.dot", 0);
            //Tree1.toFile("T1forest.dot", val);
            //Tree2.toFile("T2forest.dot", val);
            try {
                BufferedWriter out = new BufferedWriter(new FileWriter("maf.txt"));
                for (List<Tree> component : forest) {
                    for (Tree x : component) {
                        out.write(x.label + " ");
                    }
                    out.write("\n");
                }
                out.close();
            } catch (IOException e) {
                return;
            }
            System.out.println("Solution saved to file");
        } else {
            System.out.println("Error: incorrect agreement forest!");
        }
        Tree1.clearComponents();
        Tree2.clearComponents();
    }

    private static boolean checkForest(List<List<Tree>> forest, Tree Tree1, Tree Tree2) {
        // this checks if the components are edge-disjoint subtrees of T1 and T2
        // and if the subtrees are compatible
        // it also labels the each node by the component number of the incoming edge

        // first delete possible existing component numbers
        Tree1.clearComponents();
        Tree2.clearComponents();

        // first we find the edge-disjoint subtrees
        boolean out = true;
        for (List<Tree> component : forest) {
            int compnum = forest.indexOf(component) + 1;
            out = out & Tree1.colourComponent(component, compnum);
            out = out & Tree2.colourComponent(component, compnum);
        }
        if (!out) {
            return false;
        }

        if (DEBUGG) {
            Tree1.toFile("T1now.dot", 0);
            Tree2.toFile("T2now.dot", 0);
        }

        // now we need to check if the subtrees are compatible
        for (List<Tree> component : forest) {
            int compnum = forest.indexOf(component) + 1;
            List<List<Tree>> clustersetT1 = Tree1.getComponentClusters(compnum);
            List<List<Tree>> clustersetT2 = Tree2.getComponentClusters(compnum);
            // check if every pair of clusters is compatible
            for (List<Tree> c1 : clustersetT1) {
                if (c1.size() == 1) {
                    continue;
                }
                for (List<Tree> c2 : clustersetT2) {
                    if (c2.size() == 1) {
                        continue;
                    }
                    if (!compatible(c1, c2)) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    private static boolean compatible(List<Tree> cluster1, List<Tree> cluster2) {
        // checks if the given clusters are compatible
        // get the intersection
        List<Tree> intersection = new LinkedList();
        for (Tree c1 : cluster1) {
            for (Tree c2 : cluster2) {
                if (c1.equals(c2)) {
                    intersection.add(c1);
                    break;
                }
            }
        }
        if (intersection.isEmpty() | intersection.size() == cluster1.size() | intersection.size() == cluster2.size()) {
            return true;
        }
        return false;
    }

    private static Tree approxMAF(Tree Tree1, Tree Tree2) {
        boolean success = true;

        // copy trees
        Tree T1 = Tree1.clone();
        Tree T2 = Tree2.clone();

        while (success) {

            if (DEBUGG) {
                T1.toFile("T1now.dot", 0);
                T2.toFile("T2now.dot", 0);
            }

            // clean up
            T1 = T1.clean();
            T2 = T2.clean();

            if (DEBUGG) {
                T1.toFile("T1now.dot", 0);
                T2.toFile("T2now.dot", 0);
            }

            // check if finished
            if (T1.getCluster().size() < 2) {
                // T1 consists of at most one leaf, we can stop
                if (DEBUG) {
                    System.out.println("At most one leaf left, no more cuts needed...");
                }
                break;
            }

            // Step 0a: collapse all cherries
            success = collapseCherry(T1, T2);

            if (success) {
                continue;
            }

            // Step 0b: delete each singleton of T2 from T1
            success = success | deleteSingletons(T1, T2);

            if (success) {
                continue;
            }

            // Step 0c: if the leaves in C are all in different components of T2, cut them off in T2 and remove them from T1
            success = success | cutCherries(T1, T2);

            if (success) {
                continue;
            }

            // Step 1 and 2
            success = success | step1and2(T1, T2);

            if (!success & DEBUG) {
                System.out.println("No more cuts found...");
            }

        }

        // decollapse
        T2 = T2.decollapse();

        // clean up
        T2 = T2.clean();

        return T2;
    }

    static class StackObject {

        Tree T1;
        Tree T2;
        int components;

        public StackObject() {
            T1 = null;
            T2 = null;
            components = 0;
        }

        public StackObject(Tree t1, Tree t2, int c) {
            T1 = t1;
            T2 = t2;
            components = c;
        }
    }

    private static Tree fptMAF(Tree Tree1, Tree Tree2, int lb, int app, boolean recursive) {

        // TRY CLUSTER REDUCTION
        if (CLUSTER_REDUCTION) {
            // make working copies
            Tree T1copy = Tree1.clone();
            Tree T2copy = Tree2.clone();
            List<Tree> trees = T1copy.findCommonCluster(T2copy);
            if (!trees.isEmpty()) {
                // do cluster reduction
                
                // make four copies of each tree
                Tree subtree1 = trees.get(0);
                Tree s1copy1 = subtree1.clone();
                Tree s1copy2 = subtree1.clone();
                Tree subtree2 = trees.get(1);
                Tree s2copy1 = subtree2.clone();
                Tree s2copy2 = subtree2.clone();
                
                if (!SILENT) {
                    System.out.println("Applying cluster reduction to: " + s1copy1.getCluster().toString());
                }
                
                // FIRST SUBPROBLEM
                // add a dummy root
                // to first tree
                Tree dummyRoot1 = new Tree();
                Tree dummyLeaf1 = new Tree("dummy");
                dummyRoot1.children.add(dummyLeaf1);
                dummyRoot1.children.add(s1copy1);
                s1copy1.parent = dummyRoot1;
                dummyLeaf1.parent = dummyRoot1;
                // and to second tree
                Tree dummyRoot2 = new Tree();
                Tree dummyLeaf2 = new Tree("dummy");
                dummyRoot2.children.add(dummyLeaf2);
                dummyLeaf2.parent = dummyRoot2;
                dummyRoot2.children.add(s2copy1);
                s2copy1.parent = dummyRoot2;
                dummyLeaf2.parent = dummyRoot2;
                
                // recursively compute a solution for first subproblem
                Tree subtree_with_dummy = fptMAF(dummyRoot1, dummyRoot2, 0, 0, true);
                if (subtree_with_dummy == null) {
                    //System.out.println("No solution to subproblem");
                    return null;
                }
                boolean dummyRootAlone = subtree_with_dummy.isAlone(dummyLeaf1);
                // switch all edges leading to the dummy leaf on
                if(!dummyRootAlone) {
                    subtree_with_dummy.switchOnTo(dummyLeaf1);
                }
                // remove the dummy leaf
                subtree_with_dummy.deleteLeaf(dummyLeaf1);
                
                Tree subtree_without_dummy = null;
                Tree recursive_subtree = null;
                int with = 0;
                int without = 0;
                boolean withdummyroot;
                if (!dummyRootAlone) {
                    // we need to check if it's possible to get a better solution without the dummy root (i.e. if the dummy root could equally well be alone)
                    
                    // recursively compute a solution for subtree WITHOUT dummy root
                    subtree_without_dummy = fptMAF(s1copy2, s2copy2, 0, 0, true);
                    if (subtree_without_dummy == null) {
                        //System.out.println("No solution to subproblem");
                        return null;
                    }

                    // check if the solution without dummy root is better than the one with dummy root
                    without = subtree_without_dummy.getForest().size();
                    with = subtree_with_dummy.getForest().size();
                    if (without < with) {
                        recursive_subtree = subtree_without_dummy;
                        withdummyroot = false;
                    } else {
                        recursive_subtree = subtree_with_dummy;
                        withdummyroot = true;
                    }
                } else {
                   // we can just take the solution with the dummy root (from which the dummy root has been removed already)
                   recursive_subtree = subtree_with_dummy;
                   withdummyroot = true;
                }
                    
                // now create new subproblem, where subtree is replaced by a dummy leaf
                String dummylabel = "cluster_" + subtree1.getCluster().get(0).toString();
                subtree1.label = dummylabel;
                subtree2.label = dummylabel;
                subtree1.children = new LinkedList();
                subtree2.children = new LinkedList();
                Tree recursive_tree = null;
                if(withdummyroot) {
                    // find a solution where the subtree is replaced by a dummy leaf
                    recursive_tree = fptMAF(T1copy, T2copy, 0, 0, true);
                    if (recursive_tree == null) {
                        //System.out.println("No solution to subproblem");
                        return null;
                    }
                } else {
                    // find a solution where the subtree is removed
                    // actually, the subtree will be replaced by a leaf
                    // but its incoming edge is already switched off
                    T1copy.deleteLeaf(subtree1);
                    subtree2.edgeOn = false;
                    recursive_tree = fptMAF(T1copy, T2copy, 0, 0, true);
                    if (recursive_tree == null) {
                        //System.out.println("No solution to subproblem");
                        return null;
                    }
                }
                
                // replace dummy leaf by recursive_subtree
                recursive_tree.replaceLabelBySubtree(dummylabel, recursive_subtree);
                // clean up
                recursive_tree.clean();
                return recursive_tree;
            }
        }
        
        if(recursive) {
            
            if (!SILENT) {
                System.out.println("Recursively computing agreement forest for trees:");
                System.out.println(Tree1.toString());
                System.out.println(Tree2.toString());
                System.out.println("Running approximation algorithm...");
            }
            
            // RUN APPROXIMATION ALGORITHM
            Tree T2app = approxMAF(Tree1, Tree2);
            List<List<Tree>> app_forest = T2app.getForest();
            app = app_forest.size();
            lb = (int) Math.ceil(app * 1.0 / 4);

            if (!SILENT) {
                System.out.println("Approximation algorithm finished.");
                System.out.println("Approximate agreement forest has " + app + " components.");
                if(lb!=app) {
                    System.out.println("Optimal agreement forest for this subproblem has at least " + lb + " and at most " + app + " components.");
                } else {
                    System.out.println("Optimal agreement forest for this subproblem has " + app + " components.");
                }
            }
        }
        
        if(lb==app) {
            // finished
            return Tree2;
        }

        Tree T;
        int k;
        for (k = lb; k <= Math.min(app, UB); k++) {
            if(!SILENT) {
                System.out.println("Searching for agreement forest of size at most " + k + "...");
            }
            T = kfptMAF(Tree1, Tree2, k, k + 1, recursive);
            if (T != null) {
                return T;
            }
        }

        Tree sol = kfptMAF(Tree1, Tree2, Math.max(lb, k), app, recursive);
        if(!SILENT & recursive) {
            System.out.println("Agreement forest is optimal!");
        }
        return sol;
    }

    private static Tree kfptMAF(Tree Tree1, Tree Tree2, int lb, int ub, boolean recursive) {
        // searches for a solution of size lb <= k < ub
        // set boolean recursive to true if the program is called recursively
        // (if running recursively, intermediate solutions are not being saved to file)

        Tree bestForest = null;

        // make working copies
        Tree T1copy = Tree1.clone();
        Tree T2copy = Tree2.clone();

        Stack stack = new Stack();
        MAF.StackObject so = new MAF.StackObject(T1copy, T2copy, 1);
        stack.push(so);

        while (!stack.isEmpty()) {
            so = (MAF.StackObject) stack.peek();
            stack.pop();
            int components = so.components;

            // clean up
            so.T1 = so.T1.clean();
            so.T2 = so.T2.clean();

            if (DEBUGG) {
                so.T1.toFile("T1now.dot", 0);
                so.T2.toFile("T2now.dot", 0);
            }

            // prune if possible
            if (components >= ub) {
                continue;
            }

            // run approximation-algorithm
            if (APP_PRUNING) {
                // the following code can be used to use the approximation algorithm for pruning
                // however, this seems to slow down the algorithm
                Tree T2approx = approxMAF(so.T1, so.T2);
                List<List<Tree>> approxF = T2approx.getForest();
                int app = approxF.size();
                int applb = (int) Math.ceil(app * 1.0 / 4);

                // prune if possible
                if (applb >= ub) {
                    // no better solution can be found in this subproblem
                    continue;
                }

                if (app == lb) {
                    // optimal solution found
                    bestForest = T2approx;
                    ub = approxF.size();
                    if (!SILENT) {
                        System.out.println("Agreement forest found with " + ub + " components");
                    }
                    break;
                }
            }

            // check if finished
            if (so.T1.getCluster().size() < 2) {
                // T1 consists of at most one leaf, agreement forest found
                // extract forest
                List<List<Tree>> forest = so.T2.getForest();
                if (forest.size() < ub | ub == 0) {
                    // best forest so far
                    bestForest = so.T2;
                    ub = forest.size();
                    if(!SILENT) {
                        System.out.println("Agreement forest found with " + ub + " components");
                    }
                    if (forest.size() == lb) {
                        break;
                    } else if (!recursive) {
                        printSolution(so.T2, Tree1, Tree2);
                    }
                }
                continue;
            }

            // collapse cherries
            boolean success = collapseCherry(so.T1, so.T2);
            if (success) {
                MAF.StackObject newSo = new MAF.StackObject(so.T1, so.T2, components);
                stack.push(newSo);
                continue;
            }

            // delete singletons
            success |= deleteSingletons(so.T1, so.T2);

            if (success) {
                MAF.StackObject newSo = new MAF.StackObject(so.T1, so.T2, components);
                stack.push(newSo);
                continue;
            }

            // Step 0c
            List<List<Tree>> result = cutCherriesFPT(so.T1, so.T1, so.T2);

            if (result != null) {
                for (List<Tree> subproblem : result) {
                    Tree T1sub = subproblem.get(0);
                    Tree T2sub = subproblem.get(1);
                    MAF.StackObject newSo = new MAF.StackObject(T1sub, T2sub, components + 1);
                    stack.push(newSo);
                }
                continue;
            }

            // Step 1 and 2
            List<List<Tree>> result2 = step1and2FPT(so.T1, so.T2);

            if (result2 != null) {
                for (List<Tree> subproblem : result2) {
                    Tree T1sub = subproblem.get(0);
                    Tree T2sub = subproblem.get(1);
                    MAF.StackObject newSo = new MAF.StackObject(T1sub, T2sub, components + 1);
                    stack.push(newSo);
                }
                continue;
            }

            if (!success & DEBUG) {
                System.out.println("No more cuts found...");
            }

            // found an agreement forest
            // extract forest
            List<List<Tree>> forest = so.T2.getForest();
            //boolean correct = checkForest(Tree1,Tree2,forest);
            if ((forest.size() < ub | ub == 0)) {
                // best forest so far
                bestForest = so.T2;
                ub = forest.size();
                if (!SILENT) {
                    System.out.println("Agreement forest found with " + ub + " components.");
                }
                if (forest.size() == lb) {
                    break;
                }
                if (!recursive) {
                    printSolution(so.T2, Tree1, Tree2);
                }
            }
        }

        if (bestForest != null) {
            // decollapse
            bestForest = bestForest.decollapse();
            // and clean up
            bestForest = bestForest.clean();
        }

        return bestForest;
    }

    private static boolean step1and2(Tree T1, Tree T2) {
        // returns true if a cut was made, false otherwise
        // find a cherry in T1
        Tree u = T1.getCherry();
        if (u == null) {
            return false;
        }
        // among its children, find c1,c2 such that c1,c2 are in the same component of T2
        // and s.t. LCA(c1,c2) is at maximum distance from the root of T2
        // among all those pairs, we pick one for which neither of c1 and c2 is a child of their lca in T2, if possible
        List<Tree> l = sameCompMaxDist(u, T2);
        if (!l.isEmpty()) {
            Tree c1 = l.get(0);
            Tree c2 = l.get(1);
            Tree lca = l.get(2);
            // we are working in T2
            c1 = T2.findLeaf(c1);
            c2 = T2.findLeaf(c2);
            // whether we're in case 1 or 2 doesn't matter
            // make the four cuts
            c1.edgeOn = false;
            c2.edgeOn = false;
            c1.parent.edgeOn = false;
            c2.parent.edgeOn = false;
            if (DEBUG) {
                if (c1.parent == lca | c2.parent == lca) {
                    System.out.println("Making a cut of type 2 for leaves: " + c1.label + " and " + c2.label);
                } else {
                    System.out.println("Making a cut of type 1 for leaves: " + c1.label + " and " + c2.label);
                }
            }
            return true;
        }
        return false;
    }

    private static List<List<Tree>> step1and2FPT(Tree T1, Tree T2) {
        // returns at most four subproblems, or null
        // find a cherry in T1
        Tree u = T1.getCherry();
        if (u == null) {
            return null;
        }
        // among its children, find c1,c2 such that c1,c2 are in the same component of T2
        // and s.t. LCA(c1,c2) is at maximum distance from the root of T2
        // among all those pairs, we pick one for which neither of c1 and c2 is a child of their lca in T2, if possible
        List<Tree> l = sameCompMaxDist(u, T2);
        if (!l.isEmpty()) {
            Tree c1 = l.get(0);
            Tree c2 = l.get(1);
            Tree lca = l.get(2);
            c1 = T2.findLeaf(c1);
            c2 = T2.findLeaf(c2);
            if (DEBUG) {
                if (c1.parent == lca | c2.parent == lca) {
                    System.out.println("Making a cut of type 2 for leaves: " + c1.label + " and " + c2.label);
                } else {
                    System.out.println("Making a cut of type 1 for leaves: " + c1.label + " and " + c2.label);
                }
            }

            // whether we're in case 1 or 2 doesn't matter
            // make at most four subproblems

            List<List<Tree>> result = new LinkedList();

            // subproblem 1
            Tree T1S1 = T1.clone();
            Tree T2S1 = T2.clone();
            c1 = T2S1.findLeaf(c1);
            if (c1.edgeOn) {
                c1.edgeOn = false;
                List<Tree> sub1 = new LinkedList();
                sub1.add(T1S1);
                sub1.add(T2S1);
                result.add(sub1);
            }

            // subproblem 2
            Tree T1S2 = T1.clone();
            Tree T2S2 = T2.clone();
            c2 = T2S2.findLeaf(c2);
            if (c2.edgeOn) {
                c2.edgeOn = false;
                List<Tree> sub2 = new LinkedList();
                sub2.add(T1S2);
                sub2.add(T2S2);
                result.add(sub2);
            }

            // subproblem 3
            // we have to be careful here
            // we can't just cut the parent, but we need to cut and refine!
            Tree T1S3 = T1.clone();
            Tree T2S3 = T2.clone();
            c1 = T2S3.findLeaf(c1);
            boolean success = T2S3.refineAndCut(c1);
            if (success) {
                List<Tree> sub3 = new LinkedList();
                sub3.add(T1S3);
                sub3.add(T2S3);
                result.add(sub3);
            }

            // subproblem 4
            Tree T1S4 = T1.clone();
            Tree T2S4 = T2.clone();
            c2 = T2S4.findLeaf(c2);
            success = T2S4.refineAndCut(c2);
            if (success) {
                List<Tree> sub4 = new LinkedList();
                sub4.add(T1S4);
                sub4.add(T2S4);
                result.add(sub4);
            }

            if (!result.isEmpty()) {
                return result;
            }
        }
        return null;
    }

    private static List<Tree> sameCompMaxDist(Tree u, Tree T2) {
        // returns children c1,c2 of u and their lca in T2, s.t.
        // c1,c2 are in the same component of T2
        // lca(c1,c2) in T2 has maximum distance from root
        // neither of c1 and c2 is a child of their lca in T2, if possible
        List<Tree> out = new LinkedList();
        Tree bestc1 = null;
        Tree bestc2 = null;
        Tree bestlca = null;
        int maxdist = -1;
        for (Tree c1 : u.children) {
            for (Tree c2 : u.children) {
                if (c1 == c2) {
                    continue;
                }
                Tree lca = T2.findLCA(c1, c2);
                if (lca == null) {
                    // c1 and c2 are in different components of T2
                    continue;
                }
                int d = T2.getDist(lca);
                if (d > maxdist) {
                    maxdist = d;
                    bestc1 = c1;
                    bestc2 = c2;
                    bestlca = lca;
                }
                if (d == maxdist) {
                    // check if either of c1 and c2 is a parent of their lca in T2
                    if (!lca.hasChild(c1) && !lca.hasChild(c2)) {
                        // we prefer this combo
                        bestc1 = c1;
                        bestc2 = c2;
                        bestlca = lca;
                    }
                }
            }
        }
        if (maxdist > -1) {
            out.add(bestc1);
            out.add(bestc2);
            out.add(bestlca);
        }
        return out;
    }

    private static boolean cutCherries(Tree T1, Tree T2) {
        // Step 0c
        // returns true if a cut was made, false otherwise
        boolean success = false;
        // check if this is a cherry in T1
        boolean cherry = true;
        if (T1.children.isEmpty()) {
            cherry = false;
        }
        for (Tree child : T1.children) {
            if (!child.children.isEmpty()) {
                cherry = false;
            }
        }
        if (cherry) {
            // this is a cherry
            // now we need to check if its leaves are all in different components of T2
            if (diffComponents(T1.children, T2)) {
                // each leaf in the cherry is in a different component of T2
                // cut all leaves in T2
                T2.cutLeaves(T1.children);
                // remove this cherry from T1
                Tree p = T1.parent;
                p.children.remove(T1);
                success = true;
                if (DEBUG) {
                    System.out.print("Making a cut of type 0c for leaves: ");
                    boolean first = true;
                    for (Tree x : T1.children) {
                        if (first) {
                            first = false;
                        } else {
                            System.out.print(", ");
                        }
                        System.out.print(x.label + " ");
                    }
                    System.out.print("\n");
                }

                // suppress
                if (p.parent != null) {
                    p.suppress();
                } else if (p.children.size() == 1) {
                    // suppress the root
                    T1 = p.children.get(0);
                }
            }
        }
        // recurse
        for (int i = 0; i < T1.children.size(); i++) {
            Tree child = T1.children.get(i);
            success = success | cutCherries(child, T2);
        }
        return success;
    }

    private static List<List<Tree>> cutCherriesFPT(Tree v, Tree T1, Tree T2) {
        // Step 0c, FPT version
        // returns a list of tree pairs
        // check if v is a cherry in T1
        boolean cherry = true;
        if (v.children.size() < 2) {
            cherry = false;
        }
        for (Tree child : v.children) {
            if (!child.children.isEmpty()) {
                cherry = false;
            }
        }
        if (cherry) {
            // this is a cherry
            // now we need to check if its leaves are all in different components of T2
            if (diffComponents(v.children, T2)) {
                // each leaf in the cherry is in a different component of T2

                // branch into two subproblems

                // the first subproblem
                Tree T1S1 = T1.clone();
                Tree T2S1 = T2.clone();
                Tree x = v.children.get(0);
                List<Tree> lx = new LinkedList();
                lx.add(x);
                T2S1.cutLeaves(lx);
                T1S1.deleteLeaf(x);
                List<Tree> sub1 = new LinkedList();
                sub1.add(T1S1);
                sub1.add(T2S1);

                // the first subproblem
                Tree T1S2 = T1.clone();
                Tree T2S2 = T2.clone();
                Tree y = v.children.get(1);
                List<Tree> ly = new LinkedList();
                ly.add(y);
                T2S2.cutLeaves(ly);
                T1S2.deleteLeaf(y);
                List<Tree> sub2 = new LinkedList();
                sub2.add(T1S2);
                sub2.add(T2S2);

                if (DEBUG) {
                    System.out.println("Making a cut of type 0c for leaves " + x.label + " and " + y.label);
                }

                // output
                List<List<Tree>> output = new LinkedList();
                output.add(sub1);
                output.add(sub2);
                return output;
            }
        }
        // recurse
        for (int i = 0; i < v.children.size(); i++) {
            Tree child = v.children.get(i);
            List<List<Tree>> recursive = cutCherriesFPT(child, T1, T2);
            if (recursive != null) {
                return recursive;
            }
        }
        return null;
    }

    private static boolean diffComponents(List<Tree> leaves, Tree T2) {
        // returns true if the specifies leaves are all in different components of T2
        int num = 0;
        for (Tree leaf : leaves) {
            if (T2.hasDecendantInComp(leaf)) {
                num++;
            }
        }
        if (num > 1) {
            return false;
        }
        // recurse
        for (Tree child : T2.children) {
            if (!diffComponents(leaves, child)) {
                return false;
            }
        }
        return true;
    }

    private static boolean collapseCherry(Tree T1, Tree T2) {
        boolean success = false;
        if (T1.isCherry()) {
            // try to collapse
            for (Tree c1 : T1.children) {
                Tree c1T2 = T2.findLeaf(c1);
                List<Tree> leaves = new LinkedList();
                for (Tree c2 : T1.children) {
                    if (c1 == c2) {
                        continue;
                    }
                    Tree c2T2 = T2.findLeaf(c2);

                    if (c1T2 == null | c2T2 == null) {
                        System.out.println("Error");
                    }

                    if (T2.hasCherryInComp(c1T2, c2T2)) {
                        leaves.add(c2);
                    }
                }
                if (!leaves.isEmpty()) {
                    // collapse
                    String label = "(" + c1.label;
                    Tree lcaT2 = T2.findLCA(c1, leaves.get(0));
                    for (Tree leaf : leaves) {
                        label += "," + leaf.label;
                        T1.deleteLeaf(leaf);
                        T2.deleteLeaf(leaf);
                    }
                    label += ")";
                    T1.deleteLeaf(c1);
                    T2.deleteLeaf(c1);
                    if (T1.children.isEmpty()) {
                        T1.label = label;
                    } else {
                        Tree newLeaf = new Tree(label);
                        newLeaf.parent = T1;
                        T1.children.add(newLeaf);
                    }
                    if (lcaT2.children.isEmpty()) {
                        lcaT2.label = label;
                    } else {
                        Tree newLeaf = new Tree(label);
                        newLeaf.parent = lcaT2;
                        lcaT2.children.add(newLeaf);
                    }
                    if (DEBUG) {
                        System.out.println("Collapsing cherry: " + label);
                    }
                    return true;
                }
            }
        } else {
            // recurse
            for (Tree child : T1.children) {
                success |= collapseCherry(child, T2);
            }
        }
        return success;
    }

    // old version
//    private static boolean collapseCherry(Tree T1, Tree T2) {
//        // collapse leaves if they have a common parent in both trees
//        // returns true if a cherry was collapsed
//        boolean success = false;
//        List<Tree> cherries = T1.getAllCherries();
//        for (Tree cherry : cherries) {
//            for (Tree c1 : cherry.children) {
//                List<Tree> leaves = new LinkedList();
//                for (Tree c2 : cherry.children) {
//                    if (c1 == c2) {
//                        continue;
//                    }
//                    Tree c1T2 = T2.findLeaf(c1);
//                    Tree c2T2 = T2.findLeaf(c2);
//
//                    if (c1T2 == null | c2T2 == null) {
//                        System.out.println("Error");
//                    }
//
//                    if (T2.hasCherryInComp(c1T2, c2T2)) {
//                        leaves.add(c2);
//                    }
//                }
//                if (!leaves.isEmpty()) {
//                    // collapse
//                    String label = "(" + c1.label;
//                    Tree lcaT2 = T2.findLCA(c1, leaves.get(0));
//                    for (Tree leaf : leaves) {
//                        label += "," + leaf.label;
//                        T1.deleteLeaf(leaf);
//                        T2.deleteLeaf(leaf);
//                    }
//                    label += ")";
//                    T1.deleteLeaf(c1);
//                    T2.deleteLeaf(c1);
//                    if (cherry.children.isEmpty()) {
//                        cherry.label = label;
//                    } else {
//                        Tree newLeaf = new Tree(label);
//                        newLeaf.parent = cherry;
//                        cherry.children.add(newLeaf);
//                    }
//                    if (lcaT2.children.isEmpty()) {
//                        lcaT2.label = label;
//                    } else {
//                        Tree newLeaf = new Tree(label);
//                        newLeaf.parent = lcaT2;
//                        lcaT2.children.add(newLeaf);
//                    }
//                    if (DEBUG) {
//                        System.out.println("Collapsing cherry: " + label);
//                    }
//                    return true;
//                }
//            }
//        }
//        return success;
//    }
    private static boolean deleteSingletons(Tree T1, Tree T2) {
        // returns true if a singleton was deleted, false otherwise
        List<Tree> cluster = T2.getClusterInComp();
        if ((!T2.edgeOn | T2.parent == null) && cluster.size() == 1) {
            // this is a singleton
            // delete the leaf from T1
            Tree leaf = cluster.get(0);
            Tree x = T1.findLeaf(leaf);
            if (x != null) {
                T1.deleteLeaf(x);
                if (DEBUG) {
                    System.out.println("Deleting singleton: " + leaf.label);
                }
                return true;
            }
        }
        // recurse
        for (Tree child : T2.children) {
            if (deleteSingletons(T1, child)) {
                return true;
            }
        }
        return false;
    }
}

class Tree {

    List<Tree> children;
    int number;
    Tree parent;
    String label;
    boolean edgeOn; // indicates if the incoming edge is switched on or off
    int component; // the number of the component that the incoming edge is in (asigned when checking an agreement forest)

    public Tree() {
        children = new LinkedList();
        parent = null;
        label = null;
        edgeOn = true;
        number = 0;
        component = 0;
    }

    public Tree(String l) {
        children = new LinkedList();
        parent = null;
        label = l;
        edgeOn = true;
        number = 0;
        component = 0;
    }

    public void clearComponents() {
        component = 0;
        for (Tree child : children) {
            child.clearComponents();
        }
    }

    public boolean contains(Tree taxon) {
        if (this.equals(taxon)) {
            return true;
        }
        for (Tree child : children) {
            if (child.contains(taxon)) {
                return true;
            }
        }
        return false;
    }

    public static Tree newickToTree(String newick) {
        Tree tree;
        int lastclosepar = newick.lastIndexOf(")");
        int lastcolon = newick.lastIndexOf(":");

        // get rid of semicolon
        if (newick.endsWith(";")) {
            newick = newick.substring(0, newick.length() - 1);
        }

        // get rid of weights
        if (lastcolon > lastclosepar) {
            return newickToTree(newick.substring(0, lastcolon));
        }

        if (newick.startsWith("(")) {
            // internal vertex
            tree = new Tree();
            // split vertex
            int openpar = 0;
            int closepar = 0;
            int start = 1;
            List<String> childrenNewick = new LinkedList();
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
                Tree child = newickToTree(childNewick);
                tree.children.add(child);
                child.parent = tree;
            }
        } else {
            // a leaf
            tree = new Tree(newick);
        }
        return tree;
    }

    public boolean equals(Tree v) {
        if (this == v) {
            return true;
        }
        if (this.label == null | v.label == null) {
            return false;
        }
        if (this.label.equals(v.label)) {
            return true;
        }
        return false;
    }

    public boolean isCherry() {
        if (children.isEmpty()) {
            return false;
        }
        for (Tree child : children) {
            if (!child.children.isEmpty()) {
                return false;
            }
        }
        return true;
    }

    @Override
    public String toString() {
        if (children.isEmpty()) {
            if (label.contains(",")) {
                return "\"" + label + "\"";
            }
            return label;
        }
        String newick = "(";
        for (Tree child : children) {
            newick += child.toString();
            newick += ",";
        }
        // remove the last comma
        newick = newick.substring(0, newick.length() - 1);
        newick += ")";
        if (parent == null) {
            newick += ";";
        }
        return newick;
    }

    public Tree findLeaf(Tree x) {
        Tree leaf;
        if (this.equals(x)) {
            return this;
        }
        for (Tree child : children) {
            leaf = child.findLeaf(x);
            if (leaf != null) {
                return leaf;
            }
        }
        return null;
    }

    public boolean hasDecendant(Tree v) {
        if (this.equals(v)) {
            return true;
        }
        for (Tree child : children) {
            if (child.hasDecendant(v)) {
                return true;
            }
        }
        return false;
    }

    public boolean hasDecendantInComp(Tree v) {
        if (this.equals(v)) {
            return true;
        }
        for (Tree child : children) {
            if (!child.edgeOn) {
                continue;
            }
            if (child.hasDecendantInComp(v)) {
                return true;
            }
        }
        return false;
    }

    public Tree findLCA(Tree v1, Tree v2) {
        for (Tree c1 : children) {
            if (c1.hasDecendant(v1) && c1.hasDecendant(v2)) {
                // both v1 and v2 are below c1
                // search below c1
                return c1.findLCA(v1, v2);
            }
            if (!c1.edgeOn) {
                continue;
            }
            for (Tree c2 : children) {
                if (c1 == c2) {
                    continue;
                }
                if (!c2.edgeOn) {
                    continue;
                }
                if (c1.hasDecendantInComp(v1) && c2.hasDecendantInComp(v2)) {
                    // this is the lca
                    return this;
                }
            }
        }
        return null;
    }

    public void deleteLeaf(Tree x) {
        if (children.isEmpty()) {
            if (this.equals(x)) {
                System.out.println("Error: cannot delete leaf because there are no other vertices");
            }
        }
        for (Tree child : children) {
            if (child.equals(x)) {
                children.remove(child);
                return;
            }
        }
        // recurse
        for (Tree child : children) {
            child.deleteLeaf(x);
        }
    }

    public void suppress() {
        if (children.size() != 1 | parent == null) {
            // do not suppress
            return;
        }
        Tree onlychild = children.get(0);
        parent.children.set(parent.children.indexOf(this), onlychild);
        onlychild.parent = parent;
        // the edge is only on if both edges were on
        onlychild.edgeOn = onlychild.edgeOn & this.edgeOn;
    }

    public List<Tree> getClusterInComp() {
        List cluster = new LinkedList();
        if (children.isEmpty()) {
            // this is a leaf
            cluster.add(this);
        }
        for (Tree child : children) {
            if (child.edgeOn) {
                cluster.addAll(child.getClusterInComp());
            }
        }
        return cluster;
    }

    public List<Tree> getCluster() {
        List cluster = new LinkedList();
        if (children.isEmpty() & label != null) {
            // this is a leaf
            cluster.add(this);
        }
        for (Tree child : children) {
            cluster.addAll(child.getCluster());
        }
        return cluster;
    }

    public void cutLeaves(List<Tree> leaves) {
        for (Tree leaf : leaves) {
            if (this.equals(leaf)) {
                edgeOn = false;
            }
        }
        // recurse
        for (Tree child : children) {
            child.cutLeaves(leaves);
        }
    }

    public List<String> toDot(int s) {
        List<String> out = new LinkedList();
        out.add("strict digraph G {");
        int[] num = new int[1];
        num[0] = 1000;
        out.addAll(this.nodes2dot(num, s));
        out.addAll(this.arcs2dot(s));
        out.add("}");
        this.clearNumbers();
        return out;
    }

    public void clearNumbers() {
        number = 0;
        for (Tree child : children) {
            child.clearNumbers();
        }
    }

    public List<String> nodes2dot(int num[], int s) {
        List<String> out = new LinkedList();
        if (children.isEmpty()) {
            // this is a leaf
            number = num[0];
            out.add(number + " [shape=none, label=\"" + label + "\"];");
        } else {
            number = num[0];
            if (edgeOn && s == 0) {
                out.add(number + " [shape=point, width=0, color=black];");
            } else {
                out.add(number + " [shape=point, width=0, color=grey];");
            }
            for (Tree child : children) {
                num[0]++;
                out.addAll(child.nodes2dot(num, s));
            }
        }
        return out;
    }

    public List<String> arcs2dot(int s) {
        // s is the number of components
        // if there are more than 12 components we do not have enough colours to colour them
        List<String> out = new LinkedList();
        for (Tree child : children) {
            if (child.edgeOn) {
                if (child.component == 0) {
                    if (s == 0) {
                        out.add(number + " -> " + child.number + "[color=black]");
                    } else {
                        out.add(number + " -> " + child.number + "[color=grey]");
                    }
                } else {
                    if (s <= 12) {
                        out.add(number + " -> " + child.number + "[colorscheme=paired12, color=" + child.component + ", label=\"" + child.component + "\"]");
                    } else {
                        out.add(number + " -> " + child.number + "[label=\"" + child.component + "\"]");
                    }
                }
            } else {
                if (s == 0) {
                    out.add(number + " -> " + child.number + "[color=white]");
                } else {
                    out.add(number + " -> " + child.number + "[color=grey]");
                }
            }
        }
        for (Tree child : children) {
            out.addAll(child.arcs2dot(s));
        }
        return out;
    }

    public void cutSomeEdges(double p) {
        if (Math.random() < p) {
            edgeOn = false;
        }
        for (Tree child : children) {
            child.cutSomeEdges(p);
        }
    }

    public void toFile(String filename, int s) {
        List<String> dot = this.toDot(s);
        // write to file
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(filename));
            for (String str : dot) {
                out.write(str + "\n");
            }
            out.close();
        } catch (IOException e) {
            return;
        }
        // transform dot to PDF
        try {
            String line;
            Process p = Runtime.getRuntime().exec("dot -Tpdf " + filename + " -O ");
            BufferedReader bre = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            while ((line = bre.readLine()) != null) {
                System.out.println(line);
            }
            bre.close();
            p.waitFor();
            //System.out.println("Converted to PDF in " + filename + ".pdf");
        } catch (Exception err) {
        }
    }

    public Tree getCherry() {
        boolean cherry = true;
        if (children.size() == 0) {
            // this is a leaf
            cherry = false;
        }
        for (Tree child : children) {
            if (!child.children.isEmpty()) {
                // this child is not a leaf
                cherry = false;
            }
        }
        if (cherry) {
            return this;
        } else {
            // recurse
            for (Tree child : children) {
                Tree u = child.getCherry();
                if (u != null) {
                    return u;
                }
            }
        }
        return null;
    }

    public List<Tree> getAllCherries() {
        List<Tree> out = new LinkedList();
        boolean cherry = true;
        if (children.size() == 0) {
            // this is a leaf
            cherry = false;
        }
        for (Tree child : children) {
            if (!child.children.isEmpty()) {
                // this child is not a leaf
                cherry = false;
            }
        }
        if (cherry) {
            out.add(this);
        }
        // recurse
        for (Tree child : children) {
            out.addAll(child.getAllCherries());
        }
        return out;
    }

    public int getDist(Tree v) {
        // returns the distance from the root to v
        if (this.equals(v)) {
            return 0;
        }
        for (Tree child : children) {
            int d = child.getDist(v);
            if (d != -1) {
                return d + 1;
            }
        }
        return -1;
    }

    public boolean hasChild(Tree c) {
        for (Tree child : children) {
            if (child.equals(c)) {
                return true;
            }
        }
        return false;
    }

    public Tree decollapse() {
        if (label != null && label.contains(",")) {
            return newickToTree(label);
        }
        for (int i = 0; i < children.size(); i++) {
            Tree child = children.get(i);
            boolean eo = child.edgeOn;
            child = child.decollapse();
            children.set(i, child);
            child.parent = this;
            child.edgeOn = eo;
        }
        return this;
    }

    public Tree suppressRoot() {
        if (children.size() == 1) {
            Tree child = children.get(0);
            child.parent = null;
            return child.suppressRoot();
        } else {
            return this;
        }
    }

    public Tree clean() {
        boolean success = true;
        while (success) {
            success = cleanUp();
        }
        return this.suppressRoot();
    }

    public boolean cleanUp() {

        boolean success = false;

        // switch unncessary edges off
        // i.e. if only one edge incident to this vertex is on, we can switch it off, unless this is a leaf
        int num = 0;
        Tree c = null;
        for (Tree child : children) {
            if (child.edgeOn) {
                c = child;
                num++;
            }
        }
        if (num == 1 & (edgeOn == false | parent == null)) {
            c.edgeOn = false;
        }
        if (!children.isEmpty() & num == 0) {
            edgeOn = false;
        }

        // recurse
        for (int i = 0; i < children.size(); i++) {
            Tree child = children.get(i);
            List<Tree> cluster = child.getCluster();
            if (cluster.isEmpty()) {
                // remove ghost subtrees
                children.remove(child);
                success = true;
                i--;
            } else {
                if (child.children.size() == 1) {
                    // suppress indegree-1 outdegree-1 vertices
                    // and recurse
                    Tree grandchild = child.children.get(0);
                    grandchild.edgeOn &= child.edgeOn;
                    success |= grandchild.cleanUp();
                    children.set(i, grandchild);
                    grandchild.parent = this;
                    success = true;
                    i--;
                } else {
                    // recurse
                    success |= child.cleanUp();
                    children.set(i, child);
                }
            }
        }

        return success;
    }

    public List<List<Tree>> getForest() {
        List<List<Tree>> output = new LinkedList();
        if (edgeOn == false | parent == null) {
            // this is the root of a component
            List<Tree> comp = getClusterInComp();
            if (!comp.isEmpty()) {
                output.add(comp);
            }
        }
        // recurse
        for (Tree child : children) {
            output.addAll(child.getForest());
        }
        return output;
    }

    @Override
    public Tree clone() {
        Tree out = new Tree();
        out.number = number;
        out.label = label;
        out.edgeOn = edgeOn;
        out.component = component;
        for (Tree child : children) {
            Tree childclone = child.clone();
            out.children.add(childclone);
            childclone.parent = out;
        }
        return out;
    }

    public boolean colourComponent(List<Tree> comp, int num) {
        List<Tree> cluster = getCluster(); // all the leaves reachable from this vertex
        // the edge entering this vertex is in the component, iff the cluster contains at least one, but not all leaves of the component
        // or if its a leaf that is in the component (leaf edges we will always put in a component)
        boolean containsone = false;
        boolean containsall = true;
        for (Tree x : comp) {
            boolean contains = false;
            for (Tree y : cluster) {
                if (x.equals(y)) {
                    contains = true;
                }
            }
            if (contains) {
                containsone = true;
            } else {
                containsall = false;
            }
        }

        if (containsone && (!containsall || children.isEmpty())) {
            // the edge entering this vertex is in this component
            // check if it's already in another compoennt
            if (component != 0) {
                return false;
            }
            component = num;
        }

        // recurse
        boolean correct = true;
        for (Tree child : children) {
            correct = correct & child.colourComponent(comp, num);
        }

        return correct;
    }

    public boolean hasCherryInComp(Tree v1, Tree v2) {
        // checks if the given vertices form a cherry in some component of the forest
        Tree lca = findLCA(v1, v2);
        if (lca == null) {
            return false;
        }
        for (Tree c1 : lca.children) {
            List<Tree> cluster1 = c1.getClusterInComp();
            if (!c1.edgeOn) {
                continue;
            }
            if (cluster1.size() != 1) {
                continue;
            }
            if (!cluster1.contains(v1)) {
                continue;
            }
            for (Tree c2 : lca.children) {
                if (c2 == c1) {
                    continue;
                }
                if (!c2.edgeOn) {
                    continue;
                }
                List<Tree> cluster2 = c2.getClusterInComp();
                if (cluster2.size() != 1) {
                    continue;
                }
                if (!cluster2.contains(v2)) {
                    continue;
                }

                // cluster1 is a singleton cluster containing only v1
                // cluster2 is a singleton cluster containing only v2
                // that means that v1 and v2 form a cherry
                return true;
            }
        }

        return false;
    }

    public List<List<Tree>> getComponentClusters(int compnum) {
        List<List<Tree>> clusterset = new LinkedList();
        List<Tree> cluster = new LinkedList();
        if (children.isEmpty()) {
            // this is a leaf
            if (component == compnum) {
                cluster.add(this);
            }
            clusterset.add(cluster);
            return clusterset;
        }
        for (Tree child : children) {
            List<List<Tree>> childclusters = child.getComponentClusters(compnum);
            List<Tree> childcluster = childclusters.get(childclusters.size() - 1);
            if (child.component == compnum) {
                cluster.addAll(childcluster);
            }
            if (childcluster.isEmpty()) {
                childclusters.remove(childclusters.size() - 1);
            }
            clusterset.addAll(childclusters);
        }
        clusterset.add(cluster);
        return clusterset;
    }

    public boolean refineAndCut(Tree c) {
        // refines off all children EXCEPT c
        // and cuts the created edge
        // returns false if there is nothing to cut off

        Tree p = c.parent;

        // first we check if there's something to cut
        boolean cut = false;
        for (Tree s : p.children) {
            if (!s.equals(c) & s.edgeOn) {
                cut = true;
            }
        }
        if (!cut) {
            return false;
        }

        // if there's only one child apart from c, we can just switch its edge off
        if (p.children.size() == 2) {
            for (Tree s : p.children) {
                if (!s.equals(c)) {
                    s.edgeOn = false;
                }
            }
            return true;
        }

        // create a new vertex
        Tree n = new Tree();
        for (Tree s : p.children) {
            if (s.equals(c)) {
                continue;
            }
            s.parent = n;
            n.children.add(s);
        }
        p.children = new LinkedList();
        p.children.add(c);
        p.children.add(n);
        n.parent = p;
        n.edgeOn = false;

        return true;
    }

    public List<Tree> findCommonCluster(Tree Tree2) {
        // finds a smallest-possible common cluster
        // returns the corresponding vertices in this tree and in Tree2
        List<Tree> output = new LinkedList();
        List<Tree> cluster = getCluster();
        int n = Tree2.getCluster().size();
        int s = n;
        if (cluster.size() < 2) {
            // we don't want singleton clusters
            return output;
        }
        boolean isCluster = true;
        if (parent == null) {
            // this is the root
            isCluster = false;
        }
        if (cluster.size() == n - 1) {
            // this might be the child of a dummy root
            isCluster = false;
        }
        Tree vertex = Tree2.findCluster(cluster);
        if (vertex == null) {
            // not a common cluster
            isCluster = false;
        }
        if (isCluster) {
            // this is a common cluster
            s = cluster.size();
            output.add(this);
            output.add(vertex);
        }
        for (Tree child : children) {
            List<Tree> recursive = child.findCommonCluster(Tree2);
            if (!recursive.isEmpty()) {
                int cs = recursive.get(0).getCluster().size(); // size of the cluster
                // check if this cluster is bigger than the one we have
                if (cs < s) {
                    // this is the best cluster found so far
                    output = recursive;
                    s = cs;
                }
            }
        }
        return output;
    }

    public Tree findCluster(List<Tree> cluster) {
        List<Tree> myCluster = getCluster();
        if (MAF.clustersEqual(cluster, myCluster)) {
            return this;
        }
        for (Tree child : children) {
            Tree rec = child.findCluster(cluster);
            if (rec != null) {
                return rec;
            }
        }
        return null;
    }

    public void replaceLabelBySubtree(String dummylabel, Tree subtree) {
        // replaces leaves labelled "dummy" by the given subtree
        for(Tree child: children) {
            if(child.label != null && child.label.equals(dummylabel)) {
                children.set(children.indexOf(child), subtree);
                subtree.edgeOn = subtree.edgeOn & child.edgeOn;
                subtree.parent = this;
            } else {
                child.replaceLabelBySubtree(dummylabel, subtree);
            }
        }
    }

    public void switchOnTo(Tree leaf) {
        if (contains(leaf)) {
            edgeOn = true;
        }
        for (Tree child : children) {
            child.switchOnTo(leaf);
        }
    }

    public boolean isAlone(Tree leaf) {
        // returns true if leaf is alone in a component (i.e. if its incoming edge is switched off)
        if (this == leaf & edgeOn) {
            return false;
        }
        if (this == leaf & !edgeOn) {
            return true;
        }
        for (Tree child : children) {
            if (child.isAlone(leaf)) {
                return true;
            }
        }
        return false;
    }
}