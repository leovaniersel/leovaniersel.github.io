
import java.io.*;
import java.util.*;

public class MAF {

    public static int UB = 0;
    public static boolean DEBUG = false;
    public static boolean DEBUGG = false;
    public static boolean SILENT = false;
    public static boolean APP_PRUNING = false;
    public static boolean CLUSTER_REDUCTION = true;
    
    public static void main(String[] args) {
        Tree T1;
        Tree T2;
        
        Long startingTime = System.currentTimeMillis();
        
        System.out.println("----------------------------------------------------------------------");
        System.out.println("MAF: Maximum Agreement Forests for nonbinary trees");
        System.out.println("----------------------------------------------------------------------");
        System.out.println("Implements an approximation and an FPT algorithm.");
        System.out.println("The approximation algorithm constructs an agreement");
        System.out.println("forest that is at most a factor 4 from optimal.");
        System.out.println("The FPT algorithm constructs a maximum agreement forest (MAF).");
        System.out.println("By Leo van Iersel (2012).");
        System.out.println("http://homepages.cwi.nl/~iersel/MAF/");

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
            } else if(i>0) {
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
        T1.toFile("T1.dot", 0);
        System.out.println("Second input tree:");
        System.out.println(T2.toString());
        T2.toFile("T2.dot", 0);

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
        if(!allDifferent(X1)) {
            System.out.println("Error: first tree has duplicate leaves.");
            return;
        }
        if(!allDifferent(X2)) {
            System.out.println("Error: second tree has duplicate leaves.");
            return;
        }
        if (clustersEqual(X1,X2)) {
            System.out.println("Trees have " + n + " leaves.");
            T1res = T1;
            T2res = T2;
        } else {
            System.out.println("Trees do not have identical leaf sets");
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
            if(n == X2res.size()) {
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
        System.out.println("Computation took " + approxTime/1000 + " seconds.");
        // write agreement forest to file
        System.out.println("Graphical representation of agreement forest saved as app_forest.dot");
        T2app.toFile("app_forest.dot", 0);
        // write T1 with component numbers to file
        System.out.println("T1 edge-coloured by agreement forest saved as app_T1forest.dot");
        T1.toFile("app_T1forest.dot", app);
        // write T2 with component numbers to file
        System.out.println("T2 edge-coloured by agreement forest saved as app_T2forest.dot");
        T2.toFile("app_T2forest.dot", app);

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
            if(lb != app) {
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
            System.out.println("Graphical representation of agreement forest saved as forest.dot");
            T2fpt.toFile("forest.dot", 0);
            // write T1 with component numbers to file
            System.out.println("T1 edge-coloured by agreement forest saved as T1forest.dot");
            T1.toFile("T1forest.dot", opt);
            // write T2 with component numbers to file
            System.out.println("T2 edge-coloured by agreement forest saved as T2forest.dot");
            T2.toFile("T2forest.dot", opt);

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
            System.out.println("Total computation time: " + fptTime/1000 + " seconds.");
            System.out.println("Maximum agreement forest has " + opt + " components.");
        }
    }
    
    public static boolean clustersEqual(List<Tree> X1, List<Tree> X2) {
        List<Tree> common = new LinkedList();
        for (Tree x1 : X1) {
            if (clusterContainsTaxon(X2,x1)) {
                common.add(x1);
            }
        }
        return (X1.size() == X2.size() & X1.size() == common.size());
    }
    
    public static boolean allDifferent(List<Tree> X) {
        for(Tree x : X) {
            for(Tree y : X) {
                if(X.indexOf(x) == X.indexOf(y)) {
                    continue;
                }
                if(x.equals(y)) {
                    return false;
                }
            }
        }
        return true;
    }
    
    public static boolean clusterContainsTaxon(List<Tree> cluster, Tree taxon) {
        for(Tree x : cluster) {
            if(x.equals(taxon)) {
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
            sol.toFile("forest.dot", 0);
            Tree1.toFile("T1forest.dot", val);
            Tree2.toFile("T2forest.dot", val);
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

        if (recursive) {
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
                // and to second tree
                Tree dummyRoot2 = new Tree();
                Tree dummyLeaf2 = new Tree("dummy");
                dummyRoot2.children.add(dummyLeaf2);
                dummyLeaf2.parent = dummyRoot2;
                dummyRoot2.children.add(s2copy1);
                s2copy1.parent = dummyRoot2;
                
                // recursively compute a solution for first subproblem
                Tree subtree_with_dummy = fptMAF(dummyRoot1, dummyRoot2, 0, 0, true);
                if (subtree_with_dummy == null) {
                    System.out.println("Error: no solution to subproblem");
                    return null;
                }
                boolean dummyRootAlone = subtree_with_dummy.isAlone(dummyLeaf1);
                // switch all edges leading to the dummy leaf on
                subtree_with_dummy.switchOnTo(dummyLeaf1);
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
                        System.out.println("Error: no solution to subproblem");
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
                        System.out.println("Error: no solution to subproblem");
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
                        System.out.println("Error: no solution to subproblem");
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
        StackObject so = new StackObject(T1copy, T2copy, 1);
        stack.push(so);

        while (!stack.isEmpty()) {
            so = (StackObject) stack.peek();
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
                    if(!SILENT) {
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
                    } else if(!recursive) {
                        printSolution(so.T2, Tree1, Tree2);
                    }
                }
                continue;
            }

            // collapse cherries
            boolean success = collapseCherry(so.T1, so.T2);
            if (success) {
                StackObject newSo = new StackObject(so.T1, so.T2, components);
                stack.push(newSo);
                continue;
            }

            // delete singletons
            success |= deleteSingletons(so.T1, so.T2);

            if (success) {
                StackObject newSo = new StackObject(so.T1, so.T2, components);
                stack.push(newSo);
                continue;
            }

            // Step 0c
            List<List<Tree>> result = cutCherriesFPT(so.T1, so.T1, so.T2);

            if (result != null) {
                for (List<Tree> subproblem : result) {
                    Tree T1sub = subproblem.get(0);
                    Tree T2sub = subproblem.get(1);
                    StackObject newSo = new StackObject(T1sub, T2sub, components + 1);
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
                    StackObject newSo = new StackObject(T1sub, T2sub, components + 1);
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
                if(!SILENT) {
                    System.out.println("Agreement forest found with " + ub + " components.");
                }
                if (forest.size() == lb) {
                    break;
                }
                if(!recursive) {
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
        if(this.equals(taxon)) {
            return true;
        }
        for(Tree child : children) {
            if(child.contains(taxon)) {
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
                if(s==0) {
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
        if(cluster.size() < 2) {
            // we don't want singleton clusters
            return output;
        }
        boolean isCluster = true;
        if(parent == null) {
            // this is the root
            isCluster = false;
        }
        if(cluster.size() == n - 1) {
            // this might be the child of a dummy root
            isCluster = false;
        }
        Tree vertex = Tree2.findCluster(cluster);
        if(vertex == null) {
            // not a common cluster
            isCluster = false;
        }
        if(isCluster) {
            // this is a common cluster
            s = cluster.size();
            output.add(this);
            output.add(vertex);
        }
        for(Tree child: children) {
            List<Tree> recursive = child.findCommonCluster(Tree2);
            if(!recursive.isEmpty()) {
                int cs = recursive.get(0).getCluster().size(); // size of the cluster
                // check if this cluster is bigger than the one we have
                if( cs < s ) {
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
        if(MAF.clustersEqual(cluster,myCluster)) {
            return this;
        }
        for(Tree child : children) {
            Tree rec = child.findCluster(cluster);
            if(rec != null) {
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
                subtree.parent = this;
            } else {
                child.replaceLabelBySubtree(dummylabel, subtree);
            }
        }
    }
    
    public void switchOnTo(Tree leaf) {
        if(contains(leaf)) {
            edgeOn = true;
        }
        for(Tree child : children) {
            child.switchOnTo(leaf);
        }
    }
    
    public boolean isAlone(Tree leaf) {
        // returns true if leaf is alone in a component (i.e. if its incoming edge is switched off)
        if(this==leaf & edgeOn) {
            return false;
        }
        if(this==leaf & !edgeOn) {
            return true;
        }
        for(Tree child : children) {
            if(child.isAlone(leaf)) {
                return true;
            }
        }
        return false;
    }
}
