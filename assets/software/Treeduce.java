package treeduce;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author Leo van Iersel (2015)
 */
public class Treeduce {

    static boolean SILENT = false;
    static boolean CHAINS = true;
    static boolean DEBUG = false;
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        System.out.println("----------------------------------------------------------------------");
        System.out.println("Treeduce: Kernelization Algorithm for Hybridization Number");
	System.out.println("on Multiple Nonbinary Trees");
        System.out.println("By Leo van Iersel, Steven Kelk and Celine Scornavacca (2015)");
        System.out.println("-----------------------------USAGE------------------------------------");
        System.out.println("java Treeduce input.tree [k] [-s]");
        System.out.println("input.tree\t input file containing trees with equal taxon sets");
        System.out.println("\t\t in Newick format");
        System.out.println("k\t\t parameter (default value is 1)");
        System.out.println("-s\t\t silent mode");
        System.out.println("----------------------------------------------------------------------");
        
        Long startingTime = System.currentTimeMillis();

        List<Tree> trees = new LinkedList();
        List<Tree> taxa = new LinkedList();
        int t = 0; // number of trees
        int n = 0; // number of taxa
        int N = 0; // original number of taxa
        
        int subtreeremoved = 0;
        int degreechainremoved = 0;
        int starchainremoved = 0;
        
        String treeFile = args[0];
        int k = 1;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-s")) {
                SILENT = true;
            }
        }
        
        if (args.length > 1) {
            try {
                k = Integer.parseInt(args[1]);
            } catch (NumberFormatException e) {
                System.out.println("Using default value for k");
            }
        } else {
            System.out.println("Using default value for k");
        }
        System.out.println("k = " + k);
        
        // parse the trees
        boolean success = true;
        File file = new File(treeFile);
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(file));
            while (success) {
                String newick = "";
                try {
                    newick = reader.readLine();
                } catch (IOException e) {
                }
                if (newick == null) {
                    success = false;
                    reader.close();
                } else {
                    t++;
                    Tree tree = Tree.newickToTree(newick);
                    tree.clean();
                    trees.add(tree);
                    if (!SILENT) {
                        System.out.println("Input tree " + t + ": ");
                    }
                    if (!SILENT) {
                        System.out.println(tree.toString());
                    }
                    if (DEBUG) {
                        tree.toFile("tree_" + t + ".dot", 0);
                    }

                    if (newick.length() == 0) {
                        System.out.println("Error: input tree is empty");
                        return;
                    }
                }
            }
            reader.close();

        } catch (FileNotFoundException e) {
            System.out.println("Error: input tree is empty");
            return;
        } catch (IOException e) {
        }

        if (t < 2) {
            System.out.println("Error: at least two input trees required");
            return;
        }

        // check if all trees have the same leaves
        taxa = trees.get(0).getCluster();
        n = taxa.size();
        if (false) {
            if (!allDifferent(taxa)) {
                System.out.println("Error: first tree has duplicate leaves.");
                return;
            }
            for (int i = 1; i < trees.size(); i++) {
                Tree T2 = trees.get(i);
                List<Tree> X2 = T2.getCluster();
                if (!allDifferent(X2)) {
                    System.out.println("Error: tree " + (i + 1) + " has duplicate leaves.");
                    return;
                }
                if (!clustersEqual(taxa, X2)) {
                    System.out.println("Error: trees do not have identical leaf sets");
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
        }
        
        // find the maximum outdegree
        int maxout = 0;
        for(Tree tree : trees) {
            int d = tree.getOutDeg();
            if(d > maxout) {
                maxout = d;
            }
        }
        System.out.println("Maximum outdegree: " + maxout);
        System.out.println("Number of trees: " + t);
        System.out.println("Number of leaves: " + n);
        N = n;

        // -----> subtree reduction
        System.out.println("Starting subtree reduction...");
        Long subtreeStartingTime = System.currentTimeMillis();
        fastSubtreeReduction(trees, taxa);
        Long subtreeTime = System.currentTimeMillis() - subtreeStartingTime;
        taxa = trees.get(0).getCluster();
        subtreeremoved += n-taxa.size();
        n = taxa.size();
        if(!SILENT) System.out.println("There are " + n + " leaves remaining after subtree reduction.");

        Long chainStartingTime = System.currentTimeMillis();
        if (CHAINS) {

            // degree-based chain reduction
            System.out.println("Starting degree-based chain reduction...");
            int max_length = 5 * k * (maxout - 1);
            if (k == 0) {
                max_length = 2;
            }
            if(!SILENT) System.out.println("Reducing all chains to length " + max_length + " ...");
            boolean reduced = true;
            while (reduced) {
                if(max_length > taxa.size()) {
                    break;
                }
                reduced = chainReduction(trees, taxa, max_length, -1);

                if (reduced) {

                    clean(taxa);
                    
                    if(!SILENT) System.out.println("Repeating subtree reduction...");
                    
                    fastSubtreeReduction(trees, taxa);
                }
            }
            degreechainremoved += n - taxa.size();
            n = taxa.size();
            if(!SILENT) System.out.println("There are " + n + " leaves remaining after degree-based chain reduction.");

            // star-based chain reduction
            System.out.println("Starting star-based chain reduction...");
            for (int q = t - 1; q >= 0; q--) {
                max_length = (int) Math.pow(5 * k, t - q);
                if (k == 0) {
                    max_length = 2;
                }
                if(!SILENT) System.out.println("Reducing all " + q + "-star chains to length " + max_length + " ...");
                if(max_length > taxa.size()) {
                    break;
                }
                reduced = chainReduction(trees, taxa, max_length, q);
                if (reduced) {
                    if(!SILENT) System.out.println("Repeating subtree reduction...");
                    clean(taxa);
                    fastSubtreeReduction(trees, taxa);
                    q = t;
                }
            }
            starchainremoved += n - taxa.size();
            n = taxa.size();
            if(!SILENT) System.out.println("There are " + n + " leaves remaining after star-based chain reduction.");

        }
        
        Long chainTime = System.currentTimeMillis() - chainStartingTime;
        Long elapsedTime = System.currentTimeMillis() - startingTime;

        for(Tree tree : trees) {
            tree = tree.suppressRoot();
            if (!SILENT) {
                System.out.println("Reduced tree " + (trees.indexOf(tree)+1) + ": ");
            }
            if (!SILENT) {
                System.out.println(tree.toString());
            }
            if (DEBUG) {
                tree.toFile("tree_" + (trees.indexOf(tree)+1) + "_reduced.dot", 0);
            }
        }
        
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter("reduced.tree"));
            for (Tree tree : trees) {
                out.write(tree.toString() + "\n");
            }
            out.close();
        } catch (IOException e) {
            return;
        }
        System.out.println("Reduced trees in Newick format saved to reduced.tree");

        System.out.println("");
        System.out.println("---------------------------SUMMARY------------------------------------");
        System.out.println("Parameter: k = " + k);
        System.out.println("Number of leaves: n = " + N);
        System.out.println("Number of trees: t = " + t);
        System.out.println("Maximum outdegree: " + maxout);
        System.out.println("Leaves removed by subtree reduction: " + subtreeremoved);
        System.out.println("Leaves removed by degree-based chain reduction: " + degreechainremoved);
        System.out.println("Leaves removed by star-based chain reduction: " + starchainremoved);
        System.out.println("Leaves remaining: " + n);
        System.out.println("Total computation time: " + elapsedTime/1000 + " seconds.");
        System.out.println("Computation time subtree reduction: " + subtreeTime/1000 + " seconds.");
        System.out.println("Computation time chain reduction: " + chainTime/1000 + " seconds.");
        System.out.println("----------------------------------------------------------------------");
    }
    
    public static void clean(List<Tree> taxa) {
        for(int i = 0; i < taxa.size(); i++) {
            Tree x = taxa.get(i);
            if(x.label == null) {
                taxa.remove(i);
                i--;
            }
        }
    }
    
    public static void fastSubtreeReduction(List<Tree> trees, List<Tree> taxa) {
        List<List<Tree>> allRoots = new LinkedList();
        List<List<List<Tree>>> childrenInSubtree = new LinkedList();
        // first create empty lists
        for (Tree tree : trees) {
            allRoots.add(new LinkedList());
            childrenInSubtree.add(new LinkedList());
        }
        // now add the taxa as singleton common pendant subtrees
        // and their parents as the roots
        for (Tree taxon : taxa) {
            for (int t = 0; t < trees.size(); t++) {
                Tree x = trees.get(t).findLeaf(taxon);
                allRoots.get(t).add(x.parent);
                if(x.parent == null) {
                    System.out.println("Null parent");
                }
                List<Tree> children = new LinkedList();
                children.add(x);
                childrenInSubtree.get(t).add(children);
            }
        }
        // keep merging subtrees
        boolean merged = true;
        while(merged) {
            merged = false;
            int num_subtrees = allRoots.get(0).size();
            for(int i = 0; i < num_subtrees-1; i++) {
                for(int j = i+1; j < num_subtrees ; j++) {
                    boolean sameroot = true;
                    for(int t = 0; t < trees.size(); t++) {
                        Tree u = allRoots.get(t).get(i);
                        Tree v = allRoots.get(t).get(j);
                        if(u!=v) {
                            sameroot = false;
                            break;
                        }
                    }
                    if (sameroot) {
                        // merge
                        for (int t = 0; t < trees.size(); t++) {
                            Tree root = allRoots.get(t).get(i);
                            List<Tree> c_i = childrenInSubtree.get(t).get(i);
                            List<Tree> c_j = childrenInSubtree.get(t).get(j);
                            if(c_i.size() + c_j.size() == root.children.size()) {
                                // move root up
                                allRoots.get(t).set(i, root.parent);
                                c_i.clear();
                                c_i.add(root);
                            } else {
                                c_i.addAll(c_j);
                            }
                            allRoots.get(t).remove(j);
                            childrenInSubtree.get(t).set(i, c_i);
                            childrenInSubtree.get(t).remove(j);
                        }
                        merged = true;
                        break;
                    }
                }
                if(merged) break;
            }
        }
        // replace each common subtree by a single leaf
        int num_subtrees = allRoots.get(0).size();
        int stnum = 1;
        for (int i = 0; i < num_subtrees; i++) {
            List<Tree> children = childrenInSubtree.get(0).get(i);
            if(children.size()==1 & children.get(0).isLeaf()) {
                // this is a trivial common pendant subtree
                continue;
            }
            String l = "ST" + stnum;
            stnum++;
            for (int t = 0; t < trees.size(); t++) {
                Tree root = allRoots.get(t).get(i);
                List<Tree> subtreeChildren = childrenInSubtree.get(t).get(i);
                for(Tree child : subtreeChildren) {
                    root.children.remove(child);
                }
                Tree subtree = new Tree(l);
                subtree.parent = root;
                root.children.add(subtree);
            }
        }
        taxa = trees.get(0).getCluster();
        for(Tree tree : trees) {
            tree.clean();
        }
    }
    
    public static void subtreeReduction(List<Tree> trees, List<Tree> taxa) {
        boolean found = true;
        Tree x = null;
        Tree y = null;
        while (found) {
            found = false;
            for (Tree taxon1 : taxa) {
                for (Tree taxon2 : taxa) {
                    if (taxon2 == taxon1) {
                        continue;
                    }
                    boolean cherry = true;
                    for (Tree tree : trees) {
                        if (!tree.gotCherry(taxon1, taxon2)) {
                            cherry = false;
                            break;
                        }
                    }
                    if (!cherry) {
                        continue;
                    } else {
                        found = true;
                        x = taxon1;
                        y = taxon2;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
            if (found) {
                // collapse
                if (!SILENT) {
                    System.out.println("Collapsing cherry on " + x.label + " and " + y.label + " ...");
                }
                for (Tree tree : trees) {
                    tree.collapseCherry(x, y);
                    tree.clean();
                }
                taxa.remove(y);
            }
        }
    }
    
    public static boolean chainReduction(List<Tree> trees, List<Tree> taxa, int max_length, int star) {
        boolean reduced = false;
        List<Tree> chain = new LinkedList();
        for (Tree x_s : taxa) {
            for (Tree x_t : taxa) {
                if (x_s == x_t) {
                    continue;
                }
                // find a common q-star s-t chain of maximum size
                List<Tree> stchain = getMaxChain(trees, x_s, x_t, star);
                if (stchain != null && stchain.size() > chain.size()) {
                    chain = stchain;
                }
            }
        }

        if (chain.size() > max_length) {
            reduced = true;
            // order the leaves in the chain
            chain = toChain(trees, chain);
            if(!SILENT) System.out.println("Chain found: " + chain.toString());
            if(!SILENT) System.out.print("Deleting leaves: ");
            for (int i = max_length; i < chain.size(); i++) {
                Tree leaf = chain.get(i);
                if(!SILENT & i > max_length) System.out.print(", ");
                if(!SILENT) System.out.print(leaf);
                for (Tree tree : trees) {
                    tree.deleteLeaf(leaf);
                }
                taxa.remove(leaf);
            }
            if(!SILENT) System.out.print("\n");
        }
        for(Tree tree : trees) {
            tree.clean();
        }
        return reduced;
    }

    public static List<Tree> getMaxChain(List<Tree> trees, Tree s, Tree t, int q) {
        List<Tree> chain = new LinkedList();
        List<Tree> mustBe = new LinkedList(); // the leaves that have to be in the chain
        mustBe.add(s);
        mustBe.add(t);
        List<Tree> canBe = new LinkedList(); // the leaves that can be in the chain
        int stars = 0;
        for (Tree tree : trees) {
            List<Tree> mustBe_t = new LinkedList();
            List<Tree> canBe_t = new LinkedList();
            Tree p_s = tree.getParentOf(s);
            Tree p_t = tree.getParentOf(t);
            if(p_s == p_t) {
                // here we could check if the chain is pendant rather than star
                stars++;
            }
            List<Tree> internal = p_s.getInternalChainVertices(s,t);
            if (internal == null) {
                return null;
            }
            for(Tree v : internal) {
                for(Tree c : v.children) {
                    if(c.isLeaf()) {
                        mustBe_t.add(c);
                    }
                }
            }
            List<Tree> all = new LinkedList();
            all.addAll(internal);
            all.add(p_s);
            all.add(p_t);
            for(Tree v : all) {
                for(Tree c : v.children) {
                    if(c.isLeaf()) {
                        canBe_t.add(c);
                    }
                }
            }
            if(trees.indexOf(tree) == 0) {
                mustBe.addAll(mustBe_t);
                canBe.addAll(canBe_t);
            } else {
                for(Tree x : mustBe_t) {
                    if(!mustBe.contains(x)) {
                        mustBe.add(x);
                    }
                }
                canBe.retainAll(canBe_t);
            }
        }
        for(Tree x : mustBe) {
            if(!canBe.contains(x)) {
                return null;
            }
        }
        if(q != -1 && stars != q) {
            return null;
        }
        
        // check if mustBe is chainable
        if(toChain(trees,mustBe) == null) {
            return null;
        }
        
        // extend mustBe with maximum number of leaves from canBe
        chain = maximize(trees,canBe,mustBe,s,t);
        
        return chain;
    }
    
    public static List<Tree> maximize(List<Tree> trees, List<Tree> canBe, List<Tree> mustBe, Tree s, Tree t) {
        List<Tree> chain = new LinkedList();

        // the vertices of the digraph
        List<Tree> vertices = new LinkedList();
        // it's easier to also have s and t in the graph
        // we clone all nodes to make sure that we don't mess up
        Tree v_s = new Tree(s.label);
        Tree v_t = new Tree(t.label);
        vertices.add(v_s);
        vertices.add(v_t);
        for(Tree x : canBe) {
            if(mustBe.contains(x)) {
                continue;
            }
            // check if this leaf can be added to the chain
            mustBe.add(x);
            if(toChain(trees, mustBe) != null) {
                return null;
            }
            mustBe.remove(x);
            Tree v = new Tree(x.label);
            vertices.add(v);
        }
        
        // the edges of the digraph
        for(Tree u : vertices) {
            for(Tree v : vertices) {
                if(u == v) {
                    continue;
                }
                // check if there should be an edge from u to v
                boolean edge = true;
                for(Tree tree : trees) {
                    Tree p_u = tree.getParentOf(u);
                    Tree p_v = tree.getParentOf(v);
                    if(!p_u.find(p_v)) {
                        edge = false;
                    }
                }
                if(edge) {
                    u.children.add(v);
                    v.parents.add(u);
                }
            }
        }
                
        // find distances from s by longest paths
        v_s.longestPath();
        List<Tree> path = v_t.backTrack();
        
        chain.addAll(mustBe);
        for(Tree v : path) {
            if(!v.equals(s) && !v.equals(t)) {
                chain.add(new Tree(v.label));
            }
        }
        
        return chain;
    }
    
    public static List<Tree> toChain(List<Tree> trees, List<Tree> leaves) {
        
        List<Tree> chain = new LinkedList();
        
        // the vertices of the digraph
        List<Tree> vertices = new LinkedList();
        for(Tree leaf : leaves) {
            Tree vertex = new Tree(leaf.label);
            vertices.add(vertex);
        }
        
        // the edges of the digraph
        for(Tree u : vertices) {
            for(Tree v : vertices) {
                if(u == v) {
                    continue;
                }
                // check if there should be an edge from u to v
                for(Tree tree : trees) {
                    Tree p_u = tree.getParentOf(u);
                    Tree p_v = tree.getParentOf(v);
                    if(p_u == p_v) {
                        continue;
                    }
                    if(p_u.find(p_v)) {
                        u.children.add(v);
                        v.parents.add(u);
                        break;
                    }
                }
            }
        }
        
        // check if the digraph is acyclic
        boolean success = true;
        while(success) {
            success = false;
            for(Tree v : vertices) {
                if(v.parents.isEmpty()) {
                    chain.add(new Tree(v.label));
                    // delete v
                    vertices.remove(v);
                    for(Tree c : v.children) {
                        c.parents.remove(v);
                    }
                    success = true;
                    break;
                }
            }
        }
        // no vertex with outdegree 0
        if (vertices.isEmpty()) {
            // graph is acyclic
            return chain;
        } else {
            // graph is NOT acyclic
            return null;
        }
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

    public static boolean clustersEqual(List<Tree> X1, List<Tree> X2) {
        List<Tree> common = new LinkedList();
        for (Tree x1 : X1) {
            if (clusterContainsTaxon(X2, x1)) {
                common.add(x1);
            }
        }
        return (X1.size() == X2.size() & X1.size() == common.size());
    }

    public static boolean clusterContainsTaxon(List<Tree> cluster, Tree taxon) {
        for (Tree x : cluster) {
            if (x.equals(taxon)) {
                return true;
            }
        }
        return false;

    }

}

class Tree {

    List<Tree> children;
    List<Tree> parents; // used for abusing this class for digraphs
    int number;
    int number2;
    Tree parent;
    String label;
    List<String> labels;

    public Tree() {
        children = new LinkedList();
        parents = new LinkedList();
        parent = null;
        label = null;
        number = 0;
        number2 = 0;
        labels = new LinkedList();
    }

    public Tree(String l) {
        children = new LinkedList();
        parents = new LinkedList();
        parent = null;
        label = l;
        number = 0;
        number2 = 0;
        labels = new LinkedList();
        labels.add(label);
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
    
    public Tree addDummyRoot() {
        Tree root = new Tree();
        root.children.add(this);
        return root;
    }
    
    public int getOutDeg() {
        int d = children.size();
        for(Tree child : children) {
            int e = child.getOutDeg();
            if(e > d) {
                d = e;
            }
        }
        return d;
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

    
    public void longestPath() {
        for(Tree child : children) {
            child.number2++;
            if(child.number < number + 1) {
                child.number = number + 1;
            }
        }
        for(Tree child : children) {
            // only go down once all incoming arcs have been traversed
            if(child.number2 == child.parents.size()) {
                child.longestPath();
            }
        }
    }

    public List<Tree> backTrack() {
        List<Tree> path = new LinkedList();
        if(number == 0) {
            // this must be s
            path.add(this);
            return path;
        }
        for(Tree p : parents) {
            if(p.number == number - 1) {
                path = p.backTrack();
                path.add(this);
                return path;
            }
        }
        return null;
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

    @Override
    public boolean equals(Object o) {
        if ( !(o instanceof Tree) ) return false;
        Tree v = (Tree) o;
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
            if(label == null) {
                return "";
            }
            if (label.contains(",")) {
                return "\"" + label + "\"";
            } else {
                return label;
            }
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
    
    public void clean() {
        removeGhosts();
        cleanRec();
    }

    public void cleanRec() {
        // recurse first
        for(Tree child : children) {
            child.cleanRec();
        }
        // now suppress any indegree-1 outdegree-1 children
        for (Tree child : children) {
            if (child.children.size() == 1) {
                // suppress indegree-1 outdegree-1 vertices
                Tree grandchild = child.children.get(0);
                children.set(children.indexOf(child), grandchild);
                grandchild.parent = this;
            }
        }
    }

    public void removeGhosts() {
        for (int i = 0; i < children.size(); i++) {
            Tree child = children.get(i);
            List<Tree> cluster = child.getCluster();
            if (cluster.isEmpty()) {
                // remove ghost subtrees
                children.remove(child);
                i--;
            }
        }
        for(Tree child : children) {
            child.removeGhosts();
        }
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

    public Tree suppressRoot() {
        if (children.size() == 1) {
            Tree child = children.get(0);
            child.parent = null;
            return child.suppressRoot();
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
            out.add(number + " [shape=point, width=0, color=black];");
            for (Tree child : children) {
                num[0]++;
                out.addAll(child.nodes2dot(num, s));
            }
        }
        return out;
    }

    public boolean find(Tree v) {
        if(this.equals(v)) {
            return true;
        }
        for(Tree child : children) {
            if(child.find(v)) {
                return true;
            }
        }
        return false;
    }
    
    public List<String> arcs2dot(int s) {
        // s is the number of components
        // if there are more than 12 components we do not have enough colours to colour them
        List<String> out = new LinkedList();
        for (Tree child : children) {
            out.add(number + " -> " + child.number + "[color=black]");
        }
        for (Tree child : children) {
            out.addAll(child.arcs2dot(s));
        }
        return out;
    }

    public boolean collapseCherry(Tree x, Tree y) {
        boolean success = false;
        for (int i = 0; i < children.size(); i++) {
            Tree child = children.get(i);
            if (child.equals(x)) {
                child.labels.add(y.label);
            }
            if (child.equals(y)) {
                children.remove(i);
                i--;
                success = true;
            }
        }
        for (Tree child : children) {
            if (!success) {
                success = child.collapseCherry(x, y);
            }
        }
        return success;
    }
    
    public List<Tree> getInternalChainVertices(Tree s, Tree t) {
        List<Tree> internal = new LinkedList();
        boolean parentofs = false;
        for (Tree child : children) {
            if (child.equals(t)) {
                return internal;
            }
            if(child.equals(s)) {
                parentofs = true;
            }
        }
        for(Tree child : children) {
            List<Tree> recursive = child.getInternalChainVertices(s, t);
            if(!child.isLeaf() && recursive == null && !parentofs) {
                return null;
            }
            if(recursive != null) {
                internal.addAll(recursive);
                if(!parentofs) {
                    internal.add(this);
                }
                return internal;
            }
        }
        return null;
    }

    public Tree getParentOf(Tree v) {
        for (Tree child : children) {
            if (child.equals(v)) {
                return this;
            }
            Tree parent = child.getParentOf(v);
            if (parent != null) {
                return parent;
            }
        }
        return null;
    }
    
    public boolean isLeaf() {
        return children.isEmpty();
    }

    public boolean gotCherry(Tree x, Tree y) {
        boolean xchild = false;
        boolean ychild = false;
        for (Tree child : children) {
            if (child.label == null) {
                continue;
            }
            if (child.equals(x)) {
                xchild = true;
            }
            if (child.equals(y)) {
                ychild = true;
            }
        }
        if (xchild && ychild) {
            return true;
        }
        if (xchild) {
            return false;
        }
        if (ychild) {
            return false;
        }
        for (Tree child : children) {
            if (child.gotCherry(x, y)) {
                return true;
            }
        }
        return false;
    }
}
