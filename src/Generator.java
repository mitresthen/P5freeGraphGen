import java.io.*;
import java.math.BigInteger;
import java.util.*;

import edu.ucla.sspace.graph.*;
import edu.ucla.sspace.graph.isomorphism.*;

/**
 * Class to generate p5-free graphs.
 * Special considerations: Using longs to represent graphs - 
 * each long is a bit representation of an adjacency matrixes lower diagonal, with a 1 as msb.
 * e.g. if the graph is 0001, it's saved as 10001.
 * @author havard
 *
 */
public class Generator {

	/**
	 * Returns a p5 in a graph given the graph and some edge with endpoints u and v
	 * @param g graph to find p5 in
	 * @param u one of the edge endpoints
	 * @param v other edge endpoint
	 * @return A list which contains the 5 nodes in the p5 if there is one, empty list otherwise
	 */
	public static List<Integer> findP5(SparseUndirectedGraph g, int u, int v){
		List<Integer> retList = new ArrayList<Integer>();
		int a = u;
		Set<Integer> Nv = g.getNeighbors(v);
		Set<Integer> B = new HashSet<Integer>(g.vertices());
		B.removeAll(Nv);
		Set<Integer> Na = g.getNeighbors(a);
		Set<Integer> Bp = new HashSet<Integer>(B);
		Bp.removeAll(Na);
		List<HashSet<Integer>> C = listComponents(g, Bp);
		Map<Integer, Integer> compMap = componentMapping(g, Bp);

		Integer c[] = new Integer[C.size()];
		for(int i = 0; i<C.size(); i++){
			c[i] = 0;
		}
		List<Integer> L = new ArrayList<Integer>();
		Set<Integer> NaB = new HashSet<Integer>(Na);
		NaB.retainAll(B);
		for(Integer b : NaB ){
			Set<Integer> NbBp = new HashSet<Integer>(g.getNeighbors(b));
			NbBp.retainAll(Bp);
			for(Integer tmpc : NbBp){
				int j = compMap.get(tmpc);
				c[j]++;
				if(c[j] == 1)
					L.add(j);
			}
			List<Integer> Ldone = new ArrayList<Integer>();
			for(Integer j : L){
				if(c[j] < C.get(j).size()){
					Set<Integer> X = new HashSet<Integer>(g.getNeighbors(b));
					X.retainAll(C.get(j));
					Set<Integer> Y = new HashSet<Integer>(C.get(j));
					Y.removeAll(X);
					for(Integer x : X){
						for(Integer xOpp : g.getNeighbors(x)){
							if(Y.contains(xOpp)){
								retList.add(v);
								retList.add(a);
								retList.add(b);
								retList.add(x);
								retList.add(xOpp);
								return retList;
							}
						}
					}
				}
				c[j] = 0;
				Ldone.add(j);
			}
			L.removeAll(Ldone);
		}
		return retList;
	}

	/**
	 * Method to find and return an induced P5 in a given graph
	 * @param g Graph to test in
	 * @return empty list if P5-free, otherwise returns a P5.
	 */
	public static List<Integer> findP5(SparseUndirectedGraph g){
		List<Integer> retList = new ArrayList<Integer>();
		for(Edge edge : g.edges()){
			retList = findP5(g, edge.from(), edge.to()); 
			if(retList.isEmpty()){
				retList = findP5(g, edge.to(), edge.from());
			}
			if(!retList.isEmpty()){
				return retList;
			}
		}
		return retList;
	}


	/**
	 * Method to map vertexes to components in a subgraph.
	 * @param g the whole graph
	 * @param Bp subgraph to find components within
	 * @return a map of which component each vertex is in.
	 */
	public static Map<Integer, Integer> componentMapping(SparseUndirectedGraph g, Set<Integer> Bp){
		int[] comps = new int[g.vertices().size()];
		int currComp = 0;
		for(int i : g.vertices()){
			comps[i] = -2;
		}
		for(int i : Bp){
			comps[i] = -1;
		}
		Queue<Integer> q = null;
		for(int i : Bp){
			if(comps[i] != -1)
				continue;
			q = new LinkedList<Integer>();
			comps[i] = currComp;
			q.add(i);
			while(!q.isEmpty()){
				int k = q.poll();
				for(int j : g.getNeighbors(k)){
					if(comps[j] == -1){
						comps[j] = currComp;
						q.add(j);
					}
				}
			}
			currComp++;
		}
		Map<Integer, Integer> comp = new HashMap<Integer, Integer>();
		for(int i : Bp){
			comp.put(i, comps[i]);
		}
		return comp;

	}


	/**
	 * Method to get a list of components in a subgraph
	 * @param g main graph
	 * @param Bp subgraph to find components within
	 * @return a list of sets, where each set is a connected component
	 */
	public static List<HashSet<Integer>> listComponents(SparseUndirectedGraph g, Set<Integer> Bp){
		int[] comps = new int[g.vertices().size()];
		int currComp = 0;
		for(int i : g.vertices()){
			comps[i] = -2;
		}
		for(int i : Bp){
			comps[i] = -1;
		}
		Queue<Integer> q = null;
		for(int i : Bp){
			if(comps[i] != -1)
				continue;
			q = new LinkedList<Integer>();
			comps[i] = currComp;
			q.add(i);
			while(!q.isEmpty()){
				int k = q.poll();
				for(int j : g.getNeighbors(k)){
					if(comps[j] == -1){
						comps[j] = currComp;
						q.add(j);
					}
				}
			}
			currComp++;
		}
		List<HashSet<Integer>> C = new ArrayList<HashSet<Integer>>();
		for(int i = 0; i < currComp; i++){
			C.add(new HashSet<Integer>());
		}
		for(int i : Bp){
			C.get(comps[i]).add(i);
		}
		return C;

	}

	/**
	 * Method to test if a graph is P5 free.
	 * @param g graph to test
	 * @return true if P5-free, false otherwise.
	 */
	public static boolean isP5Free(SparseUndirectedGraph g){
		List<Integer> P5 = findP5(g);
		if(P5.isEmpty()){
			return true;
		}
		return false;
	}

	/**
	 * Method to find all graphs 1 node larger than input, that are non-isomorphic and p5-free
	 * @param graph graph to test
	 * @return set of non-isomorphic p5-free graphs that are 1 node larger than input graph.
	 */
	public static IsomorphicSet<SparseUndirectedGraph> p5Extender(SparseUndirectedGraph graph){
		IsomorphicSet<SparseUndirectedGraph> isoSet = new IsomorphicSet<SparseUndirectedGraph>();
		int origGraphSize = graph.vertices().size();
		//Iterate over possible subsets of vertices.
		double DoublePow2N = Math.pow(2.0f, graph.vertices().size());
		int pow2 = (int) DoublePow2N;
		for(int i = 0; i < pow2; i++){
			SparseUndirectedGraph graphClone = new SparseUndirectedGraph(graph);
			graphClone.add(origGraphSize);
			for(int j = 0; j < origGraphSize; j++){
				if((i & ( 1 << j )) >> j == 1){
					graphClone.add(new SimpleEdge(origGraphSize, j));
				}
			}
			if(isP5Free(graphClone)){
				isoSet.add(graphClone);
			}
		}
		return isoSet;
	}


	/**
	 * Method to get a random p5-free graph 1 node larger than input
	 * @param graph graph to extend
	 * @return One p5-free graph, 1 node larger than input
	 */
	public static SparseUndirectedGraph p5RandomExtender(SparseUndirectedGraph graph){
		int origGraphSize = graph.vertices().size();
		double DoublePow2N = Math.pow(2.0f, graph.vertices().size());
		long pow2n = (long) DoublePow2N;
		Set<Long> tried = new HashSet<Long>();
		Random rnd = new Random();
		do{
			long k = rnd.nextLong();
			k = Math.abs(k) % pow2n;
			if(tried.contains(k))
				continue;
			SparseUndirectedGraph graphClone = new SparseUndirectedGraph(graph);
			graphClone.add(origGraphSize);
			for(int j = 0; j < origGraphSize; j++){
				if((k & ( 1 << j )) >> j == 1){
					graphClone.add(new SimpleEdge(origGraphSize, j));
				}
			}
			tried.add(k);
//			SparseUndirectedGraph graphClone2 = new SparseUndirectedGraph(graphClone);
//			if(isP5Free(graphClone) != vP5free(graphClone, origGraphSize)){
//				System.out.println("discrepancy");
//				System.out.println(isP5Free(graphClone));
//				System.out.println(vP5free(graphClone, origGraphSize));
//				System.out.println(dotGraph(graphClone));
//				System.exit(0);
//			}
			if(vP5free(graphClone, origGraphSize)){
				return graphClone;
			}
		}while(tried.size() < pow2n);
		return null;

	}

	/**
	 * Extends a graph up to a p5-free graph of size n.
	 * @param graph graph to extend
	 * @param n size of final graph
	 * @return the extended graph
	 */
	public static SparseUndirectedGraph p5RandomExtender(SparseUndirectedGraph graph, int n){
		if(n<graph.vertices().size()){
			return graph;
		}
		SparseUndirectedGraph extendedGraph = graph.copy(graph.vertices());
		while(extendedGraph.vertices().size() < n){
//			System.out.println("at " + extendedGraph.vertices().size());
			extendedGraph = p5RandomExtender(extendedGraph);
		}
		return extendedGraph;
	}



	/**
	 * Produce a Dot graph representation of the input graph
	 * @param graph graph to print
	 * @return a string containing the Dot graph representation of input
	 */
	public static String dotGraph(SparseUndirectedGraph graph){
		String out = "graph myGraph{ \n";
		for(int i : graph.vertices()){
			out += i + "\n";
		}
		for(Edge e : graph.edges()){
			out += e.from() + " -- " + e.to() + "\n";
		}
		out += "}";
		return out;

	}

	/**
	 * Method to print a graphs Dot-Graph into a file.
	 * @param graph
	 * @param filename
	 */
	public static void writeDotGraph(SparseUndirectedGraph graph, String filename){
		String out = dotGraph(graph);
		try {
			PrintWriter output = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
			output.print(out);
			output.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}


	/**
	 * Method to divide a graph into its connected components.
	 * @param graph graph to divide
	 * @return a List of SparseUndirectedGraphs where each is a connected component in graph.
	 */
	public static List<SparseUndirectedGraph> sparseGraphComponents(SparseUndirectedGraph graph){
		Queue<Integer> q = new LinkedList<Integer>();
		int[] marked = new int[graph.vertices().size()];
		for(int i = 0; i<graph.vertices().size(); i++){
			marked[i] = -1;
		}
		int comp = 0;
		for(int i : graph.vertices()){
			if(marked[i] != -1)
				continue;
			q.add(i);
			marked[i] = comp;
			while(!q.isEmpty()){
				int k = q.poll();
				for(int j : graph.getNeighbors(k)){
					if(marked[j] == -1){
						marked[j] = comp;
						q.add(j);
					}
				}
			}
			comp++;
		}

		List<SparseUndirectedGraph> returnGraphs = new ArrayList<SparseUndirectedGraph>();
		List<HashSet<Integer>> sets = new ArrayList<HashSet<Integer>>();
		for(int i = 0; i<comp; i++){
			sets.add(new HashSet<Integer>());
		}
		for(int i = 0; i<marked.length; i++){
			sets.get(marked[i]).add(i);
		}

		for(HashSet k : sets){
			SparseUndirectedGraph tmpGraph = graph.copy(k);
			returnGraphs.add(tmpGraph);
		}
		return returnGraphs;
	}


	/**
	 * Method that takes a list of graphs and creates a new graph where each pair of graphs in the 
	 * list that differ by exactly 1 edge has an edge between them.
	 * @param graphList list of graphs to connect
	 * @return a SparseUndirectedGraph produced by the format above.
	 */
	public static SparseUndirectedGraph theP5Connection(List<Long> graphList){
		SparseUndirectedGraph connect = new SparseUndirectedGraph(graphList.size());
		for(int i = 0; i<graphList.size(); i++){
			connect.add(i); 
		}
		Map<Long, Integer> mapping = new HashMap<Long, Integer>();
		for(int i = 0; i<graphList.size(); i++){
			mapping.put(graphList.get(i), i);
		}
		for(Long l : graphList){
			for(Long l2 : graphList){
				if(l == l2)
					continue;
				long anded = l ^ l2;
				if( anded != 0 && (anded & (anded-1)) == 0){
					connect.add(new SimpleEdge(mapping.get(l), mapping.get(l2)));
				}
			}
		}
		return connect;
	}


	/**
	 * Function that takes a set of graphs, creates a new graph where each graph in the set is represented as a node
	 * and connects two of the nodes in the new graph if the two graphs they represent differ by 1 edge.
	 * @param graphSet Set of graphs to connect
	 * @return a single graph representing the connection between the graphs.
	 */
	public static SparseUndirectedGraph graphP5Connection(Set<SparseUndirectedGraph> graphSet){
		SparseUndirectedGraph connect = new SparseUndirectedGraph();
		for(int i = 0; i<graphSet.size(); i++){
			connect.add(i);
		}
		List<SparseUndirectedGraph> graphList = new ArrayList<SparseUndirectedGraph>(graphSet);
		IsomorphismTester tester = new VF2IsomorphismTester();
		Map<SparseUndirectedGraph, Integer> graphMap = new HashMap<SparseUndirectedGraph, Integer>();
		int count = 0;
		for(SparseUndirectedGraph graph : graphSet){
			graphMap.put(graph, count++);
		}
		for(int i = 0; i<graphList.size(); i++){
			for(int j= 0; j<graphList.size(); j++){
				if(graphList.get(i).size()+1 != graphList.get(j).size() || i == j)
					continue;
				for(int k : graphList.get(i).vertices()){
					for(int l : graphList.get(j).vertices()){
						if(k != l && !graphList.get(i).getNeighbors(k).contains(l)){
							Edge edge = new SimpleEdge(k,l);
							graphList.get(i).add(edge);
							if(tester.areIsomorphic(graphList.get(i), graphList.get(j))){
								connect.add(new SimpleEdge(i,j));
							}
							graphList.get(i).remove(edge);

						}
					}
				}
			}
		}

		return connect;
	}

	/**
	 * Takes a list of graphs represented as longs and prints them to a file
	 * @param graphs list of graphs to print
	 * @param filename name of file to print to
	 */
	public static void writeLongGraphsToFile(List<Long> graphs, String filename){
		try{
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
			for(Long l : graphs){
				out.println(l);
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	/**
	 * Takes a list of graphs represented as longs and prints them to a file
	 * @param graphs list of graphs to print
	 * @param filename name of file to print to
	 */
	public static void writeBigIntGraphsToFile(List<BigInteger> graphs, String filename){
		try{
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
			for(BigInteger l : graphs){
				out.println(l.toString());
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	private static boolean vuP5free(SparseUndirectedGraph g, Integer v1, Integer v2){
		if(vuP5freeunFlipped(g,v1,v2) && vuP5freeunFlipped(g,v2,v1))
			return true;
		else
			return false;

	}

	/**
	 * Method to determine if graph contains a p5 with (v1,v2) as edge uv in uvwxy in time O(m)
	 * @param g Graph to test in
	 * @return true if P5-free, false otherwise
	 */
	private static boolean vuP5freeunFlipped(SparseUndirectedGraph g, Integer v1, Integer v2){
		Integer v = v1;
		Integer u = v2;
		Set<Integer> Nv = new HashSet<Integer>(g.getNeighbors(v));
		Set<Integer> B = new HashSet<Integer>(g.vertices());
		B.removeAll(Nv);
		Set<Integer> Nu = new HashSet<Integer>(g.getNeighbors(u));
		Set<Integer> Bp = new HashSet<Integer>(B);
		Bp.removeAll(Nu);
		List<HashSet<Integer>> C = listComponents(g, Bp);
		Map<Integer, Integer> compMap = componentMapping(g, Bp);
		Integer c[] = new Integer[C.size()];
		for(int i = 0; i<C.size(); i++){
			c[i] = 0;
		}
		List<Integer> L = new ArrayList<Integer>();

		Set<Integer> NaB = new HashSet<Integer>(Nu);
		NaB.retainAll(B);
		for(Integer b : NaB ){
			Set<Integer> NbBp = new HashSet<Integer>(g.getNeighbors(b));
			NbBp.retainAll(Bp);
			for(Integer tmpc : NbBp){
				int j = compMap.get(tmpc);
				c[j]++;
				if(c[j] == 1)
					L.add(j);
			}
			List<Integer> Ldone = new ArrayList<Integer>();
			for(Integer j : L){
				if(c[j] < C.get(j).size()){
					return false;
				}
				c[j] = 0;
				Ldone.add(j);
			}
			L.removeAll(Ldone);
		}
		return true;
	}

	/**
	 * Executes containsMiddle on both directions of v1 v2
	 * @param g
	 * @param v1
	 * @param v2
	 * @return
	 */
	private static boolean containsMiddleP5(SparseUndirectedGraph g, Integer v1, Integer v2){
		if(containsMiddleP5unFlipped(g, v1, v2) || containsMiddleP5unFlipped(g, v2, v1))
			return true;
		return false;
	}

	/**
	 * Test if graph contains a p5 of the form uvwxy given the edge wx in time O(nm)
	 * @param g
	 * @param v1
	 * @param v2
	 * @return true if contains specified p5, false otherwise
	 */
	private static boolean containsMiddleP5unFlipped(SparseUndirectedGraph g, Integer v1, Integer v2){
		if(!g.contains(v1, v2))
			return false;
		Integer w = v1;
		Integer x = v2;
		Set<Integer> Nx = new HashSet<Integer>(g.getNeighbors(x));
		Set<Integer> NxminNw = new HashSet<Integer>(g.getNeighbors(x));
		Set<Integer> Nw = new HashSet<Integer>(g.getNeighbors(w));
		NxminNw.removeAll(Nw); NxminNw.remove(w);
		for(Integer y : NxminNw){
			Set<Integer> Ny = new HashSet<Integer>(g.getNeighbors(y)); 
			Set<Integer> B = new HashSet<Integer>(g.vertices());
			B.removeAll(Ny);
			B.remove(y);
			Set<Integer> Bp = new HashSet<Integer>();
			Bp.addAll(B);
			Bp.removeAll(Nx);
			Bp.remove(x);
			Bp.remove(w);
			List<HashSet<Integer>> C = listComponents(g, Bp);
			Map<Integer, Integer> compMap = componentMapping(g, Bp);
			Integer c[] = new Integer[C.size()];
			for(int i = 0; i<C.size(); i++){
				c[i] = 0;
			}
			List<Integer> L = new ArrayList<Integer>();

			Set<Integer> NbBp = new HashSet<Integer>(g.getNeighbors(w));
			
			NbBp.retainAll(Bp);
			for(Integer tmpc : NbBp){
				int j = compMap.get(tmpc);
				c[j]++;
				if(c[j] == 1)
					L.add(j);
			}
			for(Integer j : L){
				if(c[j] < C.get(j).size()){
					return true;
				}
			}
		}
		return false;

	}

	/**
	 * Determine if the graph contains a p5 involving edge (v1,v2)
	 * @param g
	 * @param v1
	 * @param v2
	 * @return true if graph contains p5 including v1 and v2, false if not
	 */
	private static boolean containsUVP5(SparseUndirectedGraph g, Integer v1, Integer v2){
		if(vuP5free(g,v1,v2) && !containsMiddleP5(g, v1, v2)){
			return false;
		}
		else{
			return true;
		}
	}

	/**
	 * Determine if the graph contains a p5 involving vertex v
	 * @param g
	 * @param v
	 * @return true if graph contains p5 including v1 and v2, false if not
	 */
	private static boolean vP5free(SparseUndirectedGraph g, Integer v){
		Set<Integer> Nv = new HashSet<Integer>(g.getNeighbors(v));
		for(Integer u : Nv){
			if(!vuP5free(g,v,u) || containsMiddleP5(g, v, u)){
//				System.out.println("vuP5: " + vuP5free(g,v,u));
//				System.out.println("midP5: " + containsMiddleP5(g, v, u));
				return false;
			}
		}
		return true;
	}


	/**
	 * Reads graphs in the long format
	 * @param filename
	 * @return
	 */
	public static List<Long> readLongGraphsFromFile(String filename){
		List<Long> graphs = new ArrayList<Long>();
		try {
			Scanner sc = new Scanner(new File(filename));
			while(sc.hasNextLong()){
				graphs.add(sc.nextLong());
			}
			sc.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return graphs;
	}

	/**
	 * Reads graphs in the long format
	 * @param filename
	 * @return
	 */
	public static List<BigInteger> readBigIntGraphsFromFile(String filename){
		List<BigInteger> graphs = new ArrayList<BigInteger>();
		try {
			Scanner sc = new Scanner(new File(filename));
			while(sc.hasNext()){
				graphs.add(new BigInteger(sc.next()));
			}
			sc.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return graphs;
	}

	/**
	 * Encodes a graph as a long, where the binary representation of the long is 
	 * the lower diagonal of the adjecency matrix.
	 * Encoded in the format 1"binary representation", where 1 is a delimiter.
	 * @param graph graph to convert
	 * @return a long that in binary represents the graph.
	 */
	public static long graphToLong(SparseUndirectedGraph graph){
		int n = graph.vertices().size();
		if(n == 0)
			return 0;
		if(n == 1)
			return 1;
		if(n > 11)
			return -1;
		long longGraph = 2;
		for(int i = 0; i<n; i++){
			Set<Integer> iNeigh = new HashSet<Integer>(graph.getNeighbors(i));
			for(int j = 0; j<i; j++){
				if(iNeigh.contains(j)){
					longGraph = (1 | longGraph); 
				}
				longGraph = longGraph << 1;

			}
		}
		longGraph = longGraph >> 1;
				return longGraph;

	}

	/**
	 * Encodes a graph as a long, where the binary representation of the long is 
	 * the lower diagonal of the adjecency matrix.
	 * Encoded in the format 1"binary representation", where 1 is a delimiter.
	 * @param graph graph to convert
	 * @return a long that in binary represents the graph.
	 */
	public static BigInteger graphToBigInt(SparseUndirectedGraph graph){
		int n = graph.vertices().size();
		if(n == 0)
			return BigInteger.ZERO;
		if(n == 1)
			return BigInteger.ONE;
		BigInteger longGraph = new BigInteger("2");
		for(int i = 0; i<n; i++){
			Set<Integer> iNeigh = new HashSet<Integer>(graph.getNeighbors(i));
			for(int j = 0; j<i; j++){
				if(iNeigh.contains(j)){
					longGraph = BigInteger.ONE.or(longGraph); 
				}
				longGraph = longGraph.shiftLeft(1);// << 1;

			}
		}
		longGraph = longGraph.shiftRight(1);//longGraph >> 1;
		return longGraph;

	}


	/**
	 * Function to convert a long to a graph according to the description given elsewhere in this general class.
	 * @param l long to convert
	 * @return a graph based on the long
	 */
	public static SparseUndirectedGraph longToGraph(Long l){
		SparseUndirectedGraph g = new SparseUndirectedGraph();
		long orig = l;
		int nBits = (63-Long.numberOfLeadingZeros(l));
		int n = (int) ((1+Math.sqrt(1+(4*nBits*2)))/2);
		for(int i = 0; i<n; i++){
			g.add(i);
		}
		for(int i = n-1; i>= 0; i--){
			for(int j = i-1; j>= 0; j--){
				if((1 & l) == 1){
					g.add(new SimpleEdge(i, j));
				}
				l = l >> 1;
			}
		}
		return g;
	}

	/**
	 * Function to convert a long to a graph according to the description given elsewhere in this general class.
	 * @param l long to convert
	 * @return a graph based on the long
	 */
	public static SparseUndirectedGraph bigintToGraph(BigInteger l){
		SparseUndirectedGraph g = new SparseUndirectedGraph();
		BigInteger orig = l;
		//Log base 2
		int nBits = orig.bitLength();// (63-Long.numberOfLeadingZeros(l));
		int n = (int) ((1+Math.sqrt(1+(4*nBits*2)))/2);
		for(int i = 0; i<n; i++){
			g.add(i);
		}
		for(int i = n-1; i>= 0; i--){
			for(int j = i-1; j>= 0; j--){
				if(BigInteger.ONE.and(orig).equals(BigInteger.ONE)){
					g.add(new SimpleEdge(i, j));
				}
				orig = orig.shiftRight(1);//l >> 1;
			}
		}
		return g;
	}


	/**
	 * Function that tests if it is possible to add an edge to the graph and have it remain p5-free
	 * @param graph graph to test
	 * @return true if its possible, false otherwise
	 */
	public static boolean addEdgeP5Free(SparseUndirectedGraph graph){
		int n = graph.vertices().size();
		int edgeCapacity = (n*(n-1))/2;
		if(graph.size() == edgeCapacity)
			return true;
		for(int i : graph.vertices()){
			for(int j : graph.vertices()){
				if(i == j)
					continue;
				if(graph.getEdges(i, j).isEmpty()){
					Edge edge = new SimpleEdge(i,j);
					graph.add(edge);
					if(isP5Free(graph)){
						graph.remove(edge);
						return true;
					}
					graph.remove(edge);
				}
			}
		}
		return false;
	}

	/**
	 * Function that takes a set of graphs and tests if For all graphs in the set there exists an edge that can be added to the graph
	 * so that it remains p5-free
	 * @param graphs set of graphs to test
	 * @return true if there is such an edge for all the graphs, false otherwise
	 */
	public static boolean allAddEdgeP5Free(Set<SparseUndirectedGraph> graphs){
		for(SparseUndirectedGraph gr : graphs){
			if(!addEdgeP5Free(gr)){
				//System.out.println(dotGraph(gr));
				return false;
			}
		}
		return true;


	}


	/**
	 * Function that takes a set of long graphs and outputs a set of graphs
	 * @param graphs set of long graphs to convert
	 * @return set of graphs
	 */
	public static Set<SparseUndirectedGraph> graphToLongSet(Collection<Long> graphs){
		Set<SparseUndirectedGraph> graphSet = new HashSet<SparseUndirectedGraph>();
		for(Long l : graphs){
			graphSet.add(longToGraph(l));
		}
		return graphSet;
	}


	/**
	 * Add a new node with identical neighbourhood as node with id=k
	 * @param inGraph graph which gets a new node
	 * @param k id of node to copy
	 * @return a new graph with an extra node according to method description
	 */
	private static Set<SparseUndirectedGraph> twinNodeExt(SparseUndirectedGraph inGraph, int k){
		Set<SparseUndirectedGraph> graphSet = new HashSet<SparseUndirectedGraph>();
		SparseUndirectedGraph graph = inGraph.copy(inGraph.vertices()); //Create new grpah which is clone of old graph
		int newNode = graph.vertices().size();
		graph.add(newNode);
		for(int i : graph.getNeighbors(k)){
			graph.add(new SimpleEdge(newNode, i));
		}
		if(isP5Free(graph)) //should always be true
			graphSet.add(graph);
		return graphSet;

	}

	public static void genRandomNonIsoGraphs(int size, int numbers, String filename){
		SparseUndirectedGraph graph = new SparseUndirectedGraph();
		graph.add(0);
		int upperBound = size;
		int lim = numbers;
		List<BigInteger> bigGraphs = new ArrayList<BigInteger>();
		IsomorphicSet<SparseUndirectedGraph> bigGraphsIso = new IsomorphicSet<>();
		while(bigGraphsIso.size() < lim)
			bigGraphsIso.add(p5RandomExtender(graph, upperBound));
		for(SparseUndirectedGraph gr : bigGraphsIso ){
			bigGraphs.add(graphToBigInt(gr));
		}
		writeBigIntGraphsToFile(bigGraphs, filename);
	}

	public static void genRandomGraphs(int size, int numbers, String filename){
		SparseUndirectedGraph graph = new SparseUndirectedGraph();
		graph.add(0);
		int upperBound = size;
		int lim = numbers;
		List<BigInteger> bigGraphs = new ArrayList<BigInteger>();
		Set<SparseUndirectedGraph> bigGraphsSet = new HashSet<>();
		while(bigGraphsSet.size() < lim)
			bigGraphsSet.add(p5RandomExtender(graph, upperBound));
		for(SparseUndirectedGraph gr : bigGraphsSet ){
			bigGraphs.add(graphToBigInt(gr));
		}
		writeBigIntGraphsToFile(bigGraphs, filename);
	}


	public static void main(String[] args) {
		//Block for testing connectivity
		//		for(int i = 5; i<=9; i++){
		//			String fileName = "" + i + "graphs.txt";
		//			System.out.println("i: " + i + " size: " + readLongGraphsFromFile(fileName).size());
		//		}
		/*List<Long> graphs = readLongGraphsFromFile("8graphs.txt");
		System.out.println(graphs.size());
		Set<SparseUndirectedGraph> graphSet = graphToLongSet(graphs);
		System.out.println(graphSet.size());
		for(SparseUndirectedGraph gr : graphSet){
			if(!isP5Free(gr))
				System.out.println(dotGraph(gr));
		}
		SparseUndirectedGraph conn = graphP5Connection(graphSet);
		List<SparseUndirectedGraph> graphComps = sparseGraphComponents(conn);*/
		//writeDotGraph(conn, "ConnectedGraph.dot");
		genRandomGraphs(12, 100000, "12graphs.txt");


	}


}
