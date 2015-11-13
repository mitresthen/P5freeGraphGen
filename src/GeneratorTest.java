import static org.junit.Assert.*;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;

import edu.ucla.sspace.graph.SimpleEdge;
import edu.ucla.sspace.graph.SparseUndirectedGraph;
import edu.ucla.sspace.graph.isomorphism.IsomorphicSet;


public class GeneratorTest {
	private int graphSizes;
	private SparseUndirectedGraph graph;
	private List<SparseUndirectedGraph> p5FreeGraphs;
	private List<SparseUndirectedGraph> p5Graphs;
	
	
	@Before
	public void setUp() {
		graphSizes = 6;
		graph = new SparseUndirectedGraph();
		p5FreeGraphs = new ArrayList<SparseUndirectedGraph>();
		p5Graphs = new ArrayList<SparseUndirectedGraph>();
		for(int i = 0; i<graphSizes; i++){
			graph.add(i);
		}
		for(int i = 0; i<graphSizes-1; i++){
			graph.add(new SimpleEdge(i, i+1));
		}
		p5Graphs.add(graph.copy(graph.vertices()));
		graph.add(new SimpleEdge(2,4));
		graph.add(new SimpleEdge(1,3));
		graph.add(new SimpleEdge(1,4));
		p5FreeGraphs.add(graph.copy(graph.vertices()));
		graph.clearEdges();
		graph.add(new SimpleEdge(0,1));
		graph.add(new SimpleEdge(0,2));
		graph.add(new SimpleEdge(0,3));
		graph.add(new SimpleEdge(0,4));
		graph.add(new SimpleEdge(4,5));
		graph.add(new SimpleEdge(1,2));
		graph.add(new SimpleEdge(2,3));
		graph.add(new SimpleEdge(1,3));
		p5FreeGraphs.add(graph.copy(graph.vertices()));
		graph.clearEdges();
		graph.add(new SimpleEdge(0,1));
		graph.add(new SimpleEdge(0,2));
		graph.add(new SimpleEdge(1,2));
		graph.add(new SimpleEdge(2,4));
		graph.add(new SimpleEdge(3,4));
		graph.add(new SimpleEdge(3,5));
		graph.add(new SimpleEdge(1,5));
		graph.add(new SimpleEdge(1,3));
		p5Graphs.add(graph.copy(graph.vertices()));
		
		
		
	}

	@Test
	public void testFindP5() {
		for(SparseUndirectedGraph gr : p5Graphs){
			assertEquals(Generator.findP5(gr).size(), 5);
		}
		for(SparseUndirectedGraph gr : p5FreeGraphs){
			assertEquals(Generator.findP5(gr).size(), 0);
		}
	}

	@Test
	public void testComponentListAndMapping() {
		graph = new SparseUndirectedGraph();
		graph.add(0); graph.add(1); graph.add(2); graph.add(3);
		graph.add(new SimpleEdge(0,1));
		graph.add(new SimpleEdge(1,2));
		graph.add(new SimpleEdge(2,3));
		Set<Integer> sub = new HashSet<Integer>();
		sub.add(0);
		sub.add(2);
		sub.add(3);
		List<HashSet<Integer>> comps = Generator.listComponents(graph, sub);
		assertEquals(comps.size(), 2);
		Map<Integer, Integer> compMap = Generator.componentMapping(graph, sub);
		assertEquals(compMap.keySet().size(), sub.size());
		
	}

	@Test
	public void testIsP5Free() {
		for(SparseUndirectedGraph gr : p5Graphs){
			assertFalse(Generator.isP5Free(gr));
		}
		for(SparseUndirectedGraph gr : p5FreeGraphs){
			assertTrue(Generator.isP5Free(gr));
		}
	}

	@Test
	public void testP5Extender() {
		for(SparseUndirectedGraph gr : p5FreeGraphs){
			IsomorphicSet<SparseUndirectedGraph> isoSet = Generator.p5Extender(gr);
			for(SparseUndirectedGraph gr2: isoSet){
				assertTrue(Generator.isP5Free(gr2));
			}
		}
	}

	@Test
	public void testP5RandomExtenderSparseUndirectedGraph() {
		for(SparseUndirectedGraph gr : p5FreeGraphs){
			assertTrue(Generator.isP5Free(Generator.p5RandomExtender(gr)));
			assertEquals(Generator.p5RandomExtender(gr).vertices().size(), gr.vertices().size()+1);
		}
		for(SparseUndirectedGraph gr : p5Graphs){
			assertNull(Generator.p5RandomExtender(gr));
		}
	}


	@Test
	public void testTheP5ConnectionAndComponents() {
		long l1 = Long.parseLong("00100", 2);
		long l2 = Long.parseLong("00110", 2);
		List<Long> longs = new ArrayList<Long>();
		longs.add(l1);
		longs.add(l2);
		SparseUndirectedGraph gr = Generator.theP5Connection(longs);
		assertEquals(Generator.sparseGraphComponents(gr).size(), 1);
		long l3 = Long.parseLong("10000",2);
		long l4 = Long.parseLong("00001",2);
		longs.add(l3);
		longs.add(l4);
		gr = Generator.theP5Connection(longs);
		assertEquals(Generator.sparseGraphComponents(gr).size(), 3);
	}

	@Test
	public void testWriteReadLongGraphsToFile() {
		long bound = 10000;
		List<Long> longs = new ArrayList<Long>();
		for(long l = 0l; l<bound; l++){
			longs.add(l);
		}
		Generator.writeLongGraphsToFile(longs, "testFile.txt");
		List<Long> graphs = Generator.readLongGraphsFromFile("testFile.txt");
		assertArrayEquals(graphs.toArray(), longs.toArray());
	}
	@Test
	public void testGraphToFromLong() {
		int graphSize = 10;
		SparseUndirectedGraph gr = new SparseUndirectedGraph();
		for(int i = 0; i<graphSize;i++){
			gr.add(i);
		}
		for(int i = 0; i<graphSize; i++){
			for(int j = i+1; j<graphSize; j += 2){
				gr.add(new SimpleEdge(i,j));
			}
		}
		Long l = Generator.graphToLong(gr);
		SparseUndirectedGraph grCopy = Generator.longToGraph(l);
		assertEquals(grCopy, gr);
	}
	
	@Test
	public void testGraphToFromBigInt() {
		int graphSize = 10;
		SparseUndirectedGraph gr = new SparseUndirectedGraph();
		for(int i = 0; i<graphSize;i++){
			gr.add(i);
		}
		for(int i = 0; i<graphSize; i++){
			for(int j = i+1; j<graphSize; j += 2){
				gr.add(new SimpleEdge(i,j));
			}
		}
		BigInteger l = Generator.graphToBigInt(gr);
		SparseUndirectedGraph grCopy = Generator.bigintToGraph(l);
		assertEquals(grCopy, gr);
	}

	@Test
	public void testWriteReadBigIntGraphsToFile() {
		BigInteger bound = new BigInteger("10000");
		List<BigInteger> longs = new ArrayList<BigInteger>();
		for(BigInteger l = BigInteger.ZERO; l.compareTo(bound) < 0; l = l.add(BigInteger.ONE)){
			longs.add(l);
		}
		Generator.writeBigIntGraphsToFile(longs, "testFile.txt");
		List<BigInteger> graphs = Generator.readBigIntGraphsFromFile("testFile.txt");
		assertArrayEquals(graphs.toArray(), longs.toArray());
	}
}
