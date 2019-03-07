package org.knoesis.louvain;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import scala.tools.nsc.classpath.FileUtils;

public class Main {
  // The graph has to have the vertex ids from 0 - num_vertices - 1
  // The graph has to be undirected graph
  public static String id = "";
  public static int outIters = 5;
  public static int innerIters = 100;
  public static int mainIters = 10;
  public static boolean printAreaScore;
  public static double initWeight = 1.0;
  public static double wn = 1.0;
  public static double lambda = 0.1;
  public static String authorFile = "sorted_user_concepts.csv";
  // public static String graphFile = "/Users/bhat236/Documents/untitled
  // folder/CommunityDetection/edges"+id;
  public static String graphFile = "processed_file.txt";
  public static String extra = "temp";

  private static int getNumTopics(String graphFile) throws IOException {
    return Files.readAllLines(Paths.get(graphFile)).get(0).split(" ").length - 2;
  }

  /*
   * 0 0 0 0 ___ 0 0 0 0 1 0 3 1 ___ 1 0 2 1 2 0 5 2 ___ 2 2 3 3 4 3 5 5
   */
  private static int getNumVertices(String graphFile) throws IOException {
    int num_vertices = 0;
    List<String> lines = Files.readAllLines(Paths.get(graphFile));
    for (String line : lines) {
      String splits[] = line.split(" ");
      int start_node = Integer.parseInt(splits[0]);
      int end_node = Integer.parseInt(splits[splits.length - 1]);
      if (start_node > num_vertices) {
        num_vertices = start_node;
      }
      if (end_node > num_vertices) {
        num_vertices = end_node;
      }
    }
    return (num_vertices + 1);
  }

  public void analyzeData(Graph g, String fileName) throws Exception {
    List<String> linesGt = Files.readAllLines(Paths.get(fileName));
    int comm_no = 0;
    HashMap<Integer, ArrayList<String>> ds = new HashMap<Integer, ArrayList<String>>();
    for (String line : linesGt) {
      String splits[] = line.split("\t");
      ArrayList<String> nodesBelonging = new ArrayList<String>();
      for (int i = 1; i < splits.length; i++) {
        nodesBelonging.add((splits[i]));
      }
      ds.put(comm_no++, nodesBelonging);
    }
    g.init(g.graphFiles.get(0)._1);
    for (int c : ds.keySet()) {
      ArrayList<String> vers = ds.get(c);
      System.out.println("putting " + c + " in " + g.vertices[1]);
      for (String v : vers) {
        System.out.println(v);
        g.vertices[Integer.parseInt(v)].community_label = c;
      }
    }

    g.prepareForLouvainIters();
    double d = g.computeGainAll();
    System.out.println(d);

  }

  public void computeClusters() throws Exception {
    // NormalizeInput ni = new NormalizeInput();
    // graphFile = ni.normalizeInput(graphFile);

    int num_topics = getNumTopics(graphFile);
    int num_vertices = getNumVertices(graphFile);
    double[] comm_weights = new double[num_topics];
    for (int i = 0; i < num_topics; i++) {
      comm_weights[i] = Main.initWeight;
    }
    GradientDescent gd = new GradientDescent(num_vertices, num_topics);
    KGOptimizer kgo = new KGOptimizer();
    Graph g = new Graph(num_vertices, num_topics, comm_weights, graphFile, gd, extra);
    // analyzeData(g, "newCircles414");
    for (int i = 0; i < mainIters; i++) {
      System.out.println(i);
      g.startProcessWithPhases(i);
      System.out.println(g.community_info.size());
      kgo.optimizeKnowledgeGraph(authorFile, num_topics, graphFile, g.community_info,
          graphFile, "foundLabels" + extra);
    }
  }

  public static void main(String[] args) throws Exception {
    new Main().computeClusters();    
    // int[] mainItersVals = {1,2,3,5,6,7,8};
    // int[] innerIterVals = {2, 3, 5, 8, 10, 15, 20};
    // int[] outerIterVals = {1, 2, 3, 4};
    // double[] wnVals = {1.1,1.0,1.2};
    // double[] lambdas_v = {0.0000000001, 0.00000000001, 0.0000001, 0.00001,
    // 0.000000000001, 0.00001, 0.01};
    // double[] initWeights = {0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45,
    // 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    // for (double inWe : initWeights) {
    // for (int mainIterVal : mainItersVals) {
    // for (double wnVal : wnVals) {
    // for (double lambda_v : lambdas_v) {
    // mainIters = mainIterVal;
    // Main.initWeight = inWe;
    // innerIters = 200;
    // outIters = 4;
    // wn = wnVal;
    // graphFile = "testEdgesWeights.gr";
    // lambda = lambda_v;
    // extra = "edges_"+id+"_"+inWe+"_"+wn+"_"+lambda+"_"+mainIterVal;
    // new Main().computeClusters();
    // }
    // }
    // }
    // }
  }
}
