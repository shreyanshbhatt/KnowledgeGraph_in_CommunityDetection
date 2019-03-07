package org.knoesis.louvain;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import org.knoesis.results.CommunityFMeasureCesna;
import org.knoesis.results.CommunityJaccard;
import scala.Tuple2;
import scala.Tuple3;

public class Graph {
  Vertex vertices[] = null; // has to be ordered by id
  Community[] communities; // each index in here corresponds to its community
                           // label
  int num_vertices_backup;
  EdgeWeights total_edgeweights;
  int num_topics;
  int num_topics_backup;
  double[] initial_commweights;
  double[] initial_commweights_backup;
  GradientDescent gd;
  ArrayList<Tuple2<String, String>> graphFiles;
  String extra;
  ArrayList<Tuple2<HashMap<Integer, Integer>, HashMap<Integer, Integer>>> levelInfo;
  Community[] communities_backup;
  HashMap<Integer, Community> vertex_to_community_index = new HashMap<Integer, Community>();
  HashMap<Integer, ArrayList<Vertex>> community_to_vertices_index =
      new HashMap<Integer, ArrayList<Vertex>>();
  HashMap<Integer, Tuple2<ArrayList<Vertex>, ArrayList<Vertex>>> community_info =
      new HashMap<Integer, Tuple2<ArrayList<Vertex>, ArrayList<Vertex>>>();

  Graph(int numVertices, int num_topics, double[] comm_weights, String graphFile,
      GradientDescent gd, String extra) throws IOException {
    vertices = new Vertex[numVertices];
    communities = new Community[numVertices];
    this.num_topics = num_topics;
    total_edgeweights = new EdgeWeights();
    graphFiles = new ArrayList<Tuple2<String, String>>();
    graphFiles.add(new Tuple2<String, String>(graphFile, ""));
    total_edgeweights.feature_edge_weight_derivative_helper = new ArrayList<Double>();
    for (int i = 0; i < num_topics; i++) {
      total_edgeweights.feature_edge_weight_derivative_helper.add(1.0);
    }
    this.initial_commweights = comm_weights;
    this.gd = gd;
    this.extra = extra;
    communities_backup = new Community[numVertices];
    // removed init from here
    // init(graphFiles.get(0)._1);
    levelInfo = new ArrayList<Tuple2<HashMap<Integer, Integer>, HashMap<Integer, Integer>>>();
    vertex_to_community_index = new HashMap<Integer, Community>();
    community_to_vertices_index = new HashMap<Integer, ArrayList<Vertex>>();
    community_info = new HashMap<Integer, Tuple2<ArrayList<Vertex>, ArrayList<Vertex>>>();
  }

  private void addEdge(Vertex src, Vertex dest, ArrayList<Double> similarity) {
    Neighbor n = new Neighbor();
    n.similarity = similarity;
    n.weighted_similarity_der = new ArrayList<Double>();
    n.neighbour_label = dest.id;
    src.neighbors.add(n);
  }

  private void addCommunity(int vertex_id) {
    if (communities[vertex_id] != null) {
      return;
    }
    Community newComm = new Community();
    ArrayList<Double> topic_weights = new ArrayList<Double>();
    for (double d : initial_commweights) {
      topic_weights.add(d);
    }
    newComm.topic_weights = topic_weights;
    communities[vertex_id] = newComm;
  }

  private void addVertexInfo(int comm_label, int vertex_id) {
    if (vertices[vertex_id] != null) {
      return;
    }
    Vertex newVert = vertices[vertex_id];
    newVert = new Vertex();
    newVert.id = vertex_id;
    newVert.community_label = comm_label;
    newVert.prev_iter_community_label = comm_label;
    newVert.neighbors = new ArrayList<Neighbor>();
    newVert.feature_vertex_derivate_helper = new ArrayList<Double>();
    newVert.should_skip = false;
    for (int i = 0; i < num_topics; i++) {
      newVert.feature_vertex_derivate_helper.add(1.0);
    }
    vertices[vertex_id] = newVert;
  }

  private void normalizeInputData(String inputFile) throws IOException {
    List<String> lines = Files.readAllLines(Paths.get(inputFile), Charset.defaultCharset());
    ArrayList<Tuple3<Integer, ArrayList<Double>, Integer>> newValues =
        new ArrayList<Tuple3<Integer, ArrayList<Double>, Integer>>();
    BufferedWriter bw = new BufferedWriter(new FileWriter(inputFile));
    ArrayList<Double> curSumArr = new ArrayList<Double>();
    for (String line : lines) {
      String splits[] = line.split(" ");
      int src = Integer.parseInt(splits[0]);
      int target = Integer.parseInt(splits[splits.length - 1]);
      ArrayList<Double> scores = new ArrayList<Double>();
      for (int i = 1; i < splits.length - 1; i++) {
        scores.add(Double.parseDouble(splits[i]));
        double curSum = 0.0;
        if (curSumArr.size() >= i) {
          curSum = curSumArr.get(i - 1);
          curSum += Math.pow(Double.parseDouble(splits[i]), 2);
          curSumArr.set(i - 1, curSum);
          continue;
        }
        curSum += Math.pow(Double.parseDouble(splits[i]), 2);
        curSumArr.add(curSum);
      }
      newValues.add(new Tuple3<Integer, ArrayList<Double>, Integer>(src, scores, target));
    }

    for (int i = 0; i < curSumArr.size(); i++) {
      double updated = Math.sqrt(curSumArr.get(i));
      curSumArr.set(i, updated);
    }

    for (Tuple3<Integer, ArrayList<Double>, Integer> tuple : newValues) {
      bw.write(tuple._1() + " ");
      for (int i = 0; i < tuple._2().size(); i++) {
        bw.write(((double) (tuple._2().get(i) / curSumArr.get(i))) + " ");
      }
      bw.write("" + tuple._3());
      bw.newLine();
    }
    bw.close();
  }

  public void init(String inputFile) throws IOException {
    // normalizeInputData(inputFile);
    List<String> lines = Files.readAllLines(Paths.get(inputFile), Charset.defaultCharset());
    for (String line : lines) {
      String[] splits = line.split(" ");
      int source_vertex_id = Integer.parseInt(splits[0]);
      int dest_vertex_id = Integer.parseInt(splits[splits.length - 1]);
      ArrayList<Double> similarity = new ArrayList<Double>();

      for (int i = 1; i < num_topics + 1; i++) {
        if (splits[i].equals("NaN")) {
          similarity.add(0.0);
          continue;
        }
        similarity.add(Double.parseDouble(splits[i]));
      }

      ArrayList<Double> similaritytedge = new ArrayList<Double>();

      for (int i = 1; i < num_topics + 1; i++) {
        if (splits[i].equals("NaN")) {
          similaritytedge.add(0.0);
          continue;
        }
        similaritytedge.add(Double.parseDouble(splits[i]));
      }

      addVertexInfo(source_vertex_id, source_vertex_id);
      if (source_vertex_id != dest_vertex_id)
        addVertexInfo(dest_vertex_id, dest_vertex_id);
      addCommunity(source_vertex_id);
      if (source_vertex_id != dest_vertex_id)
        addCommunity(dest_vertex_id);
      addEdge(vertices[source_vertex_id], vertices[dest_vertex_id], similarity);
      if (source_vertex_id != dest_vertex_id)
        addEdge(vertices[dest_vertex_id], vertices[source_vertex_id], similaritytedge);
    }
  }
  // We use the following for parallel processing
  // private void mergeVertices() {
  // for (Vertex v : vertices) {
  // if (v.neighbors.size() == 1) {
  // Vertex neighboring_vertex = vertices[v.neighbors.get(0).neighbour_label];
  // addEdge(neighboring_vertex, neighboring_vertex,
  // v.neighbors.get(0).similarity);
  // v.should_skip = true;
  // }
  // }
  // }

  private double dotProduct(ArrayList<Double> similarities, ArrayList<Double> weights) {
    double res = 0.0;
    for (int i = 0; i < similarities.size(); i++) {
      res += similarities.get(i) * weights.get(i);
    }
    return res;
  }

  private double kernel_total_edge_weight(ArrayList<Double> similarities,
      ArrayList<Double> weights) {
    double dot_product = dotProduct(similarities, weights);
    // return (double)(Math.exp(dot_product));
    return (double) ((double) (Main.wn - ((double) 1.0 / (double) Math.exp(dot_product))) / 2.0);
  }

  // private double kernel(double similarity, double weight) {
  // return (double)((double)1.0 -
  // ((double)1/(double)Math.exp(similarity*weight))); // wn has to be larger
  // than 1
  // }

  // compute gradient like this. gain based on dot product
  // (kernal_total_edge_weight). similarity1 * exp(-(alpha1, alpha2)*(sim1,
  // sim2))
  // optimize gain computed based on kernal_total_edge_weight. we don't need
  // featured vertex degree at all.

  private double kernel_total_edge_weight_der(ArrayList<Double> similarities,
      ArrayList<Double> weights, int feat_num) {
    double dot_product = dotProduct(similarities, weights);
    // return (double)((double)similarities.get(feat_num) *
    // ((double)Math.exp(dot_product)));
    // System.out.println(similarities.get(feat_num));
    // System.out.println((double)(1.0/(double)Math.exp(dot_product)));

    return (double) ((double) ((double) similarities.get(feat_num)
        * ((double) 1.0 / (double) Math.exp(dot_product)) - Main.lambda) / 2.0);
  }

  // private double kernel_der(double similarity, double weight) {
  // return (double)((double) (similarity * weight) *
  // ((double)1/(double)Math.exp(similarity*weight)));
  // }

  private void updateVertexRelatedInfo(Vertex src, Neighbor n, boolean ifSelfComm) {
    ArrayList<Double> src_comm_weights = communities[src.prev_iter_community_label].topic_weights;
    ArrayList<Double> dest_comm_weights =
        communities[vertices[n.neighbour_label].prev_iter_community_label].topic_weights;
    if (ifSelfComm) {
      // this is an edge to the vertex that belong to the community we are
      // interested in
      // hence we will update the derivation information as well
      for (int i = 0; i < src_comm_weights.size(); i++) {
        double new_edge_weight_der_computed =
            kernel_total_edge_weight_der(n.similarity, src_comm_weights, i)
                + kernel_total_edge_weight_der(n.similarity, dest_comm_weights, i);

        src.feature_vertex_derivate_helper.set(i, (double) src.feature_vertex_derivate_helper.get(i)
            - n.weighted_similarity_der.get(i) + new_edge_weight_der_computed);

        total_edgeweights.feature_edge_weight_derivative_helper.set(i,
            total_edgeweights.feature_edge_weight_derivative_helper.get(i)
                - n.weighted_similarity_der.get(i) + new_edge_weight_der_computed);
        n.weighted_similarity_der.set(i, new_edge_weight_der_computed);
      }
    }
    double new_edge_weight = kernel_total_edge_weight(n.similarity, src_comm_weights)
        + kernel_total_edge_weight(n.similarity, dest_comm_weights);
    communities[src.community_label].total_community_degree -= n.edge_weight;
    communities[src.community_label].total_community_degree += new_edge_weight;
    src.total_vertex_degree -= n.edge_weight;
    src.total_vertex_degree += new_edge_weight;
    total_edgeweights.total_edge_weight -= n.edge_weight;
    total_edgeweights.total_edge_weight += new_edge_weight;
    n.edge_weight = new_edge_weight;
  }

  private void handleEdge(Vertex src, Neighbor n) {
    updateVertexRelatedInfo(src, n, true);
    if (vertices[n.neighbour_label].community_label != src.community_label) {
      // we have a neighbor that belongs to a different community.
      // so we need to update that community information as well
      for (Neighbor neighOfNeigh : vertices[n.neighbour_label].neighbors) {
        if (neighOfNeigh.neighbour_label != src.id) {
          continue;
        }
        updateVertexRelatedInfo(vertices[n.neighbour_label], neighOfNeigh, false);
      }
    }
  }

  public void updateWeightsDegreeGradient(int comm) {
    for (Vertex v : community_to_vertices_index.get(comm)) {
      for (Neighbor n : v.neighbors) {
        handleEdge(v, n);
      }
    }
  }

  public void prepareForLouvainIters() {
    // fill m=total_edge_weight, ml=feature_edge_weight
    // ki=total_vertex_degree, kl=feature_vertex_degree,
    // aC=total_community_degree

    // process every vertex. for every neighboring vertex, increment the current
    // degree according to neighboring vertex's weight.
    // also increment neighboring vertex's weight calculated according to each
    // feature

    // now look at this vertex's community label and in that community update
    // the community degree
    // also update the m and ml according to this vertex's degree

    // look at its community
    // and update the community degree

    // for each vertex
    total_edgeweights.total_edge_weight = 0.0;
    for (int i = 0; i < num_topics; i++) {
      total_edgeweights.feature_edge_weight_derivative_helper.set(i, 0.0);
    }
    for (Community c : communities) {
      c.total_community_degree = 0.0;
      c.total_vertices = 0;
    }
    for (Vertex v : vertices) {
      // for each neighbor
      v.total_vertex_degree = 0;
      for (int i = 0; i < num_topics; i++) {
        v.feature_vertex_derivate_helper.set(i, 0.0);
      }
      for (Neighbor n : v.neighbors) {
        ArrayList<Double> community_weight_vector =
            communities[(vertices[(n.neighbour_label)].prev_iter_community_label)].topic_weights;
        ArrayList<Double> own_community_weight_vector =
            communities[(v.prev_iter_community_label)].topic_weights;
        ArrayList<Double> similarity_vector = n.similarity;
        n.edge_weight = kernel_total_edge_weight(similarity_vector, community_weight_vector)
            + kernel_total_edge_weight(similarity_vector, own_community_weight_vector);
        total_edgeweights.total_edge_weight += n.edge_weight;
        v.total_vertex_degree += n.edge_weight;
        for (int i = 0; i < community_weight_vector.size(); i++) {
          double new_edge_weight_der_computed =
              kernel_total_edge_weight_der(similarity_vector, community_weight_vector, i)
                  + kernel_total_edge_weight_der(similarity_vector, own_community_weight_vector, i);

          if ((Double.isNaN(new_edge_weight_der_computed))) {
            System.out.println(dotProduct(similarity_vector, community_weight_vector));
            System.out.println(dotProduct(similarity_vector, own_community_weight_vector));
            System.out.println("issueeeee");
            System.exit(1);
          }
          n.weighted_similarity_der.add(new_edge_weight_der_computed);

          v.feature_vertex_derivate_helper.set(i,
              (double) v.feature_vertex_derivate_helper.get(i) + new_edge_weight_der_computed);

          total_edgeweights.feature_edge_weight_derivative_helper.set(i,
              total_edgeweights.feature_edge_weight_derivative_helper.get(i)
                  + new_edge_weight_der_computed);
        }
      }
      communities[(v.community_label)].total_community_degree += v.total_vertex_degree;
      communities[(v.community_label)].total_vertices++;
    }
  }

  // keep it just in case
  // private void debugInter() {
  // System.out.println("\n\n\n");
  // System.out.println("------------------------");
  // for (Vertex v : vertices) {
  // System.out.println("Now printing for "+v.id);
  // System.out.println("\tdegree = "+v.total_vertex_degree);
  // for (Neighbor n : v.neighbors) {
  // System.out.println("\tedge "+v.id+" and "+n.neighbour_label+" =
  // "+n.edge_weight);
  // }
  // }
  // System.out.println("Now printing for communities");
  // for (int i = 0; i < communities.length; i++) {
  // System.out.println("\t for community "+i+" it is
  // "+communities[i].total_community_degree);
  // }
  // System.out.println("-----------------------");
  // System.out.println("\n\n\n");
  // }

  private boolean updateCommunityLabel() {
    // for each vertex, list its neighboring communities.
    // for each community compute the gain and assign new_label for the maximum
    // gain
    HashMap<Integer, Double> communities_seen = null;
    // Since we are doing a sequential processing instead of a parallel
    // processing.
    // we will randomize the vertex order and start our processing
    int[] vertexOrder =
        new Random().ints(0, vertices.length).distinct().limit(vertices.length).toArray();
    boolean should_stop = true;
    for (int vertexId : vertexOrder) {
      Vertex v = vertices[vertexId];
      // FIXME: we are not using the following since we are operating algorithm
      // as it is sequential instead of parallel
      if (v.should_skip) {
        continue;
      }
      communities_seen = new HashMap<Integer, Double>();
      // communities_seen.put(v.old_community_label, 0.0);
      communities_seen.put(v.community_label, 0.0);
      for (Neighbor n : v.neighbors) {
        int possible_community_label = vertices[n.neighbour_label].community_label;
        double existing_comm_weight = 0.0;
        if (communities_seen.containsKey(possible_community_label))
          existing_comm_weight = communities_seen.get(possible_community_label);
        communities_seen.put(possible_community_label, existing_comm_weight + n.edge_weight);
      }

      double max_gain = 0.0;
      int max_comm_label = v.community_label;
      for (int comm_label : communities_seen.keySet()) {
        double gain = (double) ((double) (((double) communities_seen.get(comm_label)
            - (double) communities_seen.get(v.community_label))
            / (double) total_edgeweights.total_edge_weight)
            + (double) ((double) ((double) v.total_vertex_degree
                * (((double) communities[v.community_label].total_community_degree
                    - (double) communities[(comm_label)].total_community_degree)))
                / (double) ((double) 2 * ((double) total_edgeweights.total_edge_weight
                    * (double) total_edgeweights.total_edge_weight))));

        // TODO: not braking ties randomly
        if (Double.compare(gain, max_gain) > 0) {
          max_gain = gain;
          max_comm_label = comm_label;
        }
      }

      // Parallel heuristic
      // FIXME: we are operating in sequential mode. so skipping the following
      // if (communities[v.old_community_label].total_vertices == 1 &&
      // communities[max_comm_label].total_vertices == 1 && max_comm_label >
      // v.old_community_label) {
      // continue;
      // }

      communities[max_comm_label].total_vertices++;
      communities[v.community_label].total_vertices--;
      communities[max_comm_label].total_community_degree += v.total_vertex_degree;
      communities[v.community_label].total_community_degree -= v.total_vertex_degree;
      if (v.community_label != max_comm_label) {
        should_stop = false;
      }
      v.community_label = max_comm_label;
    }
    return should_stop;
  }

  // FIXME: we are not doing parallel processing so we have the updated
  // community labels with us.
  // hence we don't need following for now.
  // public void prepareForNextIteration() {
  // for (Vertex v : vertices) {
  // v.old_community_label = v.new_community_label;
  // }
  // }

  public void printCommunityAffiliations() {
    for (Vertex v : vertices) {
      System.out.println(v.community_label);
    }
  }

  // keep it just incase. replace this with sophisticated logging mechanism but
  // I'm in too hurry to do all that by myself for now. (need to graduate soon
  // and nobody cares about coding)
  // private void printIntterWeights() {
  // for (Vertex v : vertices) {
  // for (Neighbor n : v.neighbors) {
  // System.out.println("weighted sim between "+v.id +" and
  // "+n.neighbour_label+" = "+n.edge_weight);
  // }
  // }
  // }

  public boolean louvainIter(int iterNum) {
    // Step 1 : By this time you should have the vertices filled vertex and
    // neighbors
    // where the similarity is given from the previous step and
    // weighted_similarity
    // is calculated based on a given weight vector in community.

    // debugInter();
    // Step 2 : Now we have weights, m. ki, and aC. we will update community
    // labels
    return updateCommunityLabel();
    // FIXME: do the following if you're doing the parallel processing
    // prepareForNextIteration();
  }

  private void initCommLabels() {
    for (Vertex v : vertices) {
      v.community_label = v.id;
    }
  }

  private void mainIterWithPhases(int level_id) throws IOException {
    deactivateComm();
    initCommLabels();
    prepareForLouvainIters();
    for (int i = 0; i < Main.innerIters; i++) {
      boolean should_stop = louvainIter(i);
      if (should_stop) {
        break;
      }
      
    }
    double gain = computeModularity();
    System.out.println("gain = "+gain);
    prepareForLouvainIters();
    LevelUp lu = new LevelUp();
    HashMap<Integer, Integer> vertices_to_comm = new HashMap<Integer, Integer>();
    HashMap<Integer, Integer> comm_to_id_map = new HashMap<Integer, Integer>();
    String newGraphFile = "new_level_graph_" + level_id;
    String newMapFile = "new_level_map_" + level_id;
    int new_num_vertices = lu.prepareForNextLevel(vertices, communities, newGraphFile,
        comm_to_id_map, vertices_to_comm);
    graphFiles.add(new Tuple2<String, String>(newGraphFile, newMapFile));
    levelInfo.add(new Tuple2<HashMap<Integer, Integer>, HashMap<Integer, Integer>>(vertices_to_comm,
        comm_to_id_map));
    if (level_id == 0) {
      this.communities_backup = communities;
      this.num_topics_backup = num_topics;
      this.initial_commweights_backup = initial_commweights;
      this.num_vertices_backup = vertices.length;
    }
    this.num_topics = 1;
    this.initial_commweights = new double[1];
    this.initial_commweights[0] = Main.initWeight;
    vertices = new Vertex[new_num_vertices];
    communities = new Community[new_num_vertices];
  }

  // Spread information across all the communities associated with a vertex id
  // based on community that this vertex belong to
  // Not using it for now. since we decided to change the algorithm to not
  // spread
  // community weights rather update the edge weights based on the new community
  // weights
  // computed
  // private void spreadCommunityInfo(int vertexIdComm, int belongingComm) {
  // for (int i = 0; i < communities[belongingComm].topic_weights.size(); i++) {
  // communities[vertexIdComm].topic_weights.set(i,
  // communities[belongingComm].topic_weights.get(i));
  // }
  // }

  private void deactivateComm() {
    for (int cNum = 0; cNum < communities.length; cNum++) {
      Community c = communities[cNum];
      c.isActive = false;
      for (int i = 0; i < c.topic_weights.size(); i++) {
        c.topic_weights.set(i, Main.initWeight);
      }
    }
    for (Vertex v : vertices) {
      // spreadCommunityInfo(v.id, v.community_label);
      v.community_label = v.id;
    }
  }

  private void mainIter(int level_id) throws IOException {
    deactivateComm();
    for (int i = 0; i < 40; i++) {
      louvainIter(i);
    }
    customizedFinalize();
  }

  // should be according to a feature/nodeattribute
  public double computeGainAll() {
    prepareForLouvainIters();
    double gain = 0.0;
    for (Vertex v : vertices) {
      // look at all the target vertices and update the community data structure
      // if the target vertex is in same
      for (Neighbor n : v.neighbors) {
        if (v.community_label != vertices[n.neighbour_label].community_label) {
          continue;
        }
        // compute gain and gradient at the same time
        double gain_now = (double) n.edge_weight - (double) ((double) (v.total_vertex_degree
            * vertices[n.neighbour_label].total_vertex_degree)
            / (double) ((double) 2 * total_edgeweights.total_edge_weight));
        // gain_now = (double)(gain_now/total_edgeweights.total_edge_weight);
        gain += gain_now;
      }
    }
    double reg_term = 0.0;
    for (Community c : communities) {
      if (!c.isActive)
        continue;
      for (double weight : c.topic_weights) {
        reg_term += weight * weight;
      }
    }
    reg_term = Math.sqrt(reg_term);
    gain -= Main.lambda * reg_term;
    return gain;
  }

  public ArrayList<Double> computeGradientOfCommunity(int comm) {
    // prepareForLouvainIters();
    ArrayList<Double> gradients = new ArrayList<Double>();
    for (int k = 0; k < num_topics; k++) {
      gradients.add(0.0);
    }

    for (Vertex v : community_to_vertices_index.get(comm)) {
      // look at all the target vertices and update the community data structure
      // if the target vertex is in same
      for (Neighbor n : v.neighbors) {
        if (v.community_label != vertices[n.neighbour_label].community_label) {
          continue;
        }

        if (v.community_label != comm) {
          continue;
        }

        for (int k = 0; k < num_topics; k++) {
          double gradient = gradients.get(k);
          // System.out.println("weigh sim = "+n.weighted_similarity.get(k)+",
          // feat_vert_deg = "+v.feature_vertex_degree.get(k)+
          // "total_edgeweights =
          // "+total_edgeweights.feature_edge_weight.get(k));

          // compute gain and gradient at the same time. whatever topic you
          // consider, you will have a similar gain.

          if (total_edgeweights.feature_edge_weight_derivative_helper.get(k) == 0.0) {
            gradient += 0.0;
            continue;
          }

          gradient += n.weighted_similarity_der.get(k)
              - (double) ((double) (v.feature_vertex_derivate_helper.get(k)
                  * vertices[n.neighbour_label].feature_vertex_derivate_helper.get(k))
                  / (double) ((double) 2.0
                      * total_edgeweights.feature_edge_weight_derivative_helper.get(k)));
          gradients.set(k, gradient);
        }
      }
    }

    for (int k = 0; k < gradients.size(); k++) {
      double gradient = gradients.get(k);
      if (total_edgeweights.feature_edge_weight_derivative_helper.get(k) == 0.0) {
        continue;
      }
      // gradients.set(k,
      // (double)(gradient/total_edgeweights.feature_edge_weight_derivative_helper.get(k)));
      gradients.set(k, (double) (gradient));
    }

    return gradients;
  }

  public double computeGainOfCommunity(int comm) {
    // prepareForLouvainIters();
    double gain = 0.0;
    for (Vertex v : vertices) {
      // look at all the target vertices and update the community data structure
      // if the target vertex is in same
      for (Neighbor n : v.neighbors) {
        if (v.community_label != vertices[n.neighbour_label].community_label) {
          continue;
        }

        gain += (double) n.edge_weight - (double) ((double) (v.total_vertex_degree
            * vertices[n.neighbour_label].total_vertex_degree)
            / (double) ((double) 2 * total_edgeweights.total_edge_weight));

      }
    }
    double reg_term = 0.0;
    for (int i = 0; i < communities[comm].topic_weights.size(); i++) {
      reg_term += (communities[comm].topic_weights.get(i) * communities[comm].topic_weights.get(i));
    }
    reg_term = Math.sqrt(reg_term);
    gain -= Main.lambda * reg_term;
    return gain;
    // return (gain/total_edgeweights.total_edge_weight);
  }

  private void regenerateGraph() throws IOException {
    // Now replace the communities with community_backup and recreate
    // community_backup
    // assign correct community labels.
    this.num_topics = num_topics_backup;
    vertices = new Vertex[this.num_vertices_backup];
    communities = new Community[this.num_vertices_backup];
    this.initial_commweights = initial_commweights_backup;
    init(graphFiles.get(0)._1);
    this.communities = this.communities_backup;
    // Now re-assign the communities

    // 1. get the first community label. find associated mapped digit which is
    // input to step 2
    // 2. get the second level community label and associated mapped digit and
    // input to step 3..
    //
    // at the end you will have a label
    // TODO: Check the concern. It should be according to the the actual
    // community label of first
    // iteration instead of the mapped community label of last iteration. Not
    // sure it's a
    // big concern since these mapped labels are the labels which will be used
    // next iteration
    // onwards. so we should be good to go but this requires deeper thinking.
    for (Vertex v : vertices) {
      int comm_label_for_me = v.community_label;
      for (Tuple2<HashMap<Integer, Integer>, HashMap<Integer, Integer>> t : levelInfo) {
        comm_label_for_me = t._1.get(comm_label_for_me);
        comm_label_for_me = t._2.get(comm_label_for_me);
      }
      v.community_label = comm_label_for_me;
    }
    this.levelInfo = new ArrayList<Tuple2<HashMap<Integer, Integer>, HashMap<Integer, Integer>>>();
    Tuple2<String, String> firstGraphFile = this.graphFiles.get(0);
    this.graphFiles = new ArrayList<Tuple2<String, String>>();
    this.graphFiles.add(firstGraphFile);
    customizedFinalize();
  }

  public void startProcessWithPhases(int phase) throws IOException {
    // Check the logic for stopping with modularity computation
    for (int outIter = 0; outIter < Main.outIters; outIter++) {
      if (!(phase > 0 && outIter == 0)) {
        if (outIter != 0) {
          normalizeInputData(graphFiles.get(outIter)._1);
        }
        // init(graphFiles.get(outIter)._1);
      }

      // FIXME: Careful. this change is specifically for processing different
      // optimization
      init(graphFiles.get(outIter)._1);
      // mergeVertices(); FIXME: Going sequential
      mainIterWithPhases(outIter);

    }

    regenerateGraph();
  }

  public void startProcess() throws IOException {
    mainIter(0);
  }

  private void normalizeEdgeWeights() {
    double divisor = 0.0;
    for (Vertex v : vertices) {
      for (Neighbor n : v.neighbors) {
        for (int i = 0; i < n.similarity.size(); i++) {
          divisor += n.similarity.get(i) * n.similarity.get(i);
        }
      }
    }
    divisor = Math.sqrt(divisor);
    if (Double.compare(divisor, 0.0) == 0) {
      return;
    }
    for (Vertex v : vertices) {
      for (Neighbor n : v.neighbors) {
        for (int i = 0; i < n.similarity.size(); i++) {
          double curVal = n.similarity.get(i);
          n.similarity.set(i, (double) (curVal / divisor));
        }
      }
    }
  }

  private double maxGain = 0.0;

  private double computeModularity() {
    double gain = 0.0;
    for (Vertex v : vertices) {
      // look at all the target vertices and update the community data structure
      // if the target vertex is in same
      for (Neighbor n : v.neighbors) {

        if (v.community_label != vertices[n.neighbour_label].community_label) {
          continue;
        }
        gain += (double) n.edge_weight - (double) ((double) (v.total_vertex_degree
            * vertices[n.neighbour_label].total_vertex_degree)
            / (double) ((double) 2 * total_edgeweights.total_edge_weight));
        
      }
    }
    if (Double.compare(gain, maxGain) > 0) {
      maxGain = gain/total_edgeweights.total_edge_weight;
    }
    //System.out.println(">>>>>>"+total_edgeweights.total_edge_weight);
    return gain/total_edgeweights.total_edge_weight;
  }

  public void printCommAffiliations() throws Exception {
    double gain = computeModularity();
    HashMap<Integer, ArrayList<Integer>> commToVertexMap =
        new HashMap<Integer, ArrayList<Integer>>();
    for (Vertex v : vertices) {
      System.out.println(v.id + "\t" + v.community_label);
      for (Neighbor n : v.neighbors) {
        for (int i = 0; i < n.similarity.size(); i++) {
          Community c = communities[v.community_label];
          Community target = communities[n.neighbour_label];
          n.similarity.set(i,
              (c.topic_weights.get(i) + target.topic_weights.get(i)) * n.similarity.get(i));
        }
      }
      ArrayList<Integer> verticestemp;
      if (commToVertexMap.containsKey(v.prev_iter_community_label)) {
        verticestemp = commToVertexMap.get(v.prev_iter_community_label);
      } else {
        verticestemp = new ArrayList<Integer>();
      }
      verticestemp.add(v.id);
      commToVertexMap.put(v.prev_iter_community_label, verticestemp);
    }
    normalizeEdgeWeights();
    // for (Community c : communities) {
    // if (!c.isActive) {
    // continue;
    // }
    // System.out.println("for community ");
    // for (double d : c.topic_weights) {
    // System.out.println("\t"+d);
    // }
    // }
    BufferedWriter bw = new BufferedWriter(new FileWriter("foundLabels" + extra));
    // BufferedWriter bw2 = new BufferedWriter(new
    // FileWriter("for_nmi_ver"+extra));
    for (int i : commToVertexMap.keySet()) {
      bw.write("" + i + ":");
      double areaScores[] = new double[num_topics];
      for (int vertx : commToVertexMap.get(i)) {
        bw.write("\t" + vertx);
        for (Neighbor n : vertices[vertx].neighbors) {
          if (vertices[vertx].community_label != vertices[n.neighbour_label].community_label)
            continue;
          for (int simno = 0; simno < n.similarity.size(); simno++) {
            areaScores[simno] += n.similarity.get(simno);
          }
        }
        // bw2.write(""+i+" ");
      }
      if (Main.printAreaScore) {
        for (double areaScore : areaScores) {
          bw.write("\t" + areaScore);
        }
      }
      double max = 0.0;
      for (int iter = 0; iter < communities[i].topic_weights.size(); iter++) {
        // bw.write("\t"+communities[i].topic_weights.get(iter));
        if (Double.compare(communities[i].topic_weights.get(iter), max) > 0) {
          max = communities[i].topic_weights.get(iter);
        }
      }
      bw.newLine();
    }
    bw.close();
  }

  public void customizedFinalize() throws IOException {
    vertex_to_community_index = new HashMap<Integer, Community>();
    community_to_vertices_index = new HashMap<Integer, ArrayList<Vertex>>();
    community_info = new HashMap<Integer, Tuple2<ArrayList<Vertex>, ArrayList<Vertex>>>();
    for (Vertex v : vertices) {
      v.prev_iter_community_label = v.community_label;
      vertex_to_community_index.put(v.community_label, communities[v.community_label]);
      if (community_to_vertices_index.containsKey(v.community_label)) {
        community_to_vertices_index.get(v.community_label).add(v);
        community_info.get(v.community_label)._1.add(v);
      } else {
        ArrayList<Vertex> averts = new ArrayList<Vertex>();
        ArrayList<Vertex> in_vertex = new ArrayList<Vertex>();
        ArrayList<Vertex> out_vertex = new ArrayList<Vertex>();
        in_vertex.add(v);
        averts.add(v);
        community_to_vertices_index.put(v.community_label, averts);
        community_info.put(v.community_label, new Tuple2<ArrayList<Vertex>, ArrayList<Vertex>>(in_vertex, out_vertex));
      }
      for (Neighbor neigh : v.neighbors) {
        int n_c_label = vertices[neigh.neighbour_label].community_label;
        if (n_c_label == v.community_label)
          continue;
        community_info.get(v.community_label)._2.add(vertices[neigh.neighbour_label]);
      }
      communities[v.prev_iter_community_label].isActive = true;
    }

    // printCommAffiliations();
  }
}
