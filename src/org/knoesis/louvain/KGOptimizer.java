package org.knoesis.louvain;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import scala.Tuple2;

class KOptimizer {
  ArrayList<ArrayList<Tuple2<String, Integer>>> authorDesc;
  ArrayList<HashMap<String, Tuple2<HashSet<String>, Integer>>> index;
  HashMap<Integer, ArrayList<Integer>> graphData;
  public static String[] domains = //{"<http://dbpedia.org/resource/Category:Sports_in_the_United_States>"};//,
    {"<http://dbpedia.org/resource/Category:American_music>"};
  private static String[] files = //{"sports_subgraph.ttl"};//,
  {"music_subgraph.ttl"};

  //public static String[] domains = {"<USA>"};
  //private static String[] files = {"test_kg.ttl"};


  //  public static String[] domains =
  //      {"<http://dbpedia.org/resource/Category:Companies_by_country_and_industry>",
  //          "<http://dbpedia.org/resource/Category:Occupations>",
  //          "<http://dbpedia.org/resource/Category:Countries>",
  //          "<http://dbpedia.org/resource/Category:Universities_and_colleges_by_country>"};
  //
  //  private static String[] files = {"gplus_domains/institutions.ttl",
  //      "gplus_domains/occupations.ttl", "gplus_domains/place.ttl", "gplus_domains/university.ttl"};

  private double getConceptScore(int parentLevel, int depth, int no_of_mentions) {
    double plevel = (-1.0)* (1.0/(Math.log(parentLevel) / Math.log(2.0)));
    //return (double) ((double) (no_of_mentions) / ((double) parentLevel));
    return (double) ((double) (no_of_mentions) * ((double) plevel));
  }

  public void readLines(String fileName) throws Exception {
    authorDesc = new ArrayList<ArrayList<Tuple2<String, Integer>>>();
    List<String> lines = Files.readAllLines(Paths.get(fileName));
    for (String line : lines) {
      String[] splits = line.split(",");
      ArrayList<Tuple2<String, Integer>> concepts = new ArrayList<Tuple2<String, Integer>>();
      for (int i = 1; i < splits.length; i++) {
        String inSp[] = splits[i].split(" ");
        if (!inSp[0].endsWith(">"))
          continue;
        boolean to_consider = false;
        for (HashMap<String, Tuple2<HashSet<String>, Integer>> ind : index) {
          if (ind.containsKey(inSp[0])) {
            to_consider = true;
            break;
          }
        }
        if (to_consider)
          concepts.add(new Tuple2<String, Integer>(inSp[0], Integer.parseInt(inSp[1])));
      }
      if (concepts.size() == 0) {
        for (String domain : domains) {
          concepts.add(new Tuple2<String, Integer>(domain, 1));
        }
      }
      authorDesc.add(concepts);
    }
  }

  public double getThisConceptScore(int concept_count, int expected, int total) {
    return (double)((double)concept_count - (double)((double)1.0*(double)expected/(double)total));
  }

  public HashMap<String, Double> scoreConcepts(HashMap<String, Integer> existingList,
      HashMap<String, Tuple2<HashSet<String>, Integer>> index, int total) {
    HashMap<String, Double> concepts_and_score = new HashMap<String, Double>();
    HashMap<String, Integer> this_level_list = existingList;
    HashMap<String, Integer> next_level_list = new HashMap<String, Integer>();

    for (int level = 0; level < ConceptIndexerAuthorScorer.maxParentLevel; level++) {
      for (String concept_for_this : this_level_list.keySet()) {
        if (!index.containsKey(concept_for_this)) {
          //System.out.println(concept_for_this);
          continue;
        }
        double this_concept_score = getThisConceptScore(this_level_list.get(concept_for_this), 
            index.get(concept_for_this)._2, total);
        if (concepts_and_score.containsKey(concept_for_this)) {
          double existing_score = concepts_and_score.get(concept_for_this);
          this_concept_score += existing_score; 
        }
        concepts_and_score.put(concept_for_this, this_concept_score);
        if (index.get(concept_for_this)._1 == null)
          continue;
        for (String parent : index.get(concept_for_this)._1) {
          int mention_ct = this_level_list.get(concept_for_this);
          if (next_level_list.containsKey(parent)) {
            mention_ct += next_level_list.get(parent);
          }
          next_level_list.put(parent, mention_ct);
        }
      }
      this_level_list = next_level_list;
      next_level_list = new HashMap<String, Integer>();
    }
    return concepts_and_score;
  }

  public HashMap<String, Double> getAuthorMultipleLevelList(
      ArrayList<HashMap<String, Double>> newListHierarchy,
      ArrayList<Tuple2<String, Integer>> existingList,
      HashMap<String, Tuple2<HashSet<String>, Integer>> index, int indexId) {
    int currentLevel = 0;
    HashMap<String, Double> finalList = new HashMap<String, Double>();
    for (Tuple2<String, Integer> existing : existingList) {
      currentLevel = 0;
      HashMap<String, Double> actualConceptsList = newListHierarchy.get(currentLevel++);
      if (existing == null) {
        continue;
      }
      if (existing._1 == null) {
        continue;
      }
      if (index.get(existing._1) == null) {
        continue;
      }
      // if (existing._2 != indexId) {
      //  continue;
      // }
      LinkedList<String> checkParents = new LinkedList<String>();
      int no_of_mentions = existing._2;
      checkParents.add(existing._1);
      // add current concept
      double score_for_me = getConceptScore(1, index.get(existing._1)._2, no_of_mentions);
      if (actualConceptsList.containsKey(existing._1)) {
        // The most specific concept defining this author already exists in the
        // hierachy.
        // we add the weight for that concept.
        score_for_me += actualConceptsList.get(existing._1);
      }
      actualConceptsList.put(existing._1, score_for_me);
      // now start parents
      for (int parentLevel =
          2; parentLevel < ConceptIndexerAuthorScorer.maxParentLevel; parentLevel++) {
        if (checkParents.isEmpty()) {
          break;
        }
        LinkedList<String> checkParentTemp = new LinkedList<String>();
        while (!checkParents.isEmpty()) {
          String checkParent = checkParents.removeFirst();
          if (index.get(checkParent)._1 == null) {
            boolean should_allow = false;
            for (String allowedNullParent : ConceptIndexerAuthorScorer.domains) {
              if (checkParent.equals(allowedNullParent)) {
                should_allow = true;
                break;
              }
            }
            if (!should_allow) {
              System.out.println(
                  checkParent + " caused error since it wasn't found while it should in index");
              System.exit(1);
            }
            continue;
          }
          for (String par1 : index.get(checkParent)._1) {
            double score_for_this =
                getConceptScore(parentLevel, index.get(par1)._2, no_of_mentions);
            if (finalList.containsKey(par1)) {
              // if a concept for this author already exist in parent hierachy
              // and we find it again.
              // then we add the score for that concept.
              score_for_this += finalList.get(par1);
            }
            finalList.put(par1, score_for_this);
            checkParentTemp.add(par1);
          }
        }
        checkParents = checkParentTemp;
      }
      // authorCache.put(authorId, newList);
    }
    return finalList;
  }

  public void readAuthorDesc(String authorFile) throws Exception {
    index = new ArrayList<HashMap<String, Tuple2<HashSet<String>, Integer>>>();
    for (int i = 0; i < domains.length; i++) {
      HashMap<String, Tuple2<HashSet<String>, Integer>> index_temp = 
          new HashMap<String, Tuple2<HashSet<String>, Integer>>();
      Indexer indexer = new Indexer();
      indexer.readFile(files[i]);
      indexer.addExpectedConcepts(domains[i], index_temp);
      // System.out.println(index_temp.size() + " " + files[i]);
      index.add(index_temp);
    }
    readLines(authorFile);
  }

  public Tuple2<String, Double> getMax(HashMap<String, Double> ranking, HashMap<String, Double> nRanking) {
    double max = Double.MIN_VALUE;
    String maxConcept = null;
    for (String concept : ranking.keySet()) {
      double myScore = ranking.get(concept);
      myScore = (double)((double)myScore/(double)ranking.size());
      if (nRanking.containsKey(concept)) {
        myScore -= (double)(nRanking.get(concept)/(double)nRanking.size());
      }
      if (Double.compare(myScore, max) > 0) {
        max = (myScore);
        maxConcept = concept;
      }
    }
    if (maxConcept == null) {
      // System.out.println("This community did not have anything related to this domain");
      return new Tuple2<String, Double>("null", 0.0);
    }
    return new Tuple2<String, Double>(maxConcept, max);
  }

  public Tuple2<ArrayList<String>, ArrayList<Double>> processCommunity(ArrayList<Integer> aids, ArrayList<Integer> naids, int total_indexes,
      BufferedWriter bw) throws Exception {
    String[] newDomains = new String[total_indexes];
    for (int i = 0; i < domains.length; i++) {
      newDomains[i] = domains[i];
    }
    ArrayList<Double> community_fit_scores = new ArrayList<Double>();
    ArrayList<String> community_domains = new ArrayList<String>();
    for (int indexId = 0; indexId < total_indexes; indexId++) {
      // All the concepts and their counts in this community
      HashMap<String, Integer> commList = new HashMap<String, Integer>();
      HashMap<String, Integer> nCommList = new HashMap<String, Integer>();
      for (int aid : aids) {
        for (Tuple2<String, Integer> t : authorDesc.get(aid)) {
          if (commList.containsKey(t._1)) {
            int existing = commList.get(t._1);
            existing += t._2;
            commList.put(t._1, existing);
          } else {
            commList.put(t._1, t._2);
          }
        }
      }

      for (int aid : naids) {
        for (Tuple2<String, Integer> t : authorDesc.get(aid)) {
          if (nCommList.containsKey(t._1)) {
            int existing = nCommList.get(t._1);
            existing += t._2;
            nCommList.put(t._1, existing);
          } else {
            nCommList.put(t._1, t._2);
          }
        }
      }

      ArrayList<Tuple2<String, Integer>> temp = new ArrayList<Tuple2<String, Integer>>();
      for (String key : commList.keySet()) {
        temp.add(new Tuple2<String, Integer>(key, commList.get(key)));
      }

      // public void getAuthorMultipleLevelList(ArrayList<HashMap<String,
      // Double>> newListHierarchy, ArrayList<Tuple2<String, Integer>>
      // existingList,
      // HashMap<String, Tuple2<HashSet<String>, Integer>> index, String
      // authorId, int indexId) {
      // ArrayList<HashMap<String, Double>>
      ArrayList<HashMap<String, Double>> authorTotalList = new ArrayList<HashMap<String, Double>>();
      for (int i = 0; i < ConceptIndexerAuthorScorer.maxParentLevel; i++) {
        authorTotalList.add(new HashMap<String, Double>());
      }
      int total = index.get(indexId).get(newDomains[indexId])._2;
      HashMap<String, Double> maxFind =
          scoreConcepts(commList, index.get(indexId), total);
      HashMap<String, Double> minFind = scoreConcepts(nCommList, index.get(indexId), total);
      Tuple2<String, Double> newLevel = getMax(maxFind, minFind);
      newDomains[indexId] = newLevel._1;
      community_fit_scores.add(newLevel._2);
      community_domains.add(newLevel._1);
    }

    ArrayList<HashMap<String, Tuple2<HashSet<String>, Integer>>> newIndex =
        new ArrayList<HashMap<String, Tuple2<HashSet<String>, Integer>>>();
    for (int i = 0; i < newDomains.length; i++) {
      HashMap<String, Tuple2<HashSet<String>, Integer>> index_temp =
          new HashMap<String, Tuple2<HashSet<String>, Integer>>();
      Indexer indexer = new Indexer();
      indexer.readFile(files[i]);
      indexer.addExpectedConcepts(newDomains[i], index_temp);
      // System.out.println("\t\t>>>>"+newDomains[i]);
      // System.out.println(index_temp.size() + " " + files[i]);
      newIndex.add(index_temp);
    }

    // You will be having two directed edges instead of one undirected edge
    AuthorComparator ac = new AuthorComparator();
    for (int aid : aids) {
      for (int target : graphData.get(aid)) {
        ArrayList<Double> four_area_interaction_scores = ac.getAuthorInteraction(newIndex,
            authorDesc.get(aid), authorDesc.get(target), "" + aid, "" + target);
        bw.write("" + aid);
        for (double d : four_area_interaction_scores) {
          bw.write(" " + (double) (d));
        }
        bw.write(" " + target);
        bw.newLine();
      }
    }
    return new Tuple2<ArrayList<String>, ArrayList<Double>>(community_domains, community_fit_scores);
  }

  private void addInfo(HashMap<Integer, ArrayList<Integer>> graphData, int sid, int tid) {
    if (graphData.containsKey(sid)) {
      graphData.get(sid).add(tid);
    } else {
      ArrayList<Integer> edges = new ArrayList<Integer>();
      edges.add(tid);
      graphData.put(sid, edges);
    }
  }

  public void readAdjacencyList(String graphFile) throws Exception {
    List<String> lines = Files.readAllLines(Paths.get(graphFile));
    graphData = new HashMap<Integer, ArrayList<Integer>>();
    for (String line : lines) {
      String splits[] = line.split(" ");
      int sid = Integer.parseInt(splits[0]);
      int tid = Integer.parseInt(splits[splits.length - 1]);
      addInfo(graphData, sid, tid);
      //if (sid != tid)
      //  addInfo(graphData, tid, sid);
    }
  }

  public void readAndProcessAuthorComm(HashMap<Integer, Tuple2<ArrayList<Vertex>, ArrayList<Vertex>>> commToVertex,
      int num_topics, String newGraphFile, String finalOutput) throws Exception {
    BufferedWriter bw = new BufferedWriter(new FileWriter(newGraphFile));
    BufferedWriter bwFinal = new BufferedWriter(new FileWriter(finalOutput));
    int commId = 0;
    for (int c : commToVertex.keySet()) {
      Tuple2<ArrayList<Vertex>, ArrayList<Vertex>> t = commToVertex.get(c);
      ArrayList<Integer> authorsInComm = new ArrayList<Integer>();
      ArrayList<Integer> authorsInNeighbors = new ArrayList<Integer>();
      bwFinal.write("" + commId++);
      for (Vertex v : t._1) {
        authorsInComm.add(v.id);
        bwFinal.write("\t" + v.id);
      }
      for (Vertex v : t._2) {
        authorsInNeighbors.add(v.id);
      }
      Tuple2<ArrayList<String>, ArrayList<Double>> communityFitScore = processCommunity(authorsInComm, authorsInNeighbors, num_topics, bw);
      for (int d = 0; d < communityFitScore._1.size(); d++) {
        bwFinal.write("\t"+communityFitScore._1.get(d)+" = "+communityFitScore._2.get(d));
      }
      bwFinal.newLine();
    }
    bwFinal.close();
    bw.close();
  }

}


public class KGOptimizer {

  public void optimizeKnowledgeGraph(String authorFile, int num_topics, String graphFile,
      HashMap<Integer, Tuple2<ArrayList<Vertex>, ArrayList<Vertex>>> authorCommFile, String newGraphFile, String finalOutput)
          throws Exception {
    KOptimizer ko = new KOptimizer();
    ko.readAdjacencyList(graphFile);
    ko.readAuthorDesc(authorFile);

    ko.readAndProcessAuthorComm(authorCommFile, num_topics, newGraphFile, finalOutput);
  }

  public static void main(String[] args) throws Exception {
    String authorFile = "sorted_user_concepts.csv";
    String graphFile = "conv_net.txt";
    String finalOutput = "resfile.txt";
    int num_topics = 1;
    HashMap<Integer, Tuple2<ArrayList<Vertex>, ArrayList<Vertex>>> authorCommFile =
        new HashMap<Integer, Tuple2<ArrayList<Vertex>, ArrayList<Vertex>>>();
    authorCommFile.put(0, new Tuple2<ArrayList<Vertex>, ArrayList<Vertex>>(new ArrayList<Vertex>(), new ArrayList<Vertex>()));
    //authorCommFile.put(1, new Tuple2<ArrayList<Vertex>, ArrayList<Vertex>>(new ArrayList<Vertex>(), new ArrayList<Vertex>()));
    for (int i = 0; i < 386; i++) {
      Vertex v1 = new Vertex();
      v1.id = i;
      authorCommFile.get(0)._1.add(v1);
      //if (i == 2) {
      //  authorCommFile.get(1)._2.add(v1);
      //}
    }
    //    for (int i = 3; i < 6; i++) {
    //      Vertex v1 = new Vertex();
    //      v1.id = i;
    //      authorCommFile.get(1)._1.add(v1);
    //      if (i == 3) {
    //        authorCommFile.get(0)._2.add(v1);
    //      }
    //    }
    String newGraphFile = "processed_file.txt";
    KOptimizer ko = new KOptimizer();
    ko.readAdjacencyList(graphFile);
    ko.readAuthorDesc(authorFile);
    ko.readAndProcessAuthorComm(authorCommFile, num_topics, newGraphFile, finalOutput);
  }
}