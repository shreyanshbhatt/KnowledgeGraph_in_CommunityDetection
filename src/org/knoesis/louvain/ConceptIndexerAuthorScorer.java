package org.knoesis.louvain;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import scala.Tuple2;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

class Indexer {
  HashMap<String, HashSet<String>> conceptParentIndex;
  HashMap<String, HashSet<String>> reverseParentIndex;
  public static HashSet<String> appearedConcepts = new HashSet<String>();
  public int getConceptScore(String concept, HashMap<String, Tuple2<HashSet<String>, Integer>> index) {
    if (index.containsKey(concept)) {
      return index.get(concept)._2;
    }
    HashSet<String> children = reverseParentIndex.get(concept);
    if (children == null) {
      return 0;
    }
    // myscore = all_children_score addition
    int my_score = 0;
    for (String child : children) {
      if (appearedConcepts.contains(child)) {
        continue;
      }
      int child_score = 0;
      appearedConcepts.add(child);
      child_score = getConceptScore(child, index);
      index.put(child, new Tuple2<HashSet<String>, Integer>(conceptParentIndex.get(child), child_score + 1));
      my_score += child_score + 1;
    }
    return my_score;
  }


  public void addExpectedConcepts(String inputConcept,
      HashMap<String, Tuple2<HashSet<String>, Integer>> index) {
    int root_score = getConceptScore(inputConcept, index);
    index.put(inputConcept, new Tuple2<HashSet<String>, Integer>(null, root_score + 1));
    conceptParentIndex.clear();
    reverseParentIndex.clear();
    appearedConcepts.clear();
    appearedConcepts = new HashSet<String>();
    // System.out.println("index size == "+index.size());
  }


  public void addDepth(String inputConcept,
      HashMap<String, Tuple2<HashSet<String>, Integer>> index) {
    LinkedList<String> queue = new LinkedList<String>();
    LinkedList<String> newQueue = new LinkedList<String>();
    queue.add(inputConcept);
    int level = 2;
    index.put(inputConcept, new Tuple2<HashSet<String>, Integer>(null, 1));
    while (!queue.isEmpty()) {
      String curConcept = queue.removeFirst();
      HashSet<String> children = reverseParentIndex.get(curConcept);
      if (children == null) {
        if (queue.isEmpty()) {
          queue = newQueue;
          newQueue = new LinkedList<String>();
          level++;
        }
        continue;
      }
      for (String child : children) {
        if (index.containsKey(child) && index.get(child)._2() <= level) {
          continue;
        }
        index.put(child,
            new Tuple2<HashSet<String>, Integer>(null, level));
        newQueue.addLast(child);
      }
      if (queue.isEmpty()) {
        queue = newQueue;
        newQueue = new LinkedList<String>();
        level++;
      }
    }
    conceptParentIndex.clear();
    reverseParentIndex.clear();
  }

  public void print(HashMap<String, Tuple2<HashSet<String>, Integer>> index) {
    int total = 0;
    for (String s : index.keySet()) {
      if (index.get(s)._2 == 5) {
        System.out.println(s);
        total++;
      }
    }
    System.out.println(total);
  }

  private void addForwardParent(String sub, String pred, String obj,
      HashMap<String, HashSet<String>> conceptParentIndex) {
    if (pred.equals("<http://dbpedia.org/ontology/wikiPageWikiLink>")) {
      HashSet<String> parentsForMe = null;
      if (conceptParentIndex.containsKey(obj)) {
        // we have already seen this object and we will be adding to its
        // parents.
        parentsForMe = conceptParentIndex.get(obj);
      } else {
        // it's the first time we encountered this object
        parentsForMe = new HashSet<String>();
      }

      // It's the <sub> here that should be the parent but since we consider the
      // article mentioned as the same level as article we will include the
      // artcile's parents
      // as the parents.
      if (conceptParentIndex.get(sub) == null) {
        // This is a very special case. our root itself had some articles listed
        // in itself so we won't be
        // able to find parents for this. in this case just add root.
        parentsForMe.add(sub);
      } else {
        parentsForMe.addAll(conceptParentIndex.get(sub));
      }
      conceptParentIndex.put(obj, parentsForMe);
    } else if (pred.equals("<http://purl.org/dc/terms/subject>")
        || pred.equals("<http://www.w3.org/2004/02/skos/core#broader>")) {
      HashSet<String> existing = null;
      if (conceptParentIndex.containsKey(sub)) {
        existing = conceptParentIndex.get(sub);
      } else {
        existing = new HashSet<String>();
      }
      existing.add(obj);
      conceptParentIndex.put(sub, existing);
    } else {
      System.out.println("Issue...");
      System.exit(1);
    }
  }

  private void addReverseParent(String sub, String pred, String obj,
      HashMap<String, HashSet<String>> reverseParentIndex) {
    if (pred.equals("<http://dbpedia.org/ontology/wikiPageWikiLink>")) {
      HashSet<String> children = null;
      if (reverseParentIndex.containsKey(sub)) {
        children = reverseParentIndex.get(sub);
      } else {
        children = new HashSet<String>();
      }
      children.add(obj);
      reverseParentIndex.put(sub, children);
    } else if (pred.equals("<http://purl.org/dc/terms/subject>")
        || pred.equals("<http://www.w3.org/2004/02/skos/core#broader>")) {
      HashSet<String> existing_children = null;
      if (reverseParentIndex.containsKey(obj)) {
        existing_children = reverseParentIndex.get(obj);
      } else {
        existing_children = new HashSet<String>();
      }
      existing_children.add(sub);
      reverseParentIndex.put(obj, existing_children);
    } else {
      System.out.println("Error");
      System.exit(1);
    }
  }

  public void readFile(String file) throws Exception {
    conceptParentIndex = new HashMap<String, HashSet<String>>();
    reverseParentIndex = new HashMap<String, HashSet<String>>();
    HashSet<String> encountered = new HashSet<String>();
    List<String> lines = Files.readAllLines(Paths.get(file));
    for (String line : lines) {
      String[] splits = line.split(" ");
      String sub = splits[0];
      String pred = splits[1];
      String obj = splits[2];
      encountered.add(sub);
      encountered.add(obj);
      addForwardParent(sub, pred, obj, conceptParentIndex);
      addReverseParent(sub, pred, obj, reverseParentIndex);
    }
  }

}


class AuthorComparator {
  HashMap<String, HashMap<String, Double>> authorCache;

  public AuthorComparator() {
    authorCache = new HashMap<String, HashMap<String, Double>>();
  }

  public void readLines(String fileName, ArrayList<ArrayList<Tuple2<String, Integer>>> authorDesc)
      throws Exception {
    List<String> lines = Files.readAllLines(Paths.get(fileName));
    for (String line : lines) {
      String[] splits = line.split("\t");
      ArrayList<Tuple2<String, Integer>> concepts = new ArrayList<Tuple2<String, Integer>>();

      String conceptsStr = splits[1];
      String innerSplits[] = conceptsStr.split("\\*");
      for (int i = 0; i < innerSplits.length; i++) {
        String concetSplits[] = innerSplits[i].split(" ");
        for (String conceptSplit : concetSplits) {
          concepts.add(new Tuple2<String, Integer>(conceptSplit, i));
        }
      }

      authorDesc.add(concepts);
    }
  }

  private double getJaccardWeightDistane(double d1, double d2) {
    return (double) (Math.min(d1, d2) * (1.0 - (Math.abs(d1 - d2)) / (d1 + d2)));
  }

  private Tuple2<Double, Double> compute_wighted_jaccard(ArrayList<HashMap<String, Double>> author1,
      ArrayList<HashMap<String, Double>> author2) {
    double denom = 1.0;
    double numer = 0.0;
    if (author1.size() == 0 || author2.size() == 0) {
      // Authors were not able to be described in this domain
      return new Tuple2<Double, Double>(0.0, 0.0);
    }
    for (HashMap<String, Double> author1Matched : author1) {
      for (HashMap<String, Double> author2Matched : author2) {
        for (String concept : author1Matched.keySet()) {
          if (author2Matched.containsKey(concept)) {
            numer +=
                getJaccardWeightDistane(author1Matched.get(concept), author2Matched.get(concept));
          }
        }
      }
    }

    for (String concept1 : author1.get(0).keySet()) {
      for (String concept2 : author2.get(0).keySet()) {
        denom +=
            getJaccardWeightDistane(author1.get(0).get(concept1), author2.get(0).get(concept2));
      }
    }
    return new Tuple2<Double, Double>(numer, denom);
  }

  private double getConceptScore(int parentLevel, int depth, int no_of_mentions) {
    return (double) ((double) (no_of_mentions) / ((double) parentLevel * (double) depth));
  }

  private void getAuthorMultipleLevelList(ArrayList<HashMap<String, Double>> newListHierarchy,
      ArrayList<Tuple2<String, Integer>> existingList,
      HashMap<String, Tuple2<HashSet<String>, Integer>> index, String authorId, int indexId) {
    int currentLevel = 0;
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
      //if (existing._2 != indexId) {
      //  continue;
      //}
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
        HashMap<String, Double> thisLevelList = newListHierarchy.get(parentLevel);
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
              //System.out.println(
              //    checkParent + " caused error since it wasn't found while it should in index");
              //System.exit(1);
            }
            continue;
          }
          for (String par1 : index.get(checkParent)._1) {
            if (!index.containsKey(par1))
              continue;
            double score_for_this =
                getConceptScore(parentLevel, index.get(par1)._2, no_of_mentions);
            if (thisLevelList.containsKey(par1)) {
              // if a concept for this author already exist in parent hierachy
              // and we find it again.
              // then we add the score for that concept.
              score_for_this += thisLevelList.get(par1);
            }
            thisLevelList.put(par1, score_for_this);
            checkParentTemp.add(par1);
          }
        }
        checkParents = checkParentTemp;
      }
      // authorCache.put(authorId, newList);
    }
  }

  private Tuple2<Double, Double> getAreaScore(
      HashMap<String, Tuple2<HashSet<String>, Integer>> index,
      ArrayList<Tuple2<String, Integer>> author1, ArrayList<Tuple2<String, Integer>> author2,
      String aid1, String aid2, int indexNum) {
    ArrayList<HashMap<String, Double>> authorTotalList1 = new ArrayList<HashMap<String, Double>>();
    ArrayList<HashMap<String, Double>> authorTotalList2 = new ArrayList<HashMap<String, Double>>();
    for (int i = 0; i < ConceptIndexerAuthorScorer.maxParentLevel; i++) {
      authorTotalList1.add(new HashMap<String, Double>());
      authorTotalList2.add(new HashMap<String, Double>());
    }
    getAuthorMultipleLevelList(authorTotalList1, author1, index, aid1, indexNum);
    getAuthorMultipleLevelList(authorTotalList2, author2, index, aid2, indexNum);
    return compute_wighted_jaccard(authorTotalList1, authorTotalList2);
  }

  public ArrayList<Double> getAuthorInteraction(
      ArrayList<HashMap<String, Tuple2<HashSet<String>, Integer>>> index,
      ArrayList<Tuple2<String, Integer>> author1, ArrayList<Tuple2<String, Integer>> author2,
      String aid1, String aid2) {
    ArrayList<Double> four_area_score = new ArrayList<Double>();
    ArrayList<Double> numers = new ArrayList<Double>();
    ArrayList<Double> denoms = new ArrayList<Double>();
    double denom = 0.0;
    for (int indexId = 0; indexId < index.size(); indexId++) {
      Tuple2<Double, Double> numer_denom =
          getAreaScore(index.get(indexId), author1, author2, aid1, aid2, indexId);
      numers.add(numer_denom._1);
      denoms.add(numer_denom._2);
      denom += numer_denom._2;
    }
    for (int i = 0; i < numers.size(); i++) {
      if (denom == 0.0) {
        four_area_score.add(0.0);
        continue;
      }
      four_area_score.add((double) (numers.get(i) / denoms.get(i)));
    }
    return four_area_score;
  }

  public void computeAuthorDistance(String authorEdgeFile,
      ArrayList<HashMap<String, Tuple2<HashSet<String>, Integer>>> index,
      ArrayList<ArrayList<Tuple2<String, Integer>>> authorDesc) throws Exception {
    List<String> lines = Files.readAllLines(Paths.get(authorEdgeFile));
    BufferedWriter bw = new BufferedWriter(new FileWriter(authorEdgeFile + "_annotated"));
    for (String line : lines) {
      String splits[] = line.split(" ");
      int author1 = Integer.parseInt(splits[0]);
      int author2 = Integer.parseInt(splits[1]);
      ArrayList<Double> four_area_interaction_scores = getAuthorInteraction(index,
          authorDesc.get(author1), authorDesc.get(author2), "" + author1, "" + author2);
      bw.write("" + author1);
      for (double d : four_area_interaction_scores) {
        bw.write(" " + (double) (d));
      }
      bw.write(" " + author2);
      bw.newLine();
      System.out.println("done " + author1 + " " + author2);
    }
    bw.close();
  }

}


public class ConceptIndexerAuthorScorer {
  public static String[] domains =
    {"<http://dbpedia.org/resource/Category:Companies_by_country_and_industry>",
        "<http://dbpedia.org/resource/Category:Occupations>",
        "<http://dbpedia.org/resource/Category:Countries>",
    "<http://dbpedia.org/resource/Category:Universities_and_colleges_by_country>"};
  private static String[] files = {"gplus_domains/institutions.ttl",
      "gplus_domains/occupations.ttl", "gplus_domains/place.ttl", "gplus_domains/university.ttl"};
  public static int maxParentLevel = 4;

  public static void mainHelper(String[] newDomains) throws Exception {
    domains = newDomains;
    ArrayList<HashMap<String, Tuple2<HashSet<String>, Integer>>> index =
        new ArrayList<HashMap<String, Tuple2<HashSet<String>, Integer>>>();
    for (int i = 0; i < domains.length; i++) {
      HashMap<String, Tuple2<HashSet<String>, Integer>> index_temp =
          new HashMap<String, Tuple2<HashSet<String>, Integer>>();
      Indexer indexer = new Indexer();
      indexer.readFile(files[i]);
      indexer.addDepth(domains[i], index_temp);
      System.out.println(index_temp.size() + " " + files[i]);
      index.add(index_temp);
    }
    for (final File fileEntry : new File("./").listFiles()) {
      if (!fileEntry.getName().endsWith(".edges")) {
        continue;
      }
      String authorId = fileEntry.getName().split("\\.")[0];
      AuthorComparator ap = new AuthorComparator();
      ArrayList<ArrayList<Tuple2<String, Integer>>> authorDescriptions =
          new ArrayList<ArrayList<Tuple2<String, Integer>>>();
      ap.readLines(authorId + ".featuresdb", authorDescriptions);
      // ap.readLines(authorDescbi, authorDescriptions);
      ap.computeAuthorDistance(fileEntry.getAbsolutePath(), index, authorDescriptions);
    }
  }

  public static void main(String[] args) throws Exception {
    mainHelper(domains);
  }
}
