This repository consists of the implementation of knowledge graph enhanced community detection algorithm.

The details of the algorithm can be found in the [publication](https://dl.acm.org/citation.cfm?id=3291031).

The implementation requires a domain-specific knowledge graph as an input. For a given domain-specific knowledge graph,
it first computes an updated edge list and then runs the community detection algorithm.

The current implementation includes an example that can be run by,

ant run-kg-eval
ant run-community-evaluation

Please change the input files in the KGOptimizer.java and Main.java to run the algorithm with a different dataset.
