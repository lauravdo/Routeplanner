package graphs;

import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;

public class DirectedGraph<V extends DGVertex<E>, E extends DGEdge<V>> {

    private Map<String, V> vertices = new HashMap<>();


    /**
     * representation invariants:
     * 1.  all vertices in the graph are unique by their implementation of the getId() method
     * 2.  all edges in the graph reference vertices from and to which are true members of the vertices map
     * (i.e. by true object instance equality == and not just by identity equality from the getId() method)
     * 3.  all edges of a vertex are outgoing edges, i.e. FOR ALL e in v.edges: e.from == v
     **/

    public DirectedGraph() {
    }

    public Collection<V> getVertices() {
        return this.vertices.values();
    }

    /**
     * finds the vertex in the graph identified by the given id
     *
     * @param id
     * @return the vertex that matches the given id
     * return null if none of the vertices matches the id
     */
    public V getVertexById(String id) {
        return this.vertices.get(id);
    }


    /**
     * Adds newVertex to the graph, if not yet present and in a way that maintains the representation invariants.
     * If (a duplicate of) newVertex (with the same id) already exists in the graph,
     * nothing will be added, and the existing duplicate will be kept and returned.
     *
     * @param newVertex
     * @return the duplicate of newVertex with the same id that already existed in the graph,
     * or newVertex itself if it has been added.
     */
    public V addOrGetVertex(V newVertex) {
        // TODO add and return the newVertex, or return the existing duplicate vertex
        //create new vertex, if there is already one, put it in the already excisiting vertices
        V newV;
        if (!vertices.containsKey(newVertex.getId())) {
            this.vertices.put(newVertex.getId(), newVertex);
        }
        return newVertex;
    }

    /**
     * Adds all newVertices to the graph, which are not present yet and and in a way that maintains the representation invariants.
     *
     * @param newVertices an array of vertices to be added, provided as variable length argument list
     * @return the number of vertices that actually have been added.
     */
    public int addVertices(V... newVertices) {
        int count = 0;
        for (V v : newVertices) {
            if (v == this.addOrGetVertex(v)) {
                count++;
            }
        }

        return count;
    }

    /**
     * Adds newEdge to the graph, if not yet present and in a way that maintains the representation invariants:
     * If any of the newEdge.from or newEdge.to vertices does not yet exist in the graph, it is added now.
     * If newEdge does not exist yet in the edges list of the newEdge.from vertex, it is added now,
     * otherwise no change is made to that list.
     *
     * @param newEdge the new edge to be added in the edges list of newEdge.from
     * @return the duplicate of newEdge that already existed in the graph
     * or newEdge itselves if it just has been added.
     * @throws IllegalArgumentException if newEdge.from or newEdge.to are duplicate vertices that have not
     *                                  been added to the graph yet have the same id as another vertex in the graph
     */
    public E addOrGetEdge(E newEdge) {
        // TODO add and return the newEdge, or return the existing duplicate edge or throw an exception
        //check if vertice already has this edge & throws an error
        if ((vertices.containsValue(newEdge.getFrom()) || vertices.containsValue(newEdge.getTo()))
                && vertices.containsValue(newEdge.getFrom().getId()) || vertices.containsValue(newEdge.getTo().getId())) {
            throw new IllegalArgumentException("has same id as other vertex in graph?(" + newEdge + ")");
        }
        //check if vertice doesnt have the edge and adds it
        if (!vertices.containsValue(newEdge.getFrom())) {
            vertices.put(newEdge.getFrom().getId(), newEdge.getFrom());
        }
        if (!vertices.containsValue(newEdge.getTo())) {
            vertices.put(newEdge.getTo().getId(), newEdge.getTo());
        }
        if (!newEdge.getFrom().getEdges().contains(newEdge)) {
            newEdge.getFrom().getEdges().add(newEdge);

        }
        return newEdge;


        // a proper edge shall be returned at all times

    }

    /**
     * Adds all newEdges to the graph, which are not present yet and in a way that maintains the representation invariants.
     *
     * @param newEdges an array of vertices to be added, provides as variable length argument list
     * @return the number of edges that actually have been added.
     */
    public int addEdges(E... newEdges) {
        int count = 0;
        for (E e : newEdges) {
            if (e == this.addOrGetEdge(e)) {
                count++;
            }
        }

        return count;
    }

    /**
     * @return the total number of vertices in the graph
     */
    public int getNumVertices() {
        return this.vertices.size();
    }

    /**
     * @return the total number of edges in the graph
     */
    public int getNumEdges() {
        // TODO calculate and return the total number of edges in the graph
        //counter starts at 0
        int count = 0;
        //loop through all the edges in the vertices and adds every single one to count
        for (var entry : vertices.entrySet()) {
            count = entry.getValue().getEdges().size();
        }
        return count;
    }

    /**
     * Clean-up unconnected vertices in the graph
     */
    public void removeUnconnectedVertices() {
        Set<V> unconnected = new HashSet<>();

        this.getVertices().stream().filter(v -> v.getEdges().size() == 0).forEach(unconnected::add);
        this.getVertices().stream().flatMap(v -> v.getEdges().stream().map(E::getTo)).forEach(unconnected::remove);
        unconnected.stream().map(V::getId).forEach(this.vertices::remove);
    }

    /**
     * represents a path of connected vertices and edges in the graph
     */
    public class DGPath {
        private V start = null;
        private LinkedList<E> edges = new LinkedList<>();
        private double totalWeight = 0.0;
        private Set<V> visited = new HashSet<>();

        /**
         * representation invariants:
         * 1. The edges are connected by vertices, i.e. FOR ALL i: 0 < i < edges.length: edges[i].from == edges[i-1].to
         * 2. The path begins at vertex == start
         * 3. if edges is empty, the path also ends at vertex == start
         * otherwise edges[0].from == start and the path continues along edges[i].to for all 0 <= i < edges.length
         **/

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder(
                    String.format("Weight=%f Length=%d Visited=%d (",
                            this.totalWeight, 1 + this.edges.size(), this.visited.size()));
            sb.append(start.getId());
            for (E e : edges) {
                sb.append(", " + e.getTo().getId());
            }
            sb.append(")");
            return sb.toString();
        }

        public V getStart() {
            return start;
        }

        public LinkedList<E> getEdges() {
            return edges;
        }

        public double getTotalWeight() {
            return totalWeight;
        }

        public Set<V> getVisited() {
            return visited;
        }
    }

    /**
     * Uses a depth-first search algorithm to find a path from the start vertex to the target vertex in the graph
     * The path.totalWeight should indicate the number of edges in the result path
     * All vertices that are being visited by the search should also be registered in path.visited
     *
     * @param startId
     * @param targetId
     * @return the path from start to target
     * returns null if either start or target cannot be matched with a vertex in the graph
     * or no path can be found from start to target
     */
    public DGPath depthFirstSearch(String startId, String targetId) {

        V start = this.getVertexById(startId);
        V target = this.getVertexById(targetId);
        if (start == null || target == null) return null;

        DGPath path = new DGPath();
        //start of the path
        path.start = start;


        // easy target, found path
        if (start == target) {
            path.visited.add(start);
            return path;
        }
        //recursive call
        return DFSRecursion(start, target, path);


    }

    private DGPath DFSRecursion(V currentVertice, V targetVertice, DGPath path) {

        //checks if the path has already visited this vertice
        if (path.getVisited().contains(currentVertice)) {
            return null;
        }
        //if it did not visit, then add to the path
        path.getVisited().add(currentVertice);

        //if the current vertice is the target vertice, we found our path
        if (currentVertice.equals(targetVertice)) {
            return path;
        }
        //check all edges in the vertice, and add the edge to the path edges and return the path until there are no edges left
        for (E e : currentVertice.getEdges()) {
            if (DFSRecursion(e.getTo(), targetVertice, path) != null) {
                path.getEdges().addFirst(e);
                return path;
            }
        }
        return null;

    }

    /**
     * Uses a breadth-first search algorithm to find a path from the start vertex to the target vertex in the graph
     * The path.totalWeight should indicate the number of edges in the result path
     * All vertices that are being visited by the search should also be registered in path.visited
     *
     * @param startId
     * @param targetId
     * @return the path from start to target
     * returns null if either start or target cannot be matched with a vertex in the graph
     * or no path can be found from start to target
     */
    public DGPath breadthFirstSearch(String startId, String targetId) {

        V start = this.getVertexById(startId);
        V target = this.getVertexById(targetId);
        if (start == null || target == null) return null;

        DGPath path = new DGPath();
        path.start = start;
        path.visited.add(start);

        // easy target
        if (start == target) return path;
        //create a queue for bfs
        LinkedList<V> queue = new LinkedList<>();
        //create a map to keep track of visited nodes
        Map<V, V> visited = new HashMap<>();
        //put start node in queue
        queue.offer(start);
        //put start node in the visited nodes
        visited.put(start, null);

        //keep searching until the queue is empty
        while (!queue.isEmpty()) {
            //get the current node by removing the first node from the queue
            V current = queue.poll();
            //get all the adjacent vertices of dequeued vertex
            for (E e : current.getEdges()) {
                V neighbour = e.getTo();
                path.visited.add(neighbour);

                //if the node hits target,
                if (neighbour == target) {
                    path.getEdges().addFirst(e);
                    while (current != null) {
                        path.getEdges().addFirst(e);
                        //   path.getVisited().add()
                        current = visited.get(current);
                    }
                    return path;
                    //push neighbours into queue if they havent been explored already
                } else if (!visited.containsKey(neighbour)) {
                    visited.put(neighbour, current);
                    queue.offer(neighbour);
                }
            }
        }
        //no path found
        return null;
    }

    // helper class to register the state of a vertex in dijkstra shortest path algorithm
    // your may change this class or delete it altogether follow a different approach in your implementation
    private class DSPNode implements Comparable<DSPNode> {
        public V vertex;                // the graph vertex that is concerned with this DSPNode
        public E fromEdge = null;        // the edge from the predecessor's vertex to this node's vertex
        public boolean marked = false;  // indicates DSP processing has been marked complete
        public double weightSumTo = Double.MAX_VALUE;   // sum of weights of current shortest path to this node's vertex

        public DSPNode(V vertex) {
            this.vertex = vertex;
        }

        // comparable interface helps to find a node with the shortest current path, sofar
        @Override
        public int compareTo(DSPNode dspv) {
            return Double.compare(this.weightSumTo, dspv.weightSumTo);
        }
    }


    /**
     * Calculates the edge-weighted shortest path from start to target
     * Uses a minimum distance heuristic from any vertex to the target
     * in order to reduce the number of visited vertices during the search
     *
     * @param startId
     * @param targetId
     * @param weightMapper provides a function, by which the weight of an edge can be retrieved or calculated
     * @return the shortest path from start to target
     * returns null if either start or target cannot be matched with a vertex in the graph
     * or no path can be found from start to target
     */
    public DGPath dijkstraShortestPath(String startId, String targetId,
                                       Function<E, Double> weightMapper) {

        V start = this.getVertexById(startId);
        V target = this.getVertexById(targetId);

        if (start == null || target == null) return null;

        // initialise the result path of the search
        DGPath path = new DGPath();
        path.start = start;
        path.visited.add(start);

        // easy target
        if (start == target) return path;
        // keep track of the DSP status of all visited nodes
        // you may choose a different approach of tracking progress of the algorith, if you wish

        Map<V, DSPNode> progressData = new HashMap<>();

        // initialise the progress of the start node
        DSPNode nextDspNode = new DSPNode(start);
        //mark starting node 0.0 as per rules
        nextDspNode.weightSumTo = 0.0;
        //add it to visited nodes & mark
        progressData.put(start, nextDspNode);

        nextDspNode.marked = true;
        //while we have a next node
        while (nextDspNode != null) {
            //for all the adjacent edges
            for (E edge : nextDspNode.vertex.getEdges()) {
                // get the neighbour vertice and create a new node from it
                // then add it to visited vertices
                V neighbourV = edge.getTo();
                DSPNode neighbourN;
                path.getEdges().addFirst(edge);

                //if the visited nodes contains this new node , update
                if (progressData.containsKey(neighbourV)) {
                    neighbourN = progressData.get(neighbourV);
                    //has not been visited so we create the new node
                } else {
                    neighbourN = new DSPNode(neighbourV);
                }
                //add to path
                path.visited.add(neighbourV);
                //apply weight to the edge
                double weight = weightMapper.apply(edge);
                double costNextNode = nextDspNode.weightSumTo + weight;
                //check if the cost of the next is less
                //assign its weight and its neighbour
                // add it to visited nodes.
                if (costNextNode < neighbourN.weightSumTo) {
                    neighbourN.weightSumTo = costNextNode;
                    neighbourN.fromEdge = edge;
                    progressData.put(neighbourV, neighbourN);
                }
            }

            DSPNode current = progressData.get(nextDspNode.vertex);
            current.marked = true;
            if (nextDspNode.vertex == target) break; // found our path
            //get the next node from the visited nodes that have not been marked yet
            nextDspNode = progressData.values().stream().filter(dspNode -> !dspNode.marked).
                    findFirst().orElse(null);
            // no path found, graph was not connected ???
            //get the totalweight of all the nodes by summing all the node weight
            path.totalWeight = progressData.values().stream().mapToDouble(value -> value.weightSumTo).sum();
            return path;

        }
        return null;
    }
    // helper class to register the state of a vertex in A* shortest path algorithm
//    private class ASNode extends DSPNode {
//        // TODO add and handle information for the minimumWeightEstimator
//        double minimumWeight = 0.0;
//        V locationX;
//        V locationY;
//        // TODO enhance this constructor as required
//        private ASNode(V vertex) {
//            super(vertex);
//        }
//        @Override
//        public int compareTo ASNode(ASNode asNode){
//
//
//            return Double.compare(this.minimumWeight, asNode.weightSumTo);
//
//        }
//        // TODO override the compareTo
//    }


    /**
     * Calculates the edge-weighted shortest path from start to target
     * Uses a minimum distance heuristic from any vertex to the target
     * in order to reduce the number of visited vertices during the search
     *
     * @param startId
     * @param targetId
     * @param weightMapper           provides a function, by which the weight of an edge can be retrieved or calculated
     * @param minimumWeightEstimator provides a function, by which a lower bound of the cumulative weight
     *                               between two vertices can be calculated.
     * @return the shortest path from start to target
     * returns null if either start or target cannot be matched with a vertex in the graph
     * or no path can be found from start to target
     */
    public DGPath aStarShortestPath(String startId, String targetId,
                                    Function<E, Double> weightMapper,
                                    BiFunction<V, V, Double> minimumWeightEstimator) {

        V start = this.getVertexById(startId);
        V target = this.getVertexById(targetId);
        if (start == null || target == null) return null;

        DGPath path = new DGPath();
        path.start = start;
        path.visited.add(start);

        // easy target
        if (start == target) return path;

        // TODO apply the A* algorithm to find shortest path from start to target.
        //  take dijkstra's solution as the starting point and enhance with heuristic functionality
        //  register all visited vertices while going, for statistical purposes


        // TODO END
        // no path found, graph was not connected ???
        return null;
    }

    /**
     * Calculates the edge-weighted shortest path from start to target
     *
     * @param startId
     * @param targetId
     * @param weightMapper provides a function by which the weight of an edge can be retrieved or calculated
     * @return the shortest path from start to target
     * returns null if either start or target cannot be matched with a vertex in the graph
     * or no path can be found from start to target
     */
    public DGPath dijkstraShortestPathByAStar(String startId, String targetId,
                                              Function<E, Double> weightMapper) {
        return aStarShortestPath(startId, targetId,
                weightMapper,
                // TODO provide a minimumWeightEstimator that makes A* run like regular Dijkstra
                getMinimumWeightEstimator()
        );
    }

    private BiFunction<V, V, Double> getMinimumWeightEstimator() {
        return null;
    }

    @Override
    public String toString() {
        return this.getVertices().stream()
                .map(Object::toString)
                .collect(Collectors.joining(",\n  ", "{ ", "\n}"));
    }
}

