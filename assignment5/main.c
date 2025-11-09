#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define infinity 1000000000

typedef struct Vertex {
  struct Edge *edge;
  void *data;
} vertex_t;

typedef struct Edge {
  struct Edge *next;
  vertex_t *to;
} edge_t;

typedef struct Graph {
  int numVertices, numEdges;
  vertex_t *vertex;
} graph_t;

typedef struct FileExtractResult {
  int numberOfReadData;
  int data[2];
} fileExtractResult_t;

typedef struct Predecessor {
  struct Vertex* predecessor;
  int distance;
} predecessor_t;

typedef struct BFSQueue {
  vertex_t **elements;
  int front;
  int back;
  int capacity;
} bfs_queue_t;

typedef struct TopologicalListEdge {
  int found;
  struct TopologicalListEdge* next;
} topological_list_edge_t;



graph_t *createGraph(const int numVertices) {
  graph_t *graph = malloc(sizeof(graph_t));
  graph->numVertices = numVertices;
  graph->numEdges = 0;
  graph->vertex = malloc(numVertices * sizeof(vertex_t));
  for (int i = 0; i < numVertices; i++) {
    graph->vertex[i].edge = NULL;
  }
  return graph;
}

bool freeGraph(graph_t *graph) {
  if (graph == NULL) return false;
  if (graph->vertex == NULL) {
    free(graph);
    return true;
  }

  for (int i = 0; i < graph->numVertices; i++) {
    edge_t *current = graph->vertex[i].edge;
    while (current != NULL) {
      edge_t *temp = current;
      current = current->next;
      free(temp);
    }
  }

  free(graph->vertex);
  free(graph);
  return true;
}

void addEdge(graph_t *graph, const int from, const int to) {
  if (graph == NULL) return;
  if (from < 0 || from >= graph->numVertices) return;
  if (to < 0 || to >= graph->numVertices) return;

  edge_t *newEdge = malloc(sizeof(edge_t));
  if (newEdge == NULL) return;

  newEdge->to = &graph->vertex[to];
  newEdge->next = graph->vertex[from].edge;
  graph->vertex[from].edge = newEdge;
  graph->numEdges++;
}

void remove_edge(graph_t *graph, const int from, const int to) {
  if (graph == NULL) return;
  if (from < 0 || from >= graph->numVertices) return;
  if (to < 0 || to >= graph->numVertices) return;

  edge_t *current = graph->vertex[from].edge;
  edge_t *previous = NULL;

  while (current != NULL) {
    if (current->to == &graph->vertex[to]) {
      if (previous == NULL) {
        graph->vertex[from].edge = current->next;
      } else {
        previous->next = current->next;
      }
      free(current);
      graph->numEdges--;
      return;
    }

    previous = current;
    current = current->next;
  }
}

void print_graph(const graph_t graph) {
  printf("%-10s %-10s %-10s\n", "Node", "Forgj", "Dist");
  printf("--------------------------------------------------------\n");
  for (int i = 0; i < graph.numVertices; i++) {
    if (graph.vertex[i].data != NULL) {
      const predecessor_t *pred = (predecessor_t*) graph.vertex[i].data;
      if (pred->predecessor != NULL) {
        const int predIndex = (int) (pred->predecessor - graph.vertex);
        printf("%-10d %-10d %-10d\n", i, predIndex, pred->distance);
      } else {
        printf("%-10d %-10s %-10s\n", i, " ", "∞");
      }
    }
  }
}


fileExtractResult_t extract_data_from_line(const char *line) {
  fileExtractResult_t result;
  result.numberOfReadData = 0;
  const char *ptr = line;
  while (*ptr != '\0') {
    while (*ptr == ' ' || *ptr == '\t') ptr++;
    if (*ptr == '\n' || *ptr == '\0') break;

    char *endPtr;
    int value = (int) strtol(ptr, &endPtr, 10);
    if (ptr == endPtr) {
      break;
    }
    result.data[result.numberOfReadData++] = value;
    ptr = endPtr;
  }
  return result;
}


graph_t* create_graph_from_file(const char *filename) {
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    perror("Failed to open file");
    return NULL;
  }
  char *line = malloc(256 * sizeof(char));
  fgets(line, 256, file);
  const fileExtractResult_t firstLineData = extract_data_from_line(line);
  printf("number of vertices %d, number of edges %d\n", firstLineData.data[0], firstLineData.data[1]);
  const int numVertices = firstLineData.data[0];
  const int numEdges = firstLineData.data[1];
  graph_t *graph = createGraph(numVertices);

  while (fgets(line, 256, file)) {
    const fileExtractResult_t edgeData = extract_data_from_line(line);
    if (edgeData.numberOfReadData == 2) {
      addEdge(graph, edgeData.data[0], edgeData.data[1]);
    }
  }

  // validere antall kanter
  if (graph->numEdges != numEdges) {
    printf("Warning: Expected %d edges, but found %d edges in the file.\n", numEdges, graph->numEdges);
  }
  fclose(file);
  free(line);
  return graph;
}

predecessor_t* new_predecessor() {
  predecessor_t *predecessor = malloc(sizeof(predecessor_t));
  predecessor->predecessor = NULL;
  predecessor->distance = infinity;
  return predecessor;
}

void init_predecessors(const graph_t *graph, const vertex_t* vertex) {
  for (int i = 0; i < graph->numVertices; ++i) {
    graph->vertex[i].data = new_predecessor();
  }
  ((predecessor_t*)vertex->data)->distance = 0;
}

bfs_queue_t* create_bfs_queue(int capacity) {
  bfs_queue_t *queue = malloc(sizeof(bfs_queue_t));
  queue->elements = malloc(capacity * sizeof(vertex_t*));
  queue->front = 0;
  queue->back = 0;
  queue->capacity = capacity;
  return queue;
}

void bfs(const graph_t *graph, const vertex_t* start) {
  if (graph == NULL || start == NULL) return;
  init_predecessors(graph, start);
  bfs_queue_t *queue = create_bfs_queue(graph->numVertices);
  queue->elements[queue->back++] = (vertex_t*)start;
  while (queue->front < queue->back) {
    vertex_t *current = queue->elements[queue->front++];
    const predecessor_t *currentPred = (predecessor_t*) current->data;
    const edge_t *edge = current->edge;
    while (edge != NULL) {
      predecessor_t *neighborPred = (predecessor_t*) edge->to->data;
      if (neighborPred->distance == infinity) {
        neighborPred->distance = currentPred->distance + 1;
        neighborPred->predecessor = current;
        queue->elements[queue->back++] = edge->to;
      }
      edge = edge->next;
    }
  }
  free(queue->elements);
  free(queue);
}

vertex_t* df_topo(vertex_t* start_vertex, vertex_t* vertex_list) {
  topological_list_edge_t *current = start_vertex->data;
  if (current->found) return vertex_list;
  current->found = 1;

  for (const edge_t* edge = start_vertex->edge; edge; edge = edge->next) {
    vertex_list = df_topo(edge->to, vertex_list);
  }

  if (vertex_list != NULL) {
    current->next = vertex_list->data;
  } else {
    current->next = NULL;
  }

  return start_vertex;
}

vertex_t* topological_sort(const graph_t* graph) {
  if (graph == NULL) return NULL;

  for (int i = 0; i < graph->numVertices; ++i) {
    graph->vertex[i].data = calloc(1, sizeof(topological_list_edge_t));;
  }

  vertex_t* vertex_list = NULL;
  for (int i = 0; i < graph->numVertices; ++i) {
    vertex_list = df_topo(&graph->vertex[i], vertex_list);
  }

  return vertex_list;
}


void print_topological_list(const vertex_t* list, const graph_t* graph) {
  if (list == NULL || graph == NULL || graph->vertex == NULL) {
    printf("hva skjedde her!?.\n");
    return;
  }
  printf("Topological order: ");
  for (const vertex_t* v = list; v != NULL; ) {
    printf("%d ", (int) (v - graph->vertex));
    const topological_list_edge_t* node = v->data;
    if (node->next == NULL) break;
    const vertex_t* next = NULL;
    for (int i = 0; i < graph->numVertices; i++) {
      if (graph->vertex[i].data == node->next) {
        next = &graph->vertex[i];
        break;
      }
    }
    v = next;
  }
  printf("\n");

}


void print_separator(char* title) {
  printf("\n\n================= %s =================\n", title);
}

int main() {
  // ø5g1
  print_separator("ø5g1");
  graph_t *graph = create_graph_from_file("ø5g1.txt");
  bfs(graph, &graph->vertex[0]);
  print_graph(*graph);
  freeGraph(graph);

  // ø5g2
  print_separator("ø5g2");
  graph = create_graph_from_file("ø5g2.txt");
  bfs(graph, &graph->vertex[1]);
  print_graph(*graph);
  freeGraph(graph);

  // ø5g3
  print_separator("ø5g3");
  graph = create_graph_from_file("ø5g3.txt");
  bfs(graph, &graph->vertex[2]);
  print_graph(*graph);
  freeGraph(graph);

  // ø5g5
  print_separator("ø5g5");
  graph = create_graph_from_file("ø5g5.txt");
  bfs(graph, &graph->vertex[3]);
  print_graph(*graph);
  const vertex_t* vertex_list = topological_sort(graph);
  print_topological_list(vertex_list, graph);
  freeGraph(graph);

  // ø5g7
  print_separator("ø5g7");
  graph = create_graph_from_file("ø5g7.txt");
  bfs(graph, &graph->vertex[4]);
  print_graph(*graph);
  const vertex_t* topo_list = topological_sort(graph);
  print_topological_list(topo_list, graph);
  freeGraph(graph);


  // ø5Skandinavia
  print_separator("ø5Skandinavia");
  graph = create_graph_from_file("ø5Skandinavia.txt");
  bfs(graph, &graph->vertex[37774]);
  print_graph(*graph);
  freeGraph(graph);

  return 0;
}
