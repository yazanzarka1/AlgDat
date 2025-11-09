#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#define INFINITY 1000000000

typedef struct PriorityQueueNode {
  uint32_t weight;
  struct PriorityQueueNode *left, *right;
} pq_node_t;

typedef struct PriorityQueue {
  pq_node_t **data;
  size_t size;
  size_t capacity;
} pq_t;

// Heap og prioritetskø implementasjon for dijkstra
void pq_swap_node(pq_node_t **a, pq_node_t **b) {
  pq_node_t *t = *a;
  *a = *b;
  *b = t;
}
pq_node_t *create_pq_node(int16_t byte, uint32_t weight, pq_node_t *left,
                          pq_node_t *right) {
  pq_node_t *node = (pq_node_t *)malloc(sizeof(pq_node_t));
  node->weight = weight;
  node->left = left;
  node->right = right;
  return node;
}

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

pq_t *create_pq(size_t capacity) {
  pq_t *pq = (pq_t *)malloc(sizeof(pq_t));
  pq->data = (pq_node_t **)malloc(capacity * sizeof(pq_node_t *));
  pq->size = 0;
  pq->capacity = capacity;
  return pq;
}

void heapify(pq_t *pq, size_t index) {
  while (index > 0) {
    size_t parent = (index - 1) / 2;
    if (pq->data[parent]->weight <= pq->data[index]->weight) {
      break;
    }
    pq_swap_node(&pq->data[parent], &pq->data[index]);
    index = parent;
  }
}

void free_pq(pq_t *pq) {
  for (size_t i = 0; i < pq->size; i++) {
    free(pq->data[i]);
  }
  free(pq->data);
  free(pq);
}

void pq_push(pq_t *pq, pq_node_t *node) {
  if (pq->size >= pq->capacity) {
    pq->capacity *= 2;
    pq_node_t **temp = pq->data;
    pq->data =
        (pq_node_t **)realloc(pq->data, pq->capacity * sizeof(pq_node_t *));
    if (pq->data == NULL) {
      free(temp);
      fprintf(stderr, "Memory reallocation failed\n");
      exit(EXIT_FAILURE);
    }
  }
  pq->data[pq->size] = node;
  heapify(pq, pq->size++);
}

pq_node_t *pq_pop(pq_t *pq) {
  if (pq->size == 0) {
    return NULL;
  }
  pq_node_t *min_node = pq->data[0];
  pq->data[0] = pq->data[--pq->size];
  size_t index = 0;
  while (1) {
    size_t left = 2 * index + 1;
    size_t right = 2 * index + 2;
    size_t smallest = index;
    if (left < pq->size && pq->data[left]->weight < pq->data[smallest]->weight) {
      smallest = left;
    }
    if (right < pq->size && pq->data[right]->weight < pq->data[smallest]->weight) {
      smallest = right;
    }
    if (smallest == index) {
      break;
    }
    pq_swap_node(&pq->data[index], &pq->data[smallest]);
    index = smallest;
  }
  return min_node;
}

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
  predecessor->distance = INFINITY;
  return predecessor;
}

void init_predecessors(const graph_t *graph, const vertex_t* vertex) {
  for (int i = 0; i < graph->numVertices; ++i) {
    graph->vertex[i].data = new_predecessor();
  }
  ((predecessor_t*)vertex->data)->distance = 0;
}

int main(int argc, char *argv[]) { return EXIT_SUCCESS; }
