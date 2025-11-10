#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define INFINITY 1000000000
#define NEARBY_MODE 0
#define FASTEST_MODE 1

// bit-flagger
#define PLACENAME_FLAG 0x01
#define BENSINSTATION_FLAG 0x02
#define CHARINGSTATION_FLAG 0x04
#define DININGPLACE_FLAG 0x08
#define DRINKPLACE_FLAG 0x10
#define HOTEL_FLAG 0x20

typedef struct Vertex {
  struct Edge *edge;
  float latitude, longitude;
  void *data;
} vertex_t;

typedef struct PreprocessedLandmarks {
  int64_t num_landmarks;
  int64_t num_vertices;
  uint64_t *distances;
} preprocessed_landmarks_t;

typedef struct Edge {
  struct Edge *next;
  vertex_t *to;
  uint32_t drive_time;
  uint32_t length;
  uint8_t speed;
} edge_t;

typedef struct Graph {
  int numVertices, numEdges;
  vertex_t *vertex;
} graph_t;

typedef struct Predecessor {
  struct Vertex *predecessor;
  uint32_t distance;
  uint32_t total_drive_time;
} predecessor_t;

typedef struct PriorityQueueNode {
  uint32_t weight;
  void *data;
} pq_node_t;

typedef struct PriorityQueue {
  pq_node_t **data;
  size_t size;
  size_t capacity;
} pq_t;

typedef struct Token {
  union {
    double number;
    char *string;
  } value;
  bool is_string;
} token_t;

typedef struct FileExtractResult {
  int numberOfReadData;
  void *data;
} fileExtractResult_t;

typedef struct POI {
  char *name;
  uint8_t flags;
  uint32_t node_id;
} poi_t;

bool is_flag_set(uint8_t flags, uint8_t flag_to_check) {
  return (flags & flag_to_check) != 0;
}

void toggle_flag(uint8_t *flags, uint8_t flag_to_toggle) {
  *flags ^= flag_to_toggle;
}

// Heap og prioritetskø implementasjon for dijkstra
void pq_swap_node(pq_node_t **a, pq_node_t **b) {
  pq_node_t *t = *a;
  *a = *b;
  *b = t;
}
pq_node_t *create_pq_node(uint32_t weight, void *data) {
  pq_node_t *node = (pq_node_t *) malloc(sizeof(pq_node_t));
  node->weight = weight;
  node->data = data;
  return node;
}

pq_t *create_pq(size_t capacity) {
  pq_t *pq = (pq_t *) malloc(sizeof(pq_t));
  pq->data = (pq_node_t **) malloc(capacity * sizeof(pq_node_t *));
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
        (pq_node_t **) realloc(pq->data, pq->capacity * sizeof(pq_node_t *));
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
  if (graph == NULL)
    return false;
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

void addEdge(graph_t *graph,
             const int from,
             const int to,
             uint32_t drive_time,
             uint32_t length,
             uint8_t speed) {
  if (graph == NULL)
    return;
  if (from < 0 || from >= graph->numVertices)
    return;
  if (to < 0 || to >= graph->numVertices)
    return;

  edge_t *newEdge = malloc(sizeof(edge_t));
  if (newEdge == NULL)
    return;

  newEdge->to = &graph->vertex[to];
  newEdge->drive_time = drive_time;
  newEdge->length = length;
  newEdge->speed = speed;
  newEdge->next = graph->vertex[from].edge;
  graph->vertex[from].edge = newEdge;
  graph->numEdges++;
}

void remove_edge(graph_t *graph, const int from, const int to) {
  if (graph == NULL)
    return;
  if (from < 0 || from >= graph->numVertices)
    return;
  if (to < 0 || to >= graph->numVertices)
    return;

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

fileExtractResult_t extract_tokens_from_line(const char *line, int max_count) {
  fileExtractResult_t result;
  result.numberOfReadData = 0;
  token_t *data = malloc(max_count * sizeof(token_t));
  const char *ptr = line;

  while (*ptr != '\0' && result.numberOfReadData < max_count) {
    while (*ptr == ' ' || *ptr == '\t')
      ptr++;
    if (*ptr == '\n' || *ptr == '\0')
      break;

    if (*ptr == '"') {
      ptr++;
      const char *start = ptr;
      while (*ptr != '"' && *ptr != '\0')
        ptr++;
      size_t len = ptr - start;
      data[result.numberOfReadData].value.string = malloc(len + 1);
      memcpy(data[result.numberOfReadData].value.string, start, len);
      data[result.numberOfReadData].value.string[len] = '\0';
      data[result.numberOfReadData].is_string = true;
      result.numberOfReadData++;
      if (*ptr == '"')
        ptr++;
    } else {
      char *endPtr;
      double num = strtod(ptr, &endPtr);
      if (ptr != endPtr) {
        data[result.numberOfReadData].value.number = num;
        data[result.numberOfReadData].is_string = false;
        result.numberOfReadData++;
        ptr = endPtr;
      } else {
        break;
      }
    }
  }
  result.data = data;
  return result;
}

void free_tokens(token_t *data, int count) {
  for (int i = 0; i < count; i++) {
    if (data[i].is_string) {
      free(data[i].value.string);
    }
  }
  free(data);
}

graph_t *initialize_graph(const char *nodes_filename,
                          const char *edges_filename,
                          const char *poi_filename) {
  FILE *edgesFile = fopen(edges_filename, "r");
  FILE *nodesFile = fopen(nodes_filename, "r");
  FILE *poiFile = fopen(poi_filename, "r");

  if (edgesFile == NULL || nodesFile == NULL || poiFile == NULL) {
    perror("Failed to open file");
    if (edgesFile)
      fclose(edgesFile);
    if (nodesFile)
      fclose(nodesFile);
    if (poiFile)
      fclose(poiFile);
    return NULL;
  }
  char *edgesLine = malloc(256 * sizeof(char));
  char *nodesLine = malloc(256 * sizeof(char));
  char *poiLine = malloc(256 * sizeof(char));
  fgets(edgesLine, 256, edgesFile);
  fgets(nodesLine, 256, nodesFile);
  fgets(poiLine, 256, poiFile);
  // edges
  const fileExtractResult_t firstLineDataEdges = extract_tokens_from_line(edgesLine, 5);
  printf("numEges: %d\n", (int) ((token_t *) firstLineDataEdges.data)[0].value.number);

  // nodes
  const fileExtractResult_t firstLineDataNodes = extract_tokens_from_line(nodesLine, 5);
  printf("numNodes: %d\n", (int) ((token_t *) firstLineDataNodes.data)[0].value.number);

  // poi
  const fileExtractResult_t firstLineDataPOI = extract_tokens_from_line(poiLine, 3);
  printf("numPOI: %d\n", (int) ((token_t *) firstLineDataPOI.data)[0].value.number);

  const int numEdges = (int) ((token_t *) firstLineDataEdges.data)[0].value.number;
  const int numNodes = (int) ((token_t *) firstLineDataNodes.data)[0].value.number;
  const int numPOI = (int) ((token_t *) firstLineDataPOI.data)[0].value.number;
  graph_t *graph = createGraph(numNodes);

  free_tokens(firstLineDataEdges.data, firstLineDataEdges.numberOfReadData);
  free_tokens(firstLineDataNodes.data, firstLineDataNodes.numberOfReadData);
  free_tokens(firstLineDataPOI.data, firstLineDataPOI.numberOfReadData);
  fclose(poiFile);

  printf("Processing nodes\n");
  while (fgets(nodesLine, 256, nodesFile)) {
    const fileExtractResult_t nodeData = extract_tokens_from_line(nodesLine, 3);
    if (nodeData.numberOfReadData == 3) {
      token_t *data = nodeData.data;
      const int nodeIndex = (int) data[0].value.number;
      const float lat = (float) data[1].value.number;
      const float lon = (float) data[2].value.number;
      free_tokens(data, nodeData.numberOfReadData);

      graph->vertex[nodeIndex].latitude = lat;
      graph->vertex[nodeIndex].longitude = lon;
    }
  }

  printf("Processing edges\n");
  while (fgets(edgesLine, 256, edgesFile)) {
    const fileExtractResult_t edgeData = extract_tokens_from_line(edgesLine, 5);
    if (edgeData.numberOfReadData == 5) {
      const token_t *data = edgeData.data;
      addEdge(graph,
              (int) data[0].value.number,
              (int) data[1].value.number,
              (uint32_t) data[2].value.number,
              (uint32_t) data[3].value.number,
              (uint8_t) data[4].value.number);
      free(edgeData.data);
    }
  }

  // validere antall kanter
  if (graph->numEdges != numEdges) {
    printf("Warning: Expected %d edges, but found %d edges in the file.\n",
           numEdges,
           graph->numEdges);
  }
  fclose(edgesFile);
  fclose(nodesFile);
  free(nodesLine);
  free(edgesLine);
  return graph;
}

predecessor_t *new_predecessor() {
  predecessor_t *predecessor = malloc(sizeof(predecessor_t));
  predecessor->predecessor = NULL;
  predecessor->distance = INFINITY;
  predecessor->total_drive_time = INFINITY;
  return predecessor;
}

void init_predecessors(const graph_t *graph, const vertex_t *start) {
  for (int i = 0; i < graph->numVertices; ++i) {
    if (graph->vertex[i].data == NULL) {
      graph->vertex[i].data = new_predecessor();
    } else {
      predecessor_t *p = graph->vertex[i].data;
      p->predecessor = NULL;
      p->distance = INFINITY;
      p->total_drive_time = INFINITY;
    }
  }
  ((predecessor_t *) start->data)->distance = 0;
  ((predecessor_t *) start->data)->total_drive_time = 0;
}

void dijkstra(const graph_t *graph, const vertex_t *start, const vertex_t *end) {
  init_predecessors(graph, start);
  pq_t *pq = create_pq(graph->numVertices);
  pq_push(pq, create_pq_node(0, (void *) start));
  bool *visited = calloc(graph->numVertices, sizeof(bool));
  while (pq->size) {
    pq_node_t *current = pq_pop(pq);
    vertex_t *current_vertex = current->data;
    free(current);
    int vertex_index = (int) (current_vertex - graph->vertex);
    if (visited[vertex_index]) {
      continue;
    }
    visited[vertex_index] = true;

    if (current_vertex == end) {
      free(visited);
      free_pq(pq);
      return;
    }
    const predecessor_t *current_pred = current_vertex->data;
    edge_t *current_edge = current_vertex->edge;
    while (current_edge != NULL) {
      predecessor_t *neighbor_pred = current_edge->to->data;
      const uint32_t new_drivetime = current_pred->total_drive_time + current_edge->drive_time;
      if (new_drivetime < neighbor_pred->total_drive_time) {
        neighbor_pred->total_drive_time = new_drivetime;
        neighbor_pred->predecessor = current_vertex;
        pq_push(pq, create_pq_node(new_drivetime, (void *) current_edge->to));
      }
      current_edge = current_edge->next;
    }
  }

  free(visited);
  free_pq(pq);
}

uint64_t ALT_heuristic(const graph_t *graph,
                       vertex_t *current,
                       vertex_t *end,
                       preprocessed_landmarks_t **landmarks,
                       int *landmark_indices,
                       int num_landmarks) {
  int64_t best = 0;
  int current_index = (int) (current - graph->vertex);

  int end_index = (int) (end - graph->vertex);

  for (int i = 0; i < num_landmarks; i++) {
    int64_t fwd = (int64_t) landmarks[0]->distances[i * graph->numVertices + end_index]
        - (int64_t) landmarks[0]->distances[i * graph->numVertices + current_index];
    int64_t bwd = (int64_t) landmarks[1]->distances[i * graph->numVertices + current_index]
        - (int64_t) landmarks[1]->distances[i * graph->numVertices + end_index];
    if (fwd > best) best = fwd;
    if (bwd > best) best = bwd;
  }
  if (best < 0) best = 0;
  return (uint64_t) best;
}


void ALT(graph_t *graph,
         vertex_t *start,
         vertex_t *end,
         preprocessed_landmarks_t **landmarks,
         int *landmark_indices,
         int num_landmarks) {
  init_predecessors(graph, start);
  pq_t *pq = create_pq(graph->numVertices);
  pq_push(pq, create_pq_node(0, (void *) start));
  bool *visited = calloc(graph->numVertices, sizeof(bool));
  while (pq->size) {
    pq_node_t *current = pq_pop(pq);
    vertex_t *current_vertex = current->data;
    free(current);
    int vertex_index = (int) (current_vertex - graph->vertex);
    if (visited[vertex_index]) {
      continue;
    }
    visited[vertex_index] = true;

    if (current_vertex == end) {
      free(visited);
      free_pq(pq);
      return;
    }
    const predecessor_t *current_pred = current_vertex->data;
    edge_t *current_edge = current_vertex->edge;
    while (current_edge != NULL) {
      predecessor_t *neighbor_pred = current_edge->to->data;

      uint64_t heuristic = ALT_heuristic(graph,
                                         current_edge->to,
                                         end,
                                         landmarks,
                                         landmark_indices,
                                         num_landmarks);
      const uint32_t new_drivetime = current_pred->total_drive_time + current_edge->drive_time;
      const uint32_t estimated_drive_time = new_drivetime + heuristic;
      if (new_drivetime < neighbor_pred->total_drive_time) {
        neighbor_pred->total_drive_time = new_drivetime;
        neighbor_pred->predecessor = current_vertex;
        pq_push(pq, create_pq_node(estimated_drive_time, current_edge->to));
      }
      current_edge = current_edge->next;
    }
  }

  free(visited);
  free_pq(pq);
}


graph_t *reverse_graph(graph_t *graph) {
  graph_t *reversed_graph = createGraph(graph->numVertices);
  for (int i = 0; i < graph->numVertices; i++) {
    edge_t *current_edge = graph->vertex[i].edge;
    while (current_edge != NULL) {
      int to_index = (int) (current_edge->to - graph->vertex);
      addEdge(reversed_graph,
              to_index,
              i,
              current_edge->drive_time,
              current_edge->length,
              current_edge->speed);
      current_edge = current_edge->next;
    }
  }
  return reversed_graph;
}

preprocessed_landmarks_t **preprocess_graph(const graph_t *graph,
                                            int *landmarks_indices,
                                            int num_landmarks) {
  FILE *file = fopen("distances.yz", "rb");

  if (file != NULL) {
    printf("Landmark file exists, loading landmarks.\n");
    char *line = malloc(256 * sizeof(char));
    fgets(line, 256, file);
    fileExtractResult_t firstLineData = extract_tokens_from_line(line, 2);
    const int numLandmarks = (int) ((token_t *) firstLineData.data)[0].value.number;
    const int numVertices = (int) ((token_t *) firstLineData.data)[1].value.number;

    if (numLandmarks != num_landmarks || numVertices != graph->numVertices) {
      printf("Landmark data in file does not match expected values.\n");
      free_tokens(firstLineData.data, firstLineData.numberOfReadData);
      fclose(file);
      free(line);
      return NULL;
    }
    fgets(line, 256, file);
    fileExtractResult_t indicesData = extract_tokens_from_line(line, numLandmarks);
    bool indices_match = true;
    for (int i = 0; i < numLandmarks; i++) {
      if ((int) ((token_t *) indicesData.data)[i].value.number != landmarks_indices[i]) {
        indices_match = false;
        break;
      }
    }

    if (!indices_match) {
      printf("Landmark indices in file do not match.\n");
      free_tokens(firstLineData.data, firstLineData.numberOfReadData);
      free_tokens(indicesData.data, indicesData.numberOfReadData);
      fclose(file);
      free(line);
      return NULL;
    }

    free_tokens(indicesData.data, indicesData.numberOfReadData);

    preprocessed_landmarks_t **landmarks = malloc(2 * sizeof(preprocessed_landmarks_t *));

    // distanser for fremover
    landmarks[0] = malloc(sizeof(preprocessed_landmarks_t));
    landmarks[0]->num_landmarks = numLandmarks;
    landmarks[0]->num_vertices = numVertices;
    landmarks[0]->distances = malloc(numLandmarks * numVertices * sizeof(uint64_t));
    fread(landmarks[0]->distances, sizeof(uint64_t), numLandmarks * numVertices, file);

    // distanser for bakover
    landmarks[1] = malloc(sizeof(preprocessed_landmarks_t));
    landmarks[1]->num_landmarks = numLandmarks;
    landmarks[1]->num_vertices = numVertices;
    landmarks[1]->distances = malloc(numLandmarks * numVertices * sizeof(uint64_t));
    fread(landmarks[1]->distances, sizeof(uint64_t), numLandmarks * numVertices, file);
    free_tokens(firstLineData.data, firstLineData.numberOfReadData);
    fclose(file);
    free(line);
    return landmarks;
  }

  printf("Landmark file does not exist, preprocessing landmarks.\n");
  file = fopen("distances.yz", "wb");
  if (file == NULL) {
    perror("Failed to open file for writing");
    return NULL;
  }
  fprintf(file, "%d %d\n", num_landmarks, graph->numVertices);
  for (int i = 0; i < num_landmarks; i++) {
    fprintf(file, "%d ", landmarks_indices[i]);
  }
  fprintf(file, "\n");
  preprocessed_landmarks_t **landmarks = malloc(2 * sizeof(preprocessed_landmarks_t *));

  // distanser for fremover
  landmarks[0] = malloc(sizeof(preprocessed_landmarks_t));
  landmarks[0]->num_landmarks = num_landmarks;
  landmarks[0]->num_vertices = graph->numVertices;
  landmarks[0]->distances = malloc(num_landmarks * graph->numVertices * sizeof(uint64_t));
  // regne ut distanser mellom landmarks og alle noder
  for (int i = 0; i < num_landmarks; i++) {
    dijkstra(graph, &graph->vertex[landmarks_indices[i]], NULL);
    for (int j = 0; j < graph->numVertices; j++) {
      predecessor_t *pred = graph->vertex[j].data;
      landmarks[0]->distances[i * graph->numVertices + j] = pred->total_drive_time;
    }
  }

  // distanser for bakover
  landmarks[1] = malloc(sizeof(preprocessed_landmarks_t));
  landmarks[1]->num_landmarks = num_landmarks;
  landmarks[1]->num_vertices = graph->numVertices;
  landmarks[1]->distances = malloc(num_landmarks * graph->numVertices * sizeof(uint64_t));
  graph_t *reversed_graph = reverse_graph((graph_t *) graph);
  for (int i = 0; i < num_landmarks; i++) {
    dijkstra(reversed_graph, &reversed_graph->vertex[landmarks_indices[i]], NULL);
    for (int j = 0; j < reversed_graph->numVertices; j++) {
      predecessor_t *pred = reversed_graph->vertex[j].data;
      landmarks[1]->distances[i * reversed_graph->numVertices + j] = pred->total_drive_time;
    }
  }
  freeGraph(reversed_graph);
  // skrive til fil
  fwrite(landmarks[0]->distances, 1, num_landmarks * graph->numVertices * sizeof(uint64_t), file);
  fwrite(landmarks[1]->distances, 1, num_landmarks * graph->numVertices * sizeof(uint64_t), file);
  fclose(file);
  return landmarks;
}

void free_preprocessed_landmarks(preprocessed_landmarks_t **landmarks, int num_landmarks) {
  if (landmarks == NULL)
    return;
  for (int i = 0; i < num_landmarks; i++) {
    if (landmarks[i] != NULL) {
      if (landmarks[i]->distances != NULL) {
        free(landmarks[i]->distances);
      }
      free(landmarks[i]);
    }
  }
  free(landmarks);
}

void print_graph(const graph_t graph) {
  printf("%-10s %-10s %-10s\n", "Node", "Forgj", "Dist");
  printf("--------------------------------------------------------\n");
  for (int i = 0; i < graph.numVertices; i++) {
    if (graph.vertex[i].data != NULL) {
      const predecessor_t *pred = (predecessor_t *) graph.vertex[i].data;
      if (pred->predecessor != NULL) {
        const int predIndex = (int) (pred->predecessor - graph.vertex);
        printf("%-10d %-10d %-10d\n", i, predIndex, pred->distance);
      } else {
        printf("%-10d %-10s %-10s\n", i, " ", "∞");
      }
    }
  }
}

void print_path(const graph_t *graph, const vertex_t *start, const vertex_t *end) {
  FILE *file = fopen("path.txt", "w");
  if (graph == NULL || start == NULL || end == NULL) return;

  predecessor_t *current_pred = end->data;
  if (current_pred->total_drive_time == INFINITY) {
    printf("No path from start to end.\n");
    if (file) fclose(file);
    return;
  }

  printf("Path from node %lld to node %lld (Drive time: %u):\n",
         (long long) (start - graph->vertex),
         (long long) (end - graph->vertex),
         current_pred->total_drive_time);

  vertex_t *current_vertex = (vertex_t *) end;
  while (current_vertex != NULL) {
    fprintf(file, "%f,%f\n", current_vertex->latitude, current_vertex->longitude);
    predecessor_t *pred = current_vertex->data;
    current_vertex = pred->predecessor;
  }
  fclose(file);
}

int main(int argc, char *argv[]) {
  int landmarks_indices[] = {200, 1000, 2000, 3000, 4000};
  graph_t *graph = initialize_graph("noder.txt", "kanter.txt", "interessepkt.txt");
  preprocessed_landmarks_t **landmarks = preprocess_graph(graph, landmarks_indices, 5);
  vertex_t *start = &graph->vertex[2918249];
  vertex_t *end = &graph->vertex[7916010];
  ALT(graph, start, end, landmarks, landmarks_indices, 5);
  print_path(graph, start, end);
  freeGraph(graph);
  free_preprocessed_landmarks(landmarks, sizeof(landmarks_indices) / sizeof(landmarks_indices[0]));
  return EXIT_SUCCESS;
}
