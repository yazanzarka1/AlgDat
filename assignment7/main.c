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

typedef struct Vertex
{
  struct Edge *edge;
  float latitude, longitude;
  void *data;
} vertex_t;

typedef struct Edge
{
  struct Edge *next;
  vertex_t *to;
  uint32_t drive_time;
  uint32_t length;
  uint8_t speed;
} edge_t;

typedef struct Graph
{
  int numVertices, numEdges;
  vertex_t *vertex;
} graph_t;

typedef struct Predecessor
{
  struct Vertex *predecessor;
  uint32_t distance;
  uint32_t total_drive_time;
  uint32_t total_length;
} predecessor_t;

typedef struct PriorityQueueNode
{
  uint32_t weight;
  void *data;
} pq_node_t;

typedef struct PriorityQueue
{
  pq_node_t **data;
  size_t size;
  size_t capacity;
} pq_t;

typedef struct Token
{
  union
  {
    double number;
    char *string;
  } value;
  bool is_string;
} token_t;

typedef struct FileExtractResult
{
  int numberOfReadData;
  void *data;
} fileExtractResult_t;

bool is_flag_set(uint8_t flags, uint8_t flag_to_check)
{
  return (flags & flag_to_check) != 0;
}

void toggle_flag(uint8_t *flags, uint8_t flag_to_toggle)
{
  *flags ^= flag_to_toggle;
}

// Heap og prioritetskø implementasjon for dijkstra
void pq_swap_node(pq_node_t **a, pq_node_t **b)
{
  pq_node_t *t = *a;
  *a = *b;
  *b = t;
}
pq_node_t *create_pq_node(uint32_t weight, void *data)
{
  pq_node_t *node = (pq_node_t *)malloc(sizeof(pq_node_t));
  node->weight = weight;
  node->data = data;
  return node;
}

pq_t *create_pq(size_t capacity)
{
  pq_t *pq = (pq_t *)malloc(sizeof(pq_t));
  pq->data = (pq_node_t **)malloc(capacity * sizeof(pq_node_t *));
  pq->size = 0;
  pq->capacity = capacity;
  return pq;
}

void heapify(pq_t *pq, size_t index)
{
  while (index > 0)
  {
    size_t parent = (index - 1) / 2;
    if (pq->data[parent]->weight <= pq->data[index]->weight)
    {
      break;
    }
    pq_swap_node(&pq->data[parent], &pq->data[index]);
    index = parent;
  }
}

void free_pq(pq_t *pq)
{
  for (size_t i = 0; i < pq->size; i++)
  {
    free(pq->data[i]);
  }
  free(pq->data);
  free(pq);
}

void pq_push(pq_t *pq, pq_node_t *node)
{
  if (pq->size >= pq->capacity)
  {
    pq->capacity *= 2;
    pq_node_t **temp = pq->data;
    pq->data =
        (pq_node_t **)realloc(pq->data, pq->capacity * sizeof(pq_node_t *));
    if (pq->data == NULL)
    {
      free(temp);
      fprintf(stderr, "Memory reallocation failed\n");
      exit(EXIT_FAILURE);
    }
  }
  pq->data[pq->size] = node;
  heapify(pq, pq->size++);
}

pq_node_t *pq_pop(pq_t *pq)
{
  if (pq->size == 0)
  {
    return NULL;
  }
  pq_node_t *min_node = pq->data[0];
  pq->data[0] = pq->data[--pq->size];
  size_t index = 0;
  while (1)
  {
    size_t left = 2 * index + 1;
    size_t right = 2 * index + 2;
    size_t smallest = index;
    if (left < pq->size && pq->data[left]->weight < pq->data[smallest]->weight)
    {
      smallest = left;
    }
    if (right < pq->size && pq->data[right]->weight < pq->data[smallest]->weight)
    {
      smallest = right;
    }
    if (smallest == index)
    {
      break;
    }
    pq_swap_node(&pq->data[index], &pq->data[smallest]);
    index = smallest;
  }
  return min_node;
}

graph_t *createGraph(const int numVertices)
{
  graph_t *graph = malloc(sizeof(graph_t));
  graph->numVertices = numVertices;
  graph->numEdges = 0;
  graph->vertex = malloc(numVertices * sizeof(vertex_t));
  for (int i = 0; i < numVertices; i++)
  {
    graph->vertex[i].edge = NULL;
  }
  return graph;
}

bool freeGraph(graph_t *graph)
{
  if (graph == NULL)
    return false;
  if (graph->vertex == NULL)
  {
    free(graph);
    return true;
  }

  for (int i = 0; i < graph->numVertices; i++)
  {
    edge_t *current = graph->vertex[i].edge;
    while (current != NULL)
    {
      edge_t *temp = current;
      current = current->next;
      free(temp);
    }
  }

  free(graph->vertex);
  free(graph);
  return true;
}

void addEdge(graph_t *graph, const int from, const int to, uint32_t drive_time, uint32_t length, uint8_t speed)
{
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

void remove_edge(graph_t *graph, const int from, const int to)
{
  if (graph == NULL)
    return;
  if (from < 0 || from >= graph->numVertices)
    return;
  if (to < 0 || to >= graph->numVertices)
    return;

  edge_t *current = graph->vertex[from].edge;
  edge_t *previous = NULL;

  while (current != NULL)
  {
    if (current->to == &graph->vertex[to])
    {
      if (previous == NULL)
      {
        graph->vertex[from].edge = current->next;
      }
      else
      {
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

fileExtractResult_t extract_tokens_from_line(const char *line, int max_count)
{
  fileExtractResult_t result;
  result.numberOfReadData = 0;
  token_t *data = malloc(max_count * sizeof(token_t));
  const char *ptr = line;

  while (*ptr != '\0' && result.numberOfReadData < max_count)
  {
    while (*ptr == ' ' || *ptr == '\t')
      ptr++;
    if (*ptr == '\n' || *ptr == '\0')
      break;

    if (*ptr == '"')
    {
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
    }
    else
    {
      char *endPtr;
      double num = strtod(ptr, &endPtr);
      if (ptr != endPtr)
      {
        data[result.numberOfReadData].value.number = num;
        data[result.numberOfReadData].is_string = false;
        result.numberOfReadData++;
        ptr = endPtr;
      }
      else
      {
        break;
      }
    }
  }
  result.data = data;
  return result;
}

void free_tokens(token_t *data, int count)
{
  for (int i = 0; i < count; i++)
  {
    if (data[i].is_string)
    {
      free(data[i].value.string);
    }
  }
  free(data);
}

graph_t *create_graph_from_file(const char *nodes_filename, const char *edges_filename)
{
  FILE *edgesFile = fopen(edges_filename, "r");
  FILE *nodesFile = fopen(nodes_filename, "r");

  if (edgesFile == NULL || nodesFile == NULL)
  {
    perror("Failed to open file");
    return NULL;
  }
  char *edgesLine = malloc(256 * sizeof(char));
  char *nodesLine = malloc(256 * sizeof(char));
  fgets(edgesLine, 256, edgesFile);
  fgets(nodesLine, 256, nodesFile);
  const fileExtractResult_t firstLineDataEdges = extract_tokens_from_line(edgesLine, 5);
  printf("numEges: %d\n", (int)((token_t *)firstLineDataEdges.data)[0].value.number);
  const fileExtractResult_t firstLineDataNodes = extract_tokens_from_line(nodesLine, 5);
  printf("numNodes: %d\n", (int)((token_t *)firstLineDataNodes.data)[0].value.number);
  const int numEdges = (int)((token_t *)firstLineDataEdges.data)[0].value.number;
  const int numNodes = (int)((token_t *)firstLineDataNodes.data)[0].value.number;
  graph_t *graph = createGraph(numNodes);
  free(firstLineDataEdges.data);
  free(firstLineDataNodes.data);

  printf("Processing nodes");
  while (fgets(nodesLine, 256, nodesFile))
  {
    const fileExtractResult_t nodeData = extract_tokens_from_line(nodesLine, 3);
    if (nodeData.numberOfReadData == 3)
    {
      token_t *data = nodeData.data;
      graph->vertex[(int)data[0].value.number].latitude = (float)data[1].value.number;
      graph->vertex[(int)data[0].value.number].longitude = (float)data[2].value.number;
      free_tokens(data, nodeData.numberOfReadData);
    }
  }

  printf("Processing edges");
  while (fgets(edgesLine, 256, edgesFile))
  {
    const fileExtractResult_t edgeData = extract_tokens_from_line(edgesLine, 5);
    if (edgeData.numberOfReadData == 5)
    {
      const token_t *data = edgeData.data;
      addEdge(graph, (int)data[0].value.number, (int)data[1].value.number,
              (uint32_t)data[2].value.number, (uint32_t)data[3].value.number,
              (uint8_t)data[4].value.number);
      free(edgeData.data);
    }
  }

  // validere antall kanter
  if (graph->numEdges != numEdges)
  {
    printf("Warning: Expected %d edges, but found %d edges in the file.\n", numEdges, graph->numEdges);
  }
  fclose(edgesFile);
  fclose(nodesFile);
  free(nodesLine);
  free(edgesLine);
  return graph;
}

predecessor_t *new_predecessor()
{
  predecessor_t *predecessor = malloc(sizeof(predecessor_t));
  predecessor->predecessor = NULL;
  predecessor->distance = INFINITY;
  return predecessor;
}

void init_predecessors(const graph_t *graph, const vertex_t *vertex)
{
  for (int i = 0; i < graph->numVertices; ++i)
  {
    graph->vertex[i].data = new_predecessor();
  }
  ((predecessor_t *)vertex->data)->distance = 0;
}

void dijkstra(const graph_t *graph, const vertex_t *start, const vertex_t *end)
{
  init_predecessors(graph, start);
  pq_t *pq = create_pq(graph->numVertices);
  pq_push(pq, create_pq_node(0, (void *)start));
  bool *visited = calloc(graph->numVertices, sizeof(bool));
  while (pq->size)
  {
    pq_node_t *current = pq_pop(pq);
    vertex_t *current_vertex = current->data;
    free(current);
    int vertex_index = (int)(current_vertex - graph->vertex);
    if (visited[vertex_index])
    {
      continue;
    }
    visited[vertex_index] = true;

    if (current_vertex == end)
    {
      free(visited);
      free_pq(pq);
      return;
    }
    const predecessor_t *current_pred = current_vertex->data;
    edge_t *current_edge = current_vertex->edge;
    while (current_edge != NULL)
    {
      predecessor_t *neighbor_pred = current_edge->to->data;
      const uint32_t new_distance = current_pred->distance + current_edge->length;
      if (new_distance < neighbor_pred->distance)
      {
        neighbor_pred->distance = new_distance;
        neighbor_pred->predecessor = current_vertex;
        pq_push(pq, create_pq_node(new_distance, (void *)current_edge->to));
      }
      current_edge = current_edge->next;
    }
  }

  free(visited);
  free_pq(pq);
}

void preprocess_graph(graph_t *graph)
{
  // sjekke om vi har en preprosessert graf lagret på disk
  // hvis ja, last den inn
  // hvis nei, fortsett til neste steg
  // kjør dijkstra med interessepunktene som startnoder
  // lagre avstand fra hvert interessepunkt til alle noder i grafen
}

void print_graph(const graph_t graph)
{
  printf("%-10s %-10s %-10s\n", "Node", "Forgj", "Dist");
  printf("--------------------------------------------------------\n");
  for (int i = 0; i < graph.numVertices; i++)
  {
    if (graph.vertex[i].data != NULL)
    {
      const predecessor_t *pred = (predecessor_t *)graph.vertex[i].data;
      if (pred->predecessor != NULL)
      {
        const int predIndex = (int)(pred->predecessor - graph.vertex);
        printf("%-10d %-10d %-10d\n", i, predIndex, pred->distance);
      }
      else
      {
        printf("%-10d %-10s %-10s\n", i, " ", "∞");
      }
    }
  }
}

void print_path(const graph_t *graph, const vertex_t *start, const vertex_t *end)
{
  FILE* file = fopen("path.txt", "w");
  if (graph == NULL || start == NULL || end == NULL)
    return;
  predecessor_t *current_pred = end->data;
  if (current_pred->distance == INFINITY)
  {
    printf("No path from start to end.\n");
    return;
  }
  printf("Path from node %lld to node %lld (Distance: %d):\n", start - graph->vertex, end - graph->vertex, current_pred->distance);
  vertex_t *current_vertex = (vertex_t *)end;
  while (current_vertex != NULL)
  {

    fprintf(file, "%f,%f\n", current_vertex->latitude, current_vertex->longitude);
    predecessor_t *pred = current_vertex->data;
    current_vertex = pred->predecessor;
  }
  printf("\n");
}

int main(int argc, char *argv[])
{

  graph_t *graph = create_graph_from_file("noder.txt", "kanter.txt");
  dijkstra(graph, &graph->vertex[1955665], &graph->vertex[4141533]);
  print_path(graph, &graph->vertex[1955665], &graph->vertex[4141533]);

  freeGraph(graph);
  return EXIT_SUCCESS;
}
