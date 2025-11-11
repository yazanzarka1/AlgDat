#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
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
    double latitude, longitude;
    void *data;
    char *poi_name;
    uint8_t poi_flags;
} vertex_t;

typedef struct PreprocessedLandmarks {
    int64_t num_landmarks;
    int64_t num_vertices;
    uint64_t *weights;
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

typedef struct NearestPOI {
    vertex_t *vertex;
    uint32_t drive_time;
} nearest_poi_t;

void calculate_time_difference(struct timespec start, struct timespec end, int *sec, long *nsec) {
    *sec = (int) (end.tv_sec - start.tv_sec);
    *nsec = end.tv_nsec - start.tv_nsec;
    if (*nsec < 0) {
        --(*sec);
        *nsec += 1000000000L;
    }
}

void format_time(char *buf, size_t size, int sec, long nsec) {
    if (sec >= 1.0) {
        snprintf(buf, size, "%.6f s", sec + nsec / 1000000000.0);
    } else if (nsec > 50000000) {
        snprintf(buf, size, "%.3f ms", nsec / 1000000.0);
    } else if (nsec > 50000) {
        snprintf(buf, size, "%.3f µs", nsec / 1000.0);
    } else {
        snprintf(buf, size, "%ld ns", nsec);
    }
}

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
        graph->vertex[i].poi_name = NULL;
        graph->vertex[i].poi_flags = 0;
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
        if (graph->vertex[i].poi_name != NULL) {
            free(graph->vertex[i].poi_name);
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

    int numEdges = 0;
    int numNodes = 0;
    int numPOI = 0;
    char line[1024];

    if (fgets(line, sizeof(line), edgesFile)) {
        fileExtractResult_t result = extract_tokens_from_line(line, 1);
        if (result.numberOfReadData >= 1) {
            token_t *tokens = (token_t *)result.data;
            if (!tokens[0].is_string) {
                numEdges = (int)tokens[0].value.number;
            }
            free_tokens(tokens, result.numberOfReadData);
        }
    }
    
    if (fgets(line, sizeof(line), nodesFile)) {
        fileExtractResult_t result = extract_tokens_from_line(line, 1);
        if (result.numberOfReadData >= 1) {
            token_t *tokens = (token_t *)result.data;
            if (!tokens[0].is_string) {
                numNodes = (int)tokens[0].value.number;
            }
            free_tokens(tokens, result.numberOfReadData);
        }
    }
    
    if (fgets(line, sizeof(line), poiFile)) {
        fileExtractResult_t result = extract_tokens_from_line(line, 1);
        if (result.numberOfReadData >= 1) {
            token_t *tokens = (token_t *)result.data;
            if (!tokens[0].is_string) {
                numPOI = (int)tokens[0].value.number;
            }
            free_tokens(tokens, result.numberOfReadData);
        }
    }

    printf("%d %d %d\n", numEdges, numNodes, numPOI);

    graph_t *graph = createGraph(numNodes);

    printf("Processing POI\n");
    while (fgets(line, sizeof(line), poiFile)) {
        fileExtractResult_t result = extract_tokens_from_line(line, 3);
        if (result.numberOfReadData >= 2) {
            token_t *tokens = (token_t *)result.data;
            if (!tokens[0].is_string && !tokens[1].is_string) {
                int node_id = (int)tokens[0].value.number;
                uint8_t bit_flag = (uint8_t)tokens[1].value.number;
                
                if (node_id >= 0 && node_id < graph->numVertices) {
                    graph->vertex[node_id].poi_flags = bit_flag;
                    
                    if (result.numberOfReadData >= 3 && tokens[2].is_string) {
                        graph->vertex[node_id].poi_name = malloc(strlen(tokens[2].value.string) + 1);
                        if (graph->vertex[node_id].poi_name) {
                            strcpy(graph->vertex[node_id].poi_name, tokens[2].value.string);
                        }
                    }
                }
            }
            free_tokens(tokens, result.numberOfReadData);
        }
    }
    fclose(poiFile);

    printf("Processing nodes\n");
    while (fgets(line, sizeof(line), nodesFile)) {
        fileExtractResult_t result = extract_tokens_from_line(line, 3);
        if (result.numberOfReadData >= 3) {
            token_t *tokens = (token_t *)result.data;
            if (!tokens[0].is_string && !tokens[1].is_string && !tokens[2].is_string) {
                int node_id = (int)tokens[0].value.number;
                double latitude = tokens[1].value.number;
                double longitude = tokens[2].value.number;
                
                if (node_id >= 0 && node_id < graph->numVertices) {
                    graph->vertex[node_id].latitude = latitude;
                    graph->vertex[node_id].longitude = longitude;
                }
            }
            free_tokens(tokens, result.numberOfReadData);
        }
    }
    fclose(nodesFile);

    printf("Processing edges\n");
    while (fgets(line, sizeof(line), edgesFile)) {
        fileExtractResult_t result = extract_tokens_from_line(line, 5);
        if (result.numberOfReadData >= 5) {
            token_t *tokens = (token_t *)result.data;
            if (!tokens[0].is_string && !tokens[1].is_string && 
                !tokens[2].is_string && !tokens[3].is_string && !tokens[4].is_string) {
                int from_node = (int)tokens[0].value.number;
                int to_node = (int)tokens[1].value.number;
                uint32_t drive_time = (uint32_t)tokens[2].value.number;
                uint32_t length = (uint32_t)tokens[3].value.number;
                uint8_t speed = (uint8_t)tokens[4].value.number;
                
                if (from_node >= 0 && from_node < graph->numVertices &&
                    to_node >= 0 && to_node < graph->numVertices) {
                    addEdge(graph, from_node, to_node, drive_time, length, speed);
                }
            }
            free_tokens(tokens, result.numberOfReadData);
        }
    }
    fclose(edgesFile);

    if (graph->numEdges != numEdges) {
        printf("Warning: Expected %d edges, but found %d edges in the file.\n",
               numEdges,
               graph->numEdges);
    }
    
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

void dijkstra(const graph_t *graph, const vertex_t *start, const vertex_t *end, int *nodes_visited) {
    init_predecessors(graph, start);
    int count = 0;
    pq_t *pq = create_pq(graph->numVertices);
    bool *visited = calloc(graph->numVertices, sizeof(bool));
    pq_push(pq, create_pq_node(0, (void *) start));
    while (pq->size) {
        pq_node_t *current = pq_pop(pq);
        vertex_t *current_vertex = current->data;
        free(current);
        int vertex_index = (int) (current_vertex - graph->vertex);
        if (visited[vertex_index]) {
            continue;
        }
        visited[vertex_index] = true;
        count++;

        if (end != NULL && current_vertex == end) {
            free(visited);
            free_pq(pq);
            if (nodes_visited)
                *nodes_visited = count;
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
    if (nodes_visited)
        *nodes_visited = count;
}

nearest_poi_t *dijkstra_nearest_pois(const graph_t *graph, const vertex_t *start, uint8_t poi_type_flags, int n,
                                     int *nodes_visited) {
    init_predecessors(graph, start);
    int count = 0;
    pq_t *pq = create_pq(graph->numVertices);
    bool *visited = calloc(graph->numVertices, sizeof(bool));

    nearest_poi_t *nearest = malloc(n * sizeof(nearest_poi_t));
    for (int i = 0; i < n; i++) {
        nearest[i].vertex = NULL;
        nearest[i].drive_time = INFINITY;
    }
    int found_count = 0;

    pq_push(pq, create_pq_node(0, (void *) start));

    while (pq->size && found_count < n) {
        pq_node_t *current = pq_pop(pq);
        vertex_t *current_vertex = current->data;
        free(current);
        int vertex_index = (int) (current_vertex - graph->vertex);
        if (visited[vertex_index]) {
            continue;
        }
        visited[vertex_index] = true;
        count++;

        if (is_flag_set(current_vertex->poi_flags, poi_type_flags)) {
            predecessor_t *pred = current_vertex->data;
            nearest[found_count].vertex = current_vertex;
            nearest[found_count].drive_time = pred->total_drive_time;
            found_count++;
            if (found_count >= n) {
                break;
            }
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
    if (nodes_visited)
        *nodes_visited = count;

    return nearest;
}

int compare_nearest_pois(const void *a, const void *b) {
    nearest_poi_t *poi_a = (nearest_poi_t *) a;
    nearest_poi_t *poi_b = (nearest_poi_t *) b;
    if (poi_a->drive_time < poi_b->drive_time)
        return -1;
    if (poi_a->drive_time > poi_b->drive_time)
        return 1;
    return 0;
}

nearest_poi_t *find_nearest_pois(graph_t *graph, vertex_t *start, uint8_t poi_type_flags, int n, int *nodes_visited) {
    return dijkstra_nearest_pois(graph, start, poi_type_flags, n, nodes_visited);
}

void print_nearest_pois(graph_t *graph, nearest_poi_t *nearest, int n, const char *filename, vertex_t *start) {
    FILE *file = fopen(filename, "w");

    if (graph == NULL)
        return;

    if (file && start != NULL) {
        fprintf(file, "%.6f,%.6f\n", start->latitude, start->longitude);
    }

    printf("\nNearest %d POIs:\n", n);
    printf("%-40s %-10s %-15s %-15s\n", "Name", "Node ID", "Drive Time", "Lat,Lon");
    printf("-----------------------------------------------------------------------------------\n");

    for (int i = 0; i < n; i++) {
        if (nearest[i].vertex == NULL)
            break;

        int node_id = (int) (nearest[i].vertex - graph->vertex);
        const char *name = nearest[i].vertex->poi_name ? nearest[i].vertex->poi_name : "N/A";
        uint32_t drive_time = nearest[i].drive_time;

        uint32_t hours = drive_time / 360000;
        uint32_t minutes = (drive_time % 360000) / 6000;
        uint32_t seconds = (drive_time % 6000) / 100;

        printf("%-40s %-10d %02u:%02u:%02u       %.6f,%.6f\n",
               name, node_id, hours, minutes, seconds,
               nearest[i].vertex->latitude, nearest[i].vertex->longitude);

        if (file) {
            fprintf(file, "%.6f,%.6f\n",
                    nearest[i].vertex->latitude, nearest[i].vertex->longitude);
        }
    }

    if (file)
        fclose(file);
}

uint64_t ALT_heuristic(const graph_t *graph,
                       vertex_t *current,
                       vertex_t *end,
                       preprocessed_landmarks_t **landmarks,
                       int num_landmarks) {
    int64_t best = 0;
    int current_index = (int) (current - graph->vertex);

    int end_index = (int) (end - graph->vertex);

    for (int i = 0; i < num_landmarks; i++) {
        int64_t fwd = (int64_t) landmarks[0]->weights[i * graph->numVertices + end_index] - (int64_t) landmarks[0]->
                      weights[i * graph->numVertices + current_index];
        int64_t bwd = (int64_t) landmarks[1]->weights[i * graph->numVertices + current_index] - (int64_t) landmarks[1]->
                      weights[i * graph->numVertices + end_index];
        if (fwd > best)
            best = fwd;
        if (bwd > best)
            best = bwd;
    }
    if (best < 0)
        best = 0;
    return (uint64_t) best;
}

void ALT(graph_t *graph,
         vertex_t *start,
         vertex_t *end,
         preprocessed_landmarks_t **landmarks,
         int *landmark_indices,
         int num_landmarks,
         int *nodes_visited) {
    init_predecessors(graph, start);
    pq_t *pq = create_pq(graph->numVertices);
    pq_push(pq, create_pq_node(0, (void *) start));
    bool *visited = calloc(graph->numVertices, sizeof(bool));
    int count = 0;
    while (pq->size) {
        pq_node_t *current = pq_pop(pq);
        vertex_t *current_vertex = current->data;
        free(current);
        int vertex_index = (int) (current_vertex - graph->vertex);
        if (visited[vertex_index]) {
            continue;
        }
        visited[vertex_index] = true;
        count++;

        if (current_vertex == end) {
            free(visited);
            free_pq(pq);
            if (nodes_visited)
                *nodes_visited = count;
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
    if (nodes_visited)
        *nodes_visited = count;
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
        landmarks[0]->weights = malloc(numLandmarks * numVertices * sizeof(uint64_t));
        fread(landmarks[0]->weights, sizeof(uint64_t), numLandmarks * numVertices, file);

        // distanser for bakover
        landmarks[1] = malloc(sizeof(preprocessed_landmarks_t));
        landmarks[1]->num_landmarks = numLandmarks;
        landmarks[1]->num_vertices = numVertices;
        landmarks[1]->weights = malloc(numLandmarks * numVertices * sizeof(uint64_t));
        fread(landmarks[1]->weights, sizeof(uint64_t), numLandmarks * numVertices, file);
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
    landmarks[0]->weights = malloc(num_landmarks * graph->numVertices * sizeof(uint64_t));
    // regne ut distanser mellom landmarks og alle noder
    for (int i = 0; i < num_landmarks; i++) {
        dijkstra(graph, &graph->vertex[landmarks_indices[i]], NULL, NULL);
        for (int j = 0; j < graph->numVertices; j++) {
            predecessor_t *pred = graph->vertex[j].data;
            landmarks[0]->weights[i * graph->numVertices + j] = pred->total_drive_time;
        }
    }

    // distanser for bakover
    landmarks[1] = malloc(sizeof(preprocessed_landmarks_t));
    landmarks[1]->num_landmarks = num_landmarks;
    landmarks[1]->num_vertices = graph->numVertices;
    landmarks[1]->weights = malloc(num_landmarks * graph->numVertices * sizeof(uint64_t));
    graph_t *reversed_graph = reverse_graph((graph_t *) graph);
    for (int i = 0; i < num_landmarks; i++) {
        dijkstra(reversed_graph, &reversed_graph->vertex[landmarks_indices[i]], NULL, NULL);
        for (int j = 0; j < reversed_graph->numVertices; j++) {
            predecessor_t *pred = reversed_graph->vertex[j].data;
            landmarks[1]->weights[i * reversed_graph->numVertices + j] = pred->total_drive_time;
        }
    }
    freeGraph(reversed_graph);
    // skrive til fil
    fwrite(landmarks[0]->weights, 1, num_landmarks * graph->numVertices * sizeof(uint64_t), file);
    fwrite(landmarks[1]->weights, 1, num_landmarks * graph->numVertices * sizeof(uint64_t), file);
    fclose(file);
    return landmarks;
}

void free_preprocessed_landmarks(preprocessed_landmarks_t **landmarks, int num_landmarks) {
    if (landmarks == NULL)
        return;
    for (int i = 0; i < num_landmarks; i++) {
        if (landmarks[i] != NULL) {
            if (landmarks[i]->weights != NULL) {
                free(landmarks[i]->weights);
            }
            free(landmarks[i]);
        }
    }
    free(landmarks);
}

void free_coordinates(char **coords, int32_t n) {
    if (!coords) return;
    for (int32_t i = 0; i < n; ++i) free(coords[i]);
    free(coords);
}

void print_path(const graph_t *graph, const vertex_t *start, const vertex_t *end) {
    if (!graph || !start || !end) return;

    predecessor_t *end_pred = end->data;
    if (!end_pred || end_pred->total_drive_time == INFINITY) {
        printf("No path from start to end.\n");
        return;
    }

    printf("Path from node %lld to node %lld (Drive time: %u):\n",
           (long long) (start - graph->vertex),
           (long long) (end - graph->vertex),
           end_pred->total_drive_time);

    int32_t number_of_vertices = 0;
    const vertex_t *cv = end;
    while (cv) {
        ++number_of_vertices;
        const predecessor_t *pred = cv->data;
        cv = pred ? pred->predecessor : NULL;
    }
    if (number_of_vertices <= 0) {
        printf("Empty path.\n");
        return;
    }

    char **coordinates = (char **) malloc((size_t) number_of_vertices * sizeof *coordinates);
    if (!coordinates) {
        fprintf(stderr, "OOM: coordinates array\n");
        return;
    }

    vertex_t *current_vertex = (vertex_t *) end;
    for (int i = number_of_vertices - 1; i >= 0; --i) {
        coordinates[i] = (char *) malloc(256);
        if (!coordinates[i]) {
            fprintf(stderr, "OOM: coordinates[%d]\n", i);
            free_coordinates(coordinates, number_of_vertices);
            return;
        }

        snprintf(coordinates[i], 256, "%.6f, %.6f",
                 current_vertex->latitude, current_vertex->longitude);

        predecessor_t *pred = current_vertex->data;
        current_vertex = pred ? pred->predecessor : NULL;
    }

    FILE *file = fopen("path.txt", "w");
    if (!file) {
        perror("fopen");
        free_coordinates(coordinates, number_of_vertices);
        return;
    }

    for (int i = 0; i < number_of_vertices; ++i) {
        fputs(coordinates[i], file);
        fputc('\n', file);
    }

    fclose(file);
    free_coordinates(coordinates, number_of_vertices);
}

int main(int argc, char *argv[]) {
    printf("\n========== GRAPH ROUTING SYSTEM ==========\n\n");

    if (argc < 2) {
        printf("Usage:\n");
        printf("  %s nearest <node_id> <poi_type>       - Find 5 nearest POIs\n", argv[0]);
        printf("  %s compare <start> <end>              - Compare Dijkstra vs ALT\n", argv[0]);
        printf("  %s <algorithm> <start> <end>          - Run specific algorithm\n\n", argv[0]);
        printf("Algorithms: dijkstra, alt\n");
        printf("POI types: 1=place, 2=gas, 4=charging, 8=dining, 16=drink, 32=hotel\n");
        printf("\n==========================================\n\n");
        return EXIT_FAILURE;
    }

    if (strcmp(argv[1], "nearest") == 0 && argc < 4) {
        printf("Usage: %s nearest <node_id> <poi_type>\n", argv[0]);
        printf("POI types: 1=place, 2=gas, 4=charging, 8=dining, 16=drink, 32=hotel\n");
        return EXIT_FAILURE;
    }

    if (strcmp(argv[1], "compare") == 0 && argc < 4) {
        printf("Usage: %s compare <start_node> <end_node>\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (argc == 4 && strcmp(argv[1], "dijkstra") != 0 && strcmp(argv[1], "alt") != 0 && strcmp(argv[1], "nearest") != 0
        && strcmp(argv[1], "compare") != 0) {
        printf("Unknown algorithm: %s\n", argv[1]);
        printf("Available algorithms: dijkstra, alt\n");
        printf("Or use: nearest, compare\n");
        return EXIT_FAILURE;
    }

    printf("Loading graph data...\n");
    int landmarks_indices[] = {
        7826348, //Trondheim
        5519970, // Agder
        5518023, // Moxy Bergen
        2246647, // Kirkenes
        5459754, // Horsnes
        925129, // Helsinki,
        545911 // København
    };
    graph_t *graph = initialize_graph("noder.txt", "kanter.txt", "interessepkt.txt");
    preprocessed_landmarks_t **landmarks = preprocess_graph(graph, landmarks_indices, 5);
    printf("\n==========================================\n\n");

    if (strcmp(argv[1], "nearest") == 0) {
        int node_id = atoi(argv[2]);
        int poi_type = atoi(argv[3]);

        if (node_id < 0 || node_id >= graph->numVertices) {
            printf("Invalid node ID\n");
            freeGraph(graph);
            free_preprocessed_landmarks(landmarks, 2);
            return EXIT_FAILURE;
        }

        vertex_t *start = &graph->vertex[node_id];

        printf("\nSearching for 5 nearest POIs (type %d) from node %d...\n\n", poi_type, node_id);

        int nodes = 0;
        nearest_poi_t *nearest = find_nearest_pois(graph, start, (uint8_t) poi_type, 5, &nodes);
        print_nearest_pois(graph, nearest, 5, "nearest_pois.txt", start);

        printf("\nNodes visited: %d\n", nodes);
        printf("Result saved to: nearest_pois.txt\n");

        free(nearest);
    } else if (strcmp(argv[1], "compare") == 0) {
        int start_id = atoi(argv[2]);
        int end_id = atoi(argv[3]);

        if (start_id < 0 || start_id >= graph->numVertices || end_id < 0 || end_id >= graph->numVertices) {
            printf("Invalid node IDs\n");
            freeGraph(graph);
            free_preprocessed_landmarks(landmarks, 2);
            return EXIT_FAILURE;
        }

        vertex_t *start = &graph->vertex[start_id];
        vertex_t *end = &graph->vertex[end_id];

        printf("\n=== DIJKSTRA ===\n");
        struct timespec dijk_start, dijk_end;
        int dijk_nodes = 0;

        init_predecessors(graph, start);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &dijk_start);
        dijkstra(graph, start, end, &dijk_nodes);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &dijk_end);

        int dijk_sec;
        long dijk_nsec;
        calculate_time_difference(dijk_start, dijk_end, &dijk_sec, &dijk_nsec);
        char dijk_time_buf[64];
        format_time(dijk_time_buf, sizeof(dijk_time_buf), dijk_sec, dijk_nsec);

        predecessor_t *pred = end->data;
        uint32_t dijk_distance = pred->total_drive_time;
        uint32_t hours = dijk_distance / 360000;
        uint32_t minutes = (dijk_distance % 360000) / 6000;
        uint32_t seconds = (dijk_distance % 6000) / 100;

        printf("Nodes visited: %d\n", dijk_nodes);
        printf("Time: %s\n", dijk_time_buf);
        printf("Drive time: %02u:%02u:%02u\n", hours, minutes, seconds);
        print_path(graph, start, end);

        printf("\n=== ALT ===\n");
        struct timespec alt_start, alt_end;
        int alt_nodes = 0;

        init_predecessors(graph, start);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &alt_start);
        ALT(graph, start, end, landmarks, landmarks_indices, 5, &alt_nodes);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &alt_end);

        int alt_sec;
        long alt_nsec;
        calculate_time_difference(alt_start, alt_end, &alt_sec, &alt_nsec);
        char alt_time_buf[64];
        format_time(alt_time_buf, sizeof(alt_time_buf), alt_sec, alt_nsec);

        pred = end->data;
        uint32_t alt_distance = pred->total_drive_time;
        hours = alt_distance / 360000;
        minutes = (alt_distance % 360000) / 6000;
        seconds = (alt_distance % 6000) / 100;

        printf("Nodes visited: %d\n", alt_nodes);
        printf("Time: %s\n", alt_time_buf);
        printf("Drive time: %02u:%02u:%02u\n", hours, minutes, seconds);
        print_path(graph, start, end);

        printf("\n=== COMPARISON ===\n");
        printf("Dijkstra nodes: %d, ALT nodes: %d (%.2f%% reduction)\n",
               dijk_nodes, alt_nodes, 100.0 * (dijk_nodes - alt_nodes) / dijk_nodes);
        double dijk_time_sec = dijk_sec + dijk_nsec / 1000000000.0;
        double alt_time_sec = alt_sec + alt_nsec / 1000000000.0;
        printf("Dijkstra time: %s, ALT time: %s (%.2fx speedup)\n",
               dijk_time_buf, alt_time_buf, dijk_time_sec / alt_time_sec);
    } else {
        int start_id = atoi(argv[2]);
        int end_id = atoi(argv[3]);

        if (start_id < 0 || start_id >= graph->numVertices || end_id < 0 || end_id >= graph->numVertices) {
            printf("Invalid node IDs\n");
            freeGraph(graph);
            free_preprocessed_landmarks(landmarks, 2);
            return EXIT_FAILURE;
        }

        vertex_t *start = &graph->vertex[start_id];
        vertex_t *end = &graph->vertex[end_id];
        const char *alg = argv[1];

        printf("\nRoute from node %d to node %d\n", start_id, end_id);
        printf("Algorithm: %s\n\n", alg);

        struct timespec ts_start, ts_end;
        int nodes = 0;

        if (strcmp(alg, "dijkstra") == 0) {
            init_predecessors(graph, start);
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_start);
            dijkstra(graph, start, end, &nodes);
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_end);
        } else {
            init_predecessors(graph, start);
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_start);
            ALT(graph, start, end, landmarks, landmarks_indices, 5, &nodes);
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_end);
        }

        int sec;
        long nsec;
        calculate_time_difference(ts_start, ts_end, &sec, &nsec);
        char time_buf[64];
        format_time(time_buf, sizeof(time_buf), sec, nsec);

        predecessor_t *pred = end->data;
        uint32_t hours = pred->total_drive_time / 360000;
        uint32_t minutes = (pred->total_drive_time % 360000) / 6000;
        uint32_t seconds = (pred->total_drive_time % 6000) / 100;

        printf("Nodes visited: %d\n", nodes);
        printf("Time: %s\n", time_buf);
        printf("Drive time: %02u:%02u:%02u\n", hours, minutes, seconds);
        print_path(graph, start, end);
    }

    freeGraph(graph);
    free_preprocessed_landmarks(landmarks, 2);
    return EXIT_SUCCESS;
}
