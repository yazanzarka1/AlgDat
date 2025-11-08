
typedef struct PriorityQueueNode {
  int16_t byte;
  uint32_t freq;
  struct PriorityQueueNode *left, *right;
} pq_node_t;

typedef struct PriorityQueue {
  pq_node_t **data;
  size_t size;
  size_t capacity;
} pq_t;

// Heap og prioritetskø implementasjon for Huffman-koding
void pq_swap_node(pq_node_t **a, pq_node_t **b) {
  pq_node_t *t = *a;
  *a = *b;
  *b = t;
}
pq_node_t *create_pq_node(int16_t byte, uint32_t freq, pq_node_t *left, pq_node_t *right) {
  pq_node_t *node = (pq_node_t *) malloc(sizeof(pq_node_t));
  node->byte = byte;
  node->freq = freq;
  node->left = left;
  node->right = right;
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
    if (pq->data[parent]->freq <= pq->data[index]->freq) {
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
    pq->data = (pq_node_t **) realloc(pq->data, pq->capacity * sizeof(pq_node_t *));
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
    if (left < pq->size && pq->data[left]->freq < pq->data[smallest]->freq) {
      smallest = left;
    }
    if (right < pq->size && pq->data[right]->freq < pq->data[smallest]->freq) {
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



int int main(int argc, char *argv[])
{
  return EXIT_SUCCESS;
}
