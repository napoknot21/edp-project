
from __future__ import annotations

import numpy as np
from collections import defaultdict, deque


def connected_components_boundary_edges(bnd: np.ndarray) -> np.ndarray:
    """
    Label connected components of the boundary-edge graph.
    """
    vertex_to_edges = defaultdict(list)

    for k, (a, b) in enumerate(bnd) :
        
        vertex_to_edges[int(a)].append(k)
        vertex_to_edges[int(b)].append(k)

    labels = -np.ones(len(bnd), dtype=int)
    component = 0

    for start in range(len(bnd)) :

        if labels[start] != -1:
            continue

        q = deque([start])
        labels[start] = component

        while q :

            e_id = q.popleft()
            a, b = bnd[e_id]

            for neigh in vertex_to_edges[int(a)] + vertex_to_edges[int(b)]:
            
                if labels[neigh] == -1:
            
                    labels[neigh] = component
                    q.append(neigh)

        component += 1

    return labels


def boundary_vertices_by_component(bnd: np.ndarray, labels: np.ndarray) -> dict[int, np.ndarray] :
    """
    Group boundary vertices by connected-component label.
    """
    result = {}

    for lab in np.unique(labels) :
    
        edges = bnd[labels == lab]
        result[int(lab)] = np.unique(edges.reshape(-1))
    
    return result
