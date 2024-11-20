#ifndef LINEARHEAP_H
#define LINEARHEAP_H

#include "Utility.h"

/**
 * @brief used for degeneracy order, which can be done in O(m).
 * based on linked list
 */
class LinearHeap
{
    struct Node
    {
        ui pre, ne;
    };
    ui *h;
    Node *nodes;

public:
    ui max_range;
    ui min_key;
    ui sz, n; // sz will decrease , but n won't
    LinearHeap() : max_range(0), n(0), sz(0), min_key(0), h(nullptr), nodes(nullptr)
    {
    }
    LinearHeap(ui _range, ui _n) : max_range(_range), n(_n), sz(0), h(nullptr), nodes(nullptr)
    {
        h = new ui[max_range];
        for (ui i = 0; i < max_range; i++)
            h[i] = n;
        nodes = new Node[n];
        // for (ui i = 0; i < n; i++)
        //     nodes[i] = {n, n};
        min_key = max_range;
    }
    LinearHeap(ui _range, ui _n, int pd[]) : max_range(_range), n(_n), sz(0), h(nullptr), nodes(nullptr)
    {
        h = new ui[max_range];
        for (ui i = 0; i < max_range; i++)
            h[i] = n;
        nodes = new Node[n];
        // for (ui i = 0; i < n; i++)
        //     nodes[i] = {n, n};
        min_key = max_range;
        for (ui i = 0; i < n; i++)
            insert(pd[i], i);
    }
    LinearHeap(ui _range, ui _n, vector<ui> &pd) : max_range(_range), n(_n), sz(0), h(nullptr), nodes(nullptr)
    {
        h = new ui[max_range];
        for (ui i = 0; i < max_range; i++)
            h[i] = n;
        nodes = new Node[n];
        // for (ui i = 0; i < n; i++)
        //     nodes[i] = {n, n};
        min_key = max_range;
        for (ui i = 0; i < n; i++)
            insert(pd[i], i);
    }
    ~LinearHeap()
    {
        if (h != nullptr)
            delete[] h;
        if (nodes != nullptr)
            delete[] nodes;
    }
    LinearHeap &operator=(const LinearHeap &other)
    {
        max_range = other.max_range;
        n = other.n;
        sz = other.sz;
        min_key = other.min_key;
        if (h != nullptr)
        {
            delete[] h;
        }
        h = new ui[max_range];
        memcpy(h, other.h, sizeof(ui) * max_range);
        if (nodes != nullptr)
        {
            delete[] nodes;
        }
        nodes = new Node[n];
        memcpy(nodes, other.nodes, sizeof(Node) * n);
        return *this;
    }
    // key is the degree of a vertex; val is the index of a vertex
    inline void insert(ui key, ui val)
    {
        min_key = min(min_key, key);
        sz++;
        nodes[val].pre = key + n;
        nodes[val].ne = h[key];
        if (h[key] < n)
            nodes[h[key]].pre = val;
        h[key] = val;
    }
    ui get_min_key()
    {
        assert(sz > 0);
        while (h[min_key] >= n)
            min_key++;
        return min_key;
    }
    inline void delete_node(ui id)
    {
        sz--;
        auto &cur = nodes[id];
        if (cur.pre >= n) // key = cur.pre-n
            h[cur.pre - n] = cur.ne;
        else
            nodes[cur.pre].ne = cur.ne;
        if (cur.ne < n)
            nodes[cur.ne].pre = cur.pre;
    }
    inline void decrease(ui key_now, ui id)
    {
        // delete_node(id);
        // insert(key_now, id);
        int pre = nodes[id].pre, ne = nodes[id].ne;
        if (pre >= n)
        {
            h[pre - n] = ne;
        }
        else
            nodes[pre].ne = ne;
        if (ne < n)
        {
            nodes[ne].pre = pre;
        }
        min_key = min(min_key, key_now);
        ui nxt = h[key_now];
        nodes[id] = {key_now + n, nxt};
        if (nxt < n)
        {
            nodes[nxt].pre = id;
        }
        h[key_now] = id;
    }
    ui get_min_node()
    {
        return h[get_min_key()];
    }
};

#endif