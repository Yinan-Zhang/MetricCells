//
//  priority_queue.h
//  MetricCells
//
//  Created by Yinan Zhang on 4/20/15.
//  Copyright (c) 2015 Yinan Zhang. All rights reserved.
//

#ifndef MetricCells_priority_queue_h
#define MetricCells_priority_queue_h

#include <queue>    // priority queue
#include <vector>
#include <unordered_set>

namespace data_structure {
    /***********************
     * For storing
     ***********************/
    template <class T>
    struct pqnode
    {
        explicit pqnode<T>(const T& arg_task, double arg_priority) : task(arg_task), priority(arg_priority) {}
        
        bool operator<(const pqnode<T>& rhs) const
        {
            return priority > rhs.priority;
        }
        
        T task;
        double priority;
    };
    
    /***********************
     * priority queue
     ***********************/
    template <class T>
    struct priority_queue {
        std::priority_queue<pqnode<T>, std::vector<pqnode<T>>> pq;
        
        //Insert an element into the prioity queue.
        void push(const T& task, double priority )
        {
            this->pq.emplace(task, priority);
        }
        
        void update(const T& task, double new_priority )
        {
            this->pq.emplace(task, new_priority);
        }
        
        // Return (without removing it) a highest priority element from the priority queue.
        T top() const
        {
            pqnode<T> nd = this->pq.top();
            return nd.task;
        }
        
        // Remove a highest priority element from the priority queue.
        T pop()
        {
            pqnode<T> nd = this->pq.top();
            this->pq.pop();
            return nd.task;
        }
        
        int size() const
        {
            return (int)(this->pq.size());
        }
        
        bool is_empty() const
        {
            return this->pq.empty();
        }
        
    };
}


#endif
