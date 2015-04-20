#ifndef CSPACE_DECOMPOSER_H
#define CSPACE_DECOMPOSER_H

#include <vector>
#include <list>
#include <array>

#include "../math/geoND.h"
#include "../robotics/obstacles.h"

namespace algorithms {
    template <typename IROBOT> class Planner;
    static constexpr int NOT_FOUND = -1;
    
    template <typename IROBOT>
    class Boundary {
    public:
        static constexpr int DIM = IROBOT::DIM;
        using CONFIG = typename IROBOT::CONFIG;
        
        Boundary(int u, int v) : side1(u), side2(v) {}
        
        int otherside(int side) const {
            if (side == side1)
                return side2;
            if (side == side2)
                return side1;
            return NOT_FOUND;
        }
        
        
        /**
         * Generate all configs on this boundary
         */
        void get_boundary_configs(Planner<IROBOT>& planner, double interval, int& begin, int& end) const {
            if (begin_ != NOT_FOUND) {
                begin = begin_;
                end = end_;
                return;
            }
            
            int num_each_dim[DIM], number[DIM];
            for (int i = 0; i < DIM; ++i) {
                num_each_dim[i] = (upper_bound[i] - lower_bound[i]) / interval;
                number[i] = 0;
            }
            ND::vec<DIM> coord;
            for (int i = 0; i < DIM; ++i) {
                coord[i] = 0.0;
            }
            while (true) {
                for (int i = 0; i < DIM; ++i) {
                    coord[i] = number[i] * interval + lower_bound[i];
                }
                int node_index = planner.get_new_info(CONFIG(coord));
                if (begin_ == NOT_FOUND)
                    begin_ = node_index;
                end_ = node_index + 1;
                int j;
                for (j = DIM - 1; j >= 0; --j) {
                    if (number[j] + 1 <= num_each_dim[j])
                        break;
                }
                if (j == -1)
                    break;
                number[j] += 1;
                for (++j; j < DIM; ++j) {
                    //std::cout << ' ';
                    number[j] = 0;
                }
            }
            begin = begin_;
            end = end_;
        }
        void reset() {
            begin_ = end_ = NOT_FOUND;
        }
        double get_weight() const { return weight; }
        void set_weight(double w) { weight = w; }
        int get_fix_dim() const { return fix_dim; }
        int get_side1() const { return side1; }
        int get_side2() const { return side2; }
        void set_fix_dim(int dim) { fix_dim = dim; }
        void set_upper_bound(int index, double bound) { upper_bound[index] = bound; }
        void set_lower_bound(int index, double bound) { lower_bound[index] = bound; }
        const std::array<double, DIM>& get_upper_bound() const { return upper_bound; }
        const std::array<double, DIM>& get_lower_bound() const { return lower_bound; }
    private:
        std::array<double, DIM> upper_bound;
        std::array<double, DIM> lower_bound;
        double weight;
        int side1;
        int side2;
        int fix_dim = NOT_FOUND;
        mutable int begin_ = NOT_FOUND;
        mutable int end_ = NOT_FOUND;
    };
    
    
    template <typename IROBOT>
    class Cell
    {
    public:
        static constexpr int DIM = IROBOT::DIM;
        static constexpr int NUM_SPLITS = 1 << DIM;
        using CONFIG = typename IROBOT::CONFIG;

        explicit Cell( const ND::sphere<DIM>& s ) : c_(s) {boundaries.clear(); boundaries.shrink_to_fit(); }

        Cell& operator=(const Cell<IROBOT>& rhs)
        {
            c_ = rhs.c_;
            is_visited = rhs.is_visited;
            potential_value = rhs.potential_value;
            boundaries = rhs.boundaries;
            return *this;
        }
        
        const ND::vec<DIM>& center() const { return c_.center(); }
        
        double radius() const { return c_.radius(); }
        
        const ND::sphere<DIM>& sphere() const { return c_; }
        
        /**
         * Generate all corner configs
         */
        void get_corner_configs(std::vector<CONFIG>& result) const {
            ND::vec<DIM> center;
            double radius = this->radius();
            for (int i = 0; i < NUM_SPLITS; ++i) {
                center = this->center();
                for (int j = 0; j < DIM; ++j) {
                    if (((1 << j) & i) > 0)
                        center[j] += radius;
                    else
                        center[j] -= radius;
                }
                result.emplace_back(center);
            }
        }
        
        void add_boundary(int boundary) {
            boundaries.emplace_back(boundary);
        }
        
        const std::vector<int>& get_boundaries() const { return boundaries; }
        
        bool contains( const ND::vec<DIM>& point ) const { return c_.contains(point); }
        
        bool on_boundary( const ND::vec<DIM>& point, double tolerance ) const { return c_.on_boundary(point, tolerance); }
        
        bool visited() const { return is_visited; }
        
        void set_visited() { is_visited = true; }
        
        void shrink_to_fit() { boundaries.shrink_to_fit(); }
        
        double get_potential() const { return potential_value; }
        void set_potential(double potential) { potential_value = potential; }

    private:
        ND::sphere<DIM> c_;
        std::vector<int> boundaries;
        bool is_visited = false;
        double potential_value = std::numeric_limits<double>::infinity();
    };
    
    enum class SamplingMethod {UNIFORM, KD, GRID};
    
    template<size_t DIM>
    class TreeNode {
    public:
        static constexpr int NUM_SPLITS = 1 << DIM;
        static constexpr int FREE_MASK = 1;
        static constexpr int COVERED_MASK = 2;
        TreeNode(const ND::vec<DIM>& c, int arg_depth) : center_{c}, children{NOT_FOUND}, depth_(arg_depth) {}
        const ND::vec<DIM>& center() const { return center_; }
        int depth() const { return depth_; }
        int get_children() const { return children; }
        void set_children(int first_child) { children = first_child; }
        int get_cell() const { return cell; }
        void set_cell(int cell_number) { cell = cell_number; }
        bool is_free() const { return info & FREE_MASK; }
        void set_free() { info |= FREE_MASK; }
        bool is_covered() const { return info & COVERED_MASK; }
        void set_covered() { info |= COVERED_MASK; }
    private:
        const ND::vec<DIM> center_;
        const int depth_;
        int       children = NOT_FOUND;
        int       cell = NOT_FOUND;
        unsigned char info = 0;
    };
    
    
    
    template<size_t DIM>
    struct Subspace
    {
        Subspace( std::array<robotics::Range, DIM>& cfg_ranges )
        {
            this->ranges = cfg_ranges;
        }
        
        std::array<robotics::Range, DIM> ranges;
    };
    
    
    
    
    template<typename IROBOT>
    class KDDecomposer final {
    public:
        static constexpr int DIM = IROBOT::DIM;
        static constexpr int NUM_SPLITS = 1 << DIM;
        static constexpr int PENETRATION_CONSTANT = IROBOT::PENETRATION_CONSTANT;
        static constexpr bool MERGE_CELLS = IROBOT::MERGE_CELLS;
        using CONFIG = typename IROBOT::CONFIG;
        KDDecomposer(IROBOT& robot,
                     robotics::ObstManager& obstacle_manager,
                     double epsilon);
        
        // @return a vector of indices of boundary cells. ( cells that partially contain obstacles )
        std::vector<int>  ShallowDecompose(double min_radius);
        
        // Decompose the subspace.
        void  DecomposeSubspace(int subspace_index);
        
        void DecomposeSpace();
        //void ScanCollisionSpace();
        std::vector<Cell<IROBOT>>&& move_free_cells() { return std::move(cells); }
        std::vector<Boundary<IROBOT>>&& move_boundaries() { return std::move(all_boundaries); }
        const std::vector<Cell<IROBOT>>& get_free_cells() const { return cells; }

    private:
        void build_edges();
        bool delete_unsafe_cells(int node);
        bool merge_free_cells(int node);
        void create_free_cells(int node);
        void clean_tree();
        
        /**
         * Given a point on the corner, find all cells that contain this point
         */
        void ContainedCells(const ND::vec<DIM>& point, std::vector<int>& result) const {
            double epsilon = 1e-14;
            int stack[1024]; // Set to the maximum depth * NUM_SPLIT to avoid overflow
            stack[0] = root_index;
            int top = 0;
            while (top >= 0) {
                int node_index = stack[top--];
                const TreeNode<DIM>& node = get_node(node_index);
                const ND::vec<DIM>& center = node.center();
                double radius = radius_array[node.depth()];
                
                bool contained = true;
                // Unroll loop for efficiency
                if (DIM == 3) {
                    contained = fabs(point[0] - center[0]) <= radius + epsilon
                    && fabs(point[1] - center[1]) <= radius + epsilon
                    && fabs(point[2] - center[2]) <= radius + epsilon;
                }
                else if (DIM == 2) {
                    contained = fabs(point[0] - center[0]) <= radius + epsilon
                    && fabs(point[1] - center[1]) <= radius + epsilon;
                }
                else if (DIM == 4) {
                    contained = fabs(point[0] - center[0]) <= radius + epsilon
                    && fabs(point[1] - center[1]) <= radius + epsilon
                    && fabs(point[2] - center[2]) <= radius + epsilon
                    && fabs(point[3] - center[3]) <= radius + epsilon;
                }
                else {
                    for (int d = 0; contained && d < DIM; ++d)
                        contained = fabs(point[d] - center[d]) <= radius + epsilon;
                }
                if (contained) {
                    if (node.is_free()) {
                        result.emplace_back(node.get_cell());
                    } else if (get_node(node_index).get_children() != NOT_FOUND) {
                        for (int i = 0; i < NUM_SPLITS; ++i)
                            stack[++top] = node.get_children() + i;
                    }
                }
            }
        }
        /**
         * Test whether two cells are adjacent
         */
        bool is_adjacent(int from, int to) const {
            for (int index : get_cell(from).get_boundaries()) {
                if (get_boundary(index).otherside(from) == to)
                    return true;
            }
            return false;
        }
        
        void get_unique_crn_cfgs(std::vector<const CONFIG*>& result, std::vector<CONFIG>& all_corners) const;
        /**
         * Split the cell into subcubes
         */
        void SplitCell(int node_index) {
            ND::vec<DIM> center;
            double radius = radius_array[get_node(node_index).depth() + 1];
            for (int i = 0; i < NUM_SPLITS; ++i) {
                center = get_node(node_index).center();
                for (int j = 0; j < DIM; ++j) {
                    if (((1 << j) & i) > 0)
                        center[j] += radius;
                    else
                        center[j] -= radius;
                }
                int node_number = get_new_node(center, get_node(node_index).depth() + 1);
                if (i == 0)
                    get_node(node_index).set_children(node_number);
            }
        }
        
        //int InCell(const CONFIG& cfg);
        
        double CalcFreeCellRadius() const {
            double min_dist = obstacle_manager_.dist_to_obsts(robot_);
            double min_time_to_collision = min_dist / robot_.get_max_speed();
            return min_time_to_collision * min_param_speed_;
        }
        
        //void CalcCollisionCellRadius(double& radius, double& dist_to_obstacle) const;
        
        /**
         * Compute the upper bound and the lower bound in each dimension for a boundary
         */
        bool compute_bound(int boundary_index) {
            Boundary<IROBOT>& boundary = all_boundaries[boundary_index];
            Cell<IROBOT>& cell1 = cells[boundary.get_side1()];
            Cell<IROBOT>& cell2 = cells[boundary.get_side2()];
            const ND::vec<DIM>& c1 = cell1.center();
            const ND::vec<DIM>& c2 = cell2.center();
            double r1 = cell1.radius();
            double r2 = cell2.radius();
            
            for (int i = 0; i < DIM; ++i) {
                double max = std::min(c1[i] + r1, c2[i] + r2);
                double min = std::max(c1[i] - r1, c2[i] - r2);
                if (i == boundary.get_fix_dim()) {
                    boundary.set_upper_bound(i, max);
                    boundary.set_lower_bound(i, max);
                    
                } else {
                    boundary.set_upper_bound(i, max);
                    boundary.set_lower_bound(i, min + 1e-14);
                }
            }
            return true;    
        }
        
        /**
         * Compute the fixed dimension for two cells
         * If there is no fixed dimension or there are more than one dimensions, return NOT_FOUND
         */
        int compute_fix_dim(int lhs, int rhs) const {
            const ND::vec<DIM>& c1 = get_cell(lhs).center();
            const ND::vec<DIM>& c2 = get_cell(rhs).center();
            double r1 = get_cell(lhs).radius();
            double r2 = get_cell(rhs).radius();
            int fixdim = NOT_FOUND;
            for (int i = 0; i < DIM; ++i) {
                if (fabs(fabs(c1[i] - c2[i]) - r1 - r2) <= 1e-15) {
                    if (fixdim == NOT_FOUND)
                        fixdim = i;
                    else {
                        // Degenerate case
                        return NOT_FOUND;
                    }
                }
            }
            return fixdim;
        }
        
        int get_new_cell(const ND::vec<DIM>& c, double arg_r) {
            cells.emplace_back(ND::sphere<DIM>(c, arg_r, ND::SPHEREMETRIC::LINFTY));
            if (cells.size() % 1000000 == 0)
                std::cout << "Number of cell is " << cells.size() << '\n';
            return (int)cells.size() - 1;
        }
        
        Cell<IROBOT>& get_cell(int index) {
            return cells[index];
        }
        const Cell<IROBOT>& get_cell(int index) const {
            return cells[index];
        }
        
        int get_new_node(const ND::vec<DIM>& c, int depth) {
            nodes.emplace_back(c, depth);
            if (nodes.size() % 1000000 == 0)
                std::cout << "Number of node is " << nodes.size() << '\n';
            return (int)nodes.size() - 1;
        }

        TreeNode<DIM>& get_node(int index) {
            return nodes[index];
        }
        const TreeNode<DIM>& get_node(int index) const {
            return nodes[index];
        }
        
        Boundary<IROBOT>& get_boundary(int index) {
            return all_boundaries[index];
        }
        const Boundary<IROBOT>& get_boundary(int index) const {
            return all_boundaries[index];
        }
        int get_new_boundary(int from, int to) {
            all_boundaries.emplace_back(from, to);
            return (int)all_boundaries.size() - 1;
        }
        

        std::vector<Cell<IROBOT>> cells;
        std::vector<Boundary<IROBOT>> all_boundaries;
        std::vector<TreeNode<DIM>> nodes;

        IROBOT& robot_;
        robotics::ObstManager& obstacle_manager_;
        
        SamplingMethod method_;
        double radius_array[32];
        int root_index;
        double epsilon_;
        double min_param_speed_;
    };
}
#include "cspace_decomposer.cpp"
#endif /* CSPACE_DECOMPOSER_H */
