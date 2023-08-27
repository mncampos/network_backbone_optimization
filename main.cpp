#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>
#include <random>
#include <thread>
#include <chrono>

using namespace std;

// Define the initial solution heuristic
const std::string heuristic = "greedy";

// Define the number of threads to use
const int NUM_THREADS = 6;

// Define filename
const std::string file = "instance_10000_19800.dat";

// Define the maximum number of iterations without improvement
const int MAX_ITERATIONS = 5;

// Define the size of the tabu list
const int TABU_SIZE = 25;

// Define the number of neighbors to be generated in each iteration
const int NUM_NEIGHBORS = 5;

// Define the random number generator
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

// Function to generate a random integer in the range [a, b]
int randint(int a, int b) {
    uniform_int_distribution<int> dist(a, b);
    return dist(rng);
}

// Function to check if a solution is feasible
bool is_feasible(const vector<vector<int>>& graph, const vector<bool>& solution) {
    // Check if the solution is connected
    set<int> visited;
    for (int i = 0; i < solution.size(); i++) {
        if (solution[i]) {
            visited.insert(i);
            break;
        }
    }
    if (visited.empty()) return false;
    set<int> to_visit;
    for (int neighbor : graph[*visited.begin()]) {
        if (solution[neighbor]) to_visit.insert(neighbor);
    }
    while (!to_visit.empty()) {
        int current = *to_visit.begin();
        to_visit.erase(current);
        visited.insert(current);
        for (int neighbor : graph[current]) {
            if (solution[neighbor] && visited.find(neighbor) == visited.end()) {
                to_visit.insert(neighbor);
            }
        }
    }
    for (int i = 0; i < solution.size(); i++) {
        if (solution[i] && visited.find(i) == visited.end()) return false;
    }

    // Check if all vertices are in the solution or have a neighbor in the solution
    for (int i = 0; i < solution.size(); i++) {
        if (!solution[i]) {
            bool has_neighbor = false;
            for (int neighbor : graph[i]) {
                if (solution[neighbor]) {
                    has_neighbor = true;
                    break;
                }
            }
            if (!has_neighbor) return false;
        }
    }

    return true;
}



// Function to generate a random feasible solution
vector<bool> generate_random_solution(const vector<vector<int>>& graph) {
    int n = graph.size();
    vector<bool> solution(n, false);

    while (!is_feasible(graph, solution)) {
        for (int i = 0; i < n; i++) {
            solution[i] = randint(0, 1);
        }
    }

    return solution;
}

vector<bool> generate_greedy_solution(const vector<vector<int>>& graph) {
    int n = graph.size();
    vector<bool> solution(n, false);


    // Choose the first vertex randomly
    int start_vertex = randint(0, n - 1);
    solution[start_vertex] = true;

    // Keep track of visited vertices
    vector<bool> visited(n, false);
    visited[start_vertex] = true;

    // Greedily select vertices based on the number of neighbors already selected
    for (int i = 1; i < n; i++) {
        int best_vertex = -1;
        int best_neighbors = -1;

        // Find the vertex with the most neighbors among visited vertices
        for (int vertex = 0; vertex < n; vertex++) {
            if (!visited[vertex]) {
                int neighbors_count = 0;
                for (int neighbor : graph[vertex]) {
                    if (solution[neighbor]) {
                        neighbors_count++;
                    }
                }
                if (neighbors_count > best_neighbors) {
                    best_neighbors = neighbors_count;
                    best_vertex = vertex;
                }
            }
        }

        // Mark the best vertex as selected
        solution[best_vertex] = true;
        visited[best_vertex] = true;
    }


    return solution;
}

// Function to calculate the cost of a solution
int calculate_cost(const vector<bool>& solution) {
    return count(solution.begin(), solution.end(), true);
}

// Function to generate a random neighbor of a solution
vector<bool> generate_random_neighbor(const vector<vector<int>>& graph, const vector<bool>& solution) {
    int n = graph.size();
    vector<bool> neighbor = solution;
    while (true) {
        int index = randint(0, n - 1);
        neighbor[index] = !neighbor[index];
        if (is_feasible(graph, neighbor)) break;
        neighbor[index] = !neighbor[index];
    }
    return neighbor;
}

vector<vector<bool>> generate_neighbors_parallel(const vector<vector<int>>& graph, const vector<bool>& solution) {
    vector<vector<bool>> neighbors(NUM_NEIGHBORS);


    // Function to generate neighbors in each thread
    auto generate_neighbors_thread = [&](int thread_id) {
        for (int i = thread_id; i < NUM_NEIGHBORS; i += NUM_THREADS) {
            neighbors[i] = generate_random_neighbor(graph, solution);
        }
        };

    // Create and launch threads
    vector<thread> threads;
    for (int i = 0; i < NUM_THREADS; i++) {
        threads.emplace_back(generate_neighbors_thread, i);
    }

    // Wait for threads to finish
    for (auto& thread : threads) {
        thread.join();
    }

    return neighbors;
}

// Function to perform the Tabu Search algorithm
vector<bool> tabu_search(const vector<vector<int>>& graph) {
    // Generate an initial random feasible solution
    vector<bool> best_solution;

    if (heuristic == "greedy") {
        best_solution = generate_greedy_solution(graph);
    }

    if (heuristic == "random") {
        best_solution = generate_random_solution(graph);
    }

    int best_cost = calculate_cost(best_solution);

    // Initialize the current solution and cost
    vector<bool> current_solution = best_solution;
    int current_cost = best_cost;

    // Initialize the tabu list
    set<vector<bool>> tabu_list;

    // Initialize the number of iterations without improvement
    int iterations_without_improvement = 0;

    // Perform the Tabu Search algorithm
    while (iterations_without_improvement < MAX_ITERATIONS) {
        // Generate neighbors of the current solution
      vector<vector<bool>> neighbors = generate_neighbors_parallel(graph, current_solution);
        // Find the best non-tabu neighbor
        bool found_non_tabu_neighbor = false;
        for (const auto& neighbor : neighbors) {
            if (tabu_list.find(neighbor) == tabu_list.end()) {
                int cost = calculate_cost(neighbor);
                if (!found_non_tabu_neighbor || cost < current_cost) {
                    current_solution = neighbor;
                    current_cost = cost;
                    found_non_tabu_neighbor = true;
                }
            }
        }

        // If no non-tabu neighbor was found, find the best neighbor regardless of tabu status
        if (!found_non_tabu_neighbor) {
            for (const auto& neighbor : neighbors) {
                int cost = calculate_cost(neighbor);
                if (cost < current_cost) {
                    current_solution = neighbor;
                    current_cost = cost;
                    break;
                }
            }
        }

        // Update the best solution and cost
        if (current_cost < best_cost) {
            best_solution = current_solution;
            best_cost = current_cost;
            iterations_without_improvement = 0;
        }
        else {
            iterations_without_improvement++;
        }

        // Update the tabu list
        tabu_list.insert(current_solution);
        if (tabu_list.size() > TABU_SIZE) {
            tabu_list.erase(tabu_list.begin());
        }
    }

    return best_solution;
}


int main() {
    // Open the input file
    ifstream input_file(file);

    // Read the number of vertices and edges
    int n, m;
    input_file >> n >> m;

    // Initialize the graph
    vector<vector<int>> graph(n);

    // Read the edges
    for (int i = 0; i < m; i++) {
        int u, v;
        input_file >> u >> v;
        u--; v--;
        graph[u].push_back(v);
        graph[v].push_back(u);
    }

    // Close the input file
    input_file.close();

    // Perform the Tabu Search algorithm
    auto start_time = chrono::high_resolution_clock::now();
    std::cout << "Calculating best solution..." << std::endl;
    vector<bool> solution = tabu_search(graph);
    // Record the end time
    auto end_time = chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);

    // Output the solution
    int backbone_size = 0;
    cout << "Solution: ";
    for (int i = 0; i < n; i++) {
        if (solution[i]) {
            cout << i + 1 << " ";
            backbone_size++;
        }
    }
    cout << endl;
    cout << "Backbone size: " << backbone_size << endl;
    cout << "Time taken: " << duration.count() << " milliseconds" << endl;


    return 0;
}
