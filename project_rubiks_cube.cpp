#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <sstream>
#include <climits>
#include <cstring>

using namespace std;

// Constants defining the cube structure
const int CUBE_SIZE = 3;        // 3x3x3 cube
const int NUM_FACES = 6;        // Six faces on a cube
const int MAX_SEARCH_DEPTH = 25; // Maximum search depth to prevent infinite loops

// Enumeration for the six faces of the cube
// These match the standard Rubik's cube notation
enum CubeFace
{
    UP = 0,     // U - Top face (usually white)
    FRONT = 1,  // F - Front face (usually red)
    LEFT = 2,   // L - Left face (usually blue)
    BACK = 3,   // B - Back face (usually orange)
    RIGHT = 4,  // R - Right face (usually green)
    DOWN = 5    // D - Bottom face (usually yellow)
};

// Direction of rotation for face turns
enum TurnDirection
{
    CLOCKWISE = 1,          // 90 degrees clockwise
    COUNTER_CLOCKWISE = -1  // 90 degrees counterclockwise
};

// Represents a single move on the cube (like R, U', F2, etc.)
class CubeMove
{
private:
    CubeFace face;        // Which face to turn (U, F, L, B, R, D)
    TurnDirection dir;    // Direction of the turn
    int id;               // Unique identifier for this move type (1-12)

public:
    // Constructor: creates a move with the specified face, direction, and ID
    CubeMove(CubeFace f = UP, TurnDirection d = CLOCKWISE, int i = 0)
        : face(f), dir(d), id(i) {}
    
    // Getters for accessing move properties
    CubeFace getFace() const { return face; }
    TurnDirection getDir() const { return dir; }
    int getId() const { return id; }
    
    // Converts the move to standard Rubik's cube notation
    // For example: R (right clockwise), R' (right counterclockwise)
    string toString() const
    {
        static string names[] = {"U", "F", "L", "B", "R", "D"};
        return names[face] + (dir == COUNTER_CLOCKWISE ? "'" : "");
    }
};

// Represents the current state of the entire cube
class CubeState
{
private:
    // 3D array storing the color/number of each sticker
    // c[face][row][col] = the sticker at that position
    int c[NUM_FACES][CUBE_SIZE][CUBE_SIZE];

public:
    // Default constructor: initializes all stickers to 0
    CubeState() { memset(c, 0, sizeof(c)); }

    // Copy constructor: creates a deep copy of another cube state
    CubeState(const CubeState &other) { memcpy(c, other.c, sizeof(c)); }

    // Assignment operator: allows cube1 = cube2
    CubeState &operator=(const CubeState &other)
    {
        if (this != &other)
            memcpy(c, other.c, sizeof(c));
        return *this;
    }

    // Access operator: allows cube(face, row, col) = value
    int &operator()(int face, int row, int col) { return c[face][row][col]; }

    // Const access operator: allows reading from const cube states
    const int &operator()(int face, int row, int col) const { return c[face][row][col]; }

    // Equality operator: checks if two cube states are identical
    bool operator==(const CubeState &other) const
    {
        return memcmp(c, other.c, sizeof(c)) == 0;
    }

    // Converts the cube state to a string for hashing and duplicate detection
    string hashString() const
    {
        string result;
        for (int face = 0; face < NUM_FACES; face++)
            for (int row = 0; row < CUBE_SIZE; row++)
                for (int col = 0; col < CUBE_SIZE; col++)
                    result += to_string(c[face][row][col]) + ",";
        return result;
    }
};

// Generates all possible moves for the cube (12 total: 6 faces × 2 directions)
class MoveGenerator
{
private:
    static vector<CubeMove> moves;
    static bool initialized;
    
    // Initialize the list of all possible moves
    static void initialize()
    {
        if (initialized) return;
        
        moves = {
            // Left face moves
            CubeMove(LEFT, CLOCKWISE, 1), CubeMove(LEFT, COUNTER_CLOCKWISE, 2),
            // Right face moves
            CubeMove(RIGHT, CLOCKWISE, 3), CubeMove(RIGHT, COUNTER_CLOCKWISE, 4),
            // Up face moves
            CubeMove(UP, CLOCKWISE, 5), CubeMove(UP, COUNTER_CLOCKWISE, 6),
            // Down face moves
            CubeMove(DOWN, CLOCKWISE, 7), CubeMove(DOWN, COUNTER_CLOCKWISE, 8),
            // Front face moves
            CubeMove(FRONT, CLOCKWISE, 9), CubeMove(FRONT, COUNTER_CLOCKWISE, 10),
            // Back face moves
            CubeMove(BACK, CLOCKWISE, 11), CubeMove(BACK, COUNTER_CLOCKWISE, 12)
        };
        initialized = true;
    }

public:
    // Returns all possible moves
    static const vector<CubeMove> &getAllMoves()
    {
        initialize();
        return moves;
    }
};

// Static member definitions
vector<CubeMove> MoveGenerator::moves;
bool MoveGenerator::initialized = false;

// Handles all cube rotations: face rotations and adjacent edge movements
class CubeRotator
{
private:
    // Rotates a face 90 degrees clockwise
    // This affects the 3x3 grid of stickers on the face itself
    static void rotateFaceClockwise(CubeState &state, CubeFace face)
    {
        // Save the middle edges and corners during rotation
        int temp = state(face, 0, 1);
        state(face, 0, 1) = state(face, 1, 0);
        state(face, 1, 0) = state(face, 2, 1);
        state(face, 2, 1) = state(face, 1, 2);
        state(face, 1, 2) = temp;
        
        // Rotate the corners
        temp = state(face, 0, 0);
        state(face, 0, 0) = state(face, 2, 0);
        state(face, 2, 0) = state(face, 2, 2);
        state(face, 2, 2) = state(face, 0, 2);
        state(face, 0, 2) = temp;
    }
    
    // Rotates a face 90 degrees counterclockwise
    static void rotateFaceCounterClockwise(CubeState &state, CubeFace face)
    {
        // Save the middle edges during rotation
        int temp = state(face, 0, 1);
        state(face, 0, 1) = state(face, 1, 2);
        state(face, 1, 2) = state(face, 2, 1);
        state(face, 2, 1) = state(face, 1, 0);
        state(face, 1, 0) = temp;
        
        // Rotate the corners
        temp = state(face, 0, 0);
        state(face, 0, 0) = state(face, 0, 2);
        state(face, 0, 2) = state(face, 2, 2);
        state(face, 2, 2) = state(face, 2, 0);
        state(face, 2, 0) = temp;
    }
    
    // Rotates the edges adjacent to the RIGHT face
    static void rotateAdjacentRight(CubeState &state, TurnDirection direction)
    {
        for (int i = 0; i < 3; i++)
        {
            int temp = state(UP, i, 2);
            if (direction == CLOCKWISE)
            {
                // Clockwise: UP -> FRONT -> DOWN -> BACK -> UP
                state(UP, i, 2) = state(FRONT, i, 2);
                state(FRONT, i, 2) = state(DOWN, i, 2);
                state(DOWN, i, 2) = state(BACK, i, 2);
                state(BACK, i, 2) = temp;
            }
            else
            {
                // Counterclockwise: UP -> BACK -> DOWN -> FRONT -> UP
                state(UP, i, 2) = state(BACK, i, 2);
                state(BACK, i, 2) = state(DOWN, i, 2);
                state(DOWN, i, 2) = state(FRONT, i, 2);
                state(FRONT, i, 2) = temp;
            }
        }
    }
    
    // Rotates the edges adjacent to the LEFT face
    static void rotateAdjacentLeft(CubeState &state, TurnDirection direction)
    {
        for (int i = 0; i < 3; i++)
        {
            int temp = state(UP, i, 0);
            if (direction == CLOCKWISE)
            {
                // Clockwise: UP -> BACK -> DOWN -> FRONT -> UP
                state(UP, i, 0) = state(BACK, i, 0);
                state(BACK, i, 0) = state(DOWN, i, 0);
                state(DOWN, i, 0) = state(FRONT, i, 0);
                state(FRONT, i, 0) = temp;
            }
            else
            {
                // Counterclockwise: UP -> FRONT -> DOWN -> BACK -> UP
                state(UP, i, 0) = state(FRONT, i, 0);
                state(FRONT, i, 0) = state(DOWN, i, 0);
                state(DOWN, i, 0) = state(BACK, i, 0);
                state(BACK, i, 0) = temp;
            }
        }
    }
    
    // Rotates the edges adjacent to the UP face
    static void rotateAdjacentUp(CubeState &state, TurnDirection direction)
    {
        for (int i = 0; i < 3; i++)
        {
            int temp = state(FRONT, 0, i);
            if (direction == CLOCKWISE)
            {
                // Clockwise: FRONT -> RIGHT -> BACK -> LEFT -> FRONT
                state(FRONT, 0, i) = state(RIGHT, 0, i);
                state(RIGHT, 0, i) = state(BACK, 0, i);
                state(BACK, 0, i) = state(LEFT, 0, i);
                state(LEFT, 0, i) = temp;
            }
            else
            {
                // Counterclockwise: FRONT -> LEFT -> BACK -> RIGHT -> FRONT
                state(FRONT, 0, i) = state(LEFT, 0, i);
                state(LEFT, 0, i) = state(BACK, 0, i);
                state(BACK, 0, i) = state(RIGHT, 0, i);
                state(RIGHT, 0, i) = temp;
            }
        }
    }
    
    // Rotates the edges adjacent to the DOWN face
    static void rotateAdjacentDown(CubeState &state, TurnDirection direction)
    {
        for (int i = 0; i < 3; i++)
        {
            int temp = state(FRONT, 2, i);
            if (direction == CLOCKWISE)
            {
                // Clockwise: FRONT -> LEFT -> BACK -> RIGHT -> FRONT
                state(FRONT, 2, i) = state(LEFT, 2, i);
                state(LEFT, 2, i) = state(BACK, 2, i);
                state(BACK, 2, i) = state(RIGHT, 2, i);
                state(RIGHT, 2, i) = temp;
            }
            else
            {
                // Counterclockwise: FRONT -> RIGHT -> BACK -> LEFT -> FRONT
                state(FRONT, 2, i) = state(RIGHT, 2, i);
                state(RIGHT, 2, i) = state(BACK, 2, i);
                state(BACK, 2, i) = state(LEFT, 2, i);
                state(LEFT, 2, i) = temp;
            }
        }
    }
    
    // Rotates the edges adjacent to the FRONT face
    static void rotateAdjacentFront(CubeState &state, TurnDirection direction)
    {
        for (int i = 0; i < 3; i++)
        {
            int temp = state(UP, 2, i);
            if (direction == CLOCKWISE)
            {
                // Clockwise: UP -> LEFT -> DOWN -> RIGHT -> UP
                state(UP, 2, i) = state(LEFT, 2 - i, 2);
                state(LEFT, 2 - i, 2) = state(DOWN, 0, 2 - i);
                state(DOWN, 0, 2 - i) = state(RIGHT, i, 0);
                state(RIGHT, i, 0) = temp;
            }
            else
            {
                // Counterclockwise: UP -> RIGHT -> DOWN -> LEFT -> UP
                state(UP, 2, i) = state(RIGHT, i, 0);
                state(RIGHT, i, 0) = state(DOWN, 0, 2 - i);
                state(DOWN, 0, 2 - i) = state(LEFT, 2 - i, 2);
                state(LEFT, 2 - i, 2) = temp;
            }
        }
    }
    
    // Rotates the edges adjacent to the BACK face
    static void rotateAdjacentBack(CubeState &state, TurnDirection direction)
    {
        for (int i = 0; i < 3; i++)
        {
            int temp = state(UP, 0, i);
            if (direction == CLOCKWISE)
            {
                // Clockwise: UP -> RIGHT -> DOWN -> LEFT -> UP
                state(UP, 0, i) = state(RIGHT, i, 2);
                state(RIGHT, i, 2) = state(DOWN, 2, 2 - i);
                state(DOWN, 2, 2 - i) = state(LEFT, 2 - i, 0);
                state(LEFT, 2 - i, 0) = temp;
            }
            else
            {
                // Counterclockwise: UP -> LEFT -> DOWN -> RIGHT -> UP
                state(UP, 0, i) = state(LEFT, 2 - i, 0);
                state(LEFT, 2 - i, 0) = state(DOWN, 2, 2 - i);
                state(DOWN, 2, 2 - i) = state(RIGHT, i, 2);
                state(RIGHT, i, 2) = temp;
            }
        }
    }

public:
    // Applies a move to a cube state and returns the new state
    static CubeState applyMove(const CubeState &originalState, const CubeMove &move)
    {
        CubeState newState = originalState;
        
        // First, rotate the face itself
        if (move.getDir() == CLOCKWISE)
            rotateFaceClockwise(newState, move.getFace());
        else
            rotateFaceCounterClockwise(newState, move.getFace());
        
        // Then, rotate the adjacent edges based on which face was turned
        switch (move.getFace())
        {
            case RIGHT:
                rotateAdjacentRight(newState, move.getDir());
                break;
            case LEFT:
                rotateAdjacentLeft(newState, move.getDir());
                break;
            case UP:
                rotateAdjacentUp(newState, move.getDir());
                break;
            case DOWN:
                rotateAdjacentDown(newState, move.getDir());
                break;
            case FRONT:
                rotateAdjacentFront(newState, move.getDir());
                break;
            case BACK:
                rotateAdjacentBack(newState, move.getDir());
                break;
        }
        
        return newState;
    }
};

// Improved heuristic function for better solving performance
class ImprovedHeuristic
{
public:
    // Calculates an improved Manhattan distance heuristic
    // This estimates the minimum number of moves needed to solve the cube
    int estimateMovesToSolve(const CubeState &currentState, const CubeState &goalState) const
    {
        int totalMisplacedStickers = 0;
        
        // Count all misplaced stickers across all faces
        for (int face = 0; face < NUM_FACES; face++)
        {
            for (int row = 0; row < CUBE_SIZE; row++)
            {
                for (int col = 0; col < CUBE_SIZE; col++)
                {
                    if (currentState(face, row, col) != goalState(face, row, col))
                    {
                        totalMisplacedStickers++;
                    }
                }
            }
        }
        
        // Each move can fix at most 8 stickers (but usually fewer)
        // We use a more conservative estimate to maintain admissibility
        return (totalMisplacedStickers + 7) / 8;  // Round up division
    }
};

// Represents a node in the search tree
class SearchNode
{
public:
    CubeState state;           // Current cube configuration
    SearchNode *parent;        // Parent node in the search tree
    CubeMove moveFromParent;   // Move that led to this state
    int depthFromStart;        // How many moves from the initial state
    int estimatedTotalCost;    // f(n) = g(n) + h(n)
    
    SearchNode(const CubeState &cubeState, SearchNode *parentNode, 
               CubeMove move, int depth)
        : state(cubeState), parent(parentNode), moveFromParent(move), 
          depthFromStart(depth), estimatedTotalCost(0) {}
    
    // Set the heuristic value and calculate total estimated cost
    void setHeuristicValue(int heuristicValue) 
    { 
        estimatedTotalCost = depthFromStart + heuristicValue; 
    }
    
    // Reconstructs the solution path from start to this node
    vector<CubeMove> getSolutionPath() const
    {
        vector<CubeMove> path;
        const SearchNode *currentNode = this;
        
        // Traverse back to the root, collecting moves
        while (currentNode->parent != nullptr)
        {
            path.push_back(currentNode->moveFromParent);
            currentNode = currentNode->parent;
        }
        
        // Reverse to get the path from start to goal
        reverse(path.begin(), path.end());
        return path;
    }
};

// IDA* (Iterative Deepening A*) solver implementation
class IDASolver
{
private:
    ImprovedHeuristic heuristic;
    unordered_set<string> visitedStates;
    int nodesExplored;
    int maxSearchDepth;

public:
    IDASolver(int maxDepth = MAX_SEARCH_DEPTH) 
        : nodesExplored(0), maxSearchDepth(maxDepth) {}
    
    // Main solving method using IDA*
    vector<CubeMove> solveCube(const CubeState &startState, const CubeState &goalState)
    {
        // Check if already solved
        if (startState == goalState)
            return vector<CubeMove>();
        
        // Start with the heuristic value of the initial state
        int currentThreshold = heuristic.estimateMovesToSolve(startState, goalState);
        
        // Iteratively deepen the search
        while (currentThreshold <= maxSearchDepth)
        {
            visitedStates.clear();
            auto result = depthLimitedSearch(startState, goalState, 0, currentThreshold, nullptr);
            
            // If we found a solution, return it
            if (!result.first.empty())
                return result.first;
            
            // If no more nodes to explore, we're done
            if (result.second == INT_MAX)
                break;
            
            // Increase threshold for next iteration
            currentThreshold = result.second;
        }
        
        // No solution found within the depth limit
        return vector<CubeMove>();
    }
    
    // Returns the number of nodes explored during the search
    int getNodesExplored() const { return nodesExplored; }

private:
    // Recursive depth-limited search with A* pruning
    pair<vector<CubeMove>, int> depthLimitedSearch(const CubeState &currentState, 
                                                   const CubeState &goalState, 
                                                   int currentDepth, 
                                                   int depthThreshold, 
                                                   SearchNode *parentNode)
    {
        nodesExplored++;
        
        // Calculate heuristic and total estimated cost
        int heuristicValue = heuristic.estimateMovesToSolve(currentState, goalState);
        int totalEstimatedCost = currentDepth + heuristicValue;
        
        // If we exceed the threshold, prune this branch
        if (totalEstimatedCost > depthThreshold)
            return {vector<CubeMove>(), totalEstimatedCost};
        
        // Check if we've reached the goal
        if (currentState == goalState)
            return {vector<CubeMove>(), -1};  // -1 indicates success
        
        // Avoid revisiting states (cycle detection)
        string stateKey = currentState.hashString();
        if (visitedStates.count(stateKey))
            return {vector<CubeMove>(), INT_MAX};
        
        visitedStates.insert(stateKey);
        
        // Try all possible moves
        int minNextThreshold = INT_MAX;
        for (const auto &move : MoveGenerator::getAllMoves())
        {
            CubeState nextState = CubeRotator::applyMove(currentState, move);
            auto result = depthLimitedSearch(nextState, goalState, currentDepth + 1, 
                                           depthThreshold, parentNode);
            
            // If we found a solution, prepend this move and return
            if (result.second == -1)
            {
                vector<CubeMove> solution = {move};
                solution.insert(solution.end(), result.first.begin(), result.first.end());
                return {solution, -1};
            }
            
            // Track the minimum threshold for the next iteration
            if (result.second < minNextThreshold)
                minNextThreshold = result.second;
        }
        
        // Remove from visited set to allow revisiting from other paths
        visitedStates.erase(stateKey);
        
        return {vector<CubeMove>(), minNextThreshold};
    }
};

// Utility functions for input/output
void readCubeStates(const string &filename, CubeState &startState, CubeState &goalState)
{
    ifstream inputFile(filename);
    if (!inputFile)
    {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }
    
    // Read the initial cube state
    for (int face = 0; face < NUM_FACES; face++)
    {
        for (int row = 0; row < CUBE_SIZE; row++)
        {
            for (int col = 0; col < CUBE_SIZE; col++)
            {
                inputFile >> startState(face, row, col);
            }
        }
    }
    
    // Read the goal cube state
    for (int face = 0; face < NUM_FACES; face++)
    {
        for (int row = 0; row < CUBE_SIZE; row++)
        {
            for (int col = 0; col < CUBE_SIZE; col++)
            {
                inputFile >> goalState(face, row, col);
            }
        }
    }
    
    inputFile.close();
}

// Main function - entry point of the program
int main()
{
    try
    {
        // Read the initial and goal cube states from the input file
        CubeState startState, goalState;
        readCubeStates("input.txt", startState, goalState);
        
        // Create the solver and attempt to solve the cube
        IDASolver solver;
        cout << "Solving the cube... This may take a while for complex scrambles." << endl;
        
        auto solution = solver.solveCube(startState, goalState);
        
        // Output the results
        if (solution.empty())
        {
            if (startState == goalState)
            {
                cout << "The cube is already solved!" << endl;
                cout << "0" << endl;
            }
            else
            {
                cout << "No solution found within the search depth limit." << endl;
                cout << "Try increasing MAX_SEARCH_DEPTH or check if the cube configuration is valid." << endl;
            }
        }
        else
        {
            cout << "Solution found!" << endl;
            cout << "Number of moves: " << solution.size() << endl;
            cout << "Move sequence: ";
            for (const auto &move : solution)
            {
                cout << move.toString() << " ";
            }
            cout << endl;
            
            cout << "Nodes explored: " << solver.getNodesExplored() << endl;
        }
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
