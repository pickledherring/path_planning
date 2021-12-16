/*
inspired by:
https://www.frontiersin.org/articles/10.3389/fnbot.2019.00015/full
and
https://openreview.net/forum?id=HyZ1CJZ_-r

*/
#include <random>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <cstring>
#include <semaphore.h>
#include <math.h>
using namespace std;

/*
Step 1: set up grid, start, target, and taboo.
build square grid, randomly place high weights
start and target are in the middle of the leftmost and rightmost columns.
*/
const int GRID_DIM = 16; // size of GRID_DIM x GRID_DIM grid
const float OBST_CHANCE = 4; // chance a square will have a penalty added out of 10
random_device dev;
mt19937 rng(dev());
uniform_int_distribution<mt19937::result_type> coords(0, GRID_DIM);
int start[2] = {GRID_DIM / 2, 0};
int target[2] = {GRID_DIM / 2, GRID_DIM - 1};
int grid[GRID_DIM][GRID_DIM] = {{0}};
float pher[GRID_DIM][GRID_DIM] = {{.1}};

// sem_t * sem_main;

const float alpha = 1.0;
const float beta_v = 1.5;
const float evap = .25;
const float Q = 100.0;
int n_ants = 10;
int epochs = 10;

void penalty_add() {
    int penalty_size = 8; // about half the grid width should be enough
    uniform_real_distribution<> penalty(0, 10.0);

    for (int i = 0; i < GRID_DIM; i++) {
        for (int j = 0; j < GRID_DIM; j++) {
            float pen_chance = penalty(rng);
            if (pen_chance < OBST_CHANCE) {
                grid[i][j] = grid[i][j] + penalty_size;
            }
        }
    }
}

// choice == 0 for grid, 1 for pher, 2 for pher/grid
void grid_print(int choice = 0) {
    switch(choice) {
        case 0:
            for (int i = 0; i < GRID_DIM; i++) {
                for (int j = 0; j < GRID_DIM; j++) {
                    if (grid[i][j] > 9) 
                        cout<<" "<<grid[i][j]<<" ";
                    else
                        cout<<" "<<grid[i][j]<<"  ";
                }
                cout<<endl;
            }
            cout<<endl;
            break;

        case 1:
            for (int i = 0; i < GRID_DIM; i++) {
                for (int j = 0; j < GRID_DIM; j++) {
                    if (round(pher[i][j]) > 9) 
                        cout<<" "<<round(pher[i][j] * 100.0) / 100<<" ";
                    else
                        cout<<" "<<round(pher[i][j] * 100.0) / 100<<"  ";
                }
                cout<<endl;
            }
            cout<<endl;
            break;

        case 2:
            for (int i = 0; i < GRID_DIM; i++) {
                for (int j = 0; j < GRID_DIM; j++) {
                    if (grid[i][j] > 9) 
                        cout<<" "<<round(pher[i][j] * 100.0) / 100<<"/"<<grid[i][j]<<" ";
                    else
                        cout<<" "<<round(pher[i][j] * 100.0) / 100<<"/"<<grid[i][j]<<"  ";
                }
                cout<<endl;
            }
            cout<<endl;
            break;
    }
}

/*
Step 5: Calculate probability for each grid, and select using rand_selector.
choose from probability dist. given by
( (p_ij)^alpha * (1 / d_ij)^beta ) / sum_all( (p_ij)^a * (1 / d_ij)^b )
*/
int* choose_square(int* pos, float rand_selector, bool taboo[][GRID_DIM]) {
    float total = 0;
    vector<float> numer_vec;
    vector<vector<int>> adj_pos;
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            int x = pos[0] + i;
            int y = pos[1] + j;
            if (x == target[0] && y == target[1]) {
                int* new_pos = new int[2];
                new_pos[0] = x;
                new_pos[1] = y;
                return new_pos;
            }
            else if ((i == 0 && j == 0) || x < 0 || y < 0 || x >= GRID_DIM || y >= GRID_DIM ||
            taboo[x][y])
                continue;
            else {
                adj_pos.push_back({x, y});
                float dist = sqrt(i * i + j * j);
                float numer = pow(pher[x][y], alpha) * pow(1 / dist, beta_v);
                numer_vec.push_back(numer);
                total += numer;
                // cout<<"x and y: "<<x<<", "<<y<<"\t";
                // cout<<"numer: "<<numer<<"\t";
                // cout<<"cumul: "<<total<<endl;
            }
        }
    }
    // for (int i = 0; i < numer_vec.size(); i++)
    //     cout<<numer_vec[i]<<"\t";
    // cout<<endl<<total<<endl;
    // cout<<"elements: "<<endl;
    int* new_pos = new int[2];
    float cumulative = 0;
    for (int i = 0; i < numer_vec.size(); i++) {
        float frac = numer_vec[i] / total;
        numer_vec[i] = cumulative + frac;
        if (rand_selector <= numer_vec[i]) {
            // cout<<rand_selector<<" <= "<<numer_vec[i]<<"?\t"<<(rand_selector <= numer_vec[i])<<endl;
            new_pos[0] = adj_pos[i][0];
            new_pos[1] = adj_pos[i][1];
            break;
        }
        cumulative += frac;
    }
    // cout<<"new_pos: "<<new_pos[0]<<", "<<new_pos[1]<<endl;
    // cout<<"cumulative: "<<cumulative<<endl;

    // int* new_pos = new int[2];
    // new_pos[0] = pos[0];
    // new_pos[1] = pos[1] + 1;
    return new_pos;
}

/*
Step 4: Process deadlock. According to the taboo table, judge whether ants are trapped.
If in deadlock, retract and add the deadlock node to the taboo table
*/
bool deadlock(int* pos, bool taboo[][GRID_DIM]) {
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            int x = pos[0] + i;
            int y = pos[1] + j;
            if (x >= 0 && y >= 0 && x < GRID_DIM && y < GRID_DIM && !taboo[x][y]) {
                return false;
            }
            else
                continue;
        }
    }
    return true;
}

int* retract(vector<int*> &path) {
    int* new_pos = new int[2];
    path.pop_back();
    // cout<<"path ";
    // for (int i = 0; i < path.size(); i++) {
    //     cout<<i<<":  "<<path[i][0]<<", "<<path[i][1]<<"\t    ";
    // }
    // cout<<endl;
    new_pos[0] = path.back()[0];
    new_pos[1] = path.back()[1];
    // cout<<"retract to "<<new_pos[0]<<", "<<new_pos[1]<<endl;
    return new_pos;
}

/*
Step 3: Place the ant k (k = 1, 2, ⋯ , m) on start. Proceed to have them randomly search
for the target, and return a path.
*/
vector<int*> release_ant(int ant) {
    bool taboo[GRID_DIM][GRID_DIM] = {{false}};
    taboo[start[0]][start[1]] = true;
    vector<int*> path;
    path.push_back(start);
    int pos[2];
    memcpy(pos, start, 2 * sizeof(int));
    uniform_real_distribution<> prob(0, 1.0);
    // int deadlock_count = 0;
    while (pos[0] != target[0] || pos[1] != target[1] /* && deadlock_count < 50*/) {
        float rand_selector = prob(rng);
        int* move_pos = choose_square(pos, rand_selector, taboo);
        // cout<<"chose "<<move_pos[0]<<", "<<move_pos[1]<<endl;
        taboo[move_pos[0]][move_pos[1]] = true;
        path.push_back(move_pos);
        while (deadlock(move_pos, taboo) /* && deadlock_count < 50*/) {
            // cout<<"deadlock at "<<move_pos[0]<<", "<<move_pos[1]<<endl;
            // cout<<"path ";
            // for (int i = 0; i < path.size(); i++) {
            //     cout<<i<<":  "<<path[i][0]<<", "<<path[i][1]<<"\t    ";
            // }
            // cout<<endl;
            taboo[move_pos[0]][move_pos[1]] = true;
            move_pos = retract(path);
            // deadlock_count++;
        }
        memcpy(pos, move_pos, 2 * sizeof(int));
    }

    // sem_post(sem_main);
    // for (int i = 0; i < GRID_DIM; i++) {
    //     for (int j = 0; j < GRID_DIM; j++) {
    //         if (taboo[i][j]) 
    //             cout<<" T  ";
    //         else
    //             cout<<" F  ";
    //     }
    //     cout<<endl;
    // }

    return path;
}

/*
Step 7: Update pheromone. pheromone update occurs by (1 - evap)ph_ij + q / len_k
where 0 < evap < 1 and q is a constant, and len_k is the cost of the path ant k traveled.
Only update nodes along the path.
*/
void pher_update(vector<int*> const &path) {
    int path_len = 0;
    for (int* v : path) {
        path_len += grid[v[0]][v[1]];
    }
    float pheromone = Q / path_len;
    for (int* v : path) {
        pher[v[0]][v[1]] *= (1 - evap);
        pher[v[0]][v[1]] += pheromone;
    }
    cout<<"path len "<<path_len<<"\tpheromone "<<pheromone<<endl;
    // grid_print(2);
}


int main(int argc, char** argv) {
    penalty_add();
    penalty_add();
    // grid_print();
    for (auto &arr : pher)
        fill(begin(arr), end(arr), 0.1);
    
    for (int i = 0; i < epochs; i++) {
        vector<int*> path = release_ant(0);

        /*
        After each iteration, if the number of iterations satisfies inequality N ≤ Nmax,
        update the pheromone grid and determine whether it meets the convergence conditions.
        If the number of iterations satisfies inequality N > Nmax, stop.
        */
        pher_update(path);
    }
    // print path
    // for (int i = 0; i < path.size(); i++) {
    //     cout<<i<<":  "<<path[i][0]<<", "<<path[i][1]<<"\t    ";
    // }
    // cout<<endl;

    // pthread_t *threads;
    // threads = (pthread_t*)malloc(n_ants * sizeof(pthread_t));
    // sem_open("/sem_main", NULL, NULL, n_ants);
    // for (int ant = 0; ant < n_ants; ants++) {
    //     pthread_create(&threads[i], NULL, release_ant, ant);
    //     release_ant(ant);
    // }
    // sem_wait(sem_main);
    // free(threads);
    return 0;
}