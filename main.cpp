#define _GLIBCXX_DEBUG

#include<bits/stdc++.h>
#include <utility>
#include <algorithm>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <chrono>

typedef long long ll;
typedef long double ld;
typedef std::vector<int> vi;
typedef std::vector<double> vd;

int const MAX = 1401;
double const eps = 2e-7;
int const MOD = 179179179;
int const pppower = 179;
int const INFint = 2e9 + 1000;

// #pragma GCC optimize("-ffast-math")

void solve();
void gen(vd& a);
void mergeSort(vd& input);
void startQuickSortLomuto(vd& input);
void multisetSort(vd& input);
void HeapSort(vd& input);
void RadixSort(vd& input);
bool CheckSorted(vd& input);
void priorityQueueSort(vd& input);

std::string typeOfSort(int& i) {
    switch(i) {
        case 0:
            return "mergeSort\n";
        case 1:
            return "quickSort\n";
        case 2:
            return "HeapSort\n";
        case 3:
            return "RadixSort\n";
        case 4:
            return "MultisetSort\n";
        case 5:
            return "PriorityQueueSort\n";
        case 6:
            return "StandartSort\n";
        default:
            return "Unknown\n";
    }
}

signed main() {
    //freopen("a.in", "r", stdin);
    //srand(time(0));

    std::vector<vd> ans(8, vd(7));
    int ind = 0;
    int start = 10;
    for (int n = start; n < 1e8; n *= 10) {
        int tests = 10;
        for (int w = 0; w < tests; ++w) {
            vd input(n);
            for (double& i : input) {
                i = rand() % 1000000 + static_cast<double>(rand()) / RAND_MAX;
            }
            vd inputCopy = input;
            auto start = std::chrono::high_resolution_clock::now();
            //mergeSort(inputCopy);
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = finish - start;
            ans[ind][0] += duration.count();

            input = inputCopy;
            start = std::chrono::high_resolution_clock::now();
            //startQuickSortLomuto(inputCopy);
            finish = std::chrono::high_resolution_clock::now();
            duration = finish - start;
            ans[ind][1] += duration.count();

            input = inputCopy;
            start = std::chrono::high_resolution_clock::now();
            //HeapSort(inputCopy);
            finish = std::chrono::high_resolution_clock::now();
            duration = finish - start;
            ans[ind][2] += duration.count();

            input = inputCopy;
            start = std::chrono::high_resolution_clock::now();
            //fastFunnySort(input);
            finish = std::chrono::high_resolution_clock::now();
            duration = finish - start;
            ans[ind][2] += duration.count();

            input = inputCopy;
            start = std::chrono::high_resolution_clock::now();
            //RadixSort(input);
            finish = std::chrono::high_resolution_clock::now();
            duration = finish - start;
            ans[ind][3] += duration.count();

            input = inputCopy;
            start = std::chrono::high_resolution_clock::now();
            //multisetSort(input);
            finish = std::chrono::high_resolution_clock::now();
            duration = finish - start;
            ans[ind][4] += duration.count();

            input = inputCopy;
            start = std::chrono::high_resolution_clock::now();
            //priorityQueueSort(input);
            finish = std::chrono::high_resolution_clock::now();
            duration = finish - start;
            ans[ind][5] += duration.count();

            input = inputCopy;
            start = std::chrono::high_resolution_clock::now();
            //std::sort(input.begin(), input.end());  // batska
            finish = std::chrono::high_resolution_clock::now();
            duration = finish - start;
            ans[ind][6] += duration.count();
        }
        for (double& i : ans[ind]) {
            i /= tests;
        }
        ++ind;
    }

    for (int t = 0; t < 7; ++t) {
        std::cout << typeOfSort(t);
        for (int i = 0; i < ind; ++i) {
            std::cout << std::fixed << std::setprecision(10) << ans[i][t] << '\n';
        }
        std::cout << "\n\n";
    }
}

/*bool Less(double& a, double& b) {
    return a < b + eps;
}
bool Equal(double& a, double& b) {
    return a - eps < b && a + eps > b;
}*/

void gen(vd& a) {
    for (double& i : a) {
        i = rand() % 20 + static_cast<double>(rand()) / RAND_MAX;
    }
}

std::vector<int> bubbleSort(std::vector<int>& input) {  // just why not)
    int len = input.size();
    std::vector<int> ret = input;
    bool notFinished = true;
    int start = 0;

    while (notFinished) {
        notFinished = false;
        for (int i = start; i < len - 1; ++i) {
            if (ret[i] > ret[i + 1]) {
                std::swap(ret[i], ret[i + 1]);
                notFinished = true;
            }
            if (! notFinished) start = i;
        }
    }

    return ret;
}

bool CheckSorted(vd& input) {  // it's funny, but I guessed to this after 3 sortings)
    for (int i = 1; i < input.size(); ++i) {
        if (input[i] < input[i - 1]) return false;
    }
    return true;
}

void mergeSort(vd& input) { // correct
    int k = 0;
    int n = input.size();
    vd in1 = input, in2 = input;
    bool inSec = true;
    while ((1 << k) < n) {
        int ind = 0;
        for (int start = 0; start < n; start += 1 << (k + 1)) {  // k + 1 because we see 2 groups by k
            int in1Index = start, in2Index = start + (1 << k);
            if (inSec) {  // in1 -> in2;
                for (int j = start; j < start + (1 << (k + 1)); ++j) {
                    int jnext = in2Index;
                    if (jnext >= start + (1 << (k + 1)) || jnext >= n) {
                        for (int i = in1Index; i < start + (1 << k); ++i) {  // there shouldn't be out of borders
                            if (i >= n) break;
                            in2[ind++] = in1[i];
                            ++in1Index;
                        }
                        break;
                    } else if (in1Index >= start + (1 << k) || in1Index >= n) {
                        if (in1Index >= n) break;
                        for (int i = in2Index; i < start + (1 << (k + 1)); ++i) {
                            if (i >= n) break;  // here can be out of borders
                            in2[ind++] = in1[i];
                            ++in2Index;
                        }
                    } else {
                        if (in1[in2Index] < in1[in1Index]) {
                            in2[ind] = in1[in2Index];
                            ++in2Index;
                        } else {
                            in2[ind] = in1[in1Index];
                            ++in1Index;
                        }
                        ++ind;
                    }
                }
            } else {  // in2 -> in1;
                for (int j = start; j < start + (1 << (k + 1)); ++j) {
                    int jnext = in2Index;
                    if (jnext >= start + (1 << (k + 1)) || jnext >= n) {
                        for (int i = in1Index; i < start + (1 << k); ++i) {  // there is no out of borders
                            if (i >= n) break;
                            in1[ind++] = in2[i];
                            ++in1Index;
                        }
                        break;
                    } else if (in1Index >= start + (1 << k)) {
                        for (int i = in2Index; i < start + (1 << (k + 1)); ++i) {
                            if (i >= n) break;  // here can be out of borders
                            in1[ind++] = in2[i];
                            ++in2Index;
                        }
                    } else {
                        if (in2[in2Index] < in2[in1Index]) {
                            in1[ind] = in2[in2Index];
                            ++in2Index;
                        } else {
                            in1[ind] = in2[in1Index];
                            ++in1Index;
                        }
                        ++ind;
                    }
                }
            }
        }
        ++k;
        inSec = !inSec;
    }
    input = (inSec ? in1 : in2);
}

void priorityQueueSort(vd& input) {
    std::priority_queue<double> yes;
    for (double& i : input) {
        yes.push(i);
    }
    for (double& i : input) {
        i = yes.top();
        yes.pop();
    }
}

void multisetSort(vd& input) {
    std::multiset<double> yes;
    for (double& i : input) {
        yes.insert(i);
    }
    int ind = 0;
    for (double i : yes) {
        input[ind++] = i;
    }
}

int iterOfQuickSortLomuto(vd& input, int start, int finish) {  // [start, finish)
    double p = input[start + (finish - start) / 2];
    int i = start;
    int left = 0, right = 0;
    for (int j = start; j < finish; ++j) {
        if (input[j] < p || input[j] == p && left < right) {
            std::swap(input[i], input[j]);
            ++left;
            ++i;
        } else {
            ++right;
        }
    }
    return i;
}

void quickSortLomuto(vd& input, int left, int right) {  // [left, right)
    if (left == right - 1) {
        return;
    } else {
        int bet = iterOfQuickSortLomuto(input, left, right);

        quickSortLomuto(input, left, bet);
        quickSortLomuto(input, bet, right);
    }
}

void startQuickSortLomuto(vd& input) {
    vd vcopy = input;
    vd inQSort = input;
    int n = input.size();

    quickSortLomuto(input, 0, n);

}


vd heap;

int findMinIndex() {  // first el in heap - min
    return 0;
}

void SiftDown(int ind) {  // index starts from 0
    int x1 = ind * 2 + 1, x2 = ind * 2 + 2;
    if (heap.size() <= x1) return;
    double f = heap[x1], s = f;
    if (heap.size() > x2) s = heap[x2];
    if (f <= s) {
        if (heap[ind] > f) {
            std::swap(heap[x1], heap[ind]);
            SiftDown(x1);
        }
    } else if (heap[ind] > s){
        std::swap(heap[x2], heap[ind]);
        SiftDown(x2);
    }
}

void SiftUp(int ind) {
    int pred = (ind - 1) / 2;
    if (heap[pred] > heap[ind]) {
        std::swap(heap[pred], heap[ind]);
        SiftUp(ind);
    }
}

void Insert(double val) {
    heap.push_back(val);
    SiftUp(static_cast<int>(heap.size()) - 1);
}

double ExtractMin() {  // return min value and remove it from heap
    double ret = heap[0];
    std::swap(heap[0], heap.back());
    heap.pop_back();
    SiftDown(0);
    return ret;
}

void MakeHeap() {
    int start = heap.size();
    for (int i = start - 1; i >= 0; --i) {
        SiftDown(i);
    }
}

void HeapSort(vd& input) {
    heap = input;
    MakeHeap();
    int n = input.size();
    for (double& i : input) {
        i = ExtractMin();
    }
}

void RadixSort(vd& input) {  // only for positive
    int n = input.size();
    vd inputCopy = input;
    std::vector<std::string> byteIdea(n);
    vi indexSorted(n);
    for (int i = 0; i < n; ++i) {
        double now = input[i];
        indexSorted[i] = i;
        // https://stackoverflow.com/questions/32259970/get-exact-bit-representation-of-a-double-in-c
        uint8_t *bytePointer = (uint8_t *)& now;
        for(size_t index = 0; index < sizeof(double); index++)
        {
            uint8_t byte = bytePointer[index];

            for(int bit = 0; bit < 8; bit++)
            {
                byteIdea[i] += '0' + static_cast<int>(byte&1);
                byte >>= 1;
            }
        }
    }

    std::vector<vi> buffer(2, vi(n));
    for (int w = 0; w < 64; ++w) {
        int zero = 0, one = 0;
        for (int i = 0; i < n; ++i) {
            if (byteIdea[indexSorted[i]][w] == '0') {
                buffer[0][zero++] = indexSorted[i];
            } else {
                buffer[1][one++] = indexSorted[i];
            }
        }
        for (int i = 0; i < zero; ++i) {
            indexSorted[i] = buffer[0][i];
        }
        for (int i = zero; i < n; ++i) {
            indexSorted[i] = buffer[1][i - zero];
        }
    }
    for (int i = 0; i < n; ++i) {
        input[i] = inputCopy[indexSorted[i]];
    }
}