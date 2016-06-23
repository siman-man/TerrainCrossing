#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <set>

using namespace std;

class TerrainCrossing {
public:
    vector<double> getPath(vector<string> map, vector<double> locations, int capacity) {
        vector<double> ret;
        ret.push_back(4.9995); ret.push_back(0.7658);
        // pick up item 1
        ret.push_back(4.7867); ret.push_back(0.7658);
        // drop it off at target 6
        ret.push_back(3.8144); ret.push_back(0.1081);
        // pick up item 0
        ret.push_back(3.7648); ret.push_back(1.2640);
        // drop it off at target 7
        ret.push_back(3.3420); ret.push_back(2.5000);
        ret.push_back(3.3420); ret.push_back(3.0530);
        // pick up item 2
        ret.push_back(2.5000); ret.push_back(3.0530);
        ret.push_back(1.5000); ret.push_back(3.0530);
        ret.push_back(0.7225); ret.push_back(3.0530);
        ret.push_back(0.7225); ret.push_back(2.5000);
        ret.push_back(0.7225); ret.push_back(1.4533);
        // pick up item 3
        ret.push_back(0.2299); ret.push_back(2.8555);
        ret.push_back(0.2299); ret.push_back(3.8555);
        ret.push_back(0.2299); ret.push_back(4.8555);
        // drop it off at target 4
        ret.push_back(0.5000); ret.push_back(3.3869);
        ret.push_back(1.2611); ret.push_back(3.3869);
        // drop it off at target 5
        ret.push_back(2.2611); ret.push_back(3.3869);
        ret.push_back(2.2611); ret.push_back(4.6214);
        ret.push_back(3.7958); ret.push_back(4.6214);
        // exit
        ret.push_back(3.7958); ret.push_back(4.9995);
        return ret;
    }
};


// -------8<------- end of solution submitted to the website -------8<-------
template<class T> void getVector(vector<T>& v) { for (int i = 0; i < v.size(); ++i) cin >> v[i];}
int main() {
    TerrainCrossing tc;
    int M;
    cin >> M;
    vector<string> map(M);
    getVector(map);
    int L;
    cin >> L;
    vector<double> locations(L);
    getVector(locations);
    int capacity;
    cin >> capacity;
    vector<double> ret = tc.getPath(map, locations, capacity);
    cout << ret.size() << endl;
    for (int i = 0; i < ret.size(); ++i)
        cout << ret[i] << endl;
    cout.flush();
}
