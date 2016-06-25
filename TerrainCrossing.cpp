#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <math.h>
#include <sstream>
#include <string.h>
#include <vector>
#include <set>
#include <queue>

using namespace std;

const int MAX_S = 50;
const int MAX_N = 250;
const int UNKNOWN = -1;
const int ITEM = 0;
const int TARGET = 1;

const int INF = 99999;

const int DY[4] = {-1, 0, 0, 1};
const int DX[4] = {0, -1, 1, 0};

const double EPS = 0.001;

struct Coord {
  double y;
  double x;

  Coord(double y = 0.0, double x = 0.0) {
    this->y = y;
    this->x = x;
  }
};

struct Node {
  int cid;
  double y;
  double x;
  double cost;
  vector<int> ids;

  Node(int cid, double cost) {
    this->cid = cid;
    this->cost = cost;
  }

  Node dup() {
    Node node(cid, cost);
    node.y = y;
    node.x = x;
    node.ids = ids;

    return node;
  }

  bool operator >(const Node &n) const{
    return cost > n.cost;
  }    
};

struct Object {
  int id;
  int nid;
  int type;
  double y;
  double x;

  Object(int id, int type, double y, double x) {
    this->id = id;
    this->type = type;
    this->y = y;
    this->x = x;
  }
};

struct Cell {
  vector<int> objIdList;
};

typedef Object Item;
typedef Object Target;

int g_fieldCost[MAX_S][MAX_S];
double g_pathCost[2*MAX_N][2*MAX_N];
bool g_existObject[MAX_S][MAX_S];
int S;
int V;
int N;
int g_capacity;

Cell g_field[MAX_S][MAX_S];
vector<Object> g_objectList;
vector<Object> g_itemList;
vector<Object> g_targetList;

class TerrainCrossing {
  public:
    void init(vector<string> map, vector<double> locations, int capacity) {
      S = map.size();
      N = locations.size() / 4;
      g_capacity = capacity;

      memset(g_existObject, false, sizeof(g_existObject));

      fprintf(stderr,"S = %d, N = %d, capacity = %d\n", S, N, capacity);

      V = S*S;

      for (int y = 0; y < S; y++) {
        for (int x = 0; x < S; x++) {
          g_fieldCost[y][x] = map[y][x] - '0';
        }
      }

      for (int i = 0; i < 2*N; i++) {
        for (int j = 0; j < 2*N; j++) {
          g_pathCost[i][j] = INF;
        }
      }

      for (int i = 0; i < 2*N; i++) {
        double y = locations[2*i+1];
        double x = locations[2*i];
        int iy = floor(y);
        int ix = floor(x);

        Object obj(i, UNKNOWN, y, x);
        obj.nid = iy * S + ix;
        g_existObject[iy][ix] = true;
        g_field[iy][ix].objIdList.push_back(i);

        if (i < N) {
          obj.type = ITEM;
          g_itemList.push_back(obj);
        } else {
          obj.type = TARGET;
          g_targetList.push_back(obj);
        }

        g_objectList.push_back(obj);
      }

      calcPathCost();
    }

    void calcPathCost() {
      fprintf(stderr,"calcPathCost =>\n");

      for (int i = 0; i < 2*N; i++) {
        calcMinPathCost(i);
      }
    }

    void calcMinPathCost(int oid) {
      Object *obj = getObject(oid);
      Node root(obj->nid, 0.0);

      priority_queue<Node, vector<Node>, greater<Node> > pque;
      pque.push(root);

      map<int, bool> checkList;

      while (!pque.empty()) {
        Node node = pque.top(); pque.pop();

        int y = node.cid / S;
        int x = node.cid % S;

        if (checkList[node.cid]) continue;
        checkList[node.cid] = true;

        int size = g_field[y][x].objIdList.size();
        for (int i = 0; i < size; i++) {
          g_pathCost[oid][g_field[y][x].objIdList[i]] = node.cost;
        }

        for (int i = 0; i < 4; i++) {
          int ny = y + DY[i];
          int nx = x + DX[i];
          if (isOutside(ny, nx)) continue;

          int nid = ny*S + nx;
          Node next = node.dup();
          next.cid = nid;
          next.cost += g_fieldCost[ny][nx] + pow(g_fieldCost[y][x] - g_fieldCost[ny][nx], 2);
          //next.ids.push_back(nid);
          pque.push(next);
        }
      }
    }

    vector<double> getPath(vector<string> map, vector<double> locations, int capacity) {
      vector<double> ret;

      init(map, locations, capacity);

      fprintf(stderr,"start =>\n");

      vector<int> path = createFirstPathNN();
      vector<Coord> result;
      Object *obj = getObject(path[0]);
      vector<Coord> startPath = leaveMap(path[0]);
      reverse(startPath.begin(), startPath.end());

      int ssize = startPath.size();
      fprintf(stderr,"start path size = %d\n", ssize);
      for (int i = 0; i < ssize; i++) {
        Coord coord = startPath[i];
        result.push_back(coord);
      }

      result.push_back(Coord(obj->y, obj->x));

      for (int i = 0; i < path.size()-1; i++) {
        int from = path[i];
        int to = path[i+1];

        vector<Coord> pa = createPath(from, to);

        for (int j = 0; j < pa.size(); j++) {
          result.push_back(pa[j]);
        }
      }

      vector<Coord> endPath = leaveMap(path[path.size()-1]);
      int esize = endPath.size();
      fprintf(stderr,"end path size = %d\n", esize);
      for (int i = 0; i < esize; i++) {
        Coord coord = endPath[i];
        result.push_back(coord);
      }

      ret = path2answer(result);

      return ret;
    }

    vector<double> path2answer(vector<Coord> path) {
      vector<double> answer;

      int psize = path.size();
      for (int i = 0; i < psize; i++) {
        Coord coord = path[i];
        answer.push_back(coord.x);
        answer.push_back(coord.y);
        //fprintf(stderr,"%d: y = %f, x = %f\n", i, coord.y, coord.x);
      }

      return answer;
    }

    vector<Coord> createPath(int from, int to) {
      vector<Coord> path;
      Object *fromObj = getObject(from);
      Object *toObj = getObject(to);

      Node root(fromObj->nid, 0.0);

      priority_queue<Node, vector<Node>, greater<Node> > pque;
      pque.push(root);
      map<int, bool> checkList;

      while (!pque.empty()) {
        Node node = pque.top(); pque.pop();
        int y = node.cid / S;
        int x = node.cid % S;

        if (checkList[node.cid]) continue;
        checkList[node.cid] = true;

        if (node.cid == toObj->nid) {
          int isize = node.ids.size();

          for (int i = 0; i < isize; i++) {
            int id = node.ids[i];
            Coord coord = nid2coord(id);
            path.push_back(coord);
          }

          path.push_back(Coord(toObj->y, toObj->x));
          break;
        }

        for (int i = 0; i < 4; i++) {
          int ny = y + DY[i];
          int nx = x + DX[i];
          if (isOutside(ny, nx)) continue;

          int nid = ny*S + nx;
          Node next = node.dup();
          next.cid = nid;
          next.cost += g_fieldCost[ny][nx] + pow(g_fieldCost[y][x] - g_fieldCost[ny][nx], 2);
          next.ids.push_back(nid);
          pque.push(next);
        }
      }

      return path;
    }

    vector<int> createFirstPathNN() {
      vector<int> ret;
      map<int, bool> checkList;

      ret.push_back(g_itemList[0].id);
      checkList[g_itemList[0].id] = true;

      for (int i = 0; i < 2*N-1; i++) {
        int from = ret[i];
        double minCost = DBL_MAX;
        int minId = -1;

        for (int j = 0; j < N; j++) {
          Object *obj = (i % 2 == 0)? getTarget(j) : getItem(j);

          if (checkList[obj->id]) continue;
          double cost = g_pathCost[from][obj->id];

          if (minCost > cost) {
            minCost = cost;
            minId = obj->id;
          }
        }

        assert(minId >= 0);
        ret.push_back(minId);
        checkList[minId] = true;
      }

      return ret;
    }

    vector<Coord> leaveMap(int oid) {
      vector<Coord> path;
      Object *obj = getObject(oid);

      fprintf(stderr,"leave from %d, y = %f, x = %f\n", oid, obj->y, obj->x);

      priority_queue<Node, vector<Node>, greater<Node> > pque;
      Node root(obj->nid, 0.0);
      pque.push(root);
      map<int, bool> checkList;

      while (!pque.empty()) {
        Node node = pque.top(); pque.pop();
        int y = node.cid / S;
        int x = node.cid % S;

        if (checkList[node.cid]) continue;
        checkList[node.cid] = true;

        if (isBorder(y, x)) {
          int isize = node.ids.size();
          for (int i = 0; i < isize; i++) {
            int id = node.ids[i];
            double cy = id / S + 0.5;
            double cx = id % S + 0.5;
            path.push_back(Coord(cy, cx));
          }

          if (y == 0) { 
            path.push_back(Coord(0.0001, x+0.5));
          } else if (x == 0) {
            path.push_back(Coord(y+0.5, 0.0001));
          } else if (y == S-1) {
            path.push_back(Coord(S-EPS/2, x+0.5));
          } else {
            path.push_back(Coord(y+0.5, S-EPS/2));
          }

          break;
        }

        for (int i = 0; i < 4; i++) {
          int ny = y + DY[i];
          int nx = x + DX[i];
          if (isOutside(ny, nx)) continue;

          int nid = ny*S + nx;
          Node next = node.dup();
          next.cid = nid;
          next.cost += g_fieldCost[ny][nx] + pow(g_fieldCost[y][x] - g_fieldCost[ny][nx], 2);
          next.ids.push_back(nid);
          pque.push(next);
        }
      }

      return path;
    }

    Coord nid2coord(int nid) {
      return Coord(nid / S + 0.5, nid % S + 0.5);
    }

    Object *getObject(int oid) {
      return &g_objectList[oid];
    }

    Item *getItem(int id) {
      return &g_itemList[id];
    }

    Target *getTarget(int id) {
      return &g_targetList[id];
    }

    bool isNearPoint(Coord c1, Coord c2) {
      return calcDist(c1.y, c1.x, c2.y, c2.x) < EPS;
    }

    double calcDist(int y1, int x1, int y2, int x2) {
      return sqrt((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2));
    }

    bool isBorder(int y, int x) {
      return (y == 0 || y == S-1 || x == 0 || x == S-1);
    }

    bool isInside(int y, int x) {
      return (0 <= y && y < S && 0 <= x && x < S);
    }

    bool isOutside(int y, int x) {
      return (y < 0 || S <= y || x < 0 || S <= x);
    }
};


// -------8<------- end of solution submitted to the website -------8<-------
template<class T> void getVector(vector<T>& v) { for (int i = 0; i < v.size(); ++i) cin >> v[i];}
int main() {
  TerrainCrossing tc; int L, M, capacity;
  cin >> M;
  vector<string> map(M); getVector(map);
  cin >> L;
  vector<double> locations(L); getVector(locations);
  cin >> capacity;
  vector<double> ret = tc.getPath(map, locations, capacity);
  cout << ret.size() << endl;
  for (int i = 0; i < ret.size(); ++i) cout << ret[i] << endl;
  cout.flush();
}
