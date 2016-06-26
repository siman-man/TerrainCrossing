#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>
#include <limits>
#include <map>
#include <math.h>
#include <sstream>
#include <string.h>
#include <vector>
#include <set>
#include <queue>

using namespace std;

typedef long long ll;

const int MAX_S = 50;
const int MAX_N = 250;
const int UNKNOWN = -1;
const int ITEM = 0;
const int TARGET = 1;

const ll CYCLE_PER_SEC = 2400000000;
double MAX_TIME = 10.0;

const int INF = 99999;

const int DY[4] = {-1, 0, 0, 1};
const int DX[4] = {0, -1, 1, 0};

const double EPS = 0.001;

unsigned long long xor128(){
  static unsigned long long rx=123456789, ry=362436069, rz=521288629, rw=88675123;
  unsigned long long rt = (rx ^ (rx<<11));
  rx=ry; ry=rz; rz=rw;
  return (rw=(rw^(rw>>19))^(rt^(rt>>8)));
}

unsigned long long int getCycle() {
  unsigned int low, high;
  __asm__ volatile ("rdtsc" : "=a" (low), "=d" (high));
  return ((unsigned long long int)low) | ((unsigned long long int)high << 32);
}

double getTime(unsigned long long int begin_cycle) {
  return (double)(getCycle() - begin_cycle) / CYCLE_PER_SEC;
}

int S;
int V;
int N;
int g_capacity;

struct Location {
  double y;
  double x;
  int yi;
  int xi;
  int nid;
  bool locked;

  Location(double y = -1.0, double x = -1.0, bool locked = false) {
    this->y = y;
    this->x = x;
    this->yi = floor(y);
    this->xi = floor(x);
    this->locked = locked;
    this->nid = yi * S + xi;
  }

  double dist(Location other) {
    return sqrt((other.y-y)*(other.y-y)+(other.x-x)*(other.x-x));
  }

  int manhattan(Location other) {
    return abs(other.yi - yi) + abs(other.xi - xi);
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

  bool isItem() {
    return type == ITEM;
  }

  bool isTarget() {
    return type == TARGET;
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
ll g_startCycle;

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

      g_startCycle = getCycle();

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

      vector<int> path1 = createFirstPathNN();
      assert(isValidPath(path1));
      vector<int> path2 = createFirstPathFI();
      assert(isValidPath(path2));

      double score1 = calcCost(path1);
      double score2 = calcCost(path2);

      vector<int> path = (score1 < score2)? path1 : path2;

      path = cleanPath(path);
      vector<Location> result;
      Object *obj = getObject(path[0]);
      vector<Location> startPath = leaveMap(path[0]);
      reverse(startPath.begin(), startPath.end());

      int ssize = startPath.size();
      fprintf(stderr,"start path size = %d\n", ssize);
      for (int i = 0; i < ssize; i++) {
        Location coord = startPath[i];
        result.push_back(coord);
      }
      result[0].locked = true;

      result.push_back(Location(obj->y, obj->x, true));
      int rsize = result.size();

      for (int i = 0; i < path.size()-1; i++) {
        int from = path[i];
        int to = path[i+1];

        vector<Location> pa = createPath(from, to);

        for (int j = 0; j < pa.size(); j++) {
          Location coord1 = result[rsize-1];
          Location coord2 = pa[j];

          if (!coord1.locked) {
            int nid1 = floor(coord1.y)*S + floor(coord1.x);
            int nid2 = floor(coord2.y)*S + floor(coord2.x);

            if (nid1 == nid2) {
              result[rsize-1] = coord2;
              continue;
            }
          }

          result.push_back(pa[j]);
          rsize++;
        }
      }

      vector<Location> endPath = leaveMap(path[path.size()-1]);
      int esize = endPath.size();
      fprintf(stderr,"end path size = %d\n", esize);
      for (int i = 0; i < esize; i++) {
        Location coord = endPath[i];
        result.push_back(coord);
      }

      result[result.size()-1].locked = true;

      assert(result[0].locked);
      fixPath(result);
      assert(isValidAnswer(result));
      ret = path2answer(result);

      return ret;
    }

    void fixPath(vector<Location> &path) {
      int psize = path.size();
      double minCost = calcCostDetail(path);

      double timeLimit = 9.5;
      double currentTime = getTime(g_startCycle);
      ll tryCount = 0;

      while (currentTime < timeLimit) {
        int index = xor128() % psize;
        Location coord = path[index];
        if (coord.locked) continue;
        Location temp = coord;

        int d1 = xor128()%10;
        int d2 = xor128()%10;
        coord.y += 0.01 * d1 - 0.05;
        coord.x += 0.01 * d2 - 0.05;
        int nid = floor(coord.y)*S + floor(coord.x);

        if (coord.y - coord.yi < EPS) continue;
        if (coord.x - coord.xi < EPS) continue;
        if ((coord.yi+1) - coord.y < EPS) continue;
        if ((coord.xi+1) - coord.x < EPS) continue;

        if (nid != coord.nid) continue;
        path[index] = coord;
        double cost = calcCostDetail(path);

        if (minCost > cost) {
          minCost = cost;
        } else {
          path[index] = temp;
        }

        currentTime = getTime(g_startCycle);
        tryCount++;
      }

      fprintf(stderr,"tryCount fix = %lld\n", tryCount);
    }

    vector<int> cleanPath(vector<int> path) {
      vector<int> bestPath = path;
      vector<int> goodPath = path;
      double timeLimit = 5.0;
      double currentTime = getTime(g_startCycle);
      double minCost = calcCost(bestPath);
      double goodCost = minCost;
      fprintf(stderr,"minCost = %f\n", minCost);
      ll tryCount = 0;
      ll invalidCount = 0;
      int c1, c2;
      int psize = path.size();

      double T = 10000.0;
      double k = 10.0;
      double alpha = 0.999;

      while (currentTime < timeLimit) {
        do {
          c1 = xor128() % psize;
          c2 = xor128() % psize;
        } while (c1 == c2);

        swapObject(path, c1, c2);
        double cost = calcCost(path);

        if (cost == DBL_MAX) {
          invalidCount++;
        }

        if (minCost > cost) {
          minCost = cost;
          bestPath = path;
        } else {
          swapObject(path, c1, c2);
        }

        /*
        double diffCost = goodCost - minCost;

        if (goodCost > cost) {
          goodCost = cost;
          goodPath = path;
        } else if (xor128()%100 < 100*exp(diffCost/(k*T))) {
          goodCost = cost;
          goodPath = path;
        } else {
          swapObject(path, c1, c2);
        }
        */

        T *= alpha;
        currentTime = getTime(g_startCycle);
        tryCount++;
      }

      fprintf(stderr,"(%lld/%lld), tryCount = %lld\n", invalidCount, tryCount, tryCount);

      return bestPath;
    }

    void swapObject(vector<int> &path, int c1, int c2) {
      int temp = path[c1];
      path[c1] = path[c2];
      path[c2] = temp;
    }

    vector<double> path2answer(vector<Location> path) {
      vector<double> answer;

      int psize = path.size();
      for (int i = 0; i < psize; i++) {
        Location coord = path[i];
        answer.push_back(coord.x);
        answer.push_back(coord.y);
        //fprintf(stderr,"%d: y = %f, x = %f\n", i, coord.y, coord.x);
      }

      return answer;
    }

    vector<Location> createPath(int from, int to) {
      vector<Location> path;
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
            Location coord = nid2coord(id);
            if (id == 0) coord.locked = true;
            path.push_back(coord);
          }

          path.push_back(Location(toObj->y, toObj->x, true));
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

    vector<int> createFirstPathFI() {
      vector<int> ret;
      int itemCount[2*N];
      memset(itemCount, 0, sizeof(itemCount));
      map<int, bool> checkList;

      double minDist = DBL_MAX;
      int minId = -1;
      for (int i = 0; i < N; i++) {
        Item *item = getItem(i);
        double dist = calcDist(item->y, item->x, S/2.0, S/2.0);

        if (minDist > dist) {
          minDist = dist;
          minId = i;
        }
      }
      assert(minId >= 0);
      ret.push_back(minId);
      checkList[minId] = true;

      for (int i = 0; i < 2*N-1; i++) {
        int cnt = 0;

        for (int j = 0; j < i+1; j++) {
          int oid = ret[j];
          Object *obj = getObject(oid);

          if (obj->isItem()) {
            cnt = min(g_capacity, cnt+1);
          } else {
            cnt--;
          }
          assert(cnt >= 0);
          itemCount[j] = cnt;
        }

        double maxDist = -1;
        int maxId = -1;

        if (i % 2 == 0) {
          for (int j = 0; j < N; j++) {
            Target *target = getTarget(j);
            if (checkList[target->id]) continue;

            double dist = calcDist(target->y, target->x, S/2.0, S/2.0);

            if (maxDist < dist) {
              maxDist = dist;
              maxId = target->id;
            }
          }

          double minCost = DBL_MAX;
          int index = -1;

          int rsize = ret.size();
          for (int k = 0; k < rsize; k++) {
            if (itemCount[k] == 0) continue;
            vector<int> temp = ret;
            temp.insert(temp.begin()+(k+1), maxId);
            if (!isValidPath(temp)) continue;
            int o1 = ret[k];
            int o2 = ret[(k+1)%rsize];
            double c1 = g_pathCost[o1][maxId];
            double c2 = g_pathCost[maxId][o2];
            double c3 = g_pathCost[o1][o2];
            double c = c1 + c2 - c3;

            if (minCost > c) {
              minCost = c;
              index = k+1;
            }
          }

          assert(index >= 0);
          checkList[maxId] = true;
          ret.insert(ret.begin()+index, maxId);
        } else {
          for (int j = 0; j < N; j++) {
            Item *item = getItem(j);
            if (checkList[item->id]) continue;

            double dist = calcDist(item->y, item->x, S/2.0, S/2.0);

            if (maxDist < dist) {
              maxDist = dist;
              maxId = item->id;
            }
          }

          double minCost = DBL_MAX;
          int index = -1;

          int rsize = ret.size();
          for (int k = 0; k < rsize; k++) {
            int o1 = ret[k];
            int o2 = ret[(k+1)%rsize];
            double c1 = g_pathCost[o1][maxId];
            double c2 = g_pathCost[maxId][o2];
            double c3 = g_pathCost[o1][o2];
            double c = c1 + c2 - c3;

            if (minCost > c) {
              minCost = c;

              if (k == rsize-1) {
                index = 0;
              } else {
                index = k+1;
              }
            }
          }

          assert(index >= 0);
          checkList[maxId] = true;
          ret.insert(ret.begin()+index, maxId);
        }
      }

      return ret;
    }

    vector<Location> leaveMap(int oid) {
      vector<Location> path;
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
            path.push_back(Location(cy, cx));
          }

          if (y == 0) { 
            path.push_back(Location(0.0001, x+0.5));
          } else if (x == 0) {
            path.push_back(Location(y+0.5, 0.0001));
          } else if (y == S-1) {
            path.push_back(Location(S-EPS/2, x+0.5));
          } else {
            path.push_back(Location(y+0.5, S-EPS/2));
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

    double costSeg(Location l1, Location l2) {
      if (l1.manhattan(l2) == 0) {
        return l1.dist(l2) * g_fieldCost[l1.yi][l1.xi];
      }

      int t1 = g_fieldCost[l1.yi][l1.xi];
      int t2 = g_fieldCost[l2.yi][l2.xi];
      double score = pow(t1 - t2, 2);

      double x0, y0;
      if (l1.yi == l2.yi) {
        x0 = max(l1.xi, l2.xi);
        y0 = l1.y + (l2.y - l1.y) * (x0 - l1.x) / (l2.x - l1.x);
      } else {
        y0 = max(l1.yi, l2.yi);
        x0 = l1.x + (l2.x - l1.x) * (y0 - l1.y) / (l2.y - l1.y);
      }

      Location intersection(y0, x0);

      return score + l1.dist(intersection) * t1 + l2.dist(intersection) * t2;
    }

    double calcCostDetail(vector<Location> &path) {
      double score = 0.0;
      int psize = path.size();

      for (int i = 1; i < psize; i++) {
        score += costSeg(path[i-1], path[i]);
      }

      return score;
    }

    double calcCost(vector<int> &path) {
      double score = 0.0;
      int psize = path.size();
      int itemCount = 0;

      for (int i = 0; i < psize; i++) {
        int oid = path[i];
        Object *obj = getObject(oid);

        itemCount = (obj->isItem())? min(g_capacity, itemCount+1) : itemCount-1;
        if (itemCount < 0) return DBL_MAX;

        score += g_pathCost[path[i]][path[(i+1)%psize]];
      }

      return score;
    }

    bool isValidAnswer(vector<Location> &path) {
      int psize = path.size();

      for (int i = 0; i < psize; i++) {
        Location l = path[i];

        if (i == 364) {
          fprintf(stderr,"y = %f, x = %f\n", l.y, l.x);
        }

        if (isNearInnerBorder(l)) return false;
      }

      return true;
    }

    bool isValidPath(vector<int> &path) {
      int psize = path.size();
      int cnt = 0;
      ///fprintf(stderr,"psize = %d\n", psize);

      for (int i = 0; i < psize; i++) {
        int oid = path[i];
        Object *obj = getObject(oid);

        if (obj->isItem()) {
          cnt = min(g_capacity, cnt+1);
        } else {
          cnt--;
        }

        if (cnt < 0) return false;
      }

      return true;
    }

    Location nid2coord(int nid) {
      return Location(nid / S + 0.5, nid % S + 0.5);
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

    inline bool isNearPoint(Location c1, Location c2) {
      return calcDist(c1.y, c1.x, c2.y, c2.x) < EPS;
    }

    bool isNearInnerBorder(Location l) {
      double inCellY = l.y - l.yi;
      double inCellX = l.x - l.xi;
      return ((l.xi > 0 && inCellX < EPS) || 
              (l.yi > 0 && inCellY < EPS) ||
              (l.xi < S-1 && 1 - inCellX < EPS) ||
              (l.yi < S-1 && 1 - inCellY < EPS));
    }

    inline double calcDist(double y1, double x1, double y2, double x2) {
      return sqrt((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2));
    }

    inline bool isBorder(int y, int x) {
      return (y == 0 || y == S-1 || x == 0 || x == S-1);
    }

    inline bool isInside(int y, int x) {
      return (0 <= y && y < S && 0 <= x && x < S);
    }

    inline bool isOutside(int y, int x) {
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
