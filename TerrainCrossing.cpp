#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>
#include <bitset>
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
const int ITEM = 1;
const int TARGET = -1;

const ll CYCLE_PER_SEC = 2400000000;
double MAX_TIME = 10.2;
bool g_debug = false;

const short DY[4] = {-1, 0, 0, 1};
const short DX[4] = {0, -1, 1, 0};

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

  void update() {
    this->yi = floor(y);
    this->xi = floor(x);
  }

  int manhattan(Location other) {
    return abs(other.yi - yi) + abs(other.xi - xi);
  }
};

struct Node {
  short cid;
  short yi;
  short xi;
  short beforeDirect;
  short length;
  double y;
  double x;
  double cost;
  double beforeY;
  double beforeX;
  short step[3];
  short fc[3];
  vector<short> ids;
  bool update;

  Node(short cid, double cost, double y, double x) {
    this->cid = cid;
    this->cost = cost;
    this->y = y;
    this->x = x;
    this->yi = floor(y);
    this->xi = floor(x);
    this->beforeDirect = -1;
    this->beforeY = y;
    this->beforeX = x;
    this->length = 0;
  }

  Node dup() {
    Node node(cid, cost, y, x);
    node.ids = ids;
    node.beforeDirect = beforeDirect;

    node.step[2] = step[1];
    node.step[1] = step[0];

    node.fc[2] = fc[1];
    node.fc[1] = fc[0];

    node.update = false;

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
ll g_startCycle;

Cell g_field[MAX_S][MAX_S];
vector<Object> g_objectList;
int g_objectType[2*MAX_N];
double g_leaveCost[2*MAX_N];

double g_costHistory[2*MAX_N];
double g_itemHistory[2*MAX_N];

class TerrainCrossing {
  public:
    void init(vector<string> map, vector<double> locations, int capacity) {
      S = map.size();
      N = locations.size() / 4;
      g_capacity = capacity;

      g_startCycle = getCycle();

      fprintf(stderr,"S = %d, N = %d, capacity = %d\n", S, N, capacity);

      V = S*S;

      for (int y = 0; y < S; y++) {
        for (int x = 0; x < S; x++) {
          g_fieldCost[y][x] = map[y][x] - '0';
        }
      }

      for (int i = 0; i < 2*N; i++) {
        g_leaveCost[i] = DBL_MAX;

        for (int j = 0; j < 2*N; j++) {
          g_pathCost[i][j] = DBL_MAX;
        }
      }

      for (int i = 0; i < 2*N; i++) {
        double y = locations[2*i+1];
        double x = locations[2*i];
        int iy = floor(y);
        int ix = floor(x);

        Object obj(i, UNKNOWN, y, x);
        obj.nid = iy * S + ix;
        obj.type = (i < N)? ITEM : TARGET;
        g_field[iy][ix].objIdList.push_back(i);
        g_objectType[i] = obj.type;

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
      Node root(obj->nid, 0.0, obj->y, obj->x);
      g_pathCost[oid][oid] = 0.0;

      priority_queue<Node, vector<Node>, greater<Node> > pque;
      pque.push(root);

      vector<bool> checkList(V);

      while (!pque.empty()) {
        Node node = pque.top(); pque.pop();

        short y = node.yi;
        short x = node.xi;

        if (checkList[node.cid]) continue;
        checkList[node.cid] = true;

        if (y == 0 || y == S-1 || x == 0 || x == S-1) {
          g_leaveCost[oid] = min(g_leaveCost[oid], node.cost + 0.5*g_fieldCost[y][x]);
        }

        Cell cell = g_field[y][x];
        int size = cell.objIdList.size();
        for (int i = 0; i < size; i++) {
          int eid = cell.objIdList[i];
          if (oid == cell.objIdList[i]) continue;
          Object *other = getObject(eid);

          if (node.length == 0) {
            double ncost = costSeg(Location(obj->y, obj->x), Location(other->y, other->x));
            g_pathCost[oid][eid] = ncost;
          } else {
            double bcost = costSeg(Location(node.beforeY, node.beforeX), Location(node.y, node.x));
            double ncost = costSeg(Location(node.beforeY, node.beforeX), Location(other->y, other->x));
            g_pathCost[oid][eid] = node.cost + (ncost - bcost);
          }
        }

        for (int i = 0; i < 4; i++) {
          short ny = y + DY[i];
          short nx = x + DX[i];
          if (isOutside(ny, nx)) continue;

          short nid = ny*S + nx;
          Node next = node.dup();
          next.cid = nid;
          next.beforeDirect = i;
          next.length = node.length + 1;
          next.yi = ny;
          next.xi = nx;
          next.y = ny + 0.5;
          next.x = nx + 0.5;

          double costB = costSeg(Location(node.y, node.x), Location(next.y, next.x));
          double costC = pow(g_fieldCost[y][x]-g_fieldCost[ny][nx], 2);
          next.cost += costB + costC;

          next.step[0] = i;
          next.fc[0] = g_fieldCost[ny][nx];

          updateNodeCost(next);
          pque.push(next);
        }
      }
    }

    vector<double> getPath(vector<string> map, vector<double> locations, int capacity) {
      vector<double> ret;

      init(map, locations, capacity);

      fprintf(stderr,"start =>\n");

      vector<short> path1 = createFirstPathNN();

      cerr.flush();
      assert(isValidPath(path1));
      vector<short> path2 = createFirstPathFI();
      assert(isValidPath(path2));

      double score1 = calcCost(path1);
      double score2 = calcCost(path2);

      vector<short> path = (score1 < score2)? path1 : path2;

      double currentTime = getTime(g_startCycle);
      fprintf(stderr,"current time = %f\n", currentTime);

      path = cleanPath(path);
      vector<Location> result;
      Object *obj = getObject(path[0]);
      vector<Location> startPath = leaveMap(path[0]);
      reverse(startPath.begin(), startPath.end());

      int ssize = startPath.size();
      fprintf(stderr,"start path size = %d\n", ssize);
      for (int i = 0; i < ssize; i++) {
        Location l = startPath[i];
        result.push_back(l);
      }
      result[0].locked = true;

      result.push_back(Location(obj->y, obj->x, true));
      int rsize = result.size();

      for (int i = 0; i < path.size()-1; i++) {
        short from = path[i];
        short to = path[i+1];

        vector<Location> pa = createPath(from, to);

        for (int j = 1; j < pa.size(); j++) {
          Location l1 = result[rsize-1];
          Location l2 = pa[j];

          if (!l1.locked) {
            if (l1.manhattan(l2) == 0) {
              result[rsize-1] = l2;
              continue;
            }
          }

          result.push_back(pa[j]);
          rsize++;
        }
      }

      vector<Location> endPath = leaveMap(path[2*N-1]);
      int esize = endPath.size();
      fprintf(stderr,"end path size = %d\n", esize);
      for (int i = 0; i < esize; i++) {
        Location l = endPath[i];
        result.push_back(l);
      }

      result[result.size()-1].locked = true;

      assert(result[0].locked);
      double bscore = calcCostDetail(result);
      fixPath(result);
      double ascore = calcCostDetail(result);
      assert(isValidAnswer(result));
      ret = path2answer(result);

      fprintf(stderr,"Before = %f, After = %f\n", bscore, ascore);

      return ret;
    }

    Location getLeaveLocation(double y, double x) {
      int yi = floor(y);
      int xi = floor(x);
      double minCost = DBL_MAX;
      Location best;
      Location base(y, x);

      if (yi == 0) {
        Location l(0.0001, x);
        double cost = base.dist(l);

        if (minCost > cost) {
          minCost = cost;
          best = l;
        }
      }
      if (xi == 0) {
        Location l(y, 0.0001);
        double cost = base.dist(l);

        if (minCost > cost) {
          minCost = cost;
          best = l;
        }
      }
      if (yi == S-1) {
        Location l(S-EPS/2, x);
        double cost = base.dist(l);

        if (minCost > cost) {
          minCost = cost;
          best = l;
        }
      }
      if (xi == S-1) {
        Location l(y, S-EPS/2);
        double cost = base.dist(l);

        if (minCost > cost) {
          minCost = cost;
          best = l;
        }
      }

      return best;
    }

    void fixPath(vector<Location> &path) {
      int psize = path.size();

      double timeLimit = MAX_TIME;
      double currentTime = getTime(g_startCycle);
      ll tryCount = 0;
      ll i = 0;

      while (currentTime < timeLimit) {
        i++;
        int index = i % psize;
        Location l = path[index];
        if (l.locked) continue;
        Location bl = path[index-1];
        Location al = path[index+1];
        Location temp = l;
        double s1 = costSeg(path[index-1], path[index]) + costSeg(path[index], path[index+1]);

        int d1 = xor128()%40;
        int d2 = xor128()%40;
        l.y += 0.0025 * d1 - 0.05;
        l.x += 0.0025 * d2 - 0.05;

        l.update();
        if (l.y < 0 || l.y >= S || l.x < 0 || l.x >= S) continue;
        if (l.y - l.yi < EPS) continue;
        if (l.x - l.xi < EPS) continue;
        if ((l.yi+1) - l.y < EPS) continue;
        if ((l.xi+1) - l.x < EPS) continue;
        if (bl.manhattan(l) >= 2 || al.manhattan(l) >= 2) continue;

        path[index] = l;

        double s2 = costSeg(path[index-1], path[index]) + costSeg(path[index], path[index+1]);

        if (s1 < s2) {
          path[index] = temp;
        }

        currentTime = getTime(g_startCycle);
        tryCount++;
      }

      fprintf(stderr,"tryCount fix = %lld\n", tryCount);
    }

    vector<short> cleanPath(vector<short> path) {
      vector<short> bestPath = path;
      vector<short> goodPath = path;
      updateHistory(goodPath);

      double timeLimit = MAX_TIME-1.0;
      double currentTime = getTime(g_startCycle);

      double minCost = calcCost(bestPath);
      double goodCost = minCost;

      assert(fabs(g_costHistory[2*N-1]-goodCost) < EPS);

      fprintf(stderr,"minCost = %f\n", minCost);

      ll tryCount = 0;

      int c1, c2;
      int psize = path.size();
      double alpha = 0.99999;
      double k = updateK(minCost);

      double T = 10000.0;
      double t = 1.0;
      double cost;
      ll R = 1000000;

      int nc = 0;

      while (currentTime < timeLimit) {
        do {
          c1 = xor128() % psize;
          c2 = xor128() % psize;
        } while (c1 == c2);

        int type = xor128()&3;

        if (type == 3 && (c1 > psize-3 || c2 > psize-3)) {
          continue;
        }

        switch (type) {
          case 0:
            swapObject(path, c1, c2);
            break;
          case 1:
            insertObject(path, c1, c2);
            break;
          case 2:
            crossObject(path, c1, c2);
            break;
          case 3:
            insertObject2(path, c1, c2);
            break;
        }

        if (c1 < c2 && c1 > 2) {
          cost = calcCost(path, g_costHistory[c1-2], g_itemHistory[c1-2], c1-1);
        } else if (c2 < c1 && c2 > 2) {
          cost = calcCost(path, g_costHistory[c2-2], g_itemHistory[c2-2], c2-1);
        } else {
          cost = calcCost(path);
        }

        if (cost < DBL_MAX && minCost > cost) {
          minCost = cost;
          bestPath = path;
          nc = 0;
        }

        double diffScore = goodCost - cost;

        if (goodCost > cost) {
          goodCost = cost;
          goodPath = path;
          updateHistory(goodPath);
          assert(fabs(g_costHistory[2*N-1]-goodCost) < EPS);
          nc = 0;
        } else if (cost < DBL_MAX && xor128()%R < R*exp(diffScore/(T*k))) {
          goodCost = cost;
          goodPath = path;
          updateHistory(goodPath);
          if (fabs(g_costHistory[2*N-1]-goodCost) > EPS) {
            fprintf(stderr,"costA = %f, costB = %f\n", g_costHistory[2*N-1], goodCost);
          }
          assert(fabs(g_costHistory[2*N-1]-goodCost) < EPS);
          nc++;
        } else {
          nc++;
          switch (type) {
            case 0:
              swapObject(path, c1, c2);
              break;
            default:
              path = goodPath;
              break;
          }
        }

        currentTime = getTime(g_startCycle);
        tryCount++;
        T *= alpha;
        T = max(t, T);

        if (tryCount % 20000 == 0) {
          k = updateK(minCost);
          t = updateT(timeLimit-currentTime);

          if (g_debug && tryCount % 500000 == 0) {
            fprintf(stderr,"goodCost = %f, minCost = %f, T = %f\n", goodCost, minCost, T);
            cerr.flush();
          }
        }
      }

      fprintf(stderr,"tryCount = %lld, minCost = %f, goodCost = %f\n", tryCount, minCost, goodCost);

      return bestPath;
    }

    double updateK(double minCost) {
      if (minCost > 1000) {
        return 1.0;
      } else if (minCost > 500) {
        return 0.9;
      } else {
        return 0.8;
      }
    }

    double updateT(double remainTime) {
      if (remainTime < 1.0) {
        return 0.1;
      } else if (remainTime < 2.0) {
        return 0.5;
      } else if (remainTime < 3.0) {
        return 1.0;
      } else if (remainTime < 4.0) {
        return 2.0;
      } else if (remainTime < 6.0) {
        return 3.0;
      } else {
        return 5.0;
      }
    }

    void swapObject(vector<short> &path, int c1, int c2) {
      short temp = path[c1];
      path[c1] = path[c2];
      path[c2] = temp;
    }

    void insertObject(vector<short> &path, int c1, int c2) {
      short temp = path[c1];
      path.erase(path.begin()+c1);
      path.insert(path.begin()+c2, temp);
    }

    void insertObject2(vector<short> &path, int c1, int c2) {
      short temp = path[c1];
      short temp2 = path[c1+1];
      path.erase(path.begin()+c1);
      path.erase(path.begin()+c1);
      path.insert(path.begin()+c2, temp);
      path.insert(path.begin()+c2, temp2);
    }

    void crossObject(vector<short> &path, int c1, int c2) {
      int i = min(c1, c2);
      int j = max(c1, c2);

      while (i < j) {
        short temp = path[i];
        path[i] = path[j];
        path[j] = temp;
        i++;
        j--;
      }
    }

    vector<double> path2answer(vector<Location> &path) {
      vector<double> answer;

      int psize = path.size();
      for (int i = 0; i < psize; i++) {
        Location l = path[i];
        answer.push_back(l.x);
        answer.push_back(l.y);
        //fprintf(stderr,"%d: y = %f, x = %f\n", i, l.y, l.x);
      }

      return answer;
    }

    vector<Location> createPath(int from, int to) {
      vector<Location> path;
      Object *fromObj = getObject(from);
      Object *toObj = getObject(to);

      Node root(fromObj->nid, 0.0, fromObj->y, fromObj->x);

      priority_queue<Node, vector<Node>, greater<Node> > pque;
      pque.push(root);
      vector<bool> checkList(V);

      while (!pque.empty()) {
        Node node = pque.top(); pque.pop();
        int y = node.yi;
        int x = node.xi;

        if (checkList[node.cid]) continue;
        checkList[node.cid] = true;

        if (node.cid == toObj->nid) {
          int isize = node.ids.size();

          path.push_back(Location(fromObj->y, fromObj->x, true));

          for (int i = 0; i < isize-1; i++) {
            short id = node.ids[i];
            Location l = nid2coord(id);
            if (id == 0) l.locked = true;
            path.push_back(l);
          }

          path.push_back(Location(toObj->y, toObj->x, true));
          break;
        }

        for (int i = 0; i < 4; i++) {
          short ny = y + DY[i];
          short nx = x + DX[i];
          if (isOutside(ny, nx)) continue;

          short nid = ny*S + nx;
          Node next = node.dup();
          next.cid = nid;
          next.beforeDirect = i;
          next.length = node.length + 1;
          next.yi = ny;
          next.xi = nx;
          next.y = ny + 0.5;
          next.x = nx + 0.5;

          double costB = costSeg(Location(node.y, node.x), Location(next.y, next.x));
          double costC = pow(g_fieldCost[y][x]-g_fieldCost[ny][nx], 2);
          next.cost += costB + costC;
          next.ids.push_back(nid);

          next.step[0] = i;
          next.fc[0] = g_fieldCost[ny][nx];

          updateNodeCost(next);
          pque.push(next);
        }
      }

      return path;
    }

    vector<short> createFirstPathNN() {
      vector<short> ret;
      vector<bool> checkList(V);

      ret.push_back(g_objectList[0].id);
      checkList[g_objectList[0].id] = true;

      for (int i = 0; i < 2*N-1; i++) {
        short from = ret[i];
        double minCost = DBL_MAX;
        short minId = -1;

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

    vector<short> createFirstPathFI() {
      vector<short> ret;
      int itemCount[2*N];
      memset(itemCount, 0, sizeof(itemCount));
      vector<bool> checkList(V);

      double minDist = DBL_MAX;
      int minId = -1;
      for (int i = 0; i < N; i++) {
        Item *item = getItem(i);
        double dist = calcDist(item->y, item->x, S/2.0, S/2.0) - g_leaveCost[i];

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

            double minCost = DBL_MAX;
            minId = -1;

            for (int k = 0; k < ret.size(); k++) {
              double cost = g_pathCost[target->id][ret[k]];

              if (minCost > cost) {
                minCost = cost;
                minId = target->id;
              }
            }

            if (maxDist < minCost) {
              maxDist = minCost;
              maxId = minId;
            }
          }

          double minCost = DBL_MAX;
          int index = -1;

          int rsize = ret.size();
          for (int k = 0; k < rsize; k++) {
            if (itemCount[k] == 0) continue;
            vector<short> temp = ret;
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

            double minCost = DBL_MAX;
            minId = -1;
            for (int k = 0; k < ret.size(); k++) {
              double cost = g_pathCost[item->id][ret[k]];

              if (minCost > cost) {
                minCost = cost;
                minId = item->id;
              }
            }

            if (maxDist < minCost) {
              maxDist = minCost;
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
      Node root(obj->nid, 0.0, obj->y, obj->x);
      pque.push(root);
      vector<bool> checkList(V);

      while (!pque.empty()) {
        Node node = pque.top(); pque.pop();
        short y = node.yi;
        short x = node.xi;

        if (checkList[node.cid]) continue;
        checkList[node.cid] = true;

        if (isBorder(y, x)) {
          int isize = node.ids.size();

          for (int i = 0; i < isize-1; i++) {
            short id = node.ids[i];
            double cy = id / S + 0.5;
            double cx = id % S + 0.5;
            path.push_back(Location(cy, cx));
          }

          path.push_back(getLeaveLocation(node.y, node.x));

          break;
        }

        for (int i = 0; i < 4; i++) {
          short ny = y + DY[i];
          short nx = x + DX[i];
          if (isOutside(ny, nx)) continue;

          short nid = ny*S + nx;
          Node next = node.dup();
          next.cid = nid;
          next.beforeDirect = i;
          next.length = node.length + 1;
          next.yi = ny;
          next.xi = nx;
          next.y = ny + 0.5;
          next.x = nx + 0.5;

          double costB = costSeg(Location(node.y, node.x), Location(next.y, next.x));
          double costC = pow(g_fieldCost[y][x]-g_fieldCost[ny][nx], 2);
          next.cost += costB + costC;
          next.ids.push_back(nid);

          next.step[0] = i;
          next.fc[0] = g_fieldCost[ny][nx];

          updateNodeCost(next);
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

    double calcCost(const vector<short> &path, double fc = 0.0, int fic = 0, int offset = 0) {
      int psize = 2*N;
      int itemCount = fic;
      double cost = (fc > 0.0)? fc : g_leaveCost[path[0]];

      for (int i = offset; i < psize-1; i++) {
        itemCount += g_objectType[path[i]];

        if (itemCount < 0 || itemCount > g_capacity) {
          return DBL_MAX;
        }

        cost += g_pathCost[path[i]][path[i+1]];
      }

      cost += g_leaveCost[path[psize-1]];

      return cost;
    }

    void updateNodeCost(Node &node) {
      if (node.length <= 2 || node.update) return;

      int m0 = node.step[0] % 2;
      int m1 = node.step[1] % 2;
      int m2 = node.step[2] % 2;

      if (m0 != m1 && m1 != m2) {
        node.cost -= 0.2*max(node.fc[1], node.fc[2]);
      }
      if (m0 == m1 && m1 != m2) {
        node.cost -= 0.2*node.fc[2];
      }
    }

    void updateHistory(const vector<short> &path) {
      int psize = 2*N;
      int itemCount = 0;
      double cost = g_leaveCost[path[0]];

      for (int i = 0; i < psize-1; i++) {
        itemCount += g_objectType[path[i]];

        if (itemCount < 0 || itemCount > g_capacity) {
          assert(false);
        }

        cost += g_pathCost[path[i]][path[i+1]];
        g_costHistory[i] = cost;
        g_itemHistory[i] = itemCount;
      }

      cost += g_leaveCost[path[psize-1]];
      g_costHistory[psize-1] = cost;
    }

    bool isValidAnswer(vector<Location> &path) {
      int psize = path.size();

      for (int i = 0; i < psize; i++) {
        Location l = path[i];

        if (isNearInnerBorder(l)) return false;
      }

      return true;
    }

    bool isValidPath(vector<short> &path) {
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
      return &g_objectList[id];
    }

    Target *getTarget(int id) {
      return &g_objectList[N+id];
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

    inline bool isOutside(int y, int x) {
      return (y < 0 || S <= y || x < 0 || S <= x);
    }
};


// -------8<------- end of solution submitted to the website -------8<-------
template<class T> void getVector(vector<T>& v) { for (int i = 0; i < v.size(); ++i) cin >> v[i];}
int main() {
  g_debug = true;
  TerrainCrossing tc; int L, M, capacity;
  cin >> M;
  vector<string> map(M); getVector(map);
  cin >> L;
  vector<double> locations(L); getVector(locations);
  cin >> capacity;
  vector<double> ret = tc.getPath(map, locations, capacity);
  cout << ret.size() << endl;
  for (int i = 0; i < ret.size(); ++i) {
    printf("%10.6f\n", ret[i]);
  }
  cout.flush();
}
