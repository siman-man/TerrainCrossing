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
const int ITEM = 0;
const int TARGET = 1;

const ll CYCLE_PER_SEC = 2400000000;
double MAX_TIME = 10.0;

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

struct Result {
  double cost;
  bool valid;

  Result(double cost = DBL_MAX, bool valid = true) {
    this->cost = cost;
    this->valid = valid;
  }
};

struct BeamNode {
  double cost;
  int cid;
  int itemCount;
  bitset<MAX_N> checkList;
  vector<int> path;

  BeamNode() {
    this->checkList = 0;
    this->cid = -1;
    this->cost = 0.0;
    this->itemCount = 0;
  }

  bool operator >(const BeamNode &bn) const {
    return cost > bn.cost;
  }
};

struct Node {
  int cid;
  double y;
  double x;
  int yi;
  int xi;
  double cost;
  int beforeY;
  int beforeX;
  int beforeDirect;
  int length;
  vector<int> ids;

  Node(int cid, double cost, double y, double x) {
    this->cid = cid;
    this->cost = cost;
    this->y = y;
    this->x = x;
    this->yi = floor(y);
    this->xi = floor(x);
    this->beforeDirect = -1;
    this->length = 0;
  }

  Node dup() {
    Node node(cid, cost, y, x);
    node.ids = ids;
    node.beforeDirect = beforeDirect;

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
  double leaveCost;

  Object(int id, int type, double y, double x) {
    this->id = id;
    this->type = type;
    this->y = y;
    this->x = x;
    this->leaveCost = DBL_MAX;
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

      priority_queue<Node, vector<Node>, greater<Node> > pque;
      pque.push(root);

      map<int, bool> checkList;

      while (!pque.empty()) {
        Node node = pque.top(); pque.pop();

        int y = node.yi;
        int x = node.xi;

        if (checkList[node.cid]) continue;
        checkList[node.cid] = true;

        if (y == 0 || y == S-1 || x == 0 || x == S-1) {
          obj->leaveCost = min(obj->leaveCost, node.cost);
        }

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
          next.beforeDirect = i;
          next.length = node.length + 1;
          next.yi = ny;
          next.xi = nx;
          next.beforeY = y;
          next.beforeX = x;
          next.cid = nid;

          double costB, costC;
          if (node.length == 0) {
            costB = costSeg(Location(obj->y, obj->x), Location(ny+0.5, nx+0.5));
          } else {
            costB = costSeg(Location(y+0.5, x+0.5), Location(ny+0.5, nx+0.5));
          }
          costC = pow(g_fieldCost[y][x] - g_fieldCost[ny][nx], 2);

          if (node.length == 0) {
            next.cost = costB + costC;
          } else if ((node.beforeDirect % 2) != (i % 2)) {
            next.cost = node.cost + 0.9*costB + costC;
          } else {
            next.cost = node.cost + costB + costC;
          }
          //next.ids.push_back(nid);
          pque.push(next);
        }
      }
    }

    vector<double> getPath(vector<string> map, vector<double> locations, int capacity) {
      vector<double> ret;

      init(map, locations, capacity);

      fprintf(stderr,"start =>\n");

      vector<int> path1;

      if (N <= 50) {
        path1 = createFirstPathBeam();
      } else {
        path1 = createFirstPathNN();
      }

      double currentTime = getTime(g_startCycle);
      fprintf(stderr,"current time = %f\n", currentTime);

      cerr.flush();
      assert(isValidPath(path1));
      vector<int> path2 = createFirstPathFI();
      assert(isValidPath(path2));

      Result result1 = calcCost(path1);
      Result result2 = calcCost(path2);
      double score1 = result1.cost;
      double score2 = result2.cost;

      vector<int> path = (score1 < score2)? path1 : path2;

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
        int from = path[i];
        int to = path[i+1];

        vector<Location> pa = createPath(from, to);

        for (int j = 0; j < pa.size(); j++) {
          Location l1 = result[rsize-1];
          Location l2 = pa[j];

          if (!l1.locked) {
            int nid1 = l1.yi*S + l1.xi;
            int nid2 = l2.yi*S + l2.xi;

            if (nid1 == nid2) {
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
      fixPath(result);
      assert(isValidAnswer(result));
      ret = path2answer(result);

      return ret;
    }

    void fixPath(vector<Location> &path) {
      int psize = path.size();

      double timeLimit = 10.0;
      double currentTime = getTime(g_startCycle);
      ll tryCount = 0;

      while (currentTime < timeLimit) {
        int index = xor128() % psize;
        Location l = path[index];
        if (l.locked) continue;
        Location temp = l;
        double s1 = costSeg(path[index-1], path[index]) + costSeg(path[index], path[index+1]);

        int d1 = xor128()%40;
        int d2 = xor128()%40;
        l.y += 0.0025 * d1 - 0.05;
        l.x += 0.0025 * d2 - 0.05;
        int nid = floor(l.y)*S + floor(l.x);

        if (l.y - l.yi < EPS) continue;
        if (l.x - l.xi < EPS) continue;
        if ((l.yi+1) - l.y < EPS) continue;
        if ((l.xi+1) - l.x < EPS) continue;

        if (nid != l.nid) continue;
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

    vector<int> cleanPath(vector<int> path) {
      vector<int> bestPath = path;
      vector<int> goodPath = path;

      double timeLimit = 8.0;
      double currentTime = getTime(g_startCycle);

      Result result = calcCost(bestPath);
      double minCost = result.cost;
      double goodCost = minCost;

      fprintf(stderr,"minCost = %f\n", minCost);

      ll tryCount = 0;
      ll invalidCount = 0;

      int c1, c2;
      int psize = path.size();

      double T = 10000.0;
      double alpha = 0.999;
      double penalty = 0.0;
      double remainTime;

      while (currentTime < timeLimit) {
        remainTime = timeLimit - currentTime;
        do {
          c1 = xor128() % psize;
          c2 = xor128() % psize;
        } while (c1 == c2);

        int type = xor128()%4;

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

        result = calcCost(path, penalty);

        if (!result.valid) {
          invalidCount++;
        }

        if (minCost > result.cost && result.valid) {
          minCost = result.cost;
          bestPath = path;
        }

        if (goodCost > result.cost) {
          goodCost = result.cost;
          goodPath = path;
        } else {
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

        if (remainTime < 2.0 && penalty < 100) {
          //fprintf(stderr,"goodCost = %f, minCost = %f\n", goodCost, minCost);
          penalty = 10000.0;
          result = calcCost(goodPath, penalty);
          goodCost = result.cost;
        } else if (remainTime < 4.0 && penalty < 10) {
          penalty = 30.0;
          result = calcCost(goodPath, penalty);
          goodCost = result.cost;
        } else if (remainTime < 5.0 && penalty < 5) {
          //fprintf(stderr,"goodCost = %f, minCost = %f\n", goodCost, minCost);
          penalty = 5.0;
          result = calcCost(goodPath, penalty);
          goodCost = result.cost;
        }
      }

      fprintf(stderr,"(%lld/%lld), tryCount = %lld, minCost = %f, goodCost = %f\n", invalidCount, tryCount, tryCount, minCost, goodCost);

      return bestPath;
    }

    void swapObject(vector<int> &path, int c1, int c2) {
      int temp = path[c1];
      path[c1] = path[c2];
      path[c2] = temp;
    }

    void insertObject(vector<int> &path, int c1, int c2) {
      int temp = path[c1];
      path.erase(path.begin()+c1);
      path.insert(path.begin()+c2, temp);
    }

    void insertObject2(vector<int> &path, int c1, int c2) {
      int temp = path[c1];
      int temp2 = path[c1+1];
      path.erase(path.begin()+c1);
      path.erase(path.begin()+c1);
      path.insert(path.begin()+c2, temp);
      path.insert(path.begin()+c2, temp2);
    }

    void crossObject(vector<int> &path, int c1, int c2) {
      int i = min(c1, c2);
      int j = max(c1, c2);

      while (i < j) {
        int temp = path[i];
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
            Location l = nid2coord(id);
            if (id == 0) l.locked = true;
            path.push_back(l);
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
          next.beforeDirect = i;
          next.length = node.length + 1;

          double costB, costC;

          if (node.length == 0) {
            costB = costSeg(Location(fromObj->y, fromObj->x), Location(ny+0.5, nx+0.5));
          } else {
            costB = costSeg(Location(y+0.5, x+0.5), Location(ny+0.5, nx+0.5));
          }
          costC = pow(g_fieldCost[y][x] - g_fieldCost[ny][nx], 2);

          next.cost += costB + costC;
          next.ids.push_back(nid);
          pque.push(next);
        }
      }

      return path;
    }

    vector<int> createFirstPathBeam(int BEAM_WIDTH = 2000) {
      BeamNode root;
      queue<BeamNode> que;
      que.push(root);
      BeamNode cand;

      for (int i = 0; i < 2*N; i++) {
        priority_queue<BeamNode, vector<BeamNode>, greater<BeamNode> > pque;

        while (!que.empty()) {
          BeamNode node = que.front(); que.pop();
          int itemCount = node.itemCount;

          for (int oid = 0; oid < 2*N; oid++) {
            if (node.checkList[oid]) continue;
            if (node.cid == oid) continue;
            if (itemCount == g_capacity && oid < N) continue;
            if (itemCount == 0 && oid >= N) continue;

            Object *obj = getObject(oid);

            cand.cid = oid;
            cand.path = node.path;
            cand.checkList = node.checkList;
            cand.checkList.set(oid);
            cand.path.push_back(oid);

            if (i == 0) {
              cand.cost = obj->leaveCost;
            } else if (i == 2*N-1) {
              cand.cost = node.cost + g_pathCost[node.cid][oid] + obj->leaveCost;
            } else {
              cand.cost = node.cost + g_pathCost[node.cid][oid];
            }
            if (obj->isItem()) {
              cand.itemCount = itemCount+1;
            } else {
              cand.itemCount = itemCount-1;
            }

            pque.push(cand);
          }
        }

        for (int k = 0; k < BEAM_WIDTH && !pque.empty(); k++) {
          BeamNode node = pque.top(); pque.pop();
          que.push(node);
        }
      }

      assert(!que.empty());

      BeamNode node = que.front();

      assert(node.path.size() == 2*N);
      return node.path;
    }

    vector<int> createFirstPathNN() {
      vector<int> ret;
      map<int, bool> checkList;

      ret.push_back(g_objectList[0].id);
      checkList[g_objectList[0].id] = true;

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
        //double dist = item->leaveCost;
        double dist = calcDist(item->y, item->x, S/2.0, S/2.0) - item->leaveCost;

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

    Result calcCost(vector<int> &path, double penalty = 0) {
      int psize = path.size();
      int itemCount = 0;
      double score = g_objectList[path[0]].leaveCost;
      Result result;

      for (int i = 0; i < psize-1; i++) {
        int oid = path[i];
        Object *obj = getObject(oid);

        itemCount = (obj->isItem())? itemCount+1 : itemCount-1;
        if (itemCount < 0 || itemCount > g_capacity) {
          score += penalty;
          result.valid = false;
        }

        score += g_pathCost[path[i]][path[i+1]];
      }

      score += g_objectList[path[psize-1]].leaveCost;
      result.cost = score;

      return result;
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
      return &g_objectList[id];
    }

    Target *getTarget(int id) {
      return &g_objectList[N+id];
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
