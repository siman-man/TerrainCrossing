Problem Statement
問題文
You are given a square map of S x S cells. Each cell is a square area of certain terrain type, with side of length 1. 
あなたにはSxSのセルで構成されたマップが与えられます。 各セルには地形が定義されています。 1辺の長さは1です。
Terrain types are denoted with digits from 0 to 9, which describe the cost of passing through this terrain, 
地形の種類は数字の0から9で表されます。 この数字はこの地形を通過するときのコストを表しています。
per unit: 0 is the easiest to cross, and 9 is the hardest. The cost of passing through the square is calculated 
各ユニットにおいて0は一番低いコストで、9が一番高いコストです。各地形を通過するコストはユークリッド距離 * その地形のコストから算出されます。
by multiplying the Euclidean distance travelled with the terrain type. 
When you cross a border between two terrain types, you incur additional cost of (difference of terrain types)2. 
もし、あなたが2つの地形間の移動を行う場合は追加のコストとして2を払います。


There are N items located on the map, and N target locations to which these items have to be delivered. 
マップ上にはN個のアイテムが配置されています。 そしてN個の目的にそれらのアイテムを運びます
Items are identical, so each item can be delivered to any location, but each location must have exactly one item delivered to it. 
各アイテムにはIDが付いており、運ぶ場所についてはどこでもかまいません。 ただし、各目的地に運ぶアイテムは重複しないようにしてください。
You can carry at most capacity items at once. You automatically pick up an item if you stop within 0.001 from it and still 
あなたは一度に可能な限りのアイテムを持つことが出来ます。 アイテムから0.001の距離にとまると自動的にアイテムが回収されます。 そして可能なかぎりそれを保有します。
have capacity to carry it, and you automatically drop off an item at a target location if you stop within 0.001 from it while 
                           また目的地の近くに到達すると自動的にアイテムがその場所に配置されます。そしてその場所にはもうアイテムが配置できなくなります。
carrying at least one item and no item has been delivered to this location yet. 

Your task is to enter the map at any place along its border, pick up all items and deliver them to target locations, 
あなたの目的は任意の場所からマップに入り、全てのアイテムを目的地に届けた後にマップから出ることです。
and exit the map at any place along its border.

Implementation
Your code must implement one method getPath(vector <string> terrain, vector <double> locations, int capacity):
あなたはgetPathメソッドを実装します
terrain gives the map of terrain types in the area. terrain[i][j] describes the type of terrain in the square with 
terrain引数にはマップの地形データが与えられます。 terrain[i][j]にはその地形の種類が記録されています。
coordinates [j, j+1] x [i, i+1]. Characters '0'..'9' represent types 0..9.
0から9の範囲で表されます。

locations gives the list of items and target locations for them. For N items and N target locations, locations will 
locationsにはアイテムと目的地の座標が与えられています。
contain 4*N elements. First 2*N elements will describe positions of items: ith item is located at 
locationsは4*Nのサイズで構成されており、最初の2*Nまではアイテムの座標、次の2*Nの部分には目的地の座標が格納されています。
coordinates (locations[2*i], locations[2*i+1]). Next 2*N elements will describe target locations: jth target location has coordinates (locations[2*N+2*j], locations[2*N+2*j+1]).

capacity gives the maximum number of items you can carry at once.
capacityはあなたが一度に持ち運べるアイテムの最大量を表します。
The return from this method will describe a path you want to take. The path is a sequence of points within the map 
返す値は、経路を返します。
connected by segments; i-th point of the path has coordinates (return[2*i], return[2*i+1]). The path must satisfy the following conditions:
各ポイントを接続します。 経路は以下の条件をみたす必要があります。

The path must have between 2 and 4 * S * S * N points, inclusive.
パスに含まれるポイントの数は2個から4 * S * S * N個の数におさめてください。
Each point of the path must be within the map, i.e. both coordinates must be between 0 and S.
各ポイントの座標はマップの内に修めるようにしてください。
The first and the last points of the path must be within 0.001 from the outer border of the map.
最初の座標と最後の座標についてはマップの端から0.001の距離に収まるポイントを指定してください。
All points of the path must be at least 0.001 away from internal borders between the cells of the map (even if cells on both sides of the border are of the same terrain type).
Consecutive points of the path must be at least 0.001 away from each other. (Euclidean distance).
連続した点については必ず0.001以上離すようにしてください。
Each segment of the path can cross at most one boundary between cells of the map, i.e. the Manhattan distance 
各経路において最大1までのセルを飛ばすことが出来ます。
between cells to which consecutive points of the path belong can be at most 1.

After the path is walked, all items must be picked up and all target locations must have an item delivered to them.
経路を全て探索し終えた後はアイテムを全て回収し、全ての目的地に届けておく必要があります。

Example
An example solution for seed 1 can be seen in the image. The items is shown with green dots and the delivery locations with red dots. The example path is shown in blue.

Scoring
For each test case we will calculate your raw and normalized scores. If your solution returned invalid path 
各テストケースにおいてあなたのスコアは正規化されます。
(some points outside the map or too close to the borders, not all items picked up and delivered etc.), raw score will be 0. 
もし不正な解答を行った場合にはスコアは0となります。
Otherwise, raw score will be the total cost of the path, calculated as sum of costs of passing through terrain done 
on each segment of the path (for each part of the path which passes through terrain of the same type its cost is the length of the part, 
multiplied by cost of this terrain type) and additional costs of changing terrain types. The normalized score for each 
test is 1,000,000.0 * BEST / YOUR, where BEST is the lowest raw score currently obtained on this test case 
(considering only the last submission from each competitor). Finally, your total score is equal to the arithmetic average of normalized scores on all test cases.

Tools
An offline tester is available here. You can use it to test/debug your solution locally. You can also check its source code for exact implementation of 
test case generation and score calculation. That page also contains links to useful information and sample solutions in several languages.

Definition
      
Class:  TerrainCrossing
Method: getPath
Parameters: vector <string>, vector <double>, int
Returns:  vector <double>
Method signature: vector <double> getPath(vector <string> terrain, vector <double> locations, int capacity)
(be sure your method is public)

Notes
- The time limit is 10 seconds per test case (this includes only the time spent in your code). The memory limit is 1024 megabytes.
制限時間は10秒数です。メモリの制限は1GBです。
- There is no explicit code size limit. The implicit source code size limit is around 1 MB (it is not advisable to submit codes of 
コードのサイズは1MBを上限とします。
size close to that or larger). Once your code is compiled, the binary size should not exceed 1 MB.
- The compilation time limit is 30 seconds. You can find information about compilers that we use and compilation options here.
コンパイル時間は30秒以内で終わるようにしてください。
- There are 10 example test cases and 100 full submission (provisional) test cases.
exampleは10ケースでfull submissionは100ケースで行われます。
- The match is rated.
このコンテストはレートが変動します。

Constraints
- The size of the map S will be between 10 and 50, inclusive.
マップの大きさは10から50です
- The number of terrain types T will be between 2 and 10, inclusive.
地形の種類は2から10の間です
- The number of items N will be between 5 and S*S/10, inclusive.
アイテムの数は5からS*S/10の範囲で与えられています
- The carrying capacity capacity will be between 1 and 10, inclusive.
持ち運べるアイテムの数は1から10の間です


# 与えられてるもの

* 地形情報
* アイテムと目的地の座標
** アイテムと目的地の座標はそれぞれランダムに与えられている（2つの距離は最低でも0.003以上は離れている

# 類似問題

* 巡回セールスマン問題


# 条件

* 各経路間の距離は1以内に収めること

# 考察

* 1つのアイテムを回収してそれを目的地に届けるまでに4*S*Sサイズの行動が可能 （実質無限
* 各セル間のコストはどうするか、その都度計算すると時間がもったいない
* アイテム -> 目的地 -> アイテム -> 目的地の順番で取得していけば最低限の動くコードが完成する

# 戦略

* 初期解を生成
* 焼きなまして経路を修正
* 最後に経路を出力用に変換して返す

## 生成する経路は以下の条件を満たす

* アイテムを持っていない状態で目的地に到達しない
* アイテムの所持数が上限の状態でアイテムに近づかない
