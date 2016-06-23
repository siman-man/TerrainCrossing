/*
Change Log:
2016-06-21 : Intial release
*/

import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.io.*;
import java.lang.Math.*;
import java.security.*;
import java.util.*;
import javax.swing.*;

class Location {
    public static final double eps = 1E-3;
    public double x, y;
    public int xi, yi;
    public Location() {};
    public Location(double x1, double y1) {
        x = x1;
        y = y1;
        xi = (int)(Math.floor(x));
        yi = (int)(Math.floor(y));
    }
    public double dist(Location other) {
        return Math.sqrt(Math.pow(x - other.x, 2) + Math.pow(y - other.y, 2));
    }
    public boolean near(Location other, double d) {
        return dist(other) <= d;
    }
    public int manhattan(Location other) {
        return Math.abs(xi - other.xi) + Math.abs(yi - other.yi);
    }
    public String toString() {
        return x + " " + y;
    }
}

public class TerrainCrossingVis {
    static int maxSize = 50, minSize = 10;
    static int maxTypes = 10, minTypes = 2;
    static int minItems = 5;
    static int maxCap = 10, minCap = 1;

    int maxItems;           // maximum number of items

    int S;                  // size of the terrain to cross
    int T;                  // number of terrain types used (0..T-1)
    String[] map;           // map of terrain types

    int N;                  // number of items to be picked up and delivered
    // locations of items: first half of the list (N) gives current locations of items, second half (N) - drop-off locations
    Location[] loc;         // x-coordinate corresponds to j-index (column), y = i (row), grow in the same direction
    boolean[] locused;      // whether item has been picked up/dropped off at this location
    int C;                  // max carrying capacity of the traveller

    Location[] path;
    // -----------------------------------------
    double getRandInsideCell(SecureRandom r1) {
        return r1.nextInt(S) + Location.eps + r1.nextDouble() * (1 - 2 * Location.eps);
    }
    // -----------------------------------------
    void generate(String seedStr) {
    try {
        // generate test case
        SecureRandom r1 = SecureRandom.getInstance("SHA1PRNG");
        long seed = Long.parseLong(seedStr);
        r1.setSeed(seed);
        // フィールドのサイズを決める
        S = r1.nextInt(maxSize - minSize + 1) + minSize;
        // アイテムの最大値を決める
        maxItems = S*S/10;
        T = r1.nextInt(maxTypes - minTypes + 1) + minTypes;
        N = r1.nextInt(maxItems - minItems + 1) + minItems;
        C = r1.nextInt(maxCap - minCap + 1) + minCap;
        // 
        if (seed == 1) {
            S = 5;
            T = 3;
            N = 4;
            C = 2;
        } else if (seed <= 3) {
            S = minSize * (int)seed;
            T = minTypes + 2 * (int)(seed - 1);
            N = minItems * (int)seed;
            C = minCap * (int)seed;
        } else if (seed == 4) {
            S = maxSize;
            T = maxTypes;
            N = S*S/10;
            C = maxCap;
        }

        // generate the map of terrain types

        // start with a completely random map
        char[][] m = new char[S][S];
        double[][] raw = new double[S+2][S+2];
        for (int i = 0; i < S+2; ++i)
        for (int j = 0; j < S+2; ++j) {
            raw[i][j] = r1.nextDouble();
        }
        // smooth terrain
        int nk = r1.nextInt(T);
        if (seed == 1) nk = 0;
        for (int k = 0; k < nk; ++k)
        for (int i = 1; i <= S; ++i)
        for (int j = 1; j <= S; ++j) {
            raw[i][j] = 0.6*raw[i][j] + 0.1*(raw[i-1][j]+raw[i+1][j]+raw[i][j-1]+raw[i][j+1]);
        }
        // Find maximum
        double rawmax = 0;
        for (int i = 1; i <= S; ++i)
            for (int j = 1; j <= S; ++j)
                rawmax = Math.max(rawmax, raw[i][j]);
        // Convert to terrain type
        for (int i = 0; i < S; ++i)
        for (int j = 0; j < S; ++j) {
            m[i][j] = (char)('0' + (int)(raw[i+1][j+1]*(T-1)/rawmax));
        }

        map = new String[S];
        for (int i = 0; i < S; ++i)
            map[i] = new String(m[i]);

        // generate pick up/drop off information
        locused = new boolean[2 * N];
        loc = new Location[2 * N];
        for (int i = 0; i < 2 * N; ++i) {
            boolean ok;
            do {
                ok = true;
                // all locations are at least eps away from border, so you can't visit location without entering the cell
                loc[i] = new Location(getRandInsideCell(r1), getRandInsideCell(r1));
                // verify that no two locations are within 3 * eps (so any other location can be near at most one of them)
                for (int j = 0; j < i; ++j)
                    // 他の座標との距離が3 * 0.001以内にある場合は設置を行わない
                    if (loc[i].near(loc[j], 3 * Location.eps)) {
                        ok = false;
                        break;
                    }
            } while (!ok);
        }

        if (debug) {
            System.out.println("Size of terrain S = " + S);
            System.out.println("Number of terrain types T = " + T);
            System.out.println("Terrain: ");
            for (String st : map) {
                System.out.println(st);
            }
            System.out.println("Capacity C = " + C);
            System.out.println("Number of items N = " + N);
            System.out.println("Locations:");
            for (Location l : loc) {
                System.out.println(l.toString());
            }
        }

        if (vis) {
            W = S * SZ + 40;
            H = S * SZ + 40;
            jf.setSize(W, H);
            jf.setVisible(true);
            v.repaint();
        }
    }
    catch (Exception e) {
        System.err.println("An exception occurred while generating test case.");
        e.printStackTrace();
    }
    }
    // 座標が画面外に出ているかどうかを調べる
    boolean isNearOuterBorder(Location l) {
        return !(l.x >= Location.eps && l.x <= S - Location.eps && l.y >= Location.eps && l.y <= S - Location.eps);
    }
    // セルの内部に座標があるかどうかを調べる
    boolean isNearInnerBorder(Location l) {
        double inCellX = l.x - l.xi;
        double inCellY = l.y - l.yi;
        return (l.xi > 0 && inCellX < Location.eps ||
                l.yi > 0 && inCellY < Location.eps ||
                l.xi < S - 1 && 1 - inCellX < Location.eps ||
                l.yi < S - 1 && 1 - inCellY < Location.eps);
    }
    // -----------------------------------------
    int terrainType(Location l) {
        return map[l.yi].charAt(l.xi) - '0';
    }
    // -----------------------------------------
    double scoreSeg(Location l1, Location l2) {
        // if two locations are in the same cell, score = length of segment * terrain type
        if (l1.manhattan(l2) == 0) {
            return l1.dist(l2) * terrainType(l1);
        }
        // otherwise two locations are in adjacent cells
        // score = length of segment in cell 1 * terrain type 1
        //       + length of segment in cell 2 * terrain type 2
        //       + square of level differences between terrain types 1 and 2 (direction of movement doesn't count)
        int t1 = terrainType(l1), t2 = terrainType(l2);
        double score = Math.pow(t1 - t2, 2);
        // figure out the line which intersects the segment
        double x0, y0;
        if (l1.yi == l2.yi) {
            // cells different along x - intersect segment with vertical line
            x0 = Math.max(l1.xi, l2.xi);
            y0 = l1.y + (l2.y - l1.y) * (x0 - l1.x) / (l2.x - l1.x);
        } else {
            // intersect segment with horizontal line
            y0 = Math.max(l1.yi, l2.yi);
            x0 = l1.x + (l2.x - l1.x) * (y0 - l1.y) / (l2.y - l1.y);
        }
        // find intersection point (as location)
        Location intersection = new Location(x0, y0);
        /*if (debug) {
            System.out.println(l1 + " x " + l2 + " = " + intersection);
        }*/
        return score + l1.dist(intersection) * t1 + l2.dist(intersection) * t2;
    }
    // -----------------------------------------
    public double runTest(String seed) {
    try {
        generate(seed);

        double[] locarg = new double[4 * N];
        for (int i = 0; i < 2 * N; ++i) {
            locarg[2 * i] = loc[i].x;
            locarg[2 * i + 1] = loc[i].y;
        }

        double[] pret = getPath(map, locarg, C);

        // convert return value to a path
        if (pret == null || pret.length == 0) {
            addFatalError("Failed to get result from getPath.");
            return 0.0;
        }
        if (pret.length % 2 == 1) {
            addFatalError("Return from getPath must have even number of elements.");
            return 0.0;
        }
        int nP = pret.length / 2;
        if (nP < 2) {
            addFatalError("The path must have at least 2 points.");
            return 0.0;
        }
        // allow the path to visit each cell for each location twice
        int maxP = 4 * N * S * S;
        if (nP > maxP) {
            addFatalError("The path can have at most " + maxP + " points.");
            return 0.0;
        }
        path = new Location[nP];
        for (int i = 0; i < nP; ++i) {
            if (pret[2 * i] < 0 || pret[2 * i] > S || pret[2 * i + 1] < 0 || pret[2 * i + 1] > S) {
                addFatalError("Each point of path must be within the terrain.");
                return 0.0;
            }
            path[i] = new Location(pret[2 * i], pret[2 * i + 1]);
        }
        if (vis) {
            v.repaint();
        }

        // validate path:
        // 1.1. start and end points must be near the outer border
        if (!isNearOuterBorder(path[0])) {
            addFatalError("The start point of the path must be within " + Location.eps + " from outer border.");
            return 0.0;
        }
        if (!isNearOuterBorder(path[nP - 1])) {
            addFatalError("The end point of the path must be within " + Location.eps + " from outer border.");
            return 0.0;
        }
        // 1.2. all coordinates must be NOT near any internal border between cells (including start and end points)
        //      but can be near outer border, that's ok
        for (int i = 0; i < nP; ++i)
            if (isNearInnerBorder(path[i])) {
                addFatalError("Point " + i + " of the path cannot be within " + Location.eps + " from any inner cell border.");
                return 0.0;
            }

        // 2. each point on the path is not near previous one
        for (int i = 1; i < nP; ++i)
            if (path[i].near(path[i-1], Location.eps)) {
                addFatalError("Consecutive points of the path cannot be within " + Location.eps + " from each other.");
                return 0.0;
            }

        // 3. each segment of the path crosses at most one cell boundary
        //    = Manhattan distance between start and end cell of each segment <= 1
        for (int i = 1; i < nP; ++i)
            if (path[i].manhattan(path[i-1]) > 1) {
                addFatalError("Each segment of the path cannot cross more than one cell boundary.");
                return 0.0;
            }

        // 4. after the path is traversed, all items must be picked up and all drop off locations must get item
        int carry = 0;
        for (int i = 0; i < nP; ++i) {
            // check whether endpoint is near one of locations
            int nearLoc = -1;
            for (int j = 0; j < 2 * N; ++j)
                if (path[i].near(loc[j], Location.eps)) {
                    nearLoc = j;
                    // can't have two locations near
                    break;
                }
            if (nearLoc == -1)
                continue;
            // if this location has already been used, ignore
            if (locused[nearLoc])
                continue;
            // if this location has an item: if traveller still has capacity, pick up item, otherwise ignore
            if (nearLoc < N && carry < C) {
                carry++;
                locused[nearLoc] = true;
                if (debug) {
                    System.out.println("Picked up item at location " + nearLoc + ". Carry = " + carry);
                }
            }
            // if this location is a drop-off: if traveller carries some items, drop off one, otherwise ignore
            if (nearLoc >= N && carry > 0) {
                carry--;
                locused[nearLoc] = true;
                if (debug) {
                    System.out.println("Dropped off item at location " + nearLoc + ". Carry = " + carry);
                }
            }
        }
        // all items have been delivered = all locations have been used
        if (carry != 0) {
            addFatalError("In the end of the path you must carry no items; you carry " + carry);
            return 0.0;
        }
        int unused = 0;
        for (int i = 0; i < 2 * N; ++i)
            unused += locused[i] ? 0 : 1;
        if (unused > 0) {
            addFatalError("You must pick up all items and deliver them to all locations; you have " + unused + " unused locations.");
            return 0.0;
        }

        // now that the path is known to be valid, score it
        double score = 0;
        for (int i = 1; i < nP; ++i)
            score += scoreSeg(path[i - 1], path[i]);
        return score;
    }
    catch (Exception e) {
        System.err.println("An exception occurred while trying to get your program's results.");
        e.printStackTrace();
        return 0;
    }
    }
// ------------- visualization part ------------
    JFrame jf;
    Vis v;
    static String exec;
    static boolean vis, debug;
    static Process proc;
    InputStream is;
    OutputStream os;
    BufferedReader br;
    static int SZ, W, H;
    // -----------------------------------------
    double[] getPath(String[] map, double[] locations, int capacity) throws IOException {
        if (proc == null) {
            return new double[0];
        }
        StringBuffer sb = new StringBuffer();
        sb.append(map.length).append("\n");
        for (String st : map) {
            sb.append(st).append("\n");
        }
        sb.append(locations.length).append("\n");
        for (double d : locations) {
            sb.append(d).append("\n");
        }
        sb.append(capacity).append("\n");;
        os.write(sb.toString().getBytes());
        os.flush();

        int retN = Integer.parseInt(br.readLine());
        double[] ret = new double[retN];
        for (int i = 0; i < retN; ++i)
            ret[i] = Double.parseDouble(br.readLine());
        return ret;
    }
    // -----------------------------------------
    int[] cs = {0xe5cf9f, 0xe6b46a, 0xe7a036, 0xdb8c15, 0xbf7c25, 0x9d6625, 0x774d1b, 0x5b3107, 0x362207, 0x180e02};
    // -----------------------------------------
    int getCoord(double c) {
        return (int)Math.floor(c * SZ);
    }
    // -----------------------------------------
    public class Vis extends JPanel implements WindowListener {
        public void paint(Graphics g) {
            BufferedImage cache = new BufferedImage(W, H, BufferedImage.TYPE_INT_RGB);
            Graphics2D g2 = (Graphics2D)cache.getGraphics();
            // background
            g2.setColor(new Color(0xDDDDDD));
            g2.fillRect(0, 0, W, H);

            // current colors of the cells of the board (draw every cell)
            for (int i = 0; i < S; ++i)
            for (int j = 0; j < S; ++j) {
                g2.setColor(new Color(cs[map[i].charAt(j) - '0']));
                g2.fillRect(j * SZ, i * SZ, SZ, SZ);
            }

            // pick up/drop off locations
            // green for pick up, red for drop off
            g2.setColor(Color.GREEN);
            for (int i = 0; i < N; ++i) {
                int x = getCoord(loc[i].x);
                int y = getCoord(loc[i].y);
                g2.fillOval(x - 2, y - 2, 4, 4);
            }

            g2.setColor(Color.RED);
            for (int i = 0; i < N; ++i) {
                int x = getCoord(loc[i + N].x);
                int y = getCoord(loc[i + N].y);
                g2.fillOval(x - 2, y - 2, 4, 4);
            }

            // path
            if (path != null) {
                g2.setColor(Color.BLUE);
                for (int i = 0; i < path.length - 1; ++i) {
                    g2.drawLine(getCoord(path[i].x), getCoord(path[i].y), getCoord(path[i + 1].x), getCoord(path[i + 1].y));
                }
            }

            // palette of colors for terrain types on the side to illustrate
            /*for (int i = 0; i < T; ++i) {
                g2.setColor(new Color(cs[i]));
                g2.fillRect((S + 1) * SZ, i * SZ, 3 * SZ, SZ);
            }*/
            g.drawImage(cache,0,0,W,H,null);
        }
        // -------------------------------------
        public Vis() {
            jf.addWindowListener(this);
        }
        // -------------------------------------
        // WindowListener
        public void windowClosing(WindowEvent e){
            if (proc != null)
                try { proc.destroy(); }
                catch (Exception ex) { ex.printStackTrace(); }
            System.exit(0);
        }
        public void windowActivated(WindowEvent e) { }
        public void windowDeactivated(WindowEvent e) { }
        public void windowOpened(WindowEvent e) { }
        public void windowClosed(WindowEvent e) { }
        public void windowIconified(WindowEvent e) { }
        public void windowDeiconified(WindowEvent e) { }
    }

    // -----------------------------------------
    public TerrainCrossingVis(String seed) {
      try {
        if (vis) {
            jf = new JFrame();
            v = new Vis();
            jf.getContentPane().add(v);
        }
        if (exec != null) {
            try {
                Runtime rt = Runtime.getRuntime();
                proc = rt.exec(exec);
                os = proc.getOutputStream();
                is = proc.getInputStream();
                br = new BufferedReader(new InputStreamReader(is));
                new ErrorReader(proc.getErrorStream()).start();
            } catch (Exception e) { e.printStackTrace(); }
        }
        System.out.println("Score = " + runTest(seed));
        if (proc != null)
            try { proc.destroy(); }
            catch (Exception e) { e.printStackTrace(); }
      }
      catch (Exception e) { e.printStackTrace(); }
    }
    // -----------------------------------------
    public static void main(String[] args) {
        String seed = "1";
        vis = true;
        SZ = 20;
        for (int i = 0; i<args.length; i++) {
            if (args[i].equals("-seed"))
                seed = args[++i];
            if (args[i].equals("-exec"))
                exec = args[++i];
            if (args[i].equals("-novis"))
                vis = false;
            if (args[i].equals("-size"))
                SZ = Integer.parseInt(args[++i]);
            if (args[i].equals("-debug"))
                debug = true;
        }
        TerrainCrossingVis f = new TerrainCrossingVis(seed);
    }
    // -----------------------------------------
    void addFatalError(String message) {
        System.out.println(message);
    }
}

class ErrorReader extends Thread {
    InputStream error;
    public ErrorReader(InputStream is) {
        error = is;
    }
    public void run() {
        try {
            byte[] ch = new byte[50000];
            int read;
            while ((read = error.read(ch)) > 0)
            {   String s = new String(ch,0,read);
                System.out.print(s);
                System.out.flush();
            }
        } catch(Exception e) { }
    }
}
