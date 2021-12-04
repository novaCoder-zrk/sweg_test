package sweg.algorithm;

import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import jdk.jshell.spi.ExecutionControl;
import sweg.SummaryGraphModule;

import java.time.Instant;

import static java.lang.Integer.min;
import static java.lang.Long.min;

public class SWeG extends SummaryGraphModule {

    private int T = 20; // number of iterations
    private double eps = 0; // error bound

    private ObjectArrayList<IntOpenHashSet> input = new ObjectArrayList<>(0);
    private IntArrayList[] adjList;

    private IntArrayList[] S;
    private int[] deg;
    private Int2IntOpenHashMap[] Super2Node, Node2Super;

    private int[] h;
    private int[] shingles;
    private ObjectArrayList<IntArrayList> groups;
    private int threshold = 500, max_step = 10;
    private int iter_max_step = 0;
    private Instant _start = Instant.now();

    private long uniq = 0;
    private long[] inGroup;
    private int[] val = new int[threshold];

    private int[] candidates; private int sn;
    private long[][] sep; private int sepn;

    public SWeG(boolean directed){
        super(directed);
        if(directed){
            try {
                throw new ExecutionControl.NotImplementedException("Directed version is NOT_IMPLEMENTED");
            } catch (ExecutionControl.NotImplementedException e) {
                e.printStackTrace();
                System.exit(-1);
            }
        }
    }

    @Override
    public void addVertex(int idx) {
        super.addVertex(idx);
        input.add(new IntOpenHashSet(1));
    }

    @Override
    public void processEdge(final int src, final int dst, final boolean add) {
        if(add){
            input.get(src).add(dst);
            input.get(dst).add(src);
        }else{
            input.get(src).remove(dst);
            input.get(dst).remove(src);
        }
    }

    private void addSuperEdge(int u, int v){
        // add directed edge u -> v
        if(u == v && getSize(u) == 1) return;
        P.addEdge(u, v);
        for(int _u: S[u]){
            for(int _v: S[v]){
                if(_u == _v) continue;
                if(!adjList[_u].contains(_v)){
                    Cm.addEdge(_u, _v);
                }
            }
        }
    }

    private int getSize(final int Sv){
        return S[Sv].size();
    }

    private long getPi(final int Su, final int Sv){
        long pi = getSize(Su); pi *= getSize(Sv);
        if(Su == Sv){
            pi -= getSize(Sv);
            pi /= 2;
        }
        return pi;
    }

    private Int2IntOpenHashMap getSuperDegree(int Sv){
        Int2IntOpenHashMap SuperDeg = new Int2IntOpenHashMap();
        for(Int2IntMap.Entry v: Super2Node[Sv].int2IntEntrySet()){
            SuperDeg.addTo(V.getInt(v.getIntKey()), v.getIntValue());
        }
        return SuperDeg;
    }

    private long getCostInner(int v, int vSize, Int2IntOpenHashMap Nv, boolean add){
        long cost = 0;
        for(Int2IntMap.Entry nbr: Nv.int2IntEntrySet()){
            long pi, edgeCount;
            if(v == nbr.getIntKey()){
                pi = vSize; pi *= (vSize - 1);
                edgeCount = nbr.getIntValue();
                pi /= 2; edgeCount /= 2;
            }else{
                pi = vSize; pi *= getSize(nbr.getIntKey());
                edgeCount = nbr.getIntValue();
            }
            cost += min(pi - edgeCount + 1, edgeCount);

            if(add && (pi - edgeCount + 1) < edgeCount) {
                addSuperEdge(v, nbr.getIntKey());
            }

        }
        return cost;
    }


    private long getCost(int v){
        return getCostInner(v, getSize(v), getSuperDegree(v), false);
    }

    private long getMergeCost(int u, int v){
        Int2IntOpenHashMap Nw = getSuperDegree(u);
        for(Int2IntMap.Entry nbr: getSuperDegree(v).int2IntEntrySet()){
            Nw.addTo(nbr.getIntKey(), nbr.getIntValue());
        }
        Nw.addTo(v, Nw.getOrDefault(u, 0));
        Nw.remove(u);
        return getCostInner(v, getSize(u) + getSize(v), Nw, false);
    }

    private double saving(int A, int B){
        long before = getCost(A) + getCost(B);
        long after = getMergeCost(A, B);
        long pi = getSize(A); pi *= getSize(B);
        long edgeCount = 0;
        for(int v: S[A]){
            edgeCount += Node2Super[v].getOrDefault(B, 0);
        }
        before -= min(edgeCount, pi - edgeCount + 1);
        //System.out.println(A + " vs " + B + " : " + before + " -> " + after);
        return 1 - (after / (double)before);
    }

    private void divideInner(int step){
        int me = step % 2;
        if(step == max_step){
            for(int ii=0;ii<sepn;ii++){
                int s = (int)(sep[me][ii] >> 32), e = (int)(sep[me][ii] & 0x7FFFFFFFL);
                for(int i=s;i<=e;i+=threshold){
                    if(i + (threshold - 1) <= e){
                        groups.add(new IntArrayList(candidates, i, threshold));
                    }else{
                        groups.add(new IntArrayList(candidates, i, (e-i+1)));
                    }
                }
            }
            return;
        }

        h[0] = 1;
        for (int i = 1; i < n; i++) {
            h[i] = i+1;
            int randIdx = randInt(0, i);
            h[i] = h[randIdx];
            h[randIdx] = i+1;
        }

        int nxt_sepn = 0;
        for(int ii=0;ii<sepn;ii++){
            int s = (int)(sep[me][ii] >> 32), e = (int)(sep[me][ii] & 0x7FFFFFFFL);
            for(int i=s;i<=e;i++){
                final int A = candidates[i];
                int minHash = 0x7FFFFFFF;
                for(int v: S[A]){
                    minHash = min(minHash, h[v]);
                    for(int u: adjList[v]){
                        minHash = min(minHash, h[u]);
                    }
                }
                shingles[A] = minHash;
            }
            IntArrays.parallelQuickSort(candidates, s, e + 1, new IntComparator() {
                @Override
                public int compare(int i, int i1) {
                    return Integer.compare(shingles[i], shingles[i1]);
                }
            });
            int prv = s;
            for(int i=s+1;i<=e;i++){
                if(shingles[candidates[i]] != shingles[candidates[i-1]]){
                    if(i - prv <= threshold){
                        groups.add(new IntArrayList(candidates, prv, i-prv));
                    }else{
                        sep[1-me][nxt_sepn++] = (((long)prv) << 32) + (i-1);
                    }
                    prv = i;
                }
            }
            if((e + 1) - prv <= threshold){
                groups.add(new IntArrayList(candidates, prv, (e+1) - prv));
            }else{
                sep[1-me][nxt_sepn++] = (((long)prv) << 32) + e;
            }
        }
        if(nxt_sepn > 0){
            sepn = nxt_sepn;
            divideInner(step + 1);
        }
    }

    private void divide(){
        // Input: input graph G = (V, E), current supernodes S
        // Output: disjoint groups of supernodes (shingle)
        sn = 0;
        for(int i=0;i<n;i++){
            if(deg[i] > 0) candidates[sn++] = i;
        }
        sepn = 1;
        sep = new long[2][(sn + threshold + threshold) / threshold];
        sep[0][0] = sn-1;
        groups = new ObjectArrayList<>(n / threshold);
        divideInner(0);
        /*
        int _sn = 0;
        IntArrayList checky = new IntArrayList(n);
        for(int i=0;i<n;i++){
            checky.add(0);
        }
        for(IntArrayList _Q: groups){
            _sn += _Q.size();
            if(_Q.size() > threshold) System.out.println("@");
            for(int v: _Q){
                if(checky.getInt(v) > 0){ System.out.println("?"); }
                checky.set(v, 1);
            }
        }
        System.out.println(sn + " " + _sn);
         */
    }

    private void merge(int iter){
        for(IntArrayList _Q: groups){
            sn = 0;
            IntArrayList Q = new IntArrayList(_Q);
            for(int q: Q) {
                Super2Node[q] = new Int2IntOpenHashMap(0);
            }
            for(int q: Q){
                for(int v: S[q]){
                    if(Node2Super[v] == null){
                        candidates[sn++] = v;
                        Node2Super[v] = new Int2IntOpenHashMap(0);
                    }
                    for(int u: adjList[v]){
                        if(Node2Super[u] == null){
                            candidates[sn++] = u;
                            Node2Super[u] = new Int2IntOpenHashMap(0);
                        }
                        Super2Node[q].addTo(u, 1);
                        Node2Super[u].addTo(q, 1);
                        if(Super2Node[V.getInt(u)] == null) Node2Super[v].addTo(V.getInt(u), 1);
                    }
                }
            }
            int sz = Q.size();
            while(sz > 1){
                uniq += 1;
                // pick and remove random supernode A from Q
                int randIdx = randInt(0, sz-1);
                int A = Q.getInt(randIdx);
                Q.set(randIdx, Q.getInt(sz-1));
                Q.popInt(); sz -= 1;
                for(int i=0;i<sz;i++){
                    inGroup[Q.getInt(i)] = uniq * (threshold + 1) + (i + 1);
                    val[i] = 0;
                }
                for(Int2IntMap.Entry v: Super2Node[A].int2IntEntrySet()){
                    for(Int2IntMap.Entry U: Node2Super[v.getIntKey()].int2IntEntrySet()){
                        if(inGroup[U.getIntKey()] / (threshold + 1) == uniq){
                            // Compute sum of min(w(A, v), w(B, v))
                            val[(int)(inGroup[U.getIntKey()] % (threshold + 1)) - 1] += min(v.getIntValue(), U.getIntValue());
                        }
                    }
                }
                double maxSuperJaccard = -1.0;
                int argMax = -1;
                for(int i=0;i<sz;i++){
                    double superJaccard = val[i] / (double) (deg[A] + deg[Q.getInt(i)] - val[i]);
                    if(maxSuperJaccard < superJaccard){
                        maxSuperJaccard = superJaccard;
                        argMax = i;
                    }
                }
                int B = Q.getInt(argMax);
                Q.set(argMax, Q.getInt(sz-1));
                Q.set(sz-1, B);

                double threshold = (iter < T) ? (1.0 / (double) (1 + iter)) : 0;
                double result = saving(A, B);
                if(result >= threshold){
                    if(getSize(A) + Super2Node[A].size() < getSize(B) + Super2Node[B].size()){
                        int tmp = A; A = B; B = tmp;
                    }
                    // merge A and B: A <- B
                    for(int v: S[B]){
                        S[A].add(v);
                        V.set(v, A);
                    }
                    deg[A] += deg[B];
                    deg[B] = 0;
                    for(Int2IntMap.Entry v: Super2Node[B].int2IntEntrySet()){
                        Node2Super[v.getIntKey()].addTo(A, v.getIntValue());
                        Node2Super[v.getIntKey()].remove(B);
                        Super2Node[A].addTo(v.getIntKey(), v.getIntValue());
                    }

                    S[B] = null;
                    Super2Node[B] = null;
                    Q.set(sz-1, A);
                }
            }
            // Clear Super2Node and Node2Super
            for(int q: _Q){
                Super2Node[q] = null;
            }
            for(int i=0;i<sn;i++){
                Node2Super[candidates[i]] = null;
            }
        }
        groups.clear();
    }

    private void encode(){
        for(int i=0;i<n;i++){
            for(int j: adjList[i]){
                if(Super2Node[V.getInt(i)] == null) Super2Node[V.getInt(i)] = new Int2IntOpenHashMap(0);
                Super2Node[V.getInt(i)].addTo(V.getInt(j), 1);
            }
        }
        // add superedges (P and Cm)
        for(int i=0;i<n;i++){
            if(deg[i] > 0) getCostInner(i, getSize(i), Super2Node[i], true);
        }
        // add correction edges (Cp)
        for(int i=0;i<n;i++){
            int Si = V.getInt(i);
            for(int v: adjList[i]){
                int Sv = V.getInt(v);
                if(!P.getNeighbors(Si).contains(Sv)){
                    Cp.addEdge(i, v);
                }
            }
        }
    }

    @Override
    public void processBatch(){
        System.out.println("|V|: " + n);

        adjList = new IntArrayList[n];
        for(int i=0;i<n;i++){
            adjList[i] = new IntArrayList(input.get(i));
            input.set(i, null);
        }
        input = null;

        // Overview of SWeG
        // Input: input graph G = (V, E), iterations T, error bound eps
        // Output: Summary graph G' = (S, P), corrections C+, C-

        S = new IntArrayList[n];

        // for divide
        h = new int[n];
        shingles = new int[n];
        candidates = new int[n+1];

        Super2Node = new Int2IntOpenHashMap[n];
        Node2Super = new Int2IntOpenHashMap[n];
        inGroup = new long[n];
        deg = new int[n];

        // initialize supernodes S to {{v}: v \in V}
        for(int i=0;i<n;i++){
            S[i] = new IntArrayList(1);
            S[i].add(i);
            Super2Node[i] = null;
            Node2Super[i] = null;
            inGroup[i] = 0;
            deg[i] = adjList[i].size();
        }

        for(int iter=1;iter<=T;iter++){
            // divide S into disjoint groups
            divide();
            // merge some supernodes within each group
            merge(iter);
        }
        // encode edges E into superedges P and corrections C
        encode();
        // We only implemented lossless version of SWeG,
        // since our proposed algorithm MoSSo only supports lossless summarization.
        // System.out.println(uniq * (threshold + 1) + (threshold + 1));
    }
}
