#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <tuple>
#include <numeric>
#include <unordered_set>
#include <set>


using namespace std;

enum class STATEGY {
    AGGERAGATE,
    DISPERSE,
    RANDOM
};

struct Stat {

    tuple<int, int, int> storageoverhead;
    vector<int> nodes;
    int out_trx;
    int in_trx;

    Stat(int p_nodesnum = 0) :
            nodes(p_nodesnum, 0),
            out_trx{0},
            in_trx{0} {
        std::iota(nodes.begin(), nodes.end(), 0);
    }

};

static vector<Stat> clusterstat;


static void layoutprint(const vector<tuple<vector<int>, vector<int>, vector<int>>> &layout,
                        int c,
                        int start,
                        int printnum) {
    cout << "   " << "data blocks location by cluster index|" << "local parity location by cluster index|"
         << "global parity location by cluster index\n";
    for (int i=0;i<printnum;++i) {
        auto[dataloc, lploc, gploc] = layout[start+i];
        cout << "s" << to_string(start+i) << ":";
        cout << "{ ";
        for (auto ci:dataloc) {
            cout << to_string(ci) << " ";
        }
        cout << "}|{ ";
        for (auto ci:lploc) {
            cout << to_string(ci) << " ";
        }
        cout << "}|{ ";
        for (auto ci:gploc) {
            cout << to_string(ci) << " ";
        }
        cout << "}\n";
    }
}

auto singlestriperesolve(
        const tuple<int, int, int> &para
) {
    vector<tuple<int, int, int>> ret;
    auto[k, l, g] = para;
    int r = k / l;

    //cases
    //case1:
    int theta1 = 0;
    int theta2 = 0;

    if (r <= g) {
        theta1 = g / r;
        int i = 0;
        for (i = 0; i + theta1 <= l; i += theta1) ret.emplace_back(theta1 * r, theta1, 0);
        if (i < l) {
            ret.emplace_back((l - i) * r, (l - i), 0);
        }
        ret.emplace_back(0, 0, g);
        ret.emplace_back(-1, theta1, l - i);
    } else if (0 == (r % (g + 1))) {
        theta2 = r / (g + 1);
        for (int i = 0; i < l; ++i) {
            for (int j = 0; j < theta2; ++j) ret.emplace_back(g + 1, 0, 0);
        }
        ret.emplace_back(0, l, g);
        ret.emplace_back(-1, -1, 0);
    } else {
        int m = r % (g + 1);
        theta2 = r / (g + 1);
        //each local group remains m data blocks and 1 lp block
        // m[<=>r] < g -> case1
        theta1 = g / m;
        for (int i = 0; i < l; ++i) {
            for (int j = 0; j < theta2; ++j) {
                ret.emplace_back(g + 1, 0, 0);
            }
        }

        ret.emplace_back(-1, -1, -1);//acts as a marker

        int i = 0;
        for (i = 0; i + theta1 <= l; i += theta1) ret.emplace_back(theta1 * m, theta1, 0);
        if (i < l) {
            ret.emplace_back((l - i) * m, (l - i), 0);
        }
        ret.emplace_back(0, 0, g);

        ret.emplace_back(-1, theta1, l - i);
    }
    return ret;
}

//TODO : load balanced sampler
static auto random_sample(vector<int> u,
                          int pick) {
    random_shuffle(u.begin(), u.end());
    return vector<int>(u.begin(), u.begin() + pick);
}

static void getcandidatecluster(vector<int> &ret,
                                unordered_set<int> whitelist,
                                unordered_set<int> blacklist,
                                int c) {

}

static auto generatelayout(const tuple<int, int, int> &para,
                           STATEGY s,
                           int c,
                           int stripenum,
                           int step = 2
) {
    vector<int> totalcluster(c, 0);
    iota(totalcluster.begin(), totalcluster.end(), 0);
    vector<tuple<vector<int>, vector<int>, vector<int>>> stripeslayout;
    auto[k, l, g] = para;
    int r = k / l;

    auto layout = singlestriperesolve(para);
    if (s == STATEGY::DISPERSE) {
        for (int j = 0; j * step < stripenum; ++j) {
            vector<vector<int>> datablock_location(step, vector<int>(k, -1));
            vector<vector<int>> lpblock_location(step, vector<int>(l, -1));
            vector<vector<int>> gpblock_location(step, vector<int>(g, -1));
            if (r <= g) {
                // case1
                // s1: D0D1L0 D2D3L1 G0G1G2
                // s2: D0D1L0 D2D3L1 G0G1G2



                vector<int> clusters(totalcluster.begin(), totalcluster.end());
                random_shuffle(clusters.begin(), clusters.end());
                int idx_l = 0;
                int idx_d = 0;
                int n = 0;
                //to pack residue
                auto[_ignore1, theta1, res_grpnum] = layout.back();

                int lim = (0 == res_grpnum ? layout.size() - 2 : layout.size() - 3);
                for (int i = 0; i < lim; ++i) {
                    auto[cluster_k, cluster_l, cluster_g]=layout[i];
                    for (int u = 0; u < step; ++u) {
                        int idx_l1 = idx_l;
                        int idx_d1 = idx_d;
                        for (int x = 0; x < cluster_l; ++x) {
                            lpblock_location[u][idx_l1] = clusters[n];
                            idx_l1++;
                            for (int m = 0; m < r; ++m) {
                                datablock_location[u][idx_d1] = clusters[n];
                                idx_d1++;
                            }
                        }
                        n++;
                    }
                    idx_d += cluster_k;
                    idx_l += cluster_l;
                }


                if (res_grpnum) {
                    int cur_res_grp = 0;
                    for (int x = 0; x < step; ++x) {
                        if (cur_res_grp + res_grpnum > theta1) {
                            n++;//next cluster
                            cur_res_grp = 0;//zero
                        }

                        int cur_res_idxl = idx_l;
                        int cur_res_idxd = idx_d;

                        for (int y = 0; y < res_grpnum; ++y) {
                            lpblock_location[x][cur_res_idxl++] = clusters[n];
                            for (int i = 0; i < r; ++i) {
                                datablock_location[x][cur_res_idxd++] = clusters[n];
                            }
                        }

                        cur_res_grp += res_grpnum;

                        //put res_grpnum residue group into cluster
                    }

                }


                //global cluster
                for (int x = 0; x < step; ++x) {
                    for (int i = 0; i < g; ++i) gpblock_location[x][i] = clusters[n];
                    n++;
                }


            } else if (0 == r % (g + 1)) {

                //case2
                //D0D1D2   D3D4D5   G0G1L0L1
                vector<int> clusters(totalcluster.begin(), totalcluster.end());
                random_shuffle(clusters.begin(), clusters.end());
                int n = 0;
                for (int i = 0; i < step; ++i) {
                    int idxd = 0;
                    for (int x = 0; x < layout.size() - 2; ++x) {
                        for (int y = 0; y < g + 1; ++y) {
                            datablock_location[i][idxd++] = clusters[n];
                        }
                        n++;
                    }


                    //g and l parity cluster
                    for (int x = 0; x < g; ++x) {
                        gpblock_location[i][x] = clusters[n];
                    }
                    for (int x = 0; x < l; ++x) {
                        lpblock_location[i][x] = clusters[n];
                    }
                    n++;
                }


            } else {
                //special case3
                vector<int> clusters(totalcluster.begin(), totalcluster.end());
                random_shuffle(clusters.begin(), clusters.end());
                int residue = find(layout.cbegin(), layout.cend(), tuple<int, int, int>{-1, -1, -1}) - layout.cbegin();
                int round = (r / (g + 1)) * (g + 1);
                // [round,r)*cursor+lp
                auto[_ignore1, pack_cluster_cap, frac_cluster_num] = layout.back();
                int s = 0;


                int packed_cluster_num = 0;
                int packed_residue = 0;
                if (frac_cluster_num) {
                    packed_cluster_num = (pack_cluster_cap / frac_cluster_num);
                    packed_residue = step % packed_cluster_num;
                    s = (0 != packed_residue) ? step / packed_cluster_num + 1 : step / packed_cluster_num;
                }//require s clusters to pack the residue groups




                int n = s;
                int m = residue + 1;
                int x = 0;
                int y = 0;
                int cursor = 1;
                int lim = (0 == frac_cluster_num ? layout.size() - 2 : layout.size() - 3);

                //TODO : cache optimization nested for loop
                for (; m < lim; ++m) {
                    auto[residuecluster_k, residuecluster_l, residuecluster_g] = layout[m];
                    int residuecluster_r = residuecluster_k / residuecluster_l;
                    for (int u = 0; u < step; ++u) {
                        int cur_cursor = cursor;
                        for (x = 0; x < residuecluster_l; ++x) {
                            for (y = 0; y < residuecluster_r; ++y) {
                                datablock_location[u][(cur_cursor - 1) * r + y + round] = clusters[n];
                            }
                            lpblock_location[u][cur_cursor - 1] = clusters[n];
                            cur_cursor++;
                        }
                        n++;
                    }
                    cursor += residuecluster_l;

                }
                if (frac_cluster_num) {
                    auto[res_k, res_l, res_g] = layout[m];
                    int res_r = res_k / res_l;
                    int cur_back = cursor;
                    for (int u1 = 0; u1 < step; ++u1) {
                        cursor = cur_back;
                        if (0 != u1 && 0 == (u1 % packed_cluster_num)) {
                            n++;
                        }
                        for (int x1 = 0; x1 < res_l; ++x1) {
                            lpblock_location[u1][cursor - 1] = clusters[n];
                            for (int y1 = 0; y1 < res_r; ++y1) {
                                datablock_location[u1][(cursor - 1) * r + round + y1] = clusters[n];
                            }
                            cursor++;
                        }
                    }
                }

                n++;
                //pack remained

                for (int u = 0; u < step; ++u) {
                    for (x = 0; x < l; ++x) {
                        for (y = 0; y < round; ++y) {
                            datablock_location[u][x * r + y] = clusters[n];
                            if (0 == ((y + 1) % (g + 1))) n++;
                        }
                    }
                }
                for (int u = 0; u < step; ++u) {
                    for (x = 0; x < g; ++x) {
                        gpblock_location[u][x] = clusters[n];
                    }
                    n++;
                }

            }
            for (int i = 0; i < step; ++i) {
                stripeslayout.emplace_back(datablock_location[i], lpblock_location[i], gpblock_location[i]);
            }

            datablock_location.assign(step, vector<int>(k, -1));
            lpblock_location.assign(step, vector<int>(l, -1));
            gpblock_location.assign(step, vector<int>(g, -1));
        }
    } else if (s == STATEGY::AGGERAGATE) {
        for (int j = 0; j * step < stripenum; ++j) {
            vector<int> clusters(totalcluster.begin(), totalcluster.end());
            random_shuffle(clusters.begin(), clusters.end());
            vector<int> datablock_location(k, -1);
            vector<int> lpblock_location(l, -1);
            vector<int> gpblock_location(g, -1);
            //ignore tuple if any element <0
            int idxd = 0;
            int idxl = 0;
            int n = 0;
            int r = k / l;
            int idxr = 0;
            int res_idxl = 0;
            int round = (r / (g + 1)) * (g + 1);
            for (auto cluster:layout) {
                if (auto[_k, _l, _g] = cluster;_k < 0 || _l < 0 || _g < 0) {
                    continue;
                } else {
                    //put _k data blocks ,_l lp blocks and _g gp blocks into this cluster

                    if (r <= g) {
                        for (int i = 0; i < _l; ++i) {
                            lpblock_location[idxl++] = clusters[n];
                        }
                        for (int i = 0; i < _k; ++i) {
                            datablock_location[idxd++] = clusters[n];
                        }
                        for (int i = 0; i < _g; ++i) {
                            gpblock_location[i] = clusters[n];
                        }

                        n++;
                    } else if (0 == (r % (g + 1))) {
                        //inner group index idxr
                        //group index idxl
                        //each cluster contain idxl group's subgroup
                        //D0D1   D2D3    G0L0
                        //D4D5   D6D7    G1L1

                        for (int i = 0; i < _k; ++i) {
                            datablock_location[idxl * r + idxr] = clusters[n];
                            idxr++;
                        }
                        if (idxr == r) {
                            idxl++;//+= theta1
                            idxr = 0;
                        }
                        for (int i = 0; i < _l; ++i) {
                            lpblock_location[i] = clusters[n];
                        }
                        for (int i = 0; i < _g; ++i) {
                            gpblock_location[i] = clusters[n];
                        }
                        n++;
                    } else {
                        //special case3
                        //D0D1D2   D3D4D5   D6L0   G1G2   D7D8D9   D10D11D12   D14D15D16   D17L2
                        //                  D13L1
                        if (0 == _l && 0 == _g) {
                            for (int i = 0; i < _k; ++i) {
                                datablock_location[idxl * r + idxr] = clusters[n];
                                idxr++;

                            }
                            if (idxr == round) {
                                idxl++;
                                idxr = 0;
                            }
                        } else if (0 == _g) {
                            idxr = round;
                            int res_idxr = idxr;
                            for (int i = 0; i < _l; ++i) {
                                for (int x = 0; x < (_k / _l); ++x) {
                                    datablock_location[res_idxl * r + res_idxr] = clusters[n];
                                    res_idxr++;
                                    if (res_idxr == r) {
                                        lpblock_location[res_idxl] = clusters[n];
                                        res_idxl++;
                                        res_idxr = round;
                                    }
                                }
                            }

                        } else {
                            //gp cluster
                            for (int i = 0; i < _g; ++i) {
                                gpblock_location[i] = clusters[n];
                            }
                        }
                        n++;

                    }
                }
            }

            for (int i = 0; i < step; ++i) {
                stripeslayout.emplace_back(datablock_location, lpblock_location, gpblock_location);
            }
        }
    } else {

        //random
        vector<int> clusters(totalcluster.begin(), totalcluster.end());
        for (int j = 0; j < stripenum; ++j) {
            random_shuffle(clusters.begin(), clusters.end());
            vector<int> datablock_location(k, -1);
            vector<int> lpblock_location(l, -1);
            vector<int> gpblock_location(g, -1);
            int n = 0;
            int idxl = 0;
            int idxd = 0;
            int idxr = 0;
            if (r <= g) {
                for (auto cluster:layout) {
                    auto[_k, _l, _g]=cluster;
                    if (_k < 0) continue;
                    for (int x = 0; x < _l; ++x) {
                        int _r = _k / _l;//_r == r
                        for (int y = 0; y < _r; ++y) {
                            datablock_location[idxd] = clusters[n];
                            idxd++;
                        }
                        lpblock_location[idxl] = clusters[n];
                        idxl++;
                    }

                    for (int x = 0; x < _g; ++x) {
                        gpblock_location[x] = clusters[n];
                    }
                    n++;
                }

            } else if (0 == (r % (g + 1))) {

                for (auto cluster:layout) {
                    auto[_k, _l, _g]= cluster;
                    if (_k < 0 || _l < 0 || _g < 0) continue;
                    for (int i = 0; i < _k; ++i) {
                        datablock_location[idxd++] = clusters[n];
                    }
                    for (int i = 0; i < _l; ++i) {

                        lpblock_location[idxl++] = clusters[n];
                    }
                    for (int i = 0; i < _g; ++i) {
                        gpblock_location[i] = clusters[n];
                    }
                    n++;
                }
            } else {
                int round = ((k / l) / (g + 1)) * (g + 1);

                for (auto cluster:layout) {
                    auto[_k, _l, _g]= cluster;
                    if (_k < 0 || _l < 0 || _g < 0) continue;
                    if (_l == 0 && _g == 0) {
                        //normal cluster
                        for (int x = 0; x < g + 1; ++x) {
                            datablock_location[idxd++] = clusters[n];
                        }
                        if (0 == (idxd % round)) {
                            idxd = (idxd / r + 1) * r;
                        }
                        n++;
                    } else if (_g == 0) {
                        //residue cluster
                        for (int x = 0; x < _l; ++x) {
                            lpblock_location[idxl] = clusters[n];
                            for (int y = 0; y < (_k / _l); ++y) {
                                datablock_location[idxl * r + round + y] = clusters[n];
                            }
                            idxl++;
                        }
                        n++;
                    } else {
                        //gp cluster
                        for (int x = 0; x < _g; ++x) {
                            gpblock_location[x] = clusters[n];
                        }
                    }

                }

            }
            stripeslayout.emplace_back(datablock_location, lpblock_location, gpblock_location);

        }

    }
    return stripeslayout;
}



static pair<int, int> calculate_crosscluster(const vector<tuple<vector<int>, vector<int>, vector<int>>> &layout,
                                             const tuple<int, int, int> &para_before,
                                             const tuple<int, int, int> &para_after,
                                             int from, int step = 2,
                                             bool partial_gp = true,
                                             bool partial_lp=  true ){


    //maintain g
    int rec = 0;
    int mig = 0;
    auto[k, l, g] = para_before;
    auto[_k, _l, _g]= para_after;
    int r = k / l;
    //to calculate how many blocks should be transfered to reconstruct gp
    unordered_set<int> union_cluster;
    unordered_set<int> residue_cluster;
    unordered_set<int> overlap_cluster;
    //to calculate how many blocks should be migrated
    int cluster_cap = 0;
    int residuecluster_cap = 0;
    unordered_map<int, int> cluster_tmp_datablock;//differentiate partial decoding or naive decoding
    unordered_map<int, int> cluster_tmp;
    unordered_map<int, int> residuecluster_tmp;
    unordered_map<int, vector<pair<int, int>>> migration_plan_asnormal;//{stripeid,blockid}
    unordered_map<int, vector<pair<int, int>>> migration_plan_asresidue;//{stripeid,blockid}
    for (int i = 0; i < step; ++i) {
        auto[dataloc, lploc, gploc] = layout[from + i];
        for (auto ci:dataloc) {
            union_cluster.insert(ci);
            migration_plan_asnormal.insert({ci, vector<pair<int, int>>{}});
        }
        for (auto ci:lploc) {
            union_cluster.insert(ci);
            residue_cluster.insert(ci);
            migration_plan_asresidue.insert({ci, vector<pair<int, int>>{}});
        }
    }

    if (r <= g) {
        partial_lp=false;
        cluster_cap = _g;
        //migration
        for (int i = 0; i < step; ++i) {
            auto[dataloc, lploc, gploc] = layout[from + i];
            for (int j = 0; j < dataloc.size(); ++j) {
                int ci = dataloc[j];
                if (0 != cluster_tmp.count(ci) && cluster_tmp[ci] < cluster_cap) {
                    cluster_tmp[ci]++;
                    cluster_tmp_datablock[ci]++;
                } else if (0 != cluster_tmp.count(ci)) {
                    cluster_tmp_datablock[ci]++;
                    migration_plan_asnormal[ci].emplace_back(from + i, j);
                } else {
                    //unseen cluster
                    cluster_tmp_datablock.insert({ci, 1});
                    cluster_tmp.insert({ci, 1});
                }
            }
        }
        //reconsx
        for (auto c:cluster_tmp_datablock) {
            if(partial_gp) rec += min(c.second, _g);
            else rec+=c.second;
        }

        for (auto ci:union_cluster) {
            mig += migration_plan_asnormal[ci].size();
        }
        if(!partial_lp) mig+=(mig/r);
        return {rec, mig};
    } else if (0 == (r % (g + 1))) {


        cluster_cap = _g + (_g / g);
        //migration
        for (int i = 0; i < step; ++i) {
            auto[dataloc, lploc, gploc] = layout[from + i];
            for (int j = 0; j < dataloc.size(); ++j) {
                int ci = dataloc[j];
                if (0 != cluster_tmp.count(ci) && cluster_tmp[ci] < cluster_cap) {
                    cluster_tmp[ci]++;
                    cluster_tmp_datablock[ci]++;
                } else if (0 != cluster_tmp.count(ci)) {
                    // data blocks >= cluster_cap
                    cluster_tmp_datablock[ci]++;
                    migration_plan_asnormal[ci].emplace_back(from + i, j);
                } else {
                    //unseen cluster
                    cluster_tmp.insert({ci, 1});
                    cluster_tmp_datablock.insert({ci, 1});
                }
            }
        }
        //reconsx
        for (auto c:cluster_tmp_datablock) {
            if(partial_gp) rec += min(c.second, _g);
            else rec+=c.second;
        }

        for (auto ci:union_cluster) {
            mig += migration_plan_asnormal[ci].size();
        }

        return {rec, mig};
    } else {
        cluster_cap = _g + (_g / g);
        int m = (k / l) % (g + 1);//residue_datablock per group
        residuecluster_cap = (_g / m) * (m + 1);//count lp
        int mix_cap = _g+1+(g/m);
        //pick up residue cluster id
        unordered_map<int, vector<int>> cluster_precount;
        for (auto ci:union_cluster) {
            cluster_precount.insert({ci, vector<int>(step, 0)});
        }
        for (int i = 0; i < step; ++i) {
            auto[dataloc, lploc, gploc] = layout[from + i];
            for (auto ci:dataloc) {
                cluster_precount[ci][i]++;
            }
            for (auto ci:lploc) {
                cluster_precount[ci][i]++;
            }
        }

        unordered_map<int,unordered_set<int>> special_cluster_mover;
        unordered_set<int> mixpack_cluster;
        for (auto stat:cluster_precount) {
            auto[c, countlist] = stat;
            set<int> dedup(countlist.cbegin(), countlist.cend());
            if(dedup.size()>1)
            {
                //mixed
                mixpack_cluster.insert(c);
            }
            int cap = *max_element(countlist.cbegin(),countlist.cend());
            if(accumulate(countlist.cbegin(),countlist.cend(),0)<=(mix_cap))
            {
                continue;
            }
            if(residue_cluster.contains(c)) {
                if (cap > g + 1) {
                    if (!dedup.contains(g + 1))// ==2 but all residue
                    {
                        continue;
                    } else {
                        // ==2 one normal one residue , move normal
                        if (!special_cluster_mover.contains(c)) {
                            special_cluster_mover.insert(
                                    {c, unordered_set<int>{}});
                        }
                        for (int i = 0; i < step; ++i) {
                            if (countlist[i] == g + 1) {
                                special_cluster_mover[c].insert(i);
                            }
                        }
                    }
                } else if (cap == (g + 1)) {
                    if (!special_cluster_mover.contains(c)) {
                        special_cluster_mover.insert(
                                {c, unordered_set<int>{}});
                    }
                    for (int i = 0; i < step; ++i) {
                        if (countlist[i] != g + 1) {
                            special_cluster_mover[c].insert(i);
                        }
                    }
                } else {
                    // can be concat
                }
            }

        }

        vector<pair<int,int>> migration_plan;
        for (int i = 0; i < step; ++i) {
            auto [dataloc, lploc, gploc] = layout[from + i];
            for (int j = 0; j < dataloc.size(); ++j) {
                int ci = dataloc[j];
                if (0 != cluster_tmp.count(ci)) {
                    cluster_tmp_datablock[ci]++;
                    if(special_cluster_mover.contains(ci)&&special_cluster_mover[ci].contains(i))
                    {
                        migration_plan.emplace_back(from+i,j);
                        continue;
                    }
                    int cap = (0 == residue_cluster.count(ci)) ? cluster_cap : residuecluster_cap;
                    if(mixpack_cluster.contains(ci)) cap=mix_cap;
                    if (cluster_tmp[ci] < cap) {
                        cluster_tmp[ci]++;
                    } else {
                        migration_plan.emplace_back(from + i, j);
                    }
                } else {
                    cluster_tmp.insert({ci, 1});
                    cluster_tmp_datablock.insert({ci, 1});
                    if(special_cluster_mover.contains(ci)&&special_cluster_mover[ci].contains(i))
                    {
                        migration_plan.emplace_back(from+i,j);
                        continue;
                    }
                }
            }

            for (int j = 0; j < lploc.size(); ++j) {
                int ci = lploc[j];
                if(special_cluster_mover.contains(ci)&&special_cluster_mover[ci].contains(i))
                {
                    migration_plan.emplace_back(from+i,j+k);
                    continue;
                }
                int cap = (0 == residue_cluster.count(ci)) ? cluster_cap : residuecluster_cap;
                if(mixpack_cluster.contains(ci)) cap=mix_cap;
                if (cluster_tmp[ci] < cap) {
                    cluster_tmp[ci]++;
                } else {
                    migration_plan.emplace_back(from + i, j + k);
                }
            }

        }
        //reconsx
        for (auto c:cluster_tmp_datablock) {
            if(partial_gp) rec += min(c.second, _g);
            else rec+=c.second;
        }
        mig = migration_plan.size();
        return {rec, mig};
    }
}


int main()
{
    /* sim1
    {
    std::vector<std::tuple<int,int,int,int,int>> confs{
            {4,2,2,6,100},
            {6,2,2,6,100},
            {8,2,2,8,100},
            {10,2,2,10,100},
            {12,2,2,10,100},
            {12,3,3,8,100},
            {16,2,2,14,100},
            {16,2,3,10,100},
            {20,2,2,16,100},
            {20,2,4,10,100},
            {20,4,2,18,100},
            {20,4,4,10,100}
    };
    for(const auto & conf:confs) {
        const auto[k, l, g, c, stripenum] = conf;
        std::cout << "g same :k-l-g-c-s:" << k<<"-"<<l<<"-"<<g<<"-"<<c<<"-"<<stripenum<<"\n";
        tuple<int, int, int> paraBefore{k, l, g};
        tuple<int, int, int> paraAfter{2*k, 2*l, g};
        cout << "random with partial decoding :\n";
        auto res1 = generatelayout(paraBefore, STATEGY::RANDOM, c, stripenum, 2);
        {
            int totalrec = 0;
            int totalmig = 0;
            for (int i = 0; i < stripenum / 2; ++i) {
                auto[rec, mig] = calculate_crosscluster(res1, paraBefore, paraAfter, i * 2, 2);
//            if(mig!=0){
//                layoutprint(res1,c,2*i,2);
//                cout << 2 * i << " with " << 2 * i + 1 << " cost: " << rec << "-" << mig << endl;
//            }

                totalrec += rec;
                totalmig += mig;

            }
            cout << totalrec << "-" << totalmig << endl;
        }

        cout << "random without partial decoding :\n";
        {
            int totalrec = 0;
            int totalmig = 0;
            for (int i = 0; i < stripenum / 2; ++i) {
                auto[rec, mig] = calculate_crosscluster(res1, paraBefore, paraAfter, i * 2, 2, false, false);
//            if(mig!=0){
//                layoutprint(res1,c,2*i,2);
//                cout << 2 * i << " with " << 2 * i + 1 << " cost: " << rec << "-" << mig << endl;
//            }

                totalrec += rec;
                totalmig += mig;

            }
            cout << totalrec << "-" << totalmig << endl;
        }

        cout << "aggeragate :\n";
        {
            auto res = generatelayout(paraBefore, STATEGY::AGGERAGATE, c, stripenum, 2);
//        layoutprint(res, c,0,10);
            int totalrec = 0;
            int totalmig = 0;
            for (int i = 0; i < stripenum / 2; ++i) {
                auto[rec, mig] = calculate_crosscluster(res, paraBefore, paraAfter, i * 2, 2);
//            if(i<5) cout << 2 * i << " with " << 2 * i + 1 << " cost: " << rec << "-" << mig << endl;
                totalrec += rec;
                totalmig += mig;
            }
            cout << totalrec << "-" << totalmig << endl;
        }
        cout << "disperse :\n";
        {
            auto peek = generatelayout(paraBefore, STATEGY::DISPERSE, 20, 2, 2);
//        layoutprint(peek, 20,0,2);
            auto res = generatelayout(paraBefore, STATEGY::DISPERSE, c, stripenum, 2);
//        layoutprint(res, c,0,10);
            int totalrec = 0;
            int totalmig = 0;
            for (int i = 0; i < stripenum / 2; ++i) {
                auto[rec, mig] = calculate_crosscluster(res, paraBefore, paraAfter, i * 2, 2);
//            if(i<5) cout << 2 * i << " with " << 2 * i + 1 << " cost: " << rec << "-" << mig << endl;
                totalrec += rec;
                totalmig += mig;
            }
            cout << totalrec << "-" << totalmig << endl;
        }
    }
    for(const auto & conf:confs) {
        const auto[k, l, g, c, stripenum] = conf;
        std::cout << "g double :k-l-g-c-s:" << k<<"-"<<l<<"-"<<g<<"-"<<c<<"-"<<stripenum<<"\n";
        tuple<int, int, int> paraBefore{k, l, g};
        tuple<int, int, int> paraAfter{2*k, 2*l, 2*g};
        cout << "random with partial decoding :\n";
        auto res1 = generatelayout(paraBefore, STATEGY::RANDOM, c, stripenum, 2);
        {
            int totalrec = 0;
            int totalmig = 0;
            for (int i = 0; i < stripenum / 2; ++i) {
                auto[rec, mig] = calculate_crosscluster(res1, paraBefore, paraAfter, i * 2, 2);
//            if(mig!=0){
//                layoutprint(res1,c,2*i,2);
//                cout << 2 * i << " with " << 2 * i + 1 << " cost: " << rec << "-" << mig << endl;
//            }

                totalrec += rec;
                totalmig += mig;

            }
            cout << totalrec << "-" << totalmig << endl;
        }

        cout << "random without partial decoding :\n";
        {
            int totalrec = 0;
            int totalmig = 0;
            for (int i = 0; i < stripenum / 2; ++i) {
                auto[rec, mig] = calculate_crosscluster(res1, paraBefore, paraAfter, i * 2, 2, false, false);
//            if(mig!=0){
//                layoutprint(res1,c,2*i,2);
//                cout << 2 * i << " with " << 2 * i + 1 << " cost: " << rec << "-" << mig << endl;
//            }

                totalrec += rec;
                totalmig += mig;

            }
            cout << totalrec << "-" << totalmig << endl;
        }

        cout << "aggeragate :\n";
        {
            auto res = generatelayout(paraBefore, STATEGY::AGGERAGATE, c, stripenum, 2);
//        layoutprint(res, c,0,10);
            int totalrec = 0;
            int totalmig = 0;
            for (int i = 0; i < stripenum / 2; ++i) {
                auto[rec, mig] = calculate_crosscluster(res, paraBefore, paraAfter, i * 2, 2);
//            if(i<5) cout << 2 * i << " with " << 2 * i + 1 << " cost: " << rec << "-" << mig << endl;
                totalrec += rec;
                totalmig += mig;
            }
            cout << totalrec << "-" << totalmig << endl;
        }
        cout << "disperse :\n";
        {
            auto peek = generatelayout(paraBefore, STATEGY::DISPERSE, 20, 2, 2);
//        layoutprint(peek, 20,0,2);
            auto res = generatelayout(paraBefore, STATEGY::DISPERSE, c, stripenum, 2);
//        layoutprint(res, c,0,10);
            int totalrec = 0;
            int totalmig = 0;
            for (int i = 0; i < stripenum / 2; ++i) {
                auto[rec, mig] = calculate_crosscluster(res, paraBefore, paraAfter, i * 2, 2);
//            if(i<5) cout << 2 * i << " with " << 2 * i + 1 << " cost: " << rec << "-" << mig << endl;
                totalrec += rec;
                totalmig += mig;
            }
            cout << totalrec << "-" << totalmig << endl;
        }
    }

    }
    */



    {
        std::vector<std::tuple<int, int, int, int, int, int>> confs{
                {6, 2, 2, 6,  100, 2},
                {6, 2, 2, 9,  100, 3},
                {6, 2, 2, 12, 100, 4},
                {6, 2, 2, 15, 100, 5}
        };
        for (const auto &conf:confs) {
            cout <<"******************************************************************************************\n";
            //case g -> g
            const auto[k, l, g, c, stripenum, x] = conf;
            std::cout << "g same :k-l-g-c-s-x:" << k << "-" << l << "-" << g << "-" << c << "-" << stripenum << "-" << x
                      << "\n";
            tuple<int, int, int> paraBefore{k, l, g};
            tuple<int, int, int> paraAfter{x * k, x * l, g};
            cout << "random with partial decoding :\n";
            auto res1 = generatelayout(paraBefore, STATEGY::RANDOM, c, stripenum, x);
            {
                int totalrec = 0;
                int totalmig = 0;
                for (int i = 0; x * i < stripenum; ++i) {
                    auto[rec, mig] = calculate_crosscluster(res1, paraBefore, paraAfter, i * x, x);
                    totalrec += rec;
                    totalmig += mig;
                }
                cout << totalrec << "-" << totalmig << endl;
            }

            cout << "random without partial decoding :\n";
            {
                int totalrec = 0;
                int totalmig = 0;
                for (int i = 0; x * i < stripenum; ++i) {
                    auto[rec, mig] = calculate_crosscluster(res1, paraBefore, paraAfter, i * x, x, false, false);
                    totalrec += rec;
                    totalmig += mig;
                }
                cout << totalrec << "-" << totalmig << endl;
            }

            //case g -> x*g
            std::cout << "g x times :k-l-g-c-s-x:" << k << "-" << l << "-" << g << "-" << c << "-" << stripenum << "-" << x
                      << "\n";
            paraAfter = make_tuple(x * k, x * l, x*g);
            cout << "random with partial decoding :\n";
            {
                int totalrec = 0;
                int totalmig = 0;
                for (int i = 0; x * i < stripenum; ++i) {
                    auto[rec, mig] = calculate_crosscluster(res1, paraBefore, paraAfter, i * x, x);
                    totalrec += rec;
                    totalmig += mig;
                }
                cout << totalrec << "-" << totalmig << endl;
            }

            cout << "random without partial decoding :\n";
            {
                int totalrec = 0;
                int totalmig = 0;
                for (int i = 0; x * i < stripenum; ++i) {
                    auto[rec, mig] = calculate_crosscluster(res1, paraBefore, paraAfter, i * x, x, false, false);
                    totalrec += rec;
                    totalmig += mig;
                }
                cout << totalrec << "-" << totalmig << endl;
            }
        }

    }

    return 0;
}
