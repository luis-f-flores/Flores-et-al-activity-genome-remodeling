#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
struct Region {string chr,id; long start,end; unsigned long long cs=0,ts=0,n=0;};
int main(int argc,char**argv){
 if(argc!=6){cerr<<"usage: count report regions upper_cut output_prefix tracks\n";return 2;}
 vector<Region> regs; map<string,vector<int>> ix;
 ifstream rr(argv[2]); string line;
 while(getline(rr,line)){istringstream s(line);Region r;s>>r.chr>>r.start>>r.end>>r.id;if(!s) return 3;ix[r.chr].push_back(regs.size());regs.push_back(r);}
 for(auto &kv:ix) sort(kv.second.begin(),kv.second.end(),[&](int a,int b){return regs[a].start<regs[b].start;});
 map<string,vector<long>> maxend;for(auto &kv:ix){long mx=0;for(int i:kv.second){mx=max(mx,regs[i].end);maxend[kv.first].push_back(mx);}}
 ifstream in(argv[1]); if(!in)return 4;double cut=atof(argv[3]);
 ofstream tr(argv[5]);tr<<"chr\tposition\tstrand\tnumCs\tnumTs\tcoverage\n";
 unsigned long long total=0,mincov=0,filtered=0;string prevchr;long prevpos=0;bool ordered=true;
 while(getline(in,line)){
  ++total; char* p=line.data();char* e=strchr(p,'\t');if(!e)return 5;*e=0;string chr(p);p=e+1;
  long pos=strtol(p,&e,10);p=e+1;char strand=*p;p=strchr(p,'\t')+1;
  unsigned long long cs=strtoull(p,&e,10);p=e+1;unsigned long long ts=strtoull(p,&e,10);auto cov=cs+ts;
  if(cov<4)continue;++mincov;if(cov>=cut)continue;++filtered;
  if(chr==prevchr&&pos<prevpos)ordered=false;prevchr=chr;prevpos=pos;
  if((chr=="chr2"&&pos>=109692150&&pos<=109693150)||(chr=="chr14"&&pos>=45349099&&pos<=45353609))tr<<chr<<'\t'<<pos<<'\t'<<strand<<'\t'<<cs<<'\t'<<ts<<'\t'<<cov<<'\n';
  auto it=ix.find(chr);if(it==ix.end())continue;auto &v=it->second;auto &m=maxend[chr];
  auto u=upper_bound(v.begin(),v.end(),pos,[&](long q,int j){return q<regs[j].start;});
  long k=u-v.begin()-1;for(;k>=0&&m[k]>=pos;--k){auto &r=regs[v[k]];if(r.end>=pos){r.cs+=cs;r.ts+=ts;++r.n;}}
 }
 ofstream out(string(argv[4])+".tsv");out<<"region_id\tnumCs\tnumTs\tcovered_CpG_records\n";
 for(auto&r:regs)out<<r.id<<'\t'<<r.cs<<'\t'<<r.ts<<'\t'<<r.n<<'\n';
 ofstream qc(string(argv[4])+".json");qc<<"{\"raw_records\":"<<total<<",\"mincov4_records\":"<<mincov<<",\"filtered_records\":"<<filtered<<",\"upper_cut\":"<<cut<<",\"within_chromosome_sorted\":"<<(ordered?"true":"false")<<"}\n";
 return 0;
}
