#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <fstream>
#include <iostream>
struct Q {long pos;int id;bool end;};
struct S {unsigned long long n=0,depth=0,ge4=0;};
int main(int argc,char**argv){
 if(argc!=4)return 2;
 std::ifstream regions(argv[2]);std::string chr;long s,e;int id;int maxid=-1;
 std::map<std::string,std::vector<Q>> qs;
 while(regions>>chr>>s>>e>>id){if(s<0||e<=s)return 3;qs[chr].push_back({s,id,false});qs[chr].push_back({e,id,true});maxid=std::max(maxid,id);}
 for(auto &v:qs)std::sort(v.second.begin(),v.second.end(),[](Q a,Q b){return a.pos<b.pos;});
 std::vector<S> left(maxid+1),right(maxid+1);
 FILE* in=fopen(argv[1],"rb");if(!in)return 4;
 std::vector<char> buf(8*1024*1024);setvbuf(in,buf.data(),_IOFBF,buf.size());
 char line[4096];std::string current="";size_t j=0;S running;long prev=-1;unsigned long long nrows=0;
 auto flush=[&](long pos){auto it=qs.find(current);if(it==qs.end())return;auto &v=it->second;while(j<v.size()&&v[j].pos<=pos){if(v[j].end)right[v[j].id]=running;else left[v[j].id]=running;++j;}};
 while(fgets(line,sizeof(line),in)){
  ++nrows;char*t=strchr(line,'\t');if(!t)continue;*t=0;
  if(current!=line){if(current!="")flush(9223372036854775807L);current=line;j=0;running=S();prev=-1;}
  if(!qs.count(current))continue;
  char*p=t+1,*z;long pos=strtol(p,&z,10)-1;
  if(pos<prev){std::cerr<<"Non-monotone CpG report\n";return 5;}prev=pos;
  flush(pos);p=strchr(z+1,'\t');if(!p)return 6;++p;
  unsigned long long m=strtoull(p,&z,10);unsigned long long u=strtoull(z+1,&z,10);auto d=m+u;
  running.n++;running.depth+=d;running.ge4+=(d>=4);
 }
 flush(9223372036854775807L);fclose(in);
 FILE*out=fopen(argv[3],"wb");if(!out)return 7;
 fprintf(out,"region_id\tcytosines\tdepth_sum\tge4\n");
 for(int i=0;i<=maxid;i++)fprintf(out,"%d\t%llu\t%llu\t%llu\n",i,right[i].n-left[i].n,right[i].depth-left[i].depth,right[i].ge4-left[i].ge4);
 fclose(out);std::cerr<<"rows="<<nrows<<" intervals="<<maxid+1<<"\n";return 0;
}
