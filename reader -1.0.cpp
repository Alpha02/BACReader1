#include <iostream>
#include <string.h>
#include <fstream>
#include <time.h>
#define  uint unsigned int
#define  uchar unsigned char

using namespace std;
class ERROR_BP_N{};
class ERROR_CannotOpenFile{};
class ERROR_NoFileToRead{};
class ERROR_EndOfFile{};
class ERROR_HashTableOutOfRange{};
class ERROR_KMerOutOfReads{};
class bpString{
public:
	uint length;
	uint str_len;
	uchar * str;
	bpString(){
		str=NULL;
		length=0;
	}
	void initStr(uint src_length=0){
		if(str){
			delete [str_len]str;
		}
		str_len=(src_length+3)/4;
		if(str_len){
			str=new uchar[str_len];	
		}else{
			str=NULL;
		}
	}
	uchar getBP(uint index){
		uint str_idx=index/4;
		uchar bit_shift=(index%4)*2;
		uchar ch=(str[str_idx]>>(bit_shift))&0x03;
		return ch;
	}
	void writeBP(uint index,uchar bp){
		uint str_idx=index/4;
		uchar bit_shift=(index%4)*2;
		str[str_idx]=str[str_idx]&(~(0x03<<bit_shift))|(bp<<bit_shift);
	}
	void write(char * src){
		initStr(strlen(src));
		for(int i=0;i<length;i++){
			switch(src[i]){
			case 'A':writeBP(i,0);break;
			case 'T':writeBP(i,1);break;
			case 'C':writeBP(i,2);break;
			case 'G':writeBP(i,3);break;
			default:throw ERROR_BP_N();break;
			}
		}
	}
	bpString(uchar * src,const uint len){
		str=src;
		length=len;
		str_len=(length+3)/4;
	}
	void copy(bpString * src){
		initStr(src->length);
		for(int i=0;i<src->str_len;i++){
			str[i]=src->str[i];
		}
	}
	/*
	int getOverlapLength(bString * target){
		int id_tail=length-1;
		int id_head=0;
		bool id_found=false;
		bool target_completed=false;
		while(id_tail>0){
			int id_tail_temp=id_tail;
			while((id_tail>=0) && (str[id_tail]!=target->str[id_head])){
				id_tail--;
			}
			if(id_tail<0)break; 
			id_tail_temp=id_tail;
			//记录暂时位点。从这里向后匹配。三种情况。
			//1.A串无法匹配到尾部，则放弃该匹配。
			//2.A串匹配到尾部。则接受该匹配。
			//3.B串被匹配至尾部，则认为B串是A串的子串。故应该丢弃B串。
			while(id_tail<length){
				if(str[id_tail]!=target->str[id_head])break;
				id_tail++;
				id_head++;
				if(id_head>=target->length){
					target_completed=true;
					break;
				}
				if(id_tail>=length){
					id_found=true;
					break;
				}
			}
			id_tail=id_tail_temp-1;
			id_head=0;
			if(id_found || target_completed){
				return length-id_tail_temp;
			}
		}
		return 0;
	}
	*/
};
#define READ_LEN_ID 50
#define READ_LEN_SEQ 120
class bpRead{
public:
	uchar seqID[READ_LEN_ID];
	bpString seq;
	uchar seq_src[READ_LEN_SEQ];
	uchar quality[READ_LEN_SEQ];

};
template <class T> 
class Lnk{
public:
	T * pdat;
	Lnk<T> *next;
	Lnk<T>(){
		next=NULL;
		pdat=NULL;
	}
	Lnk<T>(T * new_pdat){
		pdat=new_pdat;
		next=NULL;
	}
	Lnk<T>(T* new_pdat,Lnk<T> *nextnode){
		pdat=new_pdat;
		next=nextnode;
	}
};
#define STA_STOP 0
#define STA_START 1
#define STA_PAUSE 2
#define STA_FAILED 3
class readZipper{
public:
	FASTQData * data;
	int pos_current;
	int pos_last;
	int status;
	int k;
	readZipper(int k_new=4){
		pos_current=0;
		pos_last=0;
		status=STA_STOP;
		k=k_new;
	}
	readZipper(FASTQData * src_data,int k_new=4){
		data=src_data;
		pos_current=0;
		pos_last=0;
		status=STA_STOP;
		k=k_new;
	}
	bool fitKMer(K_Mer * kmer){
		//是否能够作为某条kmer的后继
		K_Mer new_kmer(data,pos_current,k);
		if(kmer->fitTo(&new_kmer))return true;
		return false;
	}
	bool fitTo(readZipper * read2){
		for(int i=1;i<k;i++){
			if(this->data->seq.str[pos_current+i-1]!=read2->data->seq.str[read2->pos_current+i])return false;
		}
		return true;
	}
	int setStatus(int pos_cur,K_Mer * kmer){                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
		if(fitKMer(kmer)){
			status=STA_START;
		}else if(pos_cur-pos_last<k){
			status=STA_PAUSE;
		}else{
			status=STA_FAILED;
		}
		return status;
	}
	int getHash(){
		int sum=0;
		for(int i=0;i<k;i++){
			sum=sum*4+data->seq.str[pos_current+i]-1;
		}
		return sum;
	}
	void pushForward(){
		pos_current++;
		pos_last=pos_current;
	}
};
class decisionTable{
	//决策表是这样一个对象
	//它保存着所有正在拼接的reads的指针，及其当前拼接位置、最后成功拼接位置、拼接状态
	//这些reads指针都是动态的。以链表组织的。
	linkedList<readZipper> * reads;
	decisionTable(){
		reads=NULL;
	}
	void addRead(readZipper * src){
		linkedList<readZipper> * new_read=new linkedList<readZipper>(src,reads);
		reads=new_read;
	}
	void readNext(){
	}
	
};
class HashTable{
public:
	int k;
	int width;
	linkedList<readZipper> ** link_headers;
	int *link_length;
	HashTable(int ex,int idx_begin=0){
		k=ex;
		width=1<<(2*ex);
		link_headers=new linkedList<readZipper> * [width];
		link_length=new int[width];
		for(int i=0;i<width;i++){
			link_headers[i]=NULL;
			link_length[i]=0;
		}
	}
	HashTable(HashTable & src){
		k=src.k;
		width=src.width;
		link_headers=new linkedList<readZipper> * [width];
		link_length=new int[width];
		for(int i=0;i<width;i++){
			link_headers[i]=src.link_headers[i];
			link_length[i]=src.link_length[i];
		}
	}
	void addData(readZipper * dat_new){
		int h=dat_new->getHash();
		if(h>=0 && h<width){
			linkedList<readZipper> * new_header=new linkedList<readZipper>(dat_new,link_headers[h]);
			link_headers[h]=new_header;
			link_length[h]++;
		}else{
			throw ERROR_HashTableOutOfRange();
		}
	}
	readZipper * getMaxRead(){
		int max_length=0;
		int max_index=0;
		for(int i=0;i<width;i++){
			if(link_length[i]>max_length){
				max_length=link_length[i];
				max_index=i;
			}
		}
		readZipper * my_read=link_headers[max_index]->content;
		return my_read;
	}
	HashTable * getNextTableOfKMer(readZipper * src){
		int next_pos=src->pos_current+1;
		HashTable * hashtable=new HashTable(k,next_pos);
		int src_hash=src->getHash();
		for(linkedList<readZipper> * rl=link_headers[src_hash];rl!=NULL;rl=rl->next){
			rl->content->pushForward();
			hashtable->addData(rl->content);
		}
		src_hash=(src_hash<<2)&((1<<(2*k))-1);
		for(int i=src_hash;i<src_hash+4;i++){
			//链表啊
			if=link_headers[i]>content->fitTo(src)){
				hashtable->addData(link_headers[i]->content);
			}
		}
		return hashtable;
	}

};
struct OLCNodeLinked;
class OLCNode{
	FASTQData * dat;
	OLCNodeLinked * in_linked;
	OLCNodeLinked * out_linked;
	OLCNode(){
		in_linked=NULL;
		out_linked=NULL;
		dat=NULL;
	}
	void addEdge(OLCNode * src,bool edge_direction,int weight);
};
struct OLCNodeLinked{
	OLCNode * node;
	OLCNodeLinked * next;
	int weight;
};
void OLCNode::addEdge(OLCNode * src,bool edge_direction,int weight=0){
	//direction=EDGE_IN src--->this
	OLCNodeLinked * newEdge=new OLCNodeLinked();
	newEdge->weight=weight;
	newEdge->node=src;
	if(edge_direction==EDGE_IN){
		newEdge->next=in_linked;
		in_linked=newEdge;
	}else{
		newEdge->next=out_linked;
		out_linked=newEdge;
	}
}
class OverlapGraph{

};
class FASTQReader{
private:
	ifstream * fin;
	FASTQLinked * dat_linked;
public:
	FASTQReader(){
		fin=NULL;
	}
	void open(char * file_name){
		fin=new ifstream;
		fin->open(file_name);
		if(fin->fail()){
			throw ERROR_CannotOpenFile();
		}
	}
	void close(){
		fin->close();
	}
	FASTQReader(char * file_name){
		open(file_name);
	}
	FASTQData * readSeq(){
		if((!fin) || (!fin->is_open()))throw ERROR_NoFileToRead();
		FASTQData *dat=new FASTQData;
		fin->getline(dat->seqID,FASTQ_LEN_ID);
		char seqstr[FASTQ_LEN_SEQ];
		fin->getline(seqstr,FASTQ_LEN_SEQ);
		dat->seq.read(seqstr);
		fin->getline(dat->quality,100);
		fin->getline(dat->quality,FASTQ_LEN_QUA);
		if(fin->eof())throw ERROR_EndOfFile();
		return dat;
	}
	FASTQLinked * readAllSeq(){
		if((!fin) || (!fin->is_open()))throw ERROR_NoFileToRead();
		int i=0;
		long long  num=0;
		clock_t t_start,t_end;
		t_start=clock();
		HashTable * hash_table=new HashTable(4);
		dat_linked=new FASTQLinked;
		dat_linked->dat=new FASTQData;
		dat_linked->next=NULL;
		strcpy(dat_linked->dat->seqID,"HEADER");
		FASTQLinked *tmp_linked=dat_linked;
		int a=0;
		int q[100];
		for(i=0;i<100;i++)q[i]=0;
		while(!fin->eof()){
			i++;
			if(i%1000==0)cout<<".";
			try{
			tmp_linked->next=new FASTQLinked;
			tmp_linked->next->dat=readSeq();
			readZipper * new_read=new readZipper(tmp_linked->next->dat);
			hash_table->addData(new_read);
			tmp_linked->next->next=NULL;
			}catch(ERROR_EndOfFile){
				cout<<"EOF"<<endl;
				break;
			}
			int a=0;
			if(dat_linked!=tmp_linked){
				a=tmp_linked->dat->seq.getOverlapLength(&tmp_linked->next->dat->seq);
			}
			q[a]+=1;
			num++;
			tmp_linked=tmp_linked->next;
		}
		t_end=clock();
		for(i=0;i<100;i++){
			cout<<endl;
			cout<<i<<":"<<q[i];
		}
		cout<<endl<<num;
		double t_diff=((double)(t_end-t_start))/CLOCKS_PER_SEC;
		cout<<"Read Complete in "<<t_diff<<" s "<<endl;
		while(1){
			cout<<"Hash table iteration begin"<<endl;
			cout<<"Max:";
			readZipper * tmp_read=hash_table->getMaxRead();
			cout<<endl;
			hash_table=hash_table->getNextTableOfKMer(tmp_read);
		}
		return dat_linked;
	}
};
#define CMD_READ 1
#define CMD_OPEN 2
class ShellManager{
public:
	char tmp_str[100];
	int tmp_int;
	int getCmd(){
		cin>>tmp_str;
		if(strcmp(tmp_str,"read")==0)return CMD_READ;
		if(strcmp(tmp_str,"open")==0)return CMD_OPEN;
		return 0;
	}
};
int main(){
	bString bstr;
	ShellManager sh;
	FASTQReader reader;
	while(1){
		int cmd=sh.getCmd();
		switch(cmd){
		case CMD_READ:{
			try{
			FASTQLinked * link=reader.readAllSeq();
			}catch(ERROR_NoFileToRead){
				cout<<"[ERROR]You should open a file!"<<endl;
			}
			break;
					  }
		case CMD_OPEN:{
			char tmpstr[100];
			cin>>tmpstr;
			try{
				reader.open(tmpstr);
			}catch(ERROR_CannotOpenFile){
				cout<<"[ERROR]Failed to open the file!"<<endl;
			}
			break;
					  }
		default:{
			cout<<"[ERROR]Invalid command!"<<endl;
				break;
				}
		}
	}
	return 0;
}