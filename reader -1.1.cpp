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
class ERROR_TreeNodeChildOutOfRange{};
class ERROR_BPStringOutOfRange{};
template <class T> 
class Lnk{
public:
	T * pdat;
	Lnk<T> *next;
	Lnk<T> *prev;
	Lnk<T>(){
		next=NULL;
		prev=NULL;
		pdat=NULL;
	}
	Lnk<T>(T * new_pdat){
		pdat=new_pdat;
		next=NULL;
		prev=NULL;
	}
	Lnk<T>(T* new_pdat,Lnk<T> *nextnode){
		pdat=new_pdat;
		next=nextnode;
		prev=NULL;
	}
	void append(T * new_pdat){
		if(pdat==NULL){
			pdat=new_pdat;
			return;
		}
		//0<=====>0
		//    0-->0
		Lnk<T> * new_Lnk=new Lnk<T>(pdat,next);
		pdat=new_pdat;
		//0------>0
		//    0<=>0
		if(next)next->prev=new_Lnk;
		//0       0
		//0-->0<=>0
		next=new_Lnk;
		//0       0
		//0<=>0<=>0
		new_Lnk->prev=this;
	}
	void appendFromLnk(Lnk<T> * src){
		while(src){
			append(src->pdat);
			src=src->next;
		}
	}
	uint getLength(){
		Lnk<T> * tmp=this;
		uint len=0;
		while(tmp){
			len++;
			tmp=tmp->next;
		}
		return len;
	}
	~Lnk(){
		if(next)delete next;
	}
	void remove(){
		if(prev)prev->next=next;
		if(next)next->prev=prev;
		next=NULL;
		delete this;
	}
};
class ERROR_PopFromEmptyLnk{};
template <class T> 
class LnkHeader{
public:
	Lnk<T> *lnk;
	uint length;
	LnkHeader(){
		length=0;
		lnk=NULL;
	}
	void append(T*new_pdat){
		if(!lnk){
			lnk=new Lnk<T>(new_pdat);
		}else{
			lnk->append(new_pdat);
		}
		length++;
	}
	void appendFromLnk(Lnk<T> *src){
		while(src){
			append(src->pdat);
			src=src->next;
		}
	}
	T * pop(){
		if(!lnk)return NULL;
		T * res=lnk->pdat;
		Lnk<T> *tmp=lnk;
		
		lnk=lnk->next;
		if(lnk){
			lnk->prev=tmp->prev;
		}
		tmp->next=NULL;
		delete tmp;
		length--;
		return res;
	}
	LnkHeader(T * src){
		length=0;
		lnk=NULL;
		append(src);
	}
	~LnkHeader(){
		delete lnk;
	}
	void clear(){
		delete lnk;
		lnk=NULL;
		length=0;
	}
};
class DBGKmer{
	uint idx;
	uint length;
	LnkHeader<bpRead>* reads;

};
//未访问的点
#define DBG_START 0
//正在访问的点（到达过，但是还有出路可寻）
#define DBG_VISITING 1
//访问结束的点（到达过，且所有出路已遍历）
#define DBG_OVER 2
template <class T>
#define 
class DBGNode{
	LnkHeader<DBGNode> * in;
	LnkHeader<DBGNode> * out;
	T * dat;
	uint flag;
};
class bpString{
public:
	uint length;
	uint str_len;
	uchar * str;
	bpString(){
		str=NULL;
		length=0;
		str_len=0;
	}
	void initStr(uint src_length=0){
		if(str){
			delete [str_len]str;
		}
		length=src_length;
		str_len=(src_length+3)/4;
		if(str_len){
			str=new uchar[str_len];	
		}else{
			str=NULL;
		}
	}
	uchar getBP(uint index){
		uint str_idx=index/4;
		if(index>length || str_idx>str_len)throw ERROR_BPStringOutOfRange();		
		uchar bit_shift=(index%4)*2;
		uchar ch=(str[str_idx]>>(bit_shift))&0x03;
		return ch;
	}
	void writeBP(uint index,uchar bp){
		uint str_idx=index/4;
		uchar bit_shift=(index%4)*2;
		str[str_idx]=str[str_idx]&(~(0x03<<bit_shift))|(bp<<bit_shift);
	}
	void write(char * src,int len=-1,uint idx_begin=0){
		if(len==-1)len=strlen(src);
		initStr(len);
		for(uint i=0;i<length;i++){
			try{
				switch(src[i+idx_begin]){
				case 'A':writeBP(i,0);break;
				case 'T':writeBP(i,1);break;
				case 'C':writeBP(i,2);break;
				case 'G':writeBP(i,3);break;
				default:throw ERROR_BP_N();break;
				}
			}catch(ERROR_BP_N){
				cout<<"N";
				writeBP(i,0);
			}
		}
	}
	bpString(char * src){
		str=NULL;
		write(src);
	}
	bpString(char * src,uint len,uint idx_begin=0){
		str=NULL;
		write(src,len,idx_begin);
	}
	void copy(bpString * src){
		initStr(src->length);
		for(uint i=0;i<src->str_len;i++){
			str[i]=src->str[i];
		}
	}
	void print(){
		for(uint i=0;i<length;i++){
			switch(getBP(i)){
			case 0:cout<<'A';break;
			case 1:cout<<'T';break;
			case 2:cout<<'C';break;
			case 3:cout<<'G';break;
			default:throw ERROR_BP_N();break;
			}
		}
	}
	bool equalTo(bpString * target){
		if(target->length!=length)return false;
		for(int i=0;i<length;i++){
			if(target->getBP(i)!=getBP(i))return false;
		}
		return true;
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
#define STA_STOP 0
#define STA_START 1
#define STA_PAUSE 2
#define STA_END 3
const uint KMER_LENGTH=30;
class bpRead{
public:
	char seqID[READ_LEN_ID];
	bpString seq;
	char seq_src[READ_LEN_SEQ];
	char quality[READ_LEN_SEQ];
	uint pos_last;
	uint pos_cur;
	uchar status;
	void compileSeq(){
		seq.write(seq_src);
	}
	bpRead(){
		pos_last=0;
		pos_cur=0;
		status=STA_STOP;
	}
	uchar getBP(int idx){
		return seq.getBP(idx);
	}
	uchar getBPOffset(int offset){
		return seq.getBP(pos_cur+offset);
	}
	void pushForword(){
		pos_cur++;
		if(status==STA_START){
			pos_last=pos_cur;
		}else{
		}
		if(pos_cur>=seq.length-KMER_LENGTH+1){
			status=STA_END;
		}
	}
	bpString * getStr(int length){
		return new bpString(seq_src,length,pos_cur);
	}
	bool jointERROR(uint k=KMER_LENGTH){
		if(pos_cur-pos_last>k)return true;
		return false;
	}
	void start(){
		status=STA_START;
		pos_last=pos_cur;
	}
	void pause(){
		status=STA_PAUSE;
	}
	void stop(){
		status=STA_STOP;
		pos_cur=0;
		pos_last=0;
	}
};
class bpTreeNode{
public:
	LnkHeader<bpRead> * tail;
	bpTreeNode * winner_child;
	bpTreeNode * child[4];
	bpTreeNode * parent;
	bpTreeNode(bpTreeNode * p=NULL,LnkHeader<bpRead> * t=NULL){
		child[0]=NULL;
		child[1]=NULL;
		child[2]=NULL;
		child[3]=NULL;
		parent=p;
		tail=t;
		winner_child=NULL;
	}
	bpTreeNode * getWinner(){
		if(tail!=NULL){
			return this;
		}
		if(winner_child!=NULL){
			return winner_child;
		}
		uint max_length=0;
		uchar max_idx=0;
		for(uchar i=0;i<4;i++){
			if(!child[i])continue;
			bpTreeNode * winner_tmp=child[i]->getWinner();
			if((!winner_child)||(winner_tmp && (winner_tmp->tail->length>winner_child->tail->length))){
				winner_child=winner_tmp;
			}
		}
		return winner_child;
	}
	void refreshWinner(){
		if(!parent)return;
			uint max_length=0;
			uchar max_idx=0;
			for(uchar i=0;i<4;i++){
				if(!parent->child[i])continue;
				bpTreeNode * winner_tmp=parent->child[i]->winner_child;
				if(winner_tmp && (winner_tmp->tail->length>parent->winner_child->tail->length)){
					parent->winner_child=winner_tmp;
				}
		}
		parent->refreshWinner();	
	}
};
class ERROR_TreeNodeNotFound{};

class bpTree{
public:
	bpTreeNode * root;
	bpTreeNode * winner;
	LnkHeader<bpRead> * Lnk_Longest;
	LnkHeader<bpRead> * Lnk_Decision;
	bpTree(){
		root=new bpTreeNode();
		Lnk_Longest=NULL;
		winner=NULL;
	}
	void hangOnTree(bpRead * src){
		//将一条read放进树里。Kmer有多长，树有多高。
		bpTreeNode * node=root;
		bpTreeNode * node_next;
		for(int i=0;i<KMER_LENGTH;i++){
			uchar child_idx=src->getBPOffset(i);
			node_next=node->child[child_idx];
			if(node_next==NULL){
				node_next=new bpTreeNode(node);
				node->child[child_idx]=node_next;
			}
			node=node_next;
		}
		if(!node->tail){
			node->tail=new LnkHeader<bpRead>(src);
		}else{
			node->tail->append(src);
		}
		//更新链表最长的链表头记录。
		if(!Lnk_Longest || (Lnk_Longest->length)<(node->tail->length))
			Lnk_Longest=node->tail;
	}
	LnkHeader<bpRead> * catchAllReadFromNode(bpTreeNode * root,LnkHeader<bpRead> * fishingRod){
		//递归找到节点下所有的Read
		for(int i=0;i<4;i++){
			if(root->child[i]==NULL)continue;
			catchAllReadFromNode(root->child[i],fishingRod);
		}
		if(root->tail!=NULL){
			fishingRod->appendFromLnk(root->tail->lnk);
			cout<<"\napp:";
			root->tail->lnk->pdat->seq.print();
			cout<<endl;
		}
		return fishingRod;
	}
	bpTreeNode * searchNode(bpString * seq,uint begin_idx=0,uint end_idx=KMER_LENGTH-1){
		if(seq==NULL)return NULL;
		bpTreeNode * node=root;
		for(uint i=begin_idx;i<=end_idx;i++){
			uchar child_idx=seq->getBP(i);
			node=node->child[child_idx];
			if(node==NULL){
				throw ERROR_TreeNodeNotFound();
			}
		}
		return node;
	}

	LnkHeader<bpRead> * popLeavesLnk(bpTreeNode * target){
		LnkHeader<bpRead> * lnk=target->tail;
		return lnk;
	}
	LnkHeader<bpRead> * search(bpString * seq,uint begin_idx=0,uint end_idx=KMER_LENGTH-1){
		//寻找read根据seq,结果为reads的链表。为所有可行的寻找结果。
		bpTreeNode * node=searchNode(seq,begin_idx,end_idx);
		LnkHeader<bpRead> * fishingRod=new LnkHeader<bpRead>;
		catchAllReadFromNode(node,fishingRod);
		cout<<"Search result:"<<fishingRod->length<<" found."<<endl;
	}
	void run(){		
		bpRead * bprd;
		//先处理上次决定开始的那些正在拼接的read
		while(Lnk_Decision->lnk){
			bprd=Lnk_Decision->pop();
		}
		//再找到当前最长read
		while(Lnk_Longest->lnk){
			//找到整条链表销毁，把节点挂到其它链表上
			//bprd=mylnk->pop();
			//bprd->pushForword();
			//lnk_decision->append(bprd);
		}
		//找到特征kmer的前k-1个bp值
		bpString * seq=bprd->getStr(KMER_LENGTH-1);
		//把这些值所在的树展开
		LnkHeader<bpRead> * newlnk=search(seq,0,KMER_LENGTH-1);
//		lnk_decision->appendFromLnk(newlnk);
		newlnk->clear();
	}
};

class bpDecisionTable{
public:
	LnkHeader<bpRead> * reads;
	int counter[5];
	uint joint_times;
	//counter的前4种决策，来源于正在拼接的节点。后1种决策，来源于最多可加入的节点（如果它们有不同的序列）
	bpDecisionTable(){
		reads=new LnkHeader<bpRead>;
		joint_times=0;
	}
	bpString * chooseNextKMER(bpTree * tree){
		//计数器清零
		int i;
		for(i=0;i<5;i++)counter[i]=0;
		//累加已开始的拼接能够继续的数量（假定所有拼接都维持原状态向前推进1位）
		Lnk<bpRead> * tmp_read=reads->lnk;
		Lnk<bpRead> * read_started=reads->lnk;
		bpString * feature_str=NULL;
		while(read_started && read_started->pdat->status!=STA_START){
			if(read_started->next && (read_started->next->prev!=read_started)){
				cout<<"Bingo!"<<endl;

			}
			read_started=read_started->next;
		}
		//如果找到了已经开始的点
		if(read_started){
			//记录特征前缀，用于寻找可以开始拼接的节点。
			feature_str=read_started->pdat->getStr(KMER_LENGTH-1);
			while(tmp_read){
				//是否开始，如果是的话
				if(tmp_read->pdat->status==STA_START){
					//比较前缀，如果相等，则进行统计
					bpString * target_str=tmp_read->pdat->getStr(KMER_LENGTH-1);
					if(feature_str->equalTo(target_str)){
						uchar tmp_bp=tmp_read->pdat->getBPOffset(KMER_LENGTH-1);
						counter[tmp_bp]++;
					}
				}
				tmp_read=tmp_read->next;
			}
		}
		//处理完已经开始的read，之后要收集可以加入的新read。收集有两种方式。
		//如果当前已经没有正在拼接的read，则寻找树中链表最长的read组来开始。
		//否则，既要寻找最长的，作为counter4，又要寻找已知的，作为counter0~3
		//累加可以开始继续的数量。
		bpTreeNode * max_root=tree->root->getWinner();
		
		//cout<<"最多可加入KMer:";
		//max_root->tail->lnk->pdat->getStr(KMER_LENGTH)->print();
		//cout<<"\n数量："<<max_root->tail->length<<"\n";
		
		if(feature_str &&(max_root->tail->length>0) && max_root->tail->lnk->pdat->getStr(KMER_LENGTH-1)->equalTo(feature_str)){
			
			//cout<<"此KMer与已开始序列相同\n";
			//如果前k-1个相同，不用计算，因为后面还有一次累加
		}else{
			//否则，新结果进入counter4
			//cout<<"此KMer与已开始序列不同\n";
			counter[4]+=max_root->tail->length;
		}
		bpTreeNode * tmp_root;
		try{
			//cout<<"寻找带有前缀";
			//if(feature_str)feature_str->print();
			//cout<<"的新序列\n";
		    tmp_root=tree->searchNode(feature_str,0,KMER_LENGTH-2);
		}catch(ERROR_TreeNodeNotFound){
			//cout<<"寻找失败。\n";
			tmp_root=NULL;
		}
		//寻找kmer前缀与feature_str相同的节点
		if(tmp_root){
			for(i=0;i<4;i++){
				//分别统计子节点
				if(tmp_root->child[i]){
					counter[i]+=tmp_root->child[i]->tail->length;
					//cout<<"发现序列组["<<i<<"],长度："<<tmp_root->child[i]->tail->length<<endl;
				}
			}
		}
		//现在可以计算哪个counter最大
		uchar idx_max=0;
		uint counter_max=0;
		//如果存在可以继续拼接的kemer，优先考虑它们。只有当不能继续拼接时，才考虑开始下一条。
		for(i=0;i<4;i++){
			if(counter[i]>counter_max){
				counter_max=counter[i];
				idx_max=i;
			}
		}
		if(counter_max==0){
			counter_max=counter[4];
			idx_max=4;
		}
		if(counter_max==0){
			cout<<"\END"<<endl;
			return NULL;
		}
		if(idx_max==4){
			//cout<<"\nbegin new DebruijnNode with "<<counter_max<<" KMers :"<<endl;
			cout<<"\nNEW CONTIG BEGIN--->>\n";
			max_root->tail->lnk->pdat->getStr(KMER_LENGTH)->print();
			//cout<<"\n";
			joint_times=0;
			return max_root->tail->lnk->pdat->getStr(KMER_LENGTH);
		}
		else{
			joint_times++;
			//行输出
			//cout<<"\ntimes:"<<counter_max<<"\n";
			//for(int j=0;j<joint_times;j++)cout<<"-";
			//单碱基输出
			bpString * bp=new bpString("A");
			bp->writeBP(0,idx_max);
			bp->print();
			bpString * new_str=read_started->pdat->getStr(KMER_LENGTH);
			new_str->writeBP(KMER_LENGTH-1,idx_max);
			return new_str;
		}
		cout<<"\END"<<endl;
		return NULL;
	}
	void pushForword(bpString * featureString,bpTree * tree){
		//这里包含4个步骤
		//1.将被选中的read的status设为STA_START其余PAUSE
		//2.计算PAUSE的节点中是否有与当前选中的kmer相同的点，如果有，将其start。如果超出了，将其删除并弄回树上。
		//3.将树中被选中的点从树中删除，加入决策表并从0开始start，更新赢者树。
		//4.将所有节点前推一位。
		if(featureString){
			//featureString->print();
			//cout<<"-";
		}
		Lnk<bpRead> * tmp_read=reads->lnk;
		while(tmp_read){
			//步骤1
			if(tmp_read->pdat->getStr(KMER_LENGTH)->equalTo(featureString))
				tmp_read->pdat->status=STA_START;
			else {
			//步骤2
				tmp_read->pdat->status=STA_PAUSE;
				if(tmp_read->pdat->jointERROR()){
					//回到树里去吧。。。
					tmp_read->pdat->stop();
				}
			}

			tmp_read=tmp_read->next;
		}
		//步骤3
		bpTreeNode * node=NULL;
		try{
			node=tree->searchNode(featureString);
		}catch(ERROR_TreeNodeNotFound){
		}
		if(node){
			LnkHeader<bpRead> * lnk_header=tree->popLeavesLnk(node);
			bpRead * tmp_read2=lnk_header->pop();
			while(tmp_read2){
				tmp_read2->start();
				reads->append(tmp_read2);
				tmp_read2=lnk_header->pop();
			}
			node->refreshWinner();
		}
		tmp_read=reads->lnk;
		//步骤4
		while(tmp_read){
			tmp_read->pdat->pushForword();
			if(tmp_read->pdat->status==STA_END){
//				cout<<"\n"<<tmp_read->pdat->seqID<<" FINISHED"<<endl;
				

				Lnk<bpRead> * newtmp=tmp_read->next;
				if(reads->lnk==tmp_read){
					reads->lnk=newtmp;
				}
				tmp_read->remove();
				
				tmp_read=newtmp;
			}else{
			tmp_read=tmp_read->next;
			}
		}
	}

};
class FASTQReader{
private:
	ifstream * fin;
	Lnk<bpRead> * dat_linked;
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
	bpRead * readSeq(){
		if((!fin) || (!fin->is_open()))throw ERROR_NoFileToRead();
		bpRead *dat=new bpRead;
		fin->getline(dat->seqID,READ_LEN_ID);
		fin->getline(dat->seq_src,READ_LEN_SEQ);
		dat->compileSeq();
		fin->getline(dat->quality,100);
		fin->getline(dat->quality,READ_LEN_SEQ);
		if(fin->eof())throw ERROR_EndOfFile();
		return dat;
	}
	Lnk<bpRead> * readAllSeq(){
		if((!fin) || (!fin->is_open()))throw ERROR_NoFileToRead();
		int i=0;
		long long  num=0;
		clock_t t_start,t_end;
		t_start=clock();
		bpTree * myTree=new bpTree;
		while(!fin->eof()){
			i++;
			if(i%1000==0)cout<<".";
			try{
				Lnk<bpRead> *tmp_read=new Lnk<bpRead>;
				tmp_read->pdat=readSeq();
				myTree->hangOnTree(tmp_read->pdat);
			}catch(ERROR_EndOfFile){
				cout<<"EOF"<<endl;
				break;
			}
			num++;
		}
		t_end=clock();
		cout<<endl<<num;
		double t_diff=((double)(t_end-t_start))/CLOCKS_PER_SEC;
		cout<<"Read Complete in "<<t_diff<<" s "<<endl;
		
		bpDecisionTable * table=new bpDecisionTable();
		num=0;
		while(1){
			bpString * seq=table->chooseNextKMER(myTree);
			if(seq==NULL)break;
			table->pushForword(seq,myTree);
			
			/*
			cout<<"Hash table iteration begin"<<endl;
			bpString str("AAAA");
			
			str.write("AAAAAAAA");
			str.print();
			LnkHeader<bpRead> * tmp_read=myTree->search(&str,0,7);
			str.write("AAAAAAAT");
			str.print();
			LnkHeader<bpRead> * tmp_read2=myTree->search(&str,0,7);
			str.write("AAAAAAAC");
			str.print();
			LnkHeader<bpRead> * tmp_read3=myTree->search(&str,0,7);
			str.write("AAAAAAAG");
			str.print();
			LnkHeader<bpRead> * tmp_read4=myTree->search(&str,0,7);
			str.write("AAAAAAAA");
			str.print();
			LnkHeader<bpRead> * tmp_read5=myTree->search(&str,0,6);

			cout<<"Max:"<<myTree->Lnk_Longest->length<<" from ";
			myTree->Lnk_Longest->lnk->pdat->seq.print();
			//readZipper * tmp_read=hash_table->getMaxRead();
			//cout<<endl;
			//hash_table=hash_table->getNextTableOfKMer(tmp_read);
			*/
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
	bpString bstr;
	ShellManager sh;
	FASTQReader reader;
	while(1){
		int cmd=sh.getCmd();
		switch(cmd){
		case CMD_READ:{
			try{
			Lnk<bpRead> * link=reader.readAllSeq();
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