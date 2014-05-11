#include"main.h"
#include "bpRead.h"
#include "bpTree.h"
#include "CONFIG.h"
extern uint KMER_LENGTH;
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
		
		if(feature_str &&(max_root->tail->getFreq()>0) && max_root->tail->getGetFirstRead()->pdat->getStr(KMER_LENGTH-1)->equalTo(feature_str)){
			
			//cout<<"此KMer与已开始序列相同\n";
			//如果前k-1个相同，不用计算，因为后面还有一次累加
		}else{
			//否则，新结果进入counter4
			//cout<<"此KMer与已开始序列不同\n";
			counter[4]+=max_root->tail->getFreq();
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
					counter[i]+=tmp_root->child[i]->tail->getFreq();
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
			cout<<"\nEND"<<endl;
			return NULL;
		}
		if(idx_max==4){
			//cout<<"\nbegin new DebruijnNode with "<<counter_max<<" KMers :"<<endl;
			cout<<"\nNEW CONTIG BEGIN--->>\n";
			max_root->tail->getGetFirstRead()->pdat->getStr(KMER_LENGTH)->print();
			//cout<<"\n";
			joint_times=0;
			return max_root->tail->getGetFirstRead()->pdat->getStr(KMER_LENGTH);
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
		cout<<"\nEND"<<endl;
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
			DBGNode * lnk_header=tree->popLeavesLnk(node);
			bpRead * tmp_read2=lnk_header->reads->pop();
			while(tmp_read2){
				tmp_read2->start();
				reads->append(tmp_read2);
				tmp_read2=lnk_header->reads->pop();
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
				//myTree->hangOnTree(tmp_read->pdat);
				myTree->appendDBGraph(tmp_read->pdat);
			}catch(ERROR_EndOfFile){
				cout<<"EOF"<<endl;
				cout<<"<DBGraph Build Completed.>\n";
				cout<<"Kmer length:\t\t\t"<<KMER_LENGTH<<endl;
				cout<<"Nodes Number:\t\t\t"<<myTree->number_DBGNodes<<endl;
				cout<<"Nodes Weight:\t\t\t"<<myTree->weight_DBGNodes<<endl;
				cout<<"Edges Number:\t\t\t"<<myTree->number_DBGEdges<<endl;
				cout<<"Edges Weight:\t\t\t"<<myTree->weight_DBGEdges<<endl;
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
		system("PAUSE");

		bpTreeNode * winnode=myTree->root->getWinner(false);
		cout<<"\nResult:W="<<winnode->tail->getFreq()<<"\t";
		winnode->tail->printKmer(0,KMER_LENGTH-1);		
		LnkHeader<DBGNode> * start_node_list=new LnkHeader<DBGNode>;
		myTree->root->findEularStartNodes(start_node_list);		
		while(1){
			//myTree->root->print();

			DBGNode * EularNode=start_node_list->pop();
			if(!EularNode)break;
			cout<<"\nEularNodeChoose ";
			
			EularNode->printKmer();
			cout<<endl;
			if(EularNode->edges[1]->length>0){
				cout<<"Broken..choose next."<<endl;
				continue;
			}
			//winnode->tail->GreedySearch();
			RingDestroyer * destroyer=new RingDestroyer();
			cout<<"Ring Destroyer Ready..."<<endl;
			system("PAUSE");
			destroyer->DFS(EularNode,0,NULL);
			cout<<"Ring Destroy Completed...Reset the DBGraph."<<endl;
			myTree->root->startAll();
 			EularNode->EulerSearch();
		}
		system("PAUSE");
		while(1){
			bpString * seq=table->chooseNextKMER(myTree);
			if(seq==NULL)break;
			table->pushForword(seq,myTree);
		}
		return dat_linked;
	}
};
#define CMD_READ 1
#define CMD_OPEN 2
#define CMD_SET_KMER_LENGTH 3
class ShellManager{
public:
	char tmp_str[100];
	int tmp_int;
	int getCmd(){
		cin>>tmp_str;
		if(strcmp(tmp_str,"read")==0)return CMD_READ;
		if(strcmp(tmp_str,"open")==0)return CMD_OPEN;
		if(strcmp(tmp_str,"kmer")==0)return CMD_SET_KMER_LENGTH;
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
		case CMD_SET_KMER_LENGTH:{
			uint len=10;
			cin>>len;
			KMER_LENGTH=len;
			cout<<"\nset Kmer-length to "<<len<<endl;
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