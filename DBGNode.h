#ifndef DBGNODE_H
#define DBGNODE_H
#include "DBGEdge.h"
#include "bpRead.h"
#define EDGE_OUT 0
#define EDGE_IN 1
class bpTreeNode;
class Kmer{
public:
	bpRead * my_read;
	uint idx;
	uint len;
	Kmer(bpRead * src_read,uint index=0,uint length=KMER_LENGTH){
		my_read=src_read;
		idx=index;
		len=length;
	}
};
class DBGNode{
public:
	LnkHeader<DBGEdge> * edges[2];
	LnkHeader<bpRead>* reads;
	LnkHeader<Kmer> * kmers;
	bpTreeNode * myTreeNode;
	uint flag;
	uint idx;
	uint getFreq(bool visit_limit=true){
		if(flag!=DBG_START)return 0;
		return kmers->length;
	}
	Lnk<bpRead>* getGetFirstRead(){
		return reads->lnk;
	}
	bpString * getStr(int length){
		return reads->lnk->pdat->getStr(length);
	}
	bpString * getSeq(){
		return &(reads->lnk->pdat->seq);
	}
	void append(bpRead * new_pdat){
		reads->append(new_pdat);
	}
	void appendKmer(bpRead * src){
		kmers->append(new Kmer(src,src->pos_cur));
	}
	DBGNode(bpTreeNode * root=NULL){
		myTreeNode=root;
		//cout<<"DBG["<<(uint)this<<"]:NodeCreate"<<endl;
		edges[0]=new LnkHeader<DBGEdge>;
		edges[1]=new LnkHeader<DBGEdge>;
		reads=new LnkHeader<bpRead>;
		kmers=new LnkHeader<Kmer>;
		flag=DBG_START;
	}
	void printKmer(int idx_begin=0,int len=KMER_LENGTH){
		for(int i=idx_begin;i<idx_begin+len;i++){
			cout<<kmers->lnk->pdat->my_read->seq_src[i+kmers->lnk->pdat->idx];
		}
	}
	void print(){
		cout<<"Node["<<idx<<"](vis="<<flag<<")";
		printKmer();
		cout<<"  addr="<<(uint)this<<"\t";
	}
	DBGNode(bpRead * src){
	//cout<<"DBG["<<(uint)this<<"]:NodeCreate With Source ";
		edges[0]=new LnkHeader<DBGEdge>;
		edges[1]=new LnkHeader<DBGEdge>;
		reads=new LnkHeader<bpRead>(src);
		kmers=new LnkHeader<Kmer>(new Kmer(src,src->pos_cur));
		//printKmer();
		//cout<<endl;
		flag=DBG_START;
	}
	DBGEdge * getEdgeTo(DBGNode * target,uchar edge_type=EDGE_OUT){
		Lnk<DBGEdge>* tmp_lnk=edges[edge_type]->lnk;
		while(tmp_lnk){
			if(tmp_lnk->pdat->nodes[1-edge_type]==target)return tmp_lnk->pdat;
			tmp_lnk=tmp_lnk->next;
		}
		return NULL;
	}
	DBGEdge * getMaxEdge(uchar edge_type=EDGE_OUT){
		if(edges[edge_type]->length==0)return NULL;
		Lnk<DBGEdge>* tmp_lnk=edges[edge_type]->lnk;
		Lnk<DBGEdge>* max_lnk=NULL;
		int max_w=0;
		//cout<<"\nMaxEdgeSearch:"<<endl;
		while(tmp_lnk){
			//tmp_lnk->pdat->nodes[1]->printKmer();
			//cout<<" "<<tmp_lnk->pdat->nodes[1]->getFreq()<<endl;
			if(tmp_lnk->pdat->W>max_w){
				max_w=tmp_lnk->pdat->W;
				max_lnk=tmp_lnk;
			}
			tmp_lnk=tmp_lnk->next;
		}
		if(max_lnk)return max_lnk->pdat;
		return NULL;
	}
	bool link(DBGNode * target,uchar link_type=EDGE_OUT,uint weight=1){
		//printKmer();
		//cout<<" --> ";
		//target->printKmer();
		DBGEdge * edge=getEdgeTo(target,link_type);
		if(edge!=NULL){
			//cout<<" "<<edge->W+1<<endl;
			edge->addWeight(weight);
			return 1;
		}else{
			//cout<<" Create"<<endl;
			if(link_type==EDGE_OUT){
				edge=new DBGEdge(this,target,weight);
			}else{
				edge=new DBGEdge(target,this,weight);
			}
			edges[link_type]->append(edge);
			target->edges[1-link_type]->append(edge);
			return 0;
		}
	}
	bool link(DBGEdge * target,uchar link_type=EDGE_OUT){
		DBGEdge * edge=getEdgeTo(target->nodes[1-link_type],link_type);	
		if(!edge){
			target->nodes[link_type]=this;
			edges[link_type]->append(target);
			return true;
		}else{
			return false;
		}
	}
	uchar unlink(DBGEdge * edge_target){
		if(edges[0]->remove(edge_target)){
			edge_target->nodes[0]=NULL;
			return 0;
		}
		if(edges[1]->remove(edge_target)){
			edge_target->nodes[1]=NULL;
			return 1;
		}
		return 2;
	}
	bool changeLinkTo(DBGEdge * edge,DBGNode * target){
		bool edge_type;
		if(edge->nodes[0]==this){
			edge_type=EDGE_IN;
		}
		else if(edge->nodes[1]==this){
			edge_type=EDGE_OUT;
		}else{
			return false;
		}
		print();
		cout<<"->ChangeLnkTo: ";
		target->print();
		edge->nodes[edge_type]->unlink(edge);
		target->link(edge,edge_type);
		return true;
	}
	DBGNode * createShadow(){
		DBGNode * shadow=new DBGNode;
		shadow->myTreeNode=myTreeNode;
		shadow->idx=idx;
		shadow->flag=DBG_START;
		shadow->reads=reads;
		shadow->kmers=kmers;
		return shadow;

	}
	void EulerSearch();
	void GreedySearch();
};
class NodeRecord{
public:
	DBGNode * node;
	DBGNode * tailNode;
	LnkHeader<DBGEdge> * seq_edges;
	uint idx_node;
	NodeRecord(DBGNode * n,uint idx){
		node=n;
		tailNode=NULL;
		idx_node=idx;
		seq_edges=new LnkHeader<DBGEdge>;
	}
	void append(DBGEdge * edge){
		seq_edges->append(edge);
	}
	void pop(){
		seq_edges->pop();
	}
	uint getLength(){
		return seq_edges->length;
	}
	uint getScore(){
		return seq_edges->length;
	}

};
class RingDestroyer{
public:
	LnkHeader<DBGNode> * shadows;
	LnkHeader<NodeRecord> * node_IN;
	LnkHeader<NodeRecord> * node_OUT;
	NodeRecord * record_current;
	uint idx_current;
	RingDestroyer(){
		shadows=NULL;
		node_IN=new LnkHeader<NodeRecord>;
		node_OUT=new LnkHeader<NodeRecord>;
		record_current=NULL;
		idx_current=0;
	}
	void makeRecordToAll(DBGNode * new_node,DBGEdge * from_new_edge,LnkHeader<NodeRecord> * record_list,bool make_tail=false){
		Lnk<NodeRecord> * record_tmp=record_list->lnk;
		while(record_tmp){
			cout<<"make record from ";
			record_tmp->pdat->node->printKmer();
			cout<<endl;
			record_tmp->pdat->append(from_new_edge);
			if(make_tail){
				cout<<"make tail"<<endl;
				record_tmp->pdat->tailNode=new_node;
			}
			record_tmp=record_tmp->next;
		}
	}
	void deleteEdge(DBGEdge * target){
		target->nodes[0]->unlink(target);
		target->nodes[1]->unlink(target);
		delete target;
	}
	void DestroyBubble(NodeRecord * n1,NodeRecord * n2){
		if(n1->getScore()>n2->getScore()){
			deleteEdge(n1->seq_edges->tail->pdat);
		}else{
			deleteEdge(n2->seq_edges->tail->pdat);
		}
	}
	void DestroyTip(NodeRecord * n){
		Lnk<DBGEdge> * edge_tmp=n->seq_edges->tail;
		while(edge_tmp){
			Lnk<DBGEdge> * edge_tmp2=edge_tmp;
			edge_tmp=edge_tmp->prev;
			delete edge_tmp2;
		}
	}
	void DestroyRepeat(NodeRecord * ring_back,NodeRecord * ring_forword,DBGEdge * edge_in,DBGEdge *edge_out){
		Lnk<DBGEdge> * lnk_edge=ring_forword->seq_edges->tail;
		DBGNode * shadow=NULL;
		DBGNode * shadow_next=NULL;
		DBGNode * shadow2=NULL;
		shadow_next=lnk_edge->pdat->nodes[0]->createShadow();
		while(lnk_edge){
			shadow=shadow_next;
			if(lnk_edge==ring_forword->seq_edges->tail){
				//第一个节点交叉建立
				//外节点接入,除了edge_in其它全部链接到shadow
				Lnk<DBGEdge> * edge_tmp=ring_forword->node->edges[EDGE_IN]->lnk;
				while(edge_tmp){
					if(edge_tmp->pdat!=edge_in){
						edge_tmp->pdat->nodes[0]->changeLinkTo(edge_tmp->pdat,shadow);
					}
					edge_tmp=edge_tmp->next;
				}
				//链接原通路下一节点
				shadow->link(lnk_edge->pdat->nodes[1]);
				//提前建立影子通路下一节点
				shadow2=lnk_edge->pdat->nodes[1]->createShadow();
				//链接影子通路下一节点。
				lnk_edge->pdat->nodes[0]->changeLinkTo(lnk_edge->pdat,shadow2);
			}else{
				shadow2=lnk_edge->pdat->nodes[1]->createShadow();
				shadow->link(shadow2);
			}
			shadow_next=shadow2;
			lnk_edge=lnk_edge->prev;
		//ring_forword->seq_edges->lnk->
			//第二个节点不用新建，因为第一个节点临时新建了。
		}
		//结尾接出
		Lnk<DBGEdge> * edge_tmp=ring_back->node->edges[EDGE_OUT]->lnk;
		while(edge_tmp){
			if(edge_tmp->pdat!=edge_out){
				edge_tmp->pdat->nodes[1]->changeLinkTo(edge_tmp->pdat,shadow_next);
			}
			edge_tmp=edge_tmp->next;
		}
	}
	bool DFSToTargetNode(DBGNode * begin_node,DBGEdge * from_the_edge,DBGNode * target_node,NodeRecord * record,uint depth,bool ignoreVisitFlag=true){
		if((!ignoreVisitFlag) && begin_node->flag!=DBG_START)return false;
		if(from_the_edge)record->append(from_the_edge);
		if(begin_node==target_node){
			record->tailNode=target_node;
			return true;
		}else{
			if(depth==0){
			}else{
				Lnk<DBGEdge>* lnk_edges=begin_node->edges[EDGE_OUT]->lnk;
				while(lnk_edges){
					if(DFSToTargetNode(lnk_edges->pdat->nodes[1],lnk_edges->pdat,target_node,record,depth-1)){
						return true;
					}
					lnk_edges=lnk_edges->next;
				}
			}
		}
		record->pop();
		return false;
	}
	void DFS(DBGNode * node,uint node_ID,DBGEdge * from_the_edge){
		//无论是否已经访问，先计算。
		cout<<"\nDFS to node["<<node_ID<<"]:\t";
		node->printKmer();
		cout<<endl;
		uint n_edge_out=node->edges[EDGE_OUT]->length;
		uint n_edge_in=node->edges[EDGE_IN]->length;
		cout<<"OUT:\t"<<n_edge_out<<"\tIN:\t"<<n_edge_in<<endl;
		bool node_IN_recorded=false;
		//是否为一个平凡点(1入1出)
		bool ordinary_node_recorded=false;
		if(n_edge_out<=1 && n_edge_in<=1){
			ordinary_node_recorded=true;
			cout<<"Ordinary node"<<endl;
		}
		//平凡点将延续当前记录,否则将视为当前记录的结尾
		if(ordinary_node_recorded){
			makeRecordToAll(NULL,from_the_edge,node_OUT,false);
		}else{
			makeRecordToAll(node,from_the_edge,node_OUT,false);
		}
		int depth=node->idx-record_current->node->idx;
		if(depth>0){
			//bubble
			NodeRecord * lastOutRecord=node_OUT->lnk->pdat;
			if(DFSToTargetNode(lastOutRecord->node,NULL,node,new_record,depth,false)){
				DestroyBubble(lastOutRecord,new_record);
			}else{

			}
		}
		if(node->flag!=DBG_START){
			
			NodeRecord * new_record=new NodeRecord(node,node->idx);
			//一条序列终结于与已访问节点的碰撞时，这里包含着很大的信息量呢。

			if(depth<0){
				//repeat
				DFSToTargetNode(node,NULL,node_OUT->lnk->pdat->node,new_record,-depth);
				DestroyRepeat(record_current,new_record,from_the_edge,record_current->seq_edges->tail->pdat);
			}
			return;
		}
		node->flag=DBG_VISITING;
		node->idx=node_ID;

		if(n_edge_out>1 && n_edge_in>1){
			//交叉点
		}else{
			if(n_edge_out>1){
				//分离点
				record_current=new NodeRecord(node,node_ID);
				node_OUT->append(record_current);
			}
			if(n_edge_in>1){
				//汇合点,如果未记录分离点，就不记录汇合点。
				if(node_OUT->length>0){
					record_current=new NodeRecord(node,node_ID);
					node_IN->append(record_current);
					node_IN_recorded=true;
				}
			}

		}

		if(node_IN_recorded){
			//如果这是一个被记录的汇合点且则不能再往前移动。
			return;
		}
		Lnk<DBGEdge>* lnk_edges=node->edges[EDGE_OUT]->lnk;
		while(lnk_edges){
			DFS(lnk_edges->pdat->nodes[1],node_ID+1,lnk_edges->pdat);
			lnk_edges=lnk_edges->next;
		}
		//结尾的时候别忘记把record pop掉。因为每次append都在头部，所以只要pop第一个。
		//我这个链表能当栈用，牛逼吧。
		if(record_current){
			record_current->pop();
		}
	}
};
#endif