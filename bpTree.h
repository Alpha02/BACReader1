#ifndef BPTREE_H
#define BPTREE_H
#include "DBGNode.h"
class bpTreeNode{
public:
	DBGNode * tail;
	bpTreeNode * winner_child;
	bpTreeNode * child[4];
	bpTreeNode * parent;
	bpTreeNode(bpTreeNode * p=NULL,DBGNode * t=NULL){
		child[0]=NULL;
		child[1]=NULL;
		child[2]=NULL;
		child[3]=NULL;
		parent=p;
		tail=t;
		winner_child=NULL;
	}
	uint getFreq(){
		return tail->getFreq();
	}
	bpTreeNode * getWinner(bool rechoose=false){
		if(tail!=NULL){
			return this;
		}
		if((!rechoose) && winner_child!=NULL){
			return winner_child;
		}
		for(uchar i=0;i<4;i++){
			if(!child[i])continue;
			bpTreeNode * winner_tmp=child[i]->getWinner(rechoose);
			if((!winner_child)||(winner_tmp && (winner_tmp->getFreq()>winner_child->getFreq()))){
				/*
				if(winner_child && winner_tmp){

					winner_tmp->tail->printKmer();
					cout<<" > ";
					winner_child->tail->printKmer();
					cout<<endl;
				}
				*/
				winner_child=winner_tmp;
			}
		}
		return winner_child;
	}
	void findEularStartNodes(LnkHeader<DBGNode> * nodelist){
		for(uchar i=0;i<4;i++){
			if(!child[i])continue;
			child[i]->findEularStartNodes(nodelist);
		}
	//如果出度-入度=1，则可以记录这个点。
		if(tail && tail->edges[EDGE_OUT]->length>tail->edges[EDGE_IN]->length){
			cout<<"\n newEularStartNode ";
			tail->printKmer();
			cout<<endl;
			nodelist->append(tail);
		}
}
	DBGNode * findEularStartNode(){
		bpTreeNode * tmpWinner=getWinner();
		if(tmpWinner->tail->edges[EDGE_OUT]->length>tmpWinner->tail->edges[EDGE_IN]->length){
			return tmpWinner->tail;
		}else{

		}
	}
	void refreshWinner(){
		if(!parent)return;
		for(uchar i=0;i<4;i++){
			if(!parent->child[i])continue;
			bpTreeNode * winner_tmp=parent->child[i]->winner_child;
			if(winner_tmp && (winner_tmp->getFreq()>parent->winner_child->getFreq())){
				parent->winner_child=winner_tmp;
			}
		}
		parent->refreshWinner();	
	}
	void print(){
		for(int i=0;i<4;i++){
			if(child[i]){
				child[i]->print();
			}
		}
		if(tail){
			tail->printKmer();
			cout<<" "<<tail->getFreq()<<endl;
		}
	}
	void startAll(){
		if(tail){
			tail->flag=DBG_START;
			
		}
		for(int i=0;i<4;i++){
			if(child[i]){
				child[i]->startAll();
			}
		}
	}
};
class ERROR_TreeNodeNotFound{};

class bpTree{
public:
	bpTreeNode * root;
	bpTreeNode * winner;
	uint number_DBGNodes;
	uint number_DBGEdges;
	uint weight_DBGEdges;
	uint weight_DBGNodes;
	DBGNode * EularStartNode;
	bpTree(){
		root=new bpTreeNode();
		winner=NULL;
		number_DBGNodes=0;
		number_DBGEdges=0;
		weight_DBGEdges=0;
		weight_DBGNodes=0;
		EularStartNode=NULL;
	}
	DBGNode * hangOnTree(bpRead * src){
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
		weight_DBGNodes++;
		if(!node->tail){
			node->tail=new DBGNode();
			node->tail->myTreeNode=node;
			//node->tail=new DBGNode(src);
			number_DBGNodes++;
		}else{
			//node->tail->append(src);
		}
		return node->tail;
	}
	void appendDBGraph(bpRead * src){
		src->pos_cur=0;
		src->pos_last=0;
		src->start();
		DBGNode * node=NULL;
		DBGNode * node_next=NULL;
		while(src->status==STA_START){
			node_next=hangOnTree(src);
			node_next->appendKmer(src);
			if(node){
				weight_DBGEdges++;
				if(node->link(node_next,EDGE_OUT,1))number_DBGEdges++;
			}
			node=node_next;
			src->pushForword();
		}
		src->pos_cur=0;
		src->pos_last=0;
		src->start();

		
	}
	void chooseEularStartNode(DBGNode * new_node){
		if(EularStartNode && EularStartNode->edges[EDGE_IN]->length==0)return;
		//如果出度-入度=1，则可以开始这个点。
		if(new_node && new_node->edges[EDGE_OUT]->length-new_node->edges[EDGE_IN]->length==1){
			cout<<"\n newEularStartNode is ";
			new_node->printKmer();
			cout<<endl;
			EularStartNode=new_node;
		}
	}


	LnkHeader<bpRead> * catchAllReadFromNode(bpTreeNode * root,LnkHeader<bpRead> * fishingRod){
		//递归找到节点下所有的Read
		for(int i=0;i<4;i++){
			if(root->child[i]==NULL)continue;
			catchAllReadFromNode(root->child[i],fishingRod);
		}
		if(root->tail!=NULL){
			fishingRod->appendFromLnk(root->tail->getGetFirstRead());
			cout<<"\napp:";
			root->tail->getSeq()->print();
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
	DBGNode * popLeavesLnk(bpTreeNode * target){
		DBGNode* lnk=target->tail;
		return lnk;
	}
	LnkHeader<bpRead> * search(bpString * seq,uint begin_idx=0,uint end_idx=KMER_LENGTH-1){
		//寻找read根据seq,结果为reads的链表。为所有可行的寻找结果。
		bpTreeNode * node=searchNode(seq,begin_idx,end_idx);
		LnkHeader<bpRead> * fishingRod=new LnkHeader<bpRead>;
		catchAllReadFromNode(node,fishingRod);
		cout<<"Search result:"<<fishingRod->length<<" found."<<endl;
	}
};
#endif