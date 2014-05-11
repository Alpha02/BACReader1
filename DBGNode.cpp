#include "main.h"
#include "DBGNode.h"
#include "bpTree.h"
extern uint KMER_LENGTH;
void DBGNode::GreedySearch(){
	if(flag!=DBG_START)return;
	flag=DBG_VISITING;
	if(myTreeNode)myTreeNode->refreshWinner();
	printKmer(KMER_LENGTH-1,1);
	DBGEdge * targetEdge=getMaxEdge(EDGE_OUT);
	if(targetEdge){
		targetEdge->nodes[1]->GreedySearch();
	}
	flag=DBG_OVER;
}
void DBGNode::EulerSearch(){
	cout<<"visit:";
	printKmer();
	cout<<endl;
	flag=DBG_VISITING;
	Lnk<DBGEdge>* tmp_lnk=edges[EDGE_OUT]->lnk;
	bool ended=true;
	while(tmp_lnk){
		
		if(tmp_lnk->pdat->nodes[1]->flag==DBG_START){
			tmp_lnk->pdat->nodes[1]->EulerSearch();
			ended=false;
		}else{
			cout<<"visitFailedTo ";
			tmp_lnk->pdat->nodes[1]->printKmer();
			cout<<endl;
		}
		tmp_lnk=tmp_lnk->next;
	}
//	if(ended==true){
	cout<<"\noutput:";
		printKmer();
		cout<<"<-";
//	}
}