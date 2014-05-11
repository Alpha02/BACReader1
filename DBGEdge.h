#ifndef DBGEDGE_H
#define DBGEDGE_H
class DBGNode;
class DBGEdge{
public:
	DBGNode * nodes[2];
	int W;
	DBGEdge (){
		nodes[0]=NULL;
		nodes[1]=NULL;
		W=0;
	}
	DBGEdge(DBGNode * node_from,DBGNode * node_to,int new_w=1){
		nodes[0]=node_from;
		nodes[1]=node_to;
		W=new_w;
	}
	void addWeight(int dW){
		W+=dW;
	}
};
#endif