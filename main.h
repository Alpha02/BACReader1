#ifndef MAIN_H
#define MAIN_H
#include <iostream>
#include <string.h>
#include <fstream>
#include <time.h>
using namespace std;
#define  uint unsigned int
#define  uchar unsigned char
#include "CONFIG.h"
#include "ERROR.h"

#define READ_LEN_ID 50
#define READ_LEN_SEQ 120
#define STA_STOP 0
#define STA_START 1
#define STA_PAUSE 2
#define STA_END 3
//未访问的点
#define DBG_START 0
//正在访问的点（到达过，但是还有出路可寻）
#define DBG_VISITING 1
//访问结束的点（到达过，且所有出路已遍历）
#define DBG_OVER 2
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
template <class T> 
class LnkHeader{
public:
	Lnk<T> *lnk;
	Lnk<T> *tail;
	uint length;
	LnkHeader(){
		length=0;
		lnk=NULL;
		tail=NULL;
	}
	void append(T*new_pdat){
		if(!lnk){
			lnk=new Lnk<T>(new_pdat);
			tail=lnk;
		}else{
			lnk->append(new_pdat);
			if(tail==lnk)tail=tail->next;
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
		}else{
			tail=NULL;
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
	Lnk<T> * get(T * target){
		Lnk<T> * tmp_lnk=lnk;
		while(tmp_lnk){
			if(tmp_lnk->pdat==target)return tmp_lnk;
			tmp_lnk=tmp_lnk->next;
		}
		return NULL;
	}
	bool remove(T * target){
		Lnk<T> * tmp_lnk=get(target);
		if(tmp_lnk){
			if(tail==tmp_lnk){
				tail=tmp_lnk->prev;
			}
			tmp_lnk->remove();
			
			length--;
			if(length==0){
				lnk=NULL;
			}
			return true;
		}
		return false;
	}
};


#endif